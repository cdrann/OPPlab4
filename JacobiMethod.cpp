//
// Created by Admin on 10.05.2020.
//

#include "JacobiMethod.h"

double phi(double x, double y, double z) {
    // [искомая функция] desired function: phi(x, y, z) = x^2 + y^2 + z^2
    return pow(x, x) + pow(y, y) + pow(z, z);
}

double ro(double x, double y, double z, const double a) {
    // [правая часть уравнения] right side of the equation: ro(x, y, z) = 6 - a * phi(x, y, z)
    return 6 - a * phi(x,y,z);
}

double updLayer(
        int base_z, int height, double *omega_part, double *tmp_omega_part,
        double hx, double hy, double hz,
        const int N, const double x0, const double y0, const double z0, const double a) {

    // модуль
    int abs_z = base_z + height;

    // если граница области
    if (abs_z == 0 || abs_z == N - 1) {
        // Копируем этот слой в новый массив на старое место, не пересчитывая
        memcpy(tmp_omega_part + height * N * N, omega_part + height * N * N, N * N * sizeof(double));
        return 0;
    }

    //Иначе пересчитываем каждый _элемент_ слоя с помощью итерационной формулы
    double max_delta = 0;
    double z = z0 + abs_z * hz;

    for (int i = 0; i < N; i++) {
        double x = x0 + i * hx;

        for (int j = 0; j < N; j++) {
            double y = y0 + j * hy;

            // номер клетки в слое
            int cell = height * N * N + i * N + j;

            //Если элемент находится на границе слоя, то не пересчитываем его
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                tmp_omega_part[cell] = omega_part[cell];
                continue;
            }

            // иначе пересчитываем, по формуле
            // Jacobi formula (2):
            // {[phi^m(i+1, j, k) - phi^m(i-1, j, k)] / hx^2
            // + [phi^m(i, j+1, k) - phi^m(i, j-1, k)] / delta_y^2
            // + [phi^m(i, j, k+1) - phi^m(i, j, k-1)] / delta_z^2
            // - ro(x, y, z)}
            // / ([2 / hx^2] + [2 / delta_y^2] + [2 / delta_z^2])

            tmp_omega_part[cell] = ((omega_part[height * N * N + (i + 1) * N + j]
                                    + omega_part[height * N * N + (i - 1) * N + j]) / (hx * hx)
                                    + (omega_part[height * N * N + i * N + (j + 1)]
                                    + omega_part[height * N * N + i * N + (j - 1)]) / (hy * hy)
                                    + (omega_part[(height + 1) * N * N + i * N + j]
                                    + omega_part[(height - 1) * N * N + i * N + j]) / (hz * hz)
                                    - ro(x, y, z, a)) /
                                            ((2 / (hx * hx)) + (2 / (hy * hy)) + (2 / (hz * hz)) + a);

            max_delta = std::max(max_delta, std::abs(tmp_omega_part[cell] - omega_part[cell]));
        }
    }

    return max_delta;
}

double JacobiMethod(const double epsilon, const double a, const int N,
                    const double x0, const double y0, const double z0,
                    const double x1, const double y1, const double z1) {

     /* Требуется решить уравнение (т.е найти ф-ю phi):
     * d^2(phi)/d^2(x) + d^2(phi)/d^2(y) + d^2(phi)/d^2(z) - a*phi = ro, [a >= 0] */

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (N % size) {
        if(rank == 0) {
            std::cout << "Invalid number of processes" << std::endl;
            return -1;
        }
    }

    double start_time = MPI_Wtime();

    //расстояния между соседними узлами - шаги сетки
    double hx = (x1 - x0) / (N - 1);
    double hy = (y1 - y0) / (N - 1);
    double hz = (z1 - z0) / (N - 1);

    // высота слоя
    int part_height = N / size;
    // координата для текущего процесса
    int part_base_z = rank * part_height - 1;

    // трехмерный массив для сохранения значений в точках слоя
    // +2 -- верхнее и нижнее внутреннее граничное
    double *omega = new double[(part_height + 2) * N * N]; //values
    double *tmp_omega = new double[(part_height + 2) * N * N]; //tmpvalues

    int iterationsCounter = 0;

    // _Шаг алгоритма 1_
    // Задать значения искомой функции на границе области omega: phi_{i,j,k} = F(x_i, y_j, z_k)
    // при i = 0, i = N_x, j = 0, j = N_y, k = 0, k = N_z
    // _Шаг алгоритма 2_
    // Задать начальное приближение во внутренней части области omega: phi_{i,j,k}^0
    // для i=1..Nx-1, j=0..Ny-1, k=0..Nz-1.
    for (int i = 0; i < part_height + 2; i++) {
        // Oz области омега
        int omega_z = i + part_base_z;
        double real_z = z0 + hz * omega_z;

        for (int j = 0; j < N; j++) {
            // Ox
            double x = x0 + hx * j;

            for (int k = 0; k < N; k++) {
                // Oy
                double y = y0 + hy * k;

                if (omega_z == 0 || omega_z == N - 1 || j == 0 || j == N - 1 || k == 0 || k == N - 1) {
                    // значение функции на границе области [1 шаг]
                    omega[i * N * N + j * N + k] = phi(x, y, real_z);
                } else {
                    // начальное приближение во внутренней части области [2 шаг]
                    omega[i * N * N + j * N + k] = 0;
                }
            }
        }
    }


    // _Шаг алгоритма 3_
    // Многократно вычисляем очередное приближение искомой функции по формуле phi_{i, j, k}^(m+1)
    // пока не достигнуто условие: max|phi_{i, j, k}^(m+1) - phi_{i, j, k}^(m)| < eps. [по i, j, k]
    double max_delta_shared; // порог сх-ти (3)
    do {
        // 1. вычисляются сеточные значение, прилегающие к границе локальной подобласти

        //
        //BOUNDARY LAYER:
        //

        // boundary top
        double max_delta = 0;
        double tmp_delta = updLayer(part_base_z, 1, omega, tmp_omega, hx, hy, hz, N, x0, y0, z0, a);
        max_delta = std::max(max_delta, tmp_delta);

        // boundary bottom
        tmp_delta = updLayer(part_base_z, part_height, omega, tmp_omega, hx, hy, hz, N, x0, y0, z0, a);
        max_delta = std::max(max_delta, tmp_delta);


        // 2. запускается асинхронный обмен граничных значений
        MPI_Request rq[4];

        // send[recv] TOP boundary [outer] layer to[from] prev process
        if (rank != 0) {
            MPI_Isend(tmp_omega + N * N, N * N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &rq[0]);
            MPI_Irecv(tmp_omega, N * N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &rq[2]);
        }

        //send[recv] BOTTOM boundary [outer] layer to[from] next process
        if (rank != size - 1) {
            MPI_Isend(tmp_omega + part_height * N * N, N * N, MPI_DOUBLE,
                    rank + 1, 0, MPI_COMM_WORLD, &rq[1]);
            MPI_Irecv(tmp_omega + (part_height + 1) * N * N, N * N, MPI_DOUBLE,
                    rank + 1, 0, MPI_COMM_WORLD, &rq[3]);
        }

        //3. выполняется вычисление остальных точек подобласти
        //
        // NON-BOUNDARY LAYER:
        //

        // re-calculating all elements of non-boundary layers
        for (int i = 2; i < part_height; i++) {
            double tmpdelta = updLayer(part_base_z, i, omega, tmp_omega, hx, hy, hz, N, x0, y0, z0, a);
            max_delta = std::max(max_delta, tmpdelta);
        }

        // 4. ожидание завершения обменов
        if (rank != 0) {
            MPI_Wait(&rq[0], MPI_STATUS_IGNORE);
            MPI_Wait(&rq[2], MPI_STATUS_IGNORE);
        }

        if (rank != size - 1) {
            MPI_Wait(&rq[1], MPI_STATUS_IGNORE);
            MPI_Wait(&rq[3], MPI_STATUS_IGNORE);
        }

        // fully recalculated area for this process
        memcpy(omega, tmp_omega, (part_height + 2) * N * N  * sizeof(double));

        // maximum delta of ALL processes & send to all
        MPI_Reduce(&max_delta, &max_delta_shared, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_delta_shared, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        iterationsCounter++;
    } while (max_delta_shared >= epsilon);
    // по достижению некоторого порога сх-ти -- завершение итерационного процесса

    delete[] tmp_omega; // results

    double *fullResult = nullptr;
    if (rank == 0) {
        fullResult = new double[N * N * N];
    }

    // gathering counted area-parts
    MPI_Gather(omega + N * N, part_height * N * N, MPI_DOUBLE, fullResult,
            part_height * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    // ищем максимальное отклонение max |phi_{i,j,k}^m - phi_{i,j,k}^*|
    if (rank == 0) {
        // counting maximum deviation from the standard response
        double max_delta = 0;

        for (int layer = 0; layer < N; layer++){
            double z = z0 + layer * hz;

            for (int j = 0; j < N; j++) {
                double x = x0 + j * hx;

                for (int k = 0; k < N; k++) {
                    double y = y0 + k * hy;

                    max_delta = std::max(max_delta, std::abs(fullResult[layer * N * N + j * N + k] - phi(x, y, z)));
                }
            }
        }

        // оценка точности полученного решения:
        std::cout << "Answer: delta = " << max_delta << std::endl;
        printf("Time: %lf\n", end_time - start_time);
        std::cout << iterationsCounter << " cycle iterations" << std::endl;

        delete[] fullResult;
    }

    delete[] omega;

    return end_time - start_time;
}

