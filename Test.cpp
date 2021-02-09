//
// Created by Admin on 10.05.2020.
//

#include "Test.h"
#include "JacobiMethod.h"

void JacobiMethodTest(const int repeats, const double epsilon, const double a, const int N,
                      const double x0, const double y0, const double z0,
                      const double x1, const double y1, const double z1) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) {
        /* std::cout << "Initial data:" << std::endl;
        std::cout << "N = " << N << std::endl;
        std::cout << "epsilon = " << epsilon << std::endl;
        std::cout << std::endl; */
    }

    double current_time, best_time = std::numeric_limits<double>::max();

    for (int i = 1; i <= repeats; i++) {
        if (rank == 0) {
            std::cout << "Try " << i << "/" << repeats << std::endl;
        }

        current_time = JacobiMethod(epsilon, a, N, x0, y0, z0, x1, y1, z1);

        if (rank == 0) {
            best_time = (current_time < best_time) ? current_time : best_time;
        }


    }
    if (rank == 0) {
        printf("Best time: %lf sec\n", best_time);
    }
}
