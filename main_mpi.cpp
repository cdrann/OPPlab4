#include <mpi.h>

#include "Test.h"

// [область моделирования]: [-1; 1] x [-1; 1] x [-1; 1]
#define x0 -1
#define y0 -1
#define z0 -1
#define x1 1
#define y1 1
#define z1 1

#define repeats 1
#define a 1e5
#define N 240 // 240
#define epsilon 1e-8

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    /*
     * ИСХОДНЫЕ ДАННЫЕ ЗАДАЧИ
     * [область моделирования] simulation area: [-1; 1] x [-1; 1] x [-1; 1]
     * [искомая функция] desired function: phi(x, y, z) = x^2 + y^2 + z^2
     * [правая часть уравнения] right side of the equation: ro(x, y, z) = 6 - a * phi(x, y, z)
     * [параметр уравнения] equation parameter: a = 10^5 (1e5)
     * [порог сходимости] convergence threshold: epsilon = 10^-8 (1e-8)
     * [начальное приближение] initial approximation: phi_{i, j, k}^0 = 0
     *
     * */

    JacobiMethodTest(repeats, epsilon, a, N, x0, y0, z0, x1, y1, z1);

    MPI_Finalize();
    return 0;
}