//
// Created by Admin on 10.05.2020.
//

#ifndef OPPLAB4_TEST_H
#define OPPLAB4_TEST_H

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <algorithm>

void JacobiMethodTest(const int repeats, const double epsilon, const double a, const int N,
                      const double x0, const double y0, const double z0,
                      const double x1, const double y1, const double z1);

#endif //OPPLAB4_TEST_H
