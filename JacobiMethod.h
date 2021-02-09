//
// Created by Admin on 10.05.2020.
//

#ifndef OPPLAB4_JACOBIMETHOD_H
#define OPPLAB4_JACOBIMETHOD_H

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <algorithm>

double phi(double x, double y, double z);

double ro(double x, double y, double z, const double a);

double updLayer(int base_z, int height, double *omega_part, double *tmp_omega_part,
                double hx, double hy, double hz,
                const int N, const double x0, const double y0, const double z0, const double a);

double JacobiMethod(const double epsilon, const double a, const int N,
                    const double x0, const double y0, const double z0,
                    const double x1, const double y1, const double z1);

#endif //OPPLAB4_JACOBIMETHOD_H
