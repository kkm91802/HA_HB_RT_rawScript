#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// -----------------------------
// Constants
// -----------------------------
#define N_ORIENT 6
#define G_CONST 9.8

void apply_calibration(const char *datafile,
                       double Sa[3][3],
                       double Ma[3][3],
                       double ba[3]);

// -----------------------------
// Matrix utilities
// -----------------------------
void matvec3(double M[3][3], double v[3], double result[3]);
void matmul3x3(double A[3][3], double B[3][3], double result[3][3]);
// bool invert3x3(double M[3][3], double invOut[3][3]);
double invert3x3(double M[3][3], double invOut[3][3]);

#endif // CALIBRATION_H
