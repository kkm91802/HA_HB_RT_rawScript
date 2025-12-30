#ifndef CALIBRATION_METRICES_H
#define CALIBRATION_METRICES_H

#include <stdio.h>

extern double Sa[3][3];
extern double Ma[3][3];
extern double ba[3];

#define N_ORIENT 6
#define G_CONST 9.8
#define MAX_LINE 256
#define MAX_COMB 30

typedef struct {
    double kx, ky, kz;
    double bx, by, bz;
    double alpha_yz, alpha_zy, alpha_zx;
} Params;

typedef struct {
    double meanAx, meanAy, meanAz;
} MeanReading;

// === Core Functions ===

// Reads a CSV and computes mean Ax, Ay, Az
int read_csv_means(const char *filename, MeanReading *out);

// Generates and solves combinations per the paper
Params calibrate(double V[6][3], double G[6][3]);

// Performs the full calibration process for a folder path
void run_calibration(const char *basepath);

#endif
