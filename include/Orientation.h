#ifndef ORIENTATION_H
#define ORIENTATION_H

#include <stdio.h>

// define constants
#ifdef LINE_SIZE
#undef LINE_SIZE
#endif
#define LINE_SIZE 4096

// Function to apply transformation based on orientation
void transform_row(double *Ax, double *Ay, double *Az, double *Gx, double *Gy, double *Gz);

// Function to transform an entire CSV file in-place
void process_csv_file(const char *filename);

#endif
