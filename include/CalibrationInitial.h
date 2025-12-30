#ifndef CALIBRATIONINITIAL_H
#define CALIBRATIONINITIAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINE_BUFFER_SIZE 256
#define CALIBRATION_SAMPLES 100   // Number of initial samples used for calibration

// IMU data structure with string timestamp
typedef struct {
    char time[32];     // Time in HH:MM:SS:MS format
    float Ax, Ay, Az;  // Accelerometer data
    float Gx, Gy, Gz;  // Gyroscope data
    int is_stationary; // 1 if used for calibration, 0 otherwise
} IMUData;

// Function prototypes
float compute_magnitude(float x, float y, float z);
int count_rows(const char *filename);
int read_IMU_data(const char *filename, IMUData **data);
void compute_bias(const IMUData *data, int row_count, float acc_bias[3], float gyro_bias[3]);
void apply_calibration2(IMUData *data, int row_count, const float acc_bias[3], const float gyro_bias[3]);
int write_calibrated_data(const char *output_filename, const IMUData *data, int row_count);

#endif
