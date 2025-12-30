#ifndef INEKF_H
#define INEKF_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NSTATE 15
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)
#define GGRAV 9.8
#define IDX(i, j, n) ((i) * (n) + (j))

typedef struct {
    double R[9];    // Rotation matrix
    double v[3];    // Velocity
    double p[3];    // Position
    double bw[3];   // Gyro bias
    double ba[3];   // Accel bias
    double P[NSTATE*NSTATE]; // Covariance
} InEKFState;

typedef struct {
    double gyro_noise_sigma;
    double accel_noise_sigma;
    double gyro_bias_rw;
    double accel_bias_rw;
    double zupt_velocity_noise;
    double zupt_accel_threshold;
    double zupt_gyro_threshold_deg;
} InEKFConfig;

typedef struct {
    InEKFState state;
    InEKFConfig config;
    double Qc[NSTATE*NSTATE]; // Process noise
    double R_zupt[9];         // ZUPT measurement noise
} InEKF;

// Function declarations
void inekf_init(InEKF *filter, const InEKFConfig *config);
void inekf_init_default(InEKF *filter);
void inekf_predict(InEKF *filter, const double omega_meas[3], const double acc_meas[3], double dt);
void inekf_zupt_update(InEKF *filter, const double velocity_meas[3]);
int inekf_detect_zupt(const InEKF *filter, const double acc_meas[3], const double gyro_meas_deg[3]);
int run_InEKF(const char* input_file, const char* output_file) ;

#endif