#ifndef AUTOCALIBRATION_H
#define AUTOCALIBRATION_H

#include <stdio.h>

// Global calibration parameters
extern double Ox, Oy, Oz; // biases
extern double Sxx, Syy, Szz, Sxy, Sxz, Syz; // symmetric scale matrix

// Function prototypes
int run_autocalibration(const char *folder_path, int start_idx, int end_idx,
                        int window_size, double std_frac);

// void apply_accel_calibration(double Ax, double Ay, double Az,
//                              double *Ax_cal, double *Ay_cal, double *Az_cal);

static void compute_tilt_angles_rad(double Ax, double Ay, double Az,
                                    double *phi_rad, double *rho_rad);

static void build_R_xyz(double phi_rad, double theta_rad, double R[9]);
static void build_R_yxz(double phi_rad, double theta_rad, double R[9]);

int apply_calibration_to_csv(const char *input_csv, const char *output_csv);

#endif // AUTOCALIBRATION_H
