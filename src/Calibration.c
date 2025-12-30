#include <stdio.h>      
#include <stdlib.h>     
#include <string.h>     
#include <math.h>       
#include "CalibrationMetrices.h"
#include <stdbool.h>
#include "Calibration.h"

// --- Helper Functions ---
void matvec3(double M[3][3], double v[3], double result[3]) {
    for (int i = 0; i < 3; i++)
        result[i] = M[i][0]*v[0] + M[i][1]*v[1] + M[i][2]*v[2];
}

void matmul3x3(double A[3][3], double B[3][3], double result[3][3]) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            result[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
                result[i][j] += A[i][k] * B[k][j];
        }
}

double invert3x3(double M[3][3], double Minv[3][3]) {
    double det = 
        M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) -
        M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) +
        M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    if (fabs(det) < 1e-9) {
        // Not invertible, return identity instead
        printf("⚠️ Warning: Singular matrix (det ≈ 0), using identity.\n");
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Minv[i][j] = (i == j) ? 1.0 : 0.0;
        return 0.0;
    }

    double invdet = 1.0 / det;

    Minv[0][0] =  (M[1][1]*M[2][2] - M[1][2]*M[2][1]) * invdet;
    Minv[0][1] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) * invdet;
    Minv[0][2] =  (M[0][1]*M[1][2] - M[0][2]*M[1][1]) * invdet;

    Minv[1][0] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) * invdet;
    Minv[1][1] =  (M[0][0]*M[2][2] - M[0][2]*M[2][0]) * invdet;
    Minv[1][2] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) * invdet;

    Minv[2][0] =  (M[1][0]*M[2][1] - M[1][1]*M[2][0]) * invdet;
    Minv[2][1] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) * invdet;
    Minv[2][2] =  (M[0][0]*M[1][1] - M[0][1]*M[1][0]) * invdet;

    return det;
}

void apply_calibration(const char *datafile,
                       double Sa[3][3],
                       double Ma[3][3],
                       double ba[3])
{
    FILE *fp = fopen(datafile, "r");
    if (!fp) {
        printf("Error: cannot open %s\n", datafile);
        return;
    }

    FILE *out = fopen("/home/kjw2kor/shared_folder/csv/Calibration.csv", "w");
    if (!out) {
        printf("Error: cannot create output file.\n");
        fclose(fp);
        return;
    }

    char line[256];
    char Time[32];
    double Ax, Ay, Az, Gx, Gy, Gz;
    double phi, theta;
    double Rxyz[3][3], Ryxz[3][3], R[3][3];

    // Read and skip CSV header
    fgets(line, sizeof(line), fp);
    fprintf(out, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");

    int line_num = 0, valid_lines = 0;

    while (fgets(line, sizeof(line), fp)) {
        line_num++;

        // Read timestamp as string (%[^,]) and 6 numeric values
        int n = sscanf(line, "%[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                       Time, &Ax, &Ay, &Az, &Gx, &Gy, &Gz);

        if (n != 7) {
            printf("Skipped line %d: %s", line_num, line);
            continue;
        }

        valid_lines++;

        // Compute phi and theta
        phi = atan2(Gy, sqrt(Gx * Gx + Gz * Gz));
        theta = atan2(-Gx, Gz);

        // Define Rxyz
        Rxyz[0][0] = cos(theta);  Rxyz[0][1] = 0;              Rxyz[0][2] = -sin(theta);
        Rxyz[1][0] = sin(phi)*sin(theta); Rxyz[1][1] = cos(phi); Rxyz[1][2] = sin(phi)*cos(theta);
        Rxyz[2][0] = cos(phi)*sin(theta); Rxyz[2][1] = -sin(phi); Rxyz[2][2] = cos(phi)*cos(theta);

        // Define Ryxz
        Ryxz[0][0] = cos(theta);  Ryxz[0][1] = sin(theta)*sin(phi);  Ryxz[0][2] = -sin(theta)*cos(phi);
        Ryxz[1][0] = 0;           Ryxz[1][1] = cos(phi);             Ryxz[1][2] = sin(phi);
        Ryxz[2][0] = sin(theta);  Ryxz[2][1] = -cos(theta)*sin(phi); Ryxz[2][2] = cos(theta)*cos(phi);

        // Choose between Rxyz and Ryxz based on condition
        double left = fabs(sqrt(Gy * Gy + Gz * Gz )-1);
        double right = fabs(sqrt(Gx * Gx + Gz * Gz)-1);
        if (left < right)
            memcpy(R, Rxyz, sizeof(R));
        else
            memcpy(R, Ryxz, sizeof(R));

        // Compute total = I + Sa + Ma
        double I[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
        double total[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                total[i][j] = I[i][j] + Sa[i][j] + Ma[i][j];

        // Prepare for calibration math
        double a_tilda[3] = {Ax, Ay, Az};
        double Ka[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        double KaAtilda[3], Rinverse[3][3], totalInv[3][3], temp2[3][3], a_cal[3];

        // Matrix operations
        matvec3(Ka, a_tilda, KaAtilda);
        for (int i = 0; i < 3; i++)
            KaAtilda[i] -= ba[i];



        // Compute determinant to check invertibility
        double detR = R[0][0]*(R[1][1]*R[2][2] - R[1][2]*R[2][1]) -
              R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0]) +
              R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0]);

        double detTotal = total[0][0]*(total[1][1]*total[2][2] - total[1][2]*total[2][1]) -
                  total[0][1]*(total[1][0]*total[2][2] - total[1][2]*total[2][0]) +
                  total[0][2]*(total[1][0]*total[2][1] - total[1][1]*total[2][0]);

        if (fabs(detR) < 1e-6 || fabs(detTotal) < 1e-6) {
                 printf("⚠️  Skipped line %d: nearly singular matrix (detR=%.6e detTotal=%.6e)\n",
                 line_num, detR, detTotal);
        continue; // skip this record to avoid NaN
        }


        invert3x3(R, Rinverse);
        invert3x3(total, totalInv);
        matmul3x3(totalInv, Rinverse, temp2);
        matvec3(temp2, KaAtilda, a_cal);

        // Write calibrated result
        fprintf(out, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                Time, a_cal[0], a_cal[1], a_cal[2], Gx, Gy, Gz);
    }

    fclose(fp);
    fclose(out);
    // remove(datafile);
    // rename("data_calibrated.csv", datafile);
    printf("Calibration applied successfully.\n");
    printf("Processed %d lines, valid %d lines.\n", line_num, valid_lines);
}
