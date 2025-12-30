#include "CalibrationMetrices.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// === Global Calibration Matrices ===
double Sa[3][3];
double Ma[3][3];
double ba[3];

// ---------- CSV Mean Computation ----------
int read_csv_means(const char *filename, MeanReading *out) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error opening %s\n", filename);
        return 0;
    }

    char line[MAX_LINE];
    double sumAx = 0, sumAy = 0, sumAz = 0;
    long count = 0;

    fgets(line, sizeof(line), fp);  // skip header

    while (fgets(line, sizeof(line), fp)) {
        char *token;
        double values[7];
        int i = 0;
        token = strtok(line, ",");
        while (token && i < 7) {
            values[i++] = atof(token);
            token = strtok(NULL, ",");
        }
        if (i >= 4) {
            sumAx += values[1];
            sumAy += values[2];
            sumAz += values[3];
            count++;
        }
    }
    fclose(fp);

    if (count == 0) return 0;

    out->meanAx = sumAx / count;
    out->meanAy = sumAy / count;
    out->meanAz = sumAz / count;
    return 1;
}

// ---------- Combination Generators ----------
int comb2(int (*pairs)[2]) {
    int c = 0;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++)
            pairs[c][0] = i, pairs[c++][1] = j;
    return c;
}
int comb3(int (*triplets)[3]) {
    int c = 0;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++)
            for (int k = j + 1; k < 6; k++)
                triplets[c][0] = i, triplets[c][1] = j, triplets[c++][2] = k;
    return c;
}
int comb4(int (*quads)[4]) {
    int c = 0;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++)
            for (int k = j + 1; k < 6; k++)
                for (int l = k + 1; l < 6; l++)
                    quads[c][0] = i, quads[c][1] = j, quads[c][2] = k, quads[c++][3] = l;
    return c;
}

// ---------- Matrix Solvers ----------
int solve_3x3(double A[3][3], double b[3], double x[3]) {
    double M[3][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) M[i][j] = A[i][j];
        M[i][3] = b[i];
    }
    for (int i = 0; i < 3; i++) {
        double pivot = M[i][i];
        if (fabs(pivot) < 1e-12) return 0;
        for (int j = i; j < 4; j++) M[i][j] /= pivot;
        for (int k = 0; k < 3; k++) {
            if (k == i) continue;
            double f = M[k][i];
            for (int j = i; j < 4; j++) M[k][j] -= f * M[i][j];
        }
    }
    for (int i = 0; i < 3; i++) x[i] = M[i][3];
    return 1;
}
int solve_4x4(double A[4][4], double b[4], double x[4]) {
    double M[4][5];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) M[i][j] = A[i][j];
        M[i][4] = b[i];
    }
    for (int i = 0; i < 4; i++) {
        double pivot = M[i][i];
        if (fabs(pivot) < 1e-12) return 0;
        for (int j = i; j < 5; j++) M[i][j] /= pivot;
        for (int k = 0; k < 4; k++) {
            if (k == i) continue;
            double f = M[k][i];
            for (int j = i; j < 5; j++) M[k][j] -= f * M[i][j];
        }
    }
    for (int i = 0; i < 4; i++) x[i] = M[i][4];
    return 1;
}

// ---------- Calibration Solver ----------
Params calibrate(double V[6][3], double G[6][3]) {
    int nZ = 0, nY = 0, nX = 0;
    Params p = {0};
    double kz_vals[MAX_COMB], bz_vals[MAX_COMB];
    double ky_vals[MAX_COMB], by_vals[MAX_COMB], azx_vals[MAX_COMB];
    double kx_vals[MAX_COMB], bx_vals[MAX_COMB], ayz_vals[MAX_COMB], azy_vals[MAX_COMB];printf("\nValid combinations: Z=%d, Y=%d, X=%d\n", nZ, nY, nX);
    

    int pairs[15][2]; comb2(pairs);
    for (int m = 0; m < 15; m++) {
        int i = pairs[m][0], j = pairs[m][1];
        double gz1 = G[i][2], gz2 = G[j][2];
        double vz1 = V[i][2], vz2 = V[j][2];
        if (fabs(gz1 - gz2) < 1e-12) continue;
        double kz = (vz1 - vz2) / (gz1 - gz2);
        double bz = vz1 - kz * gz1;
        kz_vals[nZ] = kz; bz_vals[nZ] = bz; nZ++;
    }

    int triplets[20][3]; comb3(triplets);
    for (int m = 0; m < 20; m++) {
        int i = triplets[m][0], j = triplets[m][1], k = triplets[m][2];
        double A[3][3] = {
            {1, G[i][1], G[i][0]},
            {1, G[j][1], G[j][0]},
            {1, G[k][1], G[k][0]}
        };
        double b[3] = {V[i][1], V[j][1], V[k][1]};
        double x[3];
        if (!solve_3x3(A, b, x)) continue;
        double by = x[0], ky = x[1], ky_azx = x[2];
        if (fabs(ky) < 1e-12) continue;
        ky_vals[nY] = ky; by_vals[nY] = by; azx_vals[nY] = ky_azx / ky; nY++;
    }

    int quads[15][4]; comb4(quads);
    for (int m = 0; m < 15; m++) {
        int i = quads[m][0], j = quads[m][1], k = quads[m][2], l = quads[m][3];
        double A[4][4] = {
            {1, G[i][0], G[i][2], G[i][1]},
            {1, G[j][0], G[j][2], G[j][1]},
            {1, G[k][0], G[k][2], G[k][1]},
            {1, G[l][0], G[l][2], G[l][1]}
        };
        double b[4] = {V[i][0], V[j][0], V[k][0], V[l][0]};
        double x[4];
        if (!solve_4x4(A, b, x)) continue;
        double bx = x[0], kx = x[1];
        if (fabs(kx) < 1e-12) continue;
        kx_vals[nX] = kx; bx_vals[nX] = bx; ayz_vals[nX] = x[2] / kx; azy_vals[nX] = x[3] / kx; nX++;
    }

    for (int i = 0; i < nZ; i++) { p.kz += kz_vals[i]; p.bz += bz_vals[i]; }
    if (nZ > 0) { p.kz /= nZ; p.bz /= nZ; }

    for (int i = 0; i < nY; i++) { p.ky += ky_vals[i]; p.by += by_vals[i]; p.alpha_zx += azx_vals[i]; }
    if (nY > 0) { p.ky /= nY; p.by /= nY; p.alpha_zx /= nY; }

    for (int i = 0; i < nX; i++) {
        p.kx += kx_vals[i]; p.bx += bx_vals[i];
        p.alpha_yz += ayz_vals[i]; p.alpha_zy += azy_vals[i];
    }
    if (nX > 0) {
        p.kx /= nX; p.bx /= nX;
        p.alpha_yz /= nX; p.alpha_zy /= nX;
    }

    printf("\nValid combinations: Z=%d, Y=%d, X=%d\n", nZ, nY, nX);
    return p;
}

// ---------- Top-level runner ----------
void run_calibration(const char *basepath) {
    MeanReading means[N_ORIENT];
    char fname[600];
    double V[6][3];

    double G[6][3] = {
        {0, -G_CONST, 0},
        { G_CONST, 0, 0},
        {0,  G_CONST, 0},
        {-G_CONST, 0, 0},
        {0, 0,  G_CONST},
        {0, 0, -G_CONST}
    };

    printf("=== Combination-based Accelerometer Calibration (C version) ===\n");

    // strcpy(basepath, "/home/kjw2kor/shared_folder/CalibrationMetrices/");

    for (int i = 0; i < N_ORIENT; i++) {
        snprintf(fname, sizeof(fname), "%s/idle%d_M0.csv", basepath, i + 1);
        if (!read_csv_means(fname, &means[i])) {
            fprintf(stderr, "Failed to read %s\n", fname);
            return;
        }
        V[i][0] = means[i].meanAx;
        V[i][1] = means[i].meanAy;
        V[i][2] = means[i].meanAz;
        printf("Orientation %d (%s): Mean Ax,Ay,Az = %.6f, %.6f, %.6f\n",
               i + 1, fname, V[i][0], V[i][1], V[i][2]);
    }

    Params p = calibrate(V, G);
    
    printf("\n=== Final Calibration Results ===\n");
    printf("kx = %.8f, ky = %.8f, kz = %.8f\n", p.kx, p.ky, p.kz);
    printf("bx = %.8f, by = %.8f, bz = %.8f\n", p.bx, p.by, p.bz);
    printf("alpha_yz = %.8f, alpha_zy = %.8f, alpha_zx = %.8f\n",
           p.alpha_yz, p.alpha_zy, p.alpha_zx);

    Sa[0][0] = p.kx; Sa[0][1] = 0; Sa[0][2] = 0;
    Sa[1][0] = 0; Sa[1][1] = p.ky; Sa[1][2] = 0;
    Sa[2][0] = 0; Sa[2][1] = 0; Sa[2][2] = p.kz;

    Ma[0][0] = 1; Ma[0][1] = 0; Ma[0][2] = 0;
    Ma[1][0] = 0; Ma[1][1] = 1; Ma[1][2] = p.alpha_yz;
    Ma[2][0] = p.alpha_zx; Ma[2][1] = p.alpha_zy; Ma[2][2] = 1;

    ba[0] = p.bx; ba[1] = p.by; ba[2] = p.bz;

    printf("\nK_a (scale factors):\n");
    for (int i = 0; i < 3; i++)
        printf("%12.8f %12.8f %12.8f\n", Sa[i][0], Sa[i][1], Sa[i][2]);

    printf("\nM_a (misalignment):\n");
    for (int i = 0; i < 3; i++)
        printf("%12.8f %12.8f %12.8f\n", Ma[i][0], Ma[i][1], Ma[i][2]);

    printf("\nb_a (bias vector):\n");
    printf("[%12.8f, %12.8f, %12.8f]\n", ba[0], ba[1], ba[2]);
}

