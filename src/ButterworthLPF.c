#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ButterworthLPF.h"

// =======================================================
// 4th Order Butterworth LPF (fc = 5 Hz, fs = 50 Hz)
// =======================================================
void butterworth_init(IIRFilter *f) {
    // Coefficients computed via bilinear transform (fc=5Hz, fs=50Hz)
    f->b[0] = 0.00041655;
    f->b[1] = 0.00166619;
    f->b[2] = 0.00249928;
    f->b[3] = 0.00166619;
    f->b[4] = 0.00041655;

    f->a[0] = 1.00000000;
    f->a[1] = -3.18063855;
    f->a[2] = 3.86119435;
    f->a[3] = -2.11215536;
    f->a[4] = 0.43826514;

    memset(f->x, 0, sizeof(f->x));
    memset(f->y, 0, sizeof(f->y));
}

double butterworth_process(IIRFilter *f, double input) {
    for (int i = FILTER_ORDER; i > 0; i--) {
        f->x[i] = f->x[i - 1];
        f->y[i] = f->y[i - 1];
    }

    f->x[0] = input;
    f->y[0] = f->b[0]*f->x[0] + f->b[1]*f->x[1] + f->b[2]*f->x[2] +
              f->b[3]*f->x[3] + f->b[4]*f->x[4]
              - f->a[1]*f->y[1] - f->a[2]*f->y[2]
              - f->a[3]*f->y[3] - f->a[4]*f->y[4];

    return f->y[0];
}

// =======================================================
// CSV I/O helpers
// =======================================================
int read_multicol_csv(const char *filename, char time[][32],
                      double data[][MAX_COLS], int max_rows, int *num_cols) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening CSV");
        return -1;
    }

    char line[MAX_LINE_LEN];
    int row = 0;

    fgets(line, sizeof(line), fp); // skip header

    while (fgets(line, sizeof(line), fp) && row < max_rows) {
        char *token = strtok(line, ",");
        if (!token) continue;
        strncpy(time[row], token, 31);
        time[row][31] = '\0';

        int col = 0;
        while ((token = strtok(NULL, ",")) && col < MAX_COLS) {
            data[row][col] = atof(token);
            col++;
        }
        *num_cols = col;
        row++;
    }

    fclose(fp);
    return row;
}

int write_multicol_csv(const char *filename, char time[][32],
                       double data[][MAX_COLS], int rows, int cols) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error writing CSV");
        return -1;
    }

    fprintf(fp, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");
    for (int i = 0; i < rows; i++) {
        fprintf(fp, "%s", time[i]);
        for (int j = 0; j < cols; j++)
            fprintf(fp, ",%.6f", data[i][j]);
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}
