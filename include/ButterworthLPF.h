#ifndef BUTTERWORTH_LPF_H
#define BUTTERWORTH_LPF_H

#define MAX_COLS 8
#define MAX_ROWS 10000
#define MAX_LINE_LEN 256
#define FILTER_ORDER 4

typedef struct {
    double a[FILTER_ORDER + 1];
    double b[FILTER_ORDER + 1];
    double x[FILTER_ORDER + 1];
    double y[FILTER_ORDER + 1];
} IIRFilter;

void butterworth_init(IIRFilter *f);
double butterworth_process(IIRFilter *f, double input);
int read_multicol_csv(const char *filename, char time[][32],
                      double data[][MAX_COLS], int max_rows, int *num_cols);
int write_multicol_csv(const char *filename, char time[][32],
                       double data[][MAX_COLS], int rows, int cols);

#endif


