#include "Orientation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void transform_row(double *Ax, double *Ay, double *Az, double *Gx, double *Gy, double *Gz) {
    double ax = *Ax, ay = *Ay, az = *Az;
    double gx = *Gx, gy = *Gy, gz = *Gz;

    double absAx = fabs(ax), absAy = fabs(ay), absAz = fabs(az);
    int orientation = 0;

    if (absAx >= absAy && absAx >= absAz) {
        orientation = (ax > 0) ? 2 : 4;
    } else if (absAy >= absAx && absAy >= absAz) {
        orientation = (ay > 0) ? 3 : 1;
    } else {
        orientation = (az < 0) ? 6 : 5;
    }

    switch (orientation) {
        case 1:
            *Ax = -az; *Ay = ax; *Az = -ay;
            *Gx = -gz; *Gy = gx; *Gz = -gy;
            break;
        case 2:
            *Ax = -az; *Ay = ay; *Az = ax;
            *Gx = -gz; *Gy = gy; *Gz = gx;
            break;
        case 3:
            *Ax = -az; *Ay = -ax; *Az = ay;
            *Gx = -gz; *Gy = -gx; *Gz = gy;
            break;
        case 4:
            *Ax = -az; *Ay = -ay; *Az = -ax;
            *Gx = -gz; *Gy = -gy; *Gz = -gx;
            break;
        case 5:
            *Ax = ax; *Ay = ay; *Az = az;
            *Gx = gx; *Gy = gy; *Gz = gz;
            break;
        case 6:
            *Ax = -ax; *Ay = ay; *Az = -az;
            *Gx = -gx; *Gy = gy; *Gz = -gz;
            break;
        default:
            break;
    }
}

void process_csv_file(const char *filename) {
    char tempfile[512];
    snprintf(tempfile, sizeof(tempfile), "%s.temp", filename);

    FILE *fp = fopen(filename, "r");
    FILE *out = fopen(tempfile, "w");
    if (!fp || !out) {
        perror("File open failed");
        if (fp) fclose(fp);
        if (out) fclose(out);
        return;
    }

    char line[LINE_SIZE];
    int firstLine = 1;

    while (fgets(line, sizeof(line), fp)) {
        if (firstLine) {
            fputs(line, out);
            firstLine = 0;
            continue;
        }

        size_t L = strlen(line);
        if (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = '\0';
        if (L && (line[L - 1] == '\r')) line[--L] = '\0';

        char buf[LINE_SIZE];
        strncpy(buf, line, LINE_SIZE);
        buf[LINE_SIZE - 1] = '\0';

        char *saveptr = NULL;
        char *tokens[64];
        int n = 0;
        char *tok = strtok_r(buf, ",", &saveptr);
        while (tok && n < (int)(sizeof(tokens) / sizeof(tokens[0]))) {
            tokens[n++] = tok;
            tok = strtok_r(NULL, ",", &saveptr);
        }

        if (n >= 7) {
            char *time_str = tokens[0];
            double Ax = strtod(tokens[1], NULL);
            double Ay = strtod(tokens[2], NULL);
            double Az = strtod(tokens[3], NULL);
            double Gx = strtod(tokens[4], NULL);
            double Gy = strtod(tokens[5], NULL);
            double Gz = strtod(tokens[6], NULL);

            transform_row(&Ax, &Ay, &Az, &Gx, &Gy, &Gz);

            fprintf(out, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f", time_str, Ax, Ay, Az, Gx, Gy, Gz);
            for (int i = 7; i < n; ++i)
                fprintf(out, ",%s", tokens[i]);
            fprintf(out, "\n");
        } else {
            fputs(line, out);
            fputc('\n', out);
        }
    }

    fclose(fp);
    fclose(out);

    if (remove(filename) != 0)
        perror("remove original");
    if (rename(tempfile, filename) != 0)
        perror("rename temp");

    printf("CSV file updated in-place: %s\n", filename);
}
