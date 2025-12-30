#include "normalization.h"
#include <stdio.h>
#include <stdlib.h>

#define LINE_SIZE 512

int normalize_csv(const char *filename) {
    const char *tempfile = "temp_normalized.csv";

    FILE *infile = fopen(filename, "r");
    FILE *outfile = fopen(tempfile, "w");

    if (!infile || !outfile) {
        perror("File open error");
        if (infile) fclose(infile);
        if (outfile) fclose(outfile);
        return 1;
    }

    char line[LINE_SIZE];
    int is_header = 1;

    while (fgets(line, sizeof(line), infile)) {
        if (is_header) {
            // Write header line directly
            fputs(line, outfile);
            is_header = 0;
            continue;
        }

        double Ax, Ay, Az, Gx, Gy, Gz;
        char Time[64];

        // Expect: Time,Ax,Ay,Az,Gx,Gy,Gz
        if (sscanf(line, "%63[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                   Time, &Ax, &Ay, &Az, &Gx, &Gy, &Gz) == 7) {
            // Normalize acceleration values
            Ax /= 9.8;
            Ay /= 9.8;
            Az /= 9.8;

            fprintf(outfile, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    Time, Ax, Ay, Az, Gx, Gy, Gz);
        }
    }

    fclose(infile);
    fclose(outfile);

    // Replace original file
    remove(filename);
    rename(tempfile, filename);

    printf("✅ Normalized accelerometer columns (Ax, Ay, Az) in %s successfully.\n", filename);
    return 0;
}
