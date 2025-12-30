#ifndef NORMALIZATION_H
#define NORMALIZATION_H

// Normalize accelerometer columns (Ax, Ay, Az) by dividing by 9.8
// in the given CSV file. Returns 0 on success, 1 on failure.
int normalize_csv(const char *filename);

#endif
