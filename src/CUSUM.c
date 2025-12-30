#include "CUSUM.h"
// #include "Polling.h" 
#include <string.h>
  // contains SpikeEvent definition

// Add a new sample to the circular reference window
void add_to_window(Window *w, double value) {
    w->data[w->idx] = value;
    w->idx = (w->idx + 1) % REF_SIZE;
    if (w->idx == 0) w->full = 1;
}

// Compute mean of full or partial window
double window_mean(Window *w) {
    int n = w->full ? REF_SIZE : w->idx;
    if (n == 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += w->data[i];
    return sum / n;
}

// Compute variance
double window_var(Window *w, double mu) {
    int n = w->full ? REF_SIZE : w->idx;
    if (n < 2) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n; i++)
        s += (w->data[i] - mu) * (w->data[i] - mu);
    return s / (n - 1);
}

// Compute mean excluding last 'exclude' samples
double window_mean_exclude_recent(Window *w, int exclude) {
    int n = w->full ? REF_SIZE : w->idx;
    if (n == 0) return 0.0;
    int limit = n - exclude;
    if (limit < 1) limit = n;
    double sum = 0.0;
    for (int i = 0; i < limit; i++) sum += w->data[i];
    return sum / limit;
}

// Main CUSUM logic with adaptive alpha and centered smoothing
void run_cusum_on_csv(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("Error opening file"); return; }

    // ------------------- First pass: read all data -------------------
    const int MAX_N = 2000000;
    char **timebuf = malloc(sizeof(char*) * MAX_N);
    double *magraw = malloc(sizeof(double) * MAX_N);
    int n = 0;
    char line[LINE_LEN];

    fgets(line, LINE_LEN, fp); // skip header
    while (fgets(line, LINE_LEN, fp)) {
        char t[32];
        double Ax, Ay, Az, Gx, Gy, Gz;
        if (sscanf(line, "%31[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                   t, &Ax, &Ay, &Az, &Gx, &Gy, &Gz) != 7) continue;
        timebuf[n] = malloc(32);
        strncpy(timebuf[n], t, 31);
        timebuf[n][31] = 0;
        magraw[n] = Ax;  // currently using Ax only
        n++;
        if (n >= MAX_N) break;
    }
    fclose(fp);
    if (n == 0) return;

    // ------------------- Centered smoothing -------------------
    int W = SMOOTH_WIN_CENTERED;
    int half = W / 2;
    double *mags = malloc(sizeof(double) * n);
    for (int i = 0; i < n; ++i) {
        int lo = i - half; if (lo < 0) lo = 0;
        int hi = i + half; if (hi >= n) hi = n - 1;
        double s = 0.0; int cnt = 0;
        for (int j = lo; j <= hi; ++j) { s += magraw[j]; cnt++; }
        mags[i] = s / (double)cnt;
    }

    // ------------------- Initialize reference window -------------------
    Window ref = { .idx = 0, .full = 0 };
    for (int i = 0; i < REF_SIZE && i < n; ++i)
        add_to_window(&ref, mags[i]);

    double S_pos = 0.0, S_neg = 0.0;
    int cooldown = 0;

    // ------------------- Main CUSUM loop -------------------
    for (int i = REF_SIZE; i < n; ++i) {
        add_to_window(&ref, mags[i]);
        double mu0 = window_mean_exclude_recent(&ref, EXCLUDE_RECENT);
        double sigma2 = window_var(&ref, mu0);
        if (sigma2 < CUSUM_EPS) sigma2 = CUSUM_EPS;

        double delta_in = DELTA_MAG;
        double mu1_in = mu0 + delta_in;
        double s_in = (delta_in / sigma2) * (mags[i] - (mu1_in + mu0) / 2.0);

        double delta_de = -DELTA_MAG;
        double mu1_de = mu0 + delta_de;
        double s_de = (delta_de / sigma2) * (mags[i] - (mu1_de + mu0) / 2.0);

        // Adaptive threshold
        double expected_per_sample = (delta_in * delta_in) / (2.0 * sigma2);
        double expected_S = expected_per_sample * (double)EVENT_LEN;
        double alpha = ALPHA_BASE + ALPHA_FACTOR * expected_S;

        // Update CUSUM stats
        S_pos = fmax(0.0, S_pos + s_in);
        S_neg = fmax(0.0, S_neg + s_de);
        if (cooldown > 0) cooldown--;

        // -------- Positive spike detection --------
        if (S_pos > alpha && cooldown == 0) {
            int lo = i - PEAK_PRE; if (lo < 0) lo = 0;
            int hi = i + PEAK_POST; if (hi >= n) hi = n - 1;
            double best = -1e9; int bi = i;
            for (int k = lo; k <= hi; ++k)
                if (mags[k] > best) { best = mags[k]; bi = k; }
            printf("Positive spike at %s (%.3f m/s²)\n", timebuf[bi], mags[bi]);
            // char timeStr[32];
            // sprintf(timeStr, "%02d:%02d:%02d:%03d", hour, min, sec, ms);
            // addCUSUMDetection("Positive", timeStr, amplitude);
            addCUSUMDetection("Positive", timebuf[bi], mags[bi]);
            S_pos = 0.0; S_neg = 0.0; cooldown = RESET_DELAY;
        }

        // -------- Negative spike detection --------
        if (S_neg > alpha && cooldown == 0) {
            int lo = i - PEAK_PRE; if (lo < 0) lo = 0;
            int hi = i + PEAK_POST; if (hi >= n) hi = n - 1;
            double best = 1e9; int bi = i;
            for (int k = lo; k <= hi; ++k)
                if (mags[k] < best) { best = mags[k]; bi = k; }
            printf("Negative spike at %s (%.3f m/s²)\n", timebuf[bi], mags[bi]);
            // char timeStr[32];
            // sprintf(timeStr, "%02d:%02d:%02d:%03d", hour, min, sec, ms);
            // addCUSUMDetection("Negetive", timeStr, amplitude);
            addCUSUMDetection("Negative", timebuf[bi], mags[bi]);
            S_pos = 0.0; S_neg = 0.0; cooldown = RESET_DELAY;
        }
    }

    // ------------------- Cleanup -------------------
    for (int i = 0; i < n; ++i) free(timebuf[i]);
    free(timebuf); free(magraw); free(mags);

    printf("\nAdaptive CUSUM thresholding applied.\n");
    printf("----------------------------------------------------------\n");
}

static Detection cusumDet[MAX_DETECTIONS];
static int cusumCount = 0;

void addCUSUMDetection(const char *type, const char *timeStr, float amp) {
    if (cusumCount >= MAX_DETECTIONS) return;
    Detection *d = &cusumDet[cusumCount++];
    strcpy(d->polarity, type);
    int h,m,s,ms;
    sscanf(timeStr, "%d:%d:%d:%d", &h,&m,&s,&ms);
    d->hour=h; d->min=m; d->sec=s; d->ms=ms; d->amplitude=amp;
}

// void runCUSUM(const char *filename) {
//     // your detection logic here → call addCUSUMDetection()
// }

int getCUSUMDetections(Detection **array) {
    *array = cusumDet;
    return cusumCount;
}