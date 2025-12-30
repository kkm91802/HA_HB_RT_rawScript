#include "Polling.h"
#include "DerivativeDetect.h"
#include "CUSUM.h"
#include "CWT.h"
#include "EnergyEnvelope.h"
#include "Threshold.h"

// int run_cusum_detection(const char *filename, SpikeEvent *events, size_t max_events) {
//     FILE *fp = fopen(filename, "r");
//     if (!fp) { perror("Error opening file"); return 0; }

//     const int MAX_N = 2000000;
//     char **timebuf = malloc(sizeof(char*) * MAX_N);
//     double *magraw = malloc(sizeof(double) * MAX_N);
//     double *t_seconds = malloc(sizeof(double) * MAX_N);  // Declare t_seconds array
//     int n = 0;
//     char line[LINE_LEN];

//     fgets(line, LINE_LEN, fp);  // Skip header
//     while (fgets(line, LINE_LEN, fp)) {
//         char t[32];
//         double Ax, Ay, Az, Gx, Gy, Gz;
//         if (sscanf(line, "%31[^,],%lf,%lf,%lf,%lf,%lf,%lf",
//                    t, &Ax, &Ay, &Az, &Gx, &Gy, &Gz) != 7) continue;
//         timebuf[n] = malloc(32);
//         strncpy(timebuf[n], t, 31);
//         timebuf[n][31] = 0;
//         magraw[n] = Ax;

//         // Convert the timestamp string to seconds and store it in t_seconds[]
//         int hh, mm, ss, ms;
//         if (sscanf(t, "%d:%d:%d:%d", &hh, &mm, &ss, &ms) == 4) {
//             t_seconds[n] = hh * 3600.0 + mm * 60.0 + ss + ms / 1000.0;
//         }

//         n++;
//         if (n >= MAX_N) break;
//     }
//     fclose(fp);
//     if (n == 0) return 0;

//     // Smoothing step (same as before)
//     int W = SMOOTH_WIN_CENTERED;
//     int half = W / 2;
//     double *mags = malloc(sizeof(double) * n);
//     for (int i = 0; i < n; ++i) {
//         int lo = i - half; if (lo < 0) lo = 0;
//         int hi = i + half; if (hi >= n) hi = n - 1;
//         double s = 0.0; int cnt = 0;
//         for (int j = lo; j <= hi; ++j) { s += magraw[j]; cnt++; }
//         mags[i] = s / (double)cnt;
//     }

//     Window ref = { .idx = 0, .full = 0 };
//     for (int i = 0; i < REF_SIZE && i < n; ++i)
//         add_to_window(&ref, mags[i]);

//     double S_pos = 0.0, S_neg = 0.0;
//     int cooldown = 0;
//     int event_count = 0;  // Counter for detected events

//     for (int i = REF_SIZE; i < n; ++i) {
//         add_to_window(&ref, mags[i]);
//         double mu0 = window_mean_exclude_recent(&ref, EXCLUDE_RECENT);
//         double sigma2 = window_var(&ref, mu0);
//         if (sigma2 < CUSUM_EPS) sigma2 = CUSUM_EPS;

//         double delta_in = DELTA_MAG;
//         double mu1_in = mu0 + delta_in;
//         double s_in = (delta_in / sigma2) * (mags[i] - (mu1_in + mu0) / 2.0);

//         double delta_de = -DELTA_MAG;
//         double mu1_de = mu0 + delta_de;
//         double s_de = (delta_de / sigma2) * (mags[i] - (mu1_de + mu0) / 2.0);

//         double expected_per_sample = (delta_in * delta_in) / (2.0 * sigma2);
//         double expected_S = expected_per_sample * (double)EVENT_LEN;
//         double alpha = ALPHA_BASE + ALPHA_FACTOR * expected_S;

//         S_pos = fmax(0.0, S_pos + s_in);
//         S_neg = fmax(0.0, S_neg + s_de);
//         if (cooldown > 0) cooldown--;

//         // ---- Positive spike ----
//         if (S_pos > alpha && cooldown == 0) {
//             int lo = i - PEAK_PRE; if (lo < 0) lo = 0;
//             int hi = i + PEAK_POST; if (hi >= n) hi = n - 1;
//             double best = -1e9; int bi = i;
//             for (int k = lo; k <= hi; ++k)
//                 if (mags[k] > best) { best = mags[k]; bi = k; }

//             // Round the timestamp and store the spike event
//             int hh, mm, ss, ms;
//             if (sscanf(timebuf[bi], "%d:%d:%d:%d", &hh, &mm, &ss, &ms) == 4) {
//                 // Round milliseconds to nearest second
//                 if (ms >= 500) {
//                     ss += 1;  // Round up if ms >= 500
//                 }

//                 // Ensure seconds are within bounds (0-59)
//                 if (ss == 60) {
//                     ss = 0;
//                     mm += 1;  // Increment minute if seconds roll over
//                 }

//                 // Ensure minutes are within bounds (0-59)
//                 if (mm == 60) {
//                     mm = 0;
//                     hh += 1;  // Increment hour if minutes roll over
//                 }

//                 // Format rounded time as HH:MM:SS
//                 snprintf(events[event_count].time_str, sizeof(events[event_count].time_str),
//                          "%02d:%02d:%02d", hh, mm, ss);
//             }

//             // Store spike details
//             events[event_count].time_ms = t_seconds[bi] * 1000.0;  // Convert to milliseconds
//             events[event_count].magnitude = mags[bi];
//             events[event_count].type = +1;  // Positive spike
//             event_count++;

//             S_pos = 0.0; S_neg = 0.0; cooldown = RESET_DELAY;
//         }

//         // ---- Negative spike ----
//         if (S_neg > alpha && cooldown == 0) {
//             int lo = i - PEAK_PRE; if (lo < 0) lo = 0;
//             int hi = i + PEAK_POST; if (hi >= n) hi = n - 1;
//             double best = 1e9; int bi = i;
//             for (int k = lo; k <= hi; ++k)
//                 if (mags[k] < best) { best = mags[k]; bi = k; }

//             // Round the timestamp and store the negative spike event
//             int hh, mm, ss, ms;
//             if (sscanf(timebuf[bi], "%d:%d:%d:%d", &hh, &mm, &ss, &ms) == 4) {
//                 // Round milliseconds to nearest second
//                 if (ms >= 500) {
//                     ss += 1;  // Round up if ms >= 500
//                 }

//                 // Ensure seconds are within bounds (0-59)
//                 if (ss == 60) {
//                     ss = 0;
//                     mm += 1;  // Increment minute if seconds roll over
//                 }

//                 // Ensure minutes are within bounds (0-59)
//                 if (mm == 60) {
//                     mm = 0;
//                     hh += 1;  // Increment hour if minutes roll over
//                 }

//                 // Format rounded time as HH:MM:SS
//                 snprintf(events[event_count].time_str, sizeof(events[event_count].time_str),
//                          "%02d:%02d:%02d", hh, mm, ss);
//             }

//             // Store spike details
//             events[event_count].time_ms = t_seconds[bi] * 1000.0;
//             events[event_count].magnitude = mags[bi];
//             events[event_count].type = -1;  // Negative spike
//             event_count++;

//             S_pos = 0.0; S_neg = 0.0; cooldown = RESET_DELAY;
//         }
//     }

//     for (int i = 0; i < n; ++i) free(timebuf[i]);
//     free(timebuf); free(magraw); free(mags);

//     return event_count;  // Return the number of detected events
// }


// int run_energy_detection(const char *filename, SpikeEvent *events, size_t max_events) {
//     // Load the CSV data
//     EnvelopeSample data[MAX_SAMPLES_ENVELOPE];
//     int totalSamples = loadEnvelopeCSV(filename, data);
//     if (totalSamples <= 0) {
//         fprintf(stderr, "No valid data loaded.\n");
//         return 0;
//     }

//     // Declare t_seconds array to store timestamp in seconds
//     double t_seconds[MAX_SAMPLES_ENVELOPE];

//     // Sampling rate (fs) - replace with your actual sampling rate if needed
//     double fs = 1000.0;  // Example: 1000 Hz, replace with actual sampling rate if available

//     // Populate t_seconds[] with time in seconds
//     for (int i = 0; i < totalSamples; i++) {
//         int hh, mm, ss, ms;
//         if (sscanf(data[i].time, "%d:%d:%d:%d", &hh, &mm, &ss, &ms) == 4) {
//             t_seconds[i] = hh * 3600.0 + mm * 60.0 + ss + ms / 1000.0;  // Convert to seconds
//         } else {
//             t_seconds[i] = 0.0;  // If the time parsing fails, set it to 0
//         }
//     }

//     // Event counter for detected spikes
//     int event_count = 0;

//     // Call the sliding window detection method
//     runImprovedEnvelopeSliding(data, totalSamples, fs);

//     // Loop over detected events and fill the SpikeEvent array
//     // Example: Store detected positive spikes in events[]
//     for (int i = 0; i < totalSamples && event_count < max_events; i++) {
//         if (fabs(data[i].Ax) >= 0.5) {  // Adjust this condition based on your detection logic
//             strncpy(events[event_count].time_str, data[i].time, sizeof(events[event_count].time_str) - 1);
//             events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//             events[event_count].time_ms = t_seconds[i] * 1000.0;  // Convert to milliseconds
//             events[event_count].magnitude = data[i].Ax;
//             events[event_count].type = +1;  // Positive spike
//             event_count++;
//         }
//     }

//     return event_count;  // Return the number of detected events
// }

// int run_cwt_detection(const char *filename, SpikeEvent *events, size_t max_events) {
//     FILE *fp = fopen(filename, "r");
//     if (!fp) {
//         perror("fopen");
//         return 1;
//     }

//     Sample *samples = malloc(MAX_LINES * sizeof(Sample));
//     int N = 0;
//     char line[512];
//     fgets(line, sizeof(line), fp); // header
//     while (fgets(line, sizeof(line), fp)) {
//         if (N >= MAX_LINES) break;
//         char *tok = strtok(line, ",");
//         if (!tok) continue;
//         strncpy(samples[N].time_str, tok, MAX_TIMESTR_LEN);
//         samples[N].t_seconds = parse_time_to_seconds(tok);

//         tok = strtok(NULL, ",");
//         if (!tok) continue;
//         samples[N].ax = atof(tok);
//         tok = strtok(NULL, ",");
//         if (!tok) continue;
//         samples[N].ay = atof(tok);
//         tok = strtok(NULL, ",");
//         if (!tok) continue;
//         samples[N].az = atof(tok);

//         samples[N].amp = samples[N].ax;  // Only using ax here
//         N++;
//     }
//     fclose(fp);
//     if (N < 10) {
//         fprintf(stderr, "Not enough data\n");
//         return 1;
//     }

//     double *tsec = malloc(N * sizeof(double));
//     double *amp = malloc(N * sizeof(double));
//     for (int i = 0; i < N; i++) {
//         tsec[i] = samples[i].t_seconds;
//         amp[i] = samples[i].amp;
//     }

//     // Estimate sampling rate
//     double *dts = malloc((N - 1) * sizeof(double));
//     for (int i = 0; i < N - 1; i++) dts[i] = tsec[i + 1] - tsec[i];
//     double med_dt = median_double(dts, N - 1);
//     double fs = (med_dt > 0) ? 1.0 / med_dt : 100.0;

//     // Define scales
//     double min_s = MIN_DUR_SEC * fs, max_s = MAX_DUR_SEC * fs;
//     double scales[N_SCALES];
//     for (int s = 0; s < N_SCALES; s++) scales[s] = min_s + (max_s - min_s) * s / (N_SCALES - 1);

//     // Compute coefficients
//     double **coeffs = malloc(N_SCALES * sizeof(double *));
//     for (int s = 0; s < N_SCALES; s++) {
//         coeffs[s] = calloc(N, sizeof(double));
//         int L = (int)ceil(10.0 * scales[s]);
//         if (L < 21) L = 21;
//         if (L > MAX_WAVELET_LEN) L = MAX_WAVELET_LEN;
//         if (L % 2 == 0) L++;
//         double *w = malloc(L * sizeof(double));
//         make_morlet(w, L, scales[s]);
//         convolve(amp, N, w, L, coeffs[s]);
//         free(w);
//     }

//     // Thresholds per scale
//     double sigma[N_SCALES], mu[N_SCALES];
//     for (int s = 0; s < N_SCALES; s++) {
//         double *absvals = malloc(N * sizeof(double));
//         for (int i = 0; i < N; i++) absvals[i] = fabs(coeffs[s][i]);
//         mu[s] = median_double(absvals, N);
//         sigma[s] = compute_mad(coeffs[s], N);
//         free(absvals);
//     }

//     unsigned char *pos_mask = calloc(N, sizeof(unsigned char));
//     unsigned char *neg_mask = calloc(N, sizeof(unsigned char));
//     for (int i = 0; i < N; i++) {
//         for (int s = 0; s < N_SCALES; s++) {
//             double thr = 1.2 * mu[s] + THRESH_K * sigma[s];
//             if (coeffs[s][i] > thr) pos_mask[i] = 1;
//             if (coeffs[s][i] < -thr) neg_mask[i] = 1;
//         }
//     }

//     double *avg = malloc(N * sizeof(double));
//     for (int i = 0; i < N; i++) {
//         avg[i] = 0;
//         for (int s = 0; s < N_SCALES; s++) avg[i] += coeffs[s][i];
//         avg[i] /= N_SCALES;
//     }

//     int *pos_idx = malloc(N * sizeof(int));
//     int *neg_idx = malloc(N * sizeof(int));
//     int pos_n = detect_regions(avg, pos_mask, N, pos_idx);
//     int neg_n = detect_regions(avg, neg_mask, N, neg_idx);

//     double *pos_t = malloc(pos_n * sizeof(double));
//     double *pos_m = malloc(pos_n * sizeof(double));
//     for (int i = 0; i < pos_n; i++) {
//         pos_t[i] = tsec[pos_idx[i]] * 1000.0;
//         pos_m[i] = amp[pos_idx[i]];
//     }

//     double *neg_t = malloc(neg_n * sizeof(double));
//     double *neg_m = malloc(neg_n * sizeof(double));
//     for (int i = 0; i < neg_n; i++) {
//         neg_t[i] = tsec[neg_idx[i]] * 1000.0;
//         neg_m[i] = amp[neg_idx[i]];
//     }

//     double *ptm = malloc(pos_n * sizeof(double)), *pmm = malloc(pos_n * sizeof(double));
//     double *ntm = malloc(neg_n * sizeof(double)), *nmm = malloc(neg_n * sizeof(double));
//     int pos_final = merge_events(pos_t, pos_m, pos_n, MERGE_WINDOW_MS, ptm, pmm);
//     int neg_final = merge_events(neg_t, neg_m, neg_n, MERGE_WINDOW_MS, ntm, nmm);

//     // Store detected spikes in SpikeEvent[] array
//     int event_count = 0;
//     for (int i = 0; i < pos_final && event_count < max_events; i++) {
//         double target = ptm[i] / 1000.0;
//         int idx = 0;
//         for (int j = 1; j < N; j++) {
//             if (fabs(tsec[j] - target) < fabs(tsec[idx] - target)) idx = j;
//         }
//         if (pmm[i] < 0.3) continue;  // Skip weak spikes
//         strncpy(events[event_count].time_str, samples[idx].time_str, sizeof(events[event_count].time_str) - 1);
//         events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//         events[event_count].time_ms = ptm[i];
//         events[event_count].magnitude = pmm[i];
//         events[event_count].type = 1;  // Positive spike
//         event_count++;
//     }

//     for (int i = 0; i < neg_final && event_count < max_events; i++) {
//         double target = ntm[i] / 1000.0;
//         int idx = 0;
//         for (int j = 1; j < N; j++) {
//             if (fabs(tsec[j] - target) < fabs(tsec[idx] - target)) idx = j;
//         }
//         if (fabs(nmm[i]) < 0.3) continue;  // Skip weak spikes
//         strncpy(events[event_count].time_str, samples[idx].time_str, sizeof(events[event_count].time_str) - 1);
//         events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//         events[event_count].time_ms = ntm[i];
//         events[event_count].magnitude = nmm[i];
//         events[event_count].type = -1;  // Negative spike
//         event_count++;
//     }

//     // Cleanup
//     free(samples);
//     free(tsec);
//     free(amp);
//     free(dts);
//     for (int s = 0; s < N_SCALES; s++) free(coeffs[s]);
//     free(coeffs);
//     free(pos_mask);
//     free(neg_mask);
//     free(avg);
//     free(pos_idx);
//     free(neg_idx);
//     free(pos_t);
//     free(pos_m);
//     free(neg_t);
//     free(neg_m);
//     free(ptm);
//     free(pmm);
//     free(ntm);
//     free(nmm);

//     return event_count;  // Return the number of detected events
// }



// int run_derivative_detection(const char *filename, SpikeEvent *events, size_t max_events) {
//     // Load the CSV data
//     DerivSample data[MAX_SAMPLES_DERIV];
//     int totalSamples = loadCSV(filename, data);
//     if (totalSamples <= 0) {
//         fprintf(stderr, "No valid data loaded.\n");
//         return 0;
//     }

//     // Declare t_seconds array to store timestamp in seconds
//     double t_seconds[MAX_SAMPLES_DERIV];

//     // Sampling rate (fs) - replace with your actual sampling rate if needed
//     double fs = 1000.0;  // Example: 1000 Hz, replace with actual sampling rate if available

//     // Populate t_seconds[] with time in seconds
//     for (int i = 0; i < totalSamples; i++) {
//         int hh, mm, ss, ms;
//         if (sscanf(data[i].time, "%d:%d:%d:%d", &hh, &mm, &ss, &ms) == 4) {
//             t_seconds[i] = hh * 3600.0 + mm * 60.0 + ss + ms / 1000.0;  // Convert to seconds
//         } else {
//             t_seconds[i] = 0.0;  // If the time parsing fails, set it to 0
//         }
//     }

//     // Event counter for detected spikes
//     int event_count = 0;

//     // Call the spike detection method
//     detectSpikes(data, totalSamples);

//     // Loop over detected events and fill the SpikeEvent array
//     // Example: Store detected positive spikes in events[]
//     for (int i = 0; i < totalSamples && event_count < max_events; i++) {
//         if (fabs(data[i].Ax) >= 0.5) {  // Adjust this condition based on your detection logic
//             strncpy(events[event_count].time_str, data[i].time, sizeof(events[event_count].time_str) - 1);
//             events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//             events[event_count].time_ms = t_seconds[i] * 1000.0;  // Convert to milliseconds
//             events[event_count].magnitude = data[i].Ax;
//             events[event_count].type = +1;  // Positive spike
//             event_count++;
//         }
//     }

//     return event_count;  // Return the number of detected events
// }



// int run_threshold_detection(const char *filename, SpikeEvent *events, size_t max_events) {
//     // Load the IMU data
//     IMUSample data[MAX_SAMPLES_IMU];
//     int totalSamples = readIMUCSV(filename, data, MAX_SAMPLES_IMU);
//     if (totalSamples <= 0) {
//         fprintf(stderr, "No valid data loaded.\n");
//         return 0;
//     }

//     // Initialize threshold parameters
//     AdaptiveThreshold th;
//     initThreshold(&th, 0.5f, -0.5f, 10);  // Example initialization

//     // Compute the magnitude of the accelerometer signal (you can modify it to use all 3 axes or a custom calculation)
//     float accMag[MAX_SAMPLES_IMU];
//     computeMagnitude(data, accMag, totalSamples);

//     // Event counter for detected spikes
//     int event_count = 0;

//     // Loop through the signal and apply the threshold detection
//     for (int i = 1; i < totalSamples - 1 && event_count < max_events; i++) {
//         // Update the adaptive threshold with the current window
//         float window[WINDOW_SIZE];
//         for (int j = 0; j < WINDOW_SIZE; j++) {
//             if (i - j >= 0) {
//                 window[j] = accMag[i - j];
//             }
//         }
//         updateAdaptiveThreshold(&th, window, WINDOW_SIZE);

//         // Detect positive and negative spikes
//         if (detectPositiveSpike(accMag, i, &th)) {
//             // Store detected positive spike in events[] array
//             strncpy(events[event_count].time_str, data[i].timeStr, sizeof(events[event_count].time_str) - 1);
//             events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//             events[event_count].time_ms = data[i].timeSec * 1000.0;  // Convert to milliseconds
//             events[event_count].magnitude = accMag[i];
//             events[event_count].type = +1;  // Positive spike
//             event_count++;
//         }

//         if (detectNegativeSpike(accMag, i, &th)) {
//             // Store detected negative spike in events[] array
//             strncpy(events[event_count].time_str, data[i].timeStr, sizeof(events[event_count].time_str) - 1);
//             events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//             events[event_count].time_ms = data[i].timeSec * 1000.0;
//             events[event_count].magnitude = accMag[i];
//             events[event_count].type = -1;  // Negative spike
//             event_count++;
//         }
//     }

//     return event_count;  // Return the number of detected events
// }



