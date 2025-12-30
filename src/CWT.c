#include "CWT.h"
// #include "Polling.h" 

// --- Utility functions ---
static int cmp_double(const void *a, const void *b) {
    double va = *(const double*)a, vb = *(const double*)b;
    return (va > vb) - (va < vb);
}

static double median_double(double *arr, int n) {
    if (n <= 0) return 0.0;
    qsort(arr, n, sizeof(double), cmp_double);
    return (n % 2) ? arr[n/2] : 0.5*(arr[n/2 - 1] + arr[n/2]);
}

static double compute_mad(double *x, int n) {
    if (n <= 0) return 0.0;
    double *copy = malloc(n*sizeof(double));
    for (int i=0;i<n;i++) copy[i]=x[i];
    double med = median_double(copy, n);
    for (int i=0;i<n;i++) copy[i] = fabs(x[i]-med);
    double mad = median_double(copy, n);
    free(copy);
    return mad / MAD_SCALE_FACTOR;
}

static double parse_time_to_seconds(const char *tstr) {
    int hh=0, mm=0, ss=0, ms=0;
    if (sscanf(tstr, "%d:%d:%d:%d", &hh, &mm, &ss, &ms)==4)
        return hh*3600.0 + mm*60.0 + ss + ms/1000.0;
    double sec_with_ms=0.0;
    if (sscanf(tstr, "%d:%d:%lf", &hh, &mm, &sec_with_ms)==3)
        return hh*3600.0 + mm*60.0 + sec_with_ms;
    return 0.0;
}

// --- Wavelet functions ---
static void make_morlet(double *w, int L, double scale_samples) {
    int mid=L/2;
    double sigma=scale_samples;
    if(sigma<1.0)sigma=1.0;
    double k0=5.0, sumabs=0.0;
    for(int i=0;i<L;i++){
        double t=(i-mid)/sigma;
        double val=cos(2*M_PI*k0*t)*exp(-0.5*t*t);
        w[i]=val;
        sumabs+=fabs(val);
    }
    if(sumabs==0.0)return;
    for(int i=0;i<L;i++) w[i]/=sumabs;
}

static void convolve(const double *sig, int N, const double *ker, int L, double *out) {
    int half=L/2;
    for(int i=0;i<N;i++){
        double s=0.0;
        for(int k=0;k<L;k++){
            int j=i-k+half;
            if(j>=0&&j<N) s+=sig[j]*ker[k];
        }
        out[i]=s;
    }
}

// --- Spike utilities ---
static int detect_regions(double *vals, unsigned char *mask, int N, int *indices) {
    int count=0;
    for(int i=0;i<N;){
        if(!mask[i]){i++;continue;}
        int start=i, best=i;
        double bestv=fabs(vals[i]);
        while(i<N && mask[i]){
            if(fabs(vals[i])>bestv){bestv=fabs(vals[i]);best=i;}
            i++;
        }
        indices[count++]=best;
    }
    return count;
}

static int merge_events(double *times_ms,double *mags,int count,double merge_ms,double *out_t,double *out_m){
    if(count==0)return 0;
    int n=0;double cur_t=times_ms[0],cur_m=mags[0];
    for(int i=1;i<count;i++){
        if(times_ms[i]-cur_t<=merge_ms){
            if(fabs(mags[i])>fabs(cur_m)){cur_t=times_ms[i];cur_m=mags[i];}
        }else{
            out_t[n]=cur_t;out_m[n++]=cur_m;
            cur_t=times_ms[i];cur_m=mags[i];
        }
    }
    out_t[n]=cur_t;out_m[n++]=cur_m;
    return n;
}

// --- Core pipeline ---
int run_cwt_to_csv(const char *filename) {
    FILE *fp=fopen(filename,"r");
    if(!fp){perror("fopen");return 1;}

    Sample *samples=malloc(MAX_LINES*sizeof(Sample));
    int N=0;
    char line[512];
    fgets(line,sizeof(line),fp); // header
    while(fgets(line,sizeof(line),fp)){
        if(N>=MAX_LINES)break;
        char *tok=strtok(line,",");
        if(!tok)continue;
        strncpy(samples[N].time_str,tok,MAX_TIMESTR_LEN);
        samples[N].t_seconds=parse_time_to_seconds(tok);

        tok=strtok(NULL,","); if(!tok)continue; samples[N].ax=atof(tok);
        tok=strtok(NULL,","); if(!tok)continue; samples[N].ay=atof(tok);
        tok=strtok(NULL,","); if(!tok)continue; samples[N].az=atof(tok);

        // samples[N].amp=sqrt(samples[N].ax*samples[N].ax+samples[N].ay*samples[N].ay+samples[N].az*samples[N].az);

        samples[N].amp=samples[N].ax;
        N++;
    }
    fclose(fp);
    if(N<10){fprintf(stderr,"Not enough data\n");return 1;}

    double *tsec=malloc(N*sizeof(double));
    double *amp=malloc(N*sizeof(double));
    for(int i=0;i<N;i++){tsec[i]=samples[i].t_seconds;amp[i]=samples[i].amp;}

    // Estimate sampling rate
    double *dts=malloc((N-1)*sizeof(double));
    for(int i=0;i<N-1;i++) dts[i]=tsec[i+1]-tsec[i];
    double med_dt=median_double(dts,N-1);
    double fs=(med_dt>0)?1.0/med_dt:100.0;

    // Define scales
    double min_s=MIN_DUR_SEC*fs,max_s=MAX_DUR_SEC*fs;
    double scales[N_SCALES];
    for(int s=0;s<N_SCALES;s++)scales[s]=min_s+(max_s-min_s)*s/(N_SCALES-1);

    // Compute coefficients
    double **coeffs=malloc(N_SCALES*sizeof(double*));
    for(int s=0;s<N_SCALES;s++){
        coeffs[s]=calloc(N,sizeof(double));
        int L=(int)ceil(10.0*scales[s]);
        if(L<21)L=21; if(L>MAX_WAVELET_LEN)L=MAX_WAVELET_LEN;
        if(L%2==0)L++;
        double *w=malloc(L*sizeof(double));
        make_morlet(w,L,scales[s]);
        convolve(amp,N,w,L,coeffs[s]);
        free(w);
    }

    // Thresholds per scale
    double sigma[N_SCALES],mu[N_SCALES];
    for(int s=0;s<N_SCALES;s++){
        double *absvals=malloc(N*sizeof(double));
        for(int i=0;i<N;i++)absvals[i]=fabs(coeffs[s][i]);
        mu[s]=median_double(absvals,N);
        sigma[s]=compute_mad(coeffs[s],N);
        free(absvals);
    }

    unsigned char *pos_mask=calloc(N,sizeof(unsigned char));
    unsigned char *neg_mask=calloc(N,sizeof(unsigned char));
    for(int i=0;i<N;i++){
        for(int s=0;s<N_SCALES;s++){
            double thr=1.2*mu[s]+THRESH_K*sigma[s];
            if(coeffs[s][i]>thr)pos_mask[i]=1;
            if(coeffs[s][i]<-thr)neg_mask[i]=1;
        }
    }

    double *avg=malloc(N*sizeof(double));
    for(int i=0;i<N;i++){avg[i]=0;for(int s=0;s<N_SCALES;s++)avg[i]+=coeffs[s][i];avg[i]/=N_SCALES;}

    int *pos_idx=malloc(N*sizeof(int));
    int *neg_idx=malloc(N*sizeof(int));
    int pos_n=detect_regions(avg,pos_mask,N,pos_idx);
    int neg_n=detect_regions(avg,neg_mask,N,neg_idx);

    double *pos_t=malloc(pos_n*sizeof(double));
    double *pos_m=malloc(pos_n*sizeof(double));
    for(int i=0;i<pos_n;i++){pos_t[i]=tsec[pos_idx[i]]*1000.0;pos_m[i]=amp[pos_idx[i]];}
    double *neg_t=malloc(neg_n*sizeof(double));
    double *neg_m=malloc(neg_n*sizeof(double));
    for(int i=0;i<neg_n;i++){neg_t[i]=tsec[neg_idx[i]]*1000.0;neg_m[i]=amp[neg_idx[i]];}

    double *ptm=malloc(pos_n*sizeof(double)),*pmm=malloc(pos_n*sizeof(double));
    double *ntm=malloc(neg_n*sizeof(double)),*nmm=malloc(neg_n*sizeof(double));
    int pos_final=merge_events(pos_t,pos_m,pos_n,MERGE_WINDOW_MS,ptm,pmm);
    int neg_final=merge_events(neg_t,neg_m,neg_n,MERGE_WINDOW_MS,ntm,nmm);

    // --- Output ---
    for(int i=0;i<pos_final;i++){
        double target=ptm[i]/1000.0;
        int idx=0;
        for(int j=1;j<N;j++) if(fabs(tsec[j]-target)<fabs(tsec[idx]-target)) idx=j;
        if (pmm[i] < 0.3) continue;  // Skip weak spikes (<0.3 m/s²)
        printf("Positive spike at %s (%.3f m/s²)\n",samples[idx].time_str,pmm[i]);
        addCWTDetection("Positive", samples[idx].time_str, pmm[i]);
    }
    for(int i=0;i<neg_final;i++){
        double target=ntm[i]/1000.0;
        int idx=0;
        for(int j=1;j<N;j++) if(fabs(tsec[j]-target)<fabs(tsec[idx]-target)) idx=j;
        if (fabs(nmm[i]) < 0.3) continue;
        printf("Negative spike at %s (%.3f m/s²)\n",samples[idx].time_str,nmm[i]);
        addCWTDetection("Negative", samples[idx].time_str, pmm[i]);
    }

    // cleanup
    free(samples);free(tsec);free(amp);free(dts);
    for(int s=0;s<N_SCALES;s++)free(coeffs[s]);free(coeffs);
    free(pos_mask);free(neg_mask);
    free(avg);free(pos_idx);free(neg_idx);
    free(pos_t);free(pos_m);free(neg_t);free(neg_m);
    free(ptm);free(pmm);free(ntm);free(nmm);

    printf("\nAdaptive CWT thresholding applied.\n");
    printf("----------------------------------------------------------\n");
    return 0;
}

// --- Core pipeline ---
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

static Detection cwtDet[MAX_DETECTIONS];
static int cwtCount = 0;

static void addCWTDetection(const char *type, const char *timeStr, float amplitude) {
    if (cwtCount >= MAX_DETECTIONS) return;
    Detection *d = &cwtDet[cwtCount++];
    strcpy(d->polarity, type);
    int h, m, s, ms;
    sscanf(timeStr, "%d:%d:%d:%d", &h, &m, &s, &ms);
    d->hour = h; d->min = m; d->sec = s; d->ms = ms;
    d->amplitude = amplitude;
}



int getCWTDetections(Detection **array) {
    *array = cwtDet;
    return cwtCount;
}