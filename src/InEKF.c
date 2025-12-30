#include "InEKF.h"
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdlib.h>

/* ---------- Simple logging ---------- */
static FILE *g_log_fp = NULL;
static void log_open(const char *path) {
    if (!path) return;
    g_log_fp = fopen(path, "w");
}
static void log_close(void) {
    if (g_log_fp) { fclose(g_log_fp); g_log_fp = NULL; }
}
static void log_printf(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    if (g_log_fp) {
        va_list ap2;
        va_copy(ap2, ap);
        vfprintf(g_log_fp, fmt, ap2);
        va_end(ap2);
    }
    va_end(ap);
}

/* ---------- Linear algebra helpers ---------- */
static void mat_identity(double *A, int n) {
    for (int i = 0; i < n*n; ++i) A[i] = 0.0;
    for (int i = 0; i < n; ++i) A[i*n + i] = 1.0;
}

static void mat_zero_elems(double *A, int elems) {
    for (int i = 0; i < elems; ++i) A[i] = 0.0;
}

static void mat_copy(double *dst, const double *src, int n) {
    for (int i=0; i<n*n; i++) dst[i] = src[i];
}

static void mat3_mul(double *C, const double *A, const double *B) {
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            double s = 0.0;
            for(int k = 0; k < 3; ++k) s += A[i*3+k] * B[k*3+j];
            C[i*3+j] = s;
        }
    }
}

static void mat3_vec(double *y, const double *A, const double *x) {
    for(int i = 0; i < 3; ++i) {
        y[i] = A[i*3+0]*x[0] + A[i*3+1]*x[1] + A[i*3+2]*x[2];
    }
}

static void mat_transpose(double *B, const double *A, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            B[i*n + j] = A[j*n + i];
}

static double norm3(const double *v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static void skew3(double *M, const double *w) {
    M[0] = 0;      M[1] = -w[2]; M[2] = w[1];
    M[3] = w[2];   M[4] = 0;     M[5] = -w[0];
    M[6] = -w[1];  M[7] = w[0];  M[8] = 0;
}

static void expSO3(double *R, const double *w, double dt) {
    double phi[3] = { w[0]*dt, w[1]*dt, w[2]*dt };
    double th = norm3(phi);
    double I3[9] = {1,0,0, 0,1,0, 0,0,1};
    double K[9];

    if (th < 1e-12) {
        skew3(K, phi);
        for(int i=0;i<9;i++) R[i] = I3[i] + K[i];
        return;
    }

    double u[3] = {phi[0]/th, phi[1]/th, phi[2]/th};
    skew3(K, u);
    double K2[9];
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) {
        double s=0;
        for (int k=0;k<3;k++) s += K[i*3+k]*K[k*3+j];
        K2[i*3+j]=s;
    }
    double c = cos(th), s = sin(th);
    for(int i=0;i<9;i++) R[i] = I3[i] + s * K[i] + (1.0 - c) * K2[i];
}

/* ---------- InEKF core functions ---------- */

void inekf_init(InEKF *filter, const InEKFConfig *config) {
    if (!filter || !config) return;

    // Copy config
    filter->config = *config;
    
    // Initialize state
    mat_identity(filter->state.R, 3);
    for (int i=0;i<3;i++) {
        filter->state.v[i] = 0.0;
        filter->state.p[i] = 0.0;
        filter->state.bw[i] = 0.0;
        filter->state.ba[i] = 0.0;
    }

    // Initialize covariance
    mat_zero_elems(filter->state.P, NSTATE*NSTATE);
    for (int i=0;i<3;i++) filter->state.P[i*NSTATE + i] = 1e-4;
    for (int i=3;i<6;i++) filter->state.P[i*NSTATE + i] = 1e-3;
    for (int i=6;i<9;i++) filter->state.P[i*NSTATE + i] = 1e-3;
    for (int i=9;i<12;i++) filter->state.P[i*NSTATE + i] = 1e-6;
    for (int i=12;i<15;i++) filter->state.P[i*NSTATE + i] = 1e-6;

    // Initialize process noise
    mat_zero_elems(filter->Qc, NSTATE*NSTATE);
    for (int i=0;i<3;i++) filter->Qc[i*NSTATE + i] = config->gyro_noise_sigma * config->gyro_noise_sigma;
    for (int i=3;i<6;i++) filter->Qc[i*NSTATE + i] = config->accel_noise_sigma * config->accel_noise_sigma;
    for (int i=9;i<12;i++) filter->Qc[i*NSTATE + i] = config->gyro_bias_rw * config->gyro_bias_rw;
    for (int i=12;i<15;i++) filter->Qc[i*NSTATE + i] = config->accel_bias_rw * config->accel_bias_rw;

    // Initialize measurement noise
    mat_zero_elems(filter->R_zupt, 9);
    double noise = config->zupt_velocity_noise;
    filter->R_zupt[0] = noise*noise;
    filter->R_zupt[4] = noise*noise;
    filter->R_zupt[8] = noise*noise;
}

void inekf_init_default(InEKF *filter) {
    InEKFConfig cfg = {
        .gyro_noise_sigma = 0.003,
        .accel_noise_sigma = 0.05,
        .gyro_bias_rw = 1e-5,
        .accel_bias_rw = 1e-3,
        .zupt_velocity_noise = 0.01,
        .zupt_accel_threshold = 0.1,
        .zupt_gyro_threshold_deg = 0.5
    };
    inekf_init(filter, &cfg);
}

static void build_F_G(double *F, double *G, const double *R, const double *v, const double *p) {
    mat_zero_elems(F, NSTATE*NSTATE);
    mat_zero_elems(G, NSTATE*NSTATE);

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) F[i*NSTATE + (9+j)] = -R[i*3+j];

    double gskew[9] = {0, -GGRAV, 0, GGRAV, 0, 0, 0, 0, 0};
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) F[(3+i)*NSTATE + j] = gskew[i*3+j];

    double vsk[9]; skew3(vsk, v);
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        double s=0;
        for (int k=0;k<3;k++) s += vsk[i*3+k]*R[k*3+j];
        F[(3+i)*NSTATE + (9+j)] = -s;
    }

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) F[(3+i)*NSTATE + (12+j)] = -R[i*3+j];

    for (int i=0;i<3;i++) F[(6+i)*NSTATE + (3+i)] = 1.0;

    double psk[9]; skew3(psk, p);
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        double s=0;
        for (int k=0;k<3;k++) s += psk[i*3+k]*R[k*3+j];
        F[(6+i)*NSTATE + (9+j)] = -s;
    }

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) G[i*NSTATE + (0+j)] = R[i*3+j];

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        double s = 0;
        for (int k=0;k<3;k++) s += vsk[i*3+k]*R[k*3+j];
        G[(3+i)*NSTATE + (0+j)] = s;
        G[(3+i)*NSTATE + (3+j)] = R[i*3+j];
    }

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        double s=0;
        for (int k=0;k<3;k++) s += psk[i*3+k]*R[k*3+j];
        G[(6+i)*NSTATE + (0+j)] = s;
        G[(6+i)*NSTATE + (6+j)] = R[i*3+j];
    }

    for (int i=0;i<3;i++) G[(9+i)*NSTATE + (9+i)] = -1.0;
    for (int i=0;i<3;i++) G[(12+i)*NSTATE + (12+i)] = -1.0;
}

static void compute_Phi_firstorder(double *Phi, const double *F, double dt) {
    mat_identity(Phi, NSTATE);
    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) Phi[i*NSTATE + j] += F[i*NSTATE + j] * dt;
}

static void propagate_covariance(double *Sigma, const double *F, const double *G, const double *Qc, double dt) {
     double Phi[NSTATE*NSTATE];
     double tmp1[NSTATE*NSTATE];
     double tmp2[NSTATE*NSTATE];
     double tmp3[NSTATE*NSTATE];
     double GQ[NSTATE*NSTATE];
     double GQGT[NSTATE*NSTATE];

    compute_Phi_firstorder(Phi, F, dt);

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += Phi[i*NSTATE + k] * Sigma[k*NSTATE + j];
        tmp1[i*NSTATE + j] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += tmp1[i*NSTATE + k] * Phi[j*NSTATE + k];
        tmp2[i*NSTATE + j] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += G[i*NSTATE + k] * Qc[k*NSTATE + j];
        GQ[i*NSTATE + j] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += GQ[i*NSTATE + k] * G[j*NSTATE + k];
        GQGT[i*NSTATE + j] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += Phi[i*NSTATE + k] * GQGT[k*NSTATE + j];
        tmp1[i*NSTATE + j] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += tmp1[i*NSTATE + k] * Phi[j*NSTATE + k];
        tmp3[i*NSTATE + j] = s * dt;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) 
        Sigma[i*NSTATE + j] = tmp2[i*NSTATE + j] + tmp3[i*NSTATE + j];
}

static void cov_update_joseph(double *P, const double *K, const double *H, const double *Rmeas, int m) {
    double I_N[NSTATE*NSTATE];
    double KH[NSTATE*NSTATE];
    double IminusKH[NSTATE*NSTATE];
    double tmp[NSTATE*NSTATE];
    double Pnew[NSTATE*NSTATE];
    double KR[NSTATE * 3];

    mat_identity(I_N, NSTATE);

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s = 0.0;
        for (int r = 0; r < m; ++r) {
            s += K[i*NSTATE + r] * H[r * NSTATE + j];
        }
        KH[i*NSTATE + j] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) 
        IminusKH[i*NSTATE + j] = I_N[i*NSTATE + j] - KH[i*NSTATE + j];

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += IminusKH[i*NSTATE + k] * P[k*NSTATE + j];
        tmp[i*NSTATE + j] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += tmp[i*NSTATE + k] * IminusKH[j*NSTATE + k];
        Pnew[i*NSTATE + j] = s;
    }

    mat_zero_elems(KR, NSTATE * 3);
    for (int i=0;i<NSTATE;i++) for (int r=0;r<m;r++) {
        double s=0;
        for (int c=0;c<m;c++) s += K[i*NSTATE + c] * Rmeas[c*3 + r];
        KR[i*3 + r] = s;
    }

    for (int i=0;i<NSTATE;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int r=0;r<m;r++) s += KR[i*3 + r] * K[j*NSTATE + r];
        Pnew[i*NSTATE + j] += s;
    }

    for (int i=0;i<NSTATE*NSTATE;i++) P[i] = Pnew[i];
}

void inekf_predict(InEKF *filter, const double omega_meas[3], const double acc_meas[3], double dt) {
    if (!filter || dt <= 0) return;

    double omega_corr[3], acc_corr[3];
    for (int i=0;i<3;i++) {
        omega_corr[i] = omega_meas[i] - filter->state.bw[i];
        acc_corr[i] = acc_meas[i] - filter->state.ba[i];
    }

    double Rinc[9];
    expSO3(Rinc, omega_corr, dt);

    double Rnew[9];
    mat3_mul(Rnew, filter->state.R, Rinc);
    mat_copy(filter->state.R, Rnew, 3);

    double tmpvec[3], a_world[3];
    mat3_vec(tmpvec, filter->state.R, acc_corr);
    a_world[0] = tmpvec[0];
    a_world[1] = tmpvec[1];
    // a_world[2] = tmpvec[2] + GGRAV;
    a_world[2] = tmpvec[2] - GGRAV;

    for (int i=0;i<3;i++) filter->state.v[i] += a_world[i] * dt;
    for (int i=0;i<3;i++) filter->state.p[i] += filter->state.v[i] * dt;

    double F[NSTATE*NSTATE], G[NSTATE*NSTATE];
    build_F_G(F, G, filter->state.R, filter->state.v, filter->state.p);
    propagate_covariance(filter->state.P, F, G, filter->Qc, dt);
}

static void measurement_update_velocity(InEKF *filter, const double z[3]) {
    double Rt[9];
    mat_transpose(Rt, filter->state.R, 3);
    double v_in_imu[3];
    mat3_vec(v_in_imu, Rt, filter->state.v);

    double y[3];
    for (int i=0;i<3;i++) y[i] = z[i] - v_in_imu[i];

    double H[3 * NSTATE];
    for (int i=0;i<3*NSTATE;i++) H[i] = 0.0;
    for (int r=0;r<3;r++) H[r*NSTATE + (3 + r)] = 1.0;

    double HP[3 * NSTATE];
    for (int i=0;i<3;i++) for (int j=0;j<NSTATE;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += H[i*NSTATE + k] * filter->state.P[k*NSTATE + j];
        HP[i*NSTATE + j] = s;
    }

    double S[9];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += HP[i*NSTATE + k] * H[j*NSTATE + k];
        S[i*3+j] = s + filter->R_zupt[i*3+j];
    }

    double detS = S[0]*(S[4]*S[8]-S[5]*S[7]) - S[1]*(S[3]*S[8]-S[5]*S[6]) + S[2]*(S[3]*S[7]-S[4]*S[6]);
    double Sinv[9];
    if (fabs(detS) < 1e-12) {
        for (int i=0;i<9;i++) Sinv[i] = 0.0;
        for (int d=0; d<3; d++) if (fabs(S[d*3+d])>1e-12) Sinv[d*3+d] = 1.0 / S[d*3+d];
    } else {
        Sinv[0] =  (S[4]*S[8]-S[5]*S[7]) / detS;
        Sinv[1] = -(S[1]*S[8]-S[2]*S[7]) / detS;
        Sinv[2] =  (S[1]*S[5]-S[2]*S[4]) / detS;
        Sinv[3] = -(S[3]*S[8]-S[5]*S[6]) / detS;
        Sinv[4] =  (S[0]*S[8]-S[2]*S[6]) / detS;
        Sinv[5] = -(S[0]*S[5]-S[2]*S[3]) / detS;
        Sinv[6] =  (S[3]*S[7]-S[4]*S[6]) / detS;
        Sinv[7] = -(S[0]*S[7]-S[1]*S[6]) / detS;
        Sinv[8] =  (S[0]*S[4]-S[1]*S[3]) / detS;
    }

    double PHT[NSTATE * 3];
    for (int i=0;i<NSTATE;i++) for (int j=0;j<3;j++) {
        double s=0;
        for (int k=0;k<NSTATE;k++) s += filter->state.P[i*NSTATE + k] * H[j*NSTATE + k];
        PHT[i*3 + j] = s;
    }

    double K[NSTATE * 3];
    for (int i=0;i<NSTATE;i++) for (int j=0;j<3;j++) {
        double s=0;
        for (int k=0;k<3;k++) s += PHT[i*3 + k] * Sinv[k*3 + j];
        K[i*3 + j] = s;
    }

    double delta[NSTATE];
    for (int i=0;i<NSTATE;i++) {
        double s=0;
        for (int j=0;j<3;j++) s += K[i*3 + j] * y[j];
        delta[i] = s;
    }

    double dphi[3] = { delta[0], delta[1], delta[2] };
    double Rcor[9];
    expSO3(Rcor, dphi, 1.0);
    double Rupdated[9];
    mat3_mul(Rupdated, Rcor, filter->state.R);
    mat_copy(filter->state.R, Rupdated, 3);

    for (int i=0;i<3;i++) {
        filter->state.v[i] += delta[3+i];
        filter->state.p[i] += delta[6+i];
        filter->state.bw[i] += delta[9+i];
        filter->state.ba[i] += delta[12+i];
    }

    cov_update_joseph(filter->state.P, K, H, filter->R_zupt, 3);
}

void inekf_zupt_update(InEKF *filter, const double velocity_meas[3]) {
    if (!filter) return;
    measurement_update_velocity(filter, velocity_meas);
}

int inekf_detect_zupt(const InEKF *filter, const double acc_meas[3], const double gyro_meas_deg[3]) {
    if (!filter) return 0;
    double acc_norm = norm3(acc_meas);
    double gyro_norm_deg = norm3(gyro_meas_deg);
    if (fabs(acc_norm - GGRAV) < filter->config.zupt_accel_threshold &&
        gyro_norm_deg < filter->config.zupt_gyro_threshold_deg) return 1;
    return 0;
}

/* Accessors */
void inekf_get_rotation(const InEKF *filter, double R[9]) { 
    if (filter && R) mat_copy(R, filter->state.R, 3); 
}
void inekf_get_velocity(const InEKF *filter, double v[3]) { 
    if (filter && v) for (int i=0;i<3;i++) v[i]=filter->state.v[i]; 
}
void inekf_get_position(const InEKF *filter, double p[3]) { 
    if (filter && p) for (int i=0;i<3;i++) p[i]=filter->state.p[i]; 
}
void inekf_get_biases(const InEKF *filter, double bw[3], double ba[3]) {
    if (!filter) return;
    if (bw) for (int i=0;i<3;i++) bw[i] = filter->state.bw[i];
    if (ba) for (int i=0;i<3;i++) ba[i] = filter->state.ba[i];
}

void inekf_get_bias_corrected_imu(const InEKF *filter, const double acc_meas[3], const double gyro_meas_deg[3],
                                  double acc_corrected[3], double gyro_corrected_deg[3]) {
    if (!filter) return;
    if (acc_corrected) for (int i=0;i<3;i++) acc_corrected[i] = acc_meas[i] - filter->state.ba[i];
    if (gyro_corrected_deg) {
        for (int i=0;i<3;i++) {
            double gyro_rad = gyro_meas_deg[i] * DEG2RAD;
            double corrected_rad = gyro_rad - filter->state.bw[i];
            gyro_corrected_deg[i] = corrected_rad * RAD2DEG;
        }
    }
}

int run_InEKF(const char* input_file, const char* output_file) {
    
    if (!input_file || !output_file) {
        return 1;
    }
    
    // Open files
    FILE *fin = fopen(input_file, "r");
    if (!fin) {
        return 1;
    }
    
    FILE *fout = fopen(output_file, "w");
    if (!fout) {
        fclose(fin);
        return 1;
    }
    
    // Read and copy header
    char header[1024];
    if (!fgets(header, sizeof(header), fin)) {
        fclose(fin);
        fclose(fout);
        return 1;
    }
    fprintf(fout, "%s", header);
    
    // Initialize filter
    InEKF filter;
    inekf_init_default(&filter);
    
    // Process all data lines
    char line[4096];
    int line_count = 0;
    int processed_count = 0;
    double prev_time = -1.0;
    
    while (fgets(line, sizeof(line), fin)) {
        line_count++;
        
        // Remove newline characters
        line[strcspn(line, "\r\n")] = 0;
        
        // Skip empty lines
        if (strlen(line) < 5) {
            continue;
        }
        
        // Parse CSV line
        char time_str[128];
        double ax, ay, az, gx, gy, gz;
        
        if (sscanf(line, "%127[^,],%lf,%lf,%lf,%lf,%lf,%lf", 
                   time_str, &ax, &ay, &az, &gx, &gy, &gz) != 7) {
            continue;
        }
        
        // Parse actual timestamp from time_str
        double t;
        if (sscanf(time_str, "%lf", &t) != 1) {
            // If timestamp parsing fails, use synthetic time
            t = (prev_time < 0) ? 0.0 : prev_time + 0.01;
        }
        
        // Calculate delta time
        double dt;
        if (prev_time < 0) {
            dt = 0.01; // Default for first sample
        } else {
            dt = t - prev_time;
            // Handle invalid or zero dt
            if (dt <= 0) dt = 0.01;
        }
        prev_time = t;
        
        // Convert to correct units
        double gyro_rad[3] = { gx * DEG2RAD, gy * DEG2RAD, gz * DEG2RAD };
        double acc_meas[3] = { ax, ay, az };
        double gyro_deg[3] = { gx, gy, gz };
        
        // Predict step
        inekf_predict(&filter, gyro_rad, acc_meas, dt);
        
        // ZUPT detection and update
        if (inekf_detect_zupt(&filter, acc_meas, gyro_deg)) {
            double zero_v[3] = {0.0, 0.0, 0.0};
            inekf_zupt_update(&filter, zero_v);
        }
        
        // Get bias-corrected data
        double acc_corr[3], gyro_corr[3];
        inekf_get_bias_corrected_imu(&filter, acc_meas, gyro_deg, acc_corr, gyro_corr);
        
        // Write output
        fprintf(fout, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                time_str,
                acc_corr[0], acc_corr[1], acc_corr[2],
                gyro_corr[0], gyro_corr[1], gyro_corr[2]);
        
        processed_count++;
        
        // Progress reporting and flush more frequently
        if (processed_count % 100 == 0) {
            fflush(fout);
        }
    }
    
    fclose(fin);
    fclose(fout);
    
    printf("Processed %d lines successfully\n", processed_count);
    
    return 0;
}