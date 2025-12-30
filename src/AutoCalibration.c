#include "AutoCalibration.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define G_ACCEL 9.80665
#define PARAMS 9
#define MAX_SAMPLES 20000
#define EPS_FD 1e-6
#define MAX_ITERS 200
#define TOL_REL 1e-9

// --- Global parameters (can be used anywhere in the program) ---
double Ox, Oy, Oz;
double Sxx, Syy, Szz, Sxy, Sxz, Syz;

// ---------- Utility functions ----------
static void *xmalloc(size_t s) {
    void *p = malloc(s);
    if (!p) { fprintf(stderr, "malloc failed\n"); exit(1); }
    return p;
}

// --- Simple linear algebra helpers ---
static void matmul_AT_A(double *A, int m, int n, double *out) {
    for (int i=0;i<n*n;i++) out[i]=0.0;
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            double a = A[i*n + j];
            for (int k=j;k<n;k++){
                out[j*n + k] += a * A[i*n + k];
            }
        }
    }
    for (int i=0;i<n;i++) for (int j=0;j<i;j++) out[i*n + j] = out[j*n + i];
}

static void matvec_AT_r(double *A, double *r, int m, int n, double *out) {
    for (int j=0;j<n;j++){
        double s=0;
        for (int i=0;i<m;i++) s += A[i*n + j]*r[i];
        out[j]=s;
    }
}

static int gauss_solve(double *A, double *b, int n) {
    for (int k=0;k<n;k++){
        int piv=k;
        double maxv=fabs(A[k*n+k]);
        for (int i=k+1;i<n;i++){
            double v=fabs(A[i*n+k]);
            if(v>maxv){maxv=v;piv=i;}
        }
        if(maxv<1e-18)return -1;
        if(piv!=k){
            for(int j=k;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[piv*n+j];A[piv*n+j]=tmp;}
            double tmpb=b[k];b[k]=b[piv];b[piv]=tmpb;
        }
        double akk=A[k*n+k];
        for(int j=k+1;j<n;j++)A[k*n+j]/=akk;
        b[k]/=akk;
        for(int i=k+1;i<n;i++){
            double f=A[i*n+k];
            if(f==0.0)continue;
            for(int j=k+1;j<n;j++)A[i*n+j]-=f*A[k*n+j];
            b[i]-=f*b[k];
        }
    }
    for(int i=n-1;i>=0;i--){
        double s=b[i];
        for(int j=i+1;j<n;j++)s-=A[i*n+j]*b[j];
        b[i]=s;
    }
    return 0;
}

// --- Residual evaluation ---
static void eval_residuals(const double *params, double *V, int N, double *r) {
    double Ox=params[0],Oy=params[1],Oz=params[2];
    double Sxx=params[3],Syy=params[4],Szz=params[5];
    double Sxy=params[6],Sxz=params[7],Syz=params[8];
    for(int k=0;k<N;k++){
        double vx=V[3*k],vy=V[3*k+1],vz=V[3*k+2];
        double dx=vx-Ox,dy=vy-Oy,dz=vz-Oz;
        double ax=Sxx*dx + Sxy*dy + Sxz*dz;
        double ay=Sxy*dx + Syy*dy + Syz*dz;
        double az=Sxz*dx + Syz*dy + Szz*dz;
        r[k]=ax*ax + ay*ay + az*az - G_ACCEL*G_ACCEL;
    }
}

static void build_jacobian_fd(const double *params,double *V,int N,double *J,double *base_r){
    eval_residuals(params,V,N,base_r);
    for(int p=0;p<PARAMS;p++){
        double eps=EPS_FD*(1.0+fabs(params[p]));
        if(eps==0.0)eps=EPS_FD;
        double *tmp=(double*)xmalloc(sizeof(double)*PARAMS);
        memcpy(tmp,params,sizeof(double)*PARAMS);
        tmp[p]+=eps;
        double *r_eps=(double*)xmalloc(sizeof(double)*N);
        eval_residuals(tmp,V,N,r_eps);
        for(int i=0;i<N;i++)J[i*PARAMS+p]=(r_eps[i]-base_r[i])/eps;
        free(tmp);free(r_eps);
    }
}

static double compute_cost(double *r,int N){
    double s=0;for(int i=0;i<N;i++)s+=r[i]*r[i];
    return s/(double)N;
}

// --- Levenberg-Marquardt optimization ---
static int optimize_params(double *params,double *V,int N,double *out_cost){
    double *r=xmalloc(sizeof(double)*N);
    double *J=xmalloc(sizeof(double)*N*PARAMS);
    double *base_r=xmalloc(sizeof(double)*N);
    eval_residuals(params,V,N,r);
    double cost=compute_cost(r,N);
    double lambda=1e-3;
    for(int iter=0;iter<MAX_ITERS;iter++){
        build_jacobian_fd(params,V,N,J,base_r);
        double H[PARAMS*PARAMS]; matmul_AT_A(J,N,PARAMS,H);
        double JT_r[PARAMS]; matvec_AT_r(J,r,N,PARAMS,JT_r);
        double H_lm[PARAMS*PARAMS]; memcpy(H_lm,H,sizeof(H_lm));
        for(int i=0;i<PARAMS;i++)H_lm[i*PARAMS+i]+=lambda*(H_lm[i*PARAMS+i]+1e-12);
        double rhs[PARAMS]; for(int i=0;i<PARAMS;i++)rhs[i]=-JT_r[i];
        double Hc[PARAMS*PARAMS]; memcpy(Hc,H_lm,sizeof(Hc));
        double rhsc[PARAMS]; memcpy(rhsc,rhs,sizeof(rhsc));
        if(gauss_solve(Hc,rhsc,PARAMS)!=0)return -1;
        double delta[PARAMS]; for(int i=0;i<PARAMS;i++)delta[i]=rhsc[i];
        double params_new[PARAMS]; for(int i=0;i<PARAMS;i++)params_new[i]=params[i]+delta[i];
        double *r_new=xmalloc(sizeof(double)*N); eval_residuals(params_new,V,N,r_new);
        double cost_new=compute_cost(r_new,N);
        if(cost_new<cost){ memcpy(params,params_new,sizeof(params_new)); memcpy(r,r_new,sizeof(double)*N); cost=cost_new; lambda*=0.5; if(lambda<1e-12)lambda=1e-12; }
        else lambda*=10.0;
        free(r_new);
        double max_rel=0.0;
        for(int i=0;i<PARAMS;i++){ double denom=fabs(params[i])+1e-12; double rel=fabs(delta[i])/denom; if(rel>max_rel)max_rel=rel; }
        if(max_rel<TOL_REL)break;
    }
    *out_cost=cost;
    free(r);free(J);free(base_r);
    return 0;
}

// --- Simple static detection using std of magnitude ---
static void mean_std(const double *x,int start,int len,double *mean,double *std){
    double s=0;for(int i=0;i<len;i++)s+=x[start+i];*mean=s/len;
    double ss=0;for(int i=0;i<len;i++){double d=x[start+i]-*mean;ss+=d*d;}*std=sqrt(ss/len);
}

// Process one file (Time,Ax,Ay,Az,Gx,Gy,Gz)
static int process_file(const char *fname,int window_size,double std_thresh,double *samples,int max_samples){
    FILE *f=fopen(fname,"r"); if(!f) return 0;
    int cap=4096,n=0;
    double *Ax=xmalloc(sizeof(double)*cap),*Ay=xmalloc(sizeof(double)*cap),*Az=xmalloc(sizeof(double)*cap);
    char line[512];
    while(fgets(line,sizeof(line),f)){
        if(line[0]=='#')continue;
        double t,ax,ay,az,gx,gy,gz;
        // int got=sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&t,&ax,&ay,&az,&gx,&gy,&gz);
        char time_str[32];
        int got = sscanf(line, "%[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                 time_str, &ax, &ay, &az, &gx, &gy, &gz);
        if(got<4)continue;
        if(n>=cap){cap*=2;Ax=realloc(Ax,sizeof(double)*cap);Ay=realloc(Ay,sizeof(double)*cap);Az=realloc(Az,sizeof(double)*cap);}
        Ax[n]=ax;Ay[n]=ay;Az[n]=az;n++;
    } fclose(f);
    if(n<window_size){free(Ax);free(Ay);free(Az);return 0;}
    double *mag=xmalloc(sizeof(double)*n);
    for(int i=0;i<n;i++)mag[i]=sqrt(Ax[i]*Ax[i]+Ay[i]*Ay[i]+Az[i]*Az[i]);
    char *is_static=xmalloc(n-window_size+1);
    for(int i=0;i<=n-window_size;i++){double m,s;mean_std(mag,i,window_size,&m,&s);is_static[i]=(s<std_thresh)?1:0;}
    int added=0;
    for(int i=0;i<=n-window_size && added<max_samples;i++){
        if(!is_static[i])continue;
        int j=i;while(j<=n-window_size && is_static[j])j++;
        int start=i,end=j+window_size-2;if(end>=n)end=n-1;
        int len=end-start+1;
        if(len>=window_size){
            double sx=0,sy=0,sz=0;for(int k=start;k<=end;k++){sx+=Ax[k];sy+=Ay[k];sz+=Az[k];}
            samples[3*added+0]=sx/len; samples[3*added+1]=sy/len; samples[3*added+2]=sz/len; added++;
        }
        i=j;
    }
    free(Ax);free(Ay);free(Az);free(mag);free(is_static);
    return added;
}

// --- Public API implementation ---
int run_autocalibration(const char *folder_path,int start_idx,int end_idx,
                        int window_size,double std_frac)
{
    double std_thresh=std_frac*G_ACCEL;
    double *samples=xmalloc(sizeof(double)*MAX_SAMPLES*3);
    int total=0;
    char fname[512];
    for(int i=start_idx;i<=end_idx && total<MAX_SAMPLES;i++){
        snprintf(fname,sizeof(fname),"%s/%d.csv",folder_path,i);
        int added=process_file(fname,window_size,std_thresh,samples+3*total,MAX_SAMPLES-total);
        total+=added;
    }
    if(total<9){fprintf(stderr,"Not enough static samples (%d)\n",total);free(samples);return -1;}
    double params[PARAMS]={0,0,0,1,1,1,0,0,0};
    double cost;
    if(optimize_params(params,samples,total,&cost)!=0){fprintf(stderr,"Optimization failed\n");free(samples);return -2;}
    Ox=params[0];Oy=params[1];Oz=params[2];
    Sxx=params[3];Syy=params[4];Szz=params[5];
    Sxy=params[6];Sxz=params[7];Syz=params[8];
    printf("Calibration complete: cost=%.6g, samples=%d\n",cost,total);
    printf("O=(%.6g, %.6g, %.6g)\n",Ox,Oy,Oz);
    printf("S=\n[%.6g %.6g %.6g]\n[%.6g %.6g %.6g]\n[%.6g %.6g %.6g]\n",
           Sxx,Sxy,Sxz,Sxy,Syy,Syz,Sxz,Syz,Szz);
    free(samples);
    return 0;
}

/*Helper FUnctions*/
// ---- compute phi/rho (radians) ----
static void compute_tilt_angles_rad(double Ax, double Ay, double Az,
                                    double *phi_rad, double *rho_rad)
{
    *phi_rad = atan2(Ax, sqrt(Ay*Ay + Az*Az));
    *rho_rad = atan2(Ay, sqrt(Ax*Ax + Az*Az));
}

static void build_R_xyz(double phi_rad, double rho_rad, double R[9])
{
    double cphi = cos(phi_rad);
    double sphi = sin(phi_rad);
    double cth  = cos(rho_rad);
    double sth  = sin(rho_rad);

    R[0] =  cth;
    R[1] =  0.0;
    R[2] = -sth;
    R[3] =  sphi * sth;
    R[4] =  cphi;
    R[5] =  sphi * cth;
    R[6] =  cphi * sth;
    R[7] = -sphi;
    R[8] =  cphi * cth;
}

static void build_R_yxz(double phi_rad, double rho_rad, double R[9])
{
    double cphi = cos(phi_rad);
    double sphi = sin(phi_rad);
    double cth  = cos(rho_rad);
    double sth  = sin(rho_rad);

    R[0] =  cth;
    R[1] =  sth * sphi;
    R[2] = -sth * cphi;
    R[3] =  0.0;
    R[4] =  cphi;
    R[5] =  sphi;
    R[6] =  sth;
    R[7] = -cth * sphi;
    R[8] =  cth * cphi;
}



int apply_calibration_to_csv(const char *input_csv, const char *output_csv)
{
    FILE *fin = fopen(input_csv, "r");
    if (!fin) {
        fprintf(stderr, "Error opening input file: %s\n", input_csv);
        return -1;
    }

    FILE *fout = fopen(output_csv, "w");
    if (!fout) {
        fprintf(stderr, "Error creating output file: %s\n", output_csv);
        fclose(fin);
        return -1;
    }

    char line[512];

    // Copy header if present
    if (fgets(line, sizeof(line), fin)) {
        if (strstr(line, "Ax") || strstr(line, "Time")) {
            // Standard header
            fprintf(fout, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");
        } else {
            // Not a header, rewind file to process this line
            fseek(fin, 0, SEEK_SET);
            fprintf(fout, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");
        }
    }

    while (fgets(line, sizeof(line), fin)) {
        if (line[0] == '#' || strlen(line) < 3)
            continue;

        char time_str[64];
        double Ax_raw, Ay_raw, Az_raw, Gx, Gy, Gz;

        // Expect 7 columns: Time, Ax, Ay, Az, Gx, Gy, Gz
        int got = sscanf(line, "%[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                         time_str, &Ax_raw, &Ay_raw, &Az_raw, &Gx, &Gy, &Gz);
        if (got < 7)
            continue; // Skip malformed lines

        // --- Apply calibration: A = S * (V - O) ---
        double dx = Ax_raw - Ox;
        double dy = Ay_raw + Oy;
        double dz = Az_raw + Oz;

        double Ax_cal = Sxx * dx + Sxy * dy + Sxz * dz;
        double Ay_cal = Sxy * dx + Syy * dy + Syz * dz;
        double Az_cal = Sxz * dx + Syz * dy + Szz * dz;

        // ---- Compute tilt angles (radians) ----
        double phi_rad, rho_rad;
        compute_tilt_angles_rad(Ax_cal, Ay_cal, Az_cal, &phi_rad, &rho_rad);

        // ---- Evaluate condition ----
        double lhs = fabs(sqrt(Ay_cal*Ay_cal + Az_cal*Az_cal )- 1.0);
        double rhs = fabs(sqrt(Ax_cal*Ax_cal + Az_cal*Az_cal )- 1.0);

        // ---- Choose rotation matrix ----
        double R[9];
        if (lhs < rhs) {
            build_R_xyz(phi_rad, rho_rad, R);   // Condition true
        } else {
            build_R_yxz(phi_rad, rho_rad, R);   // Condition false
        }

        // ---- Apply rotation to accelerometer vector ----
        double Ax_rot = R[0]*Ax_cal + R[1]*Ay_cal + R[2]*Az_cal;
        double Ay_rot = R[3]*Ax_cal + R[4]*Ay_cal + R[5]*Az_cal;
        double Az_rot = R[6]*Ax_cal + R[7]*Ay_cal + R[8]*Az_cal;

        // ---- Update accelerometer values ----
        Ax_cal = Ax_rot;
        Ay_cal = Ay_rot;
        Az_cal = Az_rot;

        // ---- Write updated values to CSV (same format) ----
        fprintf(fout, "%s,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",
                time_str, Ax_cal, Ay_cal, Az_cal, Gx, Gy, Gz);
    }

    fclose(fin);
    fclose(fout);
    printf("Calibrated and rotated data written to %s\n", output_csv);
    return 0;
}
