#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h> 

#define NX_IL 30                  
#define NX_HK 200   //Should also check the effect this number             
#define N_VARS (NX_IL + NX_HK + 2) 
#define DX_IL (1.0e-7 / NX_IL) 
// Values suggested
#define DX_HK (1.5e-7 / NX_HK)
#define DELTA 1.5e-8              
#define HALF_DELTA (DELTA / 2.0)
#define KB 8.617333262e-5

pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    double KF1_STRESS; double KR1_VAL;
    double KF2_VAL; double KR2_VAL;
    double DH_IL; double DH2_HK;
    double N01; double N02;
    double t_stress_total;
} PhysicsParams;

typedef struct {
    int thread_id;
    double T_celsius;
    double E_il_MV;
    double VG; // Added to pass voltage to the thread for filenames
} ThreadData;

double arrhenius(double prefactor, double Ea, double T) {
    return prefactor * exp(-Ea / (KB * T));
}

// CACHE-OPTIMIZED DENSE SOLVER (Chekc the mathematics)
int solve_dense(int n, double *A, double *b, double *x) {
    for (int i = 0; i < n; i++) x[i] = b[i];
    
    for (int i = 0; i < n; i++) {
        int max_row = i;
        double max_val = fabs(A[i*n + i]);
        
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k*n + i]) > max_val) {
                max_val = fabs(A[k*n + i]);
                max_row = k;
            }
        }
        if (max_val < 1e-20) return 0; 
        
        if (max_row != i) {
            double *Ai = A + i*n;
            double *Amax = A + max_row*n;
            for (int j = i; j < n; j++) {
                double tmp = Ai[j]; 
                Ai[j] = Amax[j]; 
                Amax[j] = tmp;
            }
            double tmp = x[i]; x[i] = x[max_row]; x[max_row] = tmp;
        }
        
        double *Ai = A + i*n;
        for (int k = i + 1; k < n; k++) {
            double *Ak = A + k*n;
            double factor = Ak[i] / Ai[i];
            
            for (int j = i; j < n; j++) Ak[j] -= factor * Ai[j];
            x[k] -= factor * x[i];
        }
    }
    
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        double *Ai = A + i*n;
        for (int j = i + 1; j < n; j++) sum += Ai[j] * x[j];
        x[i] = (x[i] - sum) / Ai[i];
    }
    return 1;
}

void compute_residuals(double *X, double *F, double *X_prev, double dt, double t, int is_recovery, PhysicsParams *p) {
    double *H = X;
    double *H2 = &X[NX_IL];
    double Nit1 = X[N_VARS - 2];
    double Nit2 = X[N_VARS - 1];

    // Both forward reactions shut off during recovery
    double KF1 = is_recovery ? 0.0 : p->KF1_STRESS;
    double KF2 = is_recovery ? 0.0 : p->KF2_VAL;
    
    double dNit1_dt = KF1 * (p->N01 - Nit1) - p->KR1_VAL * Nit1 * H[0];
    double dNit2_dt = KF2 * (p->N02 - Nit2) * H[NX_IL-1] - p->KR2_VAL * Nit2 * H2[0];

    F[N_VARS - 2] = (Nit1 - X_prev[N_VARS - 2]) / dt - dNit1_dt;
    F[N_VARS - 1] = (Nit2 - X_prev[N_VARS - 1]) / dt - dNit2_dt;

    double DH2_eff = is_recovery ? p->DH2_HK / (1.0 + 7.0 * (t / p->t_stress_total)) : p->DH2_HK;

    // Interface 1 (Silicon/IL) Boundary Condition
    F[0] = (H[0] - X_prev[0]) / dt - (p->DH_IL * (H[1] - H[0]) / DX_IL) / HALF_DELTA - dNit1_dt / HALF_DELTA;
    
    for(int i=1; i<NX_IL-1; i++) {
        F[i] = (H[i] - X_prev[i]) / dt - p->DH_IL * (H[i+1] - 2*H[i] + H[i-1]) / (DX_IL*DX_IL);
    }
    
    F[NX_IL-1] = (H[NX_IL-1] - X_prev[NX_IL-1]) / dt + (p->DH_IL * (H[NX_IL-1] - H[NX_IL-2]) / DX_IL + dNit2_dt) / HALF_DELTA;

    // Interface 2 (IL/High-K) Boundary Condition
    F[NX_IL] = (H2[0] - X_prev[NX_IL]) / dt - (DH2_eff * (H2[1] - H2[0]) / DX_HK + dNit2_dt) / HALF_DELTA;
    
    for(int i=1; i<NX_HK-1; i++) {
        F[NX_IL + i] = (H2[i] - X_prev[NX_IL + i]) / dt - DH2_eff * (H2[i+1] - 2*H2[i] + H2[i-1]) / (DX_HK*DX_HK);
    }
    
    // METAL GATE BOUNDARY CONDITION TOGGLE Between these two (Ask them what is happening there at metal gate)

    // Absorbing Sink (Vents H2 out. Required to prevent instant choking at 1.5nm thickness)
    // F[NX_IL + NX_HK - 1] = H2[NX_HK-1] - 0.0;

    //Reflective Wall (Traps H2. Only use this if you increase DX_HK to 1.5 microns)
    F[NX_IL + NX_HK - 1] = (H2[NX_HK-1] - X_prev[NX_IL + NX_HK - 1]) / dt - DH2_eff * (2.0 * H2[NX_HK-2] - 2.0 * H2[NX_HK-1]) / (DX_HK * DX_HK);
}

void compute_jacobian(double *X, double *F_base, double *J, double *X_prev, double dt, double t, int is_recovery, PhysicsParams *p, double *F_temp) {
    for (int j = 0; j < N_VARS; j++) {
        double epsilon = 1e-6 * fabs(X[j]) + 1.0; 
        X[j] += epsilon;
        
        compute_residuals(X, F_temp, X_prev, dt, t, is_recovery, p);
        
        for (int i = 0; i < N_VARS; i++) {
            J[i * N_VARS + j] = (F_temp[i] - F_base[i]) / epsilon;
        }
        X[j] -= epsilon; 
    }
}


// It takes some time compared to previous one so parallel computing is implimented
// THREAD PAYLOAD
void* run_simulation_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    double T_celsius = data->T_celsius;
    double VG = data->VG;
    double T_kelvin = T_celsius + 273.15;
    
    pthread_mutex_lock(&print_mutex);
    printf("Thread %d starting simulation for VG = %.1f V at %.0f C...\n", data->thread_id, VG, T_celsius);
    pthread_mutex_unlock(&print_mutex);

    PhysicsParams p;
    p.N01 = 5.0e12;
    p.N02 = 5.0e13;
    // STRESS TIME: Set to 10^8 seconds
    p.t_stress_total = 1.0e4;
    
    p.KF1_STRESS = arrhenius(0.22, 0.40, T_kelvin) * exp(data->E_il_MV * (0.38 + 1.2 * 0.01 / (KB * T_kelvin)));
    p.KR1_VAL = arrhenius(5.0e-6, 0.12, T_kelvin);
    p.KF2_VAL = arrhenius(5.75e3, 0.235, T_kelvin);
    p.KR2_VAL = arrhenius(7.5e-4, 0.20, T_kelvin);
    p.DH_IL = arrhenius(2.0e-2, 0.20, T_kelvin);
    p.DH2_HK = arrhenius(9.5e-8, 0.50, T_kelvin);

    double *X = (double*)calloc(N_VARS, sizeof(double));
    double *X_prev = (double*)calloc(N_VARS, sizeof(double));
    double *F = (double*)calloc(N_VARS, sizeof(double));
    double *J = (double*)calloc(N_VARS * N_VARS, sizeof(double));
    double *deltaX = (double*)calloc(N_VARS, sizeof(double));
    double *F_temp = (double*)calloc(N_VARS, sizeof(double));

    double *H = X;
    double *H2 = &X[NX_IL];

    double t = 0.0, dt = 1.0e-9;
    double t_recov = 1.0e4;
    
    char filename[256];
    // Changed filename to reflect voltage instead of temperature
    sprintf(filename, "nr_results_%.1fV.csv", fabs(VG));
    FILE *fp = fopen(filename, "w");
    
    if (fp != NULL) {
        fprintf(fp, "Phase,Time,Nit1,Nit2,Atomic_H0,Mol_H2,H2_Remaining,H2_Lost\n");

        for (int phase = 0; phase < 2; phase++) {
            int is_recovery = phase;
            t = 0.0; dt = 1.0e-9;
            double t_max = is_recovery ? t_recov : p.t_stress_total;
            
            while (t < t_max) {
                if (t + dt > t_max) dt = t_max - t;

                for (int i = 0; i < N_VARS; i++) X_prev[i] = X[i];

                for (int iter = 0; iter < 15; iter++) { 
                    compute_residuals(X, F, X_prev, dt, t, is_recovery, &p);
                    double max_res = 0.0;
                    for (int i=0; i<N_VARS; i++) if (fabs(F[i]) > max_res) max_res = fabs(F[i]);
                    if (max_res < 1e-4) break; 

                    compute_jacobian(X, F, J, X_prev, dt, t, is_recovery, &p, F_temp);
                    for (int i=0; i<N_VARS; i++) F[i] = -F[i];
                    if (!solve_dense(N_VARS, J, F, deltaX)) break; 

                    for (int i=0; i<N_VARS; i++) {
                        X[i] += deltaX[i];
                        if (X[i] < 0.0) X[i] = 0.0; 
                    }
                }
                t += dt;
                
                double H2_remaining = 0.0;
                for (int i = 0; i < NX_HK; i++) {
                    H2_remaining += H2[i] * DX_HK;
                }
                double Nit2_current = X[N_VARS-1];
                double H2_lost = Nit2_current - H2_remaining;
                
                fprintf(fp, "%s,%e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n", 
                        is_recovery ? "Recovery" : "Stress", t, X[N_VARS-2], Nit2_current, H[0], H2[0], H2_remaining, H2_lost);
                
                // this is important for long time duration
                dt *= 1.05;
                if (dt > 1.0e5) dt = 1.0e5;
            }
        }
        fclose(fp);
    }
    
    free(X); free(X_prev); free(F); free(J); free(deltaX); free(F_temp);
    
    pthread_mutex_lock(&print_mutex);
    printf("Thread %d finished VG = %.1f V. Saved to %s\n", data->thread_id, VG, filename);
    pthread_mutex_unlock(&print_mutex);
    
    pthread_exit(NULL);
    return NULL;
}


// MASTER RUN FUNCTION

void run_parallel_sweep() {
    // Replaced temperatures with voltages for the sweep
    double voltages_V[] = {-1.2, -1.5, -1.8, -2.0};
    double fixed_temp_C = 125.0;
    int num_sims = 4;

    // metal gate flatband voltage and strong inversion surface potential
    double V_fb = -0.1; 
    double psi_s = -0.8;                   

    printf("Booting parallel simulations for Voltage Sweep at %.0f C...\n", fixed_temp_C);

    pthread_t threads[4];
    ThreadData thread_data[4];

    for (int i = 0; i < num_sims; i++) {
        double VG = voltages_V[i];
        
        // Calculate absolute voltage drop across the oxide dynamically for each VG
        double V_ox = fabs(VG - V_fb - psi_s);    
        
        // Calculate field magnitude based on current VG
        double E_il_MV = (V_ox / (1.0e-7 + 1.5e-7 * (4.0 / 20.0))) / 1.0e6;      

        thread_data[i].thread_id = i;
        thread_data[i].T_celsius = fixed_temp_C;
        thread_data[i].E_il_MV = E_il_MV;
        thread_data[i].VG = VG;
        
        pthread_create(&threads[i], NULL, run_simulation_thread, (void*)&thread_data[i]);
    }

    for (int i = 0; i < num_sims; i++) {
        pthread_join(threads[i], NULL);
    }
    
    printf("All parallel simulations completed successfully!\n");
}


run_parallel_sweep();
