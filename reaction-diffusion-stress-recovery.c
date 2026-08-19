#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Physical Constants & Parameters from Paper
#define N0 5.0e12      // Trap precursor density (cm^-2) - Standard GAA NBTI value
#define KF_STRESS 1.0e-2 // Forward bond-breaking rate during stress (s^-1)
#define KF_RECOV 0.0    // Forward rate becomes zero during recovery
#define KR 1.0e-16     // Reverse annealing rate (cm^3/s for a=1, cm^1.5/s for a=2)
#define D_STRESS 1.0e-13 // Diffusion coefficient during stress (cm^2/s)

// Simulation Discretization
#define M 200          // Number of spatial grid points
#define DX 1.0e-7      // Spatial step size (1.0 nm = 10^-7 cm)
#define DT 1.0e-2      // Time step (s) - chosen to safely satisfy Fick's stability criterion (DT < DX^2 / 2D)
#define T_STRESS 1000.0 // Stress duration (s)
#define T_RECOV 1000.0  // Recovery duration (s)

// Solves the 1-D Reaction-Diffusion model for stress and recovery
// a = 1: Atomic Hydrogen (H)
// a = 2: Molecular Hydrogen (H2)
void run_rd_simulation(int a, int use_lockin, const char *filename) {
    double t_max = T_STRESS + T_RECOV;
    printf("Running C simulation for a = %d (Lock-in: %s) -> Output: %s...\n", 
           a, use_lockin ? "ON" : "OFF", filename);

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening output file");
        return;
    }

    // Write CSV Header
    fprintf(fp, "time,Nit,NH_interface,phase\n");

    // State Variables
    double *NH = (double *)calloc(M, sizeof(double));     // Hydrogen concentration array (cm^-3)
    double *NH_new = (double *)calloc(M, sizeof(double)); // Temporary array for next step
    double Nit = 0.0;                                     // Interface Trap Density (cm^-2)

    double t = 0.0;
    long long step = 0;
    long long save_interval = 100; // Periodically save data points

    while (t < t_max) {
        // Determine phase: Stress (0) vs Recovery (1)
        int is_recovery = (t >= T_STRESS);
        double kf = is_recovery ? KF_RECOV : KF_STRESS;
        
        // Handle Diffusivity Reduction (Lock-in effect shown in Eq (10) of Part 1 paper)
        double D = D_STRESS;
        if (is_recovery && use_lockin) {
            double S = 3.0; // Diffusivity reduction factor parameter from paper
            double t_rel = (t - T_STRESS);
            D = D_STRESS / (1.0 + (S * (t_rel / T_STRESS)));
        }

        double alpha = D * DT / (DX * DX);

        // 1. Calculate rate of change of interface traps (dNit/dt)
        // dNit/dt = kf * (N0 - Nit) - kr * Nit * (NH[0])^(1/a)
        double NH_interface = NH[0];
        if (NH_interface < 0.0) NH_interface = 0.0; // Physical boundary protection
        
        double dNit_dt = kf * (N0 - Nit) - KR * Nit * pow(NH_interface, 1.0 / a);
        
        // 2. Update Nit using Explicit Euler
        double Nit_new = Nit + dNit_dt * DT;
        if (Nit_new < 0.0) Nit_new = 0.0;
        if (Nit_new > N0) Nit_new = N0;

        // 3. Update spatial Hydrogen concentrations (Fick's Second Law)
        // Bulk diffusion (interior nodes 1 to M-2)
        for (int i = 1; i < M - 1; i++) {
            NH_new[i] = NH[i] + alpha * (NH[i+1] - 2.0 * NH[i] + NH[i-1]);
        }

        // Left Boundary Condition (x=0, interface):
        // Coupled to dNit/dt via flux: -D * dNH/dx = (1/a) * dNit/dt
        // Discretized using a ghost cell and central differences:
        // NH_new[0] = NH[0] + 2 * alpha * (NH[1] - NH[0]) + (2 * DT / (a * DX)) * dNit_dt
        NH_new[0] = NH[0] + 2.0 * alpha * (NH[1] - NH[0]) + (2.0 * DT / (a * DX)) * dNit_dt;

        // Right Boundary Condition (x=L, bulk oxide / gate interface):
        // Treated as an infinite sink (NH[M-1] = 0)
        NH_new[M-1] = 0.0;

        // Copy new state to current state
        for (int i = 0; i < M; i++) {
            NH[i] = NH_new[i];
        }
        Nit = Nit_new;

        // Save data points periodically
        if (step % save_interval == 0 || step < 1000) {
            fprintf(fp, "%e,%e,%e,%d\n", t, Nit, NH[0], is_recovery ? 1 : 0);
        }

        t += DT;
        step++;
    }

    fclose(fp);
    free(NH);
    free(NH_new);
    printf("Simulation completed successfully for a = %d (Lock-in: %s)!\n", 
           a, use_lockin ? "ON" : "OFF");
}

int main() {
    // Run atomic case (a=1) without and with lock-in
    run_rd_simulation(1, 0, "rd_stress_rec_atomic_classic.csv");
    run_rd_simulation(1, 1, "rd_stress_rec_atomic_lockin.csv");

    // Run molecular case (a=2) without and with lock-in
    run_rd_simulation(2, 0, "rd_stress_rec_molecular_classic.csv");
    run_rd_simulation(2, 1, "rd_stress_rec_molecular_lockin.csv");

    return 0;
}
