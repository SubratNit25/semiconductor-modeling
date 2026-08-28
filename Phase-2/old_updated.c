#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Added the missing Thomas Algorithm solver so the code compiles and runs!
int solve_tridiagonal(int n, double *a, double *b, double *c, double *rhs, double *x) {
    double *c_prime = (double *)malloc(n * sizeof(double));
    double *d_prime = (double *)malloc(n * sizeof(double));
    if (!c_prime || !d_prime) return 0;
    
    c_prime[0] = c[0] / b[0];
    d_prime[0] = rhs[0] / b[0];
    for (int i = 1; i < n; i++) {
        double m = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (rhs[i] - a[i] * d_prime[i - 1]) / m;
    }
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--) x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    
    free(c_prime); free(d_prime);
    return 1;
}

// We will do the actual process...

#define NX_IL 300
#define NX_HK 1000

int main()
{
    // define constants
    const double kB = 8.617333262e-5;
    const double T = 125.0 + 273.15; // check for different temperature
    const double kT = kB * T;

    // Check what are the exact values? (assumptions but is it valid?)
    const double L_IL = 1.0e-7;
    const double L_HK = 1.5e-7; //(AS per suggetion)
    const double delta = 1.5e-8;

    const double dx_IL = L_IL / NX_IL;
    const double dx_HK = L_HK / NX_HK;

    const double N01 = 5.0e12; // directly used from book
    const double N02 = 5.0e13;

    const double KR1 = 5.0e-6 * exp(-0.12 / kT);
    const double KF2 = 5.75e3 * exp(-0.235 / kT);
    const double KR2 = 7.5e-4 * exp(-0.20 / kT);

    // Two species here, so each has its own D0 and Ea.
    // need to note that these are NOT the final D values.
    // Temperature comes through the Arrhenius factor.
    const double DH_IL = 2.0e-2 * exp(-0.20 / kT);
    const double DH2_HK = 9.5e-8 * exp(-0.50 / kT);

    // --- PROFF UPDATE: Correct NBTI Electrostatics ---
    // Replaced hardcoded EOX = 7.0 with actual calculation based on VG
    const double VG = -1.5; 
    const double V_fb = -1.0;
    const double psi_s = -0.8;
    const double V_ox = fabs(VG - V_fb - psi_s); 
    const double EOX = (V_ox / (1.0e-7 + 1.5e-7 * (4.0 / 20.0))) / 1.0e6;
    // -------------------------------------------------

    const double Gamma0 = 0.38; //check it
    const double alpha = 1.2;

    // For planar MOSFET. Need to check why A = 7.
    const double S_lock = 7.0;

    // Calculate the field/temperature dependent part first.
    double Gamma_E = Gamma0 +(alpha * 0.01) / kT;

    // KF1 is different from the other rates because it also depends
    // on the stress electric field and D1 process parameters. 
    double KF1_stress = 0.22 * exp(EOX * Gamma_E) * exp(-0.40 / kT);

    // H  -> current atomic-H concentration profile
    // H_prev -> concentration profile from previous timestep
    double H[NX_IL] = {0.0};
    double H_prev[NX_IL] = {0.0};

    // Same idea for molecular H2 in the High-K part.
    double H2[NX_HK] = {0.0};
    double H2_prev[NX_HK] = {0.0};

    // Now we create coefficient arrays for the tridiagonal systems.
    // These will be filled before calling the solver.
    double a_H[NX_IL];
    double b_H[NX_IL];
    double c_H[NX_IL];
    double rhs_H[NX_IL];

    double a_H2[NX_HK];
    double b_H2[NX_HK];
    double c_H2[NX_HK];
    double rhs_H2[NX_HK];

    // initial generated trap density.
    // this is initial condition assumption for the stress simulation.
    // N01 and N02 are still finite precursor populations. might be associtaed with some parameters check it.
    double NIT1 = 0.0;
    double NIT2 = 0.0;

    double t = 0.0;       // current simulation time
    double dt = 1.0e-9;   // start small; later we increase dt

    const double t_stress = 100000000.0;
    const double t_recov = 100000000.0;

    // Only for printing results at roughly logarithmic time intervals.
    double print_target = 1.0e-4;

    //Print some values
    FILE *fp = fopen("rd_results.csv", "w");

    if (fp == NULL) {
        fprintf(stderr, "ERROR: Could not open rd_results.csv\n");
        return EXIT_FAILURE;
    }

    // --- PROFF UPDATE: Added H and H2 to the prints ---
    fprintf(fp, "Phase,Time,NIT1,NIT2,Total_NIT,Atomic_H_Int,Mol_H2_Gate\n");
    printf("=== PHASE 1: STRESS " "(V_G = -1.5V, T = 125.0 C) ===\n" );
    printf( "%-12s %-16s %-16s %-16s %-16s %-16s\n", "Time(s)","NIT_1(/cm2)", "NIT_2(/cm2)", "Total_NIT", "H[0]", "H2_Gate");

    // Lets START STRESS
    while (t < t_stress) {

        // have to make sure the last step ends exactly at t_stress.
        if (t + dt > t_stress)
            dt = t_stress - t;
        // save the old concentration profiles before calculating the new timestep.
        for (int i = 0; i < NX_IL; i++)
            H_prev[i] = H[i];
        for (int i = 0; i < NX_HK; i++)
            H2_prev[i] = H2[i];

        // also keep the previous trap densities for this timestep.
        double NIT1_prev = NIT1;
        double NIT2_prev = NIT2;

        // a few iterations are used because H, H2 and NIT are coupled.
        for (int iter = 0; iter < 400; iter++) {
            // first calculate the reaction terms at both interfaces.
            double Denom1 =1.0 + dt * KF1_stress + dt * KR1 * H[0];
            double Denom2 = 1.0 + dt * KF2 * H[NX_IL - 1] + dt * KR2 * H2[0];

            // Diffusion coefficient for the dimensionless
            // finite-difference parameter r_H.
            double r_H = DH_IL * dt / (dx_IL * dx_IL);
            double r_H_delta = DH_IL * dt / ((delta / 2.0) * dx_IL);

            // now build the tridiagonal system for atomic H.
            for (int i = 0; i < NX_IL; i++) {
                if (i == 0) {
                    // Interface 1 boundary.
                    a_H[0] = 0.0;
                    b_H[0] = 1.0 + r_H_delta + (dt / (delta / 2.0)) * (KR1 * NIT1_prev) / Denom1;
                    c_H[0] = -r_H_delta;
                    rhs_H[0] = H_prev[0] +(dt / (delta / 2.0)) *(KF1_stress * (N01 - NIT1_prev)) /Denom1;
                } else if (i == NX_IL - 1) {
                    // Interface 2 boundary.
                    a_H[i] = -r_H_delta;
                    b_H[i] = 1.0 + r_H_delta + (dt / (delta / 2.0)) * (KF2 * (N02 - NIT2_prev)) / Denom2;
                    c_H[i] = 0.0;
                    rhs_H[i] =H_prev[i] +(dt / (delta / 2.0)) *(KR2 * NIT2_prev * H2[0]) / Denom2;
                } else {
                    // Normal diffusion node.
                    a_H[i] = -r_H;
                    b_H[i] = 1.0 + 2.0 * r_H;
                    c_H[i] = -r_H;
                    rhs_H[i] = H_prev[i];
                }
            }
            // Now Solve for the new atomic-H concentration profile.
            if (!solve_tridiagonal(NX_IL, a_H, b_H, c_H, rhs_H,H)) {
                fprintf(stderr,"ERROR: H solver failed during stress.\n");
                fclose(fp);
                return EXIT_FAILURE;
            }
            // Now lets do the same thing for H2.
            double r_H2 = DH2_HK * dt / (dx_HK * dx_HK);
            double r_H2_delta = DH2_HK * dt / ((delta / 2.0) * dx_HK);

            // Build the H2 tridiagonal system.
            for (int i = 0; i < NX_HK; i++) {
                if (i == 0) {
                    // H2 starts from the second-interface side.
                    a_H2[0] = 0.0;
                    b_H2[0] = 1.0 + r_H2_delta +(dt / (delta / 2.0)) *(KR2 * NIT2_prev) / Denom2;
                    c_H2[0] = -r_H2_delta;
                    rhs_H2[0] = H2_prev[0] + (dt / (delta / 2.0)) * (KF2 * (N02 - NIT2_prev) * H[NX_IL - 1]) / Denom2;
                } else if (i == NX_HK - 1) {
                    // Far end of the H2 diffusion domain.
                    
                    // --- PROFF UPDATE: Activated Reflective Boundary ---
                    a_H2[i] = 0.0;
                    b_H2[i] = 1.0;
                    c_H2[i] = 0.0;
                    rhs_H2[i] = 0.0;
                    //reflective
                    // a_H2[i] = -2.0 * r_H2;
                    // b_H2[i] = 1.0 + 2.0 * r_H2;
                    // c_H2[i] = 0.0;
                    // rhs_H2[i] = H2_prev[i];
                    // ---------------------------------------------------
                    
                } else {
                    // Normal H2 diffusion node.
                    a_H2[i] = -r_H2;
                    b_H2[i] = 1.0 + 2.0 * r_H2;
                    c_H2[i] = -r_H2;
                    rhs_H2[i] = H2_prev[i];
                }
            }

            // Solve for the new H2 concentration profile.
            if (!solve_tridiagonal( NX_HK, a_H2, b_H2, c_H2, rhs_H2, H2)) {
                fprintf(stderr,"ERROR: H2 solver failed during stress.\n");
                fclose(fp);
                return EXIT_FAILURE;
            }

            // Now update the two interface trap densities.
            NIT1 =NIT1_prev + dt * (KF1_stress * (N01 - NIT1_prev) - KR1 * NIT1_prev * H[0]) /Denom1;
            NIT2 = NIT2_prev + dt * (KF2 *(N02 - NIT2_prev) * H[NX_IL - 1] - KR2 * NIT2_prev * H2[0]) / Denom2;

            // trap densities should be physically bounded.
            if (NIT1 < 0.0) NIT1 = 0.0;
            if (NIT1 > N01) NIT1 = N01;
            if (NIT2 < 0.0) NIT2 = 0.0;
            if (NIT2 > N02) NIT2 = N02;
        }

        // Move to the next stress time.
        t += dt;

        // Increase timestep slowly after the early-time region.
        // if (dt < 2.0)
        //     dt *= 1.05;
        dt *= 1.05;
        if (dt > 1.0e5) 
            dt = 1.0e5;
        // print only selected times so the output stays manageable.
        if (t >= print_target || t >= t_stress) {
            printf("%-12.4e %-16.4e %-16.4e %-16.4e %-16.4e %-16.4e\n", t, NIT1, NIT2, NIT1 + NIT2, H[0], H2[NX_HK-1]);
            fprintf(fp,"Stress,%e,%e,%e,%e,%e,%e\n",t,NIT1,NIT2, NIT1 + NIT2, H[0], H2[NX_HK-1]);
            print_target *= 10.0;
        }
    }

    //  RECOVERY 
    printf("\n=== PHASE 2: RECOVERY " "(V_GREC = 0V, T = 125.0 C) ===\n");
    printf("%-12s %-16s %-16s %-16s %-16s %-16s\n", "t_rec(s)","NIT_1(/cm2)","NIT_2(/cm2)","Total_NIT", "H[0]", "H2_Gate");

    double t_rec = 0.0;
    dt = 1.0e-9;
    print_target = 1.0e-4;

    // START RECOVERY smilar way
    while (t_rec < t_recov) {
        if (t_rec + dt > t_recov)
            dt = t_recov - t_rec;

        // Again save the previous state before solving
        // the next recovery timestep.
        for (int i = 0; i < NX_IL; i++)
            H_prev[i] = H[i];
        for (int i = 0; i < NX_HK; i++)
            H2_prev[i] = H2[i];

        double NIT1_prev = NIT1;
        double NIT2_prev = NIT2;

        // Check how H2 diffusion changes during recovery.
        double DH2_HK_rec = DH2_HK /(1.0 +S_lock *(t_rec / t_stress));

        // H, H2 and trap densities are still coupled here.
        for (int iter = 0; iter < 400; iter++) {
            double Denom1 =1.0 + dt * KR1 * H[0];
            double Denom2 =1.0 + dt * KR2 * H2[0];
            double r_H =DH_IL * dt / (dx_IL * dx_IL);

            double r_H_delta = DH_IL * dt / ((delta / 2.0) * dx_IL);

            // Build H system for recovery.
            for (int i = 0; i < NX_IL; i++) {
                if (i == 0) {
                    // interface 1- now mainly passivation.
                    a_H[0] = 0.0;
                    b_H[0] =1.0 + r_H_delta + (dt / (delta / 2.0)) *(KR1 * NIT1_prev) / Denom1;
                    c_H[0] = -r_H_delta;
                    rhs_H[0] = H_prev[0];
                } else if (i == NX_IL - 1) {
                    // interface 2 during recovery.
                    a_H[i] = -r_H_delta;
                    b_H[i] = 1.0 + r_H_delta;
                    c_H[i] = 0.0;
                    rhs_H[i] =H_prev[i] +(dt / (delta / 2.0)) *(KR2 * NIT2_prev * H2[0]) /Denom2;
                    
                } else {
                    a_H[i] = -r_H;
                    b_H[i] = 1.0 + 2.0 * r_H;
                    c_H[i] = -r_H;
                    rhs_H[i] = H_prev[i];
                }
            }

            if (!solve_tridiagonal(NX_IL,a_H,b_H,c_H,rhs_H,H)) 
            {
                fprintf(stderr, "ERROR: H solver failed during recovery.\n");
                fclose(fp);
                return EXIT_FAILURE;
            }

            // H2 diffusion during recovery.
            double r_H2 =DH2_HK_rec * dt /(dx_HK * dx_HK);
            double r_H2_delta = DH2_HK_rec * dt / ((delta / 2.0) * dx_HK);

            for (int i = 0; i < NX_HK; i++) {
                if (i == 0) {
                    a_H2[0] = 0.0;
                    b_H2[0] = 1.0 + r_H2_delta + (dt / (delta / 2.0)) * (KR2 * NIT2_prev) / Denom2;
                    c_H2[0] = -r_H2_delta;
                    rhs_H2[0] = H2_prev[0];
                } else if (i == NX_HK - 1) {
                    // Far boundary of H2 domain.
                    
                    // --- PROFF UPDATE: Activated Reflective Boundary ---
                    a_H2[i] = 0.0;
                    b_H2[i] = 1.0;
                    c_H2[i] = 0.0;
                    rhs_H2[i] = 0.0;

                    // Reflective / zero-flux far boundary.
                    // a_H2[i] = -2.0 * r_H2;
                    // b_H2[i] = 1.0 + 2.0 * r_H2;
                    // c_H2[i] = 0.0;
                    // rhs_H2[i] = H2_prev[i];
                    // ---------------------------------------------------
                    
                } else {
                    a_H2[i] = -r_H2;
                    b_H2[i] = 1.0 + 2.0 * r_H2;
                    c_H2[i] = -r_H2;
                    rhs_H2[i] = H2_prev[i];
                }
            }

            if (!solve_tridiagonal(NX_HK,a_H2,b_H2,c_H2,rhs_H2,H2)) 
            {
                fprintf(stderr,"ERROR: H2 solver failed during recovery.\n");
                fclose(fp);
                return EXIT_FAILURE;
            }

            // During recovery the reverse reactions reduce NIT.
            NIT1 = NIT1_prev + dt * (-KR1 * NIT1_prev * H[0]) /Denom1;
            NIT2 =NIT2_prev + dt *(-KR2 * NIT2_prev * H2[0]) / Denom2;

            //check  trap density should not become negative.
            if (NIT1 < 0.0) NIT1 = 0.0;
            if (NIT2 < 0.0) NIT2 = 0.0;
        }

        t_rec += dt;

        // if (dt < 2.0)
        //     dt *= 1.05;
        dt *= 1.05;
        if (dt > 1.0e5) 
            dt = 1.0e5;

        if (t_rec >= print_target || t_rec >= t_recov) {
            printf("%-12.4e %-16.4e %-16.4e %-16.4e %-16.4e %-16.4e\n",t_rec, NIT1,NIT2,NIT1 + NIT2, H[0], H2[NX_HK-1]);
            fprintf(fp,"Recovery,%e,%e,%e,%e,%e,%e\n",t_rec,NIT1,NIT2,NIT1 + NIT2, H[0], H2[NX_HK-1]);
            print_target *= 10.0;
        }
    }
    fclose(fp);
    return EXIT_SUCCESS;
}