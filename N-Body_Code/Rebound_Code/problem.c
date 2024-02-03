/**
 * Planetary migration in the GJ876 system
 *
 * This example applies dissipative forces to two
 * bodies orbiting a central object. The forces are specified
 * in terms of damping timescales for the semi-major axis and
 * eccentricity. This mimics planetary migration in a protostellar disc. 
 * The example reproduces the study of Lee & Peale (2002) on the 
 * formation of the planetary system GJ876. For a comparison, 
 * see figure 4 in their paper. The IAS15 or WHFAST integrators
 * can be used. Note that the forces are velocity dependent.
 * Special thanks goes to Willy Kley for helping me to implement
 * the damping terms as actual forces. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "spline.h"

// initialise variables and simulation feature variables
double tmax, agnmass, r_g, r_s, massscale, lenscale, lenscale_m, lenscale_rs, nondim_rs, r_isco, timescale, velscale, stefboltz, c_v, c_nbody, m_p, sigma_T;
const double G_pc = 4.3e-3, M_odot = 1.98e30, c = 3e8, pc_to_m = 3.086e16, gamma_coeff = 5./3.;
char output_folder[200] = "./OUTPUT_", orbits_filename[200], positions_filename[200], mergers_filename[200], r_filename[200], sigma_filename[200], temp_filename[200], aratio_filename[200], opacity_filename[200];
const double spin_mult      = 1e10;         // I'm storing the BH spins in the particle radius parameter (cursed, I know), so divide it by this huge multiplier so that the radius of the particle doesn't actually have any effect (e.g. with collisions)
int num_BH                  = 0;            // integer value to track how many seed BHs we've had in total so far
int MIGRATION_PRESCRIPTION  = 0;            // 0 for Pardekooper (2010?) migration prescription, 1 for Jimenez and Masset (2017)
int THERMAL_TORQUES         = 0;            // set to 1 to include thermal torques from Grishin et al (2023)
int RAND_BH_MASSES          = 0;            // 0 for 10 solar mass seed BHs, 1 for randomly sampled masses
int MERGER_KICKS            = 0;            // 0 for no kicks in mergers, 1 for kicks in random direction
int MERGER_CRITERION        = 0;            // 0 for binding_energy < KE, 1 for criterion in Li et al
double MUTUAL_HILL_PROP     = 1.;           // proportion of mutual hill radius to consider a merger
double ACCRETION            = 0.;           // 0 for no accretion, any other number to simulate accretion at that *proportion* of the eddington limit
int ADD_BH_RAND_TIME        = 0.;           // 0 for adding BHs at regular intervals, 1 for adding them at exponential randomly distributed times
double ADD_BH_INTERVAL      = 1e5;          // if ADD_BH_RAND_TIME==0, this is the interval for adding. if ==1, this is the mean of the distribution
double NEXT_ADD_TIME        = 0.;           // variable to say when to add the next BH (used when randomly choosing the interval)
int BH_SPINS                = 0;            // 1 if simulating BH spins from NR fits, 0 if BHs are non-spinning
double RETRO_PROB           = 0.2;          // probability of a black hole being initialised on a retrograde orbit

// now define our data for our disk parameter splines. start by initialising some needed arrays and constants
double temp_deriv_coeffs[3] = {0, 0, 0}, sigma_deriv_coeffs[] = {0, 0, 0}, asp_deriv_coeffs[] = {0, 0, 0};
const int n_spline_data = 100;
// disk surface density spline data
double log_sigma_spline[100];
// disk temperature spline data
double log_temp_spline[100];
// disk aspect ratio spline data
double log_aratio_spline[100];
// disk opacity spline data
double log_opacity_spline[100];

double log_radii_data[100], log_sigma_data[100], log_temps_data[100], log_aratio_data[100], log_opacity_data[100];


double uniform(double lower, double upper){
    // with thanks to https://stackoverflow.com/questions/63981013/generating-random-double-between-1-and-100-in-c?noredirect=1&lq=1
    double random = ((double) rand() / RAND_MAX) * (upper - lower) + lower;
    return random;
}
double exponential_rv(struct reb_simulation* r, double mean){
    // see #Related_Distributions on https://en.wikipedia.org/wiki/Rayleigh_distribution
    double x = reb_random_rayleigh(r, 1./sqrt(2. * 1./mean));
    return x*x;
}

void create_filenames(int mass, double f_edd, double alpha){
    char file_name_format[40];
    sprintf(file_name_format, "M%d-f%.2lf-a%.3lf", mass, f_edd, alpha);
    char folder_suffix[20] = "/";
    strcat(output_folder, file_name_format); strcat(output_folder, folder_suffix);
    strcpy(orbits_filename, output_folder); strcat(orbits_filename, "orbits.txt");
    strcpy(positions_filename, output_folder); strcat(positions_filename, "positions.txt");
    strcpy(mergers_filename, output_folder); strcat(mergers_filename, "mergers.txt");
    
    // define names of disk property files
    strcpy(r_filename, output_folder); strcat(r_filename, "log_radii_");
    strcpy(sigma_filename, output_folder); strcat(sigma_filename, "Sigma_");
    strcpy(temp_filename, output_folder); strcat(temp_filename, "temps_");
    strcpy(aratio_filename, output_folder); strcat(aratio_filename, "h_");
    strcpy(opacity_filename, output_folder); strcat(opacity_filename, "kappa_");
    strcat(r_filename, file_name_format); strcat(r_filename, ".csv");
    strcat(sigma_filename, file_name_format); strcat(sigma_filename, ".csv");
    strcat(temp_filename, file_name_format); strcat(temp_filename, ".csv");
    strcat(aratio_filename, file_name_format); strcat(aratio_filename, ".csv");
    strcat(opacity_filename, file_name_format); strcat(opacity_filename, ".csv");
}
void populate_spline_arrays(){
    FILE *log_r_file = fopen(r_filename, "r"); // x data, units of log10 schwarzschild radii
    FILE *sigma_file = fopen(sigma_filename, "r");  // y data, cgs units on a log10 scale
    FILE *temp_file = fopen(temp_filename, "r");  // y data, cgs units on a log10 scale
    FILE *aratio_file = fopen(aratio_filename, "r");  // y data, cgs units on a log10 scale
    FILE *opacity_file = fopen(opacity_filename, "r"); // y data, cgs units on a log10 scale
    for (int i = 0; i < n_spline_data; i++){
        if (fscanf(log_r_file, "%le", &log_radii_data[i]));
        if (fscanf(sigma_file, "%le", &log_sigma_data[i]));
        if (fscanf(temp_file, "%le", &log_temps_data[i]));
        if (fscanf(aratio_file, "%le", &log_aratio_data[i]));
        if (fscanf(opacity_file, "%le", &log_opacity_data[i]));
    }
    fclose(log_r_file); fclose(sigma_file); fclose(temp_file); fclose(aratio_file); fclose(opacity_file); 
}
void eval_splines(){
    // this evaluates the splines at run time (since we can't evaluate them at compile time).
    spline(log_radii_data, log_sigma_data, n_spline_data, log_sigma_spline);
    spline(log_radii_data, log_temps_data, n_spline_data, log_temp_spline);
    spline(log_radii_data, log_aratio_data, n_spline_data, log_aratio_spline);
    spline(log_radii_data, log_opacity_data, n_spline_data, log_opacity_spline);
}
double disk_surfdens(double logr){
    double sigma = pow(10., splint(log_radii_data, log_sigma_data, log_sigma_spline, n_spline_data, logr) + 1.); // +1 in power to go from cgs to SI
    return sigma * lenscale_m * lenscale_m / massscale;
}
double disk_temp(double logr){
    return pow(10., splint(log_radii_data, log_temps_data, log_temp_spline, n_spline_data, logr));
}
double disk_aspectratio(double logr){
    return pow(10., splint(log_radii_data, log_aratio_data, log_aratio_spline, n_spline_data, logr));
}
double disk_opacity(double logr){
    double kappa = pow(10., splint(log_radii_data, log_opacity_data, log_opacity_spline, n_spline_data, logr) - 1.); // -1 to go from cgs to SI
    return kappa * massscale / (lenscale_m * lenscale_m);
}
double disk_dens(double logr){
    return disk_surfdens(logr) / (2. * disk_aspectratio(logr) * pow(10., logr));
}
double disk_angvel(double logr){
    double omega = sqrt(G_pc * agnmass / pow((pow(10., logr) * 2. * r_g), 3)) * 1000;
    return omega * lenscale / velscale;
}
double disk_sigma_deriv(double logr){
    return splderiv(log_radii_data, log_sigma_data, log_sigma_spline, n_spline_data, logr);
}
double disk_temp_deriv(double logr){
    return splderiv(log_radii_data, log_temps_data, log_temp_spline, n_spline_data, logr);
}
double disk_aspratio_deriv(double logr){
    return splderiv(log_radii_data, log_aratio_data, log_aratio_spline, n_spline_data, logr);
}

double lognormal_mass(struct reb_simulation* r){
    // box-muller transform for normally distributed random numbers     https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    double mean = 1.21, sd = 0.2, lower = 0.3, upper = 1.65;
    double z0 = 0;
    while (lower > z0 || upper < z0){       // check to see if the mass is outside of our bounds
        double r1 = reb_random_uniform(r, 0.0, 1.);
        double r2 = reb_random_uniform(r, 0.0, 1.);
        z0 = sqrt(-2. * log(r1)) * sin(2. * M_PI * r2);
        z0 *= sd;
        z0 += mean;
    }
    return pow(10., z0);    // we were generating log-normal RVs, so transform back to normal
}


void check_mergers(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0];
    const int N = r->N;
    int break_check = 0;
    for (int i = 1; i < N - 1; i++){
        struct reb_particle* p1 = &(particles[i]); // get the particle
        double m1 = p1->m;
        const double dx1 = p1->x-com.x, dy1 = p1->y-com.y, dz1 = p1->z-com.z;
        const double r1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        for (int j = i + 1; j < N; j++){
            struct reb_particle* p2 = &(particles[j]); // get the particle
            double m2 = p2->m;
            const double dx2 = p2->x-com.x, dy2 = p2->y-com.y, dz2 = p2->z-com.z;
            const double r2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
            double R_mH = pow((m1 + m2) / (3 * com.m), 1./3.) * (r1 + r2) / 2.;   // mutual hill radius
            double dist = sqrt((dx1-dx2)*(dx1-dx2) + (dy1-dy2)*(dy1-dy2) + (dz1-dz2)*(dz1-dz2));
            if (dist < (MUTUAL_HILL_PROP * R_mH)){
                // now to check if their relative kinetic energy is low enough for capture
                const double dvx1 = p1->vx-com.vx, dvy1 = p1->vy-com.vy, dvz1 = p1->vz-com.vz;
                const double dvx2 = p2->vx-com.vx, dvy2 = p2->vy-com.vy, dvz2 = p2->vz-com.vz;
                double reduced_mass = 1. / (1. / m1 + 1. / m2);
                double rel_kin_energy = 0.5 * reduced_mass * ((dvx1-dvx2)*(dvx1-dvx2) + (dvy1-dvy2)*(dvy1-dvy2) + (dvz1-dvz2)*(dvz1-dvz2));
                double binding_energy = m1 * m2 / (2. * R_mH);
                bool merger = false;
                if (MERGER_CRITERION == 0){     // classic binding energy < relative kinetic energy
                    merger = binding_energy < rel_kin_energy;
                } else if (MERGER_CRITERION == 1){      // criterion based on the results of Li et al 2023
                    double logr = log10(r1 / nondim_rs), Sigma = disk_surfdens(logr);
                    double cond = (Sigma * r1*r2) / (m1 + m2);
                    double abs_Ehill = (m1 + m2) / (2. * R_mH);
                    double E_b = rel_kin_energy / reduced_mass - (m1 + m2) / (0.3 * R_mH);
                    double mu_crit = 19.1 * E_b / abs_Ehill + 25.6;
                    if (1.3 > mu_crit){
                        mu_crit = 1.3;
                    }
                    merger = cond > mu_crit;
                }
                if (merger){
                    puts("\nMerger!");
                    // set p1 as the merger remnant, with velocity and position as the mass-weighted average of the two merged particles
                    p1->x = (m1 * p1->x + m2 * p2->x) / (m1 + m2);
                    p1->y = (m1 * p1->y + m2 * p2->y) / (m1 + m2);
                    p1->z = (m1 * p1->z + m2 * p2->z) / (m1 + m2);
                    p1->vx = (m1 * p1->vx + m2 * p2->vx) / (m1 + m2);
                    p1->vy = (m1 * p1->vy + m2 * p2->vy) / (m1 + m2);
                    p1->vz = (m1 * p1->vz + m2 * p2->vz) / (m1 + m2);
                    int gen1 = p1->generation, gen2 = p2->generation;
                    p1->generation = fmax(gen1, gen2) + 1;
                    double kick_vel = 0., spin_eff = 0., nondim_spin = 0.;
                    if (MERGER_KICKS == 0){
                        double q = fmin(m1 / m2, m2 / m1), opq = 1. + q, nu = q / (opq*opq);
                        double eps_rad = nu * (1. - 4. * nu) * (1. - 2.*sqrt(2.)/3.) + 0.048 * (4.*nu)*(4.*nu);
                        p1->m = (m1 + m2) * (1. - eps_rad);
                    } else if (MERGER_KICKS == 1 && BH_SPINS == 0){     // include a kick in a random direction
                        double q = fmin(m1 / m2, m2 / m1), opq = 1. + q, nu = q / (opq*opq);
                        double A = 1.2 * 1e7 / velscale, B = -0.93;
                        double vel_kick = A * q*q*(1. - q) / (opq*opq*opq*opq*opq) * (1. + B * q / (opq*opq));
                        double angle = reb_random_uniform(r, 0.0, 2.*M_PI);
                        double xprop = sin(angle);
                        double yprop = -cos(angle);
                        double mult = sqrt(vel_kick*vel_kick / (xprop*xprop + yprop*yprop));
                        xprop *= mult; yprop *= mult;
                        p1->vx += xprop; 
                        p1->vy += yprop; 
                        kick_vel = vel_kick * velscale / 1000.;
                        double eps_rad = nu * (1. - 4. * nu) * (1. - 2.*sqrt(2.)/3.) + 0.048 * (4.*nu)*(4.*nu);
                        p1->m = (m1 + m2) * (1. - eps_rad);
                    } else if (MERGER_KICKS == 1 && BH_SPINS == 1){ // include a kick in a semi-analytic direction
                        // physics here modelled from chapter 14.3 of Maggiore, M. 2018, Gravitational Waves: Volume 2: Astrophysics and Cosmology, Gravitational Waves (Oxford University Press)
                        double q = fmin(m1 / m2, m2 / m1), opq = 1. + q, nu = q / (opq*opq);
                        double sx1 = p1->spin_x, sx2 = p2->spin_x, sy1 = p1->spin_y, sy2 = p2->spin_y, sz1 = p1->spin_z, sz2 = p2->spin_z;   
                        const double spin1_modulus = sqrt(sx1*sx1 + sy1*sy1 + sz1*sz1), spin2_modulus = sqrt(sx2*sx2 + sy2*sy2 + sz2*sz2);
                        double usx1 = 0., usy1 = 0., usz1 = 0., usx2 = 0., usy2 = 0., usz2 = 0.;
                        if (spin1_modulus != 0){
                            usx1 = sx1 / spin1_modulus, usy1 = sy1 / spin1_modulus, usz1 = sz1 / spin1_modulus;
                        }
                        if (spin2_modulus != 0){
                            usx2 = sx2 / spin2_modulus, usy2 = sy2 / spin2_modulus, usz2 = sz2 / spin2_modulus; 
                        }
                        double aa1 = fmax(spin1_modulus, spin2_modulus), aa2 = fmin(spin1_modulus, spin2_modulus);
                        const double Lx = (dy1 - (dy1 + dy2)/2.) * (dvz1 - (dvz1 + dvz2)/2.) - (dz1 - (dz1 + dz2)/2.) * (dvy1 - (dvy1 + dvy2)/2.);
                        const double Ly = (dz1 - (dz1 + dz2)/2.) * (dvx1 - (dvx1 + dvx2)/2.) - (dx1 - (dx1 + dx2)/2.) * (dvz1 - (dvz1 + dvz2)/2.);
                        const double Lz = (dx1 - (dx1 + dx2)/2.) * (dvy1 - (dvy1 + dvy2)/2.) - (dy1 - (dy1 + dy2)/2.) * (dvx1 - (dvx1 + dvx2)/2.);
                        const double ang_mom_modulus = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
                        const double unit_Lx = Lx / ang_mom_modulus, unit_Ly = Ly / ang_mom_modulus, unit_Lz = Lz / ang_mom_modulus;
                        double cos_alpha = usx1*usx2 + usy1*usy2 + usz1*usz2;
                        double cos_beta = usx1*unit_Lx + usy1*unit_Ly + usz1*unit_Lz;
                        double cos_gamma = usx2*unit_Lx + usy2*unit_Ly + usz2*unit_Lz;
                        double term1 = -0.129 / ((1. + q*q)*(1. + q*q)) * (aa1*aa1 + aa2*aa2*q*q*q*q + 2.*aa1*aa2*q*q*cos_alpha);
                        double term2 = (-0.384 * nu - 2.686 + 2)/(1. + q*q) * (aa1*cos_beta + aa2*q*q*cos_gamma);
                        double term3 = 3.464 - 3.454*nu + 2.353*nu*nu;
                        double l_mag = term1 + term2 + term3;
                        if (m1 > m2){
                            p1->spin_x = nu/q * (sx1 + q*q * sx2) + nu * l_mag * unit_Lx;
                            p1->spin_y = nu/q * (sy1 + q*q * sy2) + nu * l_mag * unit_Ly;
                            p1->spin_z = nu/q * (sz1 + q*q * sz2) + nu * l_mag * unit_Lz;
                        } else {
                            p1->spin_x = nu/q * (sx2 + q*q * sx1) + nu * l_mag * unit_Lx;
                            p1->spin_y = nu/q * (sy2 + q*q * sy1) + nu * l_mag * unit_Ly;
                            p1->spin_z = nu/q * (sz2 + q*q * sz1) + nu * l_mag * unit_Lz;
                        }
                        nondim_spin = sqrt(p1->spin_x*p1->spin_x + p1->spin_y*p1->spin_y + p1->spin_z*p1->spin_z);

                        // now calculate the new mass of the final BH
                        // equations from Barausse et al 2012
                        double const0 = 0.04826, const1 = 0.01559;
                        double a_tilde = (spin1_modulus * cos_beta + q*q* spin2_modulus * cos_gamma) / (opq*opq);
                        double Z1 = 1. + pow(1. - a_tilde*a_tilde, 1./3.) * (pow(1. + a_tilde, 1./3.) + pow(1. - a_tilde, 1./3.));
                        double Z2 = sqrt(3. * a_tilde*a_tilde + Z1*Z1);
                        double r_eq_isco = 3. + Z2 - copysignl(1.0, a_tilde) * sqrt((3. - Z1)*(3. + Z1 + 2.*Z2));
                        double E_eq_isco = sqrt(1. - 2./r_eq_isco);
                        double eps_rad = (1. - E_eq_isco)*nu + 4.*nu*nu*(4.*const0 + 16.*const1*a_tilde*(a_tilde + 1.) + E_eq_isco - 1.);
                        p1->m = (m1 + m2) * (1. - eps_rad);

                        // now to model the kick
                        double A = 1.2e7 / velscale, B = -0.93, C = 4.57e5 / velscale, D = 3.75e6 / velscale;
                        double alpha_rad = atan2(unit_Lz, unit_Lx), beta_rad = atan2(unit_Ly, unit_Lx);
                        double sina = sin(alpha_rad), sinb = sin(beta_rad), cosa = cos(alpha_rad), cosb = cos(beta_rad);

                        double s_ex1 = cosa*cosb*usx1 - sina*usy1 + cosa*sinb*usz1;
                        double s_ex2 = cosa*cosb*usx2 - sina*usy2 + cosa*sinb*usz2;
                        double s_ey1 = sina*cosb*usy1 + cosa*usy1 + sina*sinb*usz1;
                        double s_ey2 = sina*cosb*usy2 + cosa*usy2 + sina*sinb*usz2;
                        double s_ez1 = -sinb*usx1 + cosb*usz1;
                        double s_ez2 = -sinb*usx2 + cosb*usz2;

                        double vel_m = A * q*q*(1. - q) / (opq*opq*opq*opq*opq) * (1. + B * nu);
                        double vel_perp = 0.;
                        if (m1 > m2){
                            vel_perp = C * 16. * q*q / (opq*opq*opq*opq*opq) * abs(s_ez1 - q * s_ez2);
                        } else {
                            vel_perp = C * 16. * q*q / (opq*opq*opq*opq*opq) * abs(s_ez2 - q * s_ez1);
                        }
                        double theta = reb_random_uniform(r, 0.0, 2.*M_PI); // choose random direction for the mass asymmetry kick
                        double v_kick_e1 = vel_m * cos(theta) - vel_perp * sin(theta);
                        double v_kick_e2 = vel_m * sin(theta) + vel_perp * cos(theta);
                        double DEL_e1 = 0., DEL_e2 = 0., a_perp1 = 0., a_perp2 = 0.;
                        if (m1 > m2){
                            DEL_e1 = s_ex1 - q * s_ex2, DEL_e2 = s_ey1 - q * s_ey2;
                            a_perp1 = sqrt(s_ex1*s_ex1 + s_ey1*s_ey1), a_perp2 = sqrt(s_ex2*s_ex2 + s_ey2*s_ey2);
                        } else {
                            DEL_e1 = s_ex2 - q * s_ex1, DEL_e2 = s_ey2 - q * s_ey1;
                            a_perp1 = sqrt(s_ex2*s_ex2 + s_ey2*s_ey2), a_perp2 = sqrt(s_ex1*s_ex1 + s_ey1*s_ey1);
                        }
                        double v_kick_e3 = 0.;
                        if (DEL_e1 != 0 && DEL_e2 != 0){
                            double THETA = acos((DEL_e1*cos(theta) + DEL_e2*sin(theta)) / (sqrt((DEL_e1*DEL_e1 + DEL_e2*DEL_e2))));
                            double phase = M_PI / 2.;
                            v_kick_e3 = D * cos(THETA - phase) * 16*q*q/(opq*opq*opq*opq*opq) * abs(a_perp1 - q * a_perp2);
                        }
                        
                        double rotation_determinant = cosa*cosa*cosb*cosb + sina*sina*sinb*sinb + cosa*cosa*sinb*sinb + sina*sina*cosb*cosb + cosa*cosb*sina*sinb*sinb;
                        p1->vx += ((cosa*cosb)*v_kick_e1 + (sina*sinb)*v_kick_e2 - (sinb*(sina*sina + cosa*cosa))*v_kick_e3) / rotation_determinant;
                        p1->vy += (-(sina*(cosb*cosb + sinb*sinb))*v_kick_e1 + (cosa*(cosb*cosb + sinb*sinb))*v_kick_e2) / rotation_determinant;
                        p1->vz += ((cosa*sinb)*v_kick_e1 + (sina*sinb)*v_kick_e2 + (cosa*cosa*cosb + sina*sina*sinb)*v_kick_e3) / rotation_determinant;

                        kick_vel = sqrt(v_kick_e1*v_kick_e1 + v_kick_e2*v_kick_e2 + v_kick_e3*v_kick_e3) * velscale / 1000.; // total kick velocity in km/s
                        // now calculate effective spin parameter
                        double delta = (fmax(m1, m2) - fmin(m1, m2)) / (m1 + m2);
                        double chi_1 = 0., chi_2 = 0.;
                        if (m1 > m2){
                            chi_1 = sx1 * unit_Lx + sy1 * unit_Ly + sz1 * unit_Lz;
                            chi_2 = sx2 * unit_Lx + sy2 * unit_Ly + sz2 * unit_Lz;
                        } else {
                            chi_1 = sx2 * unit_Lx + sy2 * unit_Ly + sz2 * unit_Lz;
                            chi_2 = sx1 * unit_Lx + sy1 * unit_Ly + sz1 * unit_Lz;
                        }
                        spin_eff = (1. + delta)*chi_1 / 2. + (1. - delta) * chi_2 / 2.;
                    } 
                    // now output merger details to file
                    FILE *out_file = fopen(mergers_filename, "a");
                    fprintf(out_file, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%d\t%d\t%d\n", r->t, m1*agnmass, m2*agnmass, kick_vel, spin_eff, nondim_spin, gen1, gen2, p1->generation);
                    fclose(out_file);
                    // now remove particle and end loop
                    reb_simulation_remove_particle(r, j, 1);
                    break_check = 1;
                    break;
                }
            }
        }
        if (break_check == 1){
            break;
        }
    }
}

void add_BH(struct reb_simulation* r, double distance){
    num_BH += 1;
    struct reb_particle p = {0};        // initialise BH with no spin
    double theta = reb_random_uniform(r, 0.0, 2.*M_PI);
    double R = distance;
    p.x = R * cos(theta); // x
    p.y = R * sin(theta); // y
    p.z = 0.1 * R * reb_random_uniform(r, -1., 1.); // z
    double des_vel = sqrt(1./R);
    double angle = atan2(p.y, p.x); // arctan2(y, x)
    double xprop, yprop;
    double retrograde_check = reb_random_uniform(r, 0., 1.);
    if (retrograde_check <= RETRO_PROB){        // retrograde orbiter
        xprop = sin(angle);
        yprop = -cos(angle);
    } else {                                    // prograde orbiter
        xprop = -sin(angle);
        yprop = cos(angle);
    }
    double zprop = 0.;
    double mult = des_vel * sqrt(1. / (xprop*xprop + yprop*yprop + zprop*zprop));
    xprop *= mult; yprop *= mult; zprop *= mult;
    p.vx = xprop; 
    p.vy = yprop; 
    p.vz = zprop;
    if (RAND_BH_MASSES == 0){
            p.m = 1.0e-7; // 10 solar mass BH.
    }
    else if (RAND_BH_MASSES == 1){
        p.m = lognormal_mass(r) / agnmass;
    }
    p.hash = num_BH;
    p.generation = 1;
    reb_simulation_add(r, p);
}

void output_data(struct reb_simulation* r, char* filename){
    FILE *out_file = fopen(filename, "a");
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0];
    const int N = r->N;
    const double G = r->G;
    for (int i = 1; i < N; i++){
        struct reb_particle* p = &(particles[i]); // get the particle
        const double dx = p->x-com.x;
        const double dy = p->y-com.y;
        const double dz = p->z-com.z;
        const double dvx = p->vx-com.vx;
        const double dvy = p->vy-com.vy;
        const double dvz = p->vz-com.vz;
        const double radius = sqrt(dx*dx + dy*dy + dz*dz);
        const double vr = (dx*dvx + dy*dvy + dz*dvz) / radius;
        const double mass = p->m;
        double mu = G*(com.m + mass);
        double vel = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        double a = -mu / (vel*vel - 2. * mu / radius);            // semi major axis
        
        double ex = 1. / mu * ((vel*vel - mu / radius) * dx - radius * vr * dvx );
        double ey = 1. / mu * ((vel*vel - mu / radius) * dy - radius * vr * dvy );
        double ez = 1. / mu * ((vel*vel - mu / radius) * dz - radius * vr * dvz );
        double e = sqrt(ex*ex + ey*ey + ez*ez);

        // angular momentum vectors are the cross product of position and velocity
        const double Lx = dy * dvz - dz * dvy;
        const double Ly = dz * dvx - dx * dvz;
        const double Lz = dx * dvy - dy * dvx;
        double incl = acos(Lz / sqrt(Lx*Lx + Ly*Ly + Lz*Lz)) * 180. / M_PI;

        if (BH_SPINS == 0){
            fprintf(out_file, "%.8e\t%d\t%.8e\t%.8e\t%.8e\t%.8e\n", r->t, p->hash, a, e, incl, mass);
        } else {
            double spin = sqrt(p->spin_x*p->spin_x + p->spin_y*p->spin_y + p->spin_z*p->spin_z);
            fprintf(out_file, "%.8e\t%d\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", r->t, p->hash, a, e, incl, mass, spin);
        }
    }
    fclose(out_file);
}
void output_position_data(struct reb_simulation* r, char* filename){
    FILE *out_file = fopen(filename, "a");
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0];
    const int N = r->N;
    for (int i = 1; i < N; i++){
        struct reb_particle* p = &(particles[i]); // get the particle
        const double dx = p->x-com.x;
        const double dy = p->y-com.y;
        const double dz = p->z-com.z;
        fprintf(out_file, "%.8e\t%d\t%.8e\t%.8e\t%.8e\t%.8e\n", r->t, p->hash, p->m, dx, dy, dz);
    }
    fclose(out_file);
}


void heartbeat(struct reb_simulation* r){
    check_mergers(r);
    if(reb_simulation_output_check(r, 20.*M_PI)){
        reb_simulation_output_timing(r, tmax);
    }
    if(reb_simulation_output_check(r, 1.)){
        reb_simulation_synchronize(r);
        // reb_simulation_output_orbits(r, "orbits.txt");
        output_data(r, orbits_filename);
        if (r->t < 40000){
            output_position_data(r, positions_filename);
        }
        reb_simulation_move_to_com(r); 
    }
    if (ADD_BH_RAND_TIME == 0){
        if (reb_simulation_output_check(r, ADD_BH_INTERVAL) && r->t > 0.5){
            add_BH(r, 2.);
        }
    } else if (ADD_BH_RAND_TIME == 1){  // add BHs with exponential distribution
        if (r->t > NEXT_ADD_TIME){
            add_BH(r, 2.);
            NEXT_ADD_TIME += exponential_rv(r, ADD_BH_INTERVAL);
        }
    }
    
}

void disk_forces(struct reb_simulation* r){
    const double G = r->G;
    const int N = r->N;
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0]; // calculate migration forces with respect to center of mass;

    for (int i=1; i<N; i++){
        struct reb_particle* p = &(particles[i]); // get the particle
        // first calculate the radius of the particle
        const double dx = p->x-com.x;
        const double dy = p->y-com.y;
        const double dz = p->z-com.z;
        const double dvx = p->vx-com.vx;
        const double dvy = p->vy-com.vy;
        const double dvz = p->vz-com.vz;
        const double mass = p->m, q = mass;
        double radius = sqrt(dx*dx + dy*dy + dz*dz);
        const double vr = (dx*dvx + dy*dvy + dz*dvz);   // dot product of v and r
        
        // now define constants and calculate disk properties at this radius
        double logr = log10(radius / nondim_rs), Sigma = disk_surfdens(logr), angvel = disk_angvel(logr), asp_ratio = disk_aspectratio(logr);
        double kappa = disk_opacity(logr), temp = disk_temp(logr);
        double density = Sigma / (asp_ratio * radius);  // sigma / H

        double tau = kappa * Sigma / 2, tau_eff = 3 * tau / 8 + sqrt(3) / 4 + 1. / (4 * tau);    // define optical depth params
        // start with the Type I migration as in Pardekooper (?)
        double alpha = -disk_sigma_deriv(logr), beta = -disk_temp_deriv(logr), xi = beta - (gamma_coeff - 1) * alpha; // define disk gradient properties
        double Gamma_0 = (q/asp_ratio)*(q/asp_ratio) * Sigma * radius*radius*radius*radius * angvel*angvel;
        double Gamma = 0;

        const double Lz = dx * dvy - dy * dvx;
        double a, e;
        if (Lz >= 0){       // BH is on prograde orbit
            if (MIGRATION_PRESCRIPTION == 0){
                double Theta = (c_v * Sigma * angvel * tau_eff) / (12. * M_PI * stefboltz * pow(temp, 3));
                double Gamma_iso = -0.85 - alpha - 0.9 * beta;
                double Gamma_ad = (-0.85 - alpha - 1.7 * beta + 7.9 * xi / gamma_coeff) / gamma_coeff;
                Gamma = Gamma_0 * (Gamma_ad * Theta*Theta + Gamma_iso) / ((Theta + 1)*(Theta + 1));
            }
            else if (MIGRATION_PRESCRIPTION == 1){
                double R_mu = 8.3145 / (2.016 / 1000.) * massscale;     // ideal gass constant over the mean molecular weight of H2, nondimensionalised
                // below is thermal diffusivity, chi, over a critical thermal diffusivity value, chi_c
                double chi_chi_c = (16. * (gamma_coeff - 1.) * stefboltz * temp*temp*temp / (3. * density*density * R_mu * kappa)) / (radius*radius * asp_ratio*asp_ratio * angvel);
                double fx = (sqrt(chi_chi_c / 2.) + 1. / gamma_coeff) / (sqrt(chi_chi_c / 2.) + 1.);
                double Gamma_lindblad = - (2.34 - 0.1 * alpha + 1.5 * beta) * fx;
                double Gamma_simp_corot = (0.46 - 0.96 * alpha + 1.8 * beta) / gamma_coeff;
                Gamma = Gamma_0 * (Gamma_lindblad + Gamma_simp_corot);
            }

            double Gamma_mag = Gamma / (mass * radius);         // get the net acceleration on the particle
            // add migration to the acceleration total
            p->ax += -dy * Gamma_mag / radius;
            p->ay += dx * Gamma_mag / radius;

            //// now look at Evgeni's thermal torques
            if (THERMAL_TORQUES == 1){
                double visc = 1e-2;
                double H = asp_ratio * radius;
                
                double dHdr = radius / H * (disk_aspratio_deriv(logr)); // r/H * (dln(h)/dlnr)
                double dSigmadr = radius / Sigma * -alpha;
                double drhodr = (H * dSigmadr - Sigma * dHdr) / (H*H);
                double cs = angvel * H;
                double dangveldr = -1.5 * angvel / radius;
                double dcsdr = dangveldr * H + angvel * dHdr; 
                double dPdr = -cs * (2 * density * dcsdr + cs * drhodr);
                double chi = 9. * gamma_coeff * (gamma_coeff - 1.) / 2. * visc * H*H * angvel;
                double x_c = dPdr * H*H / (3 * gamma_coeff * radius);
                double L = 0, Lc = 1;      // set our bodies luminosity value to 0 so that it has no effect
                if (ACCRETION > 0.){    // update our luminosities to have a thermal effect
                    L = ACCRETION * 4. * M_PI * G * mass * m_p * c_nbody / sigma_T;     // some proportion (given by ACCRETION) of the eddington luminosity
                    Lc = 4. * M_PI * G * mass * density * chi / gamma_coeff;            // critical luminosity given in Grishin et al (2023)
                }
                double lambda = sqrt(2. * chi / (3. * gamma_coeff * angvel));
                double Gamma_thermal = 1.61 * (gamma_coeff - 1) / gamma_coeff * x_c / lambda * (L/Lc - 1.) * Gamma_0 / asp_ratio;
                p->ax += -dy * Gamma_thermal / (mass * radius*radius);
                p->ay += dx * Gamma_thermal / (mass * radius*radius);
            }
            
            // now lets add in eccentricity/inclination damping
            double mu = G*(com.m + mass);
            double vel = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            a = -mu / (vel*vel - 2. * mu / radius);            // semi major axis
            double tdamp = asp_ratio*asp_ratio*asp_ratio*asp_ratio * com.m*com.m / (mass * Sigma * a*a * angvel);
            
            double ex = 1. / mu * ((vel*vel - mu / radius) * dx - vr * dvx);
            double ey = 1. / mu * ((vel*vel - mu / radius) * dy - vr * dvy);
            double ez = 1. / mu * ((vel*vel - mu / radius) * dz - vr * dvz);
            e = sqrt(ex*ex + ey*ey + ez*ez);
            double eps = e / asp_ratio;

            // angular momentum vectors are the cross product of position and velocity
            const double Lx = dy * dvz - dz * dvy;
            const double Ly = dz * dvx - dx * dvz;
            const double Lz = dx * dvy - dy * dvx;
            double incl = acos(Lz / sqrt(Lx*Lx + Ly*Ly + Lz*Lz)) / M_PI;
            double l = incl / asp_ratio;
            double t_i = (tdamp / 0.544) * (1. - 0.3 * l*l + 0.24 * l*l*l + 0.14 * l * eps*eps);
            double t_e = (tdamp / 0.78) * (1. - 0.14*eps*eps + 0.06*eps*eps*eps + 0.18*eps*l*l);

            // printf("\n%f\n", t_i);
            // add the damping forces to the migration acceleration
            p->ax += -2. * vr * dx / (t_e * radius*radius); // eccentricity damping
            p->ay += -2. * vr * dy / (t_e * radius*radius); // eccentricity damping
            p->az += -2. * vr * dz / (t_e * radius*radius) - dvz / t_i; // eccentricity and inclination damping
        } else {        // BH is on retrograde orbit
            double disc_vel = 1. / sqrt(radius);
            double angle = atan2(dy, dx); // arctan2(y, x)
            double disc_vel_x = -sin(angle) * disc_vel;
            double disc_vel_y = cos(angle) * disc_vel;
            double vel = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            double vel_rel = vel + disc_vel;
            double Lambda = asp_ratio * radius * vel_rel*vel_rel / (G * mass);
            double drag_force = - 4. * M_PI * log(Lambda) * G*G*mass * density / (vel_rel*vel_rel*vel_rel);
            p->ax += drag_force * (dvx - disc_vel_x);
            p->ay += drag_force * (dvy - disc_vel_y);
            p->az += drag_force * dvz;

            double mu = G*(com.m + mass);
            double ex = 1. / mu * ((vel*vel - mu / radius) * dx - vr * dvx);
            double ey = 1. / mu * ((vel*vel - mu / radius) * dy - vr * dvy);
            double ez = 1. / mu * ((vel*vel - mu / radius) * dz - vr * dvz);
            e = sqrt(ex*ex + ey*ey + ez*ez);

            a = -mu / (vel*vel - 2. * mu / radius);
        }

        

        if (a*(1. - e) < r_isco){     // particle merged with SMBH
            printf("\nParticle merged with SMBH");
            com.m += mass;
            reb_simulation_remove_particle_by_hash(r, p->hash, 1);
        } else if (e > 0.99){
            printf("\nParticle ejected!");
            reb_simulation_remove_particle_by_hash(r, p->hash, 1);
        }
    }
}

void init_conds(int N, int mass, struct reb_simulation* r){
    agnmass = pow(10., mass);
    double Nr_s = 1e3;
    r_g = G_pc * agnmass / 9e10; r_s = 2. * r_g;
    lenscale = Nr_s * r_s;
    lenscale_m = lenscale * pc_to_m;
    lenscale_rs = lenscale / r_s;
    nondim_rs = r_s / lenscale;
    r_isco = 3. * nondim_rs;
    timescale = sqrt(pow(lenscale_m, 3.) / (G_pc * agnmass * pc_to_m * 1000.*1000.));
    velscale = sqrt(G_pc * agnmass / lenscale) * 1000.;
    massscale = agnmass * M_odot;
    stefboltz = 5.67e-8 * lenscale_m*lenscale_m * timescale;         // non-dimensionalised boltzmann constant
    c_v = 14304. * massscale;       // specific heat capacity of H2 gas, non-dimensionalised
    c_nbody = c / velscale;         // non-dimensionalised speed of light
    m_p = 1.67e-27 / massscale;     // proton mass
    sigma_T = 6.65e-29 / (lenscale_m*lenscale_m);       // thompson scattering cross-section of an electron

    // double theta, dist, R, des_vel, angle, xprop, yprop, zprop, mult;
    struct reb_particle p = {0}; // smbh
    p.m = 1.;
    reb_simulation_add(r, p);
    // uniformly (and randomly) distribute points in the unit disk
    for (int i = 1; i <= N; i++){   // start from i=1 because we want the SMBH to be at i=0
        double dist = reb_random_uniform(r, 0.5*0.5*0.5, 1.);
        double R = 2. * pow(dist, 1./3.);
        add_BH(r, R);
    }

    if (MERGER_CRITERION == 1){
        MUTUAL_HILL_PROP = 0.3;
    }
    if (ADD_BH_RAND_TIME == 1){
        NEXT_ADD_TIME = exponential_rv(r, ADD_BH_INTERVAL);
    }
}



int main(int argc, char* argv[]){
    // srand(time(NULL));
    int mass;
    double f_edd, alpha;
    if (scanf("%d", &mass)); 
    if (scanf("%lf", &f_edd)); 
    if (scanf("%lf", &alpha));
    create_filenames(mass, f_edd, alpha);
    populate_spline_arrays();
    eval_splines();
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup simulation features
    THERMAL_TORQUES             = 1;
    MIGRATION_PRESCRIPTION      = 1;        // set to jimenez and masset migration torques
    RAND_BH_MASSES              = 1;        // randomly sample bh masses
    MUTUAL_HILL_PROP            = 0.65;     //
    MERGER_CRITERION            = 1;        // set to Li et al (2023)
    MERGER_KICKS                = 1;        //
    ACCRETION                   = 0.1;      // proportion of eddington luminosity onto the satellite BHs
    ADD_BH_RAND_TIME            = 1;        // add BHs with exponential time distribution
    ADD_BH_INTERVAL             = 1e4;      // add BHs with a mean of 10k time steps
    BH_SPINS                    = 1;        // model spins of BHs during mergers

    // now set up integration parameters
    // r->integrator           = REB_INTEGRATOR_BS;
    r->integrator           = REB_INTEGRATOR_IAS15;
    // r->ri_ias15.epsilon     = 5e-10;
    r->dt                   = 5e-3;    
    r->additional_forces    = disk_forces;     //Set function pointer to add dissipative forces.
    r->heartbeat            = heartbeat;        // checks for mergers and outputs data
    r->force_is_velocity_dependent = 1;
    tmax                    = 80000.;       // multiply by ~1.4 to get the approximate number of years (1e8M SMBH)
    r->rand_seed            = 2399;

    // Initial conditions
    int initial_BH = 10;
    init_conds(initial_BH, mass, r);

    reb_simulation_move_to_com(r);          

    // delete previous output files
    remove(orbits_filename); 
    remove(positions_filename);
    remove(mergers_filename);

    reb_simulation_integrate(r, tmax);
}
