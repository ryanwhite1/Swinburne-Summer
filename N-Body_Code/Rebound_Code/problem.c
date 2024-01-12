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
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "spline.h"

double* tau_a;     /**< Migration timescale in years for all particles */
double* tau_e;     /**< Eccentricity damping timescale in years for all particles */
double tmax;
double agnmass, r_g, r_s, massscale, lenscale, lenscale_m, lenscale_rs, nondim_rs, timescale, velscale, stefboltz, c_v;
const double G_pc = 4.3e-3, M_odot = 1.98e30, c = 3e8, pc_to_m = 3.086e16, gamma_coeff = 5./3.; 
int num_BH;

double temp_deriv_coeffs[3] = {0, 0, 0}, sigma_deriv_coeffs[] = {0, 0, 0}, asp_deriv_coeffs[] = {0, 0, 0};
const int n_sigma = 12, n_temp = 10, n_aratio = 17, n_opacity = 16;
double sigma_data_r[] = {0.5, 1., 1.3, 1.5, 1.7, 2., 2.6, 3, 3.5, 4., 5.5, 7.};
double sigma_data[] = {3.6, 4.1, 4.5, 4.9, 5.1, 4.8, 5.9, 5.9, 5., 4., 2., 0.};
double log_sigma_spline[12];

double temp_data_r[] = {0.5, 1., 1.7, 2., 2.5, 3., 4., 5., 6., 7.};
double temp_data[] = {5.95, 5.75, 5.4, 5.1, 5., 4.75, 4., 3.2, 2.5, 1.8};
double log_temp_spline[10];

double aratio_data_r[] = {0.5, 0.7, 1., 1.4, 1.7, 1.9, 2., 2.2, 2.6, 3., 3.1, 3.25, 3.5, 4., 5., 6., 7.};
double aratio_data[] = {-0.7, -0.8, -0.92, -1.3, -1.45, -1.4, -1.4, -1.7, -2.1, -2.15, -2.1, -2., -1.8, -1.6, -1.05, -0.6, -0.1};
double log_aratio_spline[17];

double opacity_data_r[] = {0.5, 1., 1.5, 1.7, 2., 2.5, 3., 3.5, 4., 4.1, 4.2, 4.5, 5., 5.5, 6., 7.};
double opacity_data[] = {-0.4, -0.4, -0.4, -0.4, 0., -0.15, 0., -0.3, -0.4, -0.38, -0.4, -2., -3.1, -3.12, -3.15, -3.15};
double log_opacity_spline[16];


double uniform(double lower, double upper){
    // with thanks to https://stackoverflow.com/questions/63981013/generating-random-double-between-1-and-100-in-c?noredirect=1&lq=1
    double random = ((double) rand() / RAND_MAX) * (upper - lower) + lower;
    return random;
}


void eval_splines(){
    spline(sigma_data_r, sigma_data, n_sigma, log_sigma_spline);
    spline(temp_data_r, temp_data, n_temp, log_temp_spline);
    spline(aratio_data_r, aratio_data, n_aratio, log_aratio_spline);
    spline(opacity_data_r, opacity_data, n_opacity, log_opacity_spline);
}
double disk_surfdens(double logr){
    double sigma = pow(10., splint(sigma_data_r, sigma_data, log_sigma_spline, n_sigma, logr) + 1.); // +1 in power to go from cgs to SI
    return sigma * lenscale_m * lenscale_m / massscale;
}
double disk_temp(double logr){
    return pow(10., splint(temp_data_r, temp_data, log_temp_spline, n_temp, logr));
}
double disk_aspectratio(double logr){
    return pow(10., splint(aratio_data_r, aratio_data, log_aratio_spline, n_aratio, logr));
}
double disk_opacity(double logr){
    double kappa = pow(10., splint(opacity_data_r, opacity_data, log_opacity_spline, n_opacity, logr) - 1.); // -1 to go from cgs to SI
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
    return splderiv(sigma_data_r, sigma_data, log_sigma_spline, n_sigma, logr);
}
double disk_temp_deriv(double logr){
    return splderiv(temp_data_r, temp_data, log_temp_spline, n_temp, logr);
}
double disk_aspratio_deriv(double logr){
    return splderiv(aratio_data_r, aratio_data, log_aratio_spline, n_aratio, logr);
}


void check_mergers(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0];
    const int N = r->N;
    int break_check = 0;
    for (int i = 1; i < N - 1; i++){
        struct reb_particle* p1 = &(particles[i]); // get the particle
        double m1 = p1->m;
        const double dx1 = p1->x-com.x; const double dy1 = p1->y-com.y; const double dz1 = p1->z-com.z;
        const double r1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        for (int j = i + 1; j < N; j++){
            struct reb_particle* p2 = &(particles[j]); // get the particle
            double m2 = p2->m;
            const double dx2 = p2->x-com.x; const double dy2 = p2->y-com.y; const double dz2 = p2->z-com.z;
            const double r2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
            double R_mH = pow((m1 + m2) / (3 * com.m), 1./3.) * (r1 + r2) / 2.;   // mutual hill radius
            double dist = sqrt((dx1-dx2)*(dx1-dx2) + (dy1-dy2)*(dy1-dy2) + (dz1-dz2)*(dz1-dz2));
            if (dist < (1. * R_mH)){
                // now to check if their relative kinetic energy is low enough for capture
                const double dvx1 = p1->vx-com.vx;
                const double dvy1 = p1->vy-com.vy;
                const double dvz1 = p1->vz-com.vz;
                const double dvx2 = p2->vx-com.vx;
                const double dvy2 = p2->vy-com.vy;
                const double dvz2 = p2->vz-com.vz;
                double reduced_mass = 1. / (1. / m1 + 1. / m2);
                double rel_kin_energy = 0.5 * reduced_mass * ((dvx1-dvx2)*(dvx1-dvx2) + (dvy1-dvy2)*(dvy1-dvy2) + (dvz1-dvz2)*(dvz1-dvz2));
                double binding_energy = m1 * m2 / (2. * R_mH);
                if (rel_kin_energy < binding_energy){
                    puts("Merger!\n");
                    p1->x = (m1 * p1->x + m2 * p2->x) / (m1 + m2);
                    p1->y = (m1 * p1->y + m2 * p2->y) / (m1 + m2);
                    p1->z = (m1 * p1->z + m2 * p2->z) / (m1 + m2);
                    p1->vx = (m1 * p1->vx + m2 * p2->vx) / (m1 + m2);
                    p1->vy = (m1 * p1->vy + m2 * p2->vy) / (m1 + m2);
                    p1->vz = (m1 * p1->vz + m2 * p2->vz) / (m1 + m2);
                    p1->m = 0.95 * (m1 + m2);

                    // include a kick in a random direction
                    double q = fmin(m1 / m2, m2 / m1), opq = 1. + q;
                    double A = 1.2 * 1e7 / velscale, B = -0.93;
                    double vel_kick = A * q*q*(1. - q) / (opq*opq*opq*opq*opq) * (1. + B * q / (opq*opq));
                    double angle = reb_random_uniform(r, 0.0, 2.*M_PI);
                    double xprop = sin(angle);
                    double yprop = -cos(angle);
                    double mult = sqrt(vel_kick*vel_kick / (xprop*xprop + yprop*yprop));
                    xprop *= mult; yprop *= mult;
                    p1->vx += xprop; 
                    p1->vy += yprop; 
                    printf("%e %e %e %e %e \n", A, q, vel_kick, xprop, p1->vx);

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

void add_BH(struct reb_simulation* r){
    num_BH += 1;
    struct reb_particle p = {0};
    double theta = reb_random_uniform(r, 0.0, 2.*M_PI);
    double R = 1;
    p.x = R * cos(theta); // x
    p.y = R * sin(theta); // y
    p.z = 0.0; // z
    double des_vel = sqrt(1./R);
    double angle = atan2(p.y, p.x); // arctan2(y, x)
    double xprop = sin(angle);
    double yprop = -cos(angle);
    double zprop = 0.0;
    double mult = sqrt(des_vel*des_vel / (xprop*xprop + yprop*yprop + zprop*zprop));
    xprop *= mult; yprop *= mult; zprop *= mult;
    p.vx = xprop; 
    p.vy = yprop; 
    p.vz = zprop;
    p.m = 1.0e-7; // 10 solar mass BH.
    p.hash = num_BH;
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
        double incl = acos(Lz / sqrt(Lx*Lx + Ly*Ly + Lz*Lz));

        // char out_text[] = {r->t, '\t', p->hash, '\t', a, '\t', e, '\t', incl, '\t', mass, '\n'};
        // fprintf(out_file, out_text);
        fprintf(out_file, "%.8e\t%d\t%.8e\t%.8e\t%.8e\t%.8e\n", r->t, p->hash, a, e, incl, mass);
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
        output_data(r, "orbits.txt");
        reb_simulation_move_to_com(r); 
    }
    if (reb_simulation_output_check(r, 10000.) && r->t > 1){
        add_BH(r);
    }
}

void disk_forces(struct reb_simulation* r){
    const double G = r->G;
    const int N = r->N;
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0]; // calculate migration forces with respect to center of mass;

    for (int i=1; i<N; i++){
        // printf("%d", i);
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

        double tau = kappa * Sigma / 2, tau_eff = 3 * tau / 8 + sqrt(3) / 4 + 1. / (4 * tau);    // define optical depth params
        // start with the Type I migration as in Pardekooper (?)
        double alpha = -disk_sigma_deriv(logr), beta = -disk_temp_deriv(logr), xi = beta - (gamma_coeff - 1) * alpha; // define disk gradient properties
        double Theta = (c_v * Sigma * angvel * tau_eff) / (12. * M_PI * stefboltz * pow(temp, 3));
        double Gamma_0 = (q/asp_ratio)*(q/asp_ratio) * Sigma * radius*radius*radius*radius * angvel*angvel;

        double Gamma_iso = -0.85 - alpha - 0.9 * beta;
        double Gamma_ad = (-0.85 - alpha - 1.7 * beta + 7.9 * xi / gamma_coeff) / gamma_coeff;
        double Gamma = Gamma_0 * (Gamma_ad * Theta*Theta + Gamma_iso) / ((Theta + 1)*(Theta + 1));
        double Gamma_mag = Gamma / (mass * radius);         // get the net acceleration on the particle
        // add migration to the acceleration total
        p->ax += dy * Gamma_mag / radius;
        p->ay += -1. * dx * Gamma_mag / radius;

        //// now look at Evgeni's thermal torques
        double visc = 1e-2;
        double H = asp_ratio * radius;
        double density = Sigma / (asp_ratio * radius);  // sigma / H
        double dHdr = radius / H * (disk_aspratio_deriv(logr)); // r/H * (dln(h)/dlnr)
        double dSigmadr = radius / Sigma * -alpha;
        double drhodr = (H * dSigmadr - Sigma * dHdr) / (H*H);
        double cs = angvel * H;
        double dangveldr = -1.5 * angvel / radius;
        double dcsdr = dangveldr * H + angvel * dHdr; 
        double dPdr = -cs * (2 * density * dcsdr + cs * drhodr);
        double chi = 9. * gamma_coeff * (gamma_coeff - 1.) / 2. * visc * H*H * angvel;
        double x_c = dPdr * H*H / (3 * gamma_coeff * radius);
        // double Lc = 4. * M_PI * G * mass * density * chi / gamma_coeff;
        int L = 0, Lc = 1;
        double lambda = sqrt(2. * chi / (3. * gamma_coeff * angvel));
        double Gamma_thermal = 1.61 * (gamma_coeff - 1) / gamma_coeff * x_c / lambda * (L/Lc - 1.) * Gamma_0 / asp_ratio;
        // printf("%e", Gamma_thermal);
        p->ax += dy * Gamma_thermal / (mass * radius*radius);
        p->ay += -1. * dx * Gamma_thermal / (mass * radius*radius);


        // now lets add in eccentricity/inclination damping
        double mu = G*(com.m + mass);
        double vel = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        double a = -mu / (vel*vel - 2. * mu / radius);            // semi major axis
        double tdamp = asp_ratio*asp_ratio*asp_ratio*asp_ratio * com.m*com.m / (mass * Sigma * a*a * angvel);
        
        double ex = 1. / mu * ((vel*vel - mu / radius) * dx - vr * dvx);
        double ey = 1. / mu * ((vel*vel - mu / radius) * dy - vr * dvy);
        double ez = 1. / mu * ((vel*vel - mu / radius) * dz - vr * dvz);
        double e = sqrt(ex*ex + ey*ey + ez*ez);
        double eps = e / asp_ratio;


        // angular momentum vectors are the cross product of position and velocity
        const double Lx = dy * dvz - dz * dvy;
        const double Ly = dz * dvx - dx * dvz;
        const double Lz = dx * dvy - dy * dvx;
        double incl = acos(Lz / sqrt(Lx*Lx + Ly*Ly + Lz*Lz)) - M_PI;
        double l = incl / asp_ratio;
        double t_i = (tdamp / 0.544) * (1. - 0.3 * l*l + 0.24 * l*l*l + 0.14 * l * eps*eps);
        double t_e = (tdamp / 0.78) * (1. - 0.14*eps*eps + 0.06*eps*eps*eps + 0.18*eps*l*l);
        // printf("%f %f %f\n", tdamp, t_e, l);

        // add the damping forces to the migration acceleration
        p->ax += -2. * vr * dx / (t_e * radius*radius); // eccentricity damping
        p->ay += -2. * vr * dy / (t_e * radius*radius); // eccentricity damping
        p->az += -2. * vr * dz / (t_e * radius*radius) - dvz / t_i; // eccentricity and inclination damping
    }
}

void init_conds(int N, struct reb_simulation* r){
    agnmass = 1e8;
    double Nr_s = 1e3;
    r_g = G_pc * agnmass / 9e10; r_s = 2. * r_g;
    lenscale = Nr_s * r_s;
    lenscale_m = lenscale * pc_to_m;
    lenscale_rs = lenscale / r_s;
    nondim_rs = r_s / lenscale;
    timescale = sqrt(pow(lenscale_m, 3.) / (G_pc * agnmass * pc_to_m * 1000.*1000.));
    velscale = sqrt(G_pc * agnmass / lenscale) * 1000.;
    massscale = agnmass * M_odot;
    stefboltz = 5.67e-8 * lenscale_m*lenscale_m * timescale;         // non-dimensionalised boltzmann constant
    c_v = 14304. * massscale;

    // double theta, dist, R, des_vel, angle, xprop, yprop, zprop, mult;
    struct reb_particle p = {0}; // smbh
    p.m = 1.;
    reb_simulation_add(r, p);
    // uniformly (and randomly) distribute points in the unit disk
    for (int i = 1; i <= N; i++){   // start from i=1 because we want the SMBH to be at i=0
        struct reb_particle p = {0};
        double theta = reb_random_uniform(r, 0.0, 2.*M_PI);
        double dist = reb_random_uniform(r, 0.5*0.5*0.5, 1.0);
        double R = pow(dist, 1.0/3);
        p.x = R * cos(theta); // x
        p.y = R * sin(theta); // y
        p.z = 0.0; // z
        double des_vel = sqrt(1./R);
        double angle = atan2(p.y, p.x); // arctan2(y, x)
        double xprop = sin(angle);
        double yprop = -cos(angle);
        double zprop = 0.0;
        double mult = sqrt(des_vel*des_vel / (xprop*xprop + yprop*yprop + zprop*zprop));
        xprop *= mult; yprop *= mult; zprop *= mult;
        p.vx = xprop; 
        p.vy = yprop; 
        p.vz = zprop;
        p.m = 1.0e-7; // 10 solar mass BH.
        p.hash = i;
        reb_simulation_add(r, p);
    }
}



int main(int argc, char* argv[]){
    // srand(time(NULL));
    eval_splines();
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    // r->integrator           = REB_INTEGRATOR_MERCURIUS;
    // r->dt                   = 1e-2; 
    // r->ri_ias15.min_dt      = 1e-4 * r->dt;
    // r->ri_mercurius.r_crit_hill = 5;        // 5 hill radii critical integrator switch over

    // r->integrator           = REB_INTEGRATOR_BS;
    r->integrator           = REB_INTEGRATOR_IAS15;
    // r->ri_ias15.epsilon     = 5e-10;
    r->dt                   = 5e-3;    
    r->additional_forces    = disk_forces;     //Set function pointer to add dissipative forces.
    r->heartbeat            = heartbeat;        // checks for mergers and outputs data
    r->force_is_velocity_dependent = 1;
    tmax                    = 80000.;
    r->rand_seed            = 2399;

    num_BH = 10;
    // Initial conditions
    init_conds(num_BH, r);

    reb_simulation_move_to_com(r);          

    remove("orbits.txt"); // delete previous output file

    reb_simulation_integrate(r, tmax);
}


