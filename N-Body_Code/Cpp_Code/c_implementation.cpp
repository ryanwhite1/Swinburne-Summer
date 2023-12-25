#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <vector>
#include "spline.h"
#include <cmath>

double norm(std::vector<double> vec1, std::vector<double> vec2){
    double accum;
    for (int i = 0; i < vec1.size(); i++){accum += vec1[i] * vec2[i];}
    return sqrt(accum);
}
std::vector<double> cross_product(std::vector<double> vec1, std::vector<double> vec2){
    std::vector<double> result(vec1.size());
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    return result;
}
std::vector<double> vec_add(std::vector<double> vec1, std::vector<double> vec2){
    // element-wise addition of vector elements for like-size vectors
    std::vector<double> output(vec1.size());
    for (int i = 0; i < vec1.size(); i++){
        output[i] = vec1[i] + vec2[i];
    }
    return output;
}
std::vector<double> vec_mult(std::vector<double> vec1, std::vector<double> vec2){
    // element-wise multiplication of vector elements for two like-size vectors
    std::vector<double> output(vec1.size());
    for (int i = 0; i < vec1.size(); i++){
        output[i] = vec1[i] * vec2[i];
    }
    return output;
}
std::vector<double> vec_scalar(std::vector<double> vec1, double mult){
    // element-wise multiplication of vector elements for a single number multiplier
    std::vector<double> output(vec1.size());
    for (int i = 0; i < vec1.size(); i++){
        output[i] = vec1[i] * mult;
    }
    return output;
}
double eccentricity(std::vector<double> position, std::vector<double> velocity){
    double radius = norm(position, position);
    std::vector<double> e_vec = cross_product(velocity, cross_product(position, velocity));
    e_vec = vec_add(e_vec, vec_scalar(position, 1 / radius));
    return norm(e_vec, e_vec);
}
double inclination(std::vector<double> position, std::vector<double> velocity){
    std::vector<double> ang_momentum = cross_product(position, velocity);
    return acos(ang_momentum[2] / norm(ang_momentum, ang_momentum));
}



double G_pc = 4.3e-3, M_odot = 1.98e30, c = 3e8, pc_to_m = 3.086e16;

std::vector<double> sigma_data_r = {0.5, 1, 1.3, 1.5, 1.7, 2, 2.6, 3, 3.5, 4, 5.5, 7};
std::vector<double> sigma_data = {3.6, 4.1, 4.5, 4.9, 5.1, 4.8, 5.9, 5.9, 5, 4, 2, 0};
tk::spline log_sigma_spline(sigma_data_r, sigma_data);

std::vector<double> temp_data_r = {0.5, 1, 1.7, 2, 2.5, 3, 4, 5, 6, 7};
std::vector<double> temp_data = {5.95, 5.75, 5.4, 5.1, 5, 4.75, 4, 3.2, 2.5, 1.8};
tk::spline log_temp_spline(temp_data_r, temp_data);

std::vector<double> aratio_data_r = {0.5, 0.7, 1, 1.4, 1.7, 1.9, 2, 2.2, 2.6, 3, 3.1, 3.25, 3.5, 4, 5, 6, 7};
std::vector<double> aratio_data = {-0.7, -0.8, -0.92, -1.3, -1.45, -1.4, -1.4, -1.7, -2.1, -2.15, -2.1, -2, -1.8, -1.6, -1.05, -0.6, -0.1};
tk::spline log_aratio_spline(aratio_data_r, aratio_data);

std::vector<double> opacity_data_r = {0.5, 1, 1.5, 1.7, 2, 2.5, 3, 3.5, 4, 4.1, 4.2, 4.5, 5, 5.5, 6, 7};
std::vector<double> opacity_data = {-0.4, -0.4, -0.4, -0.4, 0, -0.15, 0, -0.3, -0.4, -0.38, -0.4, -2, -3.1, -3.12, -3.15, -3.15};
tk::spline log_opacity_spline(opacity_data_r, opacity_data);

class AGNDisk{
    public:
        double mass, r_g, r_s, lenscale, lenscale_m, lenscale_rs, nondim_rs, timescale, velscale, massscale;
        AGNDisk(double, double);
        double disk_temp (double logr){
            return pow(10, log_temp_spline(logr));
        }
        double disk_surfdens(double logr){
            double sigma = pow(10, log_sigma_spline(logr) + 1); // +1 in power to go from cgs to SI
            return sigma * lenscale_m * lenscale_m / massscale;
        }
        double disk_angvel(double logr){
            double omega = sqrt(G_pc * mass / pow((pow(10, logr) * r_s), 3)) * 1000;
            return omega * timescale;
        }
        double disk_aspectratio(double logr){
            return pow(10, log_aratio_spline(logr));
        }
        double disk_opacity(double logr){
            double kappa = pow(10, log_opacity_spline(logr) - 1); // -1 to go from cgs to SI
            return kappa * massscale / (lenscale_m * lenscale_m);
        }
        double disk_dens(double logr){
            return disk_surfdens(logr) / (2 * disk_aspectratio(logr) * pow(10, logr));
        }
        std::vector<double> disk_forces(double mass, std::vector<double> position, std::vector<double> vel){
            std::vector<double> accel(3);
            // first calculate the radius of the particle
            double radius = norm(position, position);
            // now define constants and calculate disk properties at this radius
            double stefboltz = 5.67e-8 * lenscale_m * lenscale_m * timescale;
            double gamma = 5/3, c_v = 14304 * massscale, q = mass;
            double logr = log10(radius / nondim_rs), Sigma = disk_surfdens(logr), angvel = disk_angvel(logr), asp_ratio = disk_aspectratio(logr);
            double kappa = disk_opacity(logr), temp = disk_temp(logr);

            double tau = kappa * Sigma / 2, tau_eff = 3 * tau / 8 + sqrt(3) / 4 + 1 / (4 * tau);    // define optical depth params
            // start with the Type I migration as in Pardekooper (?)
            double alpha = -log_sigma_spline.deriv(1, logr), beta = -log_temp_spline.deriv(1, logr), xi = beta - (gamma - 1) * alpha; // define disk gradient properties
            double Theta = (c_v * Sigma * angvel * tau_eff) / (12 * M_PI * stefboltz * pow(temp, 3));
            double Gamma_0 = pow(q / asp_ratio, 2) * Sigma * pow(radius, 4) * pow(angvel, 2);
            double Gamma_iso = -0.85 - alpha - 0.9 * beta;
            double Gamma_ad = (-0.85 - alpha - 1.7 * beta + 7.9 * xi / gamma) / gamma;
            double Gamma = Gamma_0 * (Gamma_ad * pow(Theta, 2) + Gamma_iso) / pow(Theta + 1, 2);
            double Gamma_mag = Gamma / (mass * radius);         // get the net acceleration on the particle
            std::vector<double> theta_vec = {- position[1], position[0], 0}; theta_vec = vec_scalar(theta_vec, 1 / radius);
            std::vector<double> migration = vec_scalar(theta_vec, Gamma_mag);
            accel = vec_add(accel, migration);      // add migration to the acceleration total

            // now lets add in eccentricity/inclination damping

            return accel;
        }
};
AGNDisk::AGNDisk (double smbhmass, double lenscale){
    mass = smbhmass;
    lenscale = lenscale;
    r_g = G_pc * mass / 9e10; r_s = r_g * 2;
    lenscale_m = lenscale * pc_to_m;
    lenscale_rs = lenscale / r_s;
    nondim_rs = r_s / lenscale;
    timescale = sqrt(pow(lenscale_m, 3) / (G_pc * mass * pc_to_m * 1000 * 1000));
    velscale = sqrt(G_pc * mass / lenscale) * 1000;
    massscale = mass * M_odot;
}


int main(){
    double Tmax = 20000, dt = 0.02;
    int nt = int(Tmax / dt) + 1, NsBH = 10;
    std::vector<double> NsBHMasses(NsBH, 1);
    double SMBHMass = 1e8, r_s = 2 * G_pc * SMBHMass / 9e10;  // 2GM / c^2     units of pc
    double Nr_s = 1e3;      // number of schwarzschild radii to initialise the sim with respect to
    double lenscale = Nr_s * r_s;
    AGNDisk agn (SMBHMass, lenscale);
    double h = agn.disk_dens(3);
    std::cout << h;
    return 0;
}