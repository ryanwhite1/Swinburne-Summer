#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "spline.h"
#include <cmath>
#include <random>
#include <ctime>
#include <algorithm>
#include <iterator>

std::vector<int> set_diff(std::vector<int> vec1, std::vector<int> vec2){
    std::vector<int> diff;
    std::set_difference(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                        std::inserter(diff, diff.begin()));
    return diff;
}
std::vector<int> int_arange(int start, int stop, int step){
    int nt = int((stop - start) / step);
    std::vector<int> range(nt);
    for (int i = 0; i < nt; i++){
        range[i] = start + i * step;
    }
    return range;
}
double norm(std::vector<double> vec1, std::vector<double> vec2){
    double accum;
    for (int i = 0; i < vec1.size(); i++){accum += vec1[i] * vec2[i];}
    return sqrt(accum);
}
double dist_norm(std::vector<double> vec1, std::vector<double> vec2){
    double accum;
    for (int i = 0; i < vec1.size(); i++){accum += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);}
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
double dot_product(std::vector<double> vec1, std::vector<double> vec2){
    std::vector<double> temp = vec_mult(vec1, vec2);
    double sum;
    for (int i = 0; i < temp.size(); i++){sum += temp[i];}
    return sum;
}
template<typename T>
T uniform(T range_from, T range_to) {
    // thanks to user Shoe on https://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_real_distribution<T>    distr(range_from, range_to);
    return distr(generator);
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
double semi_major_axis(std::vector<double> position, std::vector<double> velocity){
    double radius = norm(position, position), vel_mag = norm(velocity, velocity);
    return radius / (2 - radius * vel_mag*vel_mag);    // https://physics.stackexchange.com/questions/295431/how-can-i-calculate-the-semi-major-axis-from-velocity-position-and-pull
    // also found by rearranging the expression v^2 = GM(2/r - 1/a)
}
void save_1d_vec(std::vector<double> vector, std::string filename){
    std::ofstream outputfile(filename);
    outputfile << "[";
    for (int i = 0; i < vector.size() - 1; i++){
        outputfile << vector[i] << ", ";
    }
    outputfile << vector[vector.size() - 1] << "]";
    outputfile.close();
}
void save_3d_vec(std::vector<std::vector<std::vector<double>>> vector, std::string filename){
    std::ofstream outputfile(filename);
    int lvl1 = vector.size(), lvl2 = vector[0].size(), lvl3 = vector[0][0].size();
    outputfile << "[";
    for (int i = 0; i < lvl1; i++){
        outputfile << "[";
        for (int j = 0; j < lvl2; j++){
            outputfile << "[";
            for (int k = 0; k < lvl3 - 1; k++){
                outputfile << vector[i][j][k] << ", ";
            }
            outputfile << vector[i][j][lvl3 - 1];
            if (j == lvl2 - 1){
                outputfile << "]";
            } else {
                outputfile << "], ";
            }
        }
        if (i == lvl1 - 1){
            outputfile << "]";
        } else {
            outputfile << "], ";
        }
        
    }
    outputfile << "]";
    outputfile.close();
}



const double G_pc = 4.3e-3, M_odot = 1.98e30, c = 3e8, pc_to_m = 3.086e16;


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
        double agnmass, r_g, r_s, lenscale, lenscale_m, lenscale_rs, nondim_rs, timescale, velscale, massscale, stefboltz;
        AGNDisk(double, double);
        double disk_temp (double logr){
            return pow(10, log_temp_spline(logr));
        }
        double disk_surfdens(double logr){
            double sigma = pow(10, log_sigma_spline(logr) + 1); // +1 in power to go from cgs to SI
            return sigma * lenscale_m * lenscale_m / massscale;
        }
        double disk_angvel(double logr){
            double omega = sqrt(G_pc * agnmass / pow((pow(10, logr) * r_s), 3)) * 1000;
            return omega * lenscale / velscale;
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
            std::vector<double> accel(3, 0);
            // first calculate the radius of the particle
            double radius = norm(position, position);
            
            // now define constants and calculate disk properties at this radius
            
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
            double Gamma = Gamma_0 * (Gamma_ad * Theta*Theta + Gamma_iso) / ((Theta + 1)*(Theta + 1));
            double Gamma_mag = Gamma / (mass * radius);         // get the net acceleration on the particle
            std::vector<double> theta_vec = {- position[1], position[0], 0}; theta_vec = vec_scalar(theta_vec, 1 / radius);
            std::vector<double> migration = vec_scalar(theta_vec, Gamma_mag);
            accel = vec_add(accel, migration);      // add migration to the acceleration total

            // now lets add in eccentricity/inclination damping
            double a, tdamp, t_i, e, eps, i, l;
            a = semi_major_axis(position, vel);
            tdamp = pow(asp_ratio, 4) / (mass * Sigma * a*a * angvel);  // we're assuming M_* = 1 (nbody units) on the numerator, so don't need to include
            
            e = eccentricity(position, vel);        // https://astronomy.stackexchange.com/questions/29005/calculation-of-eccentricity-of-orbit-from-velocity-and-radius
            eps = e / asp_ratio;
            i = inclination(position, vel);
            l = i / asp_ratio;
            t_i = (tdamp / 0.544) * (1 - 0.3 * l*l + 0.24 * l*l*l + 0.14 * l * eps*eps);
            std::vector<double> e_damp = vec_scalar(position, -2 * dot_product(vel, position) / tdamp); // eccentricity damping
            std::vector<double> i_damp = {0, 0, -vel[2] / t_i};
            accel = vec_add(accel, vec_add(e_damp, i_damp));    // add the damping forces to the migration acceleration
            return accel;
        }
};
AGNDisk::AGNDisk (double smbhmass, double lengthscale){
    agnmass = smbhmass;
    lenscale = lengthscale;
    r_g = G_pc * agnmass / 9e10; r_s = r_g * 2;
    lenscale_m = lenscale * pc_to_m;
    lenscale_rs = lenscale / r_s;
    nondim_rs = r_s / lenscale;
    timescale = sqrt(pow(lenscale_m, 3) / (G_pc * agnmass * pc_to_m * 1000*1000));
    velscale = sqrt(G_pc * agnmass / lenscale) * 1000;
    massscale = agnmass * M_odot;
    stefboltz = 5.67e-8 * lenscale_m*lenscale_m * timescale;         // non-dimensionalised boltzmann constant
}

std::vector<std::vector<double>> nbody_accel(std::vector<double> masses, std::vector<std::vector<double>> pos){
    int N = masses.size();
    double dist, mag;
    std::vector<std::vector<double>> accel(N, std::vector<double> (3, 0));
    double softening = 0.01;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){continue;}
            else {
                dist = dist_norm(pos[i], pos[j]);
                mag = -masses[j] / pow(dist*dist + softening*softening, 3./2);
                for (int k = 0; k < 3; k++){
                    accel[i][k] += mag * (pos[i][k] - pos[j][k]);
                }
            }
        }
    }
    return accel;
}

void nbody_timestep(std::vector<std::vector<double>> &pos, std::vector<std::vector<double>> &vel, std::vector<std::vector<double>> &accel,
                        double dt, std::vector<double> masses, AGNDisk agn, std::vector<int> &captured){
    
    std::vector<int> check_inds = set_diff(int_arange(0, masses.size(), 1), captured);
    // leapfrog integration
    int size = masses.size();
    for (int i = 0; i < size; i++){
        vel[i] = vec_add(vel[i], vec_scalar(accel[i], 0.5 * dt));
        pos[i] = vec_add(pos[i], vec_scalar(vel[i], dt));
    }
    accel = nbody_accel(masses, pos);
    for (int i = 0; i < size; i++){
        if(std::find(check_inds.begin(), check_inds.end(), i) != check_inds.end()) {    // check if check_inds contains this index 'i'
            accel[i] = vec_add(accel[i], agn.disk_forces(masses[i], pos[i], vel[i]));
            vel[i] = vec_add(vel[i], vec_scalar(accel[i], 0.5 * dt));
        } else { // if not, then we don't care to add forces to this captured BH
            continue;
        }
    }
    for (int i = 0; i < check_inds.size(); i++){
        int primary = check_inds[i];
        double m1 = masses[primary], r1 = norm(pos[primary], pos[primary]);
        for (int j = i + 1; j < check_inds.size(); j++){
            int secondary = check_inds[j];
            double m2 = masses[secondary], r2 = norm(pos[secondary], pos[secondary]);
            double R_mH = pow((m1 + m2) / (3 * masses[0]), 1./3) * (r1 + r2) / 2;
            double dist = dist_norm(pos[primary], pos[secondary]);
            if (dist < R_mH){
                std::cout << "Capture! " << R_mH << " " << dist << std::endl;
                for (int k = 0; k < check_inds.size(); k++){std::cout << check_inds[k] << std::endl;}
                captured.push_back(primary);
                std::sort(std::begin(captured), std::end(captured)); // apparently captured needs to be sorted for the check_inds line to work
                masses[secondary] = 0.95 * (m1 + m2); masses[primary] = 0;
                pos[secondary] = vec_scalar(vec_add(vec_scalar(pos[primary], m1), vec_scalar(pos[secondary], m2)), 1./(m1 + m2));
                vel[secondary] = vec_scalar(vec_add(vec_scalar(vel[primary], m1), vec_scalar(vel[secondary], m2)), 1./(m1 + m2));
                break;
            }
        }
    }
    for (int i = 0; i < captured.size(); i++){
        int index = captured[i];
        vel[index] = {0, 0, 0}; pos[index] = {0, 0, 0}; accel[i] = {0, 0, 0};
    }
}

void nbody_integrator(std::vector<std::vector<double>> &pos, std::vector<std::vector<double>> &vel, std::vector<std::vector<double>> &accel,
                        double dt, double Tmax, std::vector<double> masses, AGNDisk agn){
    int nt = int(Tmax / dt) + 1, N = masses.size();
    std::vector<int> captured(1, 0);
    std::vector<std::vector<std::vector<double>>> positions(N, std::vector<std::vector<double>> (3, std::vector<double> (nt)));
    std::vector<std::vector<std::vector<double>>> velocities(N, std::vector<std::vector<double>> (3, std::vector<double> (nt)));
    std::vector<double> times(nt);
    for (int n = 0; n < N; n++){
        for (int j = 0; j < 3; j++){
            positions[n][j][0] = pos[n][j];
            positions[n][j][0] = vel[n][j];
        }
    }
    times[0] = 0;
    for (int n = 0; n < N; n++){
        for (int j = 0; j < 3; j++){
            positions[n][j][0] = pos[n][j];
            velocities[n][j][0] = vel[n][j];
        }
    }
    for (int t = 1; t < nt; t++){
        nbody_timestep(pos, vel, accel, dt, masses, agn, captured);
        for (int n = 0; n < N; n++){
            for (int j = 0; j < 3; j++){
                positions[n][j][t] = pos[n][j];
                velocities[n][j][t] = vel[n][j];
            }
        }
        times[t] = t * dt;
    }
    save_1d_vec(times, "times.txt");
    save_3d_vec(positions, "positions.txt");
    save_3d_vec(velocities, "velocities.txt");
}

void init_conds(std::vector<double> &masses, std::vector<std::vector<double>> &pos, std::vector<std::vector<double>> &vel){
    int N = masses.size();
    double theta, dist, R, des_vel, angle, xprop, yprop, zprop, mult;
    
    // uniformly (and randomly) distribute points in the unit disk
    for (int i = 1; i < N; i++){
        theta = uniform(0., 2*M_PI);
        dist = uniform(0.5*0.5*0.5, 1.);
        R = pow(dist, 1./3);
        pos[i][0] = R * cos(theta); // x
        pos[i][1] = R * sin(theta); // y
        pos[i][2] = 0; // z
        des_vel = sqrt(1./R);
        angle = atan2(pos[i][1], pos[i][0]); // arctan2(y, x)
        xprop = sin(angle);
        yprop = -cos(angle);
        zprop = 0;
        mult = sqrt(des_vel*des_vel / (xprop*xprop + yprop*yprop + zprop*zprop));
        xprop *= mult; yprop *= mult;
        vel[i][0] = xprop; vel[i][1] = yprop; vel[i][2] = zprop;
    }
    // insert the SMBH into the start of the array
    for (int i = 0; i < 3; i++){pos[0][i] = 0;} // set SMBH position to 0
    // now convert from solar masses to nbody masses
    double tot_mass;
    for (int i = 0; i < N; i++){tot_mass += masses[i];} // calculate total mass of system
    for (int i = 0; i < N; i++){masses[i] = masses[i] / tot_mass;}  // convert each mass
}




int main(){
    time_t t1 = time(NULL);
    double Tmax = 1000, dt = 0.01;
    int NBH = 10;
    double SMBHMass = 1e8, r_s = 2 * G_pc * SMBHMass / 9e10;  // 2GM / c^2     units of pc
    double Nr_s = 1e3;      // number of schwarzschild radii to initialise the sim with respect to
    double lenscale = Nr_s * r_s;
    std::vector<double> BHMasses(NBH + 1, 10);
    BHMasses[0] = SMBHMass;
    AGNDisk agn (SMBHMass, lenscale);
    std::vector<std::vector<double>> pos(NBH + 1, std::vector<double> (3)), vel(NBH + 1, std::vector<double> (3));
    init_conds(BHMasses, pos, vel);
    std::vector<std::vector<double>> accel = nbody_accel(BHMasses, pos);
    nbody_integrator(pos, vel, accel, dt, Tmax, BHMasses, agn);
    std::cout << time(NULL) - t1 << std::endl;
    return 0;
}