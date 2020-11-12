/*
 * Code to move a charged particle in a EM field
 * with Boris method.
 *
 * Authors: S. Adhikari and J. Zrake (Clemson University)
 */
#include <stdio.h>
#include <math.h>
#define num_steps_per_orbit 500

static const double tfinal = 50.0;
static const double charge_of_particle = 1.0;
static const double mass_of_particle   = 1.0;
static const double speed_of_light     = 1.0;
static const double initial_velocity[3] = {0.0, 0.866, 0.0};
static const double initial_position[3] = {0.0, 0.0, 0.0};
static const double electric_field[3]   = {0.0, 0.0, 0.0};
static const double magnetic_field[3]   = {0.0, 0.0, 1.0};




// ============================================================================
void cross(double a[3], double b[3], double c[3])
{
    c[0] = a[1] * b[2] - b[1] * a[2];
    c[1] = a[2] * b[0] - b[2] * a[0];
    c[2] = a[0] * b[1] - b[0] * a[1];
}

double dot(double a[3], double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void sum(double a[3], double b[3], double c[3])
{
    for (int i = 0; i < 3; ++i)
    {
        c[i] = a[i] + b[i];
    }
}

void boris_velocity_update(double B[3], double E[3], double vnmh[3], double q_over_m, double dt, double vnph[3])
{
    // input v_(n-1/2) and output v_(n+1/2)
    double vm[3]; // v_minus
    double t [3];
    double s [3];
    double vv[3]; // v_prime
    double vp[3]; // v_plus
    double vm_cross_t[3];
    double vv_cross_s[3];
    double t_squared;

    for (int i = 0; i < 3; ++i)
    {
        vm[i] = vnmh[i] + q_over_m * E[i] * 0.5 * dt;
        t[i] = q_over_m * B[i] * 0.5 * dt;
    }

    t_squared = dot(t, t);

    for (int i = 0; i < 3; ++i)
    {
        s[i] = 2.0 * t[i] / (1 + t_squared);
    }

    cross(vm, t, vm_cross_t);
    sum(vm, vm_cross_t, vv);
    cross(vv, s, vv_cross_s);
    sum(vm, vv_cross_s, vp);

    for (int i = 0; i < 3; ++i)
    {
        vnph[i] = vp[i] + q_over_m * E[i] * 0.5 * dt;
    }
}

void push_particle(double x[3], double vnmh[3], double vnph[3], double dt)
{
    for (int i = 0; i < 3; ++i)
    {
        x[i] += (vnmh[i] + vnph[i]) * 0.5 * dt;
    }
}

void reassign_vel(double vnph[3], double vnmh[3])
{
    for (int i = 0; i < 3; ++i)
    {
        vnmh[i] = vnph[i];
    }
}

void err_energy(double energy1, double vnph[3], double partm, double c, double energy_diff)
{
    double gamma2 = 1 / sqrt(1 - dot(vnph, vnph));
    double energy2 = (gamma2 - 1) * partm * c * c;

    energy_diff = fabs(energy1 - energy2) / energy1 * 100;
}




// ============================================================================
int main()
{
    double B[3];
    double E[3];
    double vnmh[3];
    double vnph[3];
    double x[3];
    double xin[3];
    double t = 0.0;
    double partq = charge_of_particle;
    double partm = mass_of_particle;
    double c = speed_of_light;
    double energy_diff;
    double phase; //initial phase
    double real_phase;
    double phase_diff;
    

    // initialize parameters

    for (int i = 0; i < 3; ++i)
    {
        B[i] = magnetic_field[i];
        E[i] = electric_field[i];
        vnmh[i] = initial_velocity[i];
        x[i] = initial_position[i];
    }
    double gamma1 = 1 / sqrt(1 - dot(vnmh, vnmh));
    double energy1 = (gamma1 - 1) * partm * c * c;
    double B_mag = sqrt(dot(B, B));
    double omega_B = B_mag * partq / (gamma1 * partm);
    double larmor_period = 2.0 * M_PI / omega_B;
    double dt = larmor_period / num_steps_per_orbit;

    FILE* outfile = fopen("position_out.dat", "w");

    while (t < tfinal)
    {
        //real_phase = fmod(1.0 * M_PI + omega_B * 0.5 * t, M_PI);
        //phase = atan(x[1] - x[0]);
        real_phase = fmod(omega_B * (t + dt), M_PI);
        boris_velocity_update(B, E, vnmh, partq / partm, dt, vnph);
        err_energy(energy1, vnph, partm, c, energy_diff);
        push_particle(x, vnmh, vnph, dt);
        reassign_vel(vnph, vnmh);
        phase = 0.5 * M_PI - atan(x[1] / x[0]);// phase = arctan(y/x)
        phase_diff = real_phase - phase;
        fprintf(outfile, "%f %f %f %f %f\n", x[0], x[1], x[2], t, phase_diff);//Writes the x, y, and z coordinates as well as the %error in energy and phase.
        t += dt;
    }

    fclose(outfile);
    return 0;
}