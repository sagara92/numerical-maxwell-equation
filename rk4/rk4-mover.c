/*
 * Code to move a charged particle in a EM field
 * with RK4 method.
 *
 * Authors: S. Adhikari and J. Zrake (Clemson University)
 */
#include <stdio.h>
#include <math.h>
#define num_steps_per_orbit 80

static const double speed_of_light          =  1.0;
static const double tfinal                  =  4.0 * M_PI;
static const double charge_of_particle      = -1.0;
static const double mass_of_particle        =  1.0;
static const double particle_lorentz_factor =  2.0;
static const double electric_field[4] = {0.0, 0.0, 0.0, 0.0};
static const double magnetic_field[4] = {0.0, 0.0, 0.0, 1.0};




// ============================================================================
void cross(const double a[4], const double b[4], double c[4])
{
    c[1] = a[2] * b[3] - b[2] * a[3];
    c[2] = a[3] * b[1] - b[3] * a[1];
    c[3] = a[1] * b[2] - b[1] * a[2];
}

double dot(const double a[4], const double b[4])
{
    return a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

void sum(const double a[4], const double b[4], double c[4])
{
    for (int i = 1; i <= 3; ++i)
    {
        c[i] = a[i] + b[i];
    }
}

void time_derivative_v(const double E[4], const double B[4], const double x[4], const double v_n[4], double dvdt[4])
{
    double e = charge_of_particle;
    double m = mass_of_particle;
    double q_to_m = e / m / v_n[0];
    double v_cross_B[4];

    cross(v_n, B, v_cross_B);
    for (int d = 1; d <= 3; ++d)
      {
          dvdt[d] = q_to_m * (E[d] + v_cross_B[d]);
      }
}

void time_derivative_x(const double E[4], const double B[4], const double x[4], const double v_n[4], double dxdt[4])
{
    for (int d = 1; d <= 3; ++d)
     {
         dxdt[d] = v_n[d] / v_n[0];
     }
}


void push_particle(const double E[4], const double B[4], double x[4], double v_n[4], const double dt)
{
    double dvdt[4], dxdt[4], x1[4], v_n1[4];
    v_n1[0] = v_n[0];
    double dx1[4], dx2[4], dx3[4], dx4[4], dv1[4], dv2[4], dv3[4], dv4[4];
    time_derivative_x(E, B, x, v_n, dxdt);
    for (int d = 1; d <= 3; ++d)
    {
        dx1[d] = dxdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        x1[d] = x[d] + 0.5 * dx1[d];
    }
    time_derivative_v(E, B, x1, v_n, dvdt);
    for (int d = 1; d <= 3; ++d)
    {
        dv1[d] = dvdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        v_n1[d] = v_n[d] + 0.5 * dv1[d];
    }
    time_derivative_x(E, B, x1, v_n1, dxdt);
    for (int d = 1; d <= 3; ++d)
    {
        dx2[d] = dxdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        x1[d] = x[d] + 0.5 * dx2[d];
    }
    time_derivative_v(E, B, x1, v_n1, dvdt);
    for (int d = 1; d <= 3; ++d)
    {
        dv2[d] = dvdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        v_n1[d] = v_n[d] + 0.5 * dv2[d];
    }
    time_derivative_x(E, B, x1, v_n1, dxdt);
    for (int d = 1; d <= 3; ++d)
    {
        dx3[d] = dxdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        x1[d] = x[d] + 1.0 * dx3[d];
    }
    time_derivative_v(E, B, x1, v_n1, dvdt);
    for (int d = 1; d <= 3; ++d)
    {
        dv3[d] = dvdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        v_n1[d] = v_n[d] + 1.0 * dv3[d];
    }
    time_derivative_x(E, B, x1, v_n1, dxdt);
    for (int d = 1; d <= 3; ++d)
    {
        dx4[d] = dxdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        x1[d] = x[d] + 1.0 * dx4[d];
    }
    time_derivative_v(E, B, x1, v_n1, dvdt);
    for (int d = 1; d <= 3; ++d)
    {
        dv4[d] = dvdt[d] * dt;
    }
    for (int d = 1; d <= 3; ++d)
    {
        x[d]   += (dx1[d] + dx2[d] * 2.0 + dx3[d] * 2.0 + dx4[d]) / 6.0;
        v_n[d] += (dv1[d] + dv2[d] * 2.0 + dv3[d] * 2.0 + dv4[d]) / 6.0;
    }

    x[0] += dt;
}




// ============================================================================
struct Particle
{
    double x[4];   // position at t
    double v_n[4]; // velocity at t
};




// ============================================================================
double larmor_frequency()
{
    double e = charge_of_particle;
    double b = sqrt(dot(magnetic_field, magnetic_field));
    double m = mass_of_particle;
    double gamma = particle_lorentz_factor;
    return fabs(e * b / gamma * m);
}

double particle_velocity()
{
    double gamma = particle_lorentz_factor;
    return sqrt(1.0 - 1.0 / gamma / gamma);
}

double larmor_radius()
{
    return particle_velocity() / larmor_frequency();
}

double larmor_period()
{
    return 2 * M_PI / larmor_frequency();
}




// ============================================================================
void exact_particle_position(double t, double position[4])
{
    position[0] = t;
    position[1] = larmor_radius() * cos(larmor_frequency() * t);
    position[2] = larmor_radius() * sin(larmor_frequency() * t);
    position[3] = 0.0;
}

void exact_particle_velocity(double t, double velocity[4])
{
    double g = particle_lorentz_factor;
    velocity[0] = g;
    velocity[1] = g * particle_velocity() * sin(larmor_frequency() * t) * -1.0;
    velocity[2] = g * particle_velocity() * cos(larmor_frequency() * t) * +1.0;
    velocity[3] = 0.0;
}



// ============================================================================
struct Particle initial_state()
{
    struct Particle state;

    double dt = larmor_period() / num_steps_per_orbit;

    exact_particle_position(state.x[0], state.x);
    exact_particle_velocity(state.x[0], state.v_n);

    return state;
}

struct Particle advance(struct Particle state)
{
    double dt = larmor_period() / num_steps_per_orbit;

    push_particle(electric_field, magnetic_field, state.x, state.v_n, dt);

    return state;
}




// ============================================================================
int main()
{

    struct Particle state = initial_state();

    FILE* outfile = fopen("solution.dat", "w");
    while (state.x[0] < tfinal)
    {
        state = advance(state);
        double vel2 = dot(state.v_n, state.v_n);
        double energy = 0.5 * mass_of_particle * vel2;
        double e_exact = particle_lorentz_factor * mass_of_particle * particle_velocity() * particle_velocity();
        double x_exact[4];
        exact_particle_position(state.x[0], x_exact);
        double exact_phase = atan2(x_exact[2], x_exact[1]);
        double phase       = atan2(state.x[2], state.x[1]);
        fprintf(outfile, "%f %f %f\n", state.x[0], exact_phase, phase);
        //fprintf(outfile, "%f %f %f %f %f\n", state.x[0], x_exact[1], x_exact[2], state.x[1], state.x[2]);
    }
    fclose(outfile);

    return 0;
}