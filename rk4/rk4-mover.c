/*
 * Code to move a charged particle in a EM field
 * with RK4 method.
 *
 * Authors: S. Adhikari and J. Zrake (Clemson University)
 */
#include <stdio.h>
#include <math.h>
#define num_steps_per_orbit 1000

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

void func(const double t,const double u_n[4], const double E[4], const double B[4], double force[4])
{
    double e = charge_of_particle;
    double m = mass_of_particle;
    double q_to_m = e / m;
    double force_e[4], force_m[4];
    double v_cross_B[4];

    cross(u_n, B, v_cross_B);
    for (int d = 1; d <= 3; ++d)
      {
          force_e[d] = q_to_m * E[d];
          force_m[d] = q_to_m * v_cross_B[d];
      }
    sum(force_e, force_m, force);
}


void push_particle(const double E[4], const double B[4], double x[4], double u_n[4], double dt)
{
    double k1[4], k2[4], k3[4], k4[4], ktemp[4], unpk[4], u_np1[4];
    
    func(x[0], u_n, E, B, k1); //k1 calculation
    for (int d = 1; d <= 3; ++d)
    {
        ktemp[d] = k1[d] * dt * 0.5;
    }
    
    sum(u_n, ktemp, unpk); //redefining velocity and time for k2
    x[0] += dt * 0.5;

    func(x[0], unpk, E, B, k2); //k2 calculation
    for (int d = 1; d <= 3; ++d)
    {
        ktemp[d] = k2[d] * dt * 0.5;
    }

    sum(u_n, ktemp, unpk); // redefining velocity for k3

    func(x[0], unpk, E, B, k3); //k3 calculation
    for (int d = 1; d <= 3; ++d)
    {
        ktemp[d] = k3[d] * dt;
    }

    sum(u_n, ktemp, unpk); //redefining velocity and time for k4
    x[0] += dt * 0.5;

    func(x[0], unpk, E, B, k4);//k4 calculation

    for (int d = 1; d <= 3; ++d) //updating the position using the velocity at time t
        x[d] += dt * u_n[d];

    for (int d = 0; d <= 3; ++d) // updating the velocity for time t + dt
    {
        u_n[d] += dt * (k1[d] + 2.0 * k2[d] + 2.0 * k3[d] + k4[d]) / 6.0;
    }

    x[0] += dt;
}




// ============================================================================
struct Particle
{
    double x[4];   // position at t
    double u_n[4]; // velocity at t
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

    exact_particle_position(state.x[0] * dt, state.x);
    exact_particle_velocity(state.x[0] * dt, state.u_n);

    return state;
}

struct Particle advance(struct Particle state)
{
    double dt = larmor_period() / num_steps_per_orbit;

    push_particle(electric_field, magnetic_field, state.x, state.u_n, dt);

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
