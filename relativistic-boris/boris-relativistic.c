/*
 * Code to move a charged particle in a EM field
 * with Boris method.
 *
 * Authors: S. Adhikari and J. Zrake (Clemson University)
 */
#include <stdio.h>
#include <math.h>
#define num_steps_per_orbit 100

static const double speed_of_light          =  1.0;
static const double tfinal                  =  80.0 * M_PI;
static const double charge_of_particle      = -1.0;
static const double mass_of_particle        =  1.0;
static const double particle_lorentz_factor =  3.0;
static const double electric_field[4] = {0.0, 0.0, 0.3, 0.0};
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
    for (int i = 0; i < 4; ++i)
    {
        c[i] = a[i] + b[i];
    }
}




/* --------------------------------------------------------
* Relativistic particle mover, attributed to Boris (1970).
*
* The algorithm is described in Birdsall & Langdon (1991)
* Plasma Physics via Computer Simulation, Section 15-4.
* It is based on a rotation of the (spatial) 4-velocity
* components around an axis parallel to the magnetic field
* through an angle
*
* theta = -2 arctan( (e B dt) / (2 gamma m c) )
*
* The method uses the time-centered 4-velocity, such that
* p->E and p->B are known on integral time steps t^{n},
* while p->u is computed for time steps t^{n+1/2}. Note
* that the positions are known at integral time steps:
*
* x^{n+1} = x^{n} + v^{n+1/2} * dt
*
* Author: Jonathan Zrake, zrake@nyu.edu
*
* --------------------------------------------------------
*/
void push_particle(const double E[4], const double B[4], double x[4], double u[4], double dt)
{
    double t[4], s[4];
    double u_nph[4], u_nmh[4];
    double u_plus[4], u_minus[4], u_prime[4];
    double u_prime_cross_s[4], u_minus_cross_t[4];
    double gamma_n, gamma_nph;

    double e = charge_of_particle;
    double m = mass_of_particle;
    double h = e / m * dt;

    for (int d = 1; d <= 3; ++d)
        u_nmh[d] = u[d];

    for (int d = 1; d <= 3; ++d)                        /* Step 1 */
        u_minus[d] = u_nmh[d] + E[d] * 0.5 * h;

    gamma_n = sqrt(1.0 + dot(u_minus, u_minus));       /* Step 2 */

    for (int d = 1; d <= 3; ++d)                        /* Step 3 */
        t[d] = 0.5 * h * B[d] / gamma_n;

    for (int d = 1; d <= 3; ++d)                        /* Step 4 */
        s[d] = 2.0 * t[d] / (1.0 + dot(t, t));

    cross(u_minus, t, u_minus_cross_t);

    for (int d = 1; d <= 3; ++d)                        /* Step 5 */
        u_prime[d] = u_minus[d] + u_minus_cross_t[d];

    cross(u_prime, s, u_prime_cross_s);

    for (int d = 1; d <= 3; ++d)                        /* Step 6 */
        u_plus[d] = u_minus[d] + u_prime_cross_s[d];

    for (int d = 1; d <= 3; ++d)                        /* Step 7 */
        u_nph[d] = u_plus[d] + E[d] * 0.5 * h;

    gamma_nph = sqrt( 1.0 + dot(u_nph, u_nph) );

    u[0] = gamma_nph;

    for (int d = 1; d <= 3; ++d)                        /* Update the velocity */
        u[d] = u_nph[d];

    x[0] += dt;

    for (int d = 1; d <= 3; ++d)                        /* Update the position */
        x[d] += dt * u_nph[d] / gamma_nph;
}




// ============================================================================
struct Particle
{
    double x[4]; // position at t
    double u[4]; // velocity at t - dt / 2 (four-velocity gamma-beta)
};




// ============================================================================
double larmor_frequency()
{
    double e = charge_of_particle;
    double b = sqrt(dot(magnetic_field, magnetic_field));
    double m = mass_of_particle;
    double gamma = particle_lorentz_factor;
    return fabs(e * b / gamma / m);
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
    velocity[0] = 0.0;
    velocity[1] = g * particle_velocity() * sin(larmor_frequency() * t) * -1.0;
    velocity[2] = g * particle_velocity() * cos(larmor_frequency() * t) * +1.0;
    velocity[3] = 0.0;
}




// ============================================================================
struct Particle initial_state()
{
    struct Particle state;

    double dt = larmor_period() / num_steps_per_orbit;

    exact_particle_position(state.x[0] + 0.0 * dt, state.x);
    exact_particle_velocity(state.x[0] - 0.5 * dt, state.u);

    return state;
}

struct Particle advance(struct Particle state)
{
    double dt = larmor_period() / num_steps_per_orbit;

    push_particle(electric_field, magnetic_field, state.x, state.u, dt);

    return state;
}




// ============================================================================
int main()
{

    struct Particle state = initial_state();

    //printf("%f %f\n", state.x[1], state.x[2]);
    FILE* outfile = fopen("solution.dat", "w");

    while (state.x[0] < tfinal)
    {
        state = advance(state);

        double x_exact[4];
        exact_particle_position(state.x[0], x_exact);
        //printf("%f %f\n", state.x[1], state.x[2]);
        double exact_phase = atan2(x_exact[2], x_exact[1]);
        double phase       = atan2(state.x[2], state.x[1]);
        //fprintf(outfile, "%+.09e %+.09e %+.09e\n", state.x[0], x_exact[1], state.x[1]);
        fprintf(outfile, "%+.09e %+.09e %+.09e %+.09e %+.09e\n", state.x[0], x_exact[1], x_exact[2], state.x[1], state.x[2]);

    }
    fclose(outfile);

    return 0;
}