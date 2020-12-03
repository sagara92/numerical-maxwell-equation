/*
 * Code to move a charged particle in a EM field
 * with RK4 method.
 *
 * Authors: S. Adhikari and J. Zrake (Clemson University)
 */
#include <stdio.h>
#include <math.h>
#define num_steps_per_orbit 100

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

double lorentz_factor(const double u[4]) // JJZ: ADDED THIS FUNCTION: u[0] SHOULD NEVER STORE ANYTHING,
                                         // LORENTZ FACTOR SHOULD BE DERIVED INSTEAD. AVOIDS THE NEED TO
                                         // SYNCHRONIZE LOGICALLY REDUNDANT DATA, A MAJOR SOURCE OF BUGS.
{
    return sqrt(1.0 + dot(u, u));
}

void time_derivative_u(const double E[4], const double B[4], const double x[4], const double u[4], double dudt[4])
{
    double e = charge_of_particle;
    double m = mass_of_particle;


    double q_to_m = e / m;
    double u_cross_B[4];

    cross(u, B, u_cross_B);

    for (int d = 1; d <= 3; ++d)
    {
        dudt[d] = q_to_m * (E[d] + u_cross_B[d]) / lorentz_factor(u);
    }
}

void time_derivative_x(const double E[4], const double B[4], const double x[4], const double u[4], double dxdt[4])
{
    dxdt[0] = 1.0; // JJZ: dt/dt = 1.0

    for (int d = 1; d <= 3; ++d)
    {
        dxdt[d] = u[d] / lorentz_factor(u);
    }
}

void push_particle(const double E[4], const double B[4], double x0[4], double u0[4], const double dt)
{
    double x1[4], x2[4], x3[4];
    double u1[4], u2[4], u3[4];
    double dxdt0[4], dxdt1[4], dxdt2[4], dxdt3[4];
    double dudt0[4], dudt1[4], dudt2[4], dudt3[4];


    // STEP 1
    // ----------------------------------
    time_derivative_x(E, B, x0, u0, dxdt0);
    time_derivative_u(E, B, x0, u0, dudt0);

    for (int d = 0; d <= 3; ++d) // JJZ: NOTE THE LOOP STARTS AT 0: NECESSARY TO UPDATE TIME
    {
        x1[d] = x0[d] + dxdt0[d] * dt * 0.5;
        u1[d] = u0[d] + dudt0[d] * dt * 0.5;
    }


    // STEP 2
    // ----------------------------------
    time_derivative_x(E, B, x1, u1, dxdt1);
    time_derivative_u(E, B, x1, u1, dudt1);

    for (int d = 0; d <= 3; ++d)
    {
        x2[d] = x0[d] + dxdt1[d] * dt * 0.5;
        u2[d] = u0[d] + dudt1[d] * dt * 0.5;
    }


    // STEP 3
    // ----------------------------------
    time_derivative_x(E, B, x2, u2, dxdt2);
    time_derivative_u(E, B, x2, u2, dudt2);

    for (int d = 0; d <= 3; ++d)
    {
        x3[d] = x0[d] + dxdt2[d] * dt * 1.0;
        u3[d] = u0[d] + dudt2[d] * dt * 1.0;
    }


    // STEP 4
    // ----------------------------------
    time_derivative_x(E, B, x3, u3, dxdt3);
    time_derivative_u(E, B, x3, u3, dudt3);


    // UPDATE POSITION AND GAMMA-BETA
    // ----------------------------------
    for (int d = 0; d <= 3; ++d)
    {
        x0[d] += (dxdt0[d] + dxdt1[d] * 2.0 + dxdt2[d] * 2.0 + dxdt3[d]) * dt / 6.0;
        u0[d] += (dudt0[d] + dudt1[d] * 2.0 + dudt2[d] * 2.0 + dudt3[d]) * dt / 6.0;
    }
}




// ============================================================================
struct Particle
{
    double x[4]; // particle position
    double u[4]; // particle four-velocity  JJZ -- PARTICLE STORES ITS GAMMA-BETA: CLEARER IN MY OPINION
};




// ============================================================================
double larmor_frequency()
{
    double e = charge_of_particle;
    double b = sqrt(dot(magnetic_field, magnetic_field));
    double m = mass_of_particle;
    double gamma = particle_lorentz_factor;
    //
    // JJZ -- CRITICAL BUG BELOW: ORDER OF OPERATIONS. THIS WAS PRESENT ALREADY
    // IN THE SNIPPET I SENT ON 11-12.
    //
    // return fabs(e * b / gamma * m);

    return fabs(e * b / gamma / m);
}

double particle_beta() // JJZ -- RENAMED THIS VELOCITY -> BETA
{
    double gamma = particle_lorentz_factor;
    return sqrt(1.0 - 1.0 / gamma / gamma);
}

double larmor_radius()
{
    return particle_beta() / larmor_frequency();
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

void exact_particle_gamma_beta(double t, double gamma_beta[4])
{
    double g = particle_lorentz_factor;
    gamma_beta[0] = 0.0; // JJZ: YOU MIGHT EXPECT u[0] TO BE GAMMA, BUT u[0] SHOULD NEVER BE
                         // ACCESSED (SEE ABOVE), SO SETTING IT TO ZERO MAKES MORE SENSE.
    gamma_beta[1] = g * particle_beta() * sin(larmor_frequency() * t) * -1.0;
    gamma_beta[2] = g * particle_beta() * cos(larmor_frequency() * t) * +1.0;
    gamma_beta[3] = 0.0;
}




// ============================================================================
struct Particle initial_state()
{
    struct Particle state;

    exact_particle_position  (state.x[0], state.x);
    exact_particle_gamma_beta(state.x[0], state.u);

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
    double e_exact = mass_of_particle * dot(state.u, state.u) / particle_lorentz_factor;
    

    FILE* outfile = fopen("updated-solution.dat", "w");

    while (state.x[0] < tfinal)
    {
        state = advance(state);
        double vel2 = dot(state.u, state.u);
        double energy = mass_of_particle * vel2 / lorentz_factor(state.u);
        double x_exact[4];
        exact_particle_position(state.x[0], x_exact);
        double exact_phase = atan2(x_exact[2], x_exact[1]);
        double phase       = atan2(state.x[2], state.x[1]);
        //fprintf(outfile, "%+.09e %+.09e %+.09e\n", state.x[0], exact_phase, phase);
        fprintf(outfile, "%+.09e %+.09e %+.09e\n", state.x[0], x_exact[1], state.x[1]);
    }
    fclose(outfile);
    return 0;
}