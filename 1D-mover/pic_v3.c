/*
 * Code to move a charged particle in 1d electromagnetic field. Uses staggered-in-space and staggered-in-time update scheme.
 *
 * Authors: S. Adhikari and J. Zrake (Clemson University)
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NUM_CELLS 100
#define wavenumber_k ((2.0 * M_PI) / (domain_x1 - domain_x0))
#define grid_spacing ((domain_x1 - domain_x0) / NUM_CELLS)
#define num_steps_per_orbit 200

static const double speed_of_light = 1.0;
static const double domain_x0 = 0.0;
static const double domain_x1 = 100.0;
static const double cfl_parameter = 0.1;

static const double tfinal                  =  40.0 * M_PI;
static const double charge_of_particle      =  1.0;
static const double mass_of_particle        =  1.0;
static const double particle_lorentz_factor =  2.0;
static const double electric_field[4] = {0.0, 0.0, 0.0, 0.0};
static const double magnetic_field[4] = {0.0, 0.0, 0.0, 1.0};





// ============================================================================
struct Particle
{
    double x[4]; // position at t
    double u[4]; // velocity at t - dt / 2 (four-velocity gamma-beta)           //Step 2//
};

// ============================================================================
void cross(const double a[4], const double b[4], double c[4])
{
    c[1] = a[2] * b[3] - b[2] * a[3];
    c[2] = a[3] * b[1] - b[3] * a[1];
    c[3] = a[1] * b[2] - b[1] * a[2];
}

// ============================================================================
double dot(const double a[4], const double b[4])
{
    return a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

// ============================================================================
void sum(const double a[4], const double b[4], double c[4])
{
    for (int i = 0; i < 4; ++i)
    {
        c[i] = a[i] + b[i];
    }
}

// ============================================================================
void exact_solution(double* E, double* B, double t, double dt)
{
    double c  = speed_of_light;
    double k  = wavenumber_k;
    double w  = k * c;
    double dx = grid_spacing;
    double E0 = sqrt(dot(electric_field, electric_field));
    double B0 = sqrt(dot(magnetic_field, magnetic_field));

    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double x = domain_x0 + i * dx;
        E[i] = E0 * sin(k * (x + 0.0 * dx) + w * (t + 0.0 * dt));               //Step 1 & 3//
        B[i] = B0 * sin(k * (x + 0.5 * dx) + w * (t + 0.5 * dt));
    }
}

// ============================================================================
void update_solution(double* E, double* B, double dt)
{
    // E := E^{n}
    // B := B^{n+1/2}

    double c  = speed_of_light;
    double dx = grid_spacing;
    double dEdx[NUM_CELLS];
    double dBdx[NUM_CELLS];

    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double Bplus  = B[i];
        double Bminus = B[(i - 1 + NUM_CELLS) % NUM_CELLS];
        dBdx[i] = (Bplus - Bminus) / dx;
        E[i] += c * dBdx[i] * dt;                                               //Step 4//
    } // E := E^{n+1}

    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double Eplus  = E[(i + 1) % NUM_CELLS];
        double Eminus = E[i];
        //printf("Eplus--> %f\n", Eplus);

        dEdx[i] = (Eplus - Eminus) / dx;
        B[i] += c * dEdx[i] * dt;                                               //Step 5//
    } // B := B^{n+3/2}
}

// ============================================================================
void mag_int(double* E, double* B, double dt, double* Bnew)
{
    // B := B^{n+1/2}

    double c  = speed_of_light;
    double dx = grid_spacing;
    double dEdx[NUM_CELLS];

    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double Eplus  = E[(i + 1) % NUM_CELLS];
        double Eminus = E[i];

        dEdx[i] = (Eplus - Eminus) / dx;
        Bnew[i] = B[i] - 0.5 * c * dEdx[i] * dt;                                         //Step 6//
    } // B := B^{n}
}

// ============================================================================
void collocation_E_B(double* E, double* B, double x[4], double Enew[4], double Bnew[4])       //Step 7//
{
    double dx = grid_spacing;
    int i = (int) (x[1] / dx);
    int j = (int) ((x[1] / dx) - 0.5);
    double Eplus  = E[(i + 1) % NUM_CELLS];
    double Eminus = E[i];
    
    Enew[2] = Eminus + (x[1] - domain_x0 - (i * dx)) * (Eplus - Eminus) / dx;
    Enew[0] = 0.0;
    Enew[1] = 0.0;
    Enew[3] = 0.0;

    double Bminus = B[j];
    double Bplus  = B[(j + 1) % NUM_CELLS];

    Bnew[3] = Bminus + (x[1] - domain_x0 - (j * dx)) * (Bplus - Bminus) / dx;
    Bnew[0] = 0.0;
    Bnew[1] = 0.0;
    Bnew[2] = 0.0;
}


// ============================================================================
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

    for (int d = 1; d <= 3; ++d){
        u_nmh[d] = u[d];
    }

    for (int d = 1; d <= 3; ++d){
        u_minus[d] = u_nmh[d] + E[d] * 0.5 * h;
    }

    gamma_n = sqrt(1.0 + dot(u_minus, u_minus));

    for (int d = 1; d <= 3; ++d){
        t[d] = 0.5 * h * B[d] / gamma_n;
    }

    for (int d = 1; d <= 3; ++d){
        s[d] = 2.0 * t[d] / (1.0 + dot(t, t));
    }

    cross(u_minus, t, u_minus_cross_t);

    for (int d = 1; d <= 3; ++d){
        u_prime[d] = u_minus[d] + u_minus_cross_t[d];
    }

    cross(u_prime, s, u_prime_cross_s);

    for (int d = 1; d <= 3; ++d){
        u_plus[d] = u_minus[d] + u_prime_cross_s[d];
    }

    for (int d = 1; d <= 3; ++d){
        u_nph[d] = u_plus[d] + E[d] * 0.5 * h;
    }

    gamma_nph = sqrt( 1.0 + dot(u_nph, u_nph) );

    //u[0] = gamma_nph;

    for (int d = 1; d <= 3; ++d){                        /* Update the velocity */
        u[d] = u_nph[d];
    }

    x[0] += dt;

    for (int d = 1; d <= 3; ++d){                        /* Update the position */
        x[d] += dt * u_nph[d] / gamma_nph;
    }
}

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

struct Particle advance(struct Particle state, double E[4], double B[4])
{
    double dt = larmor_period() / num_steps_per_orbit;

    push_particle(E, B, state.x, state.u, dt);

    return state;
}



// ============================================================================
int main()
{
    double E[NUM_CELLS];
    double B[NUM_CELLS];
    double Bn[NUM_CELLS];
    double E1[4];
    double B1[4];
    double dt = larmor_period() / num_steps_per_orbit;
    struct Particle state = initial_state();

    exact_solution(E, B, state.x[0], dt);   
    FILE* outE = fopen("Electric.bin", "wb");
    FILE* outB = fopen("Magnetic.bin", "wb");
    FILE* outP = fopen("Particle_solution.bin", "wb");
    FILE* outfile = fopen("solution.dat", "w");

    while (state.x[0] < tfinal)
    {
        fprintf(outfile, "%f %f %f %f\n", state.x[0], state.x[1], state.x[2], state.x[3]);
        mag_int(E, B, dt, Bn);
        collocation_E_B(E, Bn, state.x, E1, B1);
        //for (int i = 0; i < 4; ++i)
        //{
        //    printf("%f, %f, %f\n", state.x[0], E1[i], B1[i]);
        //}
        size_t elements_written_E = fwrite(&E, sizeof(E), 1, outE);
        size_t elements_written_B = fwrite(&B, sizeof(B), 1, outB);
        size_t elements_written_P = fwrite(&state, sizeof(state), 1, outP);
        state = advance(state, E1, B1);
        update_solution(E, B, dt);
    }
    
    fclose(outE);
    fclose(outB);
    fclose(outP);
    fclose(outfile);
    return 0;
}