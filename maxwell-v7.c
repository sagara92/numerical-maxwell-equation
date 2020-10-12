/*
 * Code to solve the 1d electromagnetic wave equation for a plane-polarized
 * wave. Uses staggered-in-space and staggered-in-time update scheme.
 *
 * Authors: S. Adhikari and J. Zrake (Clemson University)
 */
#include <stdio.h>
#include <math.h>
#define NUM_CELLS 2500
#define wavenumber_k ((2.0 * M_PI) / (domain_x1 - domain_x0))
#define grid_spacing ((domain_x1 - domain_x0) / NUM_CELLS)

static const double speed_of_light = 1.0;
static const double domain_x0 = 0.0;
static const double domain_x1 = 1.0;
static const double wave_field_amplitude = 1.0;
static const double cfl_parameter = 0.4;
static const double tfinal = 1.0;




// ============================================================================
void exact_solution(double* E, double* B, double t, double dt)
{
    double c  = speed_of_light;
    double k  = wavenumber_k;
    double w  = k * c;
    double dx = grid_spacing;
    double E0 = wave_field_amplitude;
    double B0 = wave_field_amplitude;

    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double x = domain_x0 + i * dx;
        E[i] = E0 * sin(k * (x + 0.0 * dx) + w * (t + 0.0 * dt));
        B[i] = B0 * sin(k * (x + 0.5 * dx) + w * (t + 0.5 * dt));
    }
}

void update_solution(double* E, double* B, double t, double dt)
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
        E[i] += c * dBdx[i] * dt;
    } // E := E^{n+1}

    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double Eplus  = E[(i + 1) % NUM_CELLS];
        double Eminus = E[i];

        dEdx[i] = (Eplus - Eminus) / dx;
        B[i] += c * dEdx[i] * dt;
    } // B := B^{n+3/2}
}

void output_solution(double* E, double* B, double t, double dt, int iteration)
{
    if (iteration % 200 == 0)
    {
        double dx = grid_spacing;
        char fname[1024];

        snprintf(fname, 1024, "output.%04d.dat", iteration / 200);
        FILE* outfile = fopen(fname, "w");

        double E_exact[NUM_CELLS];
        double B_exact[NUM_CELLS];

        exact_solution(E_exact, B_exact, t, dt);

        for (int i = 0; i < NUM_CELLS; ++i)
        {
            double xE = domain_x0 + (i + 0.0) * dx;
            double xB = domain_x0 + (i + 0.5) * dx;
            fprintf(outfile, "%f %f %f %f %f %f\n", xE, xB, E[i], B[i], E_exact[i], B_exact[i]);
        }
        fclose(outfile);
        printf("output to %s\n", fname);
    }
 }   

void print_message(int iteration, double t)
{
    printf("[%04d] t=%0.3f\n", iteration, t);
}




// ============================================================================
int main()
{
    double E[NUM_CELLS];
    double B[NUM_CELLS];
    double t = 0.0;
    double dt = cfl_parameter * grid_spacing / speed_of_light;
    int iteration = 0;

    exact_solution(E, B, t, dt);

    while (t < tfinal)
    {
        output_solution(E, B, t, dt, iteration);
        print_message(iteration, t);
        update_solution(E, B, t, dt);

        t += dt;
        iteration += 1;
    }
    return 0;
}