#include <stdio.h>
#include <math.h>
#define NUM_CELLS 1000




int main()
{
    // Declare variables
    // -------------------------------------------------------------
    double E[NUM_CELLS];
    double B[NUM_CELLS];

    double x0 = 0.0;
    double x1 = 1.0;
    double dx = (x1 - x0) / (NUM_CELLS - 1);

    double E0 = 1.0;
    double B0 = 1.0;
    double c = 1.0;
    double k = 2 * M_PI / (x1 - x0);
    double w = k * c;
    double t = 0.0;
    double tfinal = 0.2;
    double dt = 0.3 * dx / c;
    int iteration = 0;


    // Initial condition
    // ------------------------------------------------------------------------
    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double x = x0 + i * dx;

        E[i] = E0 * sin(k * x - w * t);
        B[i] = B0 * sin(k * x - w * t);
    }


    // Upate the solution in time
    // ------------------------------------------------------------------------
    while (t < tfinal)
    {
        // E -> E + dE/dt * dt
        // B -> B + dB/dt * dt
        // 
        // dE/dt = -c dB/dx
        // dB/dt = -c dE/dx
        
        double dEdx[NUM_CELLS];
        double dBdx[NUM_CELLS];

        for (int i = 0; i < NUM_CELLS; ++i)
        {
            double Em1 = i >             0 ? E[i - 1] : E[NUM_CELLS - 2];
            double Ep1 = i < NUM_CELLS - 1 ? E[i + 1] : E[1];
            double Bm1 = i >             0 ? B[i - 1] : B[NUM_CELLS - 2];
            double Bp1 = i < NUM_CELLS - 1 ? B[i + 1] : B[1];

            dEdx[i] = (Ep1 - Em1) / (2.0 * dx);
            dBdx[i] = (Bp1 - Bm1) / (2.0 * dx);
        }

        for (int i = 0; i < NUM_CELLS; ++i)
        {
            E[i] -= c * dBdx[i] * dt;
            B[i] -= c * dEdx[i] * dt;
        }

        t += dt;
        iteration += 1;

        printf("[%04d] t=%0.3f\n", iteration, t);
    }


    // Output to a file
    // ------------------------------------------------------------------------
    FILE* outfile = fopen("em-wave.dat", "w");

    for (int i = 0; i < NUM_CELLS; ++i)
    {
        double x = x0 + i * dx;
        fprintf(outfile, "%.12lf %+.12lf %+.12lf\n", x, E[i], B[i]);
    }    
    fclose(outfile);

    return 0;
}