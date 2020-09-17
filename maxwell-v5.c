#include<stdio.h>
#include <math.h>
#define NUM_CELLS 1000



int main()
{
	double x, dx, xi, xf, E0, B0, k, ti, w, c, tf, dt, t;
	double E[NUM_CELLS], B[NUM_CELLS];
	E0 = 5.0;
	B0 = 1.0;
	xi = 0.0;
	t  = 0.0;
	xf = 1.0;
	c  = 1.0;
	ti = 0.0;
	tf = 1.1;/*1 second evolution*/
	k  = 2.0 * M_PI /(xf - xi);
	w  = k * c;
	dx = (xf - xi)/(NUM_CELLS - 1);
	dt = 0.04 * (dx / c);

	 /*Initial condition at t = 0*/
	
	for (int i = 0; i < NUM_CELLS; i+=2)
	{
		x = xi + i*dx;
		E[i] = E0 * sin (k*x - w*t);
		
	}
	for (int i = 1; i < NUM_CELLS; i+=2)
	{
		x = xi + i*dx;
	
		B[i] = B0 * sin (k*x - w*t);
	}
	/*Time evolution of the field and position*/

	while(t < tf)
	{
		double dBdx[NUM_CELLS];
		
		double Eave[NUM_CELLS];
		
		for (int i = 0; i < NUM_CELLS; i+=2)
		{
			double Eminus = i > 0 			? E[i-2]:E[NUM_CELLS-2];
			double Eplus  = i < NUM_CELLS-2 ? E[i+2]:E[2];
			double Bminus = i > 0 			? B[i-2]:B[NUM_CELLS-2];
			double Bplus  = i < NUM_CELLS-2 ? B[i+2]:B[2];

			
			dBdx[i] = (Bplus - Bminus)/(2.0 *dx);

			Eave[i] = (Eplus + Eminus)/2.0;}
		for (int i = 0; i < NUM_CELLS; i+=2)
		{
			E[i] = Eave[i] - c * dBdx[i] * dt; /*Lax-Friedrichs Method*/
		}

		double dEdx[NUM_CELLS];
		
		
		double Bave[NUM_CELLS];
		for (int i = 1; i < NUM_CELLS; i+=2)
		{
			double Eminus = i > 0 			? E[i-2]:E[NUM_CELLS-3];
			double Eplus  = i < NUM_CELLS-1 ? E[i+2]:E[3];
			double Bminus = i > 0 			? B[i-2]:B[NUM_CELLS-3];
			double Bplus  = i < NUM_CELLS-1 ? B[i+2]:B[3];

			dEdx[i] = (Eplus - Eminus)/(2.0 *dx);
			
			Bave[i] = (Bplus + Bminus)/2.0;
			
		}
		for (int i = 1; i < NUM_CELLS; i+=2)
		{
			
			B[i] = Bave[i] - c * dEdx[i] * dt;
		}
		t += dt;
	}
	printf("Code ran well for %.2f second evolution\n", tf);
	// Output to a file
    // ------------------------------------------------------------------------
    FILE* outfile = fopen("em-wave-v5.dat", "w");

    for (int i = 1; i < NUM_CELLS-1; ++i)
    {
        double x = xi + 2.0*i * dx;
        double y = xi + (2.0*i +1) * dx;
        fprintf(outfile, "%f %f %f %f\n", x, E[i], y, B[i]);
    }
    
    fclose(outfile);
	return 0;
}