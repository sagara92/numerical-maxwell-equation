#include<stdio.h>
#include <math.h>
#define NUM_CELLS 2560
/*Yee Lattice and Leapfrogging*/


int main()
{	
	double x, dx, xi, y, yi, xf, E0, B0, k, ti, w, c, tf, dt, t, L2_E, L2_B;
	double E[NUM_CELLS], B[NUM_CELLS];
	double Diff_E[NUM_CELLS];
	double Diff_B[NUM_CELLS];
	//double L2_E[NUM_CELLS];
	//double L2_B[NUM_CELLS];
	double B_r[NUM_CELLS];
	double E_r[NUM_CELLS];
	E0 = 5.0;
	B0 = 1.0;
	xi = 0.0;
	yi = 0.0;
	t  = 0.0;
	xf = 1.0;
	c  = 1.0;
	ti = 0.0;
	tf = 2.0;	/* Time of evolution in seconds*/
	k  = 2.0 * M_PI /(xf - xi);
	w  = k * c;
	dx = (xf - xi)/(NUM_CELLS - 1);
	dt = 0.04 * (dx / c);



	/*Initial condition at t = 0*/
	

	for (int i = 0; i < NUM_CELLS; ++i)
	{
		x = xi + i*dx;
		y = yi + (i+0.5)*dx;
		E[i] = E0 * sin (k*x - w*t);
		B[i] = B0 * sin (k*y - w*t);
	}


	/*Time evolution of the field*/

	while(t < tf)
	{
		double dBdx[NUM_CELLS];
		double Eave[NUM_CELLS];
		double dEdx[NUM_CELLS];
		double Bave[NUM_CELLS];

		for (int i = 0; i < NUM_CELLS; ++i) 
		{
			double Eminus = i > 		  0	? E[i-1] : E[NUM_CELLS-2];
			double Eplus  = i < NUM_CELLS-1 ? E[i+1] : E[1];
			double Bminus = i >  		  1 ? B[i-1] : B[NUM_CELLS-2];
			double Bplus  = i < NUM_CELLS-1 ? B[i+1] : B[1];

			dBdx[i] = (Bplus - Bminus)/(2.0 *dx);
			Eave[i] = (0.9*Eplus + 0.1*Eminus);/*Weighted average*/

			E[i] = Eave[i] - c * dBdx[i] * dt; /*Lax-Friedrichs Method*/
		
			dEdx[i] = (Eplus - Eminus)/(2.0 *dx);
			Bave[i] = (0.9*Bplus + 0.1*Bminus);/*Weighted average*/

			B[i] = Bave[i] - c * dEdx[i] * dt; /*Lax-Friedrichs Method*/

		}
		t += dt;			
	}
	double E2=0.0;
	double B2=0.0;
	for (int i = 0; i < NUM_CELLS; ++i)
	{
		/*Target value of E and B*/
		x = xi + i*dx;
		y = yi + (i+0.5)*dx;
		E_r[i] = E0 * sin (k*x - w*dt); /*Real value of E*/
		B_r[i] = B0 * sin (k*y - w*dt);	/*Real value of B*/	

		/*Calculating L2 errors in each position*/			
		Diff_E[i] = E_r[i] - E[i];		/*Difference in the calculated value and real value*/	
		E2 += (Diff_E[i]*Diff_E[i])*dx;
					
		Diff_B[i] = B_r[i] - B[i];
		B2 += (Diff_B[i]*Diff_B[i])*dx;
		
	}
	L2_E = sqrt(E2)/(xf-xi);
	L2_B = sqrt(B2)/(xf-xi);
	printf("Integrated error on E and B are, %f \t %f\n", L2_E, L2_B);

	/*Visual verification of the code compiling correctly*/
	printf("Code ran well for %.2f second evolution\n", tf);

	

	// Output to a file
    // ------------------------------------------------------------------------
    FILE* outfile = fopen("em-wave-v5.dat", "w");

    for (int i = 1; i < NUM_CELLS-1; ++i)
    {
        double x = xi + i * dx;
        fprintf(outfile, "%f %f %f %f %f\n", x, E[i], E_r[i], B[i], B_r[i]);
    }
    
    fclose(outfile);
	return 0;
}