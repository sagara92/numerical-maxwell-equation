#include<stdio.h>
#include <math.h>
#define NUM_CELLS 2500
/*Yee Lattice and Leapfrogging*/


int main()
{	
	/*Defining main variables*/

	double x, y, t;
	double E[NUM_CELLS]; /*Cells for E field*/
	double B[NUM_CELLS]; /*Cells for E field*/
	double B_r[NUM_CELLS]; /*Cells for true B field*/
	double E_r[NUM_CELLS]; /*Cells for true E field*/
	double E0 = 1.0; /*Amplitude of E*/
	double B0 = 1.0; /*Amplitude of B*/
	double xi = 0.0; /*Initial position for the E and B field*/
	double xf = 1.0; /*Final position for the E and B field*/
	double c  = 1.0; /*Speed of light*/
	double ti = 0.0; /*Starting time*/
	double tf = 0.5; /*Evolution time in seconds*/
	double k  = 2.0 * M_PI /(xf - xi); /*Wave number*/
	double w  = k * c; 
	double dx = (xf - xi)/(NUM_CELLS); /*Spacing of the cells*/
	double dt = 0.04 * (dx / c);


	//------------------------------------------//
	/*Initial wave at t = ti = 0*/

	for (int i = 0; i < NUM_CELLS; ++i)
	{
		x = xi + i*dx;
		E[i] = E0 * sin (k*x + w*ti);
		B[i] = B0 * sin (k*(x+0.5*dx) + w*(ti+0.5*dt));
	}
	//-----------------------------------------//


	//***********************************************************************//
	/*Time evolution*/

	while(t <= tf)
	{
		/*Evolving the true solution*/

		//------------------------------------//
		for (int i = 0; i < NUM_CELLS; ++i)
		{
			x = xi + i*dx;
			E_r[i] = E0 * sin (k*x + w*t);
			B_r[i] = B0 * sin (k*(x+0.5*dx) + w*(t+0.5*dt));
		}
		//------------------------------------//

		/*Defining variables*/

		double dEdx[NUM_CELLS], dBdx[NUM_CELLS];

		/*Defining boundary conditions and updating E*/
		for (int i = 0; i < NUM_CELLS; ++i)
		{
			double Bplus  = B[i];
			double Bminus = B[(i-1+NUM_CELLS)%NUM_CELLS];
		
			dBdx[i] = (Bplus - Bminus)/dx;
			E[i] += c * dBdx[i] * dt;
			
		}

		/*Defining boundary conditions and updating E*/
		for (int i = 0; i < NUM_CELLS; ++i)
		{
			double Eplus  = E[(i+1)%NUM_CELLS];
			double Eminus = E[i];

			dEdx[i] = (Eplus - Eminus)/dx;
			B[i] += c * dEdx[i] * dt;
		}
		
		t+=dt;
	}
	//***********************************************************************//
	printf("Code ran well for %f seconds.\n", tf); /*Visual verification*/

	/*Error analysis*/

	double E2=0.0;
	double B2=0.0;
	double Diff_E[NUM_CELLS];
	double Diff_B[NUM_CELLS];
	for (int i = 0; i < NUM_CELLS; ++i)
	{
		/*Calculating L2 errors in each position*/			
		Diff_E[i] = E_r[i] - E[i];		/*Difference in the calculated value and real value*/	
		E2 += (Diff_E[i]*Diff_E[i])*dx;
					
		Diff_B[i] = B_r[i] - B[i];
		B2 += (Diff_B[i]*Diff_B[i])*dx;
		
	}
	double L2_E = sqrt(E2)/(xf-xi);
	double L2_B = sqrt(B2)/(xf-xi);
	printf("Integrated error on E and B are, %f \t %f\n", L2_E, L2_B);

	//-----------------------------------------------------------------//
	/*Define function*/

	//----------------------------------------------------------------//
	/*Output to a file*/

	FILE* outfile = fopen("maxwell-v6.dat", "w");

	for (int i = 0; i < NUM_CELLS; ++i)
	{
		x = xi + i*dx;
		fprintf(outfile, "%f %f %f %f %f\n", x, E_r[i], E[i], B_r[i], B[i]);
	}
	fclose(outfile);
	return 0;
}