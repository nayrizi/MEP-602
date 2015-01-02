/*
 * MEP 602 - Numerical Methods in Energy Science ~ Dr. Hesham Othman
 * Cairo University - Faculty of Engineering
 * Mechanical Power Engineering Dept.,
 * 
 * Assignment #1
 * By: Mohamed Ahmed Emara
 * BN: 86
 * 
 * SOR.d : Solves the heat equation using Gauss-Seidel method
 * 		   with Successive Over-Relaxation
 * 
 * */

import utilities;

void main()
{
	double MeshHeight = 1.5;
	double MeshWidth = 1;
	double delta_x = 0.025;
	double delta_y = delta_x;

	int Ymax = to!int(round(MeshHeight/delta_y)) + 1;
	int Xmax = to!int(round(MeshWidth/delta_x))  + 1;
	double[][] Mesh = new double[][](Ymax, Xmax);


	// setting boundary conditions
	foreach(ref row; Mesh) {
		foreach(ref elem; row) {
			elem = 0.0;
		}
	}
	Mesh[0][] = 100;
	Mesh[0][Xmax - 1] = 50;

	double epsilon = 1;
	double omega = 1.875;
	int iterations = 0;
	double[][] Tn = new double[][](Ymax, Xmax);
	copyMatrix(Mesh, Tn);

	while (epsilon >= 0.01)
	{
		copyMatrix(Tn, Mesh);
		for(auto j=1; j<Tn.length-1; j++)
		{
			for(auto i=0; i<Tn[0].length-1; i++)
			{
				if (i == 0)
				{
					Tn[j][i] = ((Tn[j][i+1]) + (200*delta_x)) / (1 + (20*delta_x));
					continue;
				}
				Tn[j][i] = ((omega*(Tn[j][i+1] + Tn[j][i-1] + Tn[j+1][i] + Tn[j-1][i])) / 4) + ( (1-omega)*Mesh[j][i] );
			}
		}

		iterations++;
		auto diff = divMatrix(difference(Tn, Mesh), Mesh);
		epsilon = maxValue(diff);
	}
	writefln("Optimum omega = %g", OptimumSOR(delta_x));
	writefln("Omega = %g Iterations: %d", omega,iterations);
	writeln(temperature_at_location(Tn, 0.0, 0.1, delta_x));
	writeln(temperature_at_location(Tn, 0.2, 0.3, delta_x));
	writeln(temperature_at_location(Tn, 0.4, 0.6, delta_x));
	writeln(temperature_at_location(Tn, 0.6, 0.9, delta_x));
	writeln(temperature_at_location(Tn, 0.8, 1.2, delta_x));
	writeln(epsilon);
	//print_matrix_to_file("sor_0.0125.txt", Tn);
}

double OptimumSOR(double delta) 
{

	int Ymax = to!int(round(1.5/delta)) + 1;
	int Xmax = to!int(round(1.0/delta))  + 1;
	auto zeta = pow((cos(PI/(Xmax-1)) + cos(PI/(Ymax-1)) )/2, 2);
	auto omega = 2*( (1 - sqrt(1.0-zeta)) / (zeta));

	return omega;
}