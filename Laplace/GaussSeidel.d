/*
 * MEP 602 - Numerical Methods in Energy Science ~ Dr. Hesham Othman
 * Cairo University - Faculty of Engineering
 * Mechanical Power Engineering Dept.,
 * 
 * Assignment #1
 * By: Mohamed Ahmed Emara
 * BN: 86
 * 
 * GaussSeidel.d : Solves the heat equation using Gauss-Seidel method
 * 				   using a stopping criteria of epsilon = 1%
 * 
 * */

import utilities;


void main()
{
	double MeshHeight = 1.5;
	double MeshWidth = 1;
	double delta_x = 0.0125; //change this value for different grid sizes
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
	double omega = 1;
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
	writefln("Delta = %g", delta_x);
	writefln("Iterations: %d", iterations);
	writefln("Epsilon   : %g", epsilon);
	writefln("Ny = %s ... Nx = %s", Ymax, Xmax);
	print_matrix_to_file("gauss_seidel_0.0125.txt", Tn);
}