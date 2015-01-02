/*
 * MEP 602 - Numerical Methods in Energy Science ~ Dr. Hesham Othman
 * Cairo University - Faculty of Engineering
 * Mechanical Power Engineering Dept.,
 * 
 * Assignment #1
 * By: Mohamed Ahmed Emara
 * BN: 86
 * 
 * Utilities files : Contains common functions used by all methods
 * 
 * */

module utilities;

public import std.stdio, std.file, std.conv, std.math;

double maxValue(double[][] grid)
{
	auto ret = grid[0][0];
	foreach(row; grid)
	{
		foreach(elem; row)
		{
			if (elem > ret)
				ret = elem;
		}
	}
	return ret;
}

double[][] difference(double[][] A, double[][] B)
{
	auto ret = new double[][](A.length, A[0].length);
	for(auto j=0; j<A.length; j++) {
		for(auto i=0; i<A[0].length; i++) {
			ret[j][i] = A[j][i] - B[j][i];
		}
	}
	return ret;
}

void copyMatrix(double[][] src, ref double[][] dst)
{
	for(auto j=0; j<src.length; j++) {
		for(auto i=0; i<src[0].length; i++) {
			dst[j][i] = src[j][i];
		}
	}
}

double[][] divMatrix(double[][] A, double[][] B)
{
	auto ret = new double[][](A.length, A[0].length);
	for(auto j=0; j<A.length; j++) {
		for(auto i=0; i<A[0].length; i++) {
			ret[j][i] = A[j][i] / B[j][i];
		}
	}
	return ret;
}

void print_matrix_to_file(string filename, double[][] Grid, double delta=0,bool excel=false)
{

	File file = File(filename, "w");

	if (excel)
	{
		file.write("Y/X,");
		for(auto i=0; i<Grid[0].length; i++) {
			if(i == Grid[0].length - 1){
				file.writef("%g\n", i * delta);
				continue;
			}
			file.writef("%g,", i * delta);
		}
	}
	for(auto j=0; j< Grid.length; j++) {
		if (excel)
		{
			file.writef("%g,", j * delta);
		}
		for(auto i=0; i<Grid[0].length; i++) {
			file.write(Grid[j][i]);
			if ( i == Grid[0].length-1) {
				file.write("\n");
				continue;
			}
			file.write(",");
		}
	}
}

double temperature_at_location(double[][] grid, double x, double y, double delta)
{
	auto Ymax = to!int(round((1.5 - y)/delta));
	auto Xmax = to!int(round((x)/delta));
	return grid[Ymax][Xmax];
}

double[][] get_solution_matrix(double delta)
{
	double MeshHeight = 1.5;
	double MeshWidth = 1;
	
	int Ymax = to!int(round(MeshHeight/delta)) + 1;
	int Xmax = to!int(round(MeshWidth/delta))  + 1;
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
					Tn[j][i] = ((Tn[j][i+1]) + (200*delta)) / (1 + (20*delta));
					continue;
				}
				Tn[j][i] = ((omega*(Tn[j][i+1] + Tn[j][i-1] + Tn[j+1][i] + Tn[j-1][i])) / 4) + ( (1-omega)*Mesh[j][i] );
			}
		}
		
		iterations++;
		auto diff = divMatrix(difference(Tn, Mesh), Mesh);
		epsilon = maxValue(diff);
	}
	return Tn;
}

double[][] expand_matrix(double[][] matrix, double old_size, double new_size)
{
	int Ymax = to!int(round(1.5/new_size)) + 1;
	int Xmax = to!int(round(1/new_size))  + 1;
	auto ret = new double[][](Ymax, Xmax);
	
	foreach(ref row; ret) {
		foreach(ref elem; row) {
			// Pick an arbitary value that's even higher than the highest temperature
			// to make it easy to detect values that will be skipped from Extrapolation
			elem = 111.0;
		}
		row[row.length - 1] = 0.0;
	}
	
	ret[0][] = 0;
	ret[Ymax-1][] = 0;
	
	for (auto j=1; j<matrix.length-1; j++) {
		for (auto i=0; i<matrix[0].length - 1; i++) {
			auto new_i = to!int( round((i*old_size)/new_size) );
			auto new_j = to!int( round((j*old_size)/new_size) );
			ret[new_j][new_i] = matrix[j][i];
		}
	}
	return ret;
}