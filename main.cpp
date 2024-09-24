#include <fstream>
#include <iostream>
#include<string>
#include "Array_define.h"
#include "class.h"
#include "solve.h"
#include "parameter.h"

using namespace std;

int main()
{
	string string1;
	double temp;


	parameter* par = new parameter;
	
	//================================================================
	ifstream myfile(".//input.txt");
	if (!myfile.is_open())
	{
		cout << "Fail to open './input.txt' ! !  \tT_T" << endl;
		system("pause");
		exit(1);
	}
	myfile >> string1 >> par->N;
	myfile >> string1 >> par->length;
	myfile >> string1 >> par->bx;
	myfile >> string1 >> par->CFL;
	myfile >> string1 >> par->total_time;
	myfile >> string1 >> par->gamma;
	myfile >> string1 >> par->if_NPP;

	data_define* left = new data_define(par);
	data_define* right = new data_define(par);
	left->bx = par->bx;
	right->bx = left->bx;

	myfile >> string1;
	myfile >> string1;
	myfile >> left->rho>> left->u>> left->v>> left->w>> left->by>> left->bz>> left->p;
	myfile >> string1;
	myfile >> string1;
	myfile >> right->rho >> right->u >> right->v >> right->w >> right->by >> right->bz >> right->p;
	myfile.close();

	left->calc();
	right->calc();
	//================================================================

	//========================================================================================
	// initialization
	par->dx = par->length / (par->N + 0.0);
	Array1D<data_define> mesh(par->N+6, 3);
	//From -3~par->N+2 calc From 0~N-1

	for (int i = -3; i < par->N+3; i++)
	{
		temp = -0.5* par->length+ par->length*(i+0.0)/(par->N+0.0);
		if (temp < 0.0)
		{
			mesh(i) = *left;
		}
		else
		{
			mesh(i) = *right;
		}
		mesh(i).x = temp;
		mesh(i).calc();
	}
	//========================================================================================
	//Solver begin
	solver solu;
	solu.Solve_half(mesh, par);

	//=========================================================================================
	// output the solution
	output(mesh, par);
	//=========================================================================================

	delete right;
	delete left;
	delete par;
	return 0;
}