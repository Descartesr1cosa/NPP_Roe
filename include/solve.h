#pragma once

#include "class.h"
#include <iostream>
#include "parameter.h"
#include "Array_define.h"

class eigen
{
public:
	eigen()
	{
		for (int i = 0; i < 7; i++)
		{
			lambda[i] = 0.0;
			alpha[i] = 0.0;
		}
		rr.SetSize(7, 7, 0);
		ll.SetSize(7, 7, 0);
	}

public:
	double lambda[7], alpha[7];
	double2D rr, ll;
};

class solver
{

public:
	void Solve_half(Array1D<data_define> &mesh, parameter *par);


	double find_dt(Array1D<data_define> &mesh, parameter *par);

	void rhs_Roe(Array1D<data_define> &mesh, double2D &flux_m, double2D &flux_p, parameter *par);

	void calc_eigen_system_half(data_define &left, data_define &right, eigen &my, parameter *par);
	// void calc_eigen(data_define &left, data_define &right, eigen &myeigen);

	void NPP_Modification_Roe(data_define &left, data_define &right, eigen &my, parameter *par, double &rho, double &u, double &v, double &w, double &betay, double &betaz, double &sig);
};

void output(Array1D<data_define> &mesh, parameter *par);

double delta_0(double x);

double sign(double x);
