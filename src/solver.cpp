#include "solve.h"
#include <fstream>

void solver::Solve_half(Array1D<data_define> &mesh, parameter *par)
{
	double tt = 0.0, dt = 0.0;
	double2D fluxp(par->N + 1, 7, 1);	 //-1~par->N-1==-1~mx
	double2D fluxm(par->N + 1, 7, 1);	 //-1~par->N-1==-1~mx

	int Nstep = 0;
	while (true)
	{
		Nstep++;
		if (fmod(Nstep, 1000) == 0)
		{
			output(mesh, par);
		}
		dt = find_dt(mesh, par);

		rhs_Roe(mesh, fluxm, fluxp, par);

		for (int i = 0; i < par->N; i++)
		{

			for (int j = 0; j < 7; j++)
			{
				mesh(i).q[j] = mesh(i).q[j] - dt * (fluxm(i, j) + fluxp(i - 1, j)) / par->dx;
			}
			mesh(i).recalc();
			mesh(i).calc();
		}
		std::cout << Nstep << "\t" << tt << std::endl;
		tt = tt + dt;
		if (tt > par->total_time)
			break;
	}
}

double solver::find_dt(Array1D<data_define> &mesh, parameter *par)
{
	double temp = 0.0, rr;
	for (int i = 0; i < par->N; i++)
	{
		rr = mesh(i).get_radius();
		if (std::isnan(rr))
		{
			std::cout << "NaN!!!!!!!!!!!!!!" << std::endl;
			exit(-1);
		}
		if (temp < rr)
		{
			temp = rr;
		}
	}
	temp = par->CFL * par->dx / temp;
	par->grid_vel = par->dx / temp;
	return temp;
}

void output(Array1D<data_define> &mesh, parameter *par)
{
	//=========================================================================================
	//output
	std::ofstream outfile(".//output.dat");
	outfile << "variables=\"x\",\"rho\",\"u\",\"v\",\"w\",\"Bx\",\"By\",\"Bz\",\"p\"" << std::endl;
	for (int i = 0; i < par->N; i++)
	{
		outfile << mesh(i).x << "\t" << mesh(i).rho << "\t" << mesh(i).u << "\t" << mesh(i).v << "\t" << mesh(i).w << "\t"
				<< mesh(i).bx << "\t" << mesh(i).q[4] << "\t" << mesh(i).q[5] << "\t" << mesh(i).p
				<< std::endl;
	}
	outfile.close();
	//=========================================================================================
}

double delta_0(double x)
{
	double temp;
	if (abs(x) < 0.6)
	{
		temp = 1.0;
	}
	else
	{
		temp = 0.0;
	}
	return temp;
}
double sign(double x)
{
	double temp;
	if (x >= 0.0)
	{
		temp = 1.0;
	}
	else
	{
		temp = -1.0;
	}
	return temp;
}
