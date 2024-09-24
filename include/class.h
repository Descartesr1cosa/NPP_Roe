#include <math.h>
#include "parameter.h"
#pragma once
class data_define
{
public:
	double x, u, v, w, p, rho;
	double bx;
	double by, bz;
	double gamma;

	double q[7];
	double f[7];
	double rhoe, b_sq;

	data_define() = default;
	data_define(parameter *par){
		gamma = par->gamma;
	};

	data_define(const data_define& data_in)
	{
		this->x = data_in.x;
		this->u = data_in.u;
		this->v = data_in.v;
		this->w = data_in.w;
		this->p = data_in.p;
		this->rho = data_in.rho;
		this->bx = data_in.bx;
		this->by = data_in.by;
		this->bz = data_in.bz;
		this->gamma = data_in.gamma;
		this->calc();
	}

public:

	void calc() {
		rhoe = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
		b_sq = bx * bx + by * by + bz * bz;

		q[0] = rho;
		q[1] = rho * u;
		q[2] = rho * v;
		q[3] = rho * w;
		q[4] = by;
		q[5] = bz;
		q[6] = rhoe + 0.5 * b_sq;

		f[0] = q[1];
		f[1] = q[1] * u + p + b_sq * 0.5;
		f[2] = q[1] * v - bx * by;
		f[3] = q[1] * w - bx * bz;
		f[4] = by * u - bx * v;
		f[5] = bz * u - bx * w;
		f[6] = (q[6] + p + b_sq * 0.5) * u - bx * (bx * u + by * v + bz * w);
		
	}
	void recalc() {
		double b_sq0;

		rho = q[0];
		u = q[1] / q[0];
		v = q[2] / q[0];
		w = q[3] / q[0];
		by = q[4];
		bz = q[5];

		b_sq0= bx * bx + by * by + bz * bz;

		p = q[6] - 0.5 * b_sq0 - 0.5 * rho * (u * u + v * v + w * w);
		p = p * (gamma - 1.0);

	}



public:

	double get_radius() {
		double cfp, cfm, cap, cam, csm, csp, cu;
		double rr;
		cfp = this->get_fp();
		cfm = this->get_fm();
		cap = this->get_ap();
		cam = this->get_am();
		csm = this->get_sm();
		csp = this->get_sp();
		cu = this->get_u();
		rr = fmax(fabs(cfp), fabs(cfm));
		rr = fmax(fabs(cap), rr);
		rr = fmax(fabs(cam), rr);
		rr = fmax(fabs(csp), rr);
		rr = fmax(fabs(csm), rr);
		rr = fmax(fabs(cu), rr);

		return rr;
	}

	double get_fp() {
		double temp;
		double c_sq;
		c_sq = gamma * p / rho;
		temp = b_sq / rho + c_sq;
		temp = temp * temp;
		temp = temp - 4 * c_sq * bx * bx / rho;
		temp = 0.5 * (sqrt(temp) + c_sq + b_sq / rho);
		temp = sqrt(temp);

		temp = u + temp;

		return temp;
	}
	double get_sp() {
		double temp;
		double c_sq;
		c_sq = gamma * p / rho;
		temp = b_sq / rho + c_sq;
		temp = temp * temp;
		temp = temp - 4 * c_sq * bx * bx / rho;
		temp = 0.5 * (-sqrt(temp) + c_sq + b_sq / rho);
		temp = sqrt(temp);

		temp = u + temp;


		return temp;
	}
	double get_ap() {
		double temp;

		temp = sqrt(bx * bx / rho);

		temp = u + temp;

		return temp;
	}

	double get_fm() {
		double temp;
		double c_sq;
		c_sq = gamma * p / rho;
		temp = b_sq / rho + c_sq;
		temp = temp * temp;
		temp = temp - 4 * c_sq * bx * bx / rho;
		temp = 0.5 * (sqrt(temp) + c_sq + b_sq / rho);
		temp = sqrt(temp);

		temp = u - temp;

		return temp;
	}
	double get_sm() {
		double temp;
		double c_sq;
		c_sq = gamma * p / rho;
		temp = b_sq / rho + c_sq;
		temp = temp * temp;
		temp = temp - 4 * c_sq * bx * bx / rho;
		temp = 0.5 * (-sqrt(temp) + c_sq + b_sq / rho);
		temp = sqrt(temp);

		temp = u - temp;


		return temp;
	}
	double get_am() {
		double temp;

		temp = sqrt(bx * bx / rho);

		temp = u - temp;

		return temp;
	}
	double get_u()
	{
		return u;
	}
};
