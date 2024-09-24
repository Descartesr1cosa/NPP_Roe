#include "solve.h"
void solver::NPP_Modification_Roe(data_define &left, data_define &right, eigen &my, parameter *par, double &rho, double &u, double &v, double &w, double &betay, double &betaz, double &sig)
{
	//======================================================================================================================
	double dv, dw, dBy, dBz;
	dv = right.v - left.v;
	dw = right.w - left.w;
	dBy = right.by - left.by;
	dBz = right.bz - left.bz;
	data_define uL(par), uR(par), uu(par);
	//======================================================================================================================


	//======================================================================================================================
	// NPP Modification of Left-going Alfvenic wave

	uL = left;

	// State after the wave of [left-fast: u-cf]
	for (int k = 0; k < 7; k++)
		uL.q[k] = uL.q[k] + my.alpha[0] * my.rr(k, 0);
	uL.recalc();

	// State after the wave of [left-Alfvenic: u-ca]
	for (int k = 0; k < 7; k++)
		uR.q[k] = uL.q[k] + my.alpha[1] * my.rr(k, 1);
	uR.recalc();

	double theta_L = atan2(uL.bz, uL.by);
	double theta_R = atan2(uR.bz, uR.by);
	// theta_L theta_R are between -pi to pi

	if (theta_R - theta_L > 3.1415926535)
	{
		theta_L = theta_L + 2.0 * 3.1415926535;
	}
	else if (theta_R - theta_L < -3.1415926535)
	{
		theta_R = theta_R + 2.0 * 3.1415926535;
	}

	double dtheta = fabs(theta_R - theta_L);

	double lambda = my.lambda[1] / par->grid_vel;
	double signal = (lambda >= 0.0) ? 0.5 : -0.5;
	double theta = (0.5 + signal - lambda) * theta_R + (0.5 - signal + lambda) * theta_L;
	double theta_m = 0.5 * ((0.5 + signal) * theta_R + (0.5 - signal) * theta_L + theta);
	double dtheta2 = (signal > 0.0) ? fabs(theta_R - theta) : fabs(theta - theta_L);

	double betaz1 = sin(theta_m);
	double betay1 = cos(theta_m);
	if (betaz1 * betaz + betay1 * betay < 0.0)
	{
		betaz1 = -betaz1;
		betay1 = -betay1;
	}

	//----------------------------------------------------
	// Modify the NPP eigenvector
	my.rr(2, 1) = -rho * betaz1;
	my.rr(3, 1) = rho * betay1;
	my.rr(4, 1) = -sig * sqrt(rho) * betaz1;
	my.rr(5, 1) = sig * sqrt(rho) * betay1;
	my.rr(6, 1) = -rho * (v * betaz1 - w * betay1);

	// Modify the strength of left-Alfvenic wave
	if (fabs(sin(0.5 * dtheta)) > 1e-8)
		my.alpha[1] = my.alpha[1] * fabs(sin(0.5 * dtheta2) / sin(0.5 * dtheta));

	// Modify the NPP eigenvalue
	my.lambda[1] = par->grid_vel * 2.0 * signal;
	//----------------------------------------------------

	double v_, w_, By_, Bz_;
	if (signal > 0.0)
	{
		//State after Modify u-ca wave
		for (int k = 0; k < 7; k++)
			uu.q[k] = uR.q[k] - my.alpha[1] * my.rr(k, 1);
		uu.recalc();

		v_ = 0.5 * (uu.v + uR.v);
		w_ = 0.5 * (uu.w + uR.w);
		By_ = 0.5 * (uu.by + uR.by);
		Bz_ = 0.5 * (uu.bz + uR.bz);
		dw = uR.w - uu.w;
		dv = uR.v - uu.v;
		dBz = uR.bz - uu.bz;
		dBy = uR.by - uu.by;
	}
	else
	{
		//State after Modify u-ca wave
		for (int k = 0; k < 7; k++)
			uu.q[k] = uL.q[k] + my.alpha[1] * my.rr(k, 1);
		uu.recalc();

		v_ = 0.5 * (uu.v + uL.v);
		w_ = 0.5 * (uu.w + uL.w);
		By_ = 0.5 * (uu.by + uL.by);
		Bz_ = 0.5 * (uu.bz + uL.bz);
		dw = uu.w - uL.w;
		dv = uu.v - uL.v;
		dBz = uu.bz - uL.bz;
		dBy = uu.by - uL.by;
	}
	// Conservative modification
	my.rr(6, 1) = -rho * (v_ * betaz1 - w_ * betay1) + By_ * (-sig * sqrt(rho) * betaz1) + Bz_ * (sig * sqrt(rho) * betay1);
	//======================================================================================================================

	//======================================================================================================================
	// NPP Modification of Right-going Alfvenic wave

	uR = right;
	// State after the wave of [right-fast: u+cf]
	for (int k = 0; k < 7; k++)
		uR.q[k] = uR.q[k] - my.alpha[6] * my.rr(k, 6);
	uR.recalc();

	// State after the wave of [right-Alfvenic: u+ca]
	for (int k = 0; k < 7; k++)
		uL.q[k] = uR.q[k] - my.alpha[5] * my.rr(k, 5);
	uL.recalc();

	theta_L = atan2(uL.bz, uL.by);
	theta_R = atan2(uR.bz, uR.by);
	// theta_L theta_R are between -pi to pi

	if (theta_R - theta_L > 3.1415926535)
	{
		theta_L = theta_L + 2.0 * 3.1415926535;
	}
	else if (theta_R - theta_L < -3.1415926535)
	{
		theta_R = theta_R + 2.0 * 3.1415926535;
	}

	dtheta = fabs(theta_R - theta_L);

	lambda = my.lambda[5] / par->grid_vel;
	signal = (lambda >= 0.0) ? 0.5 : -0.5;
	theta = (0.5 + signal - lambda) * theta_R + (0.5 - signal + lambda) * theta_L;
	theta_m = 0.5 * ((0.5 + signal) * theta_R + (0.5 - signal) * theta_L + theta);
	dtheta2 = (signal > 0.0) ? fabs(theta_R - theta) : fabs(theta - theta_L);

	betaz1 = sin(theta_m);
	betay1 = cos(theta_m);
	if (betaz1 * betaz + betay1 * betay < 0.0)
	{
		betaz1 = -betaz1;
		betay1 = -betay1;
	}

	//----------------------------------------------------
	// Modify the NPP eigenvector
	my.rr(2, 5) = rho * betaz1;
	my.rr(3, 5) = -rho * betay1;
	my.rr(4, 5) = -sig * sqrt(rho) * betaz1;
	my.rr(5, 5) = sig * sqrt(rho) * betay1;
	my.rr(6, 5) = rho * (v * betaz1 - w * betay1);

	// Modify the strength of left-Alfvenic wave
	if (fabs(sin(0.5 * dtheta)) > 1e-8)
		my.alpha[5] = my.alpha[5] * fabs(sin(0.5 * dtheta2) / sin(0.5 * dtheta));
	
	// Modify the NPP eigenvalue
	my.lambda[5] = par->grid_vel * 2.0 * signal;
	//----------------------------------------------------

	if (signal > 0.0)
	{
		//State after Modify u+ca wave
		for (int k = 0; k < 7; k++)
			uu.q[k] = uR.q[k] - my.alpha[5] * my.rr(k, 5);
		uu.recalc();

		v_ = 0.5 * (uu.v + uR.v);
		w_ = 0.5 * (uu.w + uR.w);
		By_ = 0.5 * (uu.by + uR.by);
		Bz_ = 0.5 * (uu.bz + uR.bz);
		dw = uR.w - uu.w;
		dv = uR.v - uu.v;
		dBz = uR.bz - uu.bz;
		dBy = uR.by - uu.by;
	}
	else
	{
		//State after Modify u+ca wave
		for (int k = 0; k < 7; k++)
			uu.q[k] = uL.q[k] + my.alpha[5] * my.rr(k, 5);
		uu.recalc();

		v_ = 0.5 * (uu.v + uL.v);
		w_ = 0.5 * (uu.w + uL.w);
		By_ = 0.5 * (uu.by + uL.by);
		Bz_ = 0.5 * (uu.bz + uL.bz);

		dw = uu.w - uL.w;
		dv = uu.v - uL.v;
		dBz = uu.bz - uL.bz;
		dBy = uu.by - uL.by;
	}
	// Conservative modification
	my.rr(6, 5) = rho * (v_ * betaz1 - w_ * betay1) + By_ * (-sig * sqrt(rho) * betaz1) + Bz_ * (sig * sqrt(rho) * betay1);
}
