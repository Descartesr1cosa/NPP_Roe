#include "solve.h"
void solver::rhs_Roe(Array1D<data_define> &mesh, double2D &flux_m, double2D &flux_p, parameter *par)
{
	eigen myeigen;
	for (int i = -1; i < par->N; i++) //flux_num[i](f_{i+1/2}) is calculated by mesh[i] and mesh[i+1]
	{

		//Calculate Roe average matrix and the eigenvalues, eigenvectors, strength of wave
		for (int j = 0; j < 7; j++)
		{
			flux_m(i, j) = 0.0; // mesh(i).f[j];
			flux_p(i, j) = 0.0; // mesh(i + 1).f[j];
		}

		double delta = 0.0;
		for (int j = 0; j < 7; j++)
			delta += fabs(mesh(i).q[j] - mesh(i + 1).q[j]);
		if (delta < 1e-8)
			continue;

		calc_eigen_system_half(mesh(i), mesh(i + 1), myeigen, par);

		// For seven characteristic wave
		for (int j = 0; j < 7; j++)
		{
			if (myeigen.lambda[j] <= 0.0)
			{
				for (int k = 0; k < 7; k++) //The k-th component of j-th characteristic wave
					flux_m(i, k) = flux_m(i, k) + myeigen.lambda[j] * myeigen.alpha[j] * myeigen.rr(k, j);
			}
			else if (myeigen.lambda[j] >= 0.0)
			{
				for (int k = 0; k < 7; k++) //The k-th component of j-th characteristic wave
					flux_p(i, k) = flux_p(i, k) + myeigen.lambda[j] * myeigen.alpha[j] * myeigen.rr(k, j);
			}
		}
	}
}
void solver::calc_eigen_system_half(data_define &left, data_define &right, eigen &my, parameter *par)
{
	double gamma = par->gamma;
	double rho, u, v, w, X, Bx_s, By_s, Bz_s, H_t;
	double b1, b2, b3;
	double betay, betaz;
	double r1, r2;

	double cf, cs, ca, cc2, alp_s, alp_f, sig;
	double drho, dp, du, dv, dw, dby, dbz;

	double temp;

	r1 = sqrt(left.rho);
	r2 = sqrt(right.rho);

	rho = sqrt(left.rho * right.rho);
	u = (r1 * left.u + r2 * right.u) / (r1 + r2);
	v = (r1 * left.v + r2 * right.v) / (r1 + r2);
	w = (r1 * left.w + r2 * right.w) / (r1 + r2);
	Bx_s = (r2 * left.bx + r1 * right.bx) / (r1 + r2);
	By_s = (r2 * left.by + r1 * right.by) / (r1 + r2);
	Bz_s = (r2 * left.bz + r1 * right.bz) / (r1 + r2);

	double theta1, theta2;
	theta1 = atan2(left.bz, left.by);
	theta2 = atan2(right.bz, right.by);

	//Best way to modify the Roe average of tangential magnetic field
	By_s = (r2 * left.by + r1 * right.by) / (r1 + r2) +  delta_0(3.1415926535 - fabs(theta1 - theta2)) * (sqrt(left.bz * left.bz + left.by * left.by) + sqrt(right.bz * right.bz + right.by * right.by)) * 0.5 * cos((theta1 + theta2) * 0.5);
	Bz_s = (r2 * left.bz + r1 * right.bz) / (r1 + r2) +  delta_0(3.1415926535 - fabs(theta1 - theta2)) * (sqrt(left.bz * left.bz + left.by * left.by) + sqrt(right.bz * right.bz + right.by * right.by)) * 0.5 * sin((theta1 + theta2) * 0.5);
	
	// By_s = (r2 * left.by + r1 * right.by) / (r1 + r2) +  sin(0.5*fabs(theta1 - theta2)) * (sqrt(left.bz * left.bz + left.by * left.by) + sqrt(right.bz * right.bz + right.by * right.by)) * 0.5 * cos((theta1 + theta2) * 0.5);
	// Bz_s = (r2 * left.bz + r1 * right.bz) / (r1 + r2) +  sin(0.5*fabs(theta1 - theta2)) * (sqrt(left.bz * left.bz + left.by * left.by) + sqrt(right.bz * right.bz + right.by * right.by)) * 0.5 * sin((theta1 + theta2) * 0.5);

	//By_s = (r2 * left.by + r1 * right.by) / (r1 + r2) + delta_0(fabs(left.by + right.by)) * delta_0(fabs(left.bz + right.bz)) * delta_0(3.1415926535 - fabs(theta1 - theta2)) * (sqrt(left.bz * left.bz + left.by * left.by) + sqrt(right.bz * right.bz + right.by * right.by)) * 0.5 * cos((theta1 + theta2) * 0.5);
	//Bz_s = (r2 * left.bz + r1 * right.bz) / (r1 + r2) + delta_0(fabs(left.by + right.by)) * delta_0(fabs(left.bz + right.bz)) * delta_0(3.1415926535 - fabs(theta1 - theta2)) * (sqrt(left.bz * left.bz + left.by * left.by) + sqrt(right.bz * right.bz + right.by * right.by)) * 0.5 * sin((theta1 + theta2) * 0.5);

	 //double bt_L,bt_R;
	 //bt_L=sqrt(left.bz * left.bz + left.by * left.by);
	 //bt_R=sqrt(right.bz * right.bz + right.by * right.by);
	 //By_s = (r2 * bt_L + r1 * bt_R) / (r1 + r2)*cos((theta1 + theta2) * 0.5);
	 //Bz_s = (r2 * bt_L + r1 * bt_R) / (r1 + r2)*sin((theta1 + theta2) * 0.5);

	X = 0.5 * ((left.by - right.by) * (left.by - right.by) + (left.bz - right.bz) * (left.bz - right.bz)) / (left.rho + right.rho + 2.0 * r1 * r2);
	temp = left.q[6] / left.rho + left.p / left.rho +
		   0.5 * (left.bx * left.bx + left.by * left.by + left.bz * left.bz) / left.rho;
	H_t = temp;
	temp = right.q[6] / right.rho + right.p / right.rho +
		   0.5 * (right.bx * right.bx + right.by * right.by + right.bz * right.bz) / right.rho;
	H_t = (H_t * r1 + temp * r2) / (r1 + r2);
	b1 = Bx_s / sqrt(rho);
	b2 = By_s / sqrt(rho);
	b3 = Bz_s / sqrt(rho);

	ca = b1;
	cc2 = (2.0 - gamma) * X + (gamma - 1.0) * (H_t - 0.5 * (u * u + v * v + w * w) - b1 * b1 - b2 * b2 - b3 * b3); // cc^2
	temp = cc2 + b1 * b1 + b2 * b2 + b3 * b3;																	   // a*^2
	cf = 0.5 * (temp + sqrt(temp * temp - 4 * cc2 * b1 * b1));
	cf = sqrt(cf);
	cs = 0.5 * (temp - sqrt(temp * temp - 4 * cc2 * b1 * b1));
	cs = sqrt(cs);

	//============output eigenvalue=====================
	my.lambda[0] = u - cf;
	my.lambda[1] = u - ca;
	my.lambda[2] = u - cs;
	my.lambda[3] = u;
	my.lambda[4] = u + cs;
	my.lambda[5] = u + ca;
	my.lambda[6] = u + cf;

	//============output wave stregth===================
	alp_f = sqrt((cc2 - cs * cs) / (cf * cf - cs * cs));
	alp_s = sqrt((cf * cf - cc2) / (cf * cf - cs * cs));

	sig = sign(left.bx);

	betay = By_s / sqrt(By_s * By_s + Bz_s * Bz_s);
	betaz = Bz_s / sqrt(By_s * By_s + Bz_s * Bz_s);

	drho = right.rho - left.rho;
	du = right.u - left.u;
	dv = right.v - left.v;
	dw = right.w - left.w;
	dp = right.p - left.p;
	dby = right.by - left.by;
	dbz = right.bz - left.bz;

	temp = alp_f * (X * drho + dp) + rho * alp_s * cs * sig * (betay * dv + betaz * dw);
	temp = temp - rho * alp_f * cf * du + sqrt(rho * cc2) * alp_s * (betay * dby + betaz * dbz);
	my.alpha[0] = 0.5 * temp; // alpha:u-cf
	temp = alp_f * (X * drho + dp) - rho * alp_s * cs * sig * (betay * dv + betaz * dw);
	temp = temp + rho * alp_f * cf * du + sqrt(rho * cc2) * alp_s * (betay * dby + betaz * dbz);
	my.alpha[6] = 0.5 * temp; // alpha:u+cf

	temp = alp_s * (X * drho + dp) - rho * alp_f * cf * sig * (betay * dv + betaz * dw);
	temp = temp - rho * alp_s * cs * du - sqrt(rho * cc2) * alp_f * (betay * dby + betaz * dbz);
	my.alpha[2] = 0.5 * temp; // alpha:u-cs
	temp = alp_s * (X * drho + dp) + rho * alp_f * cf * sig * (betay * dv + betaz * dw);
	temp = temp + rho * alp_s * cs * du - sqrt(rho * cc2) * alp_f * (betay * dby + betaz * dbz);
	my.alpha[4] = 0.5 * temp; // alpha:u+cs

	my.alpha[1] = 0.5 * (betay * dw - betaz * dv + sig / sqrt(rho) * (betay * dbz - betaz * dby));	// alpha:u-ca
	my.alpha[5] = 0.5 * (-betay * dw + betaz * dv + sig / sqrt(rho) * (betay * dbz - betaz * dby)); // alpha:u+ca

	my.alpha[3] = (cc2 - X) * drho - dp; // alpha:u

	//=================right eigenvector==============
	//------------------eigenvector-3-----------------//u
	my.rr(0, 3) = 1.0 / cc2;
	my.rr(1, 3) = u / cc2;
	my.rr(2, 3) = v / cc2;
	my.rr(3, 3) = w / cc2;
	my.rr(4, 3) = 0.0;
	my.rr(5, 3) = 0.0;
	my.rr(6, 3) = 0.5 * (u * u + v * v + w * w) / cc2 + (gamma - 2.0) / (gamma - 1.0) * X / cc2;
	//------------------eigenvector-0-----------------//u-cf
	my.rr(0, 0) = alp_f / cc2;
	my.rr(1, 0) = alp_f / cc2 * (u - cf);
	my.rr(2, 0) = (alp_f * v + alp_s * cs * betay * sig) / cc2;
	my.rr(3, 0) = (alp_f * w + alp_s * cs * betaz * sig) / cc2;
	my.rr(4, 0) = alp_s * betay / sqrt(cc2 * rho);
	my.rr(5, 0) = alp_s * betaz / sqrt(cc2 * rho);
	my.rr(6, 0) = rho * alp_f * (H_t - b1 * b1 - b2 * b2 - b3 * b3 - u * cf) + rho * alp_s * cs * sig * (v * betay + w * betaz) + sqrt(rho * cc2) * alp_s * sqrt(By_s * By_s + Bz_s * Bz_s);
	my.rr(6, 0) = my.rr(6, 0) / rho / cc2;
	//------------------eigenvector-6-----------------//u+cf
	my.rr(0, 6) = alp_f / cc2;
	my.rr(1, 6) = alp_f / cc2 * (u + cf);
	my.rr(2, 6) = (alp_f * v - alp_s * cs * betay * sig) / cc2;
	my.rr(3, 6) = (alp_f * w - alp_s * cs * betaz * sig) / cc2;
	my.rr(4, 6) = alp_s * betay / sqrt(cc2 * rho);
	my.rr(5, 6) = alp_s * betaz / sqrt(cc2 * rho);
	my.rr(6, 6) = rho * alp_f * (H_t - b1 * b1 - b2 * b2 - b3 * b3 + u * cf) - rho * alp_s * cs * sig * (v * betay + w * betaz) + sqrt(rho * cc2) * alp_s * sqrt(By_s * By_s + Bz_s * Bz_s);
	my.rr(6, 6) = my.rr(6, 6) / rho / cc2;
	//------------------eigenvector-2-----------------//u-cs
	my.rr(0, 2) = alp_s / cc2;
	my.rr(1, 2) = alp_s / cc2 * (u - cs);
	my.rr(2, 2) = (alp_s * v - alp_f * cf * betay * sig) / cc2;
	my.rr(3, 2) = (alp_s * w - alp_f * cf * betaz * sig) / cc2;
	my.rr(4, 2) = -alp_f * betay / sqrt(cc2 * rho);
	my.rr(5, 2) = -alp_f * betaz / sqrt(cc2 * rho);
	my.rr(6, 2) = rho * alp_s * (H_t - b1 * b1 - b2 * b2 - b3 * b3 - u * cs) - rho * alp_f * cf * sig * (v * betay + w * betaz) - sqrt(rho * cc2) * alp_f * sqrt(By_s * By_s + Bz_s * Bz_s);
	my.rr(6, 2) = my.rr(6, 2) / rho / cc2;
	//------------------eigenvector-4-----------------//u+cs
	my.rr(0, 4) = alp_s / cc2;
	my.rr(1, 4) = alp_s / cc2 * (u + cs);
	my.rr(2, 4) = (alp_s * v + alp_f * cf * betay * sig) / cc2;
	my.rr(3, 4) = (alp_s * w + alp_f * cf * betaz * sig) / cc2;
	my.rr(4, 4) = -alp_f * betay / sqrt(cc2 * rho);
	my.rr(5, 4) = -alp_f * betaz / sqrt(cc2 * rho);
	my.rr(6, 4) = rho * alp_s * (H_t - b1 * b1 - b2 * b2 - b3 * b3 + u * cs) + rho * alp_f * cf * sig * (v * betay + w * betaz) - sqrt(rho * cc2) * alp_f * sqrt(By_s * By_s + Bz_s * Bz_s);
	my.rr(6, 4) = my.rr(6, 4) / rho / cc2;
	//------------------eigenvector-1-----------------//u-ca
	my.rr(0, 1) = 0.0;
	my.rr(1, 1) = 0.0;
	my.rr(2, 1) = -rho * betaz;
	my.rr(3, 1) = rho * betay;
	my.rr(4, 1) = -sig * sqrt(rho) * betaz;
	my.rr(5, 1) = sig * sqrt(rho) * betay;
	my.rr(6, 1) = -rho * (v * betaz - w * betay);
	//------------------eigenvector-5-----------------//u+ca
	my.rr(0, 5) = 0.0;
	my.rr(1, 5) = 0.0;
	my.rr(2, 5) = rho * betaz;
	my.rr(3, 5) = -rho * betay;
	my.rr(4, 5) = -sig * sqrt(rho) * betaz;
	my.rr(5, 5) = sig * sqrt(rho) * betay;
	my.rr(6, 5) = rho * (v * betaz - w * betay);
	//==============================================================
	if(par->if_NPP)
		NPP_Modification_Roe(left, right, my, par, rho, u, v, w, betay, betaz, sig); //NPP-Roe Modification
	//=======================================================================================
	//=======================================================================================
}
