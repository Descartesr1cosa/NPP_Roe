#pragma once
#include "Array_define.h"
class parameter
{
public:
	int N;
	double length, total_time, dx, CFL;
	double grid_vel;
	// double Re_m_inver;
	double bx, gamma;
	bool if_NPP;
};

