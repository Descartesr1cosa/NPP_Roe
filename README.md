# NPP_Roe
1D solver for ideal MHD Riemann Problems, with our new scheme: Numerical Path Preserving Roe scheme.
//===========================================================================================================================
//Compile
//use CMakeLists.txt:
  cd build
  cmake ..
  make -j
//copy the program "NPP_Roe" to ./0_run
  cp ./NPP_Roe ../0_run
//Run and enjoy
  ./NPP_Roe

//===========================================================================================================================
//setup
//An input.txt is required for NPP_Roe to get some setup information. An example is :
//-------------------Example----------------------
  Number_of_Points:				1000
  Length_of_Simulation:			3.0
  Bx:							0.7746
  CFL:							0.5
  total_time:						0.3
  gamma:						1.6667
  if_NPP(0.or.1)					0
  ===========Initial_Left==================
  rho--------u-------------v-----------w-----------------By--------------Bz-----------------p
  1		0		0.0		0.0			0.7746		0.0			0.5
  ==========Initial_Right==================
  rho--------u------------v-----------w------------------By--------------Bz-----------------p
  0.2		0		0.0		0.0			-0.7746		0.1			0.12
//--------------End-of-Example--------------------
//change "if_NPP" to control whether Classical Roe scheme or NPP Roe scheme is used for solving MHD Riemann problem
// Warning: Bx can not be too small or zero, which is not test by developer.

//===========================================================================================================================
//Postprocess
//Three cases is porvided by developer in the  end of "./0_run/input.txt"
//run the program NPP_Roe, and the solution "./output.txt" will be obtained
//copy the solution to the corresponding file (Example1, Example2, or Example3)
//A tecplot layout file is provided.
//NOTE: Postprocesses of Example1 and Example2 require Exact MHD Rimeann solutions respect to their cases, which are also provided.
//All the Exact MHD Riemann solutions are REGULAR solution (without any intermediate shock)
//More information coule be obtained by literature:
//<Exact ideal magnetohydrodynamic Riemann solutions considering the strength of intermediate shocks DOI:10.1063/5.0185483>
//<Numerical path preserving Godunov schemes for hyperbolic systems doi.org/10.1016/j.jcp.2023.112297>
//<Numerical Path Preserving Roe Scheme for Ideal MHD Riemann Problem: Complete Elimination of Pseudo-Convergence>
//Or contact me: 1905xuke@buaa.edu.cn
