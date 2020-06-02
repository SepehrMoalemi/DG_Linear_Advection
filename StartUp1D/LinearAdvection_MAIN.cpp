#include <iostream>
#include <complex>
#include <math.h> //pow() , fabs(), ceil()
#include <cmath> //sqrt(), floor(), sin()

/*Downloaded Libraries from Eigen*/
#include <Eigen/Dense>		  //Using Vectors and Matrices
#include <Eigen/Eigenvalues>  //Using eigen value solver

#include "Globals.h"
#include "Error_Functions.h"
#include "Startup_Functions.h"

using namespace std;
using namespace Eigen;

#pragma region Math Constants
	const double pi = 3.14159265358979323846;
	const long double tol = 2.0e-20;
	Matrix<double, 5, 1> rk4a, rk4b, rk4c;
#pragma endregion

#pragma region Global Constants
	const double ALPHA = 0.0;
	const double xmin = 0.0;
	const double xmax = 2.0 * pi;
	const double FinalTime = 10.0;
	const double NODETOL = 1e-10;
#pragma endregion 

#pragma region Global Variables
	int N, K, 
		Nfp, Np, Nv, 
		Nfaces, TotalFaces,
		vmapI, vmapO, mapI, mapO;

	Vector2i mapB, vmapB;

	VectorXd VX, w, r, P, L2Error_N;

	MatrixXi EToE, EToF, nx;

	MatrixXd u, u_exact, x, J, V, 
			 vmapM, vmapP, rx, 
			 EToV, Dr, rhsu, invV, 
			 Fmask, Fx, LIFT, Fscale;
#pragma endregion 

/*~~~~~~~~~~~~~~~~~~~~~~ Linear Advection Main Program ~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int main()
{
#pragma region Variable initialization
	//Calling Numerical Solver
	N = 4;
	K = 8;
	Nfp = 1;
	Nfaces = 2;

	Linear_Advection1D();
	cout << "The numerical u for this setup using N = " << N << " and K = " << K;
	cout<< endl << u << endl << endl;

	//Checking Rate of Convergence for N
	K = 16;
	int Nmin = 3;
	int Nmax = 4;
	ConvergN(Nmin, Nmax); //UNCOMMENT for findig error
}



