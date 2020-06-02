#pragma once

#include <iostream>
#include <complex>
#include <math.h> //pow() , fabs(), ceil()
#include <cmath> //sqrt(), floor(), sin()
#include <Eigen/Dense> //Using Vectors and Matrices
#include <Eigen/Eigenvalues>  //Using eigen value solver

using namespace std;
using namespace Eigen;

#pragma region Math Constants
	extern const double pi;
	extern const long double tol;
	extern Matrix<double, 5, 1> rk4a, rk4b, rk4c;
#pragma endregion

#pragma region Global Constants
	extern const double ALPHA;
	extern const double xmin;
	extern const double xmax;
	extern const double FinalTime;
	extern const double NODETOL;
#pragma endregion 

#pragma region Global Variables
	extern int N, K,
		Nfp, Np, Nv,
		Nfaces, TotalFaces,
		vmapI, vmapO, mapI, mapO;

	extern Vector2i mapB, vmapB;

	extern VectorXd VX, w, r, P, L2Error_N;

	extern MatrixXi EToE, EToF, nx;

	extern MatrixXd u, u_exact, x, J, V,
			vmapM, vmapP, rx,
			EToV, Dr, rhsu, invV,
			Fmask, Fx, LIFT, Fscale;
#pragma endregion 