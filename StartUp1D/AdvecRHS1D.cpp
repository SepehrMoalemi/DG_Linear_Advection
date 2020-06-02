#include "Startup_Functions.h"

//function [rhsu] = AdvecRHS1D(u,time, a, alpha)
void AdvecRHS1D(double time, double a, double alpha)
{
	//Form field differences at faces
	MatrixXd du;
	du.resize(Nfp * Nfaces, K);

	int indx = 0;
	int rowM, rowP, colM, colP, rowN, colN;

	for (int j = 0;j < K;j++)
		for (int i = 0;i < Nfp * Nfaces;i++) {
			rowM = (int)vmapM(indx, 0) % Np;
			rowP = (int)vmapP(indx, 0) % Np;

			colM = (int)vmapM(indx, 0) / Np;
			colP = (int)vmapP(indx, 0) / Np;

			rowN = indx % (Nfaces * Nfp);
			colN = indx / (Nfaces * Nfp);

			//Strong form
			du(i, j) = (u(rowM, colM) - u(rowP, colP)) * (a * nx(rowN, colN) - (1.0 - alpha * fabs(a * nx(rowN, colN)))) / 2.0;

			indx++;
		}

	//Impose boundary condition at x = 0
	double uin = -sin(a * time);

	du(0, 0) = (u(0, 0) - uin) * (a * nx(0, 0) - (1.0 - alpha) * fabs(a * nx(0, 0))) / 2.0;
	du(Nfp * Nfaces - 1, K - 1) = 0.0;

	//Compute right hand sides of the semi - discrete PDE
	rhsu = -a * rx.cwiseProduct(Dr * u) + LIFT * (Fscale.cwiseProduct(du));
}