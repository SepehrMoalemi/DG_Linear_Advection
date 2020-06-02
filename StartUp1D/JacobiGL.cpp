#include "Startup_Functions.h"

//function [x] = JacobiGL(alpha,beta,N)
void JacobiGL(double alpha, double beta)
{
	x.resize(Np, 1);

	if (N == 1) {
		x(0, 0) = -1.0;
		x(1, 0) = 1.0;
	}
	else {
		JacobiGQ(alpha + 1.0, beta + 1.0);

		//Adding end points to x
		MatrixXd x_Temp;
		x_Temp.resize(N - 1, 1);
		x_Temp = x;
		x.resize(Np, 1);

		for (int i = 1;i < N;i++)
			x(i, 0) = x_Temp(i - 1, 0);

		x(0, 0) = -1.0;
		x(N, 0) = 1.0;
	}
}