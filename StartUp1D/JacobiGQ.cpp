#include "Startup_Functions.h"

//function [x,w] = JacobiGQ(alpha,beta,N)
void JacobiGQ(double alpha, double beta)
{

	int NN = N - 2;
	double cons = alpha + beta;
	x.resize(NN + 1, 1);

	if (NN == 0) {
		x(0, 0) = (alpha - beta) / (cons + 2.0);
		w(0) = 2.0;
	}
	else {
		VectorXd n, n0, n1, bn;
		n.resize(NN);
		n0.resize(NN + 1);
		n1.resize(NN);
		bn.resize(NN + 1);

		V.resize(NN + 1, NN + 1);
		J.resize(NN + 1, NN + 1);

		J = MatrixXd::Zero(NN + 1, NN + 1);

		//Initialize n
		n0(0) = 0.0;
		for (int i = 0;i < NN; i++) {
			n(i) = (double)i + 1.0;
			n0(i + 1) = 2.0 * ((double)i + 1.0);
			n1(i) = 2.0 * ((double)i + 1.0);
		}

		//Initialize bn
		for (int i = 0;i <= NN;i++)
			bn(i) = (pow(alpha, 2.0) - pow(beta, 2.0)) / ((n0(i) + cons) * (n0(i) + cons + 2.0));

		//Initialize J
		for (int i = 0; i <= NN; i++) {
			J(i, i) = (-1.0 / 2.0) * bn(i);

			if (i != NN) {
				J(i, i + 1) = 2.0 / (n1(i) + cons) * sqrt(n(i) * (n(i) + cons) * (n(i) + alpha) * (n(i) + beta) / ((n1(i) + cons - 1.0) * (n1(i) + cons + 1.0)));
				J(i + 1, i) = J(i, i + 1);
			}
		}

		if (cons < tol)
			J(0, 0) = 0.0;

		//Find Eigen Values and Vectors
		EigenSolver<MatrixXd> es(J);

		VectorXcd eigenVal;
		MatrixXcd eigenVec;

		eigenVal = es.eigenvalues();
		eigenVec = es.eigenvectors();

		for (int i = 0;i < NN + 1;i++) {
			x(i, 0) = real(eigenVal(i));

			for (int j = 0;j < NN + 1;j++)
				V(i, j) = real(eigenVec(i, j));
		}

		//Sorting of eignevalues
		int min;
		double temp;
		MatrixXd tempVec;
		tempVec.resize(NN + 1, 1);

		for (int i = 0; i < NN; i++) {
			min = i;

			for (int j = i + 1; j < NN + 1; j++)
				if (x(j, 0) < x(min, 0))
					min = j;

			temp = x(i, 0);
			x(i, 0) = x(min, 0);
			x(min, 0) = temp;

			tempVec = V.col(i);
			V.col(i) = V.col(min);
			V.col(min) = tempVec;
		}

		//Compute Quadrature Points
		w.resize(NN + 1);
		for (int i = 0;i < NN + 1;i++)
			w(i) = pow(V(0, i), 2.0) * pow(2.0, cons + 1.0) * tgamma(alpha + 1.0) * tgamma(beta + 1.0) / tgamma(cons + 2.0);
	}
}