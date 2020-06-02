#include "Startup_Functions.h"

//function [P] = JacobiP(x,alpha,beta,N)
void JacobiP(double alpha, double beta, int NN) {
	double cons = alpha + beta;
	double Gamma0;

	Gamma0 = pow(2.0, -(cons + 1.0)) * tgamma(cons + 2.0) / (tgamma(alpha + 1.0) * tgamma(beta + 1.0));

	//Initializing the Jacobi polynomial
	MatrixXd PL;
	PL.resize(Np, NN + 1);

	for (int i = 0; i < Np; i++)
		PL(i, 0) = sqrt(Gamma0);

	if (NN == 0) {
		P.resize(Np);
		P = PL.col(0);
		return;
	}

	VectorXd Gamma1, cons2;
	Gamma1.resize(Np);
	cons2 = VectorXd::Ones(Np);

	cons2 = (alpha - beta) * cons2;
	Gamma1 = sqrt((cons + 3.0) / ((alpha + 1.0) * (beta + 1.0))) * ((cons + 2.0) * r + cons2);

	PL.col(1) = (1.0 / 2.0) * sqrt(Gamma0) * Gamma1;

	if (NN == 1) {
		P.resize(Np);
		P = PL.col(1);
		return;
	}

	//Recurrence relation
	double a_n, a_n1, b_n;
	double j, i;

	MatrixXd DD;
	VectorXd EE;

	DD.resize(Np, Np);
	EE.resize(Np);

	for (int k = 1; k < NN; k++) {
		DD = r * PL.col(k).transpose();
		EE = DD.diagonal();

		i = (double)k;
		j = (double)i + 1.0;

		a_n = 2.0 / (2.0 * i + cons) * sqrt(i * (i + cons) * (i + alpha) * (i + beta) / ((2.0 * i + cons - 1.0) * (2.0 * i + cons + 1.0)));
		a_n1 = 2.0 / (2.0 * j + cons) * sqrt(j * (j + cons) * (j + alpha) * (j + beta) / ((2.0 * j + cons - 1.0) * (2.0 * j + cons + 1.0)));
		b_n = -(pow(alpha, 2.0) - pow(beta, 2.0)) / ((2.0 * i + cons) * (2.0 * i + cons + 2.0));

		PL.col(k + 1) = 1.0 / a_n1 * (EE - b_n * PL.col(k) - a_n * PL.col(k - 1));
	}

	P.resize(Np);
	P = PL.col(NN);
}