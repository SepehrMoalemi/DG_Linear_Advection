#include "Startup_Functions.h"

void Linear_Advection1D(void)
{
	Np = N + 1;
	Nv = K + 1;
	TotalFaces = Nfaces * K;

#pragma region Allocating size
	EToE.resize(K, Nfaces);
	EToF.resize(K, Nfaces);
	nx.resize(Nfaces * Nfp, K);
	u.resize(Np, K);
	rx.resize(Np, K);
	EToV.resize(K, 2);
	Dr.resize(Np, Np);
	rhsu.resize(Np, K);
	invV.resize(Np, Np);
	Fmask.resize(Nfp, Nfaces);
	Fx.resize(Nfaces * Nfp, K);
	LIFT.resize(Np, Nfaces * Nfp);
	Fscale.resize(Nfaces * Nfp, K);
#pragma endregion

#pragma region  MeshGen1D.m
	//Generate simple mesh
	VX.resize(Nv);

	for (int i = 0; i < Nv; i++)
		VX(i) = (xmax - xmin) * (double)i / ((double)(Nv - 1)) + xmin;

	for (int i = 0; i < K; i++) {
		EToV(i, 0) = (double)i;
		EToV(i, 1) = (double)i + 1.0;
	}
#pragma endregion

#pragma region StartUp1D.m
	//Compute basic Legendre Gauss Lobatto grid
	JacobiGL(0.0, 0.0);
	r.resize(Np);
	r = x.col(0);

	//Legendre polynomials obtained with alpha = beta = 0
	MatrixXd V_1D;
	V_1D.resize(Np, N + 1);

	for (int j = 1; j <= N + 1;j++) {
		JacobiP(0.0, 0.0, j - 1);
		V_1D.col(j - 1) = P;
	}

	//Build reference element matrices
	V.resize(Np, N + 1);
	V = V_1D;

	invV = V.inverse();

	double k;
	MatrixXd DVr, dp;

	DVr.resize(Np, N + 1);
	dp = MatrixXd::Zero(Np, 1);

	for (int i = 0;i <= N;i++) {
		if (i != 0) {
			k = (double)i;

			JacobiP(0.0 + 1.0, 0.0 + 1.0, i - 1);
			dp = sqrt(k * (k + 0.0 + 0.0 + 1.0)) * P;
		}
		DVr.col(i) = dp;
	}

	Dr = DVr * invV;

	//Create surface integral terms
	MatrixXd Emat;
	Emat = MatrixXd::Zero(Np, Nfaces * Nfp);

	Emat(0, 0) = 1.0;
	Emat(Np - 1, 1) = 1.0;

	LIFT = V * (V.transpose() * Emat);

	//Build coordinates of all the nodes
	MatrixXd M_ones;
	M_ones = MatrixXd::Ones(Np, 1);

	VectorXd va, vb;

	va.resize(K);
	vb.resize(K);
	x.resize(Np, K);

	for (int i = 0; i < K; i++) {
		va(i) = VX(i);
		vb(i) = VX(i + 1);
	}

	x = M_ones * va.transpose() + 0.5 * (r + M_ones) * (vb.transpose() - va.transpose());

	//Calculate geometric factors
	J.resize(Np, K);
	J = Dr * x;
	for (int i = 0;i < Np;i++)
		for (int j = 0;j < K;j++)
			rx(i, j) = 1.0 / J(i, j);

	//Physical coordinates of the edge nodes
	double fmask1, fmask2;
	for (int i = 0; i < Np;i++) {
		if (fabs(r(i) + 1.0) < NODETOL)
			fmask1 = double(i);
		if (fabs(r(i) - 1.0) < NODETOL)
			fmask2 = double(i);
	}
	Fmask << fmask1, fmask2;
	Fx.row(0) = x.row((int)fmask1);
	Fx.row(1) = x.row((int)fmask2);

	//Build surface normals and inverse metric at surface
	nx = MatrixXi::Zero(Nfaces * Nfp, K);
	for (int i = 0; i < K; i++) {
		nx(0, i) = -1;
		nx(1, i) = 1;
	}
	for (int i = 0; i < K; i++) {
		Fscale(0, i) = 1 / J((int)fmask1, i);
		Fscale(1, i) = 1 / J((int)fmask2, i);
	}

	//Build connectivity matrix
	Connect1D();

	//Build connectivity maps
	BuildMaps1D();
#pragma endregion

#pragma region AdvecDriver1D.m
	//Set initial conditions
	u.resize(Np, K);

	for (int i = 0;i < Np;i++)
		for (int j = 0;j < K;j++)
			u(i, j) = sin(x(i, j));

	//Solve Problem
	Advec1D(FinalTime, ALPHA);

#pragma endregion
}