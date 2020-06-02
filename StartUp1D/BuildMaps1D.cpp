#include "Startup_Functions.h"

//function [vmapM, vmapP, vmapB, mapB] = BuildMaps1D
void BuildMaps1D(void)
{
	//Number volume nodes consecutively
	int indx = 0;
	MatrixXd nodeids;
	nodeids.resize(Np, K);

	vmapM.resize(Nfaces, K);
	vmapP.resize(Nfaces, K);

	for (int j = 0;j < K;j++)
		for (int i = 0;i < Np;i++) {
			nodeids(i, j) = indx;
			indx++;
		}

	//Find index of face nodes with respect to volume node ordering
	int k2, f2, vidM, vidP, row, col, x1, x2;
	double D;

	for (int k = 0; k < K;k++)
		for (int f = 0;f < Nfaces;f++)
			vmapM(f, k) = nodeids(Fmask(0, f), k);

	for (int k1 = 0; k1 < K;k1++)
		for (int f1 = 0;f1 < Nfaces;f1++) {
			//find neighbor
			k2 = EToE(k1, f1);
			f2 = EToF(k1, f1);

			//find volume node numbers of left and right nodes
			vidM = vmapM(f1, k1);
			vidP = vmapM(f2, k2);

			row = vidM % Np;
			col = vidM / Np;
			x1 = x(row, col);

			row = vidP % Np;
			col = vidP / Np;
			x2 = x(row, col);

			//Compute distance matrix
			D = pow(x1 - x2, 2);
			if (D < NODETOL)
				vmapP(f1, k1) = vidP;
		}

	//Turning vmapP and vmapM into vectors (CAN BE IMPROVED)
	MatrixXd M_Temp;
	M_Temp.resize(Nfaces, K);
	M_Temp = vmapP;
	vmapP.resize(Nfaces * K, 1);

	indx = 0;
	for (int i = 0; i < K;i++)
		for (int j = 0;j < Nfaces;j++) {
			vmapP(indx, 0) = M_Temp(j, i);
			indx++;
		}

	M_Temp = vmapM;
	vmapM.resize(Nfaces * K, 1);

	indx = 0;
	for (int i = 0; i < K;i++)
		for (int j = 0;j < Nfaces;j++) {
			vmapM(indx, 0) = M_Temp(j, i);
			indx++;
		}

	//Create list of boundary nodes
	mapB << 0,
		Nfaces* K - 1;
	vmapB(0) = vmapM(mapB(0));
	vmapB(1) = vmapM(mapB(1));

	//Create specific left (inflow) and right (outflow) maps
	mapI = 0;
	vmapI = 0;
	vmapO = K * Np;
	mapO = K * Nfaces;
}