#include "Startup_Functions.h"

//function [EToE, EToF] = Connect1D(EToV)
void Connect1D(void)
{
	int sk = 0;
	Vector2d vn;
	MatrixXd SpFToV, M_Iden;
	M_Iden.resize(TotalFaces, TotalFaces);

	vn << 0,
		1;

	//Build global face to node sparse array
	SpFToV = MatrixXd::Zero(TotalFaces, Nv);
	M_Iden = MatrixXd::Identity(TotalFaces, TotalFaces);

	for (int k = 0; k < K; k++)
		for (int face = 0; face < Nfaces;face++) {
			SpFToV(sk, EToV(k, vn(face))) = 1;
			sk++;
		}

	//Build global face to global face sparse array
	M_Iden = SpFToV * (SpFToV.transpose()) - M_Iden;
	SpFToV.resize(TotalFaces, TotalFaces);
	SpFToV = M_Iden;

	//Find complete face to face connections
	int counter = 0;
	VectorXd faces1, faces2;
	faces1.resize(TotalFaces - 2);
	faces2.resize(TotalFaces - 2);

	for (int j = 0; j < TotalFaces; j++)
		for (int i = 0; i < TotalFaces; i++)
			if (SpFToV(i, j) == 1.0) {
				faces1(counter) = (double)i;
				faces2(counter) = (double)j;
				counter++;
			}

	//Convert face global number to element and face numbers
	VectorXi element1, element2, face1, face2;

	element1.resize(TotalFaces - 2);
	element2.resize(TotalFaces - 2);
	face1.resize(TotalFaces - 2);
	face2.resize(TotalFaces - 2);

	for (int i = 0; i < TotalFaces - 2;i++) {
		element1(i) = floor(faces1(i) / (double)Nfaces);
		element2(i) = floor(faces2(i) / (double)Nfaces);
		face1(i) = (int)faces1(i) % Nfaces;
		face2(i) = (int)faces2(i) % Nfaces;
	}

	//Rearrange into Nelements x Nfaces sized arrays
	for (int i = 0; i < K;i++) {
		EToE(i, 0) = i;
		EToE(i, 1) = i;
		EToF(i, 0) = 0;
		EToF(i, 1) = 1;
	}

	for (int i = 0;i < TotalFaces - 2;i++) {
		EToE(element1(i), face1(i)) = element2(i);
		EToF(element1(i), face1(i)) = face2(i);
	}
}