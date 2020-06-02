#include "Error_Functions.h"

void Exact1D(void)
{
	u_exact.resize(Np, K);

	for (int i = 0;i < Np;i++)
		for (int j = 0;j < K;j++)
			u_exact(i, j) = sin(x(i, j) - 2 * pi * FinalTime);
}