#include "Error_Functions.h"

double L2Error(void)
{
	double max_Error;
	max_Error = fabs(u_exact(0, 0) - u(0, 0));
	for (int i = 0;i < Np;i++)
		for (int j = 0;j < K;j++) {
			if (fabs(u_exact(i, j) - u(i, j)) > max_Error)
				max_Error = fabs(u_exact(i, j) - u(i, j));
		}
	return max_Error;
}