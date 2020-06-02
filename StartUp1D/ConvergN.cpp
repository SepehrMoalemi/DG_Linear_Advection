#include "Error_Functions.h"
#include "Startup_Functions.h"

void ConvergN(int Nmin, int Nmax)
{
	int indx = 0;
	int Maxindx = (Nmax - Nmin) + 1;
	L2Error_N.resize(Maxindx);

	for (int n = Nmin; n <= Nmax; n++) {
		N = n;

		//Numerical Solution
		Linear_Advection1D();

		//Exact Solution
		Exact1D();

		//Finding Error
		L2Error_N(indx) = L2Error();

		indx++;
	}
	cout << "The variation of the L2Error for N between " << Nmin << " and " << Nmax << " :";
	cout << endl << L2Error_N << endl << endl;
}
