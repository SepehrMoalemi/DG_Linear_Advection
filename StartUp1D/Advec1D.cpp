#include "Startup_Functions.h"

//function [u, rhsu] = Advec1D(u, FinalTime, alpha)
void Advec1D(double FinalTime, double alpha)
{
	//Runge - Kutta residual storage
	MatrixXd resu;
	Matrix<double, 5, 1> rk4a, rk4b, rk4c;

	resu.resize(Np, K);
	rk4a << 0.0, -567301805773.0 / 1357537059087.0, -2404267990393.0 / 2016746695238.0, -3550918686646.0 / 2091501179385.0, -1275806237668.0 / 842570457699.0;
	rk4b << 1432997174477.0 / 9575080441755.0, 5161836677717.0 / 13612068292357.0, 1720146321549.0 / 2090206949498.0, 3134564353537.0 / 4481467310338.0, 2277821191437.0 / 14882151754819.0;
	rk4c << 0.0, 1432997174477.0 / 9575080441755.0, 2526269341429.0 / 6820363962896.0, 2006345519317.0 / 3224310063776.0, 2802321613138.0 / 2924317926251.0;

	//Compute time step size
	int Nsteps;

	double timelocal;//Outer time step loop
	double a = 2.0 * pi;//Advection Speed
	double CFL = 0.75;
	double time = 0.0;
	double xmin, dif, dt;

	xmin = fabs(x(0, 0) - x(1, 0));
	for (int i = 0;i < Np;i++) {
		dif = x(0, i) - x(1, i);
		dif = fabs(dif);
		if (dif <= xmin)
			xmin = dif;
	}

	dt = 0.5 * CFL * xmin / (2.0 * pi);
	Nsteps = ceil(FinalTime / dt);
	dt = FinalTime / Nsteps;

	for (int tstep = 1;tstep <= Nsteps;tstep++) {
		for (int INTRK = 0; INTRK < 5;INTRK++) {
			timelocal = time + rk4c(INTRK) * dt;
			AdvecRHS1D(timelocal, a, alpha);
			resu = rk4a(INTRK) * resu + dt * rhsu;
			u = u + rk4b(INTRK) * resu;
		}
		time = time + dt;
	}
}