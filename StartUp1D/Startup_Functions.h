#pragma once

#include "Globals.h"

//Startup1D
void JacobiGL(double, double);
void JacobiGQ(double, double);
void JacobiP(double, double, int);
void Connect1D(void);
void BuildMaps1D(void);

//Advection
void Advec1D(double, double);
void AdvecRHS1D(double, double, double);

//Numerical Solution
void Linear_Advection1D(void);

