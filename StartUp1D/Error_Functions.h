#pragma once

#include "Globals.h"

//Exact Solution
void Exact1D(void);

//Error & convergence finder
double L2Error(void);
void ConvergN(int, int);