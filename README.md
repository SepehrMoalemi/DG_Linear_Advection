# DG_Linear_Advection

Usage:
The Main Program is : "LinearAdvection_MAIN.cpp:", under "StarUp1D"

This porgram uses the Eigen library which needs to be included in the path 

Purpose:
Solving the 1D Linear Advection probelem using the STRONG form of the Discontinous Galerkin Method with use of explicit
Runge-Kutta (RK) methods for integration in the temporal dimension.

For calculating the flux, this program uses an upwind scheme (Change ALPHA to 1.0 for Central).
ALPHA = 0.0 //Upwind
ALPHA = 1.0 //Central 

Problem specifications:
// for x ranging [0, 2*pi]
xmin = 0.0 ; xmax = 2.0 * pi  

FinalTime = 10.0 //State after 10s

a = 2.0 * pi     //Advection Speed

NODETOL = 1e-10
