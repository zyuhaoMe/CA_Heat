//
// Created by zyuhao on 2022-01-18.
//

#ifndef CA_HEAT_CA_HEAT_H
#define CA_HEAT_CA_HEAT_H

#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// SI unit
#define PI 3.141592653589

// material property
#define materialDensity 7780.
#define thermalConductivity 20.
#define powerAbsorptivity 0.3
#define specificHeat 523.
#define thermalDiffusivity 4.9152850619571685e-06
#define latentHeat 21e4
#define powerLaser 280.
#define beamRadiusX 42.5e-6 // the 3d-size of heat source is tunable
#define beamRadiusY 42.5e-6
#define beamRadiusZ 42.5e-7
#define metalts 1941 // liquidus temperature

// scanning path parameter
#define xOrigin 0.03791203e-3
#define yOrigin 0.62330961e-3
#define spiralRadius 0.63628482e-3
#define cycleTime 4.04824e-3
#define radiusDifference 0.10194418e-3
#define initialPhase 0.00086551


// define common functions
double* Allocate_1D_Double(int size, double value);
double** Allocate_2D_Double(int size1, int size2, double value);
double*** Allocate_3D_Double(int size1, int size2, int size3, double value);
double**** Allocate_4D_Double(int size1, int size2, int size3, int size4, double value);
void Free_1D_Double(double *a, int size);
void Free_2D_Double(double** a, int size1, int size2);
void Free_3D_Double(double*** a, int size1, int size2, int size3);
void Free_4D_Double(double**** a, int size1, int size2, int size3, int size4);


// define CA functions
double Gaussian_Distribution(double mean, double sigma);
void Seeding(double ***tnuc, double ns, double nb, double Ts, double Tb, double Sigs, double Sigb, int xsize, int ysize, int zsize, double dx);
void Crystal_to_Sample(double phi1, double Phi, double phi2, double *ptCrystal, double *ptSample);
void CrossProduct(double *v1, double *v2, double *v12);
void Subtraction(double *v1, double *v2, double *v12);
double VectorLength(double *v);
double Min(double a, double b);
double Max(double a, double b);
void Sample_to_Crystal(double phi1, double Phi, double phi2, double *ptSample, double *ptCrystal);
void Calculate_New_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct,
                              double rOct, double xc, double yc, double zc, double dx, double *newx, double *newy, double *newz, double *newr);
int Is_in_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct,
                     double rOct, double xc, double yc, double zc);
void Initialization(double ***phi1, double ***Phi, double ***phi2, double ***nucindex, double *xc, double *yc, double *zc,
                    double ***nucx, double ***nucy, double ***nucz, double grain_size, double dx, int xsize, int ysize, int zsize);


// define heat transfer functions
double computeThermalFieldWithGaussianLegendre(
        int numNodes,
        int xsize,
        int ysize,
        int zsize,
        double* xc,
        double* yc,
        double* zc,
        double dtca,
        double timeca,
        double ***temp,
        double phyConst,
        double laserTime);

double summationOfGaussianLegendre(
        double xVal,
        double yVal,
        double zVal,
        double timeca,
        double* abssicaList,
        double* weightList,
        int numNodes);

// given the value of instantaneous time "timeVal", and time of interest "timeca"
// and coordinates of point of interest, evaluate the integrand
double calculateIntegrandVal(double timeVal, double timeca, double xVal, double yVal, double zVal);

// given the value of instantaneous time "timeVal", calculate the absolute position of the moving heat source
// "*xHeat" "*yHeat" "*zHeat", within one individual circle
void calculateHeatSourcePositionSpiral(double timeVal, double *laserPosition);

// calculate the physical constant ahead of the formula
double calculatePhysicalConstant();


// read discrete data from nist dataset, not used in this program
void import_data(double* X, double* Y, double* L, double* D, double* T, int nrows);
#endif //CA_HEAT_CA_HEAT_H
