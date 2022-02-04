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
#define metalts 1673

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
void Update_Temperature_Field(double ***temp, double ***temp_old, int xsize, int ysize, int zsize, double *yc, double timeca, double Gy, double Vy);


// read discrete data from nist dataset, not used in this program
void import_data(double* X, double* Y, double* L, double* D, double* T, int nrows);
#endif //CA_HEAT_CA_HEAT_H
