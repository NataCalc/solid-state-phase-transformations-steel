/*
 *  integration.h
 *  
 *
 *  Created by nat on 06/07/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef integrationH
#define integrationH

#include <cmath>
#include <iostream>
#include <vector>
#include "phase.hpp"


//#define EPS 1.0e-5
//#define JMAX 17


//extern float Diff[10];
//extern float taille_gradient;
//extern double correc_D;
//extern float CIN;

using namespace std;



double Wang(double r);
double dWang(double r);

/*
double trapzd(double (*func)(double,double,double,double,vector <double>&,double), double a, double b, int n, \
			  double delta, double Da,double Db,vector <double>&Y,double vrf);
double qtrap(double (*func)(double,double,double,double,vector <double>&,double), double a, double b, int n, \
			 double delta, double Da,double Db,vector <double>&Y,double vrf);

double simpson(double (*func)(double,double,double,double,vector <double>&,double),int n, double a, double b, \
			   double delta, double Da,double Db,vector <double>&Y,double vrf);

double qgaus(double (*func)(double,double,double,double,vector <double>&, double),int n, double a, double b, \
			 double delta, double Da,double Db,vector <double>&Y, double vrf);



double Wangh(double Da, double Db, double r);
double Wanghder(double Da, double Db, double r);
double fp(double x, double y);
void Affiche(double *t,int m,double xi,double xf);
void Equadif1(double *t,double xi,double xf, double yi,int m,int fi,vector <double> &concentration, \
			  vector <double> &Ysol,interface *iface,double (interface::*fp)(double , double ,vector <double> &, vector <double> &));

void equadiffn(double xi,double xf,double *yi,int m,int n,int fi,double (*fp)(double, double*));
double fp1(double x, double *y);


*/



void runge_kutta(double *h,double *x, double *y,vector <double> &concentration, \
				 vector <double> &Ysol,interface *iface,double (interface::*fp)(double , double ,vector <double> &, vector <double> &)); 

double equadiff_pc(double *tx,double *ty,double xi,double xf,double yi,int m,int n,int fi,vector <double> &concentration, \
				 vector <double> &Ysol,interface *iface,double (interface::*fp)(double , double ,vector <double> &, vector <double> &), \
				 std::vector < std::vector<double> > &P, \
				 void (interface::*fp1)(double *tx, double *ty,vector <double> &, \
										vector <double> &,int,int, \
										std::vector < std::vector<double> > &));

void Affiche_pc(double *tx,double *ty,int n);
void Affiche_pot(double *tx,int m, int n, std::vector < std::vector<double> > &P);
void Affichefinal_1(int n, vector <double> &Y, \
					int m, vector <double> &Yder, \
					double G1, double G2, \
					double dGm2aC, double dGm2aVA, double dGm2bC,double dGm2bVA, \
					double ddGm2aCC, double ddGm2aCVA, double ddGm2aVAVA, double ddGm2aVAC, \
					double ddGm2bCC, double ddGm2bCVA, double ddGm2bVAVA, double ddGm2bVAC, \
					double x);


//double fp2(int k, double x, double *y);
//double fp3(int k, double x, double *y);
void Display(double *t1,double *t2,double *t3,double *t4,int m,int p,double xi,double xf);
void Eqdifp(double *t1,double *t2,double *t3,double *t4,double xi,double xf,double *Yi, \
			int m,int p,int fi,vector <double> &concentration, \
			vector <double> &Ysol,interface *iface,double (interface::*fp)(int k, double x, double *y,vector <double> &, vector <double> &), \
			std::vector < std::vector<double> > &P,void (interface::*fp1)(int, double, double,double *,double *, \
																		  std::vector < std::vector<double> > &P));
void Affiche_pot3d(int m,double xi,double xf,std::vector < std::vector<double> > &P);
void RK4n(int n, double x, double *y, double h, double *y1, \
		  vector <double> &concentration,vector <double> &Ysol, \
		  interface *iface,void (interface::*fp)(double x, double *y,double *yp,vector <double> &, vector <double> &)); 

#endif