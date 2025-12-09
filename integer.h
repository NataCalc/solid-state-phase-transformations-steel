//---------------------------------------------------------------------------
#ifndef integerH
#define integerH

#include <vector>
#include "phase.hpp"

class integer
{
public:
	interface *funcin;
	interface *funcpot;
	interface *funcin2d;	
	interface *funcpot2d;
	interface *funcderiv;
	interface *funcFeC;
	interface *funcderivFeCrC;
	//we declare index type:
  double (interface::*nrfuncin)(int k, double x, double *y,vector <double> &, vector <double> &);
  void (interface::*nrfuncpot)(int, double, double,double *,double *,std::vector < std::vector<double> > &P);
	
  double (interface::*nrfuncin2d)(double , double ,vector <double> &, vector <double> &);
  void (interface::*nrfuncpot2d)(double *tx, double *ty,vector <double> &, \
								   vector <double> &,int,int, \
								   std::vector < std::vector<double> > &);

   void (interface::*nrfuncderiv)(double *tx, double *ty,vector <double> &, \
								   vector <double> &,int,int, \
								   std::vector < std::vector<double> > &);
	
   void (interface::*nrfuncderivFeCrC)(int m, double xi, double xf,double *yc, \
									 double *ycr,vector <double> &concentration,vector <double> &Ysol, \
									std::vector < std::vector<double> > &P);
	
   void (interface::*nrfuncFeC)(double y, vector <double> &Ysol, std::vector < std::vector<double> > &P);	

public:
 integer() //Konstruktor
  {
  funcin = NULL;
  funcpot = NULL;
  funcin2d=NULL;
  funcpot2d=NULL;
  funcderiv=NULL; 
  funcFeC=NULL;	 
  funcderivFeCrC=NULL;	  
  nrfuncin=0;
  nrfuncpot=0;
  nrfuncin2d=0;
  nrfuncpot2d=0;
  nrfuncderiv=0;	
  nrfuncFeC=0;
  nrfuncderivFeCrC=0;
  }
	
void Eqdifp(double *t1,double *t2,double *t3,double *t4,double xi,double xf,double *Yi,
				int m,int p,int fi,vector <double> &concentration, \
				vector <double> &Ysol,std::vector < std::vector<double> > &P,std::vector < std::vector<double> > &P1);
void Display(double *t1,double *t2,double *t3,double *t4,int m,int p,double xi,double xf);
void Affiche(int m,double xi,double xf,std::vector < std::vector<double> > &P);

void runge_kutta(double *h,double *x, double *y,vector <double> &concentration, \
							  vector <double> &Ysol); 
double equadiff_pc(double *tx,double *ty,double xi,double xf,double yi,int m,int n,int fi,vector <double> &concentration, \
				   vector <double> &Ysol, std::vector < std::vector<double> > &P,std::vector < std::vector<double> > &P1);	

void Affiche_pc(double *tx,double *ty,int n); 	
void Affiche_pot(double *tx,int m, int n, std::vector < std::vector<double> > &P);	
void Affiche_potderiv(double *tx,int m, int n, std::vector < std::vector<double> > &P1);	

void psevdo_FeC(double y, vector <double> &Ysol, std::vector < std::vector<double> > &P);	
void AfficheFeCrC(int m,double xi,double xf,std::vector < std::vector<double> > &P);
	
double Euler_Romberg(int NC, double H, double ER,double X0, double Y0, \
								double *X, double *Y, int *NL, \
								vector <double> &concentration,vector <double> &Ysol, \
								std::vector < std::vector<double> > &P, \
								std::vector < std::vector<double> > &P1);	

double Equadiff_Rung(double *t,double xi,double xf, double yi,int m,int fi, \
								  vector <double> &concentration,vector <double> &Ysol, \
								  std::vector < std::vector<double> > &P, \
								  std::vector < std::vector<double> > &P1);
	
};

//---------------------------------------------------------------------------
#endif
