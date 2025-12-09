//------------------------------------------------------------------------------------------
#ifndef simplexH
#define simplexH

#include "phase.hpp"

#define ZEPS 1.e-10
#define tol  1.e-7

class simplex
{
public:
   interface *func;
   double (interface::*nrfuncv)(int nn, std::vector <double> &Y,std::vector <double> &concentration);
public:
  simplex()
    {
     func = NULL;
     nrfuncv=0;
    }
int nn;
void amoeba(int nn,std::vector < std::vector<double> > &simplex,std::vector <double> &concentration,std::vector <double> &Y);
virtual void evaluate_simplex(int nn,std::vector < std::vector<double> > &simplex,std::vector <double> &concentration, std::vector <double> &fx);
virtual void simplex_extremes(int nn,std::vector <double> &fx,int *ihi, int *ilo, int *ihni);
virtual void simplex_bearings(int nn,std::vector < std::vector<double> > &simplex,std::vector <double> &midpoint,int ihi);
virtual int check_tol(double fmax, double fmin, double ftol);
virtual double update_simplex(int nn,int ihi,std::vector < std::vector<double> > &simplex, std::vector <double> &midpoint,std::vector <double> &next,double scale, std::vector <double> &concentration);
virtual double update_stretch(int nn,std::vector <double> &next,std::vector <double> &midpoint,std::vector <double> &next_e,double scale, std::vector <double> &concentration);
virtual void contract_simplex(int nn, int ilo,std::vector < std::vector<double> > &simplex,std::vector <double> &fx, std::vector <double> &oncentration);
virtual int control(int nn, std::vector <double> fx);

};

//---------------------------------------------------------------------------
#endif
