//---------------------------------------------------------------------------
#ifndef newtH
#define newtH

#include <vector>
#include "phase.hpp"

#define ALFnewt 1.0e-4
#define epsilonnewt 1.0e-9

#define TOLXnewt epsilonnewt
#define FMAXnewt std::max
#define SQRnewt(a) ((a)*(a))
//#define STPMXnewt 1.
#define STPMXnewt 1.0e-2
#define TOLMIN 1.0e-6

class newt
{
private:
  double fmin(std::vector <double> &x,std::vector <double> &concentration);
  void lnsrch(int n, std::vector <double> &xold, double fold, std::vector <double> &g, std::vector <double> &p, std::vector <double> &x,double *f, double stpmax, int *check,std::vector <double> &concentration, int s); 
  void lnsrch_old(int n, std::vector <double> &xold, double fold, std::vector <double> &g, std::vector <double> &p, std::vector <double> &x,double *f, double stpmax, int *check,std::vector <double> &concentration, int s); 
  void echange(double &a,double &b) const;
public:

  interface *func;
  interface *jacobian;
  //we declare index type:
  void (interface::*nrfuncv)(int nn, std::vector <double> &v, std::vector <double> &f,std::vector <double> &concentration);
  void (interface::*jacobfunc)(std::vector <double> &v, std::vector <double> &f1, std::vector < std::vector<double> > &J,std::vector <double> &concentration);

public:
  int nn;
  std::vector <double> fvec;
  std::vector <double> fvec1;
  double correc_STPMXnewt;
  newt() //Konstruktor
  {
  func = NULL;
  jacobian = NULL;
  nn=0;
  nrfuncv=0;
  jacobfunc=0;
  correc_STPMXnewt = 1. ;
  }
  unsigned int GetN() const {return nn;};
  double norme2(std::vector <double> &V);
  void SwapRows(int i1, int i2,std::vector < std::vector<double> > &A);
  void IdentityMatrix(int n, std::vector < std::vector<double> > &A);
  bool LUP(std::vector < std::vector<double> > &A,std::vector < std::vector<double> > &C, std::vector < std::vector<double> > &P);
  bool Inverse(std::vector < std::vector<double> > &A,std::vector < std::vector<double> > &K);
  bool calc(std::vector <double> &x, int *check,std::vector <double> &concentration, int s,bool &conv); //Newton
};

//---------------------------------------------------------------------------
#endif
