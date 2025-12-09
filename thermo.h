//------------------------------------------------------------------------------------------
#ifndef thermoH
#define thermoH

#include <vector>
#define BZ_THREADSAFE 
#include <blitz/array.h>

extern int def;
using namespace blitz;

//=====================================Class Gibbse=========================================================
class Gibbs
{  
	FILE *fp;
private:
	double polynom (double T1,Array <double,1> &a,Array <double,1> &Gpower);
	Array <double,1> L1,L2,L3,Ltetra,TC1,BMAG1;
	Array <double,2> GHSERfinal,TCfinal,BMAGfinal;
	Array <int,1> iL1,jL1,kL1,lL1,iL2,jL2,kL2,lL2,iL3,jL3,kL3,lL3,mL3,iLtetra,jLtetra,kLtetra,lLtetra,mLtetra,iTC1,jTC1,
	kTC1,lTC1,iBMAG1,jBMAG1,kBMAG1,lBMAG1;
	Array <double,1> derivBMAG1,derivBMAG2;
	Array <double,1> derivTC1,derivTC2;
	double Tc, Bmag; 
	double p,A;
	double n; 
public:
	std::vector <double> dGm1, dGm2;
	Array <double,2> dG_1 , dG_2 , dG_3 , dG_4 ;
	double T1;
	Gibbs();
	Gibbs(const char *);
	const char *datafile;
	
	std::vector <double> m;
	std::vector <double> Dex;
	std::vector <double> Qex;
	
	
  	
	//===================================================================
	
	int calc_spaces(char *);
	void read(double); 
	
	//===================================================================
	virtual void rewriteGHSERp(double,Array <double,1> &, Array <double,4> &, Array <double,3> &); 
	virtual void rewriteL(double, Array <double,1> &, Array <double,6> &, Array <double,5> &);
	virtual void rewriteLtetra(double, Array <double,1> &, Array <double,7> &, Array <double,6> &);
	
	virtual void rewriteGHSERfinal(double, Array <double,3> &, Array <double,1> &Gpower,Array <double,2> &);
	virtual void rewriteLfinal(double, Array <double,5> &,Array <double,1> &, Array <double,4> &);
	virtual void rewriteL_tetra_final(double, Array <double,6> &,Array <double,1> &,Array <double,5> &);
	
	virtual void scanmatrixL(Array <double,4> &, Array <double,1> &, Array <int,1> &, Array <int,1> &, Array <int,1> &, Array <int,1> &);
	virtual void scanmatrixLtetra(Array <double,5> &,Array <double,1> &,Array <int,1> &,Array <int,1> &,Array <int,1> &,Array <int,1> &,Array <int,1> &);
	
	//===================================================================
	double ylogy(double y);
	double logy(double y);	
	
	
	virtual double functionGibbs(double,std::vector <double> &y1, std::vector <double> &y2);
	
	virtual void function_derivation_Gibbs(double,std::vector <double> &y1, std::vector <double> &y2);
	
	virtual void function_derivation_Gibbs2(double T,std::vector <double> &y1, std::vector <double> &y2);
	
	virtual double functionTC(double,std::vector <double> &y1, std::vector <double> &y2);
	
	virtual double functionBMAG(double,std::vector <double> &y1, std::vector <double> &y2);
	
	virtual double magnetic(double,double,double);
	
	virtual std::vector <double> magnetic_derivation(double T,double Tc,double Bmag,Array <double,1> &TC,Array <double,1> &BMAG); 
	
	virtual vector < std::vector<double> > magnetic_derivation2(double T,double TC,double BMAG,Array <double,1> &dTC,Array <double,1> &dBMAG,	\
															   Array <double,1> &dTC1,Array <double,1> &dBMAG1,	\
															   Array <double,2> &dTC2,Array <double,2> &dBMAG2, int kk);
	
	std::vector <double> chemical_potential(std::vector <double> &y1,std::vector <double> &y2);
	
	std::vector < std::vector<double> > chemical_potential2(std::vector <double> &y1,std::vector <double> &y2);
	
	virtual double fun_tc(double Tc);
	virtual double fun_bmag(double Bmag);	
	
};

//---------------------------------------------------------------------------
#endif
