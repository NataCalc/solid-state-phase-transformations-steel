/*
La classe phase est abstraite et ne peut être utilisée. Il faut utiliser une des classes heritieres
de phase pour definir une nouvelle phase.

Chaque classe heritiere contient la definition de parametres et fonctions additionnels
et, si necessaire, la redefinition de fonctions. Les fonctions pouvant etre redefinies ont le mot
cle virtual.

Fonctions virtuelles
les classes qui heritent directement de "phase" doivent contenir une definiton, meme vide,
des fonctions virtuelles, afin de ne pas etre considerees comme abstraites

le mot cle virtual est facultatif dans toutes les classes heritieres de la classe abstraite "phase"

Ne pas redefinir dans les classes heritieres les parametres membres des classes de base.


Arborescence des classes:
phase (abstraite)
	- matrice
	- sphere
		- sphere_germ_homogene
			- sphere_germ_homogene_elasticite
*/

#ifndef PHASE_HPP
#define PHASE_HPP

#include "thermo.h"


class Gibbs;            //thermodynamyque
class interface;

class distribution;
template <class Type> class distribution_moyenne;
template <class Type> class distribution_euler;

class phase;
class matrice;
class sphere;

class siteprecipitation;
class site_homogene_constant;

class siteprecipitation
{
public:
  siteprecipitation();
  siteprecipitation(phase *);

  phase *precipiteparent;
  double densite_initiale[10];
  double densite_disponible[10];

  virtual void evolution(double)=0;
};

class site_homogene_constant : public siteprecipitation
{
public:
  site_homogene_constant();
  site_homogene_constant(phase *);

  virtual void evolution(double);
};


class distribution
{
public:
  distribution();
  distribution(siteprecipitation *);
  
  phase *phasemere;
  siteprecipitation *siteparent;

  double fluxgermination;

  virtual void evolution(double)=0;
  virtual void germination_d(double)=0;
};


template <class Type>
class distribution_moyenne : public distribution
{
public:
  distribution_moyenne();
  distribution_moyenne(siteprecipitation *);
  ~distribution_moyenne();

  Type *classe0;
  
  virtual void evolution(double);
  virtual void germination_d(double);
};


template <class Type>
class distribution_euler : public distribution
{
public:
  distribution_euler();
  distribution_euler(siteprecipitation *);
  ~distribution_euler();
  double deltaR;
  double Rcritique;
  
  Type *classe[10];
  
  virtual void evolution(double);
  virtual void germination_d(double);
};


class interface
{
public:
  interface();
  interface(phase *,phase *); //Constructeur
  phase *ph1;
  phase *ph2;
  double xeq1[10];
  double xeq2[10];
  double param_cinetique;
  //int nalfa, nbetta;  		

void SB(int,std::vector <double> &Y,std::vector <double> &F,std::vector <double> &concentration);//system_binaire
void JB(std::vector <double> &Y,std::vector <double> &F,std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration);//Jacobian_binaire
void KTB(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_binaire
void JKB(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void KTB_prof(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_binaire
void JKB_prof(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
	
void KTB1(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_binaire
void JKB1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration);//Jacobian_kinetic_binaire

	
void ST(int,std::vector <double> &Y,std::vector <double> &F,std::vector <double> &concentration);//system_ternaire
void JT(std::vector <double> &Y, std::vector <double> &F, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_ternaire
void KTT(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_sphere(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof1(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof2(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof3(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof4(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof5(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof6(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void KTT_prof7(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire

	void KTT_prof_sphere(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void JKT(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_sphere(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof2(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof3(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof4(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof5(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof6(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void JKT_prof7(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire

	void JKT_prof_sphere(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire

void KTT1(int, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);//kinetic_transitions_ternaire
void JKT1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration);//Jacobian_kinetic_binaire
void Jacobian3D(std::vector <double> &Y, std::vector <double> &F, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration);
std::vector <double> Zenerplane(std::vector <double> &Y,std::vector <double> &concentration);
std::vector <double> Zenerplane_cin(std::vector <double> &Y,std::vector <double> &concentration);
std::vector <double> Zenerplane2d(std::vector <double> &Y,std::vector <double> &concentration);
std::vector <double> Zenerplane2d_prof(std::vector <double> &Y,std::vector <double> &concentration);
	
	
void SQ4(int,std::vector <double> &Y,std::vector <double> &F,std::vector <double> &concentration);//system_quaternaire
void JQ4(std::vector <double> &Y, std::vector <double> &F, std::vector < std::vector<double> > &Jacobian4D,std::vector <double> &concentration);//Jacobian_quaternaire
void KT4(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);
void JK4(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian4D,std::vector <double> &concentration);
void KT4_1(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);
void JK4_1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian4D,std::vector <double> &concentration);
void Jacobian4D(std::vector <double> &Y, std::vector <double> &F, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration);

void SQ5(int,std::vector <double> &Y,std::vector <double> &F,std::vector <double> &concentration);
void JQ5(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration);	
	
double phase1_binaire(double c,double fe,double cr, double mo, double co, double al);
double phase2_binaire(double c,double fe,double cr, double mo, double co, double al);

		
double phase1_ternaire(double c,double fe,double cr, double mo, double co, double al);
double phase2_ternaire(double c,double fe,double cr, double mo, double co, double al);
std::vector < std::vector<double> > JacobianQhull3Dalfa(double c,double fe,double cr, double mo, double co, double al);
std::vector < std::vector<double> > JacobianQhull3Dbetta(double c,double fe,double cr, double mo, double co, double al);

	
double phase1_quaternaire(double c,double fe,double cr, double mo, double co, double al);
double phase2_quaternaire(double c,double fe,double cr, double mo, double co, double al);
std::vector < std::vector<double> > JacobianQhull4Dalfa(double c,double fe,double cr, double mo, double co, double al);
std::vector < std::vector<double> > JacobianQhull4Dbetta(double c,double fe,double cr, double mo, double co, double al);

double phase1_quinary(double c,double fe,double cr, double mo, double co, double al);
double phase2_quinary(double c,double fe,double cr, double mo, double co, double al);
	
std::vector <double> calc_diffusion_alpha(std::vector <double> &concentration);
std::vector <double> calc_diffusion_betta(std::vector <double> &concentration);
	
vector <double> integration(int n,vector <double> &Y,vector <double> &concentration);	
//std::vector <double> FeCfrmolaire(int n,vector <double> &Y);	
double funcFeC(double x, double y,vector <double> &concentration, vector <double> &Ysol);
std::vector <double> derivFeC(double y);
void funconeFeC(double *x, double *y,vector <double> &concentration, \
			 vector <double> &Ysol, int m, int n, \
			 std::vector < std::vector<double> > &P);
void funcderivmuFeC(double *x, double *y,vector <double> &concentration, \
								   vector <double> &Ysol, int m, int n, \
					std::vector < std::vector<double> > &P);
void funcderivmuFeC_(double y,vector <double> &Ysol, \
									std::vector < std::vector<double> > &P);	
	
std::vector < std::vector<double> > derivFeCrC(double yc, double ycr);	
std::vector < std::vector<double> > derivFeCrC_U(double yc, double ycr);	

void funconeFeCrC(int m, double xi, double xf,double *yc,double *ycr,std::vector < std::vector<double> > &P);

double funcFeCrC(int k, double x, double *y,vector <double> &concentration, vector <double> &Ysol);
void funcderivmuFeCrC(int m, double xi, double xf,double *yc,double *ycr,vector <double> &concentration,vector <double> &Ysol, \
					  std::vector < std::vector<double> > &P);

void JKB_prof_U(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D, \
					  std::vector <double> &concentration);
void KTB_prof_U(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration);

double deriv_mu(double ya,double yb);
std::vector <double> derivFeC_U(double y);	
std::vector <double> deriv_mu_FeCCr(double yc,double ycr,double x,vector <double> &concentration);
std::vector <double> deriv_mu_FeCCr1(double yc,double ycr,double x,vector <double> &concentration);

std::vector <double> deriv_mu_FeCCr_X(double yc,double ycr,double x,vector <double> &concentration);
std::vector <double> deriv_mu_FeCCr_KTT_prof5(double yc,double ycr,double x,vector <double> &concentration);



};


class phase
{
public:
  phase();
  phase(distribution *);
  ~phase();

  siteprecipitation *site[10];
  distribution *distri[10];
  phase *produit[10];

  phase *voisin[10];

  phase *phasemere;
  distribution *distributionmere;

  interface *frontiere[10];
  Gibbs *Thermo;
	

  double xintm[10];	//Concentration d'equilibre a l'interface dans la matrice
  double xintp[10];     //Concentration d'equilibre a l'interface dans le precipite
  double xini[10];	//Concentration initiale
  double xmoy[10];	//Concentration moyenne
  double omega[10];	//Sursaturation
  double D0[10];	//Coefficient de diffusion
  double E0[10];
  double D[10];
  double A[10];         //Coefficient cinetique
  double (*potentiel[10]) (double,double); //potentiels chimiques

  double Taille[10];
  double dTaille[10];

  double fraction_apparente;	//densite apparente * volume individuel des precipites
  double fraction;	//Fraction volumique reelle
  double volmolaire;
  double parametre_maille;
  double diffusion;
  double module_Young;
  double coefficient_Poisson;
  
  void ini(double);
  void calc_diffusion(double);
  
  virtual void calc_composition()=0;
  virtual void calc_fraction()=0;
  virtual void evolution(double)=0;
  virtual void calc_germination(double)=0;
  virtual void calc_croissance(double)=0;
};


class matrice : public phase
{
public:
  //Constructeur
  matrice();
  //Fonctions
  virtual void calc_composition();
  virtual void calc_fraction();
  virtual void evolution(double);
  virtual void calc_germination(double);
  virtual void calc_croissance(double);
};


class sphere : public phase
{
public:
  //Constructeurs
  sphere();
  sphere(distribution *);

  //Donnee d'entree
  double eninterface;
  double K0;//Constante d'equilibre
  double Q;
  double Keq;
  //Parametres microstructuraux
  //  double rayon;
  double densite_apparente;
  double densite;
  //Parametre intermediaire
  double effet_Gibbs_Thomson;		
  //Parametre d'evolution
  // double drayon;//vitesse de croissance en m.s-1
  //Fonctions
  virtual void calc_composition();
  virtual void calc_fraction();		
  virtual void evolution(double);
  virtual void calc_germination(double);		
  virtual void calc_croissance(double);
  void copie(sphere *);
};


class sphere_voisin : public sphere
{
public:
  //Constructeurs
  sphere_voisin();
  sphere_voisin(phase *);

  double drayon2;
  double rayon2;

  virtual void evolution(double);
  virtual void calc_croissance(double);
};

class sphere_germ_homogene : public sphere
{
public:
  //Constructeurs
  sphere_germ_homogene();
  sphere_germ_homogene(distribution *);

  //Donnees d'entree
  double eninterfacegerm;
  double vatomique;
  //Parametre de microstructure
  double nsites;
  //Parametres intermediaires pour le calcul du flux de germination
  double forcemotrice;//par unite de volume
  double zeldovich;
  double transport;
  double barriere;
  double rayoncritique;
  //Parametre d'evolution
  double fluxgermination;
  //Fonctions		
  virtual void evolution(double);
  virtual void calc_germination(double);
};

class sphere_germ_homogene_elasticite : public sphere_germ_homogene
{
public:
  //Constructeurs
  sphere_germ_homogene_elasticite();
  sphere_germ_homogene_elasticite(distribution *);

  //Donnee d'entree
  double misfit;
  //Parametre intermediaire
  double effet_energie_elastique;
  //Fonctions
  virtual void calc_germination(double);
};

//------------------Classe Distribution Moyenne---------------------------------------

template <class Type>
distribution_moyenne<Type>::distribution_moyenne() : distribution()
{
  classe0=0;
}

template <class Type>
distribution_moyenne<Type>::distribution_moyenne(siteprecipitation *site) : distribution(site)
{
  classe0=new Type(this);
}

template <class Type>
distribution_moyenne<Type>::~distribution_moyenne()
{
  if(classe0) delete classe0;
}

template <class Type>
void distribution_moyenne<Type>::evolution(double dt)
{
  int i=1;
  siteparent -> densite_disponible[1] -= fluxgermination;
  classe0 -> densite_apparente += fluxgermination;

  classe0 -> Taille[1] += classe0 -> dTaille[1] * dt;
}

template <class Type>
void distribution_moyenne<Type>::germination_d(double T) 
{
  classe0->calc_germination(T);
}

//------------------Classe Distribution Eulerienne-------------------------------------

template <class Type>
distribution_euler<Type>::distribution_euler() : distribution()
{
  deltaR=0;
  Rcritique=0;
  for(int i=0;i<10;i++)
    classe[i]=0;
}

template <class Type>
distribution_euler<Type>::distribution_euler(siteprecipitation *site) : distribution(site)
{
  deltaR=0;
  Rcritique=0;
  for(int i=1;i<=10;i++)
    {
      classe[i]=new Type(site);
    }
}

template <class Type>
distribution_euler<Type>::~distribution_euler()
{
  for(int i=0;i<10;i++)
    if(classe[i])
    {
      delete classe[i];
      classe[i]=0;
    }
}

template <class Type>
void distribution_euler<Type>::evolution(double dt)
{
  //  for(int i=1;i<=10;i++)
  //    {
  //      distri[i]->densite+=dt/deltaR*distri[i-1]->drayon*distri[i-1]->densite;
  //    }
  double rcritique;
  rcritique=classe[1]->rayoncritique;
}

template <class Type>
void distribution_euler<Type>::germination_d(double T) 
{
  classe[1]->calc_germination(T);
}

#endif
