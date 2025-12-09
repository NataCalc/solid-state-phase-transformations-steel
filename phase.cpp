#include "phase.hpp"
#include "donnees.h"

//extern float Diff[10];
//extern float taille_gradient;

//------------------Classe Distribution-----------------------------------------------

distribution::distribution()
{
  phasemere=0;
  siteparent=0;
  fluxgermination=0;
}

distribution::distribution(siteprecipitation *site) : phasemere(site->precipiteparent) , siteparent(site)
{
  int i=1;
  while(site->precipiteparent->distri[i]!=NULL) i++;
  site->precipiteparent->distri[i]=this;
}


//------------------Classe siteprecipitation------------------------------------------

siteprecipitation::siteprecipitation()
{
  precipiteparent=0;
  for(int i=0;i<10;i++)
  {
    densite_initiale[i]=0;
    densite_disponible[i]=0;
  }
}

siteprecipitation::siteprecipitation(phase *phmere) : precipiteparent(phmere)
{
  int i=1;
  while(precipiteparent->site[i]!=NULL) i++;
  precipiteparent->site[i]=this;
}
//------------------Classe site_homogene_constant-------------------------------------

site_homogene_constant::site_homogene_constant() : siteprecipitation() {}

site_homogene_constant::site_homogene_constant(phase *phmere) : siteprecipitation(phmere) {}

void site_homogene_constant::evolution(double dt) {}




//------------------Classe Interface--------------------------------------------------- 
interface::interface()
{
  ph1=ph2=0;
  //nalfa=nbetta=0;
  param_cinetique=0;
  for(int i=0;i<10;i++)
  {
    xeq1[i]=0;
    xeq2[i]=0;
  }
}

interface::interface(phase *phase1,phase *phase2) : ph1(phase1) , ph2(phase2) 
{
  int i=1;
  while(ph2->frontiere[i]!=NULL) i++;
  ph2->frontiere[i]=this;
  ph2->frontiere[i+1]=NULL;
 
  //nalfa=nbetta=0;
  param_cinetique=0;
  for(int i=0;i<10;i++)
  {
    xeq1[i]=0;
    xeq2[i]=0;
  }
}

//------------------Classe phase-------------------------------------------------------

phase::phase()
{
  phasemere=0;
  distributionmere=0;
  Thermo=0;
  for(int i=0;i<10;i++)
  {
    site[i]=0;
    distri[i]=0;
    produit[i]=0;
    voisin[i]=0;
    frontiere[i]=0;
    xintm[i]=0;	//Concentration d'equilibre a l'interface dans la matrice
    xintp[i]=0;     //Concentration d'equilibre a l'interface dans le precipite
    xini[i]=0;	//Concentration initiale
    xmoy[i]=0;	//Concentration moyenne
    omega[i]=0;	//Sursaturation
    D0[i]=0;	//Coefficient de diffusion
    E0[i]=0;
    D[i]=0;
    A[i]=0;         //Coefficient cinetique
    Taille[i]=0;
    dTaille[i]=0;
    potentiel[i]=0;
  }
  fraction_apparente=0;	//densite apparente * volume individuel des precipites
  fraction=0;	//Fraction volumique reelle
  volmolaire=0;
  parametre_maille=0;
  diffusion=0;
  module_Young=0;
  coefficient_Poisson=0;
}

phase::phase(distribution *distrimere)
{
  phasemere=0;
  distributionmere=0;
  Thermo=0;
  for(int i=0;i<10;i++)
  {
    site[i]=0;
    distri[i]=0;
    produit[i]=0;
    voisin[i]=0;
    frontiere[i]=0;
    xintm[i]=0;	//Concentration d'equilibre a l'interface dans la matrice
    xintp[i]=0;     //Concentration d'equilibre a l'interface dans le precipite
    xini[i]=0;	//Concentration initiale
    xmoy[i]=0;	//Concentration moyenne
    omega[i]=0;	//Sursaturation
    D0[i]=0;	//Coefficient de diffusion
    E0[i]=0;
    D[i]=0;
    A[i]=0;         //Coefficient cinetique
    Taille[i]=0;
    dTaille[i]=0;
    potentiel[i]=0;
  }
  fraction_apparente=0;	//densite apparente * volume individuel des precipites
  fraction=0;	//Fraction volumique reelle
  volmolaire=0;
  parametre_maille=0;
  diffusion=0;
  module_Young=0;
  coefficient_Poisson=0;
  distributionmere=distrimere;
  phase *phm=distributionmere->phasemere;
  phasemere=phm;
  int i=1; 
  while(phm->produit[i]!=NULL) i++;
  phm->produit[i]=this;
  phm->produit[i+1]=NULL;
  frontiere[1]=new interface(phm,this);
  frontiere[2]=NULL;
}

phase::~phase()
{
    if(frontiere[1]) delete frontiere[1];
	if(frontiere[2]) delete frontiere[2];
	if(Thermo) delete Thermo;
}

void phase::ini(double temperature)
{
  calc_diffusion(temperature);
  calc_germination(temperature);
  calc_croissance(temperature);
  calc_composition();
  calc_fraction();
  int i=1;
  while(produit[i]!=NULL)
    {
      produit[i]->ini(temperature);
      i++;
    }
}

void phase::calc_diffusion(double temperature)
{
	for(int i=1;i<=10;i++) D[i]=D0[i]*exp(-E0[i]/R/temperature);
}

//------------------Classe matrice------------------------------------------------------

matrice::matrice() : phase() {}

void matrice::calc_composition()
{
	double temp1=0.,temp2=0.;
	int i=1;
	while(produit[i]!=NULL)
	{ 
		temp1+=produit[i]->fraction_apparente * produit[i]->xmoy[2] / produit[i]->volmolaire;
		temp2+=produit[i]->fraction_apparente / produit[i]->volmolaire;
		i++;
	}
	xmoy[2]=(xini[2]/volmolaire-temp1)/(1./volmolaire-temp2);
	//	if(xmoy[2]<xint[2]) xmoy[2]=xint[2];
	xmoy[1]=1.-xmoy[2];
}

void matrice::calc_fraction()
{
	int i=1;
	fraction=1.;
	while(produit[i]!=NULL)	{fraction -= produit[i]->fraction;i++;}
}

void matrice::evolution(double dtemps)
{
  int i=1;int j=1;
  while(distri[i]!=NULL) {distri[i]->evolution(dtemps);i++;}
  while(produit[j]!=NULL) {produit[j]->evolution(dtemps);j++;}
}

void matrice::calc_germination(double temperature)
{
}

void matrice::calc_croissance(double temperature)
{
}

//------------------Classe sphere----------------------------------------------------

sphere::sphere() : phase()
{
  eninterface=0;
  K0=0;//Constante d'equilibre
  Q=0;
  Keq=0;
  //Parametres microstructuraux
  //  double rayon;
  densite_apparente=0;
  densite=0;
  //Parametre intermediaire
  effet_Gibbs_Thomson=0;	
}

sphere::sphere(distribution *distri) : phase(distri)
{
    eninterface=0;
  K0=0;//Constante d'equilibre
  Q=0;
  Keq=0;
  //Parametres microstructuraux
  //  double rayon;
  densite_apparente=0;
  densite=0;
  //Parametre intermediaire
  effet_Gibbs_Thomson=0;
}

void sphere::calc_composition()
{
  for(int i=1;i<=10;i++) xmoy[i]=frontiere[1]->xeq2[i];
}

void sphere::calc_fraction()
{
  fraction_apparente=densite_apparente*4./3.*PI*pow(Taille[1],3.);
  
  double temp=0.;
  int i=1;
  phase *tmp;
  
  while(phasemere->produit[i]!=NULL)	
    {
      tmp=phasemere->produit[i];
      temp += tmp->fraction_apparente *
	(phasemere->xmoy[2] - tmp->xmoy[2] * phasemere->volmolaire / tmp->volmolaire);
      
      i++;
    }
  fraction=fraction_apparente * phasemere->xmoy[2] / (phasemere->xini[2] + temp);
  
  densite=fraction / fraction_apparente * densite_apparente;
}

void sphere::evolution(double dtemps)
{
  int i=1;int j=1;
  while(distri[i]!=NULL) {distri[i]->evolution(dtemps);i++;}
  while(produit[j]!=NULL) {produit[j]->evolution(dtemps);j++;}
}

void sphere::calc_germination(double temperature)
{
}

void sphere::calc_croissance(double temperature)
{
  dTaille[1]=phasemere->D[2]/Taille[1]*frontiere[1]->param_cinetique;
}

void sphere::copie(sphere *sph)
{
  phasemere=sph->phasemere;
  eninterface=sph->eninterface;
  xmoy[2]=sph->xmoy[2];
  volmolaire=sph->volmolaire;
}


//--------------------Classe sphere_germ_homogene-----------------------------------------

sphere_germ_homogene::sphere_germ_homogene() : sphere() {}

sphere_germ_homogene::sphere_germ_homogene(distribution *distri) : sphere(distri) {}

void sphere_germ_homogene::evolution(double dtemps)
{
	int i=1;
	while(produit[i]!=NULL) {produit[i]->evolution(dtemps);i++;}
	
	Taille[1]+=dTaille[1]*dtemps;
	densite_apparente+=fluxgermination*dtemps;
	nsites-=fluxgermination*dtemps;
}

void sphere_germ_homogene::calc_germination(double temperature)
{
zeldovich=vatomique/sqrt(kB*temperature)*pow(forcemotrice,2.)/8./PI/pow(eninterfacegerm,1.5);

rayoncritique=2.*eninterfacegerm/forcemotrice;

transport=4.*PI*pow(rayoncritique,2.)/pow(phasemere->parametre_maille,4.)*phasemere->D[2]*phasemere->xmoy[2];

barriere=16./3*PI*pow(eninterfacegerm,3.)/pow(forcemotrice,2.);

fluxgermination=nsites*zeldovich*transport*exp(-barriere/kB/temperature);

fluxgermination=eninterfacegerm*1.e16;
}

//--------------------Classe sphere_germ_homogene_elasticite-----------------------------------------

sphere_germ_homogene_elasticite::sphere_germ_homogene_elasticite() : sphere_germ_homogene() {}

sphere_germ_homogene_elasticite::sphere_germ_homogene_elasticite(distribution *distri) : sphere_germ_homogene(distri) {}


void sphere_germ_homogene_elasticite::calc_germination(double temperature)
{
zeldovich=vatomique/sqrt(kB*temperature)*pow(forcemotrice,2.)/8./PI/pow(eninterfacegerm,1.5);

rayoncritique=2.*eninterfacegerm/forcemotrice;

transport=4.*PI*pow(rayoncritique,2.)/pow(phasemere->parametre_maille,4.)*phasemere->D[2]*phasemere->xmoy[2];

barriere=16./3*PI*pow(eninterfacegerm,3.)/pow(forcemotrice,2.);

fluxgermination=nsites*zeldovich*transport*exp(-barriere/kB/temperature);

fluxgermination=eninterfacegerm*1.e16*2.;
}
