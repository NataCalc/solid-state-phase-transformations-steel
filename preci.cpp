#include "donnees.h"
#include "phase.hpp"
#include "newt.h"
#include "simplex.h"
#include "thermo.h"
#include "qhull.h"
#include "integration.h"
#include "integer.h"
#include <fstream>


float Diff[10];
float taille_gradient;
double dtemps;
float taille_gradient_alpha;
float taille_gradient_betta;
double a,b,Longeur,Ngrain;
extern double Rzin;
extern double dc , dcr;



char fichiersortie[20];
segment seg[20];
int nbseg;

chargement ch; 
chargement ch1;

sortie so;
sortie so1;
sortie so2;
sortie so3;
sortie so4;
sortie so5;
sortie so6;
sortie so7;
sortie so8;
sortie so9;


matrice mat;

site_homogene_constant site_cem ( &mat );
site_homogene_constant site_eps ( &mat );

distribution_moyenne<sphere> cem ( &site_cem );
distribution_moyenne<sphere> eps ( &site_eps );

double temps,dt;
double T,dT;

std::vector <double> CONCENTRATION_MASStoCONCENTRATION_MOLE ( std::vector <double> &C,std::vector <double> &M );
std::vector <double> CONVERT_CONCENTRATION_MOLEtoSITE_FRACTION ( int n,std::vector <double> &m,std::vector <double> &CONCENTRATION_MOLE );
std::vector <double> CONVERT_SITE_FRACTIONtoCONCENTRATION_MOLE ( std::vector <double> &m,std::vector <double> &y1,std::vector <double> &y2 );
void ecriture2D(double temps,std::vector <double> &diff_d, \
				std::vector <double> &diff, \
				std::vector <double> &omega, \
				std::vector <double> &Gnr, \
				std::vector <double> &Gnra, \
				std::vector <double> &Gv, \
				std::vector <double> &Gtr, \
				std::vector <int> &Gcompt, \
				std::vector <double> &Crk, \
				std::vector <double> &Grka, \
				std::vector <double> &del, \
				std::vector <double> &Ccompt_);
void ecriture2D_(double temps, double Rzin,std::vector <double> &diff_d, \
				 std::vector <double> &diff,std::vector <double> &del, \
				 std::vector <double> &omega, \
				 std::vector <double> &Grka, \
				 std::vector <double> &Gnra,std::vector <double> &Gnr, \
				 std::vector <double> &Grk, std::vector <double> &Gv, \
				 std::vector <double> &Gtr, \
				 std::vector <int> &Gcompt_);
void ecriture_growth_dissolution_FeC(std::vector <double> &t,std::vector <double> &Rlong, \
									 std::vector <double> &Xca, \
									 std::vector <double> &Xfea, \
									 std::vector <double> &Xcb, \
									 std::vector <double> &Xfeb, \
									 std::vector <double> &XV, \
									 std::vector <double> &Gtrf, \
									 std::vector <double> &Gfricf);
void ecriture_growth_dissolution_FeCCr(std::vector <double> &t, \
											  std::vector <double> &Rlong, \
											  std::vector <double> &Xca, \
											  std::vector <double> &Xfea, \
											  std::vector <double> &Xcra, \
											  std::vector <double> &Xcb, \
											  std::vector <double> &Xfeb, \
											  std::vector <double> &Xcrb, \
											  std::vector <double> &XV, \
											  std::vector <double> &Gtrf, \
											  std::vector <double> &Gfricf, \
											  std::vector <double> &G); 

void ecriture3D(double temps,std::vector <double> &Gtr, \
				std::vector <double> &Gmuc, \
				std::vector <double> &Gmucr, \
				std::vector <double> &Gfric, \
				std::vector <double> &Gx_c_bcc, \
				std::vector <double> &Gx_cr_bcc, \
				std::vector <double> &Gx_c_fcc, \
				std::vector <double> &Gx_cr_fcc, \
				std::vector <double> &GtrV, \
				std::vector <double> &taille_t, \
				std::vector <double> &Gcomp, \
				std::vector <double> &Gx_c_bcc_eq, \
				std::vector <double> &Gx_cr_bcc_eq, \
				std::vector <double> &Gx_c_fcc_eq, \
				std::vector <double> &Gx_cr_fcc_eq, \
				std::vector <double> &GVeq, \
				std::vector <double> &G_c_nom, \
				std::vector <double> &G_cr_nom);

void ecriture_profileC_FeCCr(double temps,double Rzin,double SURSAT_c,double SURSAT_cr,int p, double *tx, \
							 double *ty, double xi, double xf);

void ecriture_profilePOT_FeCCr(double temps,double Rzin,int m,std::vector < std::vector<double> > &P, \
							   double xi, double xf);


void ecriture_profileC_FeC(double temps,double Rzin,double SURSAT,int m, double *tx, \
						   double *ty);
void ecriture_profilePOT_FeC(double temps,double Rzin,int m, double *tx, std::vector < std::vector<double> > &P);
void cryptograph(std::vector <double> &L,std::vector <double> &B, \
				 std::vector <double> &Fx_c_bcc,std::vector <double> &Fx_cr_bcc, \
				 std::vector <double> &Fx_c_fcc,std::vector <double> &Fx_cr_fcc);

std::vector <double> BILAN_MASSE_PLANE(std::vector <double> &U,double dtempsFeCCr);


extern void calc3dequilibre();
extern bool calc3cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
extern bool calc3cinetique2(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
extern bool calc3cinetique_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);


extern void calc3dequilibreFeNiC();
extern bool calc3cinetique1FeNiC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
extern bool calc3cinetique2FeNiC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
extern bool calc3cinetique1FeNiC_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetique2FeNiC_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetique3FeNiC_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetique4FeNiC_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetique5FeNiC_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetique6FeNiC_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);


extern bool calc3cinetiqueFeNiC_para_sphere(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetique1FeNiC_sphere(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);


extern void calc3dequilibreFeMnC();
extern bool calc3cinetique1FeMnC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
extern bool calc3cinetiqueFeMnC_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetiqueFeMnC_para_sphere(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern bool calc3cinetique1FeMnC_sphere(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);




extern void calc3dequilibreFeMoC();
extern bool calc3cinetique1FeMoC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
extern bool calc3cinetique2FeMoC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);


extern void calc2dequilibre();

extern bool calc2cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin);
extern bool calc2cinetique2(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
extern bool calc2cinetique_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
extern double f_vitess_x1(double V, double x_c_gamma,std::vector <double> CONCENTRATION_NOMINAL);
extern double f_vitess_x2(double V, double x_c_gamma,std::vector <double> CONCENTRATION_NOMINAL);
extern double regulafalsi(double a,double b,double V,std::vector <double> CONCENTRATION_NOMINAL);

extern void calc2dequilibreFeNi();
extern void calc2dequilibreFeC();


extern bool calc4dequilibre();
extern bool calc4cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);

extern void calc5dequilibre();
int calc_spaces3(char *str); 
void selection_sort(std::vector <double> &a , int length, std::vector <int> &index);
void swap(double& first, double& second);
void swapint(int& first, int& second);
int minimum_from(std::vector <double> &a, int position, int length);
void roll(std::vector <double> &a,std::vector <double> &atmp);
int roll_shifted(std::vector <double> &a);


int main ( int argc, char *argv[] )
{
    if ( argc<2 )
    {
        cout<<"Usage: preci dimension (number only)"<<endl;
        return 0;
    }
    int dim=atoi ( argv[1] );
    string bcc;
    string fcc;
    string fail;
	string fail1;
	string fail2;
	string fail3;
	string fail4;
	string fail5;
	string fail6;
	string fail7;
	string fail8;
	string fail9;

    switch ( dim )
    {
    case 2:
        bcc="_bccFeC.txt";
        fcc="_fccFeC.txt";
        fail="sortieFeCeq_Y.txt";
		fail1="sortieFeCcin_Y.txt";	
		fail2="sortieFeCeq_X.txt";
		fail3="sortieFeCcin_X.txt";	
		fail4="sortieFeC_para_U.txt";
		fail5="sortieFeC_para_X.txt";	
        break;
    case 3:
		bcc="_bccFeCrC.txt";
		fcc="_fccFeCrC.txt";
		fail="sortieFeCrCeq_Y.txt";
		fail1="sortieFeCrCcin_Y.txt";
		fail2="sortieFeCrCeq_X.txt";	
		fail3="sortieFeCrCcin_X.txt";
		fail4="sortieFeCCr_para_U.txt";
		fail5="sortieFeCCr_para_X.txt";	
		fail6="sortieFeCrCcin_U.txt";	
		fail7="V.txt";	
		fail8="V1.txt";
        break;
    case 4:
		bcc="_bccFeCrMoC.txt";
		fcc="_fccFeCrMoC.txt";
		fail="sortieFeCrMoCeq_Y.txt";
		fail1="sortieFeCrMoCcin_Y.txt";
		fail2="sortieFeCrMoCeq_X.txt";	
		fail3="sortieFeCrMoCcin_X.txt";
		fail4="sortieFeCrMoC_para_U.txt";
		fail5="sortieFeCrMoC_para_X.txt";	
		fail6="sortieFeCrMoCcin_U.txt";	
		fail7="V.txt";	
		fail8="V1.txt";
        break;
	case 5:
		bcc="_bccFeNiC.txt";
		fcc="_fccFeNiC.txt";
		fail="sortieFeNiCeq_Y.txt";
		fail1="sortieFeNiCcin_Y.txt";
		fail2="sortieFeNiCeq_X.txt";	
		fail3="sortieFeNiCcin_X.txt";
		fail4="sortieFeCNi_para_U.txt";
		fail5="sortieFeCNi_para_X.txt";	
		fail6="sortieFeNiCcin_U.txt";	
		fail7="V.txt";	
		fail8="V1.txt";	
		break;	
	case 6:
		bcc="_bccFeMoC.txt";
		fcc="_fccFeMoC.txt";
		fail="sortieFeMoCeq_Y.txt";
		fail1="sortieFeMoCcin_Y.txt";
		fail2="sortieFeMoCeq_X.txt";	
		fail3="sortieFeMoCcin_X.txt";	
		break;
	case 7:
		bcc="_bccFeMnC.txt";
		fcc="_fccFeMnC.txt";
		fail="sortieFeMnCeq_Y.txt";
		fail1="sortieFeMnCcin_Y.txt";
		fail2="sortieFeMnCeq_X.txt";	
		fail3="sortieFeMnCcin_X.txt";
		fail4="sortieFeCMn_para_U.txt";
		fail5="sortieFeCMn_para_X.txt";	
		fail6="sortieFeMnCcin_U.txt";	
		fail7="V.txt";	
		fail8="V1.txt";	
		break;
	case 8:
		bcc="_bccFeNi.txt";
		fcc="_fccFeNi.txt";
		fail="sortieFeNieq_Y.txt";
		fail1="sortieFeNicin_Y.txt";	
		fail2="sortieFeNieq_X.txt";
		fail3="sortieFeNicin_X.txt";	
		fail4="sortieFeNi_para_U.txt";
		fail5="sortieFeNi_para_X.txt";	
		fail6="sortieFeNicin_U.txt";	
		fail7="V.txt";	
		fail8="V1.txt";
		fail9="T_vs_u_massive.txt";
		break;
	case 9:
		bcc="_bccFeC.txt";
		fcc="_fccFeC.txt";
		fail="sortieFeCeq_Y.txt";
		fail1="sortieFeCcin_Y.txt";	
		fail2="sortieFeCeq_X.txt";
		fail3="sortieFeCcin_X.txt";	
		fail4="sortieFeC_para_U.txt";
		fail5="sortieFeC_para_X.txt";	
		fail6="sortieFeCcin_U.txt";	
		fail7="V.txt";	
		fail8="V1.txt";
		fail9="T_vs_u_massive.txt";
		break;	
		//break;
	default:
        cout<<"Unknown dimension."<<endl<<"Usage: preci dimension (number only)"<<endl<<"Dimension must be 2-7"<<endl;
        return 0;
    }
    mat.Thermo = new Gibbs ( bcc.c_str() );//Fe-Cr-C / Fe-Ni-C
    cem.classe0->Thermo = new Gibbs ( fcc.c_str() );//Fe-Cr-C / Fe-Ni-C
	//return (0);

    ch.ouverture ( "param.txt" );
    ch.lecture();
    ch.fermeture();
	
	ch1.ouverture ( "parametres.txt" );
    ch1.lecture1();
    ch1.fermeture();

	so.ouverture(fail.c_str());
	so1.ouverture(fail1.c_str());
	so2.ouverture(fail2.c_str());
	so3.ouverture(fail3.c_str());
	so4.ouverture(fail4.c_str());
	so5.ouverture(fail5.c_str());
	so6.ouverture(fail6.c_str());
	so7.ouverture(fail7.c_str());
	so8.ouverture(fail8.c_str());
	so9.ouverture(fail9.c_str());

	
    if ( dim==2 )
    {
        calc2dequilibre();
		return ( 0 );
    }
		
    if ( dim==3 )
    {
		calc3dequilibre();
		return ( 0 );
    }

    if ( dim==4 )
    {
		calc4dequilibre();	
        return ( 0 );
	}
	
    if ( dim==5)
    {
		calc3dequilibreFeNiC();
		return ( 0 );
    }
	
	
    if ( dim==6)
    {
        calc3dequilibreFeMoC();
		return ( 0 );
    }
	
    if ( dim==7)
    {
        calc3dequilibreFeMnC();
		return ( 0 );
    }
	if ( dim==8)
    {
        calc2dequilibreFeNi();
		return ( 0 );
    }
	if ( dim==9)
    {
        calc2dequilibreFeC();
		return ( 0 );
    }

}

void chargement::lecture()
{
    sphere *tempcem,*tempeps;
    tempcem=cem.classe0;
    tempeps=eps.classe0;

    char buf[80];
    lignesuivante();
    lignesuivante();
    lignesuivante();

    fscanf ( fe,"%lf\n",&mat.D0[2] );
    lignesuivante();
    fscanf ( fe,"%lf\n",&mat.E0[2] );
    lignesuivante();
    fscanf ( fe,"%lf\n",&mat.xini[2] );
    mat.xmoy[2]=mat.xini[2];
    lignesuivante();
    fscanf ( fe,"%lf\n",&mat.volmolaire );
    lignesuivante();
    fscanf ( fe,"%lf\n",&site_cem.densite_initiale[1] );
    cem.classe0->densite_apparente=site_cem.densite_initiale[1];
    site_cem.densite_disponible[1]=1.;
    lignesuivante();
    fscanf ( fe,"%lf\n",&tempcem->Taille[1] );
    lignesuivante();
    fscanf ( fe,"%lf\n",&tempcem->xmoy[2] );
    lignesuivante();
    fscanf ( fe,"%lf %lf\n",&tempcem->K0,&tempcem->Q );
    lignesuivante();
    fscanf ( fe,"%lf\n",&tempcem->eninterface );
    lignesuivante();
    fscanf ( fe,"%lf\n",&tempcem->volmolaire );
    lignesuivante();
    fscanf ( fe,"%i\n",&nbseg );
    lignesuivante();
    for ( int i=1;i<=1+nbseg;i++ )
    {
        fscanf ( fe,"%lf %lf %lf\n",&seg[i].date,&seg[i].temperature,&seg[i].pas );
    }
	lignesuivante();
    lignesuivante();
    fscanf ( fe,"%s\n",fichiersortie );
	
    //	tempeps->densite_apparente=tempcem->densite_apparente*2.;
    site_eps.densite_initiale[1]=site_cem.densite_initiale[1]*2.;
    eps.classe0->densite_apparente=site_eps.densite_initiale[1];
    site_eps.densite_disponible[1]=1.;
    tempeps->Taille[1]=1.e-8;
    tempeps->xmoy[2]=tempcem->xmoy[2]*2.;
    tempeps->eninterface=tempcem->eninterface*2.;
    tempeps->volmolaire=tempcem->volmolaire*2.;
    tempeps->potentiel[1]=tempcem->potentiel[1];
    tempeps->potentiel[2]=tempcem->potentiel[2];
}

void chargement::lecture1()
{
	char buf[100];
    lignesuivante();
    fscanf ( fe,"%e \n",&Diff[0]);
    lignesuivante();
    fscanf ( fe,"%e \n",&taille_gradient);
	lignesuivante();
	fscanf ( fe,"%e \n",&taille_gradient_alpha);
	lignesuivante();
	fscanf ( fe,"%e \n",&taille_gradient_betta);
	lignesuivante();
	fscanf ( fe,"%lf \n",&dtemps);
}

void sortie::ecriture ( std::vector <double> &Y )
{
    int n = Y.size();
    int i;
    int lignenb = 1;
	if(n==1)
	{
		fprintf ( fs,"%e %e \n",T,Y[0]);
	}
	if(n==5)
	{
		fprintf ( fs,"%e %e %e %e %e %e \n",T,Y[0],Y[1],Y[2],Y[3],Y[4] );
	}	
	if ( n==6 )
    {
        fprintf ( fs,"%e %e %e %e %e %e %e\n",T,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5] );
    }
	if (n==7)
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e\n",T,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6]);
	}		
    if ( n==9 )
    {
		fprintf ( fs,"#T Yx_bcc Yc_bcc Yva_bcc Yx_fcc Yc_fcc Yva_fcc v \n");
        fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e\n",T,Y[1],Y[0],Y[2],Y[3],Y[5],Y[4],Y[6],Y[7],Y[8] );
    }
    if ( n==11 )
    {
        fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e\n",T,Y[1],Y[0],Y[2],Y[3],Y[4],Y[6],Y[5],Y[7],Y[8],Y[9],Y[10] );
    }
    if (n==13)
    {
        fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",T,Y[1],Y[0],Y[2],Y[3],Y[4],Y[5],Y[7],Y[6],Y[8],Y[9],Y[10],Y[11] );
    }

}

void sortie::ecritureX ( std::vector <double> &X )
{
	int n = X.size();
	if ( n==4 )
	{
		fprintf ( fs,"%e %e %e %e %e\n",T,X[0],X[1],X[2],X[3]);
	}
	if ( n==5 )
	{
		fprintf ( fs,"%e %e %e %e %e %e\n",T,X[0],X[1],X[2],X[3],X[4]);
	}
	if ( n ==7 )
	{
		//fprintf ( fs,"#T Xfe_bcc Xx_bcc Xc_bcc Xfe_fcc Xx_fcc Xc_fcc v\n");
		fprintf ( fs,"%e %e %e %e %e %e %e %e\n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6]);
	}	
	if ( n ==8 )
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e\n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7]);
	}
	if ( n==9 )
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e\n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8] );
	}
	if ( n==10 )
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e\n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8],X[9]);
	}
	if ( n==12 )
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8],X[9],X[10],X[11]);
	}
	if ( n==13 )
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8],X[9],X[10],X[11],X[12]);
	}
	if ( n==14 )
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8],X[9],X[10],X[11],X[12],X[13]);
	}
	if ( n==16 )
	{
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",T,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8],X[9],X[10],X[11],X[12],X[13], \
					X[14],X[15]);
	}
}



void sortie::ecritureU(std::vector <double> &U)
{
	int n = U.size();
	if ( n==5 )
	{
		//fprintf ( fs,"#T Ux_fcc Ux_bcc Uc_fcc Uc_bcc v r Delta_x Delta_c time U0x U0c\n");
		fprintf ( fs,"%e %e %e %e %e %e \n",T,U[0],U[1],U[2],U[3],U[4]);
	}
	if ( n==7 )
	{
		//fprintf ( fs,"#T Ux_fcc Ux_bcc Uc_fcc Uc_bcc v r Delta_x Delta_c time U0x U0c\n");
		fprintf ( fs,"%e %e %e %e %e %e %e %e \n",T,U[0],U[1],U[2],U[3],U[4],U[5],U[6]);
	}
	if ( n==11 )
	{
		//fprintf ( fs,"#T Ux_fcc Ux_bcc Uc_fcc Uc_bcc v r Delta_x Delta_c time U0x U0c\n");
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e\n",T,U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10]);
	}
	if(n==15)
	{
		//fprintf ( fs,"#T Ux_fcc Ux_bcc Uc_fcc Uc_bcc muc mux gtr v r Delta_x Delta_c time U0x U0c Dni\n");
		fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",T,U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10], \
								U[11],U[12],U[13],U[14]);
	}	

}


void sortie::ecritureV(std::vector <double> &U,int compt)
{
int size = U.size();
if(size ==11)
{
fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %d\n",T,U[0],U[1],U[2],U[3],U[4],U[5], \
			 U[6],U[7],U[8],U[9],U[10],compt);	
}
if(size ==12)
{
fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %d %e \n",T,U[0],U[1],U[2],U[3],U[4],U[5], \
				 U[6],U[7],U[8],U[9],U[10],compt,U[11]);	
}	
else{	
fprintf ( fs,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",T,U[0],U[1],U[2],U[3],U[4],U[5], \
		 U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14],U[15],compt);}	
}




void ecriture2D(double temps,std::vector <double> &diff_d, \
				std::vector <double> &diff, \
				std::vector <double> &omega, \
				std::vector <double> &Gnr, \
				std::vector <double> &Gnra, \
				std::vector <double> &Gv, \
				std::vector <double> &Gtr, \
				std::vector <int> &Gcompt, \
				std::vector <double> &Crk, \
				std::vector <double> &Grka, \
				std::vector <double> &del, \
				std::vector <double> &Ccompt_)
{
	fstream file;
    char filename[30];
	sprintf ( filename,"t_converge/%e.txt",temps );
	
	int size = Gcompt.size();
	
	file.open ( filename,fstream::out );
    if ( file.fail() )
	cout<<"Error opening file"<<endl;

	file<<"#Temps Rayon moyenne"<<endl;
	file<<"#"<<temps<<"	"<<endl;
	file<<"#U(C,BCC)-(NR)	U(C,FCC)-(NR)	V	Gtr	 mu_c  D_c	taille	Niter U(C,BCC,eq)	U(C,FCC,eq)	V,eq	U0"<<endl;
	for(int i=0; i< size; i++)
	{
		file<<Gnr[i]<<" "<<Gnra[i]<<" "<<Gv[i]<<" "<<Gtr[i]<<" "<<omega[i]<<" "<<diff[i]<<" "<<diff_d[i]<<" "<<Gcompt[i]<<" "<<Crk[i]<<" "<< \
		Grka[i]<<" "<<del[i]<<" "<<Ccompt_[i]<<endl;
	}	
	
	//for(int i=0; i< size; i++)
	//{fprintf(fs,"%e		%e		%e		%i\n",Gnr[i],Grk[i],Gv[i],Gcompt[i]);}
	
	file.close();
	
}


void ecriture2D_(double temps, double Rzin,std::vector <double> &diff_d, \
				std::vector <double> &diff,std::vector <double> &del, \
				std::vector <double> &omega, \
				std::vector <double> &Grka, \
				std::vector <double> &Gnra,std::vector <double> &Gnr, \
				std::vector <double> &Grk, std::vector <double> &Gv, \
				std::vector <double> &Gtr, \
				std::vector <int> &Gcompt_)
{
	fstream file;
    char filename[30];
	sprintf ( filename,"t_converge/%e.txt",temps );
	
	int size = Gcompt_.size();
	
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	
	file<<"Temps Rayon moyenne"<<endl;
	file<<temps<<"	"<<Rzin<<endl;
	file<<"#U(C,BCC)-(NR)	U(C,FCC)-(NR)	U(C,BCC)-(RK)	U(C,FCC)-(RK)	V	Gtr	 OMEGA	delta	D_c		Dcfaux/Dcvrai	Niter"<<endl;
	for(int i=0; i< size; i++)
	{
		file<<Gnra[i]<<" "<<Gnr[i]<<" "<<Grka[i]<<" "<<Grk[i]<<" "<<Gv[i]<<" "<<Gtr[i]<<" "<<omega[i]<<" "<<del[i]<<" "<<diff[i]<<" "<<diff_d[i]<<" "<<Gcompt_[i]<<endl;
	}	
	
	//for(int i=0; i< size; i++)
	//{fprintf(fs,"%e		%e		%e		%i\n",Gnr[i],Grk[i],Gv[i],Gcompt[i]);}
	
	file.close();
	
}

void ecriture3D(double temps,std::vector <double> &Gtr, \
				std::vector <double> &Gmuc, \
				std::vector <double> &Gmucr, \
				std::vector <double> &Gfric, \
				std::vector <double> &Gx_c_bcc, \
				std::vector <double> &Gx_cr_bcc, \
				std::vector <double> &Gx_c_fcc, \
				std::vector <double> &Gx_cr_fcc, \
				std::vector <double> &GtrV, \
				std::vector <double> &taille_t, \
				std::vector <double> &Gcomp, \
				std::vector <double> &Gx_c_bcc_eq, \
				std::vector <double> &Gx_cr_bcc_eq, \
				std::vector <double> &Gx_c_fcc_eq, \
				std::vector <double> &Gx_cr_fcc_eq, \
				std::vector <double> &GVeq, \
				std::vector <double> &G_c_nom, \
				std::vector <double> &G_cr_nom)
{
	fstream file;
    char filename[30];
	sprintf ( filename,"t_converge/%e.txt",temps );
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	
	int size = Gtr.size();

	file<<"#Gtr	Muc	Mucr Gfric U(C,BCC)	U(CR,BCC) U(C,FCC) U(CR,FCC) V  taille ITER U(C,BCC,eq) U(CR,BCC,eq) U(C,FCC,eq) U(CR,FCC,eq) Veq U0_c U0_cr"<<endl;
	for(int i=0; i< size; i++)
	{
		file<<Gtr[i]<<" "<<Gmuc[i]<<" "<<Gmucr[i]<<" "<<Gfric[i]<<" "<<Gx_c_bcc[i]<<" "<<Gx_cr_bcc[i]<<" "<<Gx_c_fcc[i]<<" "<<Gx_cr_fcc[i]<<" "<<GtrV[i]<<" "<<taille_t[i]<<" "<<Gcomp[i]<<" "<<Gx_c_bcc_eq[i]<<" "<<Gx_cr_bcc_eq[i]<<" "<<Gx_c_fcc_eq[i]<<" "<<Gx_cr_fcc_eq[i]<<" "<<GVeq[i]<<" "<<G_c_nom[i]<<" "<<G_cr_nom[i]<<endl;
	}	
	
	//for(int i=0; i< size; i++)
	//{fprintf(fs,"%e		%e		%e		%i\n",Gnr[i],Grk[i],Gv[i],Gcompt[i]);}
	
	file.close();
	
}


void ecriture_growth_dissolution_FeC(std::vector <double> &t,std::vector <double> &Rlong, \
											 std::vector <double> &Xca, \
											 std::vector <double> &Xfea, \
											 std::vector <double> &Xcb, \
											 std::vector <double> &Xfeb, \
											 std::vector <double> &XV, \
											 std::vector <double> &Gtrf, \
											 std::vector <double> &Gfricf) 
{
	fstream file;
    char filename[30];
	

	
	int size = Rlong.size();

	sprintf ( filename,"t_growth_coarsening/%e.txt",t[size-1] );
	
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;

	file<<"#temps	Rayon moyenne	U(C,BCC)	U(Fe,BCC)	U(C,FCC)	U(Fe,FCC)	V"<<endl;
	
	for(int i=0; i<size;i++)
	{
		
file<<t[i]<<"	"<<Rlong[i]<<"	"<<Xca[i]<<"	"<<Xfea[i]<<"	"<<Xcb[i]<<"	"<<Xfeb[i]<<"	"<<XV[i]<<"		"<<Gtrf[i]<<"	"<<Gfricf[i]<<endl;
	}	
	
	file.close();
}


void ecriture_growth_dissolution_FeCCr(std::vector <double> &t, \
											  std::vector <double> &Rlong, \
											  std::vector <double> &Xca, \
											  std::vector <double> &Xfea, \
											  std::vector <double> &Xcra, \
											  std::vector <double> &Xcb, \
											  std::vector <double> &Xfeb, \
											  std::vector <double> &Xcrb, \
											  std::vector <double> &XV, \
											  std::vector <double> &Gtrf, \
											  std::vector <double> &Gfricf, \
												std::vector <double> &G)
{
	fstream file;
    char filename[30];
	
	
	
	int size = Rlong.size();
	
	sprintf ( filename,"t_growth_coarsening/%e.txt",t[size-1] );
	
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	
	file<<"#temps	Rayon moyenne	U(C,BCC)	U(Fe,BCC)   U(Cr,BCC)	U(C,FCC)	U(Fe,FCC)	U(Cr,BCC)  V"<<endl;
	
	for(int i=0; i<size;i++)
	{
		
		file<<t[i]<<"	"<<Rlong[i]<<"	"<<Xca[i]<<"	"<<Xfea[i]<<"	"<<Xcra[i]<<" "<<Xcb[i]<<"	"<<Xfeb[i]<<" "<<Xcrb[i]<<" "<<XV[i]<<"	"<<Gtrf[i]<<"	"<<Gfricf[i]<<" "<<G[i]<<endl;
	}	
	
	file.close();
}



void ecriture_profileC_FeC(double temps,double Rzin,double SURSAT,int m, double *tx, \
						   double *ty)
{
	fstream file;
    char filename[30];
	sprintf ( filename,"t_profileC_FeC/%e.txt",temps);
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	
	file<<"#temps	Rayon moyenne	Delta_Zener"<<endl;
	file<<"#"<<temps<<"	"<<Rzin<<"	"<<SURSAT<<endl;
	file<<"#X		U_c(X)"<<endl;
	for(int i=0;i<=m;i++)
	{
		file<<tx[i]<<"			"<<ty[i]<<endl;
	}	
   file.close();
}

void ecriture_profileC_FeCCr(double temps,double Rzin,double SURSAT_c,double SURSAT_cr,int p, double *tx, \
						   double *ty, double xi, double xf)
{
	fstream file;
    char filename[30];
	sprintf ( filename,"t_profileC_FeCCr/%e.txt",temps);
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	
	
	file<<"#temps	Rayon moyenne	Delta_Zener_C	Delta_Zener_CR"<<endl;
	file<<"#"<<temps<<"	"<<Rzin<<"	"<<SURSAT_c<<" "<<SURSAT_cr<<endl;
	file<<"#X		U_c(X)		U_cr(X)"<<endl;
	int    i;
	double h,x;
	h=(xf-xi)/(p-1); 
	x=xi-h;
	
	for(i=1;i<=p;i++)
	{
		x += h;
		file<<x<<"		"<<tx[i]<<"		"<<ty[i]<<endl;
	}	
	file.close();
}



void ecriture_profilePOT_FeC(double temps,double Rzin,int m, double *tx, std::vector < std::vector<double> > &P)
{
	fstream file;
    char filename[30];
	sprintf ( filename,"t_profilePOT_FeC/%e.txt",temps);
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	file<<"#temps	Rayon moyenne "<<endl;
	file<<"#"<<temps<<"	"<<Rzin<<endl;	
	file<<"#X	C	Fe"<<endl;

	for (int i=0; i<=m; i++)
	{
		file<<tx[i]<<"	"<<P[i][0]<<"	"<<P[i][1]<<endl;
	}
	file.close();

}

void ecriture_profilePOT_FeCCr(double temps,double Rzin,int m,std::vector < std::vector<double> > &P, \
							   double xi, double xf)
{
	fstream file;
    char filename[30];
	sprintf ( filename,"t_profilePOT_FeCCr/%e.txt",temps);
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	file<<"#temps	Rayon moyenne "<<endl;
	file<<"#"<<temps<<"	"<<Rzin<<endl;	
	file<<"#X	C	Fe	CR	CR-Fe"<<endl;
	
	
	int    i;
	double h,x;
	h=(xf-xi)/(m-1); 
	x=xi-h;
	
	for (i=1; i<=m; i++)
	{
		x += h;
		file<<x<<"	"<<P[i][0]<<"	"<<P[i][1]<<"	"<<P[i][2]<<"	"<<P[i][3]<<endl;
	}
	file.close();
	
}

void cryptograph(std::vector <double> &L,std::vector <double> &B, \
				 std::vector <double> &Fx_c_bcc,std::vector <double> &Fx_cr_bcc, \
				 std::vector <double> &Fx_c_fcc,std::vector <double> &Fx_cr_fcc)
{
	
	fstream file;
    char filename[30];
	sprintf ( filename,"t_converge/cryptograph.txt",temps );
	file.open ( filename,fstream::out );
    if ( file.fail() )
		cout<<"Error opening file"<<endl;
	
	int n = L.size();
	
	file<<"#L	B	U(C,BCC)	U(CR,BCC)	U(C,FCC)	U(CR,FCC)"<<endl;
	for(int i=0; i< n; i++)
	{
		file<<L[i]<<" "<<B[i]<<" "<<Fx_c_bcc[i]<<" "<<Fx_cr_bcc[i]<<" "<<Fx_c_fcc[i]<<" "<<Fx_cr_fcc[i]<<endl;
	}	
	

	
	file.close();
	
	
	
}

//------------------------------------------------------------------------------
//              Convert concentration mass to concentration mole
//______________________________________________________________________________
std::vector <double> CONCENTRATION_MASStoCONCENTRATION_MOLE ( std::vector <double> &C,std::vector <double> &M )
{
    int n = C.size();
    int i;
    std::vector <double> CONCENTRATION_MOLE ( n );
	
    double dummy = 0.;
	
    for ( i = 0; i < n; i++ )
    {
        dummy += C[i]/M[i];
    }

	
    for ( i = 0; i < n; i++ )
    {
        CONCENTRATION_MOLE[i] = (C[i]/M[i])/dummy;
    }

	
    return CONCENTRATION_MOLE;
}

//------------------------------------------------------------------------------
//              Convert concentration mole to concentration mass
//______________________________________________________________________________

std::vector <double> CONCENTRATION_MOLEtoCONCENTRATION_MASS ( std::vector <double> &X, std::vector <double> &M )
{
    int n = X.size();
    int i;
    std::vector <double> CONCENTRATION_MASS ( n );
    double dummy = 0.;
	
    for ( i = 0; i < n; i++ )
    {
        dummy += X[i]*M[i];
    }
	
    for ( i = 0; i < n; i++ )
    {
        CONCENTRATION_MASS[i] = X[i]*M[i]/dummy;
    }
    return CONCENTRATION_MASS;
}

//------------------------------------------------------------------------------
//              Convert concentration mole to site fraction
//______________________________________________________________________________
std::vector <double> CONVERT_CONCENTRATION_MOLEtoSITE_FRACTION ( int n,std::vector <double> &m,std::vector <double> &CONCENTRATION_MOLE )
{
    int i;
    std::vector <double> SITE_FRACTION ( n );
	
    SITE_FRACTION[0] = m[0]/m[1]*CONCENTRATION_MOLE[0]/ ( 1.-CONCENTRATION_MOLE[0] );
	
    for ( i = 1; i < n; i++ )
    {
        SITE_FRACTION[i] = CONCENTRATION_MOLE[i]* ( 1.+m[1]/m[0]*SITE_FRACTION[0] );
    }
    return SITE_FRACTION;
}

//---------------------------------------------------------------------------------
//               Convert site fraction to concentration mole Fe-C
//---------------------------------------------------------------------------------
std::vector <double> CONVERT_SITE_FRACTIONtoCONCENTRATION_MOLE ( std::vector <double> &m,std::vector <double> &y1,std::vector <double> &y2 )
{
    int size1 = y1.size();
    int size2 = y2.size();
    int n = size1+size2-1;
    std::vector <double> CONCENTRATION_MOLE ( n );
	
    CONCENTRATION_MOLE[0] = y2[0]*m[1]/ ( m[0]+m[1]-m[1]*y2[1] ); //C
    CONCENTRATION_MOLE[1] = y1[0]*m[0]/ ( m[0]+m[1]-m[1]*y2[1] ); //FE
	
    return CONCENTRATION_MOLE;
}





std::vector <double> BILAN_MASSE_PLANE(std::vector <double> &U,double dtempsFeCCr)
{
	double surface,Volume_fcc,Volume_bcc;
	int n=U.size();
	std::vector <double> delta(2), omega(2),X(2);
	double Uc0_fcc, Ux0_fcc;
	
	surface=a*b;
	Volume_fcc=a*b*(Longeur-Rzin);
	//if(n==7)
	//{
		//U[0]=Ucbcc U[1]=Uxbcc U[2]=Ucfcc U[3]=Uxfcc U[4]=Uc0fcc U[5]=Ux0fcc U[6]=vitesse
		omega[0] = (U[2]-U[4])/(U[2]-U[0]);
		omega[1] = (U[3]-U[5])/(U[3]-U[1]);
		for (int i =0; i < 2; i++)
			delta[i] = 2.*Rzin*(1. - omega[i])/omega[i];	
		
		Uc0_fcc=( surface*dc/(delta[0]*Volume_fcc) ) * dtempsFeCCr \
		*(U[2]-U[4]) \
		-( U[6]*surface/Volume_fcc ) * dtempsFeCCr \
		*(U[2]-U[4])+U[4];
		Ux0_fcc=( surface*dcr/(delta[1]*Volume_fcc) ) * dtempsFeCCr \
		*(U[3]-U[5]) \
		-( U[6]*surface/Volume_fcc ) * dtempsFeCCr \
		*(U[3]-U[5])+U[5];
		X[0] = Uc0_fcc/(1.+Uc0_fcc); //c0_fcc
		X[1] = Ux0_fcc/(1.+Uc0_fcc); //x0_fcc
	//}
	return X;
}


int calc_spaces3(char *str) 
{
char *tmp;
int ans = 0;
bool sp = false;
while (*str == ' ')
	++str;
for (tmp = str; *tmp; ++tmp)
{
if (isdigit(*tmp) && sp)
{
++ans;
sp = false;
}
if (*tmp == ' ')
	sp = true;
}
return ans;
}	

/* DEFINITION OF FUNCTION "selection_sort" */
void selection_sort(std::vector <double> &a , int length, std::vector <int> &index)
{
for (int count = 0 ; count < length - 1 ; count++)
{
swapint(index[count],index[minimum_from(a,count,length)]);
swap(a[count],a[minimum_from(a,count,length)]);
}
}


void swap(double& first, double& second)
{
double temp = first;
first = second;
second = temp;
}

void swapint(int& first, int& second)
{
int temp = first;
first = second;
second = temp;
}


int minimum_from(std::vector <double> &a, int position, int length)
{
int min_index = position;
for (int count = position + 1 ; count < length ; count ++)
if (a[count] < a[min_index])
	min_index = count;
return min_index;
}


void roll(std::vector <double> &a,std::vector <double> &atmp)
{
int length= a.size();
int count_change= roll_shifted(a);
//printf("a[%d]=%e a[0]=%e \n",count_change,a[count_change],a[0])
atmp.resize(length);
atmp[0]=a[count_change];
for(int count=1; count<length;count++)
	atmp[count]=a[count-1];
}

int roll_shifted(std::vector <double> &a)
{
int length=a.size();
int min_index = length-1;
return min_index;
}

