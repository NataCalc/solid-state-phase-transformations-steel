#include "donnees.h"
#include "phase.hpp"
#include "newt.h"
#include "thermo.h"
#include "qhull.h"
#include "integer.h"
#include <fstream>
#include <time.h>
#include <iomanip>
#include <iostream>


extern double temps,dt,T,dT;
extern segment seg[20];
extern int nbseg;
extern double x_i,x_j;
extern float CIN;
extern float taille_gradient;
extern double M, Vm, coeff_n, coeff_k, Lcoef, Bcoef;
extern double M0;
extern double VIT_TEMPORAIRE;
extern double gtr1,gtr, muc, mucr;
extern sortie so;
extern sortie so1;
extern sortie so2;
extern sortie so3;
extern sortie so4;
extern sortie so5;
extern sortie so6;
extern sortie so7;
extern sortie so8;
extern double x_c0;
extern double x_cr0;

extern matrice mat;
extern distribution_moyenne<sphere> cem;
extern distribution_moyenne<sphere> eps;
extern double dc,dcr,Rzin,correc_D;
extern double y_fraction_c_bcc;
extern double y_fraction_x_bcc; 
extern double y_fraction_x_fcc;
extern double y_fraction_c_fcc;
extern double dtempsFeCCr,tempsFeCCr;
extern int a1,a2,def;//lattice 


void calc3dequilibreFeMnC();
bool cinetiqueLE_FeCMn(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc3cinetiqueFeMnC_growth(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
bool calc3cinetiqueFeMnC_complet(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
bool cinetiqueLE_ch_FeCMn(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc3cinetiqueFeMnC_growth_ch(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
bool cinetiqueLE_qh_FeCMn(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &Ycinetique_para_ch, bool &flage_ch);
bool calc3cinetiqueFeMnC_BalanceC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);
bool calc3cinetiqueFeMnC_para_or_LE(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage_ch);
bool calc3cinetiqueFeMnC_para_or_LE_T(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage_ch);

extern void print_vect(std::vector <double> &Y);
extern void C_YtoX( std::vector <double> &Y,std::vector <double> &X);
extern void C_YatoXa( std::vector <double> &Y,std::vector <double> &X);
extern void C_XtoU( std::vector <double> &X,std::vector <double> &U);
extern void C_XatoUa( std::vector <double> &X,std::vector <double> &U);
extern void C_X0toU0(std::vector <double> &X0, std::vector <double> &U0);
extern double Sursaturation_c(std::vector <double> &U,std::vector <double> &U0);
extern double Sursaturation_ni(std::vector <double> &U,std::vector <double> &U0);
extern int calc_spaces3(char *str); 
extern void selection_sort(std::vector <double> &a , int length, std::vector <int> &index);
extern void swap(double& first, double& second);
extern void swapint(int& first, int& second);
extern void roll(std::vector <double> &a,std::vector <double> &atmp);
extern int roll_shifted(std::vector <double> &a);
extern std::vector <double> BILAN_MASSE_PLANE_U(std::vector <double> &U,double dtempsFeCCr);



extern fstream file_v;
extern clock_t t_debut_main , t_debut_qhull , t_fin_qhull , t_debut_NR , t_fin_NR , t_fin_main ;


void calc3dequilibreFeMnC()
{
//The difinition dimension of alloy
int dim=3;
def=5;
//Definition Newton-Raphson and convex hull	
newt nt;
qhull qh(dim);
//Temporaire parameters	
t_debut_main = clock();
double counter=time(0);
double prevT=0;
bool conv = false;	
std::vector<double> prevVector(9);		
std::vector<double> prevVectorcin(9);
for(int i=0;i<9;i++){prevVectorcin[i]=0;}		
for(int i=0;i<9;i++){prevVector[i]=0;}	
//Definition of nominal composition of alloy	
std::vector <double> CONCENTRATION_NOMINAL ( dim );
std::vector <double> CONCENTRATION_MOLE ( dim );
CONCENTRATION_NOMINAL[0] = 0.107e-2;// carbon -C (fraction mole)
CONCENTRATION_NOMINAL[2] = 0.1726e-2;// nickel -Mn (fraction mole)
CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0]-CONCENTRATION_NOMINAL[2];// iron -FE
x_i=CONCENTRATION_NOMINAL[0]; //convex hull
x_j=CONCENTRATION_NOMINAL[2];	

printf ( "Begin calcul FE-Mn-C system : \n" );
for ( int iterseg=1;iterseg<=nbseg;iterseg++ )
{
//Definition of temperature range	
temps=seg[iterseg].date;
dt= ( seg[iterseg+1].date-temps ) /seg[iterseg+1].pas;
		
T=seg[iterseg].temperature;
dT= dt * ( seg[iterseg+1].temperature-seg[iterseg].temperature )
/ ( seg[iterseg+1].date-seg[iterseg].date );
		
while ( temps<seg[iterseg+1].date )
{
bool prevExts=false;
bool good;
bool flage = false;

qh=qhull(dim);
temps+=dt;
T+=dT;
//if (( T>=1109. ) && ( T<=1771.)) continue;
if (prevT==0)
{prevT=T;}
else if ((T-prevT)<=dT)
{
prevExts=true;
prevT=T;
}
	
//Input parameters:
dc =  2.5e-5 * exp( -1.4421e5 /R/ T )*pow(1e6,2.);//true
//dc =  6.8*1.e-13*pow(1e6,2.);	//odqvist
dcr = 1.6e-5 * exp( -2.6125e5 /R/ T ) *pow(1e6,2.); //ni
Rzin = 1e-2;//micro meters
tempsFeCCr=0.0;///delete
taille_gradient=1.e-3;	
a1 = 1;//lattice
a2 = 3;//lattice 
dtempsFeCCr=0.1;
M=0.035*exp(-17700/T)*pow(1e6,4);
Vm=1e13;	
coeff_n=2.2e-4 * exp(-1.225e5/R/T)*pow(1e6,2.);//1e13; //(Dc)
coeff_k=1e-4;//1e6; //(Dni)
Lcoef=0.0;Bcoef=0.0;
nt=newt();
		
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func=mat.produit[1]->frontiere[1];
nt.nrfuncv= &interface::ST;
		
nt.jacobian=mat.produit[1]->frontiere[1];
nt.jacobfunc = &interface::JT;
		
qh.func1    =mat.produit[1]->frontiere[1];
qh.nrfuncv1 =&interface::phase1_ternaire;
		
qh.func2    =mat.produit[1]->frontiere[1];
qh.nrfuncv2 =&interface::phase2_ternaire;
		
qh.Jacob1=mat.produit[1]->frontiere[1];
qh.nrJacob1=&interface::JacobianQhull3Dalfa;
		
qh.Jacob2=mat.produit[1]->frontiere[1];
qh.nrJacob2=&interface::JacobianQhull3Dbetta;
i++;
}
	
if ( nt.func == NULL ) printf ( "newton.func = null\n" );
if ( qh.func1== NULL ) printf ( "qhull.func = null\n" );
				
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 9;
std::vector <double> Y ( nt.nn );
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn ); //Jacobian
	
printf ( "T = %lf \n", T );
mat.Thermo->read ( T );
cem.classe0->Thermo->read ( T );
			
/*if ((prevExts) && (good))Y=prevVector;
else*/
{
t_debut_qhull=clock();
qh.discretization (T,25);
qh.build_hull();
vector<point> ext=qh.extremums(/*2.e-5*/);
t_fin_qhull = clock();

if(ext.empty()) continue;
if ( ( !ext.empty() ) )
{
printf ( "extremums:" );
for ( int i=0;i<ext.size();i++ )
	cout<<i+1<<":"<<ext[i]<<endl;
}
else
{
printf ( "No extremums T = %lf\n",T );
continue;
}
cout<<"Current extremums:"<<ext[0]<<";"<<ext[1]<<endl;
Y[1] = ( ext[0].func()?ext[1][1]:ext[0][1] );  //Mn
if ( Y[1]==0. )
	Y[1]=1.e-4;
Y[0] = 1.-Y[1];       //FE
Y[2] = ( ext[0].func()?ext[1][0]:ext[0][0] );  //C
if ( Y[2]==0. )
	Y[2]=1.e-4;
Y[3] = 1.-Y[2];       //VA
Y[5] =( ext[0].func()?ext[0][1]:ext[1][1] );   //Mn
if ( Y[5]==0. )
	Y[5]=1.e-4;
Y[4] = 1.0-Y[5];      //FE
Y[6] = ( ext[0].func()?ext[0][0]:ext[1][0] );  //C
if ( Y[6]==0.0 )
	Y[6]=1.e-4;
Y[7] = 1.-Y[6];       //VA
Y[8] = 0.45;			
}
printf("Input vector from Convex Hull");	
print_vect(Y);
t_debut_NR = clock();			
int check=0;

//char filename[30];
//sprintf ( filename,"newt55.txt");
//file_v.open ( filename,fstream::out );
//if ( file_v.fail() )
//cout<<"Error opening file"<<endl;
	
	
good=nt.calc ( Y, &check,CONCENTRATION_NOMINAL,3,conv);
	
//file_v.close();

	
t_fin_NR = clock() ;
t_fin_main = clock();
double ff= CLOCKS_PER_SEC;
			
for ( int l = 0; l < nt.nn; l++ )
{
if ( Y[l]<0.0 )
	break;
}
printf("Output of NR");	
print_vect(Y);

			
//convergance FRACTION_SITE dans FRACTION_MOLAIRE
std::vector <double> X(nt.nn-2);
C_YtoX(Y,X);	

if (good)
{
so.ecriture ( Y );
so2.ecritureX ( X );
prevVector[0]=Y[0];
prevVector[1]=Y[1];
prevVector[2]=Y[2];
prevVector[3]=Y[3];
prevVector[4]=Y[4];
prevVector[5]=Y[5];
prevVector[6]=Y[6];
prevVector[7]=Y[7];
prevVector[8]=Y[8];
prevT=T;
}
printf("\n\nFIN\n\n");
printf("\n\nTemps CPU : temps (s)	\
\n-QuickHull       : %e  %e s		\
\n-NR              : %e  %e s		\
\n-Total           : %e  %e s		\
\n",							\
(double) t_fin_qhull - (double) t_debut_qhull , ((double) t_fin_qhull - (double) t_debut_qhull) / ff ,	\
(double) t_fin_NR - (double) t_debut_NR , ((double) t_fin_NR - (double) t_debut_NR) / ff ,  \
(double) t_fin_main - (double) t_debut_main , ((double) t_fin_main - (double) t_debut_main) / ff ) ;
			
	

/***************************************Cinetique********************************/

if(good)
{
printf("Calcul cinetique plane \n");
if((T>=400.) && (T<=1500.))
cinetiqueLE_FeCMn(CONCENTRATION_NOMINAL,Y,prevVectorcin,flage);
}

}
}
so1.fermeture();
so2.fermeture();
so3.fermeture();
so4.fermeture();
so5.fermeture();
so6.fermeture();
	
so.fermeture();
}


bool cinetiqueLE_FeCMn(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage)
{
newt nt;
nt=newt();
integer itt2;
itt2=integer();
int i=1;
while (mat.produit[i]!=NULL)
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT;
		
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT;
	
itt2.funcin	 =mat.produit[1]->frontiere[1];
itt2.nrfuncin =&interface::funcFeCrC;
itt2.funcpot  =mat.produit[1]->frontiere[1];
itt2.nrfuncpot=&interface::funconeFeCrC;
itt2.funcderivFeCrC = mat.produit[1]->frontiere[1];
itt2.nrfuncderivFeCrC = &interface::funcderivmuFeCrC;	
	
i++;
}
if (nt.func == NULL) printf ("newton.func = null\n");
if (nt.jacobian==NULL) printf ("jacobian.func=null\n");	

int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn=9;
std::vector <double> Ycinetique(nt.nn);
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 
double before=0.0; double after=0.0; float v_ini;
int j=0; int diff_step; int compteur=0;
bool conv = false;

FILE *fpp=fopen("coeff_cin.txt","r") ;
fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
fscanf(fpp , "%i  \n" , &diff_step) ;	
fclose(fpp);	
double raison=pow(dc/dcr,1./((double)diff_step)) ;	
	
//Principle
printf("Kinetic Fe-C-Mn with local equilibrium at interface \n\n\n");	
//while(tempsFeCCr<=400)
{
if(j==0)
{	
Ycinetique[0] = Y[0];//FE
Ycinetique[1] = 1.-Y[0];//CR
Ycinetique[2] = Y[2];//C
Ycinetique[3] = 1.-Y[2];//VA
Ycinetique[4] = Y[4];//FE
Ycinetique[5] =	1.0-Y[4];//CR
Ycinetique[6] = Y[6];//C
Ycinetique[7] = 1.-Y[6];//VA
Ycinetique[8] = 0.0;
}
else
{
dtempsFeCCr=0.1;
Ycinetique[0]=prevVectorcin[0];
Ycinetique[1]=prevVectorcin[1];
Ycinetique[2]=prevVectorcin[2];
Ycinetique[3]=prevVectorcin[3];
Ycinetique[4]=prevVectorcin[4];
Ycinetique[5]=prevVectorcin[5];
Ycinetique[6]=prevVectorcin[6];
Ycinetique[7]=prevVectorcin[7];
Ycinetique[8]=prevVectorcin[8];
before=Rzin;
Rzin=Rzin+prevVectorcin[8]*dtempsFeCCr;
after=Rzin; //printf("after=%e \n",after);
while((after-before)/before>1e-2)
{
dtempsFeCCr=dtempsFeCCr/10.;
Rzin=before+prevVectorcin[8]*dtempsFeCCr;
after=Rzin;
}
tempsFeCCr=tempsFeCCr+dtempsFeCCr;
}
correc_D=0;
if(flage)
{	
correc_D=dc/dcr;
printf("Input Y in NR:\n");	
print_vect(Ycinetique);	
nt.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
printf("Output Y from NR:\n");	
print_vect(Ycinetique);	
printf("tempsFeCCr=%e \n",tempsFeCCr);
}		
else
{	
compteur=0;
while(compteur<= diff_step)
{			
correc_D = pow(raison,compteur);
printf("Input Y in NR:\n");	
print_vect(Ycinetique);
nt.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
printf("Output Y from NR:\n");	
print_vect(Ycinetique);
printf("diff_step= %i \n",diff_step);
printf("compt=%i \n",compteur);
printf("correc_D=%e \n", correc_D);	
printf("dc=%e dni=%e dni=%e \n",dc,dcr,correc_D*dcr);
compteur ++;
}}

if(conv)
{
std::vector <double> Xcinetique(7);
std::vector <double> Ucinetique(4),U(11),Unom(2);

C_YtoX(Ycinetique,Xcinetique);	
C_XtoU(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("\n");	
printf("Solution in mole fraction\n");
printf("X(C,BCC) = %e  X(C,FCC)= %e \n",Xcinetique[2],Xcinetique[5]);
printf("X(Mn,BCC) = %e  X(Mn,FCC)= %e \n",Xcinetique[1],Xcinetique[4]);
printf("XMn = %e  XC= %e \n",CONCENTRATION_NOMINAL[2],CONCENTRATION_NOMINAL[0]);
printf("v = %e \n",Xcinetique[6]);
printf("\n\n");
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",Xcinetique[6]);

	
U[0]=Ucinetique[0];U[1]=Ucinetique[1];
U[2]=Ucinetique[2];U[3]=Ucinetique[3];
U[4]=Ycinetique[8];U[5]=Rzin;U[6]=Sursaturation_ni(Ucinetique,Unom);
U[7]=Sursaturation_c(Ucinetique,Unom);U[8]=tempsFeCCr;
U[9]=Unom[1];U[10]=Unom[0];	
so1.ecriture(Ycinetique);
so3.ecritureX(Xcinetique);
so6.ecritureU(U);
flage = true;
}
else 
{
printf("No convergence \n temps= %e dt = %e \n\n",tempsFeCCr,dtempsFeCCr);
flage=false;
}
j++;	
}	
if((T>=600.) && (T<=1500.)){
printf("\n\n");	
printf("Kinetic of ferrite formation with dissipation of Gibbs energy at the interface\n");
printf("Simplified model in Velocity - > Complet model\n");
bool flage_para = false;
std::vector<double> prevVectorcin_para(5);
prevVectorcin_para[0] = Y[2];//C in bcc
prevVectorcin_para[1] =	Y[1];//X in bcc
prevVectorcin_para[2] = Y[6];//C in fcc
prevVectorcin_para[3] = Y[5];//X in fcc
prevVectorcin_para[4] = 1.e-9;//v

	
	
/*//=========================Profile===============//
taille_gradient=1e-7;
int    fi,i,ndata,p;
double xi,xf;
double yi[10];
double v1[50],v2[50],v3[50],v4[50];
double hdiff,sdiff,temporaire;
double tdiff[12],tdiff1[12];
p=2;
xi=0.0;
xf=1.0;
ndata=11;
fi=30;
std::vector < std::vector<double> > P(ndata+1);
for(i = 1; i < P.size(); i++)
{
P[i].resize(4);
}
std::vector < std::vector<double> > P1(ndata+1);
for(i = 1; i < P1.size(); i++)
{
P1[i].resize(4);
}
std::vector <double> Ycinetique_para_profile ( 9 );
//======================================================//
	
printf("\n      DIFFERENTIAL EQUATIONS WITH P VARIABLES OF ORDER 1\n");
printf("            of type yi' = f(y1,y2,...,yn), i=1..n\n\n");
	
VIT_TEMPORAIRE=Ycinetique[8];
double xc=3. * Ycinetique[2] / (1. + 3. * Ycinetique[2]);
double xcr=1. * Ycinetique[1] / (1. + 3. * Ycinetique[2]);
yi[0]=xc/(1.-xc);
yi[1]=xcr/(1.-xc);
x_c0	= yi[0];
x_cr0	= yi[1];		
printf("X(C,BCC) = %e  X(CR,BCC) = %e V = %e correct_D=%e taille= %e\n",xc,xcr,VIT_TEMPORAIRE,correc_D, \
		   taille_gradient);
printf("U(C,BCC) = %e  U(CR,BCC) = %e \n",yi[0],yi[1]);
itt2.Eqdifp(v1,v2,v3,v4,xi,xf,yi,ndata,p,fi,CONCENTRATION_NOMINAL,Ycinetique_para_profile,P,P1);
printf ( "\n" );printf ("End Runga-Kutta\n");
printf("\n*******Convergence growth diffusion\n");
*/	
	
//calc3cinetiqueFeMnC_growth(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin_para,flage_para);
//calc3cinetiqueFeMnC_BalanceC(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin_para,flage_para);
//calc3cinetiqueFeMnC_para_or_LE(CONCENTRATION_NOMINAL,Y,prevVectorcin_para,flage_para);
}

/*printf("              Calcul LE or PE\n");
calc3cinetiqueFeMnC_para_or_LE(CONCENTRATION_NOMINAL,Y,prevVectorcin_para,flage_para);}*/
	
/*printf("              Calcul paraequilibre avec la precipitation du plane dans la zone problematique\n");
printf("              Calcul la croissance avec le bilan du masse\n");
printf("              Balance in C varied Ucbcc\n");	
bool flage_para = false;
std::vector<double> prevVectorcin_para(4);	
calc3cinetiqueFeMnC_BalanceC(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin_para,flage_para);
*/


	
}


bool calc3cinetiqueFeMnC_growth(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage)
{
newt nt;
nt=newt();
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT_prof5;
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT_prof5;
i++;
}


bool conv_para = false;
int j=0;
int compteur, diff_step;
double VIT_TEMPORAIRE0,raison,Delta_c,Delta_ni; 	
bool flage_growth_plus, flage_growth_minus;	
correc_D=1.;

	
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 4;
std::vector <double> Ycinetique_veloc(nt.nn);
std::vector <double> Xcinetique(nt.nn);	
std::vector <double> Ucinetique(nt.nn);		
std::vector <double> Unom(2);
int compt = 0;	
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 

Ycinetique_veloc[0]=Y[2]; //C in ferrite
Ycinetique_veloc[1]=Y[1];//X in ferrite
Ycinetique_veloc[2]=Y[6];//C in austenite
Ycinetique_veloc[3]=Y[5];//X in austenite
VIT_TEMPORAIRE=Y[8];//Velocity fixe


diff_step=1000;
raison = pow( pow (10. , 8) , 1. / (double) diff_step ) ;	
compteur = 0 ;	
VIT_TEMPORAIRE0=VIT_TEMPORAIRE;
flage_growth_plus=false;
flage_growth_minus=false;	

	
while (compteur <= diff_step)
{
VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;
	
nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,1,conv_para);
if(!conv_para) {flage_growth_minus=false;break ;}
else flage_growth_minus=true;


	
C_YatoXa(Ycinetique_veloc,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);

	
if(conv_para){	
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr;Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
	
so7.ecritureV ( Ua,compt );
	
}

	
compteur ++;
compt++;	
}

printf("\nBoucle dans l'autre sens : augmentation de la vitesse\n") ;
	
Ycinetique_veloc[0]=Y[2]; //C in ferrite
Ycinetique_veloc[1]=Y[1];//X in ferrite
Ycinetique_veloc[2]=Y[6];//C in austenite
Ycinetique_veloc[3]=Y[5];//X in austenite
VIT_TEMPORAIRE=Y[8];//Velocity fixe


compteur = 0 ;

while (compteur >= -diff_step)
{		
	
VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;

nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,1,conv_para);
if(!conv_para){ flage_growth_plus=false;break ;}
if(VIT_TEMPORAIRE>200.){break;}
else flage_growth_plus=true;

C_YatoXa(Ycinetique_veloc,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);
//printf(" %e  %e %e %e \n",Ucinetique[0],Ucinetique[1],Ucinetique[2],Ucinetique[3]);	
if(conv_para){
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr;Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
so7.ecritureV ( Ua,compt );	
}
compteur --;
compt++;	
}
	
so7.fermeture();
//--------------------------------
//--------------------------------
//--------------------------------
//Open the file with all solution
FILE *fee;
char mystring [100];
//sprintf(mystring,"LENP/%eV.txt",T);
sprintf(mystring,"V.txt");
cout<<mystring<<endl;	
fee=fopen(mystring,"rt");
printf("open \n");	
if(!fee){
printf("Error open file \n");
return 0;}
else
{ 
long lSize;
fseek(fee , 0 , SEEK_END);
lSize = ftell(fee);
rewind(fee);
if(lSize==0)
{fclose (fee);
printf("File is empty\n");
return 0;
}	
}
if(fee) printf("File is opened\n");

char str[256];
int nb_lignes = 0;
int nm = 0;
bool wasNumbers = false;

int nx=1000000;
std::vector <float> xTs(nx),xUxfcc(nx),xUxbcc(nx), \
	xUcbcc(nx),xUcfcc(nx),xmusc(nx),xmusx(nx),xgstr(nx), \
	xfvs(nx),xts(nx), \
	xU0sx(nx),xU0sc(nx),xcompts(nx), \
	xrs(nx),xDeltasx(nx), \
	xvs(nx),xDeltasc(nx),xDsni(nx);
float Ts,Uxfcc,Uxbcc, \
	Ucbcc,Ucfcc,musc,musx,gstr, \
	vs,fvs,ts, \
	U0sx,U0sc,rs,Deltasx,Deltasc,Dsni;
int compts;

	
//While file with results is opened we search the solution	
while (!feof(fee))
{
str[0] = 0;
char *tmp = fgets(str, sizeof(str), fee);
if (!tmp || str[0] == 0) //If have found nothing,
break; //we leave
		
if (str[0] == '/') //Comment
{
if (wasNumbers)//If we considered numbers from a matrix
{
wasNumbers = false;
++nm; //Means we pass to the second matrix
}continue;
}
//printf("File name: %s\n",str);
int sp = calc_spaces3(str); //We consider quantity of blanks

if(nm==0)
{
//"#T  U(Mn,fcc) U(Mn,bcc) U(C,fcc) U(C,bcc) muc mucx gtr v Fv r Delta_x time Delta_c U0x U0c Dx/Dx"<<endl;	
sscanf(str, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",&Ts,&Uxfcc,&Uxbcc, \
				&Ucfcc,&Ucbcc,&musc,&musx,&gstr, \
				&vs,&fvs,&rs,&Deltasx,&ts,&Deltasc, \
				&U0sx,&U0sc,&Dsni,&compts);
//printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",Ts,Uxfcc,Uxbcc, \
	Ucfcc,Ucbcc,musc,musx,gstr, \
	vs,fvs,rs,Deltasx,ts,Deltasc, \
	   U0sx,U0sc,Dsni,compts);
	
xTs[compts]		= Ts;
xUxfcc[compts]	= Uxfcc;
xUxbcc[compts]	= Uxbcc;
xUcbcc[compts]	= Ucbcc;
xUcfcc[compts]	= Ucfcc;
xmusc[compts]	= musc;
xmusx[compts]	= musx;
xgstr[compts]	= gstr;
xvs[compts]		= vs;
xfvs[compts]	= fvs;
xrs[compts]		= rs;
xDeltasx[compts]= Deltasx;
xts[compts]		= ts;
xDeltasc[compts]= Deltasc;
xU0sx[compts]	= U0sx;
xU0sc[compts]	= U0sc;
xDsni[compts]	= Dsni;
xcompts[compts] = compts;	
wasNumbers = true;
}
}
int fi = compts+1;
xTs.resize(fi);xUxfcc.resize(fi);xUxbcc.resize(fi);
xUcbcc.resize(fi);xUcfcc.resize(fi);xmusc.resize(fi);
xmusx.resize(fi);xgstr.resize(fi);
xvs.resize(fi);xfvs.resize(fi);xrs.resize(fi);
xDeltasx.resize(fi);xts.resize(fi);xDeltasc.resize(fi);
xU0sx.resize(fi);xU0sc.resize(fi);xDsni.resize(fi);
xcompts.resize(fi);
fclose(fee);

std::vector <int> index(fi);
std::vector <double> xVx_sorted(fi),xUcbcc_sorted(fi),xUxbcc_sorted(fi), \
xUcfcc_sorted(fi), xUxfcc_sorted(fi), xUxc_sorted(fi), \
xUxx_sorted(fi), xgxtr_sorted(fi), \
xfm_sorted(fi), xmuxc_sorted(fi), \
xmuxx_sorted(fi);
	
for(int count=0;count<fi;count++)
{
xVx_sorted[count]=xvs[count];
index[count]	 =count;
}
	
selection_sort(xVx_sorted,fi,index);
	
for(int count=0;count<fi;count++)
{
xUcbcc_sorted[count]=xUcbcc[index[count]];
xUxbcc_sorted[count]=xUxbcc[index[count]];
xUcfcc_sorted[count]=xUcfcc[index[count]];
xUxfcc_sorted[count]=xUxfcc[index[count]];
xgxtr_sorted[count]	=xgstr[index[count]];
xfm_sorted[count]	=xfvs[index[count]];
xmuxc_sorted[count]	=xmusx[index[count]];
xmuxx_sorted[count]	=xmusx[index[count]];
}	
std::vector <double> xgxtr_shifted(fi), xfm_shifted(fi);
for(int count=0;count<fi;count++)
{
xgxtr_shifted[count]	=xgxtr_sorted[count];
xfm_shifted[count]		=xfm_sorted[count];
}
	
roll(xgxtr_sorted,xgxtr_shifted);
roll(xfm_sorted,xfm_shifted);

std::vector <int> cross(fi), icross(fi);
double expr;
for(int count=0; count<fi; count++)
{
expr=(xfm_sorted[count]-xgxtr_sorted[count])* \
(xfm_shifted[count]-xgxtr_shifted[count]);
if(expr<=0.0){cross[count]=1;}
else cross[count]=0;
}

double temp,temp_cbcc,temp_cfcc,Vtmp,temp_nibcc, temp_nifcc;
std::vector <double> Vcross(fi),Ucbcc_cross(fi),Ucfcc_cross(fi), \
Unibcc_cross(fi),Unifcc_cross(fi);
int counter=0;
for(int count=1; count<fi; count++)
{
if(cross[count]==1)
{
temp = (xgxtr_sorted[count+1]-xfm_sorted[count+1]) \
	-(xgxtr_sorted[count]-xfm_sorted[count]);
if(temp==0.)
	temp=1.e-15;
Vtmp=xVx_sorted[count] \
	+(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(xgxtr_sorted[count]-xfm_sorted[count]) \
	/temp;
Vcross[count]=xVx_sorted[count];///Vtmp;
			
temp_cbcc=xUcbcc_sorted[count+1]-xUcbcc_sorted[count];
if(temp_cbcc==0.)
	temp_cbcc=0.;
else
temp_cbcc=(xUcbcc_sorted[count+1]-xUcbcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
Ucbcc_cross[count]=xUcbcc_sorted[count]+temp_cbcc;
			
temp_cfcc=xUcfcc_sorted[count+1]-xUcfcc_sorted[count];
if (temp_cfcc==0.)
	temp_cfcc=0.;
else
temp_cfcc=(xUcfcc_sorted[count+1]-xUcfcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
Ucfcc_cross[count]=xUcfcc_sorted[count]+temp_cfcc;
			
temp_nibcc=(xUxbcc_sorted[count+1]-xUxbcc_sorted[count]);
if(temp_nibcc==0.)
	temp_nibcc=0.;
else
temp_nibcc=(xUxbcc_sorted[count+1]-xUxbcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
Unibcc_cross[count]=xUxbcc_sorted[count]+temp_nibcc;
			
temp_nifcc=(xUxfcc_sorted[count+1]-xUxfcc_sorted[count]);
if(temp_nifcc==0.)
	temp_nifcc=0.;
else
temp_nifcc=(xUxfcc_sorted[count+1]-xUxfcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
Unifcc_cross[count]=xUxfcc_sorted[count]+temp_nifcc;
counter++;
}}

std::vector <double> Vcrossi(counter),Ucbcc_crossi(counter),Ucfcc_crossi(counter), \
Unibcc_crossi(counter),Unifcc_crossi(counter);
int counti=0;
for(int count=0; count<fi;count++)
{
if(Ucbcc_cross[count]!=0.)
{
Ucbcc_crossi[counti]=Ucbcc_cross[count];
Ucfcc_crossi[counti]=Ucfcc_cross[count];
Unibcc_crossi[counti]=Unibcc_cross[count];
Unifcc_crossi[counti]=Unifcc_cross[count];
Vcrossi[counti]=Vcross[count];				
counti++;
}	
}
printf("Solutions:\n");	
for(int count=0; count<counter; count++)
{
printf("%e %e %e %e %e\n",Ucbcc_crossi[count], \
Ucfcc_crossi[count],Unibcc_crossi[count], \
Unifcc_crossi[count],Vcrossi[count]);	
}
//-------------------------------------
//-------------------------------------
//-------------------------------------	
//Complet model	
std::vector <double> Ycinetique_comp(5);
//if(counter==0)
//{printf("No solutions \n");continue;}		
if(counter==3)
{
Ycinetique_comp[0]=Ucbcc_crossi[1]/3.;//C in ferrite
Ycinetique_comp[1]=Unibcc_crossi[1];//X in ferrite
Ycinetique_comp[2]=Ucfcc_crossi[1];//C in austenite
Ycinetique_comp[3]=Unifcc_crossi[1];//X in austenite
Ycinetique_comp[4]=Vcrossi[1];//Velocity fixe
}
else{	
Ycinetique_comp[0]=Ucbcc_crossi[2]/3.;//C in ferrite
Ycinetique_comp[1]=Unibcc_crossi[2];//X in ferrite
Ycinetique_comp[2]=Ucfcc_crossi[2];//C in austenite
Ycinetique_comp[3]=Unifcc_crossi[2];//X in austenite
Ycinetique_comp[4]=Vcrossi[2];//Velocity fixe
}
//exit(0);
bool flage_ch = false;	

	
//calc3cinetiqueFeMnC_complet(CONCENTRATION_NOMINAL,Ycinetique_comp,prevVectorcin,flage_ch);
//calc3cinetiqueFeMnC_para_or_LE_T(CONCENTRATION_NOMINAL,Ycinetique_comp,prevVectorcin,flage_ch);
	
}



/*bool calc3cinetiqueFeMnC_growth(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage)
{
newt nt;
nt=newt();
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT_prof5;
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT_prof5;
i++;
}
	
bool conv_para = false;
int j=0;
int compteur, diff_step;
double VIT_TEMPORAIRE0,raison,Delta_c,Delta_ni; 	
bool flage_growth_plus, flage_growth_minus;	
correc_D=1.;
	
int check=0;
nt.fvec.clear();	
nt.fvec1.clear();
nt.nn = 4;
std::vector <double> Ycinetique_veloc(nt.nn);
std::vector <double> Xcinetique(nt.nn);	
std::vector <double> Ucinetique(nt.nn);		
std::vector <double> Unom(2);
int compt = 0;	
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 
	
Ycinetique_veloc[0]=Y[2]; //C in ferrite
Ycinetique_veloc[1]=Y[1];//X in ferrite
Ycinetique_veloc[2]=Y[6];//C in austenite
Ycinetique_veloc[3]=Y[5];//X in austenite
VIT_TEMPORAIRE=Y[8];//Velocity fixe
	
	
diff_step=1000;
raison = pow( pow (10. , 8) , 1. / (double) diff_step ) ;	
compteur = 0 ;	
VIT_TEMPORAIRE0=VIT_TEMPORAIRE;
flage_growth_plus=false;
flage_growth_minus=false;	
	
std::vector < std::vector <double> > solutions ( 5 );	
double gtr_old , gtr2_old ;
double div ;	
for (int ii=0 ; ii<=4 ; ii++) solutions[ii] . resize(5) ;
std::vector <double> Y_para_temp ( 5 );
std::vector <double> Y_dernier ( 5 );
int 	nbsol = 0 ;	

	
while (compteur <= diff_step)
{
VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;
		
nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,1,conv_para);
if(!conv_para) {flage_growth_minus=false;break ;}
else flage_growth_minus=true;
if(conv_para == 1) Y_dernier = Y_para_temp ; ;
		
		
C_YatoXa(Ycinetique_veloc,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);
		
		
if(conv_para){	
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr;Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
			
so7.ecritureV ( Ua,compt );
	
if ( (gtr1 - gtr) * (gtr2_old - gtr_old) < 0. )
{
solutions[nbsol][0] = Ycinetique_veloc[0] ;
solutions[nbsol][1] = Ycinetique_veloc[1] ;
solutions[nbsol][2] = Ycinetique_veloc[2] ;
solutions[nbsol][3] = Ycinetique_veloc[3] ;
solutions[nbsol][4] = VIT_TEMPORAIRE ;
nbsol ++ ;
}
gtr2_old = gtr1 ; gtr_old = gtr ;	
			
}
Y_para_temp = Ycinetique_veloc ;
compteur ++;
compt++;	
}
	
printf("\nBoucle dans l'autre sens : augmentation de la vitesse\n") ;
	
Ycinetique_veloc[0]=Y[2]; //C in ferrite
Ycinetique_veloc[1]=Y[1];//X in ferrite	
Ycinetique_veloc[2]=Y[6];//C in austenite
Ycinetique_veloc[3]=Y[5];//X in austenite
VIT_TEMPORAIRE=Y[8];//Velocity fixe
	
	
compteur = 0 ;
	
while (compteur >= -diff_step)
{		
		
VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;
		
nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,1,conv_para);
if(!conv_para){ flage_growth_plus=false;break ;}
if(VIT_TEMPORAIRE>200.){break;}
else flage_growth_plus=true;
		
C_YatoXa(Ycinetique_veloc,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);
//printf(" %e  %e %e %e \n",Ucinetique[0],Ucinetique[1],Ucinetique[2],Ucinetique[3]);	
if(conv_para){
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr;Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
so7.ecritureV ( Ua,compt );
	
if ( (gtr1 - gtr) * (gtr2_old - gtr_old) < 0. )
{
solutions[nbsol][0] = Ycinetique_veloc[0] ;
solutions[nbsol][1] = Ycinetique_veloc[1] ;
solutions[nbsol][2] = Ycinetique_veloc[2] ;
solutions[nbsol][3] = Ycinetique_veloc[3] ;
solutions[nbsol][4] = VIT_TEMPORAIRE ;
nbsol ++ ;
}
gtr2_old = gtr1 ; gtr_old = gtr ;	
}
	
compteur --;
compt++;	
}
	
so7.fermeture();
//--------------------------------
//-------------------------------------
//-------------------------------------
//-------------------------------------
printf("\nNombre de solutions = %i",nbsol);
for(int ii = 1 ; ii < nbsol ; ii++) printf("\nVitesse no %i = %e" , ii , solutions[ii][4] ) ;

	//Classement des solutions par vitesse croissante
	printf("\nnbsol=%i",nbsol) ;
	if(nbsol == 3)
	{
		double vmin = 1.e30 ;
		double vmax = 0. ;
		std::vector <double> Ytemp1 ( 5 ) ;
		std::vector <double> Ytemp2 ( 5 ) ;
		std::vector <double> Ytemp3 ( 5 ) ;
		
		for(int ii=0 ; ii<=2 ; ii++)
		{
			if(solutions[ii][4] < vmin) { vmin = solutions[ii][4] ; Ytemp1 = solutions[ii] ; }
			if(solutions[ii][4] > vmax) { vmax = solutions[ii][4] ; Ytemp3 = solutions[ii] ; }		 
		}
		
		for(int ii=0 ; ii<=2 ; ii++)
			if(solutions[ii][4] > vmin && solutions[ii][4] < vmax) Ytemp2 = solutions[ii] ;
		
		solutions[0] = Ytemp1 ;
		solutions[1] = Ytemp2 ;
		solutions[2] = Ytemp3 ;
	}

	
Ycinetique_veloc = solutions[1] ;
printf("Solutions:\n");
printf("%e %e %e %e %e \n",solutions[0][4],solutions[1][4],solutions[2][4],solutions[3][4],solutions[4][4]);
printf("\n%e %e %e %e %e \n",Ycinetique_veloc[0],Ycinetique_veloc[1],Ycinetique_veloc[2],Ycinetique_veloc[3],Ycinetique_veloc[4]);	
//exit(0);
bool flage_ch = false;	
calc3cinetiqueFeMnC_complet(CONCENTRATION_NOMINAL,Ycinetique_veloc,prevVectorcin,flage_ch);	
}*/
	

bool calc3cinetiqueFeMnC_complet(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage_ch)
{	
//std::vector <double> interface::deriv_mu_FeCCr(double yc,double ycr,double x,vector <double> &concentration)
//	Diffbetta[0] = dc*coeff_n;
//	Diffbetta[1] = dcr*coeff_k*correct_D; and correct_D = 1
//void interface::KTT_prof(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
//	F[3] = (Uc_gamma-Uc_alpha)*Y[4]-(Diffbetta[0]/( Rzin*delta[0]))* \
(Uc_gamma-U0c); //C (gamma->alpha)
//F[4] = (Ucr_gamma-Ucr_alpha)*Y[4]-(Diffbetta[1]/(Rzin*delta[1]))* \
(Ucr_gamma-U0cr); //CR (gamma ->alpha)
printf("Complete model : \n\n\n");	
newt nt2;
nt2=newt();

int i = 1;
while ( mat.produit[i]!=NULL )
{
nt2.func     =mat.produit[1]->frontiere[1];
nt2.nrfuncv  =&interface::KTT_prof;
nt2.jacobian =mat.produit[1]->frontiere[1];
nt2.jacobfunc=&interface::JKT_prof;
i++;
}
if ( nt2.func == NULL ) printf ( "newton.func = null\n" );
if ( nt2.jacobian==NULL) printf ("jacobian.func=null\n");
	
int check=0;
nt2.fvec.clear();
nt2.fvec1.clear();
nt2.nn = 5;
std::vector <double> Ycinetique_para(nt2.nn);
std::vector <double> Ycinetique_para_temp(nt2.nn);	
std::vector <double> prevVectorcin_para(nt2.nn);
	
nt2.fvec.resize(nt2.nn);
nt2.fvec1.resize(nt2.nn); 
int j=0;
double before=0.;
double after=0.;
double before_v = 0.;
prevVectorcin_para[0] = Y[0];//c_alpha	
prevVectorcin_para[1] = Y[1];//cr_alpha	
prevVectorcin_para[2] = Y[2];//c_gamma
prevVectorcin_para[3] = Y[3];//cr_gamma
prevVectorcin_para[4] = Y[4];

//coeff_k=1.e11;//1e6; //(Dni)
//correc_D= 1.e-1;
//coeff_n = 1.e13;


while(tempsFeCCr<= 1.e30)	
{
if(!flage_ch)
{	
Ycinetique_para[0]=prevVectorcin_para[0];//c_alpha	
Ycinetique_para[1]=prevVectorcin_para[1];//cr_alpha	
Ycinetique_para[2]=prevVectorcin_para[2];//c_gamma
Ycinetique_para[3]=prevVectorcin_para[3];//cr_gamma
Ycinetique_para[4]=prevVectorcin_para[4];
VIT_TEMPORAIRE	= Ycinetique_para[4];
}
else
{
//dtempsFeCCr=0.1;			
Ycinetique_para[0]=prevVectorcin_para[0];//c_alpha	
Ycinetique_para[1]=prevVectorcin_para[1];//cr_alpha	
Ycinetique_para[2]=prevVectorcin_para[2];//c_gamma
Ycinetique_para[3]=prevVectorcin_para[3];//cr_gamma
Ycinetique_para[4]=prevVectorcin_para[4];
VIT_TEMPORAIRE	= Ycinetique_para[4];

	



	
	
double incr ;
incr = 1.e-1 ;
if(Ycinetique_para[4]<=1.e-9)	
{
incr =Ycinetique_para[4]*1e-3;//
}
	
dtempsFeCCr = Rzin * incr / Ycinetique_para[4] ;


before=Rzin;		
Rzin=before+Ycinetique_para[4]*dtempsFeCCr;
tempsFeCCr += dtempsFeCCr ;
	
//std::vector <double> Xcinetique_t(5),Unom_fin(2);
//std::vector <double> Ucinetique_t(4),Unom_t(2),U_fin(7);
//C_YatoXa(Ycinetique_para,Xcinetique_t);	
//C_XatoUa(Xcinetique_t,Ucinetique_t);
//C_X0toU0(CONCENTRATION_NOMINAL,Unom_t);	
//U_fin[0]=Ucinetique_t[0];U_fin[1]=Ucinetique_t[1];	
//U_fin[2]=Ucinetique_t[2];U_fin[3]=Ucinetique_t[3];
//U_fin[4]=Ycinetique_para[4];U_fin[5]=Unom_t[1];
//U_fin[6]=Unom_t[0];
//printf("Bilan de soluté \n");	
//Unom_fin = BILAN_MASSE_PLANE_U(U_fin,dtempsFeCCr);	
//CONCENTRATION_NOMINAL[0] = Unom_fin[1]/(1.+Unom_fin[1]);
//CONCENTRATION_NOMINAL[2] = Unom_fin[0]/(1.+Unom_fin[1]);
//printf("Unom_ni =%e Unom_c = %e \n",Unom_fin[0],Unom_fin[1]);
//printf("Xnom_ni =%e Cnom_c = %e \n",Unom_fin[0]/(1.+Unom_fin[1]),Unom_fin[1]/(1.+Unom_fin[1]));

//exit(0);	
//CONCENTRATION_NOMINAL[0]+=1e-6;
}

bool conv_para = false;
std::vector <double> Xcinetique(5);
std::vector <double> Ucinetique(4),U(15),Unom(2);

printf("\n\n");
printf("U(C,BCC) = %e U(C,FCC)=%e \n",Ycinetique_para[0]*3.,Ycinetique_para[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)=%e \n",Ycinetique_para[1],Ycinetique_para[3]);
printf("U(C0) = %e  U(X0) = %e",CONCENTRATION_NOMINAL[0]/(1.-CONCENTRATION_NOMINAL[0]), \
						CONCENTRATION_NOMINAL[2]/(1.-CONCENTRATION_NOMINAL[0]));	
printf("v = %e \n",Ycinetique_para[4]);
printf("correc_D =%e coeff_n = %e coeff_k = %e VITTEMPORAIRE = %e \n",correc_D,coeff_n,coeff_k,VIT_TEMPORAIRE);
printf("dc = %e dni = %e \n",dc,dcr);
printf("T = %e K \n",T);
printf("Rzin = %e temps = %e \n",Rzin,tempsFeCCr);	
muc = 0.0; gtr = 0.0; mucr = 0.0;
	
printf("Input Y in NR Dc/Dcr:\n");	
print_vect(Ycinetique_para);
//exit(0);	
before_v = 	Ycinetique_para[4];
nt2.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);

printf("Output Y from NR:\n");	
print_vect(Ycinetique_para);	

C_YatoXa(Ycinetique_para,Xcinetique);	
Xcinetique[4] = Ycinetique_para[4];	
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
	
printf("\n");	
printf("Solution in mole fraction\n");
printf("X(C,BCC) = %e  X(C,FCC)= %e \n",Xcinetique[0],Xcinetique[1]);
printf("X(Mn,BCC) = %e  X(Mn,FCC)= %e \n",Xcinetique[2],Xcinetique[3]);
printf("XMn = %e  XC= %e \n",CONCENTRATION_NOMINAL[2],CONCENTRATION_NOMINAL[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("\n\n");
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",Ycinetique_para[4]);


	
if(conv_para)
{
//std::vector <double> Xcinetique(5);
//std::vector <double> Ucinetique(4),U(15),Unom(2);
	
C_YatoXa(Ycinetique_para,Xcinetique);	
Xcinetique[4] = Ycinetique_para[4];	
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);

printf("\n");	
printf("Solution in mole fraction\n");
printf("X(C,BCC) = %e  X(C,FCC)= %e \n",Xcinetique[0],Xcinetique[1]);
printf("X(Mn,BCC) = %e  X(Mn,FCC)= %e \n",Xcinetique[2],Xcinetique[3]);
printf("XMn = %e  XC= %e \n",CONCENTRATION_NOMINAL[2],CONCENTRATION_NOMINAL[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("\n\n");
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("muc = %e mux= %e\n",muc/(R*T),mucr/(R*T),gtr/(R*T));
printf("R = %e temps = %e \n",Rzin,tempsFeCCr);	
prevVectorcin_para[0]=Ycinetique_para[0];//c_alpha	
prevVectorcin_para[1]=Ycinetique_para[1];//cr_alpha	
prevVectorcin_para[2]=Ycinetique_para[2];//c_gamma
prevVectorcin_para[3]=Ycinetique_para[3];//cr_gamma
prevVectorcin_para[4]=Ycinetique_para[4];
	
	
	
	
	
// Uni_fcc				//Uni_bcc
U[0]=Ucinetique[0];U[1]=Ucinetique[1];
//Uc_fcc				//Uc_bcc
U[2]=Ucinetique[2];U[3]=Ucinetique[3];
//muc			//mux			//gtr
U[4]=muc/(R*T);U[5]=mucr/(R*T);U[6]=gtr/(R*T);
//v						//r			//Delta_ni
U[7]=Ycinetique_para[4];U[8]=Rzin;U[9]=Sursaturation_ni(Ucinetique,Unom);
//Delta_c							//time
U[10]=Sursaturation_c(Ucinetique,Unom);U[11]=tempsFeCCr;	
//U0ni			//U0c			//Dni/Dni
U[12]=Unom[1];U[13]=Unom[0];U[14]=coeff_k;	
so5.ecritureX(Xcinetique);
so4.ecritureU(U);
flage_ch = true;	
}
else{
flage_ch = false;	
printf("Searching solution for new r= %e time=%e\n",Rzin,tempsFeCCr);
std::vector <double> Ycinetique_para_ch(5);	
cinetiqueLE_ch_FeCMn(CONCENTRATION_NOMINAL,prevVectorcin,Ycinetique_para_ch,flage_ch);	
//cinetiqueLE_qh_FeCMn(CONCENTRATION_NOMINAL,prevVectorcin,Ycinetique_para_ch,flage_ch);	
prevVectorcin_para[0]=Ycinetique_para_ch[0];//c_alpha	
prevVectorcin_para[1]=Ycinetique_para_ch[1];//cr_alpha	
prevVectorcin_para[2]=Ycinetique_para_ch[2];//c_gamma
prevVectorcin_para[3]=Ycinetique_para_ch[3];//cr_gamma
prevVectorcin_para[4]=Ycinetique_para_ch[4];
printf("Results:\n");
for(int na=0; na<5;na++)printf("%e \n",prevVectorcin_para[na]);
}
j++;
}

}



bool cinetiqueLE_ch_FeCMn(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &Ycinetique_para_ch, bool &flage_ch)
{
newt nt;
nt=newt();
int i=1;
while (mat.produit[i]!=NULL)
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT;
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT;
i++;
}
if (nt.func == NULL) printf ("newton.func = null\n");
if (nt.jacobian==NULL) printf ("jacobian.func=null\n");	
	
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn=9;
std::vector <double> Ycinetique(nt.nn);
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 
double before=0.0; double after=0.0; float v_ini;
int j=0; int diff_step; int compteur=0;
bool conv = false;
FILE *fpp=fopen("coeff_cin.txt","r") ;
fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
fscanf(fpp , "%i  \n" , &diff_step) ;	
fclose(fpp);	
double raison=pow(dc/dcr,1./((double)diff_step)) ;	
//Principle
printf("Kinetic Fe-C-Mn with local equilibrium at interface \n\n\n");	
	
	//prevVectorcin_para[0] = Y[2];//C in bcc
	//prevVectorcin_para[1] =	Y[1];//X in bcc
	//prevVectorcin_para[2] = Y[6];//C in fcc
	//prevVectorcin_para[3] = Y[5];//X in fcc
	//prevVectorcin_para[4] = 1.e-9;//v
	
Ycinetique[0] = 1.-Y[1];//FE
Ycinetique[1] = Y[1];//CR
Ycinetique[2] = Y[0];//C
Ycinetique[3] = 1.-Y[0];//VA
Ycinetique[4] = 1.-Y[3];//FE
Ycinetique[5] =	Y[3];//CR
Ycinetique[6] = Y[2];//C
Ycinetique[7] = 1.-Y[2];//VA
Ycinetique[8] = Y[4];
correc_D=0;
compteur=0;
while(compteur<= diff_step)
{			
correc_D = pow(raison,compteur);
printf("Input Y in NR:\n");	
print_vect(Ycinetique);
nt.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
printf("Output Y from NR:\n");	
print_vect(Ycinetique);
printf("diff_step= %i \n",diff_step);
printf("compt=%i \n",compteur);
printf("correc_D=%e \n", correc_D);	
printf("dc=%e dni=%e dni=%e \n",dc,dcr,correc_D*dcr);
compteur ++;
}
if(!conv) 
{
printf("Model LE : No convergence \n temps= %e dt = %e \n\n",tempsFeCCr,dtempsFeCCr);
exit(0);
}
else{
calc3cinetiqueFeMnC_growth_ch(CONCENTRATION_NOMINAL,Ycinetique,Ycinetique_para_ch,flage_ch);	
}	

}

bool calc3cinetiqueFeMnC_growth_ch(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &Ycinetique_para_ch,bool &flage)
{
printf("V\n");	
newt nt;
nt=newt();
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT_prof5;
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT_prof5;
i++;
}
bool conv_para = false;
int j=0;
int compteur, diff_step;
double VIT_TEMPORAIRE0,raison,Delta_c,Delta_ni; 	
bool flage_growth_plus, flage_growth_minus;	
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 4;
std::vector <double> Ycinetique_veloc(nt.nn);
std::vector <double> Xcinetique(nt.nn);	
std::vector <double> Ucinetique(nt.nn);		
std::vector <double> Unom(2);
int compt = 0;	
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 
Ycinetique_veloc[0]=Y[2]; //C in ferrite
Ycinetique_veloc[1]=Y[1];//X in ferrite
Ycinetique_veloc[2]=Y[6];//C in austenite
Ycinetique_veloc[3]=Y[5];//X in austenite
VIT_TEMPORAIRE=Y[8];//Velocity fixe
diff_step=1000;
raison = pow( pow (10. , 8) , 1. / (double) diff_step ) ;	
compteur = 0 ;	
VIT_TEMPORAIRE0=VIT_TEMPORAIRE;
flage_growth_plus=false;
flage_growth_minus=false;	
correc_D = 1.0;		
while (compteur <= diff_step)
{
VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;
		
//nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,7,conv_para);
nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,1,conv_para);
	
if(!conv_para) {flage_growth_minus=false;break ;}
else flage_growth_minus=true;
		
C_YatoXa(Ycinetique_veloc,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
	
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);
		
if(conv_para){	
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr;Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
so8.ecritureV ( Ua,compt );
}
compteur ++;
compt++;	
}
printf("\nBoucle dans l'autre sens : augmentation de la vitesse\n") ;
Ycinetique_veloc[0]=Y[2]; //C in ferrite
Ycinetique_veloc[1]=Y[1];//X in ferrite
Ycinetique_veloc[2]=Y[6];//C in austenite
Ycinetique_veloc[3]=Y[5];//X in austenite
VIT_TEMPORAIRE=Y[8];//Velocity fixe
compteur = 0 ;
while (compteur >= -diff_step)
{		
VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;
		
nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,1,conv_para);
if(!conv_para){ flage_growth_plus=false;break ;}
if(VIT_TEMPORAIRE>200.){break;}
else flage_growth_plus=true;
		
C_YatoXa(Ycinetique_veloc,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);
//printf(" %e  %e %e %e \n",Ucinetique[0],Ucinetique[1],Ucinetique[2],Ucinetique[3]);	
if(conv_para){
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr;Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
so8.ecritureV ( Ua,compt );	
}
compteur --;
compt++;	
}
so8.fermeture();
//--------------------------------
//--------------------------------
//--------------------------------
//Open the file with all solution
FILE *fee;
char mystring [100];
//sprintf(mystring,"LENP/%eV.txt",T);
sprintf(mystring,"V1.txt");
cout<<mystring<<endl;	
fee=fopen(mystring,"rt");
printf("open \n");	
if(!fee){
printf("Error open file \n");
return 0;}
else
{ 
long lSize;
fseek(fee , 0 , SEEK_END);
lSize = ftell(fee);
rewind(fee);
if(lSize==0)
{fclose (fee);
printf("File is empty\n");
exit(0);	
return 0;
}}if(fee) printf("File is opened\n");
char str[256];
int nb_lignes = 0;
int nm = 0;
bool wasNumbers = false;	
int nx=1000000;
std::vector <float> xTs(nx),xUxfcc(nx),xUxbcc(nx), \
xUcbcc(nx),xUcfcc(nx),xmusc(nx),xmusx(nx),xgstr(nx), \
xfvs(nx),xts(nx), \
xU0sx(nx),xU0sc(nx),xcompts(nx), \
xrs(nx),xDeltasx(nx), \
xvs(nx),xDeltasc(nx),xDsni(nx);
float Ts,Uxfcc,Uxbcc, \
Ucbcc,Ucfcc,musc,musx,gstr, \
vs,fvs,ts, \
U0sx,U0sc,rs,Deltasx,Deltasc,Dsni;
int compts;
	
//While file with results is opened we search the solution	
while (!feof(fee))
{
str[0] = 0;
char *tmp = fgets(str, sizeof(str), fee);
if (!tmp || str[0] == 0) //If have found nothing,
 break; //we leave
		
if (str[0] == '/') //Comment
{
if (wasNumbers)//If we considered numbers from a matrix
{
wasNumbers = false;
++nm; //Means we pass to the second matrix
}continue;
}
//printf("File name: %s\n",str);
int sp = calc_spaces3(str); //We consider quantity of blanks
	
if(nm==0)
{
//"#T  U(Mn,fcc) U(Mn,bcc) U(C,fcc) U(C,bcc) muc mucx gtr v Fv r Delta_x time Delta_c U0x U0c Dx/Dx"<<endl;	
sscanf(str, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",&Ts,&Uxfcc,&Uxbcc, \
		&Ucfcc,&Ucbcc,&musc,&musx,&gstr, \
		&vs,&fvs,&rs,&Deltasx,&ts,&Deltasc, \
		&U0sx,&U0sc,&Dsni,&compts);
//printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d\n",Ts,Uxfcc,Uxbcc, \
Ucfcc,Ucbcc,musc,musx,gstr, \
vs,fvs,rs,Deltasx,ts,Deltasc, \
U0sx,U0sc,Dsni,compts);
			
xTs[compts]		= Ts;
xUxfcc[compts]	= Uxfcc;
xUxbcc[compts]	= Uxbcc;
xUcbcc[compts]	= Ucbcc;
xUcfcc[compts]	= Ucfcc;
xmusc[compts]	= musc;
xmusx[compts]	= musx;
xgstr[compts]	= gstr;
xvs[compts]		= vs;
xfvs[compts]	= fvs;
xrs[compts]		= rs;
xDeltasx[compts]= Deltasx;
xts[compts]		= ts;
xDeltasc[compts]= Deltasc;
xU0sx[compts]	= U0sx;
xU0sc[compts]	= U0sc;
xDsni[compts]	= Dsni;
xcompts[compts] = compts;	
wasNumbers = true;
}
}
int fi = compts+1;
xTs.resize(fi);xUxfcc.resize(fi);xUxbcc.resize(fi);
xUcbcc.resize(fi);xUcfcc.resize(fi);xmusc.resize(fi);
xmusx.resize(fi);xgstr.resize(fi);
xvs.resize(fi);xfvs.resize(fi);xrs.resize(fi);
xDeltasx.resize(fi);xts.resize(fi);xDeltasc.resize(fi);
xU0sx.resize(fi);xU0sc.resize(fi);xDsni.resize(fi);
xcompts.resize(fi);
fclose(fee);
std::vector <int> index(fi);
std::vector <double> xVx_sorted(fi),xUcbcc_sorted(fi),xUxbcc_sorted(fi), \
xUcfcc_sorted(fi), xUxfcc_sorted(fi), xUxc_sorted(fi), \
xUxx_sorted(fi), xgxtr_sorted(fi), \
xfm_sorted(fi), xmuxc_sorted(fi), \
xmuxx_sorted(fi);
	
for(int count=0;count<fi;count++)
{
xVx_sorted[count]=xvs[count];
index[count]	 =count;
}
	
selection_sort(xVx_sorted,fi,index);

for(int count=0;count<fi;count++)
{
xUcbcc_sorted[count]=xUcbcc[index[count]];
xUxbcc_sorted[count]=xUxbcc[index[count]];
xUcfcc_sorted[count]=xUcfcc[index[count]];
xUxfcc_sorted[count]=xUxfcc[index[count]];
xgxtr_sorted[count]	=xgstr[index[count]];
xfm_sorted[count]	=xfvs[index[count]];
xmuxc_sorted[count]	=xmusx[index[count]];
xmuxx_sorted[count]	=xmusx[index[count]];
}	
std::vector <double> xgxtr_shifted(fi), xfm_shifted(fi);
for(int count=0;count<fi;count++)
{
xgxtr_shifted[count]	=xgxtr_sorted[count];
xfm_shifted[count]		=xfm_sorted[count];
}
	
roll(xgxtr_sorted,xgxtr_shifted);	
roll(xfm_sorted,xfm_shifted);
	
std::vector <int> cross(fi), icross(fi);
double expr;
for(int count=0; count<fi; count++)
{
expr=(xfm_sorted[count]-xgxtr_sorted[count])* \
(xfm_shifted[count]-xgxtr_shifted[count]);
if(expr<=0.0){cross[count]=1;}
else cross[count]=0;
}
	
double temp,temp_cbcc,temp_cfcc,Vtmp,temp_nibcc, temp_nifcc;
std::vector <double> Vcross(fi),Ucbcc_cross(fi),Ucfcc_cross(fi), \
Unibcc_cross(fi),Unifcc_cross(fi);
int counter=0;
for(int count=1; count<fi; count++)
{
if(cross[count]==1)
{
temp = (xgxtr_sorted[count+1]-xfm_sorted[count+1]) \
	-(xgxtr_sorted[count]-xfm_sorted[count]);
if(temp==0.)
 temp=1.e-15;
Vtmp=xVx_sorted[count] \
	+(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(xgxtr_sorted[count]-xfm_sorted[count]) \
	/temp;
Vcross[count]=Vtmp;
			
temp_cbcc=xUcbcc_sorted[count+1]-xUcbcc_sorted[count];
if(temp_cbcc==0.)
 temp_cbcc=0.;
else
temp_cbcc=(xUcbcc_sorted[count+1]-xUcbcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
Ucbcc_cross[count]=xUcbcc_sorted[count]+temp_cbcc;
			
temp_cfcc=xUcfcc_sorted[count+1]-xUcfcc_sorted[count];
if (temp_cfcc==0.)
	temp_cfcc=0.;
else
temp_cfcc=(xUcfcc_sorted[count+1]-xUcfcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
Ucfcc_cross[count]=xUcfcc_sorted[count]+temp_cfcc;
			
temp_nibcc=(xUxbcc_sorted[count+1]-xUxbcc_sorted[count]);
if(temp_nibcc==0.)
	temp_nibcc=0.;
else
temp_nibcc=(xUxbcc_sorted[count+1]-xUxbcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
	Unibcc_cross[count]=xUxbcc_sorted[count]+temp_nibcc;
			
temp_nifcc=(xUxfcc_sorted[count+1]-xUxfcc_sorted[count]);
if(temp_nifcc==0.)
	temp_nifcc=0.;
else
temp_nifcc=(xUxfcc_sorted[count+1]-xUxfcc_sorted[count]) \
	/(xVx_sorted[count+1]-xVx_sorted[count]) \
	*(Vtmp-xVx_sorted[count]);
Unifcc_cross[count]=xUxfcc_sorted[count]+temp_nifcc;
counter++;
}}
std::vector <double> Vcrossi(counter),Ucbcc_crossi(counter),Ucfcc_crossi(counter), \
Unibcc_crossi(counter),Unifcc_crossi(counter);
int counti=0;
for(int count=0; count<fi;count++)
{
if(Ucbcc_cross[count]!=0.)
{
Ucbcc_crossi[counti]=Ucbcc_cross[count];
Ucfcc_crossi[counti]=Ucfcc_cross[count];
Unibcc_crossi[counti]=Unibcc_cross[count];
Unifcc_crossi[counti]=Unifcc_cross[count];
Vcrossi[counti]=Vcross[count];				
counti++;
}	
}
printf("Solutions:\n");	
for(int count=0; count<counter; count++)
{printf("%e %e %e %e %e\n",Ucbcc_crossi[count], \
Ucfcc_crossi[count],Unibcc_crossi[count], \
Unifcc_crossi[count],Vcrossi[count]);	
}
//-------------------------------------
//-------------------------------------
//-------------------------------------
if(counter==1)
{
Ycinetique_para_ch[0]=Ucbcc_crossi[0]/3.;
Ycinetique_para_ch[1]=Unibcc_crossi[0];
Ycinetique_para_ch[2]=Ucfcc_crossi[0];
Ycinetique_para_ch[3]=Unifcc_crossi[0];
Ycinetique_para_ch[4]=Vcrossi[0];
}
if(counter==2)
{
Ycinetique_para_ch[0]=Ucbcc_crossi[1]/3.;
Ycinetique_para_ch[1]=Unibcc_crossi[1];
Ycinetique_para_ch[2]=Ucfcc_crossi[1];
Ycinetique_para_ch[3]=Unifcc_crossi[1];
Ycinetique_para_ch[4]=Vcrossi[1];
}	
if(counter==3)
{
Ycinetique_para_ch[0]=Ucbcc_crossi[2]/3.;
Ycinetique_para_ch[1]=Unibcc_crossi[2];
Ycinetique_para_ch[2]=Ucfcc_crossi[2];
Ycinetique_para_ch[3]=Unifcc_crossi[2];
Ycinetique_para_ch[4]=Vcrossi[2];
	
/*Ycinetique_para_ch[0]=Ucbcc_crossi[1]/3.;
Ycinetique_para_ch[1]=Unibcc_crossi[1];
Ycinetique_para_ch[2]=Ucfcc_crossi[1];
Ycinetique_para_ch[3]=Unifcc_crossi[1];
Ycinetique_para_ch[4]=Vcrossi[1];*/
}

}
	


bool cinetiqueLE_qh_FeCMn(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y1,std::vector<double> &Ycinetique_para_ch, bool &flage_ch)
{
printf("Equilibrium\n");	
//The difinition dimension of alloy
int dim=3;
def=5;
//Definition Newton-Raphson and convex hull	
newt nt;
qhull qh(dim);
//Temporaire parameters	
x_i=CONCENTRATION_NOMINAL[0]; //convex hull
x_j=CONCENTRATION_NOMINAL[2];	
	
bool prevExts=false;
bool good;
bool flage = false;
			
qh=qhull(dim);
nt=newt();
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func=mat.produit[1]->frontiere[1];
nt.nrfuncv= &interface::ST;
nt.jacobian=mat.produit[1]->frontiere[1];
nt.jacobfunc = &interface::JT;
qh.func1    =mat.produit[1]->frontiere[1];
qh.nrfuncv1 =&interface::phase1_ternaire;
qh.func2    =mat.produit[1]->frontiere[1];
qh.nrfuncv2 =&interface::phase2_ternaire;
qh.Jacob1=mat.produit[1]->frontiere[1];
qh.nrJacob1=&interface::JacobianQhull3Dalfa;
qh.Jacob2=mat.produit[1]->frontiere[1];
qh.nrJacob2=&interface::JacobianQhull3Dbetta;
i++;
}
if ( nt.func == NULL ) printf ( "newton.func = null\n" );
if ( qh.func1== NULL ) printf ( "qhull.func = null\n" );
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 9;
std::vector <double> Y ( nt.nn );
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn ); //Jacobian
printf ( "T = %lf \n", T );
mat.Thermo->read ( T );
cem.classe0->Thermo->read ( T );

qh.discretization (T,55);
qh.build_hull();
vector<point> ext=qh.extremums(/*2.e-5*/);
				
if(ext.empty()) return 0;
if ( ( !ext.empty() ) )
{
printf ( "extremums:" );
for ( int i=0;i<ext.size();i++ )
cout<<i+1<<":"<<ext[i]<<endl;
}
else
{
printf ( "No extremums T = %lf\n",T );
return 0;
}
cout<<"Current extremums:"<<ext[0]<<";"<<ext[1]<<endl;
Y[1] = ( ext[0].func()?ext[1][1]:ext[0][1] );  //Mn
if ( Y[1]==0. )
Y[1]=1.e-4;
Y[0] = 1.-Y[1];       //FE
Y[2] = ( ext[0].func()?ext[1][0]:ext[0][0] );  //C
if ( Y[2]==0. )
Y[2]=1.e-4;
Y[3] = 1.-Y[2];       //VA
Y[5] =( ext[0].func()?ext[0][1]:ext[1][1] );   //Mn
if ( Y[5]==0. )
Y[5]=1.e-4;
Y[4] = 1.0-Y[5];      //FE
Y[6] = ( ext[0].func()?ext[0][0]:ext[1][0] );  //C
if ( Y[6]==0.0 )
Y[6]=1.e-4;
Y[7] = 1.-Y[6];       //VA
Y[8] = 0.45;			



printf("Input vector from Convex Hull");	
print_vect(Y);
int check=0;
bool conv = false;	
good=nt.calc ( Y, &check,CONCENTRATION_NOMINAL,3,conv);
for ( int l = 0; l < nt.nn; l++ )
{
if ( Y[l]<0.0 )
return 0;
}
printf("Output of NR");	
print_vect(Y);
			
/***************************************Cinetique********************************/
if(good)
{
std::vector <double> Yinput(5);
Yinput[0] = Y[2];//C
Yinput[1] = Y[1];//Mn
Yinput[2] = Y[6];//C
Yinput[3] = Y[5];//Mn
Yinput[4] = 1.e-9;			
printf("Calcul cinetique plane LE\n\n\n\n");
cinetiqueLE_ch_FeCMn(CONCENTRATION_NOMINAL,Yinput,Ycinetique_para_ch,flage_ch);
}
else {
	return 0;
}
	
}



bool calc3cinetiqueFeMnC_BalanceC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage)
{
newt nt;
nt=newt();

	
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT_prof1;
		
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT_prof1;

i++;
}
if ( nt.func == NULL ) printf ( "newton.func = null\n" );
if ( nt.jacobian==NULL) printf ("jacobian.func=null\n");
	
float v_ini;	
int diff_step;	
FILE *fpp=fopen("coeff_cin.txt","r") ;
fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
fscanf(fpp , "%i  \n" , &diff_step) ;	
fclose(fpp);
	
	
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 4;
std::vector <double> Ycinetique_para ( nt.nn );
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn ); 

double raison = pow ( dc / dcr , 1. / ( (double) diff_step ) ) ;
int ib,max;
max=5000;
double compteur;
bool conv_para = false;	
std::vector <double> Xcinetique(4);	
std::vector <double> Ucinetique(4);		
std::vector <double> Unom(2);
double Delta_c,Delta_ni; 	
ib = 0;	
int compt = 0;		
printf("U(C,BCC) increase \n");
while(ib<max)
{
y_fraction_c_bcc=0.0;
y_fraction_c_bcc  = (0.00054+ib*1.e-5)/(3.*(1.-(0.00054+ib*1.e-5)));//Y[2];
	
Ycinetique_para[0]=Y[1];//cr_alpha	
Ycinetique_para[1]=Y[6];//c_gamma
Ycinetique_para[2]=Y[5];//cr_gamma
Ycinetique_para[3]=Y[8];	
VIT_TEMPORAIRE	= Ycinetique_para[3];

compteur=diff_step;	
while(compteur>=-1./diff_step)
{
correc_D = pow ( raison , compteur ) ;
printf("Input Y in NR Dc/Dcr:\n");	
print_vect(Ycinetique_para);
nt.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
printf("Output Y in NR Dc/Dcr:\n");	
print_vect(Ycinetique_para);
compteur--;
}
if(!conv_para)break;	

//Mass balance equation
newt nt2;
nt2=newt();
int k = 1;
while ( mat.produit[k]!=NULL )
{
nt2.func     =mat.produit[1]->frontiere[1];
nt2.nrfuncv  =&interface::KTT_prof2;
nt2.jacobian =mat.produit[1]->frontiere[1];
nt2.jacobfunc=&interface::JKT_prof2;				
k++;
}	
int check2=0;
nt2.fvec.clear();
nt2.fvec1.clear();
nt2.nn = 1;
std::vector <double> Ycinetique_para2 ( nt2.nn );
nt2.fvec.resize ( nt2.nn );
nt2.fvec1.resize ( nt2.nn );	

VIT_TEMPORAIRE=Ycinetique_para[3];
y_fraction_x_bcc=Ycinetique_para[0];
y_fraction_c_fcc=Ycinetique_para[1];
y_fraction_x_fcc=Ycinetique_para[2];
bool conv_para2 = false;	

Ycinetique_para2[0]=CONCENTRATION_NOMINAL[2]/(1.-CONCENTRATION_NOMINAL[0]);
nt2.calc ( Ycinetique_para2, &check2,CONCENTRATION_NOMINAL,3,conv_para2);
	
if(Ycinetique_para[3]<0.0){printf("V<0\n");conv_para=false;break;}
if(!conv_para2)break;	
else 
{
C_YatoXa(Ycinetique_para,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
Unom[0]=CONCENTRATION_NOMINAL[0]/(1.-CONCENTRATION_NOMINAL[0]);
Unom[1]=Ycinetique_para2[0];	
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr/(R*T);Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
so7.ecritureV ( Ua,compt );	
}
compt++;	
ib++;	
}
exit(0);
	
ib = 0;	
printf("U(C,BCC) decreases \n");
while(ib<max)
{
y_fraction_c_bcc=0.0;
y_fraction_c_bcc  = (0.00054-ib*1.e-5)/(3.*(1.-(0.00054-ib*1.e-5)));//Y[2];
		
Ycinetique_para[0]=Y[1];//cr_alpha	
Ycinetique_para[1]=Y[6];//c_gamma
Ycinetique_para[2]=Y[5];//cr_gamma
Ycinetique_para[3]=Y[8];	
VIT_TEMPORAIRE	= Ycinetique_para[3];
		
compteur=diff_step;	
while(compteur>=-1./diff_step)
{
correc_D = pow ( raison , compteur ) ;
printf("Input Y in NR Dc/Dcr:\n");	
print_vect(Ycinetique_para);
nt.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
printf("Output Y in NR Dc/Dcr:\n");	
print_vect(Ycinetique_para);
compteur--;
}
//Mass balance equation
newt nt2;
nt2=newt();
int k = 1;
while ( mat.produit[k]!=NULL )
{
nt2.func     =mat.produit[1]->frontiere[1];
nt2.nrfuncv  =&interface::KTT_prof2;
nt2.jacobian =mat.produit[1]->frontiere[1];
nt2.jacobfunc=&interface::JKT_prof2;				
k++;
}	
int check2=0;
nt2.fvec.clear();
nt2.fvec1.clear();
nt2.nn = 1;
std::vector <double> Ycinetique_para2 ( nt2.nn );
nt2.fvec.resize ( nt2.nn );
nt2.fvec1.resize ( nt2.nn );	
		
VIT_TEMPORAIRE=Ycinetique_para[3];
y_fraction_x_bcc=Ycinetique_para[0];
y_fraction_c_fcc=Ycinetique_para[1];
y_fraction_x_fcc=Ycinetique_para[2];
bool conv_para2 = false;	
		
Ycinetique_para2[0]=CONCENTRATION_NOMINAL[2]/(1.-CONCENTRATION_NOMINAL[0]);
nt2.calc ( Ycinetique_para2, &check2,CONCENTRATION_NOMINAL,3,conv_para2);
		
if(Ycinetique_para[3]<0.0){printf("V<0\n");conv_para=false;break;}
if(!conv_para2)break;	
if(!conv_para)break;	
else 
{
C_YatoXa(Ycinetique_para,Xcinetique);
C_XatoUa(Xcinetique,Ucinetique);
Unom[0]=CONCENTRATION_NOMINAL[0]/(1.-CONCENTRATION_NOMINAL[0]);
Unom[1]=Ycinetique_para2[0];	
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",VIT_TEMPORAIRE);
Delta_c = Sursaturation_c(Ucinetique,Unom);
Delta_ni = Sursaturation_ni(Ucinetique,Unom);
vector <double> Ua(16);
Ua[0]= Ucinetique[0];Ua[1]=Ucinetique[1];
Ua[2]=Ucinetique[2];Ua[3]=Ucinetique[3];
Ua[4]=muc/(R*T);Ua[5]=mucr/(R*T);
Ua[6]=gtr/(R*T);Ua[7]=VIT_TEMPORAIRE;
Ua[8]=gtr1;Ua[9]=Rzin;Ua[10]=Delta_c;
Ua[11]=tempsFeCCr;Ua[12]=Delta_ni;
Ua[13]=Unom[1];Ua[14]=Unom[0];Ua[15]=coeff_k;
so7.ecritureV ( Ua,compt );	
}
compt++;	
ib++;	
}
so7.fermeture();
}





bool calc3cinetiqueFeMnC_para_or_LE_T(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage_ch)
{	
//std::vector <double> interface::deriv_mu_FeCCr(double yc,double ycr,double x,vector <double> &concentration)
//	Diffbetta[0] = dc*coeff_n;
//	Diffbetta[1] = dcr*coeff_k*correct_D; and correct_D = 1
//void interface::KTT_prof(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
//	F[3] = (Uc_gamma-Uc_alpha)*Y[4]-(Diffbetta[0]/( Rzin*delta[0]))* \
(Uc_gamma-U0c); //C (gamma->alpha)
//F[4] = (Ucr_gamma-Ucr_alpha)*Y[4]-(Diffbetta[1]/(Rzin*delta[1]))* \
(Ucr_gamma-U0cr); //CR (gamma ->alpha)
printf("Complete model : \n\n\n");	
newt nt;
nt=newt();
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT_prof;
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT_prof;
i++;
}
if ( nt.func == NULL ) printf ( "newton.func = null\n" );
if ( nt.jacobian==NULL) printf ("jacobian.func=null\n");
	
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 5;
std::vector <double> Ycinetique_para(nt.nn);
std::vector <double> Ycinetique_para_temp(nt.nn);	
std::vector <double> prevVectorcin_para(nt.nn);
	
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 
int j=0;
double before=0.;
double after=0.;
double before_v = 0.;

int diff_step = 1000;	
double raison = pow ( 0.1 , 1. / ( (double) diff_step ) ) ;
int compteur = 0;
correc_D=1.;


prevVectorcin_para[0] = Y[0];//c_alpha	
prevVectorcin_para[1] = Y[1];//cr_alpha	
prevVectorcin_para[2] = Y[2];//c_gamma
prevVectorcin_para[3] = Y[3];//cr_gamma
prevVectorcin_para[4] = Y[4];

while(compteur<= diff_step)
{
Ycinetique_para[0]=prevVectorcin_para[0];//c_alpha	
Ycinetique_para[1]=prevVectorcin_para[1];//cr_alpha	
Ycinetique_para[2]=prevVectorcin_para[2];//c_gamma
Ycinetique_para[3]=prevVectorcin_para[3];//cr_gamma
Ycinetique_para[4]=prevVectorcin_para[4];
VIT_TEMPORAIRE	= Ycinetique_para[4];
	
T = 1010. - (double) compteur / (double) diff_step * 1000. ;
	
//Input parameters:
dc =  6.8*1.e-13*pow(1e6,2.);	//odqvist
dcr = 3.5e-5 * exp( -2.8600e5 /R/ T )*pow(1e6,2.); //ni
	
printf("\nBoucle augmentation de la temperature   T = %e" , T) ;
mat.Thermo->read ( T ) ;
cem.classe0->Thermo->read ( T ) ;

bool conv_para = false;
std::vector <double> Xcinetique(5);
std::vector <double> Ucinetique(4),U(15),Unom(2);
printf("\n\n");
printf("U(C,BCC) = %e U(C,FCC)=%e \n",Ycinetique_para[0]*3.,Ycinetique_para[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)=%e \n",Ycinetique_para[1],Ycinetique_para[3]);
printf("v = %e \n",Ycinetique_para[4]);
printf("correc_D =%e coeff_n = %e coeff_k = %e \n",correc_D,coeff_n,coeff_k);
printf("dc = %e dni = %e \n",dc,dcr);
printf("T = %e K \n",T);
printf("Rzin = %e temps = %e \n",Rzin,tempsFeCCr);	
muc = 0.0; gtr = 0.0; mucr = 0.0;
printf("Input Y in NR Dc/Dcr:\n");	
print_vect(Ycinetique_para);
before_v = 	Ycinetique_para[4];

	nt.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);

printf("Output Y from NR:\n");	
print_vect(Ycinetique_para);	
		
C_YatoXa(Ycinetique_para,Xcinetique);	
Xcinetique[4] = Ycinetique_para[4];	
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
		
printf("\n");	
printf("Solution in mole fraction\n");
printf("X(C,BCC) = %e  X(C,FCC)= %e \n",Xcinetique[0],Xcinetique[1]);
printf("X(Mn,BCC) = %e  X(Mn,FCC)= %e \n",Xcinetique[2],Xcinetique[3]);
printf("XMn = %e  XC= %e \n",CONCENTRATION_NOMINAL[2],CONCENTRATION_NOMINAL[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("\n\n");
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",Ycinetique_para[4]);
		
if(conv_para)
{
C_YatoXa(Ycinetique_para,Xcinetique);	
Xcinetique[4] = Ycinetique_para[4];	
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("\n");	
printf("Solution in mole fraction\n");
printf("X(C,BCC) = %e  X(C,FCC)= %e \n",Xcinetique[0],Xcinetique[1]);
printf("X(Mn,BCC) = %e  X(Mn,FCC)= %e \n",Xcinetique[2],Xcinetique[3]);
printf("XMn = %e  XC= %e \n",CONCENTRATION_NOMINAL[2],CONCENTRATION_NOMINAL[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("\n\n");
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("muc = %e mux= %e\n",muc/(R*T),mucr/(R*T),gtr/(R*T));
printf("R = %e temps = %e \n",Rzin,tempsFeCCr);	
prevVectorcin_para[0]=Ycinetique_para[0];//c_alpha	
prevVectorcin_para[1]=Ycinetique_para[1];//cr_alpha	
prevVectorcin_para[2]=Ycinetique_para[2];//c_gamma
prevVectorcin_para[3]=Ycinetique_para[3];//cr_gamma
prevVectorcin_para[4]=Ycinetique_para[4];
			
// Uni_fcc				//Uni_bcc
U[0]=Ucinetique[0];U[1]=Ucinetique[1];
//Uc_fcc				//Uc_bcc
U[2]=Ucinetique[2];U[3]=Ucinetique[3];
//muc			//mux			//gtr
U[4]=muc/(R*T);U[5]=mucr/(R*T);U[6]=gtr/(R*T);
//v						//r			//Delta_ni
U[7]=Ycinetique_para[4];U[8]=Rzin;U[9]=Sursaturation_ni(Ucinetique,Unom);
//Delta_c							//time
U[10]=Sursaturation_c(Ucinetique,Unom);U[11]=tempsFeCCr;	
//U0ni			//U0c			//Dni/Dni
U[12]=Unom[1];U[13]=Unom[0];U[14]=coeff_k;	
so5.ecritureX(Xcinetique);
so4.ecritureU(U);
flage_ch = true;	
}
else{return 0;}
compteur ++ ;
}	

}



bool calc3cinetiqueFeMnC_para_or_LE(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage_ch)
{
	//	Diffbetta[1] = dcr*coeff_k;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	//F[3] = (Uc_gamma-Uc_alpha)*Y[4]-(Diffbetta[1]*correc_D/( Rzin*delta[0]))* \
	(Uc_gamma-U0c); //C (gamma->alpha)
printf("Complete model : \n\n\n");	
newt nt;
nt=newt();
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTT_prof;
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKT_prof;
i++;
}
if ( nt.func == NULL ) printf ( "newton.func = null\n" );
if ( nt.jacobian==NULL) printf ("jacobian.func=null\n");

int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 5;
std::vector <double> Ycinetique_para(nt.nn);
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 
int j=0;
double before=0.;
double after=0.;
double before_v = 0.;
int diff_step = 100;	
int compteur = 0;
double raison = pow ( dc / dcr , 1. / ( (double) diff_step ) ) ;

	
Ycinetique_para[0]=Y[2];//c_alpha	
Ycinetique_para[1]=Y[1];//cr_alpha	
Ycinetique_para[2]=Y[6];//c_gamma
Ycinetique_para[3]=Y[5];//cr_gamma
Ycinetique_para[4]=1e-9;//Y[8];
VIT_TEMPORAIRE	= Ycinetique_para[4];

bool conv_para = false;
std::vector <double> Xcinetique(5);
std::vector <double> Ucinetique(4),U(15),Unom(2);

while(compteur<= diff_step)
{
correc_D = pow ( raison , compteur ) ;
	nt.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
compteur ++;
if(!conv_para)break;
	
}

if(conv_para)
{
C_YatoXa(Ycinetique_para,Xcinetique);	
Xcinetique[4] = Ycinetique_para[4];	
C_XatoUa(Xcinetique,Ucinetique);
C_X0toU0(CONCENTRATION_NOMINAL,Unom);
printf("\n");	
printf("Solution in mole fraction\n");
printf("X(C,BCC) = %e  X(C,FCC)= %e \n",Xcinetique[0],Xcinetique[1]);
printf("X(Mn,BCC) = %e  X(Mn,FCC)= %e \n",Xcinetique[2],Xcinetique[3]);
printf("XMn = %e  XC= %e \n",CONCENTRATION_NOMINAL[2],CONCENTRATION_NOMINAL[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("\n\n");
printf("Solution in U-fraction\n");
printf("U(C,BCC) = %e  U(C,FCC)= %e \n",Ucinetique[3],Ucinetique[2]);
printf("U(Mn,BCC) = %e U(Mn,FCC)= %e \n",Ucinetique[1],Ucinetique[0]);
printf("UMn = %e  UC= %e \n",Unom[1],Unom[0]);
printf("v = %e \n",Ycinetique_para[4]);
printf("muc = %e mux= %e\n",muc/(R*T),mucr/(R*T),gtr/(R*T));
printf("R = %e temps = %e \n",Rzin,tempsFeCCr);	

// Uni_fcc				//Uni_bcc
U[0]=Ucinetique[0];U[1]=Ucinetique[1];
//Uc_fcc				//Uc_bcc
U[2]=Ucinetique[2];U[3]=Ucinetique[3];
//muc			//mux			//gtr
U[4]=muc/(R*T);U[5]=mucr/(R*T);U[6]=gtr/(R*T);
//v						//r			//Delta_ni
U[7]=Ycinetique_para[4];U[8]=Rzin;U[9]=Sursaturation_ni(Ucinetique,Unom);
//Delta_c							//time
U[10]=Sursaturation_c(Ucinetique,Unom);U[11]=tempsFeCCr;	
//U0ni			//U0c			//Dni/Dni
U[12]=Unom[1];U[13]=Unom[0];U[14]=coeff_k;	
so5.ecritureX(Xcinetique);
so4.ecritureU(U);
}	
else{return 0;}
}
