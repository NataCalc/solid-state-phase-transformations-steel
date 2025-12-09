/*
 *  2dpreci.cpp
 *  
 *
 *  Created by nat on 23/03/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "donnees.h"
#include "phase.hpp"
#include "newt.h"
#include "simplex.h"
#include "thermo.h"
#include "qhull.h"
#include "integer.h"
#include <cmath>
#include <time.h>
#include <fstream>


extern double temps,dt;
extern double T,dT;
extern float taille_gradient;
extern double dtemps;
extern float taille_gradient_alpha;
extern float taille_gradient_betta;
extern segment seg[20];
extern int nbseg;
extern int def;
extern double compteur_c;
extern double M;
extern double Vm;
extern double diffFeC;
extern double x_c[10];
extern double x_c0;
extern double x_c01;
extern double fricFeC;
extern double VIT_TEMPORAIRE;
extern double y_c_Svoboda;
extern double gtr;
extern double gtr1;
extern double muc;
extern double dcr,Rzin,correc_D;
extern double y_fraction_x_bcc,y_fraction_x_fcc;


double tempsFeC, dtempsFeC;
extern sortie so;
extern sortie so1;
extern sortie so2;
extern sortie so3;
extern sortie so4;
extern sortie so5;
extern sortie so6;
extern sortie so7;
extern sortie so8;
extern sortie so9;

extern matrice mat;

extern distribution_moyenne<sphere> cem;
extern distribution_moyenne<sphere> eps;


extern double correc_D;
extern float CIN;

extern int calc_spaces3(char *str); 


void calc2dequilibreFeNi();
bool calc2cinetique1FeNi(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc2cinetique_FeNi(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
						 std::vector<double> &prevVectorcin,bool &flage);
bool calc2cinetique_V(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
					  std::vector<double> &prevVectorcin,bool &flage);
bool calc2cinetique_FeNi_growth(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc2cinetique_FeNi_growth_growth(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
									   std::vector<double> &prevVectorcin,bool &flage);


extern clock_t t_debut_main , t_debut_qhull , t_fin_qhull , t_debut_NR , t_fin_NR , t_fin_main ;


void calc2dequilibreFeNi()
{
def=5;
int dim = 2;
newt nt;
//qhull qh(dim);
	
double prevT=0;
bool conv = false;
bool flage= false;
	
vector<double> prevVector(6);
std::vector<double> prevVectorcin(7);	
for(int i=0;i<7;i++){prevVectorcin[i]=0;}	
for(int i=0;i<6;i++){prevVector[i]=0;}

	
std::vector <double> CONCENTRATION_NOMINAL ( 2 );
CONCENTRATION_NOMINAL[0] = 6.e-2;// carbon -Ni
CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0];// iron -FE

	
printf ( "Begin calcul the system Fe-Ni equilibre thermodynamique (2D): \n" );
for ( int iterseg=1;iterseg<=nbseg;iterseg++ )
{
		
temps=seg[iterseg].date;
dt= ( seg[iterseg+1].date-temps ) /seg[iterseg+1].pas;
T=seg[iterseg].temperature;
dT= dt * ( seg[iterseg+1].temperature-seg[iterseg].temperature )
	/ ( seg[iterseg+1].date-seg[iterseg].date );
while ( temps<seg[iterseg+1].date )
{
temps+=dt;
T+=dT;
bool good;
bool prevExts=false;
printf("T = %e \n",T);
	
int ki = 0;

{

	
//############# Parameters
Rzin=1.e-3;
dcr = 3.5e-5 * exp( -2.8600e5 /R/ T )*pow(1e6,2.); //ni-fcc
dtempsFeC=0.1;
tempsFeC=0.0;
taille_gradient=1.e-3;
M = 0.035*exp(-17700/T)*pow(1e6,4);
correc_D=1.e0;
Vm=1e13;	

/*if(ki!=0)
{	
double dTempsFeC = -1.;	
dtempsFeC=0.1;
printf("VIT = %e \n",VIT_TEMPORAIRE);	
double before=Rzin;
double Rzin=Rzin+VIT_TEMPORAIRE*dtempsFeC;
double after=Rzin; //printf("after=%e \n",after);
while((after-before)/before>1e-1)
{
dtempsFeC=dtempsFeC/10.;
Rzin=before+VIT_TEMPORAIRE*dtempsFeC;
after=Rzin;
}
printf("R = %e \n",Rzin);	
T = 1000 + dTempsFeC* dtempsFeC;
printf("T = %e \n",T);	
}*/
	
//Rzin = 1.549318e-03;
//dtempsFeC =  4.115952e+02;	
	
nt=newt();
int i = 1;
while ( mat.produit[i]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::SB_FeNi;
				
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JB_FeNi;
				
//qh.func1    =mat.produit[1]->frontiere[1]; 
//qh.nrfuncv1 =&interface::phase1_binaireFeNi;   
				
//qh.func2    =mat.produit[1]->frontiere[1]; 
//qh.nrfuncv2 =&interface::phase2_binaireFeNi;   
				i++;
}
if ( nt.func == NULL ) printf ( "newton.func = null\n" );
//if ( qh.func1== NULL ) printf ( "qhull.func = null\n" );  //
			
mat.Thermo->read (T);
cem.classe0->Thermo->read (T);
			
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 4;
std::vector <double> Y ( nt.nn );
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn ); 
int check = 0;

//Y(BCC_A2,FE)=0.96735535, Y(BCC_A2,NI)=3.2644648E-2, Y(BCC_A2,VA#2)=1,
//Y(FCC_A1,FE)=0.89639918, Y(FCC_A1,NI)=0.10360082, Y(FCC_A1,VA#2)=1
Y[1] = 1e-2/*2.5481696e-2*/;//ni
Y[0] = 1.-Y[1]; //fe
Y[3] = 1e-1/*6.7550725e-2*/;//ni
Y[2] = 1.-Y[3];//fe
//Y[1] = 1e-2;//ni
//Y[0] = 1.-Y[1]; //fe
//Y[3] = 0.1;//ni
//Y[2] = 1.-Y[3];//fe	


	
//printf ( "\n" );printf ( "Vector Y[i] before:" );
//for ( int k=0 ; k < 4 ; k++ ){printf ( "%g ",Y[k] );}printf ( "\n" );
good=nt.calc ( Y, &check,CONCENTRATION_NOMINAL,2,conv);

			
printf ( "\n" );printf ( "Vector Y[i] after:" );
for ( int k=0 ; k < 4; k++ ){printf ( "%g ",Y[k] );}printf ( "\n" );
//Convegance FRACTION_SITE dans FRACTION_MOLAIRE
std::vector <double> X(4);
double x_ni_fcc,x_ni_bcc,x_fe_fcc,x_fe_bcc;			
x_ni_bcc  = Y[1];
x_fe_bcc = Y[0];
x_ni_fcc  = Y[3];			
x_fe_fcc = Y[4];
X[0] = Y[0];
X[1] = Y[1];
X[2] = Y[2];
X[3] = Y[3];
	
if (good)
{
so2.ecritureX ( X );
}

printf("\n\n\n");	
printf("Kinetic model of Zener \n");
calc2cinetique1FeNi(CONCENTRATION_NOMINAL,Y,prevVectorcin,flage);
ki++;

}	
}
}
so1.fermeture();
so2.fermeture();
so3.fermeture();
so4.fermeture();
so5.fermeture();
so6.fermeture();
so7.fermeture();
so8.fermeture();
so9.fermeture();
so.fermeture();
}



bool calc2cinetique1FeNi(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage)
{
newt nt;
nt=newt();
	
bool good;
bool conv=false;
	
nt=newt();
int j = 1;
while ( mat.produit[j]!=NULL )
{
nt.func     =mat.produit[1]->frontiere[1];
nt.nrfuncv  =&interface::KTB_FeNi;
		
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKB_FeNi;
		
j++;
}
if ( nt.func == NULL ) printf ( "newton.func = null\n" );
if (nt.jacobian==NULL) printf ("jacobian.func=null\n");
	

float v_ini;
int diff_step;
FILE *fpp=fopen("coeff_cin.txt","r") ;
fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
fscanf(fpp , "%i  \n" , &diff_step) ;	
fclose(fpp);

//yni_bcc, yfe_bcc,yni_fcc,yfe_fcc,v	
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 5;
std::vector <double> Ycinetique ( nt.nn );
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn ); //Jacobian
int check=0;
	

Ycinetique[0] = Y[0]; //fe-bcc
Ycinetique[1] = Y[1]; //ni-bcc
Ycinetique[2] = Y[2];//fe-fcc
Ycinetique[3] = Y[3];//ni-fcc
Ycinetique[4] = 1.e-9;

	
printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
for ( int k=0 ; k < 5 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );
				
good=nt.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,2,conv);
				
printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
for ( int k=0 ; k < 5 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );		

if (conv)
{
std::vector <double> X(5);
X[0] = Ycinetique[0];
X[1] = Ycinetique[1];
X[2] = Ycinetique[2];
X[3] = Ycinetique[3];
X[4] = Ycinetique[4];
	
so1.ecriture ( Ycinetique );
so3.ecritureX ( X );
so6.ecritureU(X);
}	
	
printf("\n\n\n");	
printf("Kinetics of Odqvist\n");	
//calc2cinetique_FeNi(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin,flage);	
calc2cinetique_V(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin,flage);
//calc2cinetique_FeNi_growth(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin,flage);	
//calc2cinetique_FeNi_growth_growth(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin,flage);	

}




bool calc2cinetique_FeNi(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
						 std::vector<double> &prevVectorcin,bool &flage)
{	
	
newt nt;
nt = newt();
bool conv_para=false;
	
int f = 1;
while(mat.produit[f]!=NULL)
{
nt.func	= mat.produit[1]->frontiere[1];
nt.nrfuncv	= &interface::KTB_prof_FeNi;
	
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKB_prof_FeNi;
f++;		
}

//######  The parameters for Newton-Raphson ########
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 3;
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn );
std::vector <double> Ycinetique_para ( nt.nn );
	
Ycinetique_para[0] = Y[1]; //ni-bcc
Ycinetique_para[1] = Y[3];//ni-fcc
Ycinetique_para[2] = Y[4];	
	
printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );	

flage=nt.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,2,conv_para);

printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );		

if(conv_para)
{	
std::vector <double> X(5);
X[0] = 1.-Ycinetique_para[0];
X[1] = Ycinetique_para[0];
X[2] = 1.-Ycinetique_para[1];
X[3] = Ycinetique_para[1];
X[4] = Y[2];		
so5.ecritureX(X);
so4.ecritureU(X);
}	
}

bool calc2cinetique_V(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
					  std::vector<double> &prevVectorcin,bool &flage)
{
newt nt;
nt = newt();	
//fstream file_v;

int f = 1;
while(mat.produit[f]!=NULL)
{
nt.func	= mat.produit[1]->frontiere[1];
nt.nrfuncv	= &interface::KTB_prof_FeNi_V;
		
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKB_prof_FeNi_V;
f++;		
}

bool conv_para = false;
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 2;
std::vector <double> Ycinetique_veloc(nt.nn);
int compt = 0;	
nt.fvec.resize(nt.nn);
nt.fvec1.resize(nt.nn); 

	
newt nt1;
nt1=newt();
int w = 1;
while ( mat.produit[w]!=NULL )
{
nt1.func     =mat.produit[1]->frontiere[1];
nt1.nrfuncv  =&interface::KTB_profFeNi0;
nt1.jacobian =mat.produit[1]->frontiere[1];
nt1.jacobfunc=&interface::JKB_prof_FeNi0;				
w++;
}
int check1=0;
nt1.fvec.clear();
nt1.fvec1.clear();
nt1.nn = 1;
std::vector <double> Ycinetique_para1 ( nt1.nn );
nt1.fvec.resize ( nt1.nn );
nt1.fvec1.resize ( nt1.nn ); 
	
int j=0;
int compteur, diff_step;
double VIT_TEMPORAIRE0,raison,Delta_c,Delta_ni; 	
bool flage_growth_plus, flage_growth_minus;	
	
Ycinetique_veloc[0]=Y[1]; //C in ferrite
Ycinetique_veloc[1]=Y[3];//X in ferrite
VIT_TEMPORAIRE=Y[4];//Velocity fixe	
	
diff_step=1000;
raison = pow( pow (10. , 8) , 1. / (double) diff_step ) ;	
compteur = 0 ;	
VIT_TEMPORAIRE0=VIT_TEMPORAIRE;
flage_growth_plus=false;
flage_growth_minus=false;	
	
std::vector <double> Xcinetique(12);
int fi = 100000;	
std::vector <float> Omega(fi),DeltaZ(fi),Vv(fi),Xnibccx(fi),Xnifccx(fi),Vvv(fi);

while (compteur <= diff_step)
{
	
VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;
		
nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,2,conv_para);
if(!conv_para) {flage_growth_minus=false;break ;}
else flage_growth_minus=true;

if(conv_para)
{
y_fraction_x_bcc=Ycinetique_veloc[0];
y_fraction_x_fcc=Ycinetique_veloc[1];
bool conv_para1 = false;

Ycinetique_para1[0]=CONCENTRATION_NOMINAL[0];	
nt1.calc ( Ycinetique_para1, &check1,CONCENTRATION_NOMINAL,1,conv_para1);	
	
	
Xcinetique[0]=1.-Ycinetique_veloc[0];
Xcinetique[1]=Ycinetique_veloc[0];
Xcinetique[2]=1.-Ycinetique_veloc[1];
Xcinetique[3]=Ycinetique_veloc[1];
Xcinetique[4]=VIT_TEMPORAIRE;
Xcinetique[5]=gtr;
Xcinetique[6]=gtr1;	
Xcinetique[7]=muc;
Xcinetique[8]=VIT_TEMPORAIRE*(Vm/M)/(R*T);	
Xcinetique[9]=CONCENTRATION_NOMINAL[0];
Xcinetique[10]=correc_D;
Xcinetique[11]=Ycinetique_para1[0];	
so7.ecritureV ( Xcinetique,compt );
	

Omega[compt] = (Xcinetique[3]-CONCENTRATION_NOMINAL[0])/(Xcinetique[3]-Xcinetique[1]);
DeltaZ[compt]= 2*Rzin*(1.-Omega[compt])/Omega[compt];
Vv[compt]=dcr*Omega[compt]/DeltaZ[compt];
Xnibccx[compt]=Xcinetique[1];
Xnifccx[compt]=Xcinetique[3];	
Vvv[compt]=VIT_TEMPORAIRE;	
//printf("%d %e %e %e \n",compt,Omega[compt],DeltaZ[compt],Vv[compt]);		
	
}
compteur++;
compt++;
}
	
	
printf("\nBoucle dans l'autre sens : augmentation de la vitesse\n") ;
Ycinetique_veloc[0]=Y[1]; //C in ferrite
Ycinetique_veloc[1]=Y[3];//X in ferrite
VIT_TEMPORAIRE=Y[4];//Velocity fixe	
compteur = 0 ;
while (compteur >= -diff_step)
{		

VIT_TEMPORAIRE = VIT_TEMPORAIRE0 / pow(raison,(double) compteur) ;
nt.calc ( Ycinetique_veloc, &check,CONCENTRATION_NOMINAL,2,conv_para);
if(!conv_para){ flage_growth_plus=false;break ;}
if(VIT_TEMPORAIRE>200.){break;}
else flage_growth_plus=true;
	
if(conv_para)
{
y_fraction_x_bcc=Ycinetique_veloc[0];
y_fraction_x_fcc=Ycinetique_veloc[1];
bool conv_para1 = false;

Ycinetique_para1[0]=CONCENTRATION_NOMINAL[0];	
nt1.calc ( Ycinetique_para1, &check1,CONCENTRATION_NOMINAL,1,conv_para1);	
	
Xcinetique[0]=1.-Ycinetique_veloc[0];
Xcinetique[1]=Ycinetique_veloc[0];
Xcinetique[2]=1.-Ycinetique_veloc[1];
Xcinetique[3]=Ycinetique_veloc[1];
Xcinetique[4]=VIT_TEMPORAIRE;
Xcinetique[5]=gtr;
Xcinetique[6]=gtr1;	
Xcinetique[7]=muc;
Xcinetique[8]=VIT_TEMPORAIRE*(Vm/M)/(R*T);	
Xcinetique[9]=CONCENTRATION_NOMINAL[0];
Xcinetique[10]=correc_D;
Xcinetique[11]=Ycinetique_para1[0];
	
so7.ecritureV ( Xcinetique,compt );

/*Omega[compt] = (Xcinetique[3]-CONCENTRATION_NOMINAL[0])/(Xcinetique[3]-Xcinetique[1]);
DeltaZ[compt]= 2*Rzin*(1.-Omega[compt])/Omega[compt];
Vv[compt]=dcr*Omega[compt]/DeltaZ[compt];	
Xnibccx[compt]=Xcinetique[1];
Xnifccx[compt]=Xcinetique[3];
Vvv[compt]=VIT_TEMPORAIRE;*/	
	
}
compteur--;
compt++;
}
	
so7.fermeture();

/*Omega.resize(compt);DeltaZ.resize(compt);Vv.resize(compt);
Xnifccx.resize(compt);Xnibccx.resize(compt);Vvv.resize(compt);
fstream file_vv;
char filename[30];
sprintf ( filename,"V2.txt");
file_vv.open ( filename,fstream::out );
if ( file_vv.fail() )
	cout<<"Error opening file"<<endl;

for(int zz=0;zz<compt;zz++)
{
file_vv<<Xnibccx[zz]<<" "<<Xnifccx[zz]<<" "<<Omega[zz]<<" "<<Vv[zz]<<" "<<Vvv[zz]<<" "<<DeltaZ[zz]<<" "<<zz<<endl;
	
}*/		

//file_vv.close();

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
xUcbcc(nx),xUcfcc(nx),xvx(nx),xgtrx(nx),xgtr1x(nx), \
xmucx(nx),xfricx(nx),xnomx(nx),xcorx(nx),xnomefx(nx);	

float Ts,Uxfcc,Uxbcc, \
Ucbcc,Ucfcc,vx,gtrx,gtr1x, \
mucx,fricx,nomx,corx,nomefx;
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
//"#T  U(Ni,fcc) U(Ni,bcc) U(C,fcc) U(C,bcc) muc mucx gtr v Fv r Delta_x time Delta_c U0x U0c Dx/Dx"<<endl;	
sscanf(str, "%e %e %e %e %e %e %e %e %e %e %e %e %d %e \n",&Ts,&Ucbcc,&Uxbcc, \
	   &Ucfcc,&Uxfcc,&vx,&gtrx,&gtr1x, \
	   &mucx,&fricx,&nomx,&corx,&compts,&nomefx);
//printf("%d \n",compts);	
//printf("%e %e %e %e %e %e %e %e %e %e %e %e %d %e \n",Ts,Ucbcc,Uxbcc, \
		Ucfcc,Uxfcc,vx,gtrx,gtr1x,mucx,fricx,nomx,corx,compts,nomefx);
			
xTs[compts]=Ts;
xUxfcc[compts]=Uxfcc;
xUxbcc[compts]=Uxbcc;
xUcbcc[compts]=Ucbcc;
xUcfcc[compts]=Ucfcc;
xvx[compts]=vx;
xgtrx[compts]=gtrx;
xgtr1x[compts]=gtr1x;
xmucx[compts]=mucx;
xfricx[compts]=fricx;
xnomx[compts]=nomx;
xcorx[compts]=corx;
xnomefx[compts]=nomefx;	
wasNumbers = true;
}
}
int fik = compts+1;
xTs.resize(fik);xUxfcc.resize(fik);xUxbcc.resize(fik); \
xUcbcc.resize(fik);xUcfcc.resize(fik);xvx.resize(fik); \
xgtrx.resize(fik);xgtr1x.resize(fik); \
xmucx.resize(fik);xfricx.resize(fik); \
xnomx.resize(fik);xcorx.resize(fik);xnomefx.resize(fik);
fclose(fee);
	
float txnom;
for(int count=0;count<fik;count++)
{
	if((xUxfcc[count]-xUxbcc[count])<=1e-4)
		txnom = xnomefx[count];
	/*printf("%e %e %e \n",xUxfcc[count],xUxbcc[count],xnomefx[count]);*/
}	
printf("NOM = %e \n",txnom);	
std::vector <double> tTm(1);
tTm[0] = txnom;
so9.ecriture(tTm);
	
//for(int zz=0;zz<compt;zz++)
//			printf("%e %e %e \n",Omega[zz],DeltaZ[zz],Vv[zz]);
	
}





bool calc2cinetique_FeNi_growth(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
						 std::vector<double> &prevVectorcin,bool &flage)
{	

	
newt nt;
nt = newt();
	
int f = 1;
while(mat.produit[f]!=NULL)
{
nt.func	= mat.produit[1]->frontiere[1];
nt.nrfuncv	= &interface::KTB_prof_FeNi;
		
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKB_prof_FeNi;
f++;		
}
	
	//######  The parameters for Newton-Raphson ########
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 2;
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn );
std::vector <double> Ycinetique_para ( nt.nn );
//std::vector <double> prevVectorcin_para(2);


	
int diff_step = 1000;	
double raison = pow ( 0.1 , 1. / ( (double) diff_step ) ) ;
int compteur = 0;
flage=false;	
std::vector <double> prevVectorcin_para(2);
	
	
//while (tempsFeC<=25)	
{
if(!flage)
{	
Ycinetique_para[0] = Y[1]; //ni-bcc
Ycinetique_para[1] = Y[3];//ni-fcc
VIT_TEMPORAIRE = Y[4];	
}
/*else
{
Ycinetique_para[0] = a;//Y[1]; //ni-bcc
Ycinetique_para[1] = b;//Y[3];//ni-fccdouble incr ;
double incr = 1.e-1 ;
dtempsFeC = Rzin * incr / VIT_TEMPORAIRE ;
printf("dtempsFeC = %e \n",dtempsFeC);	
double before=Rzin;		
Rzin=before+VIT_TEMPORAIRE*dtempsFeC;
	
tempsFeC += dtempsFeC ;
double dTempsFeC = -1.;	
T = 1000 + dTempsFeC* dtempsFeC;	
mat.Thermo->read (T);
cem.classe0->Thermo->read (T);	
}*/

else {
Ycinetique_para[0] = prevVectorcin_para[0]; //ni-bcc
Ycinetique_para[1] = prevVectorcin_para[1];//ni-fcc

	
dtempsFeC=0.1;
double before=Rzin;
Rzin=Rzin+VIT_TEMPORAIRE*dtempsFeC;
double after=Rzin; //printf("after=%e \n",after);
printf("rez=%e \n",(after-before)/before);
while((after-before)/before>1e-3)
{
dtempsFeC=dtempsFeC/10.;
Rzin=before+VIT_TEMPORAIRE*dtempsFeC;
printf("dt=%e R=%e \n",dtempsFeC,Rzin);
after=Rzin;
printf("after-before/before=%e \n",(after-before)/after);
}
printf("dt= %e \n",dtempsFeC);	
tempsFeC=tempsFeC+dtempsFeC;			
}
	
printf("temps = %e Rzin = %e T = %e \n",tempsFeC,Rzin,T);
dcr = 3.5e-5 * exp( -2.8600e5 /R/ T )*pow(1e6,2.); //ni-fcc


bool conv_para = false;

printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );	

flage=nt.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);

printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );		

double Omega=(Ycinetique_para[1]-CONCENTRATION_NOMINAL[0])/(Ycinetique_para[1]-Ycinetique_para[0]);
double delta = 2.*(1.-Omega)/Omega;
if(delta<1.e-16)delta=1.e-15;
VIT_TEMPORAIRE=(dcr/(Rzin*delta))*(Ycinetique_para[1]-CONCENTRATION_NOMINAL[0])/(Ycinetique_para[1]-Ycinetique_para[0]);		
printf("VIT = %e \n",VIT_TEMPORAIRE);	
	
if(conv_para)
{		
std::vector <double> X(7);
X[0] = 1.-Ycinetique_para[0];
X[1] = Ycinetique_para[0]; // ni bcc
X[2] = 1.-Ycinetique_para[1];
X[3] = Ycinetique_para[1]; // ni fcc
X[4] = VIT_TEMPORAIRE;	
X[5] = Rzin;
X[6] = tempsFeC;
	
prevVectorcin_para[0]=Ycinetique_para[0];
prevVectorcin_para[1]=Ycinetique_para[1];
	
so5.ecritureX(X);
so4.ecritureU(X);
flage = true;
}
else{flage=false;}
compteur ++ ;
}
}



bool calc2cinetique_FeNi_growth_growth(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
								std::vector<double> &prevVectorcin,bool &flage)
{	
newt nt;
nt = newt();
	
int f = 1;
while(mat.produit[f]!=NULL)
{
nt.func	= mat.produit[1]->frontiere[1];
nt.nrfuncv	= &interface::KTB_prof_FeNigrowth;
		
nt.jacobian =mat.produit[1]->frontiere[1];
nt.jacobfunc=&interface::JKB_prof_FeNigrowth;
f++;		
}
	
//######  The parameters for Newton-Raphson ########
int check=0;
nt.fvec.clear();
nt.fvec1.clear();
nt.nn = 3;
nt.fvec.resize ( nt.nn );
nt.fvec1.resize ( nt.nn );
std::vector <double> Ycinetique_para ( nt.nn );
std::vector <double> prevVectorcin_para(nt.nn);
	
int diff_step = 1000;	
double raison = pow ( 0.1 , 1. / ( (double) diff_step ) ) ;
int compteur = 0;
flage=false;	
//Rzin = 12;	
//while (tempsFeC<=100)	
{
if(!flage)
{	
Ycinetique_para[0] = Y[1]; //ni-bcc
Ycinetique_para[1] = Y[3];//ni-fcc
Ycinetique_para[2] = Y[4];	
}
else {
Ycinetique_para[0] = prevVectorcin_para[0]; //ni-bcc
Ycinetique_para[1] = prevVectorcin_para[1];//ni-fcc
Ycinetique_para[2] = prevVectorcin_para[2];//ni-fcc

dtempsFeC=0.05;
Rzin=Rzin+Ycinetique_para[2]*dtempsFeC;
tempsFeC=tempsFeC+dtempsFeC;			

/*double before=Rzin;
Rzin=Rzin+Ycinetique_para[2]*dtempsFeC;
double after=Rzin; //printf("after=%e \n",after);
printf("rez=%e \n",(after-before)/before);
while((after-before)/before>1e1)
{
dtempsFeC=dtempsFeC/10.;
Rzin=before+Ycinetique_para[2]*dtempsFeC;
printf("dt=%e R=%e \n",dtempsFeC,Rzin);
after=Rzin;
printf("after-before/before=%e \n",(after-before)/after);
}
printf("dt= %e \n",dtempsFeC);	
tempsFeC=tempsFeC+dtempsFeC;*/		
}
		
printf("temps = %e Rzin = %e T = %e \n",tempsFeC,Rzin,T);
dcr = 3.5e-5 * exp( -2.8600e5 /R/ T )*pow(1e6,2.); //ni-fcc
bool conv_para = false;

printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );	
		
flage=nt.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
		
printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );		

		
if(conv_para)
{		
std::vector <double> X(7);
X[0] = 1.-Ycinetique_para[0];
X[1] = Ycinetique_para[0]; // ni bcc
X[2] = 1.-Ycinetique_para[1];
X[3] = Ycinetique_para[1]; // ni fcc
X[4] = Ycinetique_para[2];	
X[5] = Rzin;
X[6] = tempsFeC;
			
prevVectorcin_para[0]=Ycinetique_para[0];
prevVectorcin_para[1]=Ycinetique_para[1];
prevVectorcin_para[2]=Ycinetique_para[2];
	
			
so5.ecritureX(X);
so4.ecritureU(X);
flage = true;
}
else{flage=false;}
compteur ++ ;
}
}







