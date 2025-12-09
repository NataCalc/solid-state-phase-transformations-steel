 /*
cinetic.cpp
 */

#include "phase.hpp"
#include "donnees.h"
#include "integration.h"
#include <fstream>
#include <algorithm>

extern float Diff[10];
extern float taille_gradient;
extern double correc_D;
extern float CIN;
extern double compteur_c;
extern double compteur_cr;
extern double Pmob;
extern double x_c[10];
extern double M;
extern double diffFeC;
extern double fricFeC;
extern double x_cFeCrC[11];
extern double x_crFeCrC[11];
extern double x_c0;
extern double x_c01;
extern double x_cr0;
extern double diffFeCrC;
extern double fricFeCrC;
extern double VIT_TEMPORAIRE;
extern double Rzin;
extern double y_c_Svoboda;
double gtr;
double gtr1;
double muc;
double mucr;
extern double Vm;
extern double dc , dcr;
extern double coeff_n;
extern double coeff_k;
extern double Lcoef,Bcoef;
extern int le_para;
extern double y_fraction_c_bcc;
extern double y_fraction_x_bcc; 
extern double y_fraction_x_fcc;
extern double y_fraction_c_fcc;

//extern double taille_temporaire;

std::vector <double> interface::derivFeC(double y)
{
	vector <double> Y(8);

	//1 order
	Y[0]= -(ph1->Thermo->m[0]/ph1->Thermo->m[1])* \
			((1.+4.*y)/pow((1.+2*y),2));//C dans ferrite
	
	Y[1]= -(ph1->Thermo->m[0]/ph1->Thermo->m[1])* \
			((4.*y+5.)/pow((1.+y),2));//VA dans ferrite
	
	Y[2]= -(ph2->Thermo->m[0]/ph2->Thermo->m[1])* \
			((1.+4.*y)/pow((1.+2*y),2));	//C dans austenite
	
	Y[3]= -(ph2->Thermo->m[0]/ph2->Thermo->m[1])* \
			(1./pow((1.+y),2)); //VA dans austenite
	
	//2 order
	Y[4]=  (ph1->Thermo->m[0]/ph1->Thermo->m[1])* \
			(8.*(1.+3.*y)/pow((1.+2*y),2));	//C dans ferrite
	
	Y[5]= (ph1->Thermo->m[0]/ph1->Thermo->m[1])* \
			((12.*y+14.)/pow((y+1),3)); //VA dans ferrite
	
	
	Y[6]= (ph2->Thermo->m[0]/ph2->Thermo->m[1])* \
			(8.*(1.+3.*y)/pow((1.+2*y),2));		//C dans austenite
	
	Y[7]= (ph2->Thermo->m[0]/ph2->Thermo->m[1])* \
			(2./pow((1.+y),3)); //VA dans austenite
	
	
	return Y;
	
}

std::vector <double> interface::derivFeC_U(double y)
{
	vector <double> Y(8);
	
	//1 order
	Y[0]=ph1->Thermo->m[0]/ph1->Thermo->m[1];//C dans ferrite
	
	Y[1]= -(ph1->Thermo->m[0]/ph1->Thermo->m[1]);	//VA dans ferrite
	
	Y[2]= ph2->Thermo->m[0]/ph2->Thermo->m[1];//C dans austenite
	
	Y[3]= -(ph2->Thermo->m[0]/ph2->Thermo->m[1]);	//VA dans austenite
	
	//2 order
	Y[4]= 0.;	//C dans ferrite
	
	Y[5]= 0.;	//VA dans ferrite
	
	
	Y[6]= 0.;		//C dans austenite
	
	Y[7]= 0.;		 //VA dans austenite
	
	
	return Y;
	
}




std::vector < std::vector<double> > interface::derivFeCrC(double yc, double ycr)
{
	int i;
	std::vector < std::vector<double> > P(12);
	for(i=0;i<12;i++)
		{
			P[i].resize(4);
		}


	P[0][0]=-ycr/pow((1.-yc),2.); /*dyfe/dcc (ferrite)*/	P[0][1]=-1./(1.-yc);/*dyfe/dccr (ferrite)*/ 
	P[0][2]=-ycr/pow((1.-yc),2.);/*dyfe/dcc (austenite)*/	P[0][3]=-1./(1.-yc);/*dyfe/dccr (austenite)*/
	
	P[1][0]=ycr/pow((1.-yc),2.);/*dycr/dcc (ferrite)*/		P[1][1]=1./(1.-yc);/*dycr/dccr (ferrite)*/
	P[1][2]=ycr/pow((1.-yc),2.);/*dycr/dcc (austenite)*/	P[1][3]=1./(1.-yc);/*dycr/dccr (austenite)*/
	
	P[2][0]=ph1->Thermo->m[0]/(ph1->Thermo->m[1]*pow((1.-yc),2.));/*dyc/dcc (ferrite)*/		P[2][1]=0.0;/*dyc/dccr (ferrite)*/
	P[2][2]=ph2->Thermo->m[0]/(ph2->Thermo->m[1]*pow((1.-yc),2.));/*dyc/dcc (austenite)*/	P[2][3]=0.0;/*dyc/dccr (austenite)*/

	P[3][0]=-ph1->Thermo->m[0]/(ph1->Thermo->m[1]*pow((1.-yc),2.));/*dyva/dcc (ferrite)*/	P[3][1]=0.0;/*dyva/dccr (ferrite)*/
	P[3][2]=-ph2->Thermo->m[0]/(ph2->Thermo->m[1]*pow((1.-yc),2.));/*dyva/dcc (austenite)*/	P[3][3]=0.0;/*dyva/dccr (austenite)*/
	
	
	P[4][0]=-2.*ycr/pow((1.-yc),3.);/*d2yfe/d2cc (ferrite)*/	P[4][1]=0.0;/*d2yfe/d2ccr (ferrite)*/
	P[4][2]=-2.*ycr/pow((1.-yc),3.);/*d2yfe/d2cc (austenite)*/	P[4][3]=0.0;/*d2yfe/d2ccr (austenite)*/
	
	P[5][0]=2.*ycr/pow((1.-yc),3.);/*d2ycr/d2cc (ferrite)*/		P[5][1]=0.0;/*d2ycr/d2ccr (ferrite)*/
	P[5][2]=2.*ycr/pow((1.-yc),3.);/*d2ycr/d2cc (austenite)*/	P[5][3]=0.0;/*d2ycr/d2ccr (austenite)*/
	
	P[6][0]=2.*ph1->Thermo->m[0]/(ph1->Thermo->m[1]*pow((1.-yc),3.));/*d2yc/d2cc (ferrite)*/	P[6][1]=0.0;/*d2yc/d2yccr (ferrite)*/
	P[6][2]=2.*ph2->Thermo->m[0]/(ph2->Thermo->m[1]*pow((1.-yc),3.));/*d2yc/d2cc (austenite)*/	P[6][3]=0.0;/*d2yc/d2yccr (austenite)*/
	
	P[7][0]=-2.*ph1->Thermo->m[0]/(ph1->Thermo->m[1]*pow((1.-yc),3.));/*d2yva/d2cc (ferrite)*/	P[7][1]=0.0;/*d2yva/d2yccr (ferrite)*/
	P[7][2]=-2.*ph2->Thermo->m[0]/(ph2->Thermo->m[1]*pow((1.-yc),3.));/*d2yva/d2cc (austenite)*/P[7][3]=0.0;/*d2yva/d2yccr (austenite)*/
	
	P[8][0]=-1./pow((1.-yc),2.);/*d(dyfe/dcc)/dccr (ferrite)*/		P[8][1]=-1./pow((1.-yc),2.);/*d(dyfe/dccr)/dcc (ferrite)*/
	P[8][2]=-1./pow((1.-yc),2.);/*d(dyfe/dcc)/dccr (austenite)*/	P[8][3]=-1./pow((1.-yc),2.);/*d(dyfe/dccr)/dcc (austenite)*/
	
	P[9][0]=1./pow((1.-yc),2.);/*d(dycr/dcc)/dccr (ferrite)*/		P[9][1]=1./pow((1.-yc),2.);/*d(dycr/dccr)/dcc (ferrite)*/
	P[9][2]=1./pow((1.-yc),2.);/*d(dycr/dcc)/dccr (austenite)*/		P[9][3]=1./pow((1.-yc),2.);/*d(dycr/dccr)/dcc (austenite)*/
	
	
	P[10][0]=0.0;/*d(dyc/dcc)/dccr (ferrite)*/		P[10][1]=0.0;/*d(dyc/dccr)/dcc (ferrite)*/
	P[10][2]=0.0;/*d(dyc/dcc)/dccr (austenite)*/	P[10][3]=0.0;/*d(dyc/dccr)/dcc (austenite)*/
	
	return P;
	
}



std::vector < std::vector<double> > interface::derivFeCrC_U(double yc, double ycr)
{
	int i;
	std::vector < std::vector<double> > P(12);
	for(i=0;i<12;i++)
	{
		P[i].resize(4);
	}
	
	
	P[0][0]=0.0; /*dyfe/dcc (ferrite)*/	P[0][1]=-1.0;/*dyfe/dccr (ferrite)*/ 
	P[0][2]=0.0;/*dyfe/dcc (austenite)*/	P[0][3]=-1.0;/*dyfe/dccr (austenite)*/
	
	P[1][0]=0.0;/*dycr/dcc (ferrite)*/		P[1][1]=1.0;/*dycr/dccr (ferrite)*/
	P[1][2]=0.0;/*dycr/dcc (austenite)*/	P[1][3]=1.0;/*dycr/dccr (austenite)*/
	
	P[2][0]=ph1->Thermo->m[0]/ph1->Thermo->m[1];/*dyc/dcc (ferrite)*/		P[2][1]=0.0;/*dyc/dccr (ferrite)*/
	P[2][2]=ph2->Thermo->m[0]/ph2->Thermo->m[1];/*dyc/dcc (austenite)*/	P[2][3]=0.0;/*dyc/dccr (austenite)*/
	
	P[3][0]=-ph1->Thermo->m[0]/ph1->Thermo->m[1];/*dyva/dcc (ferrite)*/	P[3][1]=0.0;/*dyva/dccr (ferrite)*/
	P[3][2]=-ph2->Thermo->m[0]/ph2->Thermo->m[1];/*dyva/dcc (austenite)*/	P[3][3]=0.0;/*dyva/dccr (austenite)*/
	
	
	P[4][0]=0.0;/*d2yfe/d2cc (ferrite)*/	P[4][1]=0.0;/*d2yfe/d2ccr (ferrite)*/
	P[4][2]=0.0;/*d2yfe/d2cc (austenite)*/	P[4][3]=0.0;/*d2yfe/d2ccr (austenite)*/
	
	P[5][0]=0.0;/*d2ycr/d2cc (ferrite)*/		P[5][1]=0.0;/*d2ycr/d2ccr (ferrite)*/
	P[5][2]=0.0;/*d2ycr/d2cc (austenite)*/	P[5][3]=0.0;/*d2ycr/d2ccr (austenite)*/
	
	P[6][0]=0.0;/*d2yc/d2cc (ferrite)*/	P[6][1]=0.0;/*d2yc/d2yccr (ferrite)*/
	P[6][2]=0.0;/*d2yc/d2cc (austenite)*/	P[6][3]=0.0;/*d2yc/d2yccr (austenite)*/
	
	P[7][0]=0.0;/*d2yva/d2cc (ferrite)*/	P[7][1]=0.0;/*d2yva/d2yccr (ferrite)*/
	P[7][2]=0.0;/*d2yva/d2cc (austenite)*/P[7][3]=0.0;/*d2yva/d2yccr (austenite)*/
	
	P[8][0]=0.0;/*d(dyfe/dcc)/dccr (ferrite)*/		P[8][1]=0.0;/*d(dyfe/dccr)/dcc (ferrite)*/
	P[8][2]=0.0;/*d(dyfe/dcc)/dccr (austenite)*/	P[8][3]=0.0;/*d(dyfe/dccr)/dcc (austenite)*/
	
	P[9][0]=0.0;/*d(dycr/dcc)/dccr (ferrite)*/		P[9][1]=0.0;/*d(dycr/dccr)/dcc (ferrite)*/
	P[9][2]=0.0;/*d(dycr/dcc)/dccr (austenite)*/		P[9][3]=0.0;/*d(dycr/dccr)/dcc (austenite)*/
	
	
	P[10][0]=0.0;/*d(dyc/dcc)/dccr (ferrite)*/		P[10][1]=0.0;/*d(dyc/dccr)/dcc (ferrite)*/
	P[10][2]=0.0;/*d(dyc/dcc)/dccr (austenite)*/	P[10][3]=0.0;/*d(dyc/dccr)/dcc (austenite)*/

	P[11][0]=0.0;/*d(dyc/dcc)/dccr (ferrite)*/		P[11][1]=0.0;/*d(dyc/dccr)/dcc (ferrite)*/
	P[11][2]=0.0;/*d(dyc/dcc)/dccr (austenite)*/	P[11][3]=0.0;/*d(dyc/dccr)/dcc (austenite)*/
	
	return P;
	
}

void interface::funconeFeC(double *x, double *y,vector <double> &concentration, \
						   vector <double> &Ysol, int m, int n, \
						   std::vector < std::vector<double> > &P)
{
	//Diffusion
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	
	int i,j;
	std::vector <double> Y1alfa(1); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(1);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	double G1,G2;
	double ytemp,ytempU,xtemp, h, dh, dGm2aC, dGm2aVA, dGm2bC,dGm2bVA;
	double Energy_gibbs, Der_energy_cons, potchem, potchem_fe;
	vector <double> Y(8),YU(8);
	
	for(i=0; i<m+1; i++)
	{
		Y1alfa[0]  =  1.; //bcc_A2: FE
		Y2alfa[0]  =  (ph1->Thermo->m[0]/ph1->Thermo->m[1])*y[i]; //bcc_A2: C
		Y2alfa[1]  =  1.-Y2alfa[0]; //bcc_A2: VA
		Y1betta[0] =  1.; //fcc_A1: FE
		Y2betta[0] =  (ph2->Thermo->m[0]/ph2->Thermo->m[1])*y[i]; //fcc_A1: C 
		Y2betta[1] =  1.-Y2betta[0]; //fcc_A1: VA
		
		//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
		G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
		*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
		G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
		*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite

		
		//Derivation 1 order and 2 order y
		ytemp=y[i]/(1.+y[i]);	//fraction molaire x
		Y=derivFeC(ytemp);
		
		ytempU=y[i];			//U-fraction
		YU=derivFeC_U(ytempU);
		
		xtemp=x[i];	//l'espace
		h=Wang(xtemp);
		dh=dWang(xtemp);
		//Derivation energy gibbs 1 order
		ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
		dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
		dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans austenite
		
		ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
		dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
		dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite

		
		Energy_gibbs	= 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0])*G2*h \
						+1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0])*G1*(1.-h);

		
		
		Der_energy_cons	= h*(-ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.)*YU[2]*G2		\
							 +1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0])*(dGm2bC*YU[2]+dGm2bVA*YU[3]))			\
		+(1.-h)*(-ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.)*YU[0]*G1	\
				 +1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*(dGm2aC*YU[0]+dGm2aVA*YU[1]));

		
		potchem			= Energy_gibbs+(1.-ytemp)*Der_energy_cons*(1./pow((ytemp-1.),2));
		potchem_fe		= Energy_gibbs-ytemp*Der_energy_cons*(1./pow((ytemp-1.),2));
		
		P[i][0]=potchem;
		P[i][1]=potchem_fe;
	}	
}


void interface::funcderivmuFeC(double *x, double *y,vector <double> &concentration, \
							   vector <double> &Ysol, int m, int n, \
							   std::vector < std::vector<double> > &P)
{
	//Diffusion
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	
	int i,j;
	std::vector <double> Y1alfa(1); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(1);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	double G1,G2, h, dh;
	double X[2];
	double dGm2aC,dGm2aVA;
	double dGm2bC,dGm2bVA;
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	//Parametres:
	double fac_gam, fac_alf, der_gam, der_alf, \
	part1, part2, \
	der_sec_gam, der_sec_alf, \
	part3, part4;
	double a,b,res,L;
	double potchem_c;
	double ytemp, xtemp;
	vector <double> Y(8),YU(8);
	
	//double VIT= Ysol[6];
	for(i=0; i<m+1; i++)
	{
		Y1alfa[0]  =  1.; //bcc_A2: FE
		Y2alfa[0]  =  (ph1->Thermo->m[0]/ph1->Thermo->m[1])*y[i];//bcc_A2: C
		Y2alfa[1]  =  1.-Y2alfa[0]; //bcc_A2: VA
		Y1betta[0] =  1.; //fcc_A1: FE
		Y2betta[0] = (ph2->Thermo->m[0]/ph2->Thermo->m[1])*y[i];//fcc_A1: C 
		Y2betta[1] =  1.-Y2betta[0]; //fcc_A1: VA
		
		
		//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
		G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
		*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
		G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
		*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite
		//Derivation 1 order and 2 order y
		
		ytemp=y[i];	
		YU=derivFeC_U(ytemp);
		
		xtemp=x[i];	
		h=Wang(xtemp);
		dh=dWang(xtemp);
		//Fraction molaire de C
		//X[0]=(ph1->Thermo->m[1]*Y2alfa[0]) \
		/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y2alfa[1]);//carbone dans ferrite
		X[0]=x_c0;
		
		X[1]=(ph2->Thermo->m[1]*Y2betta[0]) \
		/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y2betta[1]);//carbone dans austenite
		
		ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
		dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
		dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans austenite
		
		ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
		dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
		dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
		
		
		//Derivation energy gibbs 2 order
		ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
		ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
		ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
		ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
		ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
		
		ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
		ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
		ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
		ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
		ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
		
		
		
		
		fac_gam	=	1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]);
		fac_alf	=	1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]);
		
		der_gam	=	dGm2bC*YU[2]+dGm2bVA*YU[3];
		der_alf	=	dGm2aC*YU[0]+dGm2aVA*YU[1];
		
		part1	=	-ph2->Thermo->m[1]*pow(fac_gam,2.)*YU[2]*G2+fac_gam*der_gam;
		part2	=	-ph1->Thermo->m[1]*pow(fac_alf,2.)*YU[0]*G1+fac_alf*der_alf;
		
		
		der_sec_gam	=	ddGm2bCC*pow(YU[2],2.)+ddGm2bVAVA*pow(YU[3],2.)+dGm2bC*YU[6]+dGm2bVA*YU[7] \
						+ddGm2bCVA*YU[2]*YU[3]+ddGm2bVAC*YU[2]*YU[3];
		der_sec_alf	=	ddGm2aCC*pow(YU[0],2.)+ddGm2aVAVA*pow(YU[1],2.)+dGm2aC*YU[4]+dGm2aVA*YU[5] \
						+ddGm2aCVA*YU[0]*YU[1]+ddGm2aVAC*YU[0]*YU[1];
		
		
		part3	=	2.*pow(ph2->Thermo->m[1],2.)*pow(fac_gam,3.)*pow(YU[2],2.)*G2 \
					-ph2->Thermo->m[1]*pow(fac_gam,2)*YU[6]*G2 \
					-2.*ph2->Thermo->m[1]*pow(fac_gam,2.)*YU[2]*der_gam \
					+fac_gam*der_sec_gam;
		
		part4	=	2.*pow(ph1->Thermo->m[1],2.)*pow(fac_alf,3.)*pow(YU[0],2.)*G1 \
					-ph1->Thermo->m[1]*pow(fac_alf,2)*YU[4]*G1 \
					-2.*ph1->Thermo->m[1]*pow(fac_alf,2.)*YU[0]*der_alf \
					+fac_alf*der_sec_alf;
		
		
		L = correc_D/( (2.*(h*part1+(1.-h)*part2)+(y[i]+1.)*(h*part3+(1.-h)*part4))*Vm);//pow(taille_gradient,2);
		
		double vitess=VIT_TEMPORAIRE;//taille_gradient;
		
		a		= vitess*( y[i]-X[0])+L*Vm*dh*(1./taille_gradient)*(fac_gam*G2-fac_alf*G1+(y[i]+1.)*(part1-part2));
		b		= -L*Vm*( (y[i]+1)*(h*part3+(1.-h)*part4)+2.*(h*part1+(1.-h)*part2));
		res		= taille_gradient*a/b;		
		
		
		double potchem_c_inter =taille_gradient* \
			(dh*(1./taille_gradient)*(fac_gam*G2-fac_alf*G1+(y[i]+1.)*(part1-part2))+ \
			res*(1./taille_gradient)*((y[i]+1.)*(h*part3+(1.-h)*part4)+2.*(h*part1+(1.-h)*part2)));
		

		
		potchem_c = pow(potchem_c_inter,2);

		
		P[i][0]=(1./taille_gradient)*Vm*L*potchem_c;
		
		P[i][1]=L;
		P[i][2]=(1.-ytemp)*(h*part3+(1.-h)*part4);
		//printf("P = %e \n",P[i][0]);
	}	
}

void interface::funconeFeCrC(int m, double xi, double xf,double *yc,double *ycr,std::vector < std::vector<double> > &P)
{
	double h1,x;
	
	h1=(xf-xi)/(m-1); 
	x=xi-h1;
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);//bcc_A2: Fe,Cr
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: Fe,Cr
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	double G1,G2, yctemp, ycrtemp, xtemp, h, dh;
	double dGm2aC,dGm2aVA,dGm2aFe,dGm2aCr,dGm2bC,dGm2bVA,dGm2bFe,dGm2bCr;
	double Energy_gibbs, potchem_fe, potchem_cr, potchem_c;
	double ddGm[2], dGma[2], dGmb[2];
	int i,j;
	std::vector < std::vector<double> > D(12);
	for(i=0;i<12;i++)
	{
		D[i].resize(4);
	}
	
	for(i=1; i<m+1; i++)
	{
		x += h1;
		//cout<<"x"<<x<<endl;
		
		Y1alfa[1]  = ycr[i];
		Y1alfa[0]  = 1.-Y1alfa[1]; //FE
		Y2alfa[0]  = yc[i]/ph1->Thermo->m[1]; //C
		Y2alfa[1]  = 1.-Y2alfa[0]; //VA
		Y1betta[1] = ycr[i];
		Y1betta[0] = 1.-Y1betta[1]; //FE
		Y2betta[0] = yc[i]/ph2->Thermo->m[1]; //C
		Y2betta[1] = 1.-Y2betta[0]; //VA

		//cout<<Y1alfa[1]<<Y2alfa[0]<<endl;
		
		//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
		G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
		*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
		//cout<<"G1="<<G1<<endl;
		G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
		*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite
		//cout<<"G2="<<G2<<endl;
		
		//Derivation 1 order and 2 order y
		yctemp=yc[i];	
		ycrtemp=ycr[i];
		D=derivFeCrC_U(yctemp,ycrtemp);
		
		//for(i=0;i<12;i++)
		//{
		//	for(j=0;j<4;j++)
		//	{printf("D[%i][%i]=%lf ",i,j,D[i][j]);}printf("\n");
		//}printf("\n");
		h=Wang(x);
		dh=dWang(x);
		
		//Derivation energy gibbs 1 order
		ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
		dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
		//cout<<"dGm2aC="<<dGm2aC<<endl;
		dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans ferrite
		//cout<<"dGm2aVA="<<dGm2aVA<<endl;
		dGm2aFe	=ph1->Thermo->dGm1[0];	//Fe dans ferrite
		//cout<<"dGm2aFe="<<dGm2aFe<<endl;
		dGm2aCr	=ph1->Thermo->dGm1[1];	//Cr dans ferrite
		//cout<<"dGm2aCr="<<dGm2aCr<<endl;
		ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
		dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
		//cout<<"dGm2bC="<<dGm2bC<<endl;
		dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
		//cout<<"dGm2bVA="<<dGm2bVA<<endl;
		dGm2bFe	=ph2->Thermo->dGm1[0];	//Fe dans austenite
		//cout<<"dGm2bFe="<<dGm2bFe<<endl;
		dGm2bCr	=ph2->Thermo->dGm1[1];	//Cr dans austenite
		//cout<<"dGm2bCr="<<dGm2bCr<<endl;
		
		Energy_gibbs	= 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0])*G2*h \
		+1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0])*G1*(1.-h);
		
		
		dGmb[0]=dGm2bC*D[2][2]+dGm2bVA*D[3][2]+dGm2bFe*D[0][2]+dGm2bCr*D[1][2]; //*dy/dcc (austenite)
		dGma[0]=dGm2aC*D[2][0]+dGm2aVA*D[3][0]+dGm2aFe*D[0][0]+dGm2aCr*D[1][0]; //*dy/dcc (ferrite)

		ddGm[0]=h*(-ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.)*D[2][2]*G2 \
				+1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0])*dGmb[0]) \
				+(1.-h)*(-ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.)*D[2][0]*G1 \
				+1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*dGma[0]);/*dGm/dcc*/

		dGmb[1]=dGm2bC*D[2][3]+dGm2bVA*D[3][3]+dGm2bFe*D[0][3]+dGm2bCr*D[1][3]; //*dy/dccr (austenite)
		dGma[1]=dGm2aC*D[2][1]+dGm2aVA*D[3][1]+dGm2aFe*D[0][1]+dGm2aCr*D[1][1]; //*dy/dccr (ferrite)

		ddGm[1]=h*(-ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.)*D[2][3]*G2 \
				+1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0])*dGmb[1]) \
				+(1.-h)*(-ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.)*D[2][1]*G1 \
				+1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*dGma[1]);/*dGm/dccr*/
		
		potchem_fe = Energy_gibbs-yc[i]*(yc[i]+1.)*ddGm[0]-ycr[i]*(yc[i]+1.)*ddGm[1];
		potchem_cr = Energy_gibbs-yc[i]*(yc[i]+1.)*ddGm[0]+(yc[i]-ycr[i]+1.-yc[i]*ycr[i])*ddGm[1];
		potchem_c  = Energy_gibbs+(yc[i]+1.)*ddGm[0];
		
		//printf("potchem_fe = %e \n",potchem_fe);
		//printf("potchem_cr = %e \n",potchem_cr);
		//printf("potchem_c  = %e \n",potchem_c);
		//printf("energie=%e \n",Energy_gibbs);
		
		P[i][0]=potchem_c;
		P[i][1]=potchem_fe;
		P[i][2]=potchem_cr;
		P[i][3]=potchem_cr-potchem_fe;
    }		
		
}


double interface::funcFeC(double x, double y,vector <double> &concentration, vector <double> &Ysol)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	
	int i,j;
	std::vector <double> Y1alfa(1); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(1);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	vector <double> Y(8),YU(8);
	double G1,G2, h, dh,xtemp;
	double X[2];
	double dGm2aC,dGm2aVA;
	double dGm2bC,dGm2bVA;
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	//Parametres:
	double fac_gam, fac_alf, der_gam, der_alf, \
	part1, part2, \
	der_sec_gam, der_sec_alf, \
	part3, part4;
	double a,b,res,L;
	
	h=Wang(x);
	dh=dWang(x);
	
	Y1alfa[0]  =  1.; //bcc_A2: FE
	Y2alfa[0]  =  (ph1->Thermo->m[0]/ph1->Thermo->m[1])*y;//bcc_A2: C
	Y2alfa[1]  =  1.-Y2alfa[0]; //bcc_A2: VA
	Y1betta[0] =  1.; //fcc_A1: FE
	Y2betta[0] = (ph2->Thermo->m[0]/ph2->Thermo->m[1])*y;//fcc_A1: C 
	Y2betta[1] =  1.-Y2betta[0]; //fcc_A1: VA
	
	
	//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
	G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
	*ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa);		//ferrite
	G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
	*ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta);	//austenite
	
	
	YU=derivFeC_U(y);
	
	xtemp=y/(1.+y);
	Y =derivFeC(xtemp);  	
	
	X[0]=x_c0;
	
	double duc = Vm;
	double dduc=0.0;
	
	
	X[1]=(ph2->Thermo->m[1]*Y2betta[0]) \
	/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y2betta[1]);//carbone dans austenite
	
	
	//Function Wang h(x)=x*x*(3-2*x)
	//x=x/taille_gradient;

	//Derivation function Wang dh/dx=

	
	
	//Derivation energy gibbs 1 order
	ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
	dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
	dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans austenite
	
	ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
	dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
	dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
	
	
	//Derivation energy gibbs 2 order
	ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
	ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
	ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
	ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
	ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
	
	ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
	ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
	ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
	ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
	ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
	
	
	fac_gam	=	1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]);
	fac_alf	=	1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]);
	
	der_gam	=	dGm2bC*YU[2]+dGm2bVA*YU[3];
	der_alf	=	dGm2aC*YU[0]+dGm2aVA*YU[1];

	
	
	part1	=	-ph2->Thermo->m[1]*pow(fac_gam,2.)*YU[2]*G2+fac_gam*der_gam;
	part2	=	-ph1->Thermo->m[1]*pow(fac_alf,2.)*YU[0]*G1+fac_alf*der_alf;
	
	
	der_sec_gam	=	ddGm2bCC*pow(YU[2],2.)+ddGm2bVAVA*pow(YU[3],2.)+dGm2bC*YU[6]+dGm2bVA*YU[7] \
					+ddGm2bCVA*YU[2]*YU[3]+ddGm2bVAC*YU[2]*YU[3];
	der_sec_alf	=	ddGm2aCC*pow(YU[0],2.)+ddGm2aVAVA*pow(YU[1],2.)+dGm2aC*YU[4]+dGm2aVA*YU[5] \
					+ddGm2aCVA*YU[0]*YU[1]+ddGm2aVAC*YU[0]*YU[1];
	
	
	part3	=	2.*pow(ph2->Thermo->m[1],2.)*pow(fac_gam,3.)*pow(YU[2],2.)*G2 \
					-ph2->Thermo->m[1]*pow(fac_gam,2)*YU[6]*G2 \
					-2.*ph2->Thermo->m[1]*pow(fac_gam,2.)*YU[2]*der_gam \
					+fac_gam*der_sec_gam;
	
	part4	=	2.*pow(ph1->Thermo->m[1],2.)*pow(fac_alf,3.)*pow(YU[0],2.)*G1 \
					-ph1->Thermo->m[1]*pow(fac_alf,2)*YU[4]*G1 \
					-2.*ph1->Thermo->m[1]*pow(fac_alf,2.)*YU[0]*der_alf \
					+fac_alf*der_sec_alf;

	
	
	L = correc_D/( (2.*(h*part1+(1.-h)*part2)+(y+1.)*(h*part3+(1.-h)*part4))*Vm);//pow(taille_gradient,2);
	//printf("L=%e \n",L*Vm);
	
	double vitess=VIT_TEMPORAIRE;//taille_gradient;
	
	a		= vitess*( y-X[0])+L*Vm*dh*(1./taille_gradient)*(fac_gam*G2-fac_alf*G1+(y+1.)*(part1-part2));
	b		= -L*Vm*( (y+1)*(h*part3+(1.-h)*part4)+2.*(h*part1+(1.-h)*part2));
	res		= taille_gradient*a/b;
	
	return res;
	
}




double interface::deriv_mu(double y,double x)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	
	int i,j;
	std::vector <double> Y1alfa(1); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(1);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	vector <double> Y(8),YU(8);
	double G1,G2, h, dh,xtemp;
	double X[2];
	double dGm2aC,dGm2aVA;
	double dGm2bC,dGm2bVA;
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	//Parametres:
	double fac_gam, fac_alf, der_gam, der_alf, \
	part1, part2, \
	der_sec_gam, der_sec_alf, \
	part3, part4;
	double a,b,res,L;
	
	h=Wang(x);
	dh=dWang(x);
	
	Y1alfa[0]  =  1.; //bcc_A2: FE
	Y2alfa[0]  =  (ph1->Thermo->m[0]/ph1->Thermo->m[1])*y;//bcc_A2: C
	Y2alfa[1]  =  1.-Y2alfa[0]; //bcc_A2: VA
	Y1betta[0] =  1.; //fcc_A1: FE
	Y2betta[0] = (ph2->Thermo->m[0]/ph2->Thermo->m[1])*y;//fcc_A1: C 
	Y2betta[1] =  1.-Y2betta[0]; //fcc_A1: VA
	
	
	//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
	G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
	*ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa);		//ferrite
	G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
	*ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta);	//austenite
	
	
	YU=derivFeC_U(y);
	
	xtemp=y/(1.+y);
	Y =derivFeC(xtemp);  	
	
	X[0]=x_c0;
	
	double duc = Vm;
	double dduc=0.0;
	
	
	X[1]=(ph2->Thermo->m[1]*Y2betta[0]) \
	/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y2betta[1]);//carbone dans austenite
	
	
	//Function Wang h(x)=x*x*(3-2*x)
	//x=x/taille_gradient;
	
	//Derivation function Wang dh/dx=
	
	
	
	//Derivation energy gibbs 1 order
	ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
	dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
	dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans austenite
	
	ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
	dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
	dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
	
	
	//Derivation energy gibbs 2 order
	ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
	ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
	ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
	ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
	ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
	
	ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
	ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
	ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
	ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
	ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
	
	
	fac_gam	=	1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]);
	fac_alf	=	1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]);
	
	der_gam	=	dGm2bC*YU[2]+dGm2bVA*YU[3];
	der_alf	=	dGm2aC*YU[0]+dGm2aVA*YU[1];
	
	
	
	part1	=	-ph2->Thermo->m[1]*pow(fac_gam,2.)*YU[2]*G2+fac_gam*der_gam;
	part2	=	-ph1->Thermo->m[1]*pow(fac_alf,2.)*YU[0]*G1+fac_alf*der_alf;
	
	
	der_sec_gam	=	ddGm2bCC*pow(YU[2],2.)+ddGm2bVAVA*pow(YU[3],2.)+dGm2bC*YU[6]+dGm2bVA*YU[7] \
	+ddGm2bCVA*YU[2]*YU[3]+ddGm2bVAC*YU[2]*YU[3];
	der_sec_alf	=	ddGm2aCC*pow(YU[0],2.)+ddGm2aVAVA*pow(YU[1],2.)+dGm2aC*YU[4]+dGm2aVA*YU[5] \
	+ddGm2aCVA*YU[0]*YU[1]+ddGm2aVAC*YU[0]*YU[1];
	
	
	part3	=	2.*pow(ph2->Thermo->m[1],2.)*pow(fac_gam,3.)*pow(YU[2],2.)*G2 \
	-ph2->Thermo->m[1]*pow(fac_gam,2)*YU[6]*G2 \
	-2.*ph2->Thermo->m[1]*pow(fac_gam,2.)*YU[2]*der_gam \
	+fac_gam*der_sec_gam;
	
	part4	=	2.*pow(ph1->Thermo->m[1],2.)*pow(fac_alf,3.)*pow(YU[0],2.)*G1 \
	-ph1->Thermo->m[1]*pow(fac_alf,2)*YU[4]*G1 \
	-2.*ph1->Thermo->m[1]*pow(fac_alf,2.)*YU[0]*der_alf \
	+fac_alf*der_sec_alf;
	
	
	
	L = correc_D/( (2.*(h*part1+(1.-h)*part2)+(y+1.)*(h*part3+(1.-h)*part4)));//pow(taille_gradient,2);
	//printf("L=%e \n",L);

	
	return L;
	
}


								//Ni-fcc //1 Ni-bcc	
double interface::deriv_mu_FeNi(double y,double x)
{
	
double TT2 = 0.;
TT2 = ph1->Thermo->T1;
	
int i,j;
std::vector <double> Y1alfa(2); //bcc_A2: FE
std::vector <double> Y2alfa(2); //bcc_A2: C,VA
std::vector <double> Y1betta(2);//fcc_A1: FE
std::vector <double> Y2betta(2);//fcc_A1: C,VA
vector <double> Y(8),YU(8);
double G1,G2, h, dh,xtemp;
double X[2];
double dGm2aC,dGm2aVA;
double dGm2bC,dGm2bVA;
double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
//Parametres:
double fac_gam, fac_alf, der_gam, der_alf, \
part1, part2, \
der_sec_gam, der_sec_alf, \
part3, part4;
double a,b,res,L;
	

	
//printf("%e %e\n",y,x);	
Y1alfa[0]  =  1.-y; //bcc_A2: FE
Y1alfa[1]  =  y;//bcc_A2: Ni
Y2alfa[0]  =  0.0;	
Y2alfa[1]  =  1.; //bcc_A2: VA
Y1betta[0] =  1.-y; //fcc_A1: FE
Y1betta[1] =  y;	
Y2betta[0] =  0.0;//fcc_A1: C 
Y2betta[1] =  1.; //fcc_A1: VA
	
//Derivation energy gibbs 1 order
ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
double dGm2aFe	=ph1->Thermo->dGm1[0];	//C dans ferrite
double dGm2aNi	=ph1->Thermo->dGm1[1];	//VA dans austenite
	
ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
double dGm2bFe	=ph2->Thermo->dGm1[0];	//C dans austenite
double dGm2bNi	=ph2->Thermo->dGm1[1];	//VA dans austenite
	
	
double ddGm2bFeFe, ddGm2bFeNi, ddGm2bNiNi,ddGm2bNiFe;	
ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
ddGm2bFeFe  =ph2->Thermo->dG_1(0,0);
ddGm2bFeNi  =ph2->Thermo->dG_1(0,1);
ddGm2bNiNi  = ph2->Thermo->dG_1(1,0);	
ddGm2bNiFe  =ph2->Thermo->dG_1(1,1);

double bb = 2.*ddGm2bFeFe+2.*ddGm2bNiNi-2.*ddGm2bFeNi-2.*ddGm2bNiFe;
	

L = dcr*correc_D/bb;//pow(taille_gradient,2);
//printf("L=%e \n",L);
	
return L;
	
}


std::vector <double> interface::deriv_mu_FeCCr(double yc,double ycr,double x,vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);//bcc_A2: Fe,Cr
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: Fe,Cr
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	double G1,G2,h, dh;
	double dGm2aC,dGm2aVA,dGm2aFe,dGm2aCr,dGm2bC,dGm2bVA,dGm2bFe,dGm2bCr;
	double ddGm[2], dGma[2], dGmb[2];
	double XdGm_lib, XdGmCc_lib, XdGmCc_C, XdGmCc_CR, XdGmCcr_lib, XdGmCcr_CR, XdGmCcr_C;
	double part1_gam, part2_gam, part1_alf, part2_alf, \
	part11_gam, part22_gam, part11_alf, part22_alf, \
	XdGm_C, XdGm_CR;
	
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2aFeFe, ddGm2aFeCr, ddGm2aCrFe, ddGm2aCrCr;
	double ddGm2aFeC, ddGm2aFeVA, ddGm2aCrC, ddGm2aCrVA;
	
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	double ddGm2bFeFe, ddGm2bFeCr, ddGm2bCrFe, ddGm2bCrCr;
	double ddGm2bFeC, ddGm2bFeVA, ddGm2bCrC, ddGm2bCrVA;
	double Lc, Lcrcr, Lcx;
	double Xa[2],Xb[2];
	
	
	
	int i,j;
	std::vector < std::vector<double> > D(12);
	for(i=0;i<12;i++)
	{
		D[i].resize(4);
	}
	
	
	Y1alfa[1]  = ycr; //CR
	Y1alfa[0]  = 1.-Y1alfa[1]; //FE
	Y2alfa[0]  = yc/ph1->Thermo->m[1]; //C
	Y2alfa[1]  = 1.-Y2alfa[0]; //VA
	Y1betta[1] = ycr; //CR
	Y1betta[0] = 1.-Y1betta[1]; //FE
	Y2betta[0] = yc/ph2->Thermo->m[1]; //C
	Y2betta[1] = 1.-Y2betta[0]; //VA
	
	//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
	G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
	*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
	//cout<<"G1="<<G1<<endl;
	G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
	*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite
	//cout<<"G2="<<G2<<endl;
	
	//Derivation 1 order and 2 order y
	D=derivFeCrC_U(yc,ycr);
	
	
	h=Wang(x);
	dh=dWang(x);
	
	//Derivation energy gibbs 1 order
	ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
	dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
	//cout<<"dGm2aC="<<dGm2aC<<endl;
	dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans ferrite
	//cout<<"dGm2aVA="<<dGm2aVA<<endl;
	dGm2aFe	=ph1->Thermo->dGm1[0];	//Fe dans ferrite
	//cout<<"dGm2aFe="<<dGm2aFe<<endl;
	dGm2aCr	=ph1->Thermo->dGm1[1];	//Cr dans ferrite
	//cout<<"dGm2aCr="<<dGm2aCr<<endl;
	ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
	dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
	//cout<<"dGm2bC="<<dGm2bC<<endl;
	dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
	//cout<<"dGm2bVA="<<dGm2bVA<<endl;
	dGm2bFe	=ph2->Thermo->dGm1[0];	//Fe dans austenite
	//cout<<"dGm2bFe="<<dGm2bFe<<endl;
	dGm2bCr	=ph2->Thermo->dGm1[1];	//Cr dans austenite
	//cout<<"dGm2bCr="<<dGm2bCr<<endl;
	
	//Energy_gibbs	= 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0])*G2*h \
	+1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0])*G1*(1.-h);
	
	
	dGmb[0]=dGm2bC*D[2][2]+dGm2bVA*D[3][2]+dGm2bFe*D[0][2]+dGm2bCr*D[1][2]; //*dy/dcc (austenite)
	dGma[0]=dGm2aC*D[2][0]+dGm2aVA*D[3][0]+dGm2aFe*D[0][0]+dGm2aCr*D[1][0]; //*dy/dcc (ferrite)
	
	dGmb[1]=dGm2bC*D[2][3]+dGm2bVA*D[3][3]+dGm2bFe*D[0][3]+dGm2bCr*D[1][3]; //*dy/dccr (austenite)
	dGma[1]=dGm2aC*D[2][1]+dGm2aVA*D[3][1]+dGm2aFe*D[0][1]+dGm2aCr*D[1][1]; //*dy/dccr (ferrite)
	
	ddGm[0]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*dGma[0] ) );/*dGm/dcc*/
	
	ddGm[1]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1] ) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);/*dGm/dccr*/
	
	
	XdGm_lib	= ( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*G2 \
	-( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*G1;
	
	XdGm_C		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
	+(1.-h)*(- ( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
			 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0]);
	
	XdGm_CR		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
			 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);	
	
	//Derivation energy gibbs 2 order
	ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
	ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
	ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
	ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
	ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
	
	
	
	ddGm2aFeFe	=ph1->Thermo->dG_1(0,0);	//FE:FE dans ferrite
	ddGm2aFeCr	=ph1->Thermo->dG_1(0,1);	//FE:CR dans ferrite
	ddGm2aCrFe	=ph1->Thermo->dG_1(1,0);	//CR:FE dans ferrite
	ddGm2aCrCr	=ph1->Thermo->dG_1(1,1);	//CR:CR dans ferrite
	
	ddGm2aFeC	=ph1->Thermo->dG_2(0,0);	//FE:C dans ferrite
	ddGm2aFeVA	=ph1->Thermo->dG_2(0,1);	//FE:VA dans ferrite
	ddGm2aCrC	=ph1->Thermo->dG_2(1,0);	//CR:C dans ferrite
	ddGm2aCrVA	=ph1->Thermo->dG_2(1,1);	//CR:VA dans ferrite
	
	
	ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
	ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
	ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
	ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
	ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
	
	ddGm2bFeFe	=ph2->Thermo->dG_1(0,0);	//FE:FE dans austenite
	ddGm2bFeCr	=ph2->Thermo->dG_1(0,1);	//FE:CR dans austenite
	ddGm2bCrFe	=ph2->Thermo->dG_1(1,0);	//CR:FE dans austenite
	ddGm2bCrCr	=ph2->Thermo->dG_1(1,1);	//CR:CR dans austenite
	
	ddGm2bFeC	=ph2->Thermo->dG_2(0,0);	//FE:C dans austenite
	ddGm2bFeVA	=ph2->Thermo->dG_2(0,1);	//FE:VA dans austenite
	ddGm2bCrC	=ph2->Thermo->dG_2(1,0);	//CR:C dans austenite
	ddGm2bCrVA	=ph2->Thermo->dG_2(1,1);	//CR:VA dans austenite
	
	
	
	
	
	XdGmCc_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
	+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]	\
	+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
	-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0];
	
	//printf("XdGmCc_lib=%e \n",XdGmCc_lib);
	
	part11_gam	=ddGm2bFeFe*pow(D[0][2],2.)+ddGm2bCrCr*pow(D[1][2],2.) \
	+ddGm2bCC*pow(D[2][2],2.)+ddGm2bVAVA*pow(D[3][2],2.)	\
	+dGm2bFe*D[4][2]+dGm2bCr*D[5][2]+dGm2bC*D[6][2]+dGm2bVA*D[7][2]	\
	+D[0][2]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
	+D[1][2]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
	+D[2][2]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2]) \
	+D[3][2]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
	
	//printf("part11_gam = %e \n",part11_gam);
	
	part22_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
	+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
	+dGm2bFe*D[8][2]+dGm2bCr*D[9][2]+dGm2bC*D[10][2]+dGm2bVA*D[11][2] \
	+D[0][2]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
	+D[1][2]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
	+D[2][2]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3])	\
	+D[3][2]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
	
	//printf("part22_gam = %e \n", part22_gam);
	
	part11_alf	=ddGm2aFeFe*pow(D[0][0],2.)+ddGm2aCrCr*pow(D[1][0],2.) \
	+ddGm2aCC*pow(D[2][0],2.)+ddGm2aVAVA*pow(D[3][0],2.)	\
	+dGm2aFe*D[4][0]+dGm2aCr*D[5][0]+dGm2aC*D[6][0]+dGm2aVA*D[7][0]	\
	+D[0][0]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
	+D[1][0]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
	+D[2][0]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0]) \
	+D[3][0]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	//printf("part11_alf= %e \n",part11_alf);
	
	part22_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
	+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
	+dGm2aFe*D[8][0]+dGm2aCr*D[9][0]+dGm2aC*D[10][0]+dGm2aVA*D[11][0] \
	+D[0][0]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
	+D[1][0]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
	+D[2][0]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1])	\
	+D[3][0]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	//printf("part22_alf = %e \n",part22_alf);
	
	XdGmCc_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][2],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][2]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[0]
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part11_gam)\
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][0],2.)*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][0]*G1 \
			 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[0] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part11_alf);
	
	
	//printf("XdGmCc_c =%e \n",XdGmCc_C);
	
	XdGmCc_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][2]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part22_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][0]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part22_alf);
	
	//printf("XdGmCc_Cr= %e \n",XdGmCc_CR);
	
	part1_gam	=ddGm2bFeFe*pow(D[0][3],2.)+ddGm2bCrCr*pow(D[1][3],2.)+ddGm2bCC*pow(D[2][3],2.)+ddGm2bVAVA*pow(D[3][3],2.) \
	+dGm2bFe*D[4][3]+dGm2bCr*D[5][3]+dGm2bC*D[6][3]+dGm2bVA*D[7][3] \
	+D[0][3]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
	+D[1][3]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
	+D[2][3]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3]) \
	+D[3][3]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
	
	//printf("part1_gam = %e \n", part1_gam);
	
	part2_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
	+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
	+dGm2bFe*D[8][3]+dGm2bCr*D[9][3]+dGm2bC*D[10][3]+dGm2bVA*D[11][3] \
	+D[0][3]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
	+D[1][3]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
	+D[2][3]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2])	\
	+D[3][3]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
	
	//printf("part2_gam= %e \n", part2_gam);
	
	part1_alf	=ddGm2aFeFe*pow(D[0][1],2.)+ddGm2aCrCr*pow(D[1][1],2.)+ddGm2aCC*pow(D[2][1],2.)+ddGm2aVAVA*pow(D[3][1],2.) \
	+dGm2aFe*D[4][1]+dGm2aCr*D[5][1]+dGm2aC*D[6][1]+dGm2aVA*D[7][1] \
	+D[0][1]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
	+D[1][1]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
	+D[2][1]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1]) \
	+D[3][1]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	//printf("part1_alf = %e \n",part1_alf);
	
	part2_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
	+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
	+dGm2aFe*D[8][1]+dGm2aCr*D[9][1]+dGm2aC*D[10][1]+dGm2aVA*D[11][1] \
	+D[0][1]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
	+D[1][1]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
	+D[2][1]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0])	\
	+D[3][1]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	//printf("part2_alf = %e \n", part2_alf);
	
	XdGmCcr_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
	+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]	\
	+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
	-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1];
	
	
	//printf("XdGmCcr_lib = %e\n",XdGmCcr_lib);
	
	XdGmCcr_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][3],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][3]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part1_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][1],2.)*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][1]*G1 \
			 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[1] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part1_alf);
	
	//printf("XdGmCcr_CR = %e \n",XdGmCcr_CR);
	
	XdGmCcr_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part2_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part2_alf);
	
	//printf("XdGmCcr_C=%e \n",XdGmCcr_C);
 	
	//int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(2);	
	//Diffbetta = calc_diffusion_betta(concentration);
	//Diffbetta[0] = Diffbetta[0]*pow(1.e6,2.);//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//Diffbetta[1] = Diffbetta[1]*pow(1.e6,2.);//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	Diffbetta[0] = dc*coeff_n;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	
	Diffbetta[1] = dcr*coeff_k*correc_D;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	//Ispol'zovat' etu zapis' v sluchae polnoi sistemi uravnenij
	
	//Diffbetta[1] = dcr*coeff_k;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	
	
	//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	//printf("ddGm[0] = %e ddGm[1]=%e XdGmCc_C=%e XdGmCc_CR=%e XdGmCcr_CR=%e XdGmCcr_C=%e\n",ddGm[0],ddGm[1], \
	XdGmCc_C,XdGmCc_CR,XdGmCcr_CR,XdGmCcr_C);
	double mum_c, mum_cr;	
	mum_c	=(1.+yc)*(XdGmCc_C+XdGmCc_CR)+2*ddGm[0]+ddGm[1];
	mum_cr	=(1.+yc)*(XdGmCcr_CR+XdGmCcr_C)+ddGm[1];
	
	
	
	Lc=Diffbetta[0]/mum_c/Vm; 
	Lcrcr=Diffbetta[1]/mum_cr/Vm;
	Lcx=Lcoef*(Diffbetta[0]/mum_c/Vm)+Bcoef*(Diffbetta[1]/mum_cr/Vm);
	
	
	
	std::vector <double> L(3);
	L[0]=Lc;
	L[1]=Lcrcr;
	L[2]=Lcx;	
	return L;
	
}


std::vector <double> interface::deriv_mu_FeCCr_KTT_prof5(double yc,double ycr,double x,vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);//bcc_A2: Fe,Cr
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: Fe,Cr
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	double G1,G2,h, dh;
	double dGm2aC,dGm2aVA,dGm2aFe,dGm2aCr,dGm2bC,dGm2bVA,dGm2bFe,dGm2bCr;
	double ddGm[2], dGma[2], dGmb[2];
	double XdGm_lib, XdGmCc_lib, XdGmCc_C, XdGmCc_CR, XdGmCcr_lib, XdGmCcr_CR, XdGmCcr_C;
	double part1_gam, part2_gam, part1_alf, part2_alf, \
	part11_gam, part22_gam, part11_alf, part22_alf, \
	XdGm_C, XdGm_CR;
	
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2aFeFe, ddGm2aFeCr, ddGm2aCrFe, ddGm2aCrCr;
	double ddGm2aFeC, ddGm2aFeVA, ddGm2aCrC, ddGm2aCrVA;
	
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	double ddGm2bFeFe, ddGm2bFeCr, ddGm2bCrFe, ddGm2bCrCr;
	double ddGm2bFeC, ddGm2bFeVA, ddGm2bCrC, ddGm2bCrVA;
	double Lc, Lcrcr, Lcx;
	double Xa[2],Xb[2];
	
	
	
	int i,j;
	std::vector < std::vector<double> > D(12);
	for(i=0;i<12;i++)
	{
		D[i].resize(4);
	}
	
	
	Y1alfa[1]  = ycr; //CR
	Y1alfa[0]  = 1.-Y1alfa[1]; //FE
	Y2alfa[0]  = yc/ph1->Thermo->m[1]; //C
	Y2alfa[1]  = 1.-Y2alfa[0]; //VA
	Y1betta[1] = ycr; //CR
	Y1betta[0] = 1.-Y1betta[1]; //FE
	Y2betta[0] = yc/ph2->Thermo->m[1]; //C
	Y2betta[1] = 1.-Y2betta[0]; //VA
	
	//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
	G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
	*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
	//cout<<"G1="<<G1<<endl;
	G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
	*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite
	//cout<<"G2="<<G2<<endl;
	
	//Derivation 1 order and 2 order y
	D=derivFeCrC_U(yc,ycr);
	
	
	h=Wang(x);
	dh=dWang(x);
	
	//Derivation energy gibbs 1 order
	ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
	dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
	//cout<<"dGm2aC="<<dGm2aC<<endl;
	dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans ferrite
	//cout<<"dGm2aVA="<<dGm2aVA<<endl;
	dGm2aFe	=ph1->Thermo->dGm1[0];	//Fe dans ferrite
	//cout<<"dGm2aFe="<<dGm2aFe<<endl;
	dGm2aCr	=ph1->Thermo->dGm1[1];	//Cr dans ferrite
	//cout<<"dGm2aCr="<<dGm2aCr<<endl;
	ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
	dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
	//cout<<"dGm2bC="<<dGm2bC<<endl;
	dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
	//cout<<"dGm2bVA="<<dGm2bVA<<endl;
	dGm2bFe	=ph2->Thermo->dGm1[0];	//Fe dans austenite
	//cout<<"dGm2bFe="<<dGm2bFe<<endl;
	dGm2bCr	=ph2->Thermo->dGm1[1];	//Cr dans austenite
	//cout<<"dGm2bCr="<<dGm2bCr<<endl;
	
	//Energy_gibbs	= 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0])*G2*h \
	+1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0])*G1*(1.-h);
	
	
	dGmb[0]=dGm2bC*D[2][2]+dGm2bVA*D[3][2]+dGm2bFe*D[0][2]+dGm2bCr*D[1][2]; //*dy/dcc (austenite)
	dGma[0]=dGm2aC*D[2][0]+dGm2aVA*D[3][0]+dGm2aFe*D[0][0]+dGm2aCr*D[1][0]; //*dy/dcc (ferrite)
	
	dGmb[1]=dGm2bC*D[2][3]+dGm2bVA*D[3][3]+dGm2bFe*D[0][3]+dGm2bCr*D[1][3]; //*dy/dccr (austenite)
	dGma[1]=dGm2aC*D[2][1]+dGm2aVA*D[3][1]+dGm2aFe*D[0][1]+dGm2aCr*D[1][1]; //*dy/dccr (ferrite)
	
	ddGm[0]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*dGma[0] ) );/*dGm/dcc*/
	
	ddGm[1]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1] ) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);/*dGm/dccr*/
	
	
	XdGm_lib	= ( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*G2 \
	-( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*G1;
	
	XdGm_C		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
	+(1.-h)*(- ( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
			 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0]);
	
	XdGm_CR		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
			 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);	
	
	//Derivation energy gibbs 2 order
	ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
	ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
	ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
	ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
	ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
	
	
	
	ddGm2aFeFe	=ph1->Thermo->dG_1(0,0);	//FE:FE dans ferrite
	ddGm2aFeCr	=ph1->Thermo->dG_1(0,1);	//FE:CR dans ferrite
	ddGm2aCrFe	=ph1->Thermo->dG_1(1,0);	//CR:FE dans ferrite
	ddGm2aCrCr	=ph1->Thermo->dG_1(1,1);	//CR:CR dans ferrite
	
	ddGm2aFeC	=ph1->Thermo->dG_2(0,0);	//FE:C dans ferrite
	ddGm2aFeVA	=ph1->Thermo->dG_2(0,1);	//FE:VA dans ferrite
	ddGm2aCrC	=ph1->Thermo->dG_2(1,0);	//CR:C dans ferrite
	ddGm2aCrVA	=ph1->Thermo->dG_2(1,1);	//CR:VA dans ferrite
	
	
	ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
	ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
	ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
	ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
	ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
	
	ddGm2bFeFe	=ph2->Thermo->dG_1(0,0);	//FE:FE dans austenite
	ddGm2bFeCr	=ph2->Thermo->dG_1(0,1);	//FE:CR dans austenite
	ddGm2bCrFe	=ph2->Thermo->dG_1(1,0);	//CR:FE dans austenite
	ddGm2bCrCr	=ph2->Thermo->dG_1(1,1);	//CR:CR dans austenite
	
	ddGm2bFeC	=ph2->Thermo->dG_2(0,0);	//FE:C dans austenite
	ddGm2bFeVA	=ph2->Thermo->dG_2(0,1);	//FE:VA dans austenite
	ddGm2bCrC	=ph2->Thermo->dG_2(1,0);	//CR:C dans austenite
	ddGm2bCrVA	=ph2->Thermo->dG_2(1,1);	//CR:VA dans austenite
	
	
	
	
	
	XdGmCc_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
	+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]	\
	+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
	-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0];
	
	//printf("XdGmCc_lib=%e \n",XdGmCc_lib);
	
	part11_gam	=ddGm2bFeFe*pow(D[0][2],2.)+ddGm2bCrCr*pow(D[1][2],2.) \
	+ddGm2bCC*pow(D[2][2],2.)+ddGm2bVAVA*pow(D[3][2],2.)	\
	+dGm2bFe*D[4][2]+dGm2bCr*D[5][2]+dGm2bC*D[6][2]+dGm2bVA*D[7][2]	\
	+D[0][2]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
	+D[1][2]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
	+D[2][2]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2]) \
	+D[3][2]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
	
	//printf("part11_gam = %e \n",part11_gam);
	
	part22_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
	+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
	+dGm2bFe*D[8][2]+dGm2bCr*D[9][2]+dGm2bC*D[10][2]+dGm2bVA*D[11][2] \
	+D[0][2]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
	+D[1][2]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
	+D[2][2]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3])	\
	+D[3][2]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
	
	//printf("part22_gam = %e \n", part22_gam);
	
	part11_alf	=ddGm2aFeFe*pow(D[0][0],2.)+ddGm2aCrCr*pow(D[1][0],2.) \
	+ddGm2aCC*pow(D[2][0],2.)+ddGm2aVAVA*pow(D[3][0],2.)	\
	+dGm2aFe*D[4][0]+dGm2aCr*D[5][0]+dGm2aC*D[6][0]+dGm2aVA*D[7][0]	\
	+D[0][0]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
	+D[1][0]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
	+D[2][0]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0]) \
	+D[3][0]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	//printf("part11_alf= %e \n",part11_alf);
	
	part22_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
	+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
	+dGm2aFe*D[8][0]+dGm2aCr*D[9][0]+dGm2aC*D[10][0]+dGm2aVA*D[11][0] \
	+D[0][0]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
	+D[1][0]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
	+D[2][0]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1])	\
	+D[3][0]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	//printf("part22_alf = %e \n",part22_alf);
	
	XdGmCc_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][2],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][2]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[0]
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part11_gam)\
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][0],2.)*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][0]*G1 \
			 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[0] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part11_alf);
	
	
	//printf("XdGmCc_c =%e \n",XdGmCc_C);
	
	XdGmCc_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][2]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part22_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][0]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part22_alf);
	
	//printf("XdGmCc_Cr= %e \n",XdGmCc_CR);
	
	part1_gam	=ddGm2bFeFe*pow(D[0][3],2.)+ddGm2bCrCr*pow(D[1][3],2.)+ddGm2bCC*pow(D[2][3],2.)+ddGm2bVAVA*pow(D[3][3],2.) \
	+dGm2bFe*D[4][3]+dGm2bCr*D[5][3]+dGm2bC*D[6][3]+dGm2bVA*D[7][3] \
	+D[0][3]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
	+D[1][3]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
	+D[2][3]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3]) \
	+D[3][3]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
	
	//printf("part1_gam = %e \n", part1_gam);
	
	part2_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
	+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
	+dGm2bFe*D[8][3]+dGm2bCr*D[9][3]+dGm2bC*D[10][3]+dGm2bVA*D[11][3] \
	+D[0][3]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
	+D[1][3]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
	+D[2][3]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2])	\
	+D[3][3]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
	
	//printf("part2_gam= %e \n", part2_gam);
	
	part1_alf	=ddGm2aFeFe*pow(D[0][1],2.)+ddGm2aCrCr*pow(D[1][1],2.)+ddGm2aCC*pow(D[2][1],2.)+ddGm2aVAVA*pow(D[3][1],2.) \
	+dGm2aFe*D[4][1]+dGm2aCr*D[5][1]+dGm2aC*D[6][1]+dGm2aVA*D[7][1] \
	+D[0][1]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
	+D[1][1]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
	+D[2][1]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1]) \
	+D[3][1]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	//printf("part1_alf = %e \n",part1_alf);
	
	part2_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
	+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
	+dGm2aFe*D[8][1]+dGm2aCr*D[9][1]+dGm2aC*D[10][1]+dGm2aVA*D[11][1] \
	+D[0][1]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
	+D[1][1]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
	+D[2][1]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0])	\
	+D[3][1]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	//printf("part2_alf = %e \n", part2_alf);
	
	XdGmCcr_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
	+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]	\
	+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
	-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1];
	
	
	//printf("XdGmCcr_lib = %e\n",XdGmCcr_lib);
	
	XdGmCcr_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][3],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][3]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part1_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][1],2.)*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][1]*G1 \
			 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[1] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part1_alf);
	
	//printf("XdGmCcr_CR = %e \n",XdGmCcr_CR);
	
	XdGmCcr_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part2_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part2_alf);
	
	//printf("XdGmCcr_C=%e \n",XdGmCcr_C);
 	
	//int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(2);	
	//Diffbetta = calc_diffusion_betta(concentration);
	//Diffbetta[0] = Diffbetta[0]*pow(1.e6,2.);//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//Diffbetta[1] = Diffbetta[1]*pow(1.e6,2.);//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	

	Diffbetta[0] = dc*coeff_n;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	
	Diffbetta[1] = dcr*correc_D*coeff_k;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	//Ispol'zovat' etu zapis' v sluchae polnoi sistemi uravnenij
	
	//Diffbetta[1] = dcr*coeff_k;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;


	
	
	//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	//printf("ddGm[0] = %e ddGm[1]=%e XdGmCc_C=%e XdGmCc_CR=%e XdGmCcr_CR=%e XdGmCcr_C=%e\n",ddGm[0],ddGm[1], \
		   XdGmCc_C,XdGmCc_CR,XdGmCcr_CR,XdGmCcr_C);
	double mum_c, mum_cr;	
	mum_c	=(1.+yc)*(XdGmCc_C+XdGmCc_CR)+2*ddGm[0]+ddGm[1];
	mum_cr	=(1.+yc)*(XdGmCcr_CR+XdGmCcr_C)+ddGm[1];
	
	
	Lc=Diffbetta[0]/mum_c/Vm; 
	Lcrcr=Diffbetta[1]/mum_cr/Vm;
	Lcx=Lcoef*(Diffbetta[0]/mum_c/Vm)+Bcoef*(Diffbetta[1]/mum_cr/Vm);


	
	std::vector <double> L(3);
	L[0]=Lc;
	L[1]=Lcrcr;
	L[2]=Lcx;	
	return L;

}




std::vector <double> interface::deriv_mu_FeCCr_X(double yc,double ycr,double x,vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);//bcc_A2: Fe,Cr
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: Fe,Cr
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	double G1,G2,h, dh;
	double dGm2aC,dGm2aVA,dGm2aFe,dGm2aCr,dGm2bC,dGm2bVA,dGm2bFe,dGm2bCr;
	double ddGm[2], dGma[2], dGmb[2];
	double XdGm_lib, XdGmCc_lib, XdGmCc_C, XdGmCc_CR, XdGmCcr_lib, XdGmCcr_CR, XdGmCcr_C;
	double part1_gam, part2_gam, part1_alf, part2_alf, \
	part11_gam, part22_gam, part11_alf, part22_alf, \
	XdGm_C, XdGm_CR;
	
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2aFeFe, ddGm2aFeCr, ddGm2aCrFe, ddGm2aCrCr;
	double ddGm2aFeC, ddGm2aFeVA, ddGm2aCrC, ddGm2aCrVA;
	
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	double ddGm2bFeFe, ddGm2bFeCr, ddGm2bCrFe, ddGm2bCrCr;
	double ddGm2bFeC, ddGm2bFeVA, ddGm2bCrC, ddGm2bCrVA;
	double Lc, Lcrcr, Lcx;
	double Xa[2],Xb[2];
	
	
	
	int i,j;
	std::vector < std::vector<double> > D(12);
	for(i=0;i<12;i++)
	{
		D[i].resize(4);
	}
	
	
	Y1alfa[1]  = ycr; //CR
	Y1alfa[0]  = 1.-Y1alfa[1]; //FE
	Y2alfa[0]  = yc/ph1->Thermo->m[1]; //C
	Y2alfa[1]  = 1.-Y2alfa[0]; //VA
	Y1betta[1] = ycr; //CR
	Y1betta[0] = 1.-Y1betta[1]; //FE
	Y2betta[0] = yc/ph2->Thermo->m[1]; //C
	Y2betta[1] = 1.-Y2betta[0]; //VA
	
	//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
	G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
	*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
	//cout<<"G1="<<G1<<endl;
	G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
	*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite
	//cout<<"G2="<<G2<<endl;
	
	//Derivation 1 order and 2 order y
	D=derivFeCrC_U(yc,ycr);
	
	
	h=Wang(x);
	dh=dWang(x);
	
	//Derivation energy gibbs 1 order
	ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
	dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
	//cout<<"dGm2aC="<<dGm2aC<<endl;
	dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans ferrite
	//cout<<"dGm2aVA="<<dGm2aVA<<endl;
	dGm2aFe	=ph1->Thermo->dGm1[0];	//Fe dans ferrite
	//cout<<"dGm2aFe="<<dGm2aFe<<endl;
	dGm2aCr	=ph1->Thermo->dGm1[1];	//Cr dans ferrite
	//cout<<"dGm2aCr="<<dGm2aCr<<endl;
	ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
	dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
	//cout<<"dGm2bC="<<dGm2bC<<endl;
	dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
	//cout<<"dGm2bVA="<<dGm2bVA<<endl;
	dGm2bFe	=ph2->Thermo->dGm1[0];	//Fe dans austenite
	//cout<<"dGm2bFe="<<dGm2bFe<<endl;
	dGm2bCr	=ph2->Thermo->dGm1[1];	//Cr dans austenite
	//cout<<"dGm2bCr="<<dGm2bCr<<endl;
	
	//Energy_gibbs	= 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0])*G2*h \
	+1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0])*G1*(1.-h);
	
	
	dGmb[0]=dGm2bC*D[2][2]+dGm2bVA*D[3][2]+dGm2bFe*D[0][2]+dGm2bCr*D[1][2]; //*dy/dcc (austenite)
	dGma[0]=dGm2aC*D[2][0]+dGm2aVA*D[3][0]+dGm2aFe*D[0][0]+dGm2aCr*D[1][0]; //*dy/dcc (ferrite)
	
	dGmb[1]=dGm2bC*D[2][3]+dGm2bVA*D[3][3]+dGm2bFe*D[0][3]+dGm2bCr*D[1][3]; //*dy/dccr (austenite)
	dGma[1]=dGm2aC*D[2][1]+dGm2aVA*D[3][1]+dGm2aFe*D[0][1]+dGm2aCr*D[1][1]; //*dy/dccr (ferrite)
	
	ddGm[0]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*dGma[0] ) );/*dGm/dcc*/
	
	ddGm[1]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1] ) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);/*dGm/dccr*/
	
	
	XdGm_lib	= ( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*G2 \
	-( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*G1;
	
	XdGm_C		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
	+(1.-h)*(- ( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
			 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0]);
	
	XdGm_CR		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
			 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);	
	
	//Derivation energy gibbs 2 order
	ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
	ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
	ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
	ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
	ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
	
	
	
	ddGm2aFeFe	=ph1->Thermo->dG_1(0,0);	//FE:FE dans ferrite
	ddGm2aFeCr	=ph1->Thermo->dG_1(0,1);	//FE:CR dans ferrite
	ddGm2aCrFe	=ph1->Thermo->dG_1(1,0);	//CR:FE dans ferrite
	ddGm2aCrCr	=ph1->Thermo->dG_1(1,1);	//CR:CR dans ferrite
	
	ddGm2aFeC	=ph1->Thermo->dG_2(0,0);	//FE:C dans ferrite
	ddGm2aFeVA	=ph1->Thermo->dG_2(0,1);	//FE:VA dans ferrite
	ddGm2aCrC	=ph1->Thermo->dG_2(1,0);	//CR:C dans ferrite
	ddGm2aCrVA	=ph1->Thermo->dG_2(1,1);	//CR:VA dans ferrite
	
	
	ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
	ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
	ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
	ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
	ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
	
	ddGm2bFeFe	=ph2->Thermo->dG_1(0,0);	//FE:FE dans austenite
	ddGm2bFeCr	=ph2->Thermo->dG_1(0,1);	//FE:CR dans austenite
	ddGm2bCrFe	=ph2->Thermo->dG_1(1,0);	//CR:FE dans austenite
	ddGm2bCrCr	=ph2->Thermo->dG_1(1,1);	//CR:CR dans austenite
	
	ddGm2bFeC	=ph2->Thermo->dG_2(0,0);	//FE:C dans austenite
	ddGm2bFeVA	=ph2->Thermo->dG_2(0,1);	//FE:VA dans austenite
	ddGm2bCrC	=ph2->Thermo->dG_2(1,0);	//CR:C dans austenite
	ddGm2bCrVA	=ph2->Thermo->dG_2(1,1);	//CR:VA dans austenite
	
	
	
	
	
	XdGmCc_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
	+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]	\
	+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
	-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0];
	
	part11_gam	=ddGm2bFeFe*pow(D[0][2],2.)+ddGm2bCrCr*pow(D[1][2],2.) \
	+ddGm2bCC*pow(D[2][2],2.)+ddGm2bVAVA*pow(D[3][2],2.)	\
	+dGm2bFe*D[4][2]+dGm2bCr*D[5][2]+dGm2bC*D[6][2]+dGm2bVA*D[7][2]	\
	+D[0][2]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
	+D[1][2]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
	+D[2][2]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2]) \
	+D[3][2]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
	
	
	
	part22_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
	+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
	+dGm2bFe*D[8][2]+dGm2bCr*D[9][2]+dGm2bC*D[10][2]+dGm2bVA*D[11][2] \
	+D[0][2]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
	+D[1][2]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
	+D[2][2]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3])	\
	+D[3][2]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
	
	
	part11_alf	=ddGm2aFeFe*pow(D[0][0],2.)+ddGm2aCrCr*pow(D[1][0],2.) \
	+ddGm2aCC*pow(D[2][0],2.)+ddGm2aVAVA*pow(D[3][0],2.)	\
	+dGm2aFe*D[4][0]+dGm2aCr*D[5][0]+dGm2aC*D[6][0]+dGm2aVA*D[7][0]	\
	+D[0][0]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
	+D[1][0]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
	+D[2][0]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0]) \
	+D[3][0]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	part22_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
	+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
	+dGm2aFe*D[8][0]+dGm2aCr*D[9][0]+dGm2aC*D[10][0]+dGm2aVA*D[11][0] \
	+D[0][0]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
	+D[1][0]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
	+D[2][0]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1])	\
	+D[3][0]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	XdGmCc_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][2],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][2]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[0]
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part11_gam)\
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][0],2.)*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][0]*G1 \
			 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[0] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part11_alf);
	
	
	
	XdGmCc_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][2]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part22_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][0]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part22_alf);
	
	
	part1_gam	=ddGm2bFeFe*pow(D[0][3],2.)+ddGm2bCrCr*pow(D[1][3],2.)+ddGm2bCC*pow(D[2][3],2.)+ddGm2bVAVA*pow(D[3][3],2.) \
	+dGm2bFe*D[4][3]+dGm2bCr*D[5][3]+dGm2bC*D[6][3]+dGm2bVA*D[7][3] \
	+D[0][3]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
	+D[1][3]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
	+D[2][3]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3]) \
	+D[3][3]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
	
	part2_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
	+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
	+dGm2bFe*D[8][3]+dGm2bCr*D[9][3]+dGm2bC*D[10][3]+dGm2bVA*D[11][3] \
	+D[0][3]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
	+D[1][3]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
	+D[2][3]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2])	\
	+D[3][3]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
	
	part1_alf	=ddGm2aFeFe*pow(D[0][1],2.)+ddGm2aCrCr*pow(D[1][1],2.)+ddGm2aCC*pow(D[2][1],2.)+ddGm2aVAVA*pow(D[3][1],2.) \
	+dGm2aFe*D[4][1]+dGm2aCr*D[5][1]+dGm2aC*D[6][1]+dGm2aVA*D[7][1] \
	+D[0][1]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
	+D[1][1]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
	+D[2][1]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1]) \
	+D[3][1]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	part2_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
	+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
	+dGm2aFe*D[8][1]+dGm2aCr*D[9][1]+dGm2aC*D[10][1]+dGm2aVA*D[11][1] \
	+D[0][1]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
	+D[1][1]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
	+D[2][1]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0])	\
	+D[3][1]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	XdGmCcr_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
	+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]	\
	+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
	-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1];
	
	XdGmCcr_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][3],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][3]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part1_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][1],2.)*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][1]*G1 \
			 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[1] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part1_alf);
	
	XdGmCcr_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part2_gam) \
	+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][1]*G1 \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
			 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part2_alf);
	
 	
	//int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(2);	
	//Diffbetta = calc_diffusion_betta(concentration);
	//Diffbetta[0] = Diffbetta[0]*pow(1.e6,2.);//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//Diffbetta[1] = Diffbetta[1]*pow(1.e6,2.);//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	Diffbetta[0] = dcr*correc_D*coeff_n;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr*coeff_k;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	//Diffbetta[0] = Diffbetta[1]*correc_D*10;
	//Diffbetta[1] = Diffbetta[1];
	
	
	//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	double mum_c_Uc,mum_c_Ucr,mum_cr_Uc,mum_cr_Ucr;	
	mum_c_Uc= 2*ddGm[0]+(1.+yc)*XdGmCc_C;
	mum_c_Ucr=ddGm[1]+(1.+yc)*XdGmCc_CR;
	mum_cr_Ucr=(1.+yc)*XdGmCcr_CR;
	mum_cr_Uc=ddGm[1]+(1.+yc)*XdGmCcr_C;

	
	
	std::vector <double> X(2);
	X[0]=mum_c_Ucr/mum_c_Uc;
	X[1]=mum_cr_Uc/mum_cr_Ucr;
	return X;

}


double interface::funcFeCrC(int k, double x, double *y,vector <double> &concentration, vector <double> &Ysol)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);//bcc_A2: Fe,Cr
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: Fe,Cr
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	double G1,G2,h, dh;
	double dGm2aC,dGm2aVA,dGm2aFe,dGm2aCr,dGm2bC,dGm2bVA,dGm2bFe,dGm2bCr;
	double ddGm[2], dGma[2], dGmb[2];
	double XdGm_lib, XdGmCc_lib, XdGmCc_C, XdGmCc_CR, XdGmCcr_lib, XdGmCcr_CR, XdGmCcr_C;
	double part1_gam, part2_gam, part1_alf, part2_alf, \
			part11_gam, part22_gam, part11_alf, part22_alf, \
			XdGm_C, XdGm_CR;
	
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2aFeFe, ddGm2aFeCr, ddGm2aCrFe, ddGm2aCrCr;
	double ddGm2aFeC, ddGm2aFeVA, ddGm2aCrC, ddGm2aCrVA;
	
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	double ddGm2bFeFe, ddGm2bFeCr, ddGm2bCrFe, ddGm2bCrCr;
	double ddGm2bFeC, ddGm2bFeVA, ddGm2bCrC, ddGm2bCrVA;
	double Lc, Lcrcr;
	double Xa[2],Xb[2];


	
	int i,j;
	std::vector < std::vector<double> > D(12);
	for(i=0;i<12;i++)
	{
		D[i].resize(4);
	}
	
	
	Y1alfa[1]  = y[1]; //CR
	Y1alfa[0]  = 1.-Y1alfa[1]; //FE
	Y2alfa[0]  = y[0]/ph1->Thermo->m[1]; //C
	Y2alfa[1]  = 1.-Y2alfa[0]; //VA
	Y1betta[1] = y[1]; //CR
	Y1betta[0] = 1.-Y1betta[1]; //FE
	Y2betta[0] = y[0]/ph2->Thermo->m[1]; //C
	Y2betta[1] = 1.-Y2betta[0]; //VA
	
	//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
	G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
	*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
	//cout<<"G1="<<G1<<endl;
	G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
	*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite
	//cout<<"G2="<<G2<<endl;
	
	//Derivation 1 order and 2 order y
	D=derivFeCrC_U(y[0],y[1]);
	
	
	h=Wang(x);
	dh=dWang(x);
	
	//Derivation energy gibbs 1 order
	ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
	dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
	//cout<<"dGm2aC="<<dGm2aC<<endl;
	dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans ferrite
	//cout<<"dGm2aVA="<<dGm2aVA<<endl;
	dGm2aFe	=ph1->Thermo->dGm1[0];	//Fe dans ferrite
	//cout<<"dGm2aFe="<<dGm2aFe<<endl;
	dGm2aCr	=ph1->Thermo->dGm1[1];	//Cr dans ferrite
	//cout<<"dGm2aCr="<<dGm2aCr<<endl;
	ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
	dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
	//cout<<"dGm2bC="<<dGm2bC<<endl;
	dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
	//cout<<"dGm2bVA="<<dGm2bVA<<endl;
	dGm2bFe	=ph2->Thermo->dGm1[0];	//Fe dans austenite
	//cout<<"dGm2bFe="<<dGm2bFe<<endl;
	dGm2bCr	=ph2->Thermo->dGm1[1];	//Cr dans austenite
	//cout<<"dGm2bCr="<<dGm2bCr<<endl;
	
	//Energy_gibbs	= 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0])*G2*h \
	+1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0])*G1*(1.-h);
	
	
	dGmb[0]=dGm2bC*D[2][2]+dGm2bVA*D[3][2]+dGm2bFe*D[0][2]+dGm2bCr*D[1][2]; //*dy/dcc (austenite)
	//cout<<"dGmb[0]="<<dGmb[0]<<endl;
	dGma[0]=dGm2aC*D[2][0]+dGm2aVA*D[3][0]+dGm2aFe*D[0][0]+dGm2aCr*D[1][0]; //*dy/dcc (ferrite)
	//cout<<"dGma[0]="<<dGma[0]<<endl;

	dGmb[1]=dGm2bC*D[2][3]+dGm2bVA*D[3][3]+dGm2bFe*D[0][3]+dGm2bCr*D[1][3]; //*dy/dccr (austenite)
	//cout<<"dGmb[1]="<<dGmb[1]<<endl;

	dGma[1]=dGm2aC*D[2][1]+dGm2aVA*D[3][1]+dGm2aFe*D[0][1]+dGm2aCr*D[1][1]; //*dy/dccr (ferrite)
	//cout<<"dGma[1]="<<dGma[1]<<endl;

	
	ddGm[0]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*dGma[0] ) );/*dGm/dcc*/
	//cout <<"ddGm[0]="<<ddGm[0]<<endl;
	
	ddGm[1]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
			   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1] ) \
	+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
			 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);/*dGm/dccr*/
	//cout <<"ddGm[1]="<<ddGm[1]<<endl;

	
	XdGm_lib	= ( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*G2 \
					-( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*G1;
	//cout<<"XdGm_lib="<<XdGm_lib<<endl;
	
	XdGm_C		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
				+(1.-h)*(- ( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
						 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0]);

	//cout<<"XdGm_C="<<XdGm_C<<endl;

	
	XdGm_CR		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
					 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]) \
				+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
						 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);	
	
	//cout<<"XdGm_CR="<<XdGm_CR<<endl;

	//Derivation energy gibbs 2 order
	ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
	ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
	ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
	ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
	ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
	
	//cout<<"ddGm2aCC="<<ddGm2aCC<<endl;
	//cout<<"ddGm2aCVA="<<ddGm2aCVA<<endl;
	//cout<<"ddGm2aVAC="<<ddGm2aVAC<<endl;
	//cout<<"ddGm2aVAVA="<<ddGm2aVAVA<<endl;

	
	ddGm2aFeFe	=ph1->Thermo->dG_1(0,0);	//FE:FE dans ferrite
	ddGm2aFeCr	=ph1->Thermo->dG_1(0,1);	//FE:CR dans ferrite
	ddGm2aCrFe	=ph1->Thermo->dG_1(1,0);	//CR:FE dans ferrite
	ddGm2aCrCr	=ph1->Thermo->dG_1(1,1);	//CR:CR dans ferrite

	//cout<<"ddGm2aFeFe="<<ddGm2aFeFe<<endl;
	//cout<<"ddGm2aFeCr="<<ddGm2aFeCr<<endl;
	//cout<<"ddGm2aCrFe="<<ddGm2aCrFe<<endl;
	//cout<<"ddGm2aCrCr="<<ddGm2aCrCr<<endl;
	
	ddGm2aFeC	=ph1->Thermo->dG_2(0,0);	//FE:C dans ferrite
	ddGm2aFeVA	=ph1->Thermo->dG_2(0,1);	//FE:VA dans ferrite
	ddGm2aCrC	=ph1->Thermo->dG_2(1,0);	//CR:C dans ferrite
	ddGm2aCrVA	=ph1->Thermo->dG_2(1,1);	//CR:VA dans ferrite

	//cout<<"ddGm2aFeC="<<ddGm2aFeC<<endl;
	//cout<<"ddGm2aFeVA="<<ddGm2aFeVA<<endl;
	//cout<<"ddGm2aCrC="<<ddGm2aCrC<<endl;
	//cout<<"ddGm2aCrVA="<<ddGm2aCrVA<<endl;
	
	ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
	ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
	ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
	ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
	ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
	
	/*cout<<"ddGm2bCC="<<ddGm2bCC<<endl;
	cout<<"ddGm2bCVA="<<ddGm2bCVA<<endl;
	cout<<"ddGm2bVAC="<<ddGm2bVAC<<endl;
	cout<<"ddGm2bVAVA="<<ddGm2bVAVA<<endl;
	*/
	ddGm2bFeFe	=ph2->Thermo->dG_1(0,0);	//FE:FE dans austenite
	ddGm2bFeCr	=ph2->Thermo->dG_1(0,1);	//FE:CR dans austenite
	ddGm2bCrFe	=ph2->Thermo->dG_1(1,0);	//CR:FE dans austenite
	ddGm2bCrCr	=ph2->Thermo->dG_1(1,1);	//CR:CR dans austenite
	
	/*cout<<"ddGm2bFeFe="<<ddGm2bFeFe<<endl;
	cout<<"ddGm2bFeCr="<<ddGm2bFeCr<<endl;
	cout<<"ddGm2bCrFe="<<ddGm2bCrFe<<endl;
	cout<<"ddGm2bCrCr="<<ddGm2bCrCr<<endl;
	*/
	ddGm2bFeC	=ph2->Thermo->dG_2(0,0);	//FE:C dans austenite
	ddGm2bFeVA	=ph2->Thermo->dG_2(0,1);	//FE:VA dans austenite
	ddGm2bCrC	=ph2->Thermo->dG_2(1,0);	//CR:C dans austenite
	ddGm2bCrVA	=ph2->Thermo->dG_2(1,1);	//CR:VA dans austenite
	
	/*cout<<"ddGm2bFeC="<<ddGm2bFeC<<endl;
	cout<<"ddGm2bFeVA="<<ddGm2bFeVA<<endl;
	cout<<"ddGm2bCrC="<<ddGm2bCrC<<endl;
	cout<<"ddGm2bCrVA="<<ddGm2bCrVA<<endl;
	*/
	
	
	XdGmCc_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
					+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]	\
					+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
					-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0];
	
	part11_gam	=ddGm2bFeFe*pow(D[0][2],2.)+ddGm2bCrCr*pow(D[1][2],2.) \
				+ddGm2bCC*pow(D[2][2],2.)+ddGm2bVAVA*pow(D[3][2],2.)	\
				+dGm2bFe*D[4][2]+dGm2bCr*D[5][2]+dGm2bC*D[6][2]+dGm2bVA*D[7][2]	\
				+D[0][2]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
				+D[1][2]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
				+D[2][2]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2]) \
				+D[3][2]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);

	
	
	part22_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
				 +ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
				 +dGm2bFe*D[8][2]+dGm2bCr*D[9][2]+dGm2bC*D[10][2]+dGm2bVA*D[11][2] \
				 +D[0][2]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
				 +D[1][2]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
				 +D[2][2]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3])	\
				 +D[3][2]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);

	
	part11_alf	=ddGm2aFeFe*pow(D[0][0],2.)+ddGm2aCrCr*pow(D[1][0],2.) \
				 +ddGm2aCC*pow(D[2][0],2.)+ddGm2aVAVA*pow(D[3][0],2.)	\
				 +dGm2aFe*D[4][0]+dGm2aCr*D[5][0]+dGm2aC*D[6][0]+dGm2aVA*D[7][0]	\
				 +D[0][0]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
				 +D[1][0]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
				 +D[2][0]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0]) \
				 +D[3][0]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	part22_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
				 +ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
				 +dGm2aFe*D[8][0]+dGm2aCr*D[9][0]+dGm2aC*D[10][0]+dGm2aVA*D[11][0] \
				 +D[0][0]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
				 +D[1][0]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
				 +D[2][0]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1])	\
				 +D[3][0]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	XdGmCc_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][2],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][2]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[0]
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part11_gam)\
			+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][0],2.)*G1 \
					 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][0]*G1 \
					 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[0] \
					 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part11_alf);

	
	
	XdGmCc_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][2]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part22_gam) \
			+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
					 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][0]*G1 \
					 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
					 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
					 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part22_alf);


	part1_gam	=ddGm2bFeFe*pow(D[0][3],2.)+ddGm2bCrCr*pow(D[1][3],2.)+ddGm2bCC*pow(D[2][3],2.)+ddGm2bVAVA*pow(D[3][3],2.) \
				 +dGm2bFe*D[4][3]+dGm2bCr*D[5][3]+dGm2bC*D[6][3]+dGm2bVA*D[7][3] \
				 +D[0][3]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
				 +D[1][3]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
				 +D[2][3]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3]) \
				 +D[3][3]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
	
	part2_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
				 +ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
				 +dGm2bFe*D[8][3]+dGm2bCr*D[9][3]+dGm2bC*D[10][3]+dGm2bVA*D[11][3] \
				 +D[0][3]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
				 +D[1][3]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
				 +D[2][3]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2])	\
				 +D[3][3]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
	
	part1_alf	=ddGm2aFeFe*pow(D[0][1],2.)+ddGm2aCrCr*pow(D[1][1],2.)+ddGm2aCC*pow(D[2][1],2.)+ddGm2aVAVA*pow(D[3][1],2.) \
				 +dGm2aFe*D[4][1]+dGm2aCr*D[5][1]+dGm2aC*D[6][1]+dGm2aVA*D[7][1] \
				 +D[0][1]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
				 +D[1][1]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
				 +D[2][1]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1]) \
				 +D[3][1]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
	
	part2_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
				 +ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
				 +dGm2aFe*D[8][1]+dGm2aCr*D[9][1]+dGm2aC*D[10][1]+dGm2aVA*D[11][1] \
				 +D[0][1]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
				 +D[1][1]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
				 +D[2][1]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0])	\
				 +D[3][1]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
	
	XdGmCcr_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
				  +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]	\
				  +( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
				  -( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1];
	
	XdGmCcr_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][3],2.)*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][3]*G2 \
					 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part1_gam) \
			 +(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][1],2.)*G1 \
					-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][1]*G1 \
					-2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[1] \
					+( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part1_alf);
	
	XdGmCcr_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][3]*G2 \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
					 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
					 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part2_gam) \
			+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
					-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][1]*G1 \
					-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
					-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
					+( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part2_alf);
	
	//Fraction molaire de C
	Xa[0]=x_c0;
	Xa[1]=x_cr0;
	
	Xb[0]=(ph2->Thermo->m[1]*Y2betta[0]) \
	/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y2betta[1]);//carbone dans austenite
	Xb[1]=(ph2->Thermo->m[0]*Y1betta[1]) \
	/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y2betta[1]);//chrome dans austenite

	
	double AA,BB,CC,DD,EE,GG,KK,LL,FF,PP;
	
	//int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(2);	
	//Diffbetta = calc_diffusion_betta(concentration);
	//Diffbetta[0] = Diffbetta[0]*pow(1.e6,2.);//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//Diffbetta[1] = Diffbetta[1]*pow(1.e6,2.);//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	//Diffbetta[0] = Diffbetta[1]*correc_D;
	//Diffbetta[1] = Diffbetta[1];
	
	Diffbetta[0] = dc*coeff_n;
	//Diffbetta[1] = dcr*coeff_k;
	//printf("dx = %e\n",Diffbetta[1]);
	
	//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;

	double mum_c, mum_cr;	
	mum_c	=(1.+y[0])*(XdGmCc_C+XdGmCc_CR)+2*ddGm[0]+ddGm[1];
	//printf("muc = %e y[0]=%e XdGmCc_C=%e XdGmCc_CR=%e ddGm[0]=%e ddGm[1]=%e \n",mum_c, \
		   y[0],XdGmCc_C,XdGmCc_CR,ddGm[0],ddGm[1]);
	mum_cr	=(1.+y[0])*(XdGmCcr_CR+XdGmCcr_C)+ddGm[1];
	//printf("mu_cr = %e \n",mum_cr);
	Lc=Diffbetta[0]/mum_c/pow(taille_gradient,2)/Vm; 
	Lcrcr=Diffbetta[1]/mum_cr/pow(taille_gradient,2)/Vm;
	double VIT_TEMPORAIRE_=VIT_TEMPORAIRE/taille_gradient;
	//printf("Lc=%e \n Lxx=%e \n",Lc,Lcrcr);
	
	
	
	AA	=	Vm*Lc*dh*(XdGm_lib+(1.+y[0])*XdGmCc_lib);
	//printf("Lc= %e XdGm_lib= %e XdGmCc_lib= %e XdGmCcr_lib= %e dh= %e \n",Lc,XdGm_lib,XdGmCc_lib,XdGmCcr_lib,dh);
	BB	=	Vm*Lc*(XdGm_C+ddGm[0]+(1.+y[0])*XdGmCc_C);
	//printf("Lc= %e XdGm_CR= %e ddGm[1]= %e XdGmCc_CR= %e XdGmCcr_CR= %e \n",Lc,XdGm_C,ddGm[0],XdGmCc_C,XdGmCcr_C);
	CC	=	Vm*Lc*(XdGm_CR+(1.+y[0])*XdGmCc_CR);
	//printf("Lc= %e XdGm_CR= %e ddGm[1]= %e XdGmCc_CR= %e XdGmCcr_CR= %e \n",Lc,XdGm_CR,ddGm[1],XdGmCc_CR,XdGmCcr_CR);
	DD	=	Vm*Lcrcr*(dh*(1.+y[0])*XdGmCcr_lib);
	//printf("Lcrcr= %e  dh = %e  XdGmCcr_lib = %e \n",Lcrcr,dh,XdGmCcr_lib);
	EE	=	Vm*Lcrcr*((1.+y[0])*XdGmCcr_C+ddGm[1]);
	//printf("Lcrcr= %e  XdGmCcr_C = %e \n",Lcrcr,XdGmCcr_C);
	GG	=	Vm*Lcrcr*(1.+y[0])*(XdGmCcr_CR);
	//printf("Lcrcr= %e  XdGmCcr_CR = %e \n",Lcrcr,XdGmCcr_CR);
	KK	=	-VIT_TEMPORAIRE_*y[0]+VIT_TEMPORAIRE_*Xa[0];
	//printf("KK= %e y[0]= %e Xa[0]= %e \n",KK,y[0],Xa[0]);
	LL	=	-VIT_TEMPORAIRE_*y[1]+VIT_TEMPORAIRE_*Xa[1];
	//printf("LL= %e y[1]= %e Xa[1]= %e \n",LL,y[1],Xa[1]);
	double dominateur = (GG*BB-CC*EE);
	if (dominateur==0.0) dominateur = 1.e-9;
	FF	=	((LL-DD)*BB-(KK-AA)*EE)/dominateur; 
	//printf ("FF= %e , GG*BB-CC*EE = %e\n",FF,(GG*BB-CC*EE));
	if(BB==0.0)BB=1.e-9;
	PP	=	(KK-AA-FF*CC)/BB; 	
	//printf("PP= %e BB = %e \n",PP,BB);
	
	
	
	switch(k) { 
		case 0: return PP;
		case 1: return FF;
		default: return 0;
	}
	
}

void interface::funcderivmuFeCrC(int m, double xi, double xf,double *yc,double *ycr,vector <double> &concentration,vector <double> &Ysol, \
					  std::vector < std::vector<double> > &P)
{
	double h1,x;
	h1=(xf-xi)/(m-1); 
	x=xi-h1;
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);//bcc_A2: Fe,Cr
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: Fe,Cr
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	double G1,G2, yctemp, ycrtemp, xtemp, h, dh;
	double dGm2aC,dGm2aVA,dGm2aFe,dGm2aCr,dGm2bC,dGm2bVA,dGm2bFe,dGm2bCr;
	double ddGm[2], dGma[2], dGmb[2];
	double XdGm_lib, XdGmCc_lib, XdGmCc_C, XdGmCc_CR, XdGmCcr_lib, XdGmCcr_CR, XdGmCcr_C;
	double part1_gam, part2_gam, part1_alf, part2_alf, \
	part11_gam, part22_gam, part11_alf, part22_alf, \
	XdGm_C, XdGm_CR;
	
	double ddGm2aCC, ddGm2aCVA, ddGm2aVAC, ddGm2aVAVA;
	double ddGm2aFeFe, ddGm2aFeCr, ddGm2aCrFe, ddGm2aCrCr;
	double ddGm2aFeC, ddGm2aFeVA, ddGm2aCrC, ddGm2aCrVA;
	
	double ddGm2bCC, ddGm2bCVA, ddGm2bVAC, ddGm2bVAVA;
	double ddGm2bFeFe, ddGm2bFeCr, ddGm2bCrFe, ddGm2bCrCr;
	double ddGm2bFeC, ddGm2bFeVA, ddGm2bCrC, ddGm2bCrVA;
	double Lc, Lcrcr;
	double Xa[2],Xb[2];
	
	double Energy_gibbs, potchem_fe, potchem_cr, potchem_c, potchem_cr_fe;
	int i,j;
	std::vector < std::vector<double> > D(12);
	for(i=0;i<12;i++)
	{
		D[i].resize(4);
	}
	
	for(i=1; i<m+1; i++)
	{
		x += h1;
		//cout<<"x"<<x<<endl;
		
		Y1alfa[1]  = ycr[i]; //CR
		Y1alfa[0]  = 1.-Y1alfa[1]; //FE
		Y2alfa[0]  = yc[i]/ph1->Thermo->m[1]; //C
		Y2alfa[1]  = 1.-Y2alfa[0]; //VA
		Y1betta[1] = ycr[i]; //CR
		Y1betta[0] = 1.-Y1betta[1]; //FE
		Y2betta[0] = yc[i]/ph2->Thermo->m[1]; //C
		Y2betta[1] = 1.-Y2betta[0]; //VA
		
		//Calcul energy Gibbs two phases (ferrite and austenite for system Fe-C-Va)
		G1 = (ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) \
		*(ph1->Thermo->functionGibbs(TT2,Y1alfa,Y2alfa));	//ferrite
		//cout<<"G1="<<G1<<endl;
		G2 = (ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) \
		*(ph2->Thermo->functionGibbs(TT2,Y1betta,Y2betta));	//austenite
		//cout<<"G2="<<G2<<endl;
		
		//Derivation 1 order and 2 order y
		double yc_=yc[i];
		double ycr_=ycr[i];
		D=derivFeCrC_U(yc_,ycr_);
		
		
		h=Wang(x);
		dh=dWang(x);
		
		//Derivation energy gibbs 1 order
		ph1->Thermo->function_derivation_Gibbs(TT2,Y1alfa,Y2alfa);
		dGm2aC	=ph1->Thermo->dGm2[0];	//C dans ferrite
		//cout<<"dGm2aC="<<dGm2aC<<endl;
		dGm2aVA	=ph1->Thermo->dGm2[1];	//VA dans ferrite
		//cout<<"dGm2aVA="<<dGm2aVA<<endl;
		dGm2aFe	=ph1->Thermo->dGm1[0];	//Fe dans ferrite
		//cout<<"dGm2aFe="<<dGm2aFe<<endl;
		dGm2aCr	=ph1->Thermo->dGm1[1];	//Cr dans ferrite
		//cout<<"dGm2aCr="<<dGm2aCr<<endl;
		ph2->Thermo->function_derivation_Gibbs(TT2,Y1betta,Y2betta);
		dGm2bC	=ph2->Thermo->dGm2[0];	//C dans austenite
		//cout<<"dGm2bC="<<dGm2bC<<endl;
		dGm2bVA	=ph2->Thermo->dGm2[1];	//VA dans austenite
		//cout<<"dGm2bVA="<<dGm2bVA<<endl;
		dGm2bFe	=ph2->Thermo->dGm1[0];	//Fe dans austenite
		//cout<<"dGm2bFe="<<dGm2bFe<<endl;
		dGm2bCr	=ph2->Thermo->dGm1[1];	//Cr dans austenite
		//cout<<"dGm2bCr="<<dGm2bCr<<endl;
		
		//Energy_gibbs	= 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0])*G2*h \
		+1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0])*G1*(1.-h);
		
		
		dGmb[0]=dGm2bC*D[2][2]+dGm2bVA*D[3][2]+dGm2bFe*D[0][2]+dGm2bCr*D[1][2]; //*dy/dcc (austenite)
		dGma[0]=dGm2aC*D[2][0]+dGm2aVA*D[3][0]+dGm2aFe*D[0][0]+dGm2aCr*D[1][0]; //*dy/dcc (ferrite)
		
		dGmb[1]=dGm2bC*D[2][3]+dGm2bVA*D[3][3]+dGm2bFe*D[0][3]+dGm2bCr*D[1][3]; //*dy/dccr (austenite)
		dGma[1]=dGm2aC*D[2][1]+dGm2aVA*D[3][1]+dGm2aFe*D[0][1]+dGm2aCr*D[1][1]; //*dy/dccr (ferrite)
		
		ddGm[0]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
				   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
		+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
				 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0])*dGma[0] ) );/*dGm/dcc*/
		
		ddGm[1]=h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
				   +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1] ) \
		+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
				 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);/*dGm/dccr*/
		
		
		XdGm_lib	= ( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*G2 \
		-( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*G1;
		
		XdGm_C		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
						 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]) \
		+(1.-h)*(- ( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
				 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0]);
		
		XdGm_CR		= h*(-( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
						 +( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]) \
		+(1.-h)*(-( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
				 +( 1./(ph1->Thermo->m[0] + ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1]);	
		
		//Derivation energy gibbs 2 order
		ph1->Thermo->function_derivation_Gibbs2(TT2,Y1alfa,Y2alfa);
		ddGm2aCC	=ph1->Thermo->dG_3(0,0);	//C:C dans ferrite
		ddGm2aCVA	=ph1->Thermo->dG_3(0,1);	//C:VA dans ferrite
		ddGm2aVAC	=ph1->Thermo->dG_3(1,0);	//VA:C dans ferrite
		ddGm2aVAVA	=ph1->Thermo->dG_3(1,1);	//VA:VA dans ferrite
		
		
		
		ddGm2aFeFe	=ph1->Thermo->dG_1(0,0);	//FE:FE dans ferrite
		ddGm2aFeCr	=ph1->Thermo->dG_1(0,1);	//FE:CR dans ferrite
		ddGm2aCrFe	=ph1->Thermo->dG_1(1,0);	//CR:FE dans ferrite
		ddGm2aCrCr	=ph1->Thermo->dG_1(1,1);	//CR:CR dans ferrite
		
		ddGm2aFeC	=ph1->Thermo->dG_2(0,0);	//FE:C dans ferrite
		ddGm2aFeVA	=ph1->Thermo->dG_2(0,1);	//FE:VA dans ferrite
		ddGm2aCrC	=ph1->Thermo->dG_2(1,0);	//CR:C dans ferrite
		ddGm2aCrVA	=ph1->Thermo->dG_2(1,1);	//CR:VA dans ferrite
		
		
		ph2->Thermo->function_derivation_Gibbs2(TT2,Y1betta,Y2betta);
		ddGm2bCC	=ph2->Thermo->dG_3(0,0);	//C:C dans austenite
		ddGm2bCVA	=ph2->Thermo->dG_3(0,1);	//C:VA dans austenite
		ddGm2bVAC	=ph2->Thermo->dG_3(1,0);	//VA:C dans austenite
		ddGm2bVAVA	=ph2->Thermo->dG_3(1,1);	//VA:VA dans austenite
		
		ddGm2bFeFe	=ph2->Thermo->dG_1(0,0);	//FE:FE dans austenite
		ddGm2bFeCr	=ph2->Thermo->dG_1(0,1);	//FE:CR dans austenite
		ddGm2bCrFe	=ph2->Thermo->dG_1(1,0);	//CR:FE dans austenite
		ddGm2bCrCr	=ph2->Thermo->dG_1(1,1);	//CR:CR dans austenite
		
		ddGm2bFeC	=ph2->Thermo->dG_2(0,0);	//FE:C dans austenite
		ddGm2bFeVA	=ph2->Thermo->dG_2(0,1);	//FE:VA dans austenite
		ddGm2bCrC	=ph2->Thermo->dG_2(1,0);	//CR:C dans austenite
		ddGm2bCrVA	=ph2->Thermo->dG_2(1,1);	//CR:VA dans austenite
		
		
		
		
		
		XdGmCc_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*G2 \
		+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[0]	\
		+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*G1 \
		-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[0];
		
		part11_gam	=ddGm2bFeFe*pow(D[0][2],2.)+ddGm2bCrCr*pow(D[1][2],2.) \
		+ddGm2bCC*pow(D[2][2],2.)+ddGm2bVAVA*pow(D[3][2],2.)	\
		+dGm2bFe*D[4][2]+dGm2bCr*D[5][2]+dGm2bC*D[6][2]+dGm2bVA*D[7][2]	\
		+D[0][2]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
		+D[1][2]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
		+D[2][2]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2]) \
		+D[3][2]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
		
		
		
		part22_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
		+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
		+dGm2bFe*D[8][2]+dGm2bCr*D[9][2]+dGm2bC*D[10][2]+dGm2bVA*D[11][2] \
		+D[0][2]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
		+D[1][2]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
		+D[2][2]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3])	\
		+D[3][2]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
		
		
		part11_alf	=ddGm2aFeFe*pow(D[0][0],2.)+ddGm2aCrCr*pow(D[1][0],2.) \
		+ddGm2aCC*pow(D[2][0],2.)+ddGm2aVAVA*pow(D[3][0],2.)	\
		+dGm2aFe*D[4][0]+dGm2aCr*D[5][0]+dGm2aC*D[6][0]+dGm2aVA*D[7][0]	\
		+D[0][0]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
		+D[1][0]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
		+D[2][0]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0]) \
		+D[3][0]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
		
		part22_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
		+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
		+dGm2aFe*D[8][0]+dGm2aCr*D[9][0]+dGm2aC*D[10][0]+dGm2aVA*D[11][0] \
		+D[0][0]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
		+D[1][0]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
		+D[2][0]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1])	\
		+D[3][0]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
		
		XdGmCc_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][2],2.)*G2 \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][2]*G2 \
						 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[0]
						 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part11_gam)\
		+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][0],2.)*G1 \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][0]*G1 \
				 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[0] \
				 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part11_alf);
		
		
		
		XdGmCc_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][2]*G2 \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
						 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part22_gam) \
		+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][0]*G1 \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
				 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part22_alf);
		
		
		part1_gam	=ddGm2bFeFe*pow(D[0][3],2.)+ddGm2bCrCr*pow(D[1][3],2.)+ddGm2bCC*pow(D[2][3],2.)+ddGm2bVAVA*pow(D[3][3],2.) \
		+dGm2bFe*D[4][3]+dGm2bCr*D[5][3]+dGm2bC*D[6][3]+dGm2bVA*D[7][3] \
		+D[0][3]*(ddGm2bFeCr*D[1][3]+ddGm2bFeC*D[2][3]+ddGm2bFeVA*D[3][3]) \
		+D[1][3]*(ddGm2bCrFe*D[0][3]+ddGm2bCrC*D[2][3]+ddGm2bCrVA*D[3][3]) \
		+D[2][3]*(ddGm2bFeC*D[0][3]+ddGm2bCrC*D[1][3]+ddGm2bCVA*D[3][3]) \
		+D[3][3]*(ddGm2bFeVA*D[0][3]+ddGm2bCrVA*D[1][3]+ddGm2bVAC*D[2][3]);
		
		part2_gam	= ddGm2bFeFe*D[0][3]*D[0][2]+ddGm2bCrCr*D[1][3]*D[1][2] \
		+ddGm2bCC*D[2][2]*D[2][3]+ddGm2bVAVA*D[3][3]*D[3][2] \
		+dGm2bFe*D[8][3]+dGm2bCr*D[9][3]+dGm2bC*D[10][3]+dGm2bVA*D[11][3] \
		+D[0][3]*(ddGm2bFeCr*D[1][2]+ddGm2bFeC*D[2][2]+ddGm2bFeVA*D[3][2]) \
		+D[1][3]*(ddGm2bCrFe*D[0][2]+ddGm2bCrC*D[2][2]+ddGm2bCrVA*D[3][2]) \
		+D[2][3]*(ddGm2bFeC*D[0][2]+ddGm2bCrC*D[1][2]+ddGm2bCVA*D[3][2])	\
		+D[3][3]*(ddGm2bFeVA*D[0][2]+ddGm2bCrVA*D[1][2]+ddGm2bVAC*D[2][2]);
		
		part1_alf	=ddGm2aFeFe*pow(D[0][1],2.)+ddGm2aCrCr*pow(D[1][1],2.)+ddGm2aCC*pow(D[2][1],2.)+ddGm2aVAVA*pow(D[3][1],2.) \
		+dGm2aFe*D[4][1]+dGm2aCr*D[5][1]+dGm2aC*D[6][1]+dGm2aVA*D[7][1] \
		+D[0][1]*(ddGm2aFeCr*D[1][1]+ddGm2aFeC*D[2][1]+ddGm2aFeVA*D[3][1]) \
		+D[1][1]*(ddGm2aCrFe*D[0][1]+ddGm2aCrC*D[2][1]+ddGm2aCrVA*D[3][1]) \
		+D[2][1]*(ddGm2aFeC*D[0][1]+ddGm2aCrC*D[1][1]+ddGm2aCVA*D[3][1]) \
		+D[3][1]*(ddGm2aFeVA*D[0][1]+ddGm2aCrVA*D[1][1]+ddGm2aVAC*D[2][1]);
		
		part2_alf	= ddGm2aFeFe*D[0][1]*D[0][0]+ddGm2aCrCr*D[1][1]*D[1][0] \
		+ddGm2aCC*D[2][0]*D[2][1]+ddGm2aVAVA*D[3][1]*D[3][0] \
		+dGm2aFe*D[8][1]+dGm2aCr*D[9][1]+dGm2aC*D[10][1]+dGm2aVA*D[11][1] \
		+D[0][1]*(ddGm2aFeCr*D[1][0]+ddGm2aFeC*D[2][0]+ddGm2aFeVA*D[3][0]) \
		+D[1][1]*(ddGm2aCrFe*D[0][0]+ddGm2aCrC*D[2][0]+ddGm2aCrVA*D[3][0]) \
		+D[2][1]*(ddGm2aFeC*D[0][0]+ddGm2aCrC*D[1][0]+ddGm2aCVA*D[3][0])	\
		+D[3][1]*(ddGm2aFeVA*D[0][0]+ddGm2aCrVA*D[1][0]+ddGm2aVAC*D[2][0]);
		
		XdGmCcr_lib	= -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*G2 \
		+( 1./(ph2->Thermo->m[0] + ph2->Thermo->m[1]*Y2betta[0]) )*dGmb[1]	\
		+( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*G1 \
		-( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*dGma[1];
		
		XdGmCcr_CR	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*pow(D[2][3],2.)*G2 \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[6][3]*G2 \
						 -2.*( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[1] \
						 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part1_gam) \
		+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*pow(D[2][1],2.)*G1 \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[6][1]*G1 \
				 -2.*( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[1] \
				 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part1_alf);
		
		XdGmCcr_C	= h*(2.*( pow(ph2->Thermo->m[1],2.)/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),3.) )*D[2][2]*D[2][3]*G2 \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[10][3]*G2 \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][3]*dGmb[0] \
						 -( ph2->Thermo->m[1]/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]),2.) )*D[2][2]*dGmb[1] \
						 +( 1./(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]) )*part2_gam) \
		+(1.-h)*(2.*( pow(ph1->Thermo->m[1],2.)/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),3.) )*D[2][0]*D[2][1]*G1 \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[10][1]*G1 \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][1]*dGma[0] \
				 -( ph1->Thermo->m[1]/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]),2.) )*D[2][0]*dGma[1] \
				 +( 1./(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]) )*part2_alf);
		
		//Fraction molaire de C
		Xa[0]=x_c0;//(ph1->Thermo->m[1]*Y2alfa[0]) \
		/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y2alfa[1]);//carbone dans ferrite
		Xa[1]=x_cr0;//(ph1->Thermo->m[0]*Y1alfa[1]) \
		/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y2alfa[1]);//chrome dans ferrite
		//Xa[0]=(ph1->Thermo->m[1]*y[0]) \
		/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*(1.-y[0]));//carbone dans ferrite
		//Xa[1]=(ph1->Thermo->m[0]*y[1]) \
		/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*(1.-y[1]));//chrome dans ferrite
		
		Xb[0]=(ph2->Thermo->m[1]*Y2betta[0]) \
		/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y2betta[1]);//carbone dans austenite
		Xb[1]=(ph2->Thermo->m[0]*Y1betta[1]) \
		/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y2betta[1]);//chrome dans austenite
		
		
		double AA,BB,CC,DD,EE,GG,KK,LL,FF,PP;
		
		//int dim = ph2->Thermo->Dex.size();
		std::vector <double> Diffbetta(2);	
		//Diffbetta = calc_diffusion_betta(concentration);
		//Diffbetta[0] = Diffbetta[0]*pow(1.e6,2.);//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
		//Diffbetta[1] = Diffbetta[1]*pow(1.e6,2.);//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;		
		//Diffbetta[0] = Diffbetta[1]*correc_D;
		//Diffbetta[1] = Diffbetta[1];
		Diffbetta[0] = dcr*correc_D*coeff_n;
		Diffbetta[1] = dcr*coeff_k;

		double mum_c, mum_cr;	
		mum_c	=(1.+yc[i])*(XdGmCc_C+XdGmCc_CR)+2*ddGm[0]+ddGm[1];
		mum_cr	=(1.+yc[i])*(XdGmCcr_CR+XdGmCcr_C)+ddGm[1];
		Lc=Diffbetta[0]/mum_c/pow(taille_gradient,2)/Vm; 
		Lcrcr=Diffbetta[1]/mum_cr/pow(taille_gradient,2)/Vm;
		double VIT_TEMPORAIRE_=VIT_TEMPORAIRE/taille_gradient;
		
		
		
		
		AA	=	Vm*Lc*dh*(XdGm_lib+(1.+yc[i])*XdGmCc_lib);
		//printf("Lc= %e XdGm_lib= %e XdGmCc_lib= %e XdGmCcr_lib= %e dh= %e \n",Lc,XdGm_lib,XdGmCc_lib,XdGmCcr_lib,dh);
		BB	=	Vm*Lc*(XdGm_C+ddGm[0]+(1.+yc[i])*XdGmCc_C);
		//printf("Lc= %e XdGm_CR= %e ddGm[1]= %e XdGmCc_CR= %e XdGmCcr_CR= %e \n",Lc,XdGm_C,ddGm[0],XdGmCc_C,XdGmCcr_C);
		CC	=	Vm*Lc*(XdGm_CR+(1.+yc[i])*XdGmCc_CR);
		//printf("Lc= %e XdGm_CR= %e ddGm[1]= %e XdGmCc_CR= %e XdGmCcr_CR= %e \n",Lc,XdGm_CR,ddGm[1],XdGmCc_CR,XdGmCcr_CR);
		DD	=	Vm*Lcrcr*(dh*(1.+yc[i])*XdGmCcr_lib);
		//printf("Lcrcr= %e  dh = %e  XdGmCcr_lib = %e \n",Lcrcr,dh,XdGmCcr_lib);
		EE	=	Vm*Lcrcr*((1.+yc[i])*XdGmCcr_C+ddGm[1]);
		//printf("Lcrcr= %e  XdGmCcr_C = %e \n",Lcrcr,XdGmCcr_C);
		GG	=	Vm*Lcrcr*(1.+yc[i])*(XdGmCcr_CR);
		//printf("Lcrcr= %e  XdGmCcr_CR = %e \n",Lcrcr,XdGmCcr_CR);
		KK	=	-VIT_TEMPORAIRE_*yc[i]+VIT_TEMPORAIRE_*Xa[0];
		//printf("KK= %e y[0]= %e Xa[0]= %e \n",KK,y[0],Xa[0]);
		LL	=	-VIT_TEMPORAIRE_*yc[i]+VIT_TEMPORAIRE_*Xa[1];
		//printf("LL= %e y[1]= %e Xa[1]= %e \n",LL,y[1],Xa[1]);
		FF	=	((LL-DD)*BB-(KK-AA)*EE)/(GG*BB-CC*EE); 
		//printf ("FF= %e \n",FF);
		PP	=	(KK-AA-FF*CC)/BB; 	//printf("PP= %e \n",PP);
		
		
		
		
		AA	=	dh*(XdGm_lib+(1.+yc[i])*XdGmCc_lib);
		BB	=	XdGm_C+ddGm[0]+(1.+yc[i])*XdGmCc_C;
		CC	=	XdGm_CR+(1.+yc[i])*XdGmCc_CR;
        potchem_c = AA+PP*BB+FF*CC;
		
		DD	=	dh*(1.+yc[i])*XdGmCcr_lib;
		EE	=	(1.+yc[i])*XdGmCcr_C+ddGm[1];
		GG	=   (1.+yc[i])*(XdGmCcr_CR);
		potchem_cr_fe = DD+PP*EE+FF*GG;
		
		

		
		P[i][0]=(1./taille_gradient)*pow(taille_gradient,2)*Vm*Lc*pow(potchem_c,2);
		P[i][1]=(1./taille_gradient)*pow(taille_gradient,2)*Vm*Lcrcr*pow(potchem_cr_fe,2);
		P[i][2]=pow(taille_gradient,2)*Lc*Vm;
		P[i][3]=pow(taille_gradient,2)*Lcrcr*Vm;
    }		
	
}



//___________________kinetic transition Fe-C_____________________
void interface::KTB(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(1); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(1);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	std::vector <double> Y1(7);
	
	Y1alfa[0]  =  Y[0]; //bcc_A2: FE
	Y2alfa[0]  =  Y[1]; //bcc_A2: C
	Y2alfa[1]  =  Y[2]; //bcc_A2: VAte
	Y1betta[0] =  Y[3]; //fcc_A1: FE
	Y2betta[0] =  Y[4]; //fcc_A1: C 
	Y2betta[1] =  Y[5]; //fcc_A1: VA
	Y1=Y;
	
	std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
	std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"
	
	int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(dim);	
	//Diffbetta = calc_diffusion_betta(concentration);
	//Diffbetta[0] = Diffbetta[0]*pow(1.e9,2.);
	//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	std::vector <double> delta(1);
	delta = Zenerplane2d(Y,concentration);
	
	//vector <double> res(1);
	//res=integration(n,Y1,concentration);
	//cout<<"res="<<res[0]<<endl;
	//double res=0.0;
	
	pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
	pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
	F[0] = (pot_chemicalbetta[0]-pot_chemicalalfa[0])/(R*TT2); //C_bcc - C_fcc
	F[1] = (pot_chemicalbetta[1]-pot_chemicalalfa[1])/(R*TT2); //Fe_bcc - Fe_fcc
	F[2] = Y[0] - 1.;        //Fe-1.
	F[3] = Y[3] - 1.;        //Fe-1.
	F[4] = Y[1] + Y[2] - 1.; //C+VA-1.
	F[5] = Y[4] + Y[5] - 1.; //C+VA-1.

	
	double x_c_alfa=(ph1->Thermo->m[1]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y[1]);
	double x_c_gamma=(ph2->Thermo->m[1]*Y[4])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y[4]);
	
	
	F[6] = (x_c_gamma-x_c_alfa)*Y[6]-(correc_D/(Rzin*delta[0]))*(x_c_gamma-concentration[0]);
	
	//-((ph1->Thermo->m[1]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[5]))*Y[6] \
	+((ph2->Thermo->m[1]*Y[4])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[5]))*(Y[6]-(correc_D/(Rzin*delta[0]))) \
	+(correc_D/(Rzin*delta[0]))*concentration[0];//C (gamma->alpha)
	//cout<<"correc_D="<<correc_D<<endl;
}



//___________________kinetic transition Fe-C_____________________
void interface::KTB_FeNi(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
double TT2 = 0.;
TT2 = ph1->Thermo->T1;
std::vector <double> Y1alfa(2); //bcc_A2: FE
std::vector <double> Y2alfa(2); //bcc_A2: C,VA
std::vector <double> Y1betta(2);//fcc_A1: FE
std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
Y1alfa[0]  =  Y[0]; //bcc_A2: FE
Y1alfa[1]  =  Y[1]; //bcc: Ni
Y2alfa[0]  =  0.0; //bcc_A2: C
Y2alfa[1]  =  1.; //bcc_A2: VA
Y1betta[0] =  Y[2]; //fcc_A1: FE
Y1betta[1] =  Y[3]; //fcc: Ni
Y2betta[0] =  0.0; //fcc_A1: C 
Y2betta[1] =  1.0; //fcc_A1: VA
	
std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"
	
std::vector <double> Diffbetta(1);	
Diffbetta[0] = dcr;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	
std::vector <double> delta(1);
std::vector <double> omega(1);

omega[0] = (Y[3]-concentration[0])/(Y[3]-Y[1]);
delta[0] = 	2.*(1.-omega[0])/omega[0];
	
pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
F[0] = (pot_chemicalbetta[2]-pot_chemicalbetta[1])/(R*TT2) \
		-(pot_chemicalalfa[2]-pot_chemicalalfa[1])/(R*TT2); //Ni_bcc - Ni_fcc
F[1] = Y[0]*(pot_chemicalbetta[1]-pot_chemicalalfa[1])/(R*TT2) \
		+Y[1]*(pot_chemicalbetta[2]-pot_chemicalalfa[2]); 
F[2] = Y[0] + Y[1] - 1.;        //Fe-1.
F[3] = Y[2] + Y[3] - 1.;        //Fe-1.

double x_c_alfa=(ph1->Thermo->m[1]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y[1]);
double x_c_gamma=(ph2->Thermo->m[1]*Y[4])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y[4]);
	
F[4] = (Y[3]-Y[1])*Y[4]-(dcr/(Rzin*delta[0]))*(Y[3]-concentration[0]);
	
}



//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_FeNi(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_FeNi(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_FeNi(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	



//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	

void interface::KTB_prof(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(1); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(1);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	//std::vector <double> Y1(7);
	
	Y1alfa[0]  =  1.0; //bcc_A2: FE
	Y2alfa[0]  =  Y[0]; //bcc_A2: C
	Y2alfa[1]  =  1.-Y2alfa[0]; //bcc_A2: VAte
	Y1betta[0] =  1.0; //fcc_A1: FE
	Y2betta[0] =  Y[1]; //fcc_A1: C 
	Y2betta[1] =  1.-Y2betta[0]; //fcc_A1: VA

	
	std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
	std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"


	
	pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
	pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1


	double x_c_alpha	= ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);						
	double x_fe_alpha	= 1.- x_c_alpha;
	double x_c_gamma	= ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	= 1.-x_c_gamma;	
	
	double Uc_alpha		= x_c_alpha/(1.-x_c_alpha);
	double Ufe_alpha	= x_fe_alpha/(1.-x_c_alpha);
	double Uc_gamma		= x_c_gamma/(1.-x_c_gamma);
	double Ufe_gamma	= x_fe_gamma/(1.-x_c_gamma);
	double U0			= concentration[0]/(1.-concentration[0]);

	double Om=(Uc_gamma-U0)/(Uc_gamma-Uc_alpha);
	double delta1 = 2.*(1.-Om)/Om;
	if(delta1<1.e-16)delta1=1.e-15;


	double L=deriv_mu(Uc_gamma,1.);
	//L = 1e-10;
	//L = 1e-10;
	//L = 1.e-13;
	L = 1e-4;
	//L = 1e-8;
	F[0] = (pot_chemicalbetta[0]-pot_chemicalalfa[0])/(R*TT2)+ \
	((Y[2]/L)*taille_gradient*(Uc_gamma-Uc_alpha))/(R*TT2)*1.e3;
	 //C_bcc - C_fcc
	
	F[1] = Ufe_alpha*(pot_chemicalbetta[1]-pot_chemicalalfa[1])/(R*TT2)+ \
	Uc_alpha*(pot_chemicalbetta[0]-pot_chemicalalfa[0])/(R*TT2)- \
	( (Y[2]*taille_gradient)/L)*pow((Uc_gamma-Uc_alpha),2)/(R*TT2)*1.e3-Y[2]*(Vm/M)/(R*TT2);

	gtr=( (Y[2]*taille_gradient)/L)*pow((Uc_gamma-Uc_alpha),2)/(R*TT2);
	muc=((Y[2]/(L))*taille_gradient*(Uc_gamma-Uc_alpha))/(R*TT2);
	
	printf("L = %e gtr = %e   muc =%e fric = %e \n",L,gtr,muc,Y[2]*(Vm/M)/(R*TT2));

	F[2] = (Uc_gamma-Uc_alpha)*Y[2]-(correc_D/(Rzin*delta1))*(Uc_gamma-U0);
	
	printf("F[0]=%e F[1] = %e F[2]=%e \n",F[0],F[1],F[2]);
}


void interface::KTB_prof_FeNigrowth(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
double TT2 = 0.;
TT2 = ph1->Thermo->T1;
std::vector <double> Y1alfa(2); //bcc_A2: FE
std::vector <double> Y2alfa(2); //bcc_A2: C,VA
std::vector <double> Y1betta(2);//fcc_A1: FE
std::vector <double> Y2betta(2);//fcc_A1: C,VA

	
Y1alfa[0]  =  1.-Y[0]; //bcc_A2: FE
Y1alfa[1]  =  Y[0]; //bcc : NI
Y2alfa[0]  =  0.0; //bcc_A2: C
Y2alfa[1]  =  1.0; //bcc_A2: VAte
Y1betta[0] =  1.-Y[1]; //fcc_A1: FE
Y1betta[1] =  Y[1]; //fcc : NI
Y2betta[0] =  0.0; //fcc_A1: C 
Y2betta[1] =  1.0; //fcc_A1: VA
	
std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"
	
pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	

double Omega=(Y[1]-concentration[0])/(Y[1]-Y[0]);
double delta = 2.*(1.-Omega)/Omega;
if(delta<1.e-16)delta=1.e-15;
	
double L=deriv_mu_FeNi(Y[1],1.);
L=1.e-13;
//printf("L = %e \n",L);
	
F[0] = (pot_chemicalbetta[2] - pot_chemicalbetta[1])/(R*TT2) \
	- (pot_chemicalalfa[2] - pot_chemicalalfa[1])/(R*TT2)+ \
	((Y[2]/(L))*taille_gradient*(Y[1]-Y[0]))/(R*TT2)*8.;
	//Ni_bcc - Ni_fcc


gtr= ( (Y[2]*taille_gradient)/L)*pow((Y[1]-Y[0]),2)/(R*TT2) ;	
//printf("gtr=%e L = %e \n",gtr,L);	
muc = ((Y[2]/(L))*taille_gradient*(Y[1]-Y[0]))/(R*TT2); ///(R*TT2); VRAI
//printf("muni=%e \n",muc);
	
F[1] = (1.-Y[0])*((pot_chemicalbetta[1] - pot_chemicalalfa[1])/(R*TT2))+ \
	Y[0]*((pot_chemicalbetta[2]-pot_chemicalalfa[2])/(R*TT2))- \
	( (Y[2]*taille_gradient)/L)*pow((Y[1]-Y[0]),2)/(R*TT2)*8.-Y[2]*(Vm/M)/(R*TT2)*1e6;	
	 
F[2] = (Y[1]-Y[0])*Y[2]-(dcr/(Rzin*delta))*(Y[1]-concentration[0]);	
}


void interface::KTB_prof_FeNi(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	
	Y1alfa[0]  =  1.-Y[0]; //bcc_A2: FE
	Y1alfa[1]  =  Y[0]; //bcc : NI
	Y2alfa[0]  =  0.0; //bcc_A2: C
	Y2alfa[1]  =  1.0; //bcc_A2: VAte
	Y1betta[0] =  1.-Y[1]; //fcc_A1: FE
	Y1betta[1] =  Y[1]; //fcc : NI
	Y2betta[0] =  0.0; //fcc_A1: C 
	Y2betta[1] =  1.0; //fcc_A1: VA
	
	std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
	std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"
	
	pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
	pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
	
	double Omega=(Y[1]-concentration[0])/(Y[1]-Y[0]);
	double delta = 2.*(1.-Omega)/Omega;
	if(delta<1.e-16)delta=1.e-15;
	
	double L=deriv_mu_FeNi(Y[1],1.);
	L=1e-13;
	///L=1e-15;	

	//printf("L = %e \n",L);
	
	F[0] = (pot_chemicalbetta[2] - pot_chemicalbetta[1])/(R*TT2) \
	- (pot_chemicalalfa[2] - pot_chemicalalfa[1])/(R*TT2)+ \
	((VIT_TEMPORAIRE/(L))*taille_gradient*(Y[1]-Y[0]))/(R*TT2);
	//Ni_bcc - Ni_fcc
	
	
	gtr= ( (VIT_TEMPORAIRE*taille_gradient)/L)*pow((Y[1]-Y[0]),2)/(R*TT2) ;	
	//printf("gtr=%e L = %e \n",gtr,L);	
	muc = ((VIT_TEMPORAIRE/(L))*taille_gradient*(Y[1]-Y[0]))/(R*TT2); ///(R*TT2); VRAI
	//printf("muni=%e \n",muc);
	
	F[1] = (1.-Y[0])*((pot_chemicalbetta[1] - pot_chemicalalfa[1])/(R*TT2))+ \
	Y[0]*((pot_chemicalbetta[2]-pot_chemicalalfa[2])/(R*TT2))- \
	( (VIT_TEMPORAIRE*taille_gradient)/L)*pow((Y[1]-Y[0]),2)/(R*TT2)-VIT_TEMPORAIRE*(Vm/M)/(R*TT2);	
	
	//F[2] = (Y[1]-Y[0])*Y[2]-(dcr/(Rzin*delta))*(Y[1]-concentration[0]);		
}

//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTB_profFeNi0(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);

	
	Y1alfa[0]  =  1.-y_fraction_x_bcc; //bcc_A2: FE
	Y1alfa[1]  =  y_fraction_x_bcc; //bcc : NI
	Y2alfa[0]  =  0.0; //bcc_A2: C
	Y2alfa[1]  =  1.0; //bcc_A2: VAte
	Y1betta[0] =  1.-y_fraction_x_fcc; //fcc_A1: FE
	Y1betta[1] =  y_fraction_x_fcc; //fcc : NI
	Y2betta[0] =  0.0; //fcc_A1: C 
	Y2betta[1] =  1.0; //fcc_A1: VA
	
	
	double Omega=(y_fraction_x_fcc-concentration[0])/(y_fraction_x_fcc-y_fraction_x_bcc);
	double delta = 2.*(1.-Omega)/Omega;
	if(delta<1.e-16)delta=1.e-15;
	
	//if(delta[0]<=0.0) delta[0]=abs(delta[0]);
	//if(delta[0]<1.e-16) delta[0]=1.e-15;
	//if(omega[0]<=0.0) {delta[0]=abs(delta[0]);printf("sursaturation <=0.0 FK2\n");}
	//delta[0] = 1.e-2+1./pow(cosh(omega[0]/0.3),2.);
	
	
	//printf("delta=%e omega_x=%e \n",delta[0],omega[0]);
	//printf("y=%e \n",Y[0]);
	F[0] = (y_fraction_x_fcc-y_fraction_x_bcc)*VIT_TEMPORAIRE*Rzin*delta-dcr*(y_fraction_x_fcc-Y[0]); //CR (gamma ->alpha)
	//printf("F[0]=%e Ucr_fcc=%e Y[0]=%e \n",F[0],Uc_gamma,Y[0]);
	//printf("term=%e \n",(Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE*Rzin*delta[0]);
	
}


//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTB_profFeC0(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(1);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(1);
	std::vector <double> Y2betta(2);
	
	
	Y1alfa[0]  =  1.; //bcc_A2: FE
	Y2alfa[0]  =  y_fraction_c_bcc; //bcc_A2: C
	Y2alfa[1]  =  1.0-y_fraction_c_bcc; //bcc_A2: VAte
	Y1betta[0] =  1.; //fcc_A1: FE
	Y2betta[0] =  y_fraction_c_fcc; //fcc_A1: C 
	Y2betta[1] =  1.0-y_fraction_c_fcc; //fcc_A1: VA
	
	
	double Omega=(y_fraction_c_fcc-concentration[0])/(y_fraction_c_fcc-y_fraction_c_bcc);
	double delta = 2.*(1.-Omega)/Omega;
	if(delta<1.e-16)delta=1.e-15;
	
	//if(delta[0]<=0.0) delta[0]=abs(delta[0]);
	//if(delta[0]<1.e-16) delta[0]=1.e-15;
	//if(omega[0]<=0.0) {delta[0]=abs(delta[0]);printf("sursaturation <=0.0 FK2\n");}
	//delta[0] = 1.e-2+1./pow(cosh(omega[0]/0.3),2.);
	
	
	//printf("delta=%e omega_x=%e \n",delta[0],omega[0]);
	//printf("y=%e \n",Y[0]);
	F[0] = (y_fraction_c_fcc-y_fraction_c_bcc)*VIT_TEMPORAIRE*Rzin*delta-dc*(y_fraction_c_fcc-Y[0]); //CR (gamma ->alpha)
	//printf("F[0]=%e Ucr_fcc=%e Y[0]=%e \n",F[0],Uc_gamma,Y[0]);
	//printf("term=%e \n",(Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE*Rzin*delta[0]);
	
}





void interface::KTB_prof_FeC_V(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
double TT2 = 0.;
TT2 = ph1->Thermo->T1;
std::vector <double> Y1alfa(1); //bcc_A2: FE
std::vector <double> Y2alfa(2); //bcc_A2: C,VA
std::vector <double> Y1betta(1);//fcc_A1: FE
std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
Y1alfa[0]  =  1.0; //bcc_A2: FE
Y2alfa[0]  =  Y[0]; //bcc_A2: C
Y2alfa[1]  =  1.-Y2alfa[0]; //bcc_A2: VAte
Y1betta[0] =  1.0; //fcc_A1: FE
Y2betta[0] =  Y[1]; //fcc_A1: C 
Y2betta[1] =  1.-Y2betta[0]; //fcc_A1: VA
	
	
std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"
	
pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
double x_c_alpha	= ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);						
double x_fe_alpha	= 1.- x_c_alpha;
double x_c_gamma	= ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
double x_fe_gamma	= 1.-x_c_gamma;	
	
double Uc_alpha		= x_c_alpha/(1.-x_c_alpha);
double Ufe_alpha	= x_fe_alpha/(1.-x_c_alpha);
double Uc_gamma		= x_c_gamma/(1.-x_c_gamma);
double Ufe_gamma	= x_fe_gamma/(1.-x_c_gamma);
double U0			= concentration[0]/(1.-concentration[0]);
	
double Om=(Uc_gamma-U0)/(Uc_gamma-Uc_alpha);
double delta1 = 2.*(1.-Om)/Om;
if(delta1<1.e-16)delta1=1.e-15;
	


double L=deriv_mu(Uc_gamma,1.);
L = 1e-12;
	
F[0] = (pot_chemicalbetta[0]-pot_chemicalalfa[0])/(R*TT2)+ \
	((VIT_TEMPORAIRE/L)*taille_gradient*(Uc_gamma-Uc_alpha))/(R*TT2)*1e3;//1e2
	//C_bcc - C_fcc
	
F[1] = (Ufe_alpha*(pot_chemicalbetta[1]-pot_chemicalalfa[1])/(R*TT2)+ \
	Uc_alpha*(pot_chemicalbetta[0]-pot_chemicalalfa[0])/(R*TT2))- \
	( (VIT_TEMPORAIRE*taille_gradient)/L)*pow((Uc_gamma-Uc_alpha),2)/(R*TT2)*1e4-VIT_TEMPORAIRE*(Vm/M)/(R*TT2)*1e6;
	
gtr = Ufe_alpha*(pot_chemicalbetta[1]-pot_chemicalalfa[1])/(R*TT2)+ \
	Uc_alpha*(pot_chemicalbetta[0]-pot_chemicalalfa[0])/(R*TT2);

gtr1 = 	( (VIT_TEMPORAIRE*taille_gradient)/L)*pow((Uc_gamma-Uc_alpha),2)/(R*TT2)*1e4-VIT_TEMPORAIRE*(Vm/M)/(R*TT2)*1e6;

	
muc=-((VIT_TEMPORAIRE/L)*taille_gradient*(Uc_gamma-Uc_alpha))/(R*TT2)*1e3;//1e3
	
//F[1] = (Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE-(correc_D/(Rzin*delta1))*(Uc_gamma-U0);	
	
}



void interface::KTB_prof_FeNi_V(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	
	Y1alfa[0]  =  1.-Y[0]; //bcc_A2: FE
	Y1alfa[1]  =  Y[0]; //bcc : NI
	Y2alfa[0]  =  0.0; //bcc_A2: C
	Y2alfa[1]  =  1.0; //bcc_A2: VAte
	Y1betta[0] =  1.-Y[1]; //fcc_A1: FE
	Y1betta[1] =  Y[1]; //fcc : NI
	Y2betta[0] =  0.0; //fcc_A1: C 
	Y2betta[1] =  1.0; //fcc_A1: VA
	
	std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
	std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"
	
	pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
	pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
	
double Omega=(Y[1]-concentration[0])/(Y[1]-Y[0]);
double delta = 2.*(1.-Omega)/Omega;
if(delta<1.e-16)delta=1.e-15;
	
double L=deriv_mu_FeNi(Y[1],1.);
L=1.e-17;	

F[0] = (pot_chemicalbetta[2] - pot_chemicalbetta[1])/(R*TT2) \
- (pot_chemicalalfa[2] - pot_chemicalalfa[1])/(R*TT2)+ \
((VIT_TEMPORAIRE/(L))*taille_gradient*(Y[1]-Y[0]))/(R*TT2)*2.;
//Ni_bcc - Ni_fcc
//F[0] = (pot_chemicalbetta[2] - pot_chemicalbetta[1])/(R*TT2) \
	- (pot_chemicalalfa[2] - pot_chemicalalfa[1])/(R*TT2)+ \
	((VIT_TEMPORAIRE/(L))*taille_gradient*(Y[1]-Y[0]))/(R*TT2);
	//Ni_bcc - Ni_fcc
	
muc = -((VIT_TEMPORAIRE/(L))*taille_gradient*(Y[1]-Y[0]))/(R*TT2); ///(R*TT2); VRAI

gtr = ( (VIT_TEMPORAIRE*taille_gradient)/L)*pow((Y[1]-Y[0]),2)/(R*TT2)*2.+VIT_TEMPORAIRE*(Vm/M)/(R*TT2)*1e5;
	
gtr1 = (1.-Y[0])*((pot_chemicalbetta[1] - pot_chemicalalfa[1])/(R*TT2))+ \
	Y[0]*((pot_chemicalbetta[2]-pot_chemicalalfa[2])/(R*TT2));	

F[1] = (1.-Y[0])*((pot_chemicalbetta[1] - pot_chemicalalfa[1])/(R*TT2))+ \
		Y[0]*((pot_chemicalbetta[2]-pot_chemicalalfa[2])/(R*TT2)) \
		-( (VIT_TEMPORAIRE*taille_gradient)/L)*pow((Y[1]-Y[0]),2)/(R*TT2)*2. \
		-VIT_TEMPORAIRE*(Vm/M)/(R*TT2)*1e5;
	
//F[1] = (1.-Y[0])*((pot_chemicalbetta[1] - pot_chemicalalfa[1])/(R*TT2))+ \
	Y[0]*((pot_chemicalbetta[2]-pot_chemicalalfa[2])/(R*TT2)) \
	-( (VIT_TEMPORAIRE*taille_gradient)/L)*pow((Y[1]-Y[0]),2)/(R*TT2)* \
	-VIT_TEMPORAIRE*(Vm/M)/(R*TT2)*1e5;	
//F[1]= (Y[1]-Y[0])*VIT_TEMPORAIRE-(dcr/(Rzin*delta))*(Y[1]-concentration[0]);	
}



//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_prof(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_prof(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_prof(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	


//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_prof_FeNigrowth(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_prof_FeNigrowth(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_prof_FeNigrowth(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	


//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_prof_FeNi0(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_profFeNi0(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_profFeNi0(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	

//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_prof_FeC0(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_profFeC0(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_profFeC0(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	

//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_prof_FeNi(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_prof_FeNi(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_prof_FeNi(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	
//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_prof_FeNi_V(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_prof_FeNi_V(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_prof_FeNi_V(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	

//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB_prof_FeC_V(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian2D.resize(i1);
	
	for (i = 0; i < Jacobian2D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian2D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian2D.size() ; i++)
		for(j=0 ; j < Jacobian2D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTB_prof_FeC_V(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB_prof_FeC_V(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	



//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
double TT2 = 0.;
TT2 = ph1->Thermo->T1;
std::vector <double> Y1alfa(2);
std::vector <double> Y2alfa(2);
std::vector <double> Y1betta(2);
std::vector <double> Y2betta(2);
	
Y1alfa[0]  = Y[0]; //FE
Y1alfa[1]  = Y[1]; //CR
Y2alfa[0]  = Y[2]; //C
Y2alfa[1]  = Y[3]; //VA
Y1betta[0] = Y[4]; //FE
Y1betta[1] = Y[5]; //CR
Y2betta[0] = Y[6]; //C
Y2betta[1] = Y[7]; //VA

std::vector <double> pot_chemical1alfa(3);
std::vector <double> pot_chemical1betta(3);
	
//int dim = ph2->Thermo->Dex.size();
std::vector <double> Diffbetta(2);	
Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;

double x_c_bcc = 3.*Y[2]/(1.+3.*Y[2]);
double x_c_fcc = Y[6]/(1.+Y[6]);
double x_cr_bcc= Y[1]/(1.+3.*Y[2]);
double x_cr_fcc= Y[5]/(1.+Y[6]);
double x_fe_fcc=1.-x_c_fcc-x_cr_fcc;
double x_fe_bcc=1.-x_c_bcc-x_cr_bcc;	
double Uc_bcc	=	x_c_bcc/(1.-x_c_bcc);
double Ucr_bcc	=	x_cr_bcc/(1.-x_c_bcc);
double Ufe_bcc	=	x_fe_bcc/(1.-x_c_bcc);
double Uc_fcc	=	x_c_fcc/(1.-x_c_fcc);
double Ucr_fcc	=	x_cr_fcc/(1.-x_c_fcc);
double Ufe_fcc	=	x_fe_fcc/(1.-x_c_fcc);
double U0c			=	concentration[0]/(1.-concentration[0]);
double U0cr			=	concentration[2]/(1.-concentration[0]);
	
std::vector<double> delta(2);
std::vector<double> omega(2);

omega[0] = ( Uc_fcc - U0c )/(Uc_fcc-Uc_bcc)  ;//c
omega[1] = ( Ucr_fcc - U0cr )/(Ucr_fcc-Ucr_bcc);//cr
	
delta[0] = 2.*(1.-omega[0])/omega[0];
delta[1] = 2.*(1.-omega[1])/omega[1];

for (int i = 0; i < 2; i++)
	if (delta[i]<1.e-16) delta[i]=1.e-15;
pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2); // carbone
F[1] = ((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2)); // iron
F[2] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
		- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2); //cr-fe
F[3] = Y[0] + Y[1] - 1; //printf("%i\n",F[3]);
F[4] = Y[4] + Y[5] - 1; //printf("%i\n",F[4]);
F[5] = Y[2] + Y[3] - 1; //printf("%i\n",F[5]);
F[6] = Y[6] + Y[7] - 1; //printf("%i\n",F[6]);
F[7] = (Uc_fcc-Uc_bcc)*Y[8]-Diffbetta[1]*correc_D/( Rzin*delta[0] )*(Uc_fcc-U0c);//c
F[8] = (Ucr_fcc-Ucr_bcc)*Y[8]-Diffbetta[1]/( Rzin * delta[1] )*(Ucr_fcc-U0cr); //cr
}


//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_sphere(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[0]  = Y[0]; //FE
	Y1alfa[1]  = Y[1]; //CR
	Y2alfa[0]  = Y[2]; //C
	Y2alfa[1]  = Y[3]; //VA
	Y1betta[0] = Y[4]; //FE
	Y1betta[1] = Y[5]; //CR
	Y2betta[0] = Y[6]; //C
	Y2betta[1] = Y[7]; //VA

	
	
	std::vector <double> pot_chemical1alfa(3);
	std::vector <double> pot_chemical1betta(3);
	
	//int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(2);	
	Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	
	
	double x_c_bcc = 3.*Y[2]/(1.+3.*Y[2]);
	double x_c_fcc = Y[6]/(1.+Y[6]);
	double x_cr_bcc= Y[1]/(1.+3.*Y[2]);
	double x_cr_fcc= Y[5]/(1.+Y[6]);
	double x_fe_fcc=1.-x_c_fcc-x_cr_fcc;
	double x_fe_bcc=1.-x_c_bcc-x_cr_bcc;
	
	
	double Uc_bcc	=	x_c_bcc/(1.-x_c_bcc);
	double Ucr_bcc	=	x_cr_bcc/(1.-x_c_bcc);
	double Ufe_bcc	=	x_fe_bcc/(1.-x_c_bcc);
	double Uc_fcc	=	x_c_fcc/(1.-x_c_fcc);
	double Ucr_fcc	=	x_cr_fcc/(1.-x_c_fcc);
	double Ufe_fcc	=	x_fe_fcc/(1.-x_c_fcc);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);

	
	std::vector<double> delta(2),F_omega(2);
	std::vector<double> omega(2);

	omega[0] = ( Uc_fcc - U0c )/(Uc_fcc-Uc_bcc)  ;//c
	omega[1] = ( Ucr_fcc - U0cr )/(Ucr_fcc-Ucr_bcc);//cr
	
	
	//for (int i =0; i < 2; i++)
	//{
	//	F_omega[i]=omega[i]/pow((1.-omega[i]),(1.+1./50.));
	//}
	
	for (int i=0;i<2;i++)
	{	
		delta[i] = pow((1.-omega[i]),(1.+1./50.));
	}	 

	for (int i = 0; i < 2; i++)
		if (delta[i]<1.e-16) delta[i]=1.e-15;
	
	
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2); // carbone
	F[1] = ((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2)); // iron
	F[2] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2); //cr-fe
	F[3] = Y[0] + Y[1] - 1; //printf("%i\n",F[3]);
	F[4] = Y[4] + Y[5] - 1; //printf("%i\n",F[4]);
	F[5] = Y[2] + Y[3] - 1; //printf("%i\n",F[5]);
	F[6] = Y[6] + Y[7] - 1; //printf("%i\n",F[6]);
	

	F[7] = (Uc_fcc-Uc_bcc)*Y[8]-Diffbetta[1]*correc_D/( Rzin*delta[0] )*(Uc_fcc-U0c);//c
	F[8] = (Ucr_fcc-Ucr_bcc)*Y[8]-Diffbetta[1]/( Rzin * delta[1] )*(Ucr_fcc-U0cr); //cr
}



//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
double TT2 = 0.;
TT2 = ph1->Thermo->T1;
std::vector <double> Y1alfa(2);
std::vector <double> Y2alfa(2);
std::vector <double> Y1betta(2);
std::vector <double> Y2betta(2);
	
Y1alfa[1]  = Y[1]; //CR
Y1alfa[0]  = 1.-Y[1]; //FE
Y2alfa[0]  = Y[0]; //C
Y2alfa[1]  = 1.-Y[0]; //VA
Y1betta[1] = Y[3]; //CR
Y1betta[0] = 1.-Y[3]; //FE
Y2betta[0] = Y[2]; //C
Y2betta[1] = 1.-Y[2]; //VA
//printf("Y(C,BCC)=%e Y(C,FCC)=%e \n",Y2alfa[0],Y2betta[0]);	
//printf("Y(Ni,BCC)=%e Y(Ni,FCC)=%e \n",Y1alfa[1],Y1betta[1]);	

	//########
std::vector <double> pot_chemical1alfa(3);
std::vector <double> pot_chemical1betta(3);
std::vector <double> Diffbetta(2);	
Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
//double compteur_fe;
double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
double U0c			=	concentration[0]/(1.-concentration[0]);
double U0cr			=	concentration[2]/(1.-concentration[0]);
//printf("U0c = %e U0Ni = %e \n",U0c,U0cr);	
//#####Muc
vector <double> L(3);
L=deriv_mu_FeCCr(Uc_gamma,Ucr_gamma,1.,concentration);
double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
//	muc = (Y[4]*taille_gradient/Vm)*(num_c/dem_c)/(R*TT2); //FAUX
muc = (Y[4]*taille_gradient/Vm)*(num_c/dem_c); ///(R*TT2); VRAI
	
//####MuCr
double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
double dem_cr=L[1]*L[0]-pow(L[2],2);
//	mucr= (Y[4]*taille_gradient/Vm)*(num_cr/dem_cr)/(R*TT2); //FAUX
mucr= (Y[4]*taille_gradient/Vm)*(num_cr/dem_cr); ///(R*TT2); VRAI
//#####Integrale Gtr
int i,n1;
gtr= Vm / Y[4] / taille_gradient * \
	( ( L[2] * muc + L[1] * mucr) * mucr + ( L[0] * muc + L[2] * mucr) * muc) ;

//gtr= ((Vm/vtest)*(taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
	
//printf("muc = %e mucr = %e gtr = %e \n",muc/(R*TT2),mucr/(R*TT2),gtr/(R*TT2));
	
F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc/(R*TT2); // carbone
F[1] = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2))-gtr/(R*TT2);// iron
F[2] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr/(R*TT2); //cr-fe

std::vector <double> delta(2),omega(2);
omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);
	
//printf("omegac = %e omegani=%e \n",omega[0],omega[1]);
//printf("Ucfcc=%e Ucbcc=%e\n",Uc_gamma,Uc_alpha);
//printf("Unifcc=%e Unibcc=%e\n",Ucr_gamma,Ucr_alpha);
//printf("Uc0=%e Uni0=%e\n",U0c,U0cr);

//if(omega[1]>1.)omega[1]=9.999952e-01;
	
for (i =0; i < 2; i++)
{
delta[i] = 2.*(1. - omega[i])/omega[i] ;
}	 
for (i = 0; i < 2; i++)
	if (delta[i]<1.e-16) delta[i]=1.e-15;
//Ispol'zovat' eto zapis' v sluchae LENP->PE ili LEP->PE
//F[3] = (Uc_gamma-Uc_alpha)*Y[4]-(Diffbetta[1]*correc_D/( Rzin*delta[0]))* \
(Uc_gamma-U0c); //C (gamma->alpha)
F[3] = (Uc_gamma-Uc_alpha)*Y[4]-(Diffbetta[0]/( Rzin*delta[0]))* \
	(Uc_gamma-U0c); //C (gamma->alpha)
F[4] = (Ucr_gamma-Ucr_alpha)*Y[4]-(Diffbetta[1]/(Rzin*delta[1]))* \
	(Ucr_gamma-U0cr); //CR (gamma ->alpha)
//printf("F[0] = %e F[1] = %e F[2] = %e F[3] = %e F[4] = %e delta_c = %e delta_ni=%e dc=%e dni = %e r =%e coeff_c=%e coeff_n=%e Bcoeff=%e Lcoeff = %e \n",F[0], \
	   F[1],F[2],F[3],F[4],delta[0],delta[1],dc,dcr,Rzin,coeff_k,coeff_n,Lcoef,Bcoef);	
//printf("taille =%e Vm= %e  correc_D=%e \n",taille_gradient,Vm,correc_D);	
//printf("pot_chemical1betta[0] = %e pot_chemical1alfa[0] = %e muc/(R*TT2)=%e",pot_chemical1betta[0],pot_chemical1alfa[0], \
				muc/(R*TT2));	
//exit(0);	
}


//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof1(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = Y[0]; //CR
	Y1alfa[0]  = 1.-Y[0]; //FE
	Y2alfa[0]  = y_fraction_c_bcc; //C
	Y2alfa[1]  = 1.-y_fraction_c_bcc; //VA
	Y1betta[1] = Y[2]; //CR
	Y1betta[0] = 1.-Y[2]; //FE
	Y2betta[0] = Y[1]; //C
	Y2betta[1] = 1.-Y[1]; //VA


	
	
	std::vector <double> pot_chemical1alfa(3);
	std::vector <double> pot_chemical1betta(3);
	
	std::vector <double> Diffbetta(2);	
	Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	
	//#####Muc
	vector <double> L(3);
	L=deriv_mu_FeCCr_KTT_prof5(Uc_gamma,Ucr_gamma,1.,concentration);
	double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
	double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
	muc = (Y[3]*taille_gradient/Vm)*(num_c/dem_c);
	
	//####MuCr
	double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
	double dem_cr=L[1]*L[0]-pow(L[2],2);
	mucr= (Y[3]*taille_gradient/Vm)*(num_cr/dem_cr);
	
	//#####Integrale Gtr
	int i,n1;

	//gtr= ((Vm/Y[3])*(1./taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
	
	gtr= ((Vm/Y[3])*(taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)));
	
	gtr1=Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2));
	
	F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc/(R*TT2); // carbone
	
	F[1] = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2))-gtr/(R*TT2)/*10.-Y[3]*(Vm/M)/(R*TT2)*/;// iron
	
	F[2] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr/(R*TT2); //cr-fe
	
	
	std::vector <double> delta(2),omega(2);
	
	
	omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
	omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);
	
	
	
	for (i =0; i < 2; i++)
	{
		delta[i] = 2.*(1. - omega[i])/max(omega[i],1e-10);//omega[i] ;
	}	 
	
	//for (i = 0; i < 2; i++)
	//	if (delta[i]<1.e-16) delta[i]=1.e-15;
		//if (delta[i]<1.e-16) delta[i]=abs(delta[i]);
		//if(omega[i]<0.0){delta[i]=abs(delta[i]);printf("sursaturation <=0.0 FK1\n");}
	//delta[i] = 1.e-2+1./pow(cosh(omega[i]/0.3),2.);

	
	F[3] = (Uc_gamma-Uc_alpha)*Y[3]-(Diffbetta[0]/( Rzin*delta[0]))* \
	(Uc_gamma-U0c); //C (gamma->alpha)	
	
}

//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof2(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = y_fraction_x_bcc; //CR
	Y1alfa[0]  = 1.-y_fraction_x_bcc; //FE
	Y2alfa[0]  = y_fraction_c_bcc; //C
	Y2alfa[1]  = 1.-y_fraction_c_bcc; //VA
	Y1betta[1] = y_fraction_x_fcc; //CR
	Y1betta[0] = 1.-y_fraction_x_fcc; //FE
	Y2betta[0] = y_fraction_c_fcc; //C
	Y2betta[1] = 1.-y_fraction_c_fcc; //VA


	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	std::vector <double> delta(1),omega(1);
	if(Ucr_gamma!=Ucr_alpha){
	
		omega[0] = (Ucr_gamma-Y[0])/(Ucr_gamma-Ucr_alpha);
		delta[0] = 2.*(1. - omega[0])/max(omega[0],1e-10) ;
	}
	//if(delta[0]<=0.0) delta[0]=abs(delta[0]);
	//if(delta[0]<1.e-16) delta[0]=1.e-15;
	//if(omega[0]<=0.0) {delta[0]=abs(delta[0]);printf("sursaturation <=0.0 FK2\n");}
	//delta[0] = 1.e-2+1./pow(cosh(omega[0]/0.3),2.);

	
	//printf("delta=%e omega_x=%e \n",delta[0],omega[0]);
	//printf("y=%e \n",Y[0]);
	F[0] = (Ucr_gamma-Ucr_alpha)*VIT_TEMPORAIRE*Rzin*delta[0]-dcr*(Ucr_gamma-Y[0]); //CR (gamma ->alpha)
	//printf("F[0]=%e Ucr_fcc=%e Y[0]=%e \n",F[0],Ucr_gamma,Y[0]);
	//printf("term=%e \n",(Ucr_gamma-Ucr_alpha)*VIT_TEMPORAIRE*Rzin*delta[0]);
	
}

//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof3(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = y_fraction_x_bcc; //CR
	Y1alfa[0]  = 1.-y_fraction_x_bcc; //FE
	Y2alfa[0]  = Y[0]; //C
	Y2alfa[1]  = 1.-Y[0]; //VA
	Y1betta[1] = Y[2]; //CR
	Y1betta[0] = 1.-Y[2]; //FE
	Y2betta[0] = Y[1]; //C
	Y2betta[1] = 1.-Y[1]; //VA
	
	
	
	
	std::vector <double> pot_chemical1alfa(3);
	std::vector <double> pot_chemical1betta(3);
	
	std::vector <double> Diffbetta(2);	
	Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	
	//#####Muc
	vector <double> L(3);
	L=deriv_mu_FeCCr(Uc_gamma,Ucr_gamma,1.,concentration);
	double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
	double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
	muc = (Y[3]*taille_gradient/Vm)*(num_c/dem_c);
	
	//####MuCr
	double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
	double dem_cr=L[1]*L[0]-pow(L[2],2);
	mucr= (Y[3]*taille_gradient/Vm)*(num_cr/dem_cr);
	
	//#####Integrale Gtr
	int i,n1;
	
	/*gtr= ((Vm/Y[3])*(1./taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);*/
	gtr= ((Vm/Y[3])*(taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)));
	
	
	F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc/(R*TT2); // carbone
	
	F[1] = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2))-gtr/(R*TT2);// iron
	
	F[2] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr/(R*TT2); //cr-fe
	
	
	std::vector <double> delta(2),omega(2);
	
	
	omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
	omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);
	
	
	
	for (i =0; i < 2; i++)
	{
		delta[i] = 2.*(1. - omega[i])/max(omega[i],1e-10);//omega[i] ;
	}	 
	
	//for (i = 0; i < 2; i++)
	//if (delta[i]<1.e-16) delta[i]=1.e-15;
	//if (delta[i]<1.e-16) delta[i]=abs(delta[i]);
	//if(omega[i]<0.0){delta[i]=abs(delta[i]);printf("sursaturation <=0.0 FK1\n");}
	//delta[i] = 1.e-2+1./pow(cosh(omega[i]/0.3),2.);
	
	
	F[3] = (Ucr_gamma-Ucr_alpha)*Y[3]-(Diffbetta[1]/( Rzin*delta[1]))* \
	(Ucr_gamma-U0cr); //C (gamma->alpha)	
	
}


//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof4(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = y_fraction_x_bcc; //CR
	Y1alfa[0]  = 1.-y_fraction_x_bcc; //FE
	Y2alfa[0]  = y_fraction_c_bcc; //C
	Y2alfa[1]  = 1.-y_fraction_c_bcc; //VA
	Y1betta[1] = y_fraction_x_fcc; //CR
	Y1betta[0] = 1.-y_fraction_x_fcc; //FE
	Y2betta[0] = y_fraction_c_fcc; //C
	Y2betta[1] = 1.-y_fraction_c_fcc; //VA
	
	
	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	std::vector <double> delta(1),omega(1);
	if(Uc_gamma!=Uc_alpha){
		
		omega[0] = (Uc_gamma-Y[0])/(Uc_gamma-Uc_alpha);
		delta[0] = 2.*(1. - omega[0])/max(omega[0],1e-10) ;
	}
	//if(delta[0]<=0.0) delta[0]=abs(delta[0]);
	//if(delta[0]<1.e-16) delta[0]=1.e-15;
	//if(omega[0]<=0.0) {delta[0]=abs(delta[0]);printf("sursaturation <=0.0 FK2\n");}
	//delta[0] = 1.e-2+1./pow(cosh(omega[0]/0.3),2.);
	
	
	printf("delta=%e omega_x=%e \n",delta[0],omega[0]);
	printf("y=%e \n",Y[0]);
	F[0] = (Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE*Rzin*delta[0]-dc*(Uc_gamma-Y[0]); //CR (gamma ->alpha)
	printf("F[0]=%e Ucr_fcc=%e Y[0]=%e \n",F[0],Uc_gamma,Y[0]);
	printf("term=%e \n",(Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE*Rzin*delta[0]);
	
}


//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof5(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
double TT2 = 0.;
TT2 = ph1->Thermo->T1;
std::vector <double> Y1alfa(2);
std::vector <double> Y2alfa(2);
std::vector <double> Y1betta(2);
std::vector <double> Y2betta(2);
	
Y1alfa[1]  = Y[1]; //CR
Y1alfa[0]  = 1.-Y[1]; //FE
Y2alfa[0]  = Y[0]; //C
Y2alfa[1]  = 1.-Y[0]; //VA
Y1betta[1] = Y[3]; //CR
Y1betta[0] = 1.-Y[3]; //FE
Y2betta[0] = Y[2]; //C
Y2betta[1] = 1.-Y[2]; //Va
std::vector <double> pot_chemical1alfa(3);
std::vector <double> pot_chemical1betta(3);
	
//int dim = ph2->Thermo->Dex.size();
std::vector <double> Diffbetta(2);	
Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
//double compteur_fe;
double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
double U0c			=	concentration[0]/(1.-concentration[0]);
double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	
//#####Muc
vector <double> L(3);
	
L=deriv_mu_FeCCr_KTT_prof5(Uc_gamma,Ucr_gamma,1.,concentration);

double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
muc = (VIT_TEMPORAIRE*taille_gradient/Vm)*(num_c/dem_c);
//printf("muc=%e \n",muc);
//####MuCr
double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
double dem_cr=L[1]*L[0]-pow(L[2],2);
mucr= (VIT_TEMPORAIRE*taille_gradient/Vm)*(num_cr/dem_cr);
	
	//#####Integrale Gtr
int i,n1;

gtr= ((Vm/VIT_TEMPORAIRE)*(taille_gradient)* \
	( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
	((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
	
//gtr= ((Vm/VIT_TEMPORAIRE)*(taille_gradient)* \
	( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
	((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
	
gtr1 = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2));
		
F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc/(R*TT2); // carbone
F[1] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr/(R*TT2); //cr-fe
	
	
std::vector <double> delta(2),omega(2);
	
	
omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);

	
for (i =0; i < 2; i++)
{
	delta[i] = 2.*(1. - omega[i])/omega[i] ;
}	 

for (i = 0; i < 2; i++)
	if (delta[i]<1.e-16) delta[i]=1.e-15;

F[2] = (Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE-(Diffbetta[0]/( Rzin*delta[0]))* \
	(Uc_gamma-U0c); //C (gamma->alpha)
	
F[3] = (Ucr_gamma-Ucr_alpha)*VIT_TEMPORAIRE-(Diffbetta[1]/(Rzin*delta[1]))* \
	(Ucr_gamma-U0cr); //CR (gamma ->alpha)

	
/*	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = Y[1]; //CR
	Y1alfa[0]  = 1.-Y[1]; //FE
	Y2alfa[0]  = Y[0]; //C
	Y2alfa[1]  = 1.-Y[0]; //VA
	Y1betta[1] = Y[3]; //CR
	Y1betta[0] = 1.-Y[3]; //FE
	Y2betta[0] = Y[2]; //C
	Y2betta[1] = 1.-Y[2]; //VA
	
	//########
	
	
	std::vector <double> pot_chemical1alfa(3);
	std::vector <double> pot_chemical1betta(3);
	
	//int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(2);	
	//Diffbetta = calc_diffusion_betta(concentration);
	//Diffbetta[0] = Diffbetta[0]*pow(1.e6,2.);//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	//Diffbetta[1] = Diffbetta[1]*pow(1.e6,2.);//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	
	//#####Muc
	vector <double> L(3);
	L=deriv_mu_FeCCr(Uc_gamma,Ucr_gamma,1.,concentration);
	double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
	double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
	muc = (VIT_TEMPORAIRE*taille_gradient/Vm)*(num_c/dem_c) ; // /(R*TT2); VRAI
	
	//####MuCr
	double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
	double dem_cr=L[1]*L[0]-pow(L[2],2);
	mucr= (VIT_TEMPORAIRE*taille_gradient/Vm)*(num_cr/dem_cr) ; // /(R*TT2); VRAI
	
	//#####Integrale Gtr
	int i,n1;

	gtr= ((Vm/VIT_TEMPORAIRE)*(taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
	
	
	gtr1 = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2)) ;
	
	
	//delta_mu_fe = (pot_chemical1betta[1] - pot_chemical1alfa[1]) / (R * TT2) ;
	//delta_mu_ni = (pot_chemical1betta[2] - pot_chemical1alfa[2]) / (R * TT2) ;
	
	F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc/(R*TT2); // carbone VRAI
	

	
	
	F[1] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr/(R*TT2); //cr-fe VRAI
	
	
	std::vector <double> delta(2),omega(2);
	//delta = Zenerplane_cin(Y,concentration);
	
	
	omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
	omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);
	
	for (i =0; i < 2; i++)
	{
		delta[i] = 2*(1. - omega[i])/omega[i] ;
	}	 

	for (i = 0; i < 2; i++)		
		if (delta[i]<1.e-18) delta[i]=1.e-17;
	
	
	F[2] = (
			(Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE-(Diffbetta[0]/( Rzin*delta[0]))* \
			(Uc_gamma-U0c)
			) ;/// Diffbetta[1] ;
	
	//C (gamma->alpha)
	F[3] = (
			(Ucr_gamma-Ucr_alpha)*VIT_TEMPORAIRE-(Diffbetta[1]/(Rzin*delta[1]))* \
			(Ucr_gamma-U0cr)
			) * 1.e0 ;/// Diffbetta[0] ;
	
*/
}


//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof6(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = Y[0]; //CR
	Y1alfa[0]  = 1.-Y[0]; //FE
	Y2alfa[0]  = y_fraction_c_bcc; //C
	Y2alfa[1]  = 1.-y_fraction_c_bcc; //VA
	Y1betta[1] = Y[2]; //CR
	Y1betta[0] = 1.-Y[2]; //FE
	Y2betta[0] = Y[1]; //C
	Y2betta[1] = 1.-Y[1]; //VA
	
	
	
	
	std::vector <double> pot_chemical1alfa(3);
	std::vector <double> pot_chemical1betta(3);
	
	std::vector <double> Diffbetta(2);	
	Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	
	//#####Muc
	vector <double> L(3);
	L=deriv_mu_FeCCr(Uc_gamma,Ucr_gamma,1.,concentration);
	double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
	double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
	muc = (Y[3]*taille_gradient/Vm)*(num_c/dem_c);
	
	//####MuCr
	double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
	double dem_cr=L[1]*L[0]-pow(L[2],2);
	mucr= (Y[3]*taille_gradient/Vm)*(num_cr/dem_cr);
	
	//#####Integrale Gtr
	int i,n1;
	
	//gtr= ((Vm/Y[3])*(1./taille_gradient)* \
	( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
	((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
	
	gtr= ((Vm/Y[3])*(taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)));
	
	gtr1=Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2));
	
	F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc/(R*TT2); // carbone
	
	F[1] = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2))-gtr/(R*TT2)/*10.*//*-Y[3]*(Vm/M)/(R*TT2)*/;// iron
	
	F[2] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr/(R*TT2); //cr-fe
	
	
	std::vector <double> delta(2),omega(2);
	
	
	omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
	omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);
	
	
	
	for (i =0; i < 2; i++)
	{
		delta[i] = 2.*(1. - omega[i])/max(omega[i],1e-10);//omega[i] ;
	}	 
	
	//for (i = 0; i < 2; i++)
	//if (delta[i]<1.e-16) delta[i]=1.e-15;
	//if (delta[i]<1.e-16) delta[i]=abs(delta[i]);
	//if(omega[i]<0.0){delta[i]=abs(delta[i]);printf("sursaturation <=0.0 FK1\n");}
	//delta[i] = 1.e-2+1./pow(cosh(omega[i]/0.3),2.);
	
	
	F[3] = (Ucr_gamma-Ucr_alpha)*Y[3]-(Diffbetta[1]/( Rzin*delta[1]))* \
	(Ucr_gamma-U0cr); //Ni (gamma->alpha)
	
}


//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT_prof7(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = y_fraction_x_bcc; //CR
	Y1alfa[0]  = 1.-y_fraction_x_bcc; //FE
	Y2alfa[0]  = Y[0]; //C
	Y2alfa[1]  = 1.-Y[0]; //VA
	Y1betta[1] = Y[2]; //CR
	Y1betta[0] = 1.-Y[2]; //FE
	Y2betta[0] = Y[1]; //C
	Y2betta[1] = 1.-Y[1]; //VA
	
	
	
	
	std::vector <double> pot_chemical1alfa(3);
	std::vector <double> pot_chemical1betta(3);
	
	std::vector <double> Diffbetta(2);	
	Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	
	//#####Muc
	vector <double> L(3);
	L=deriv_mu_FeCCr(Uc_gamma,Ucr_gamma,1.,concentration);
	double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
	double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
	muc = (VIT_TEMPORAIRE*taille_gradient/Vm)*(num_c/dem_c);
	
	//####MuCr
	double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
	double dem_cr=L[1]*L[0]-pow(L[2],2);
	mucr= (VIT_TEMPORAIRE*taille_gradient/Vm)*(num_cr/dem_cr);
	
	//#####Integrale Gtr
	int i,n1;
	
	//gtr= ((Vm/Y[3])*(1./taille_gradient)* \
	( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
	((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
	
	gtr= ((Vm/VIT_TEMPORAIRE)*(taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)));
	
	gtr1=Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2));
	
	F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc/(R*TT2); // carbone
	
	//F[1] = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2))-gtr/(R*TT2)/*10.*//*-Y[3]*(Vm/M)/(R*TT2)*/;// iron
	
	F[1] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr/(R*TT2); //cr-fe
	
	
	std::vector <double> delta(2),omega(2);
	
	
	omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
	omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);
	
	
	//
	for (i =0; i < 2; i++)
	{
		delta[i] = 2.*(1. - omega[i])/max(omega[i],1e-10);//omega[i] ;
	}	 
	
	//for (i = 0; i < 2; i++)
	//if (delta[i]<1.e-16) delta[i]=1.e-15;
	//if (delta[i]<1.e-16) delta[i]=abs(delta[i]);
	//if(omega[i]<0.0){delta[i]=abs(delta[i]);printf("sursaturation <=0.0 FK1\n");}
	//delta[i] = 1.e-2+1./pow(cosh(omega[i]/0.3),2.);
	
	
	F[2] = (Uc_gamma-Uc_alpha)*VIT_TEMPORAIRE-(Diffbetta[0]/( Rzin*delta[0]))* \
	(Uc_gamma-U0c); //C (gamma->alpha)	
	
	//F[2] = (Ucr_gamma-Ucr_alpha)*VIT_TEMPORAIRE-(Diffbetta[1]/( Rzin*delta[1]))* \
	(Ucr_gamma-U0cr); //C (gamma->alpha)	
}


//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT_prof2(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof2(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof2(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	/*
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = y_fraction_x_bcc; //CR
	Y1alfa[0]  = 1.-y_fraction_x_bcc; //FE
	Y2alfa[0]  = y_fraction_c_bcc; //C
	Y2alfa[1]  = 1.-y_fraction_c_bcc; //VA
	Y1betta[1] = y_fraction_x_fcc; //CR
	Y1betta[0] = 1.-y_fraction_x_fcc; //FE
	Y2betta[0] = y_fraction_c_fcc; //C
	Y2betta[1] = 1.-y_fraction_c_fcc; //VA
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	Y[0]/(1.-concentration[0]);
	
	std::vector <double> delta(1),omega(1);
	
	if(Ucr_gamma!=Ucr_alpha){
		
		omega[0] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);
	}
	//if(omega[0]<0.0)omega[0]=abs(omega[0]);
	//Jacobian3D[0][0]=dcr-2.*Rzin/pow(omega[0],2.)/(Ucr_alpha-Ucr_gamma);
	printf("omega=%e \n",omega[0]);
	Jacobian3D[0][0]=dcr+2.*Rzin*VIT_TEMPORAIRE/pow(omega[0],2.);
*/
	
}

//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT_prof4(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof4(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof4(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
}

void interface::JKT_prof5(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
int i,j;
int n;
int i1 = Y.size();
n = i1;
std::vector <double> F(i1);
std::vector <double> Ynew(i1);
std::vector < std::vector<double> > Fnew;

Fnew.resize(i1);
Jacobian3D.resize(i1);

for (i = 0; i < Jacobian3D.size(); i++)
{
Fnew[i].resize(i1);
Jacobian3D[i].resize(i1);
}
for( i=0 ; i < i1 ; i++)
{
Ynew[i] = Y[i];
}  	
for( i=0 ; i < Jacobian3D.size() ; i++)
 for(j=0 ; j < Jacobian3D[i].size() ; j++) 
 {
Ynew[j] = Y[j] + kdelta;
KTT_prof5(n,Ynew,F,concentration);
Fnew[i][j] = F[i];
KTT_prof5(n,Y,F,concentration); 
Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
Ynew[j] = Y[j]; 
}
}


void interface::JKT_prof6(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof6(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof6(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
}

void interface::JKT_prof7(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof7(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof7(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
}


void interface::KTT_prof_sphere(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	
	Y1alfa[1]  = Y[1]; //CR
	Y1alfa[0]  = 1.-Y[1]; //FE
	Y2alfa[0]  = Y[0]; //C
	Y2alfa[1]  = 1.-Y[0]; //VA
	Y1betta[1] = Y[3]; //CR
	Y1betta[0] = 1.-Y[3]; //FE
	Y2betta[0] = Y[2]; //C
	Y2betta[1] = 1.-Y[2]; //VA

	
	std::vector <double> pot_chemical1alfa(3);
	std::vector <double> pot_chemical1betta(3);
	
	std::vector <double> Diffbetta(2);	
	Diffbetta[0] = dc;//cout<<"Diffbetta[0]="<<Diffbetta[0]<<endl;
	Diffbetta[1] = dcr;//cout<<"Diffbetta[1]="<<Diffbetta[1]<<endl;
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	//double compteur_fe;
	double x_c_alfa		=	ph1->Thermo->m[1]*Y2alfa[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_cr_alfa	=	ph1->Thermo->m[0]*Y1alfa[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y2alfa[0]);
	double x_fe_alfa	=	1.-x_c_alfa-x_cr_alfa;
	double x_c_gamma	=	ph2->Thermo->m[1]*Y2betta[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_cr_gamma	=	ph2->Thermo->m[0]*Y1betta[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y2betta[0]);
	double x_fe_gamma	=	1.-x_c_gamma-x_cr_gamma;
	
	double Uc_alpha		=	x_c_alfa/(1.-x_c_alfa);
	double Ucr_alpha	=	x_cr_alfa/(1.-x_c_alfa);
	double Ufe_alpha	=	x_fe_alfa/(1.-x_c_alfa);
	double Uc_gamma		=	x_c_gamma/(1.-x_c_gamma);
	double Ucr_gamma	=	x_cr_gamma/(1.-x_c_gamma);
	double Ufe_gamma	=	x_fe_gamma/(1.-x_c_gamma);
	double U0c			=	concentration[0]/(1.-concentration[0]);
	double U0cr			=	concentration[2]/(1.-concentration[0]);
	
	
	//#####Muc
	vector <double> L(3);
	L=deriv_mu_FeCCr(Uc_gamma,Ucr_gamma,1.,concentration);
	double num_c=L[2]*L[0]*(Ucr_gamma-Ucr_alpha)-L[1]*L[0]*(Uc_gamma-Uc_alpha);
	double dem_c=L[1]*pow(L[0],2)-L[0]*pow(L[2],2);
	muc = (Y[4]*taille_gradient/Vm)*(num_c/dem_c)/(R*TT2);
	
	//####MuCr
	double num_cr=-L[0]*(Ucr_gamma-Ucr_alpha)+L[2]*(Uc_gamma-Uc_alpha);
	double dem_cr=L[1]*L[0]-pow(L[2],2);
	mucr= (Y[4]*taille_gradient/Vm)*(num_cr/dem_cr)/(R*TT2);
	
	//#####Integrale Gtr
	int i,n1;

	gtr= ((Vm/Y[4])*(1./taille_gradient)* \
		  ( ( (L[2]/taille_gradient)*muc+(L[1]/taille_gradient)*mucr)*(mucr/taille_gradient)+ \
		   ((L[0]/taille_gradient)*muc+(L[2]/taille_gradient)*mucr)*(muc/taille_gradient)))/(R*TT2);
 
	
	F[0] = (pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2)-muc; // carbone
	
	F[1] = Ufe_alpha*((pot_chemical1betta[1] - pot_chemical1alfa[1])/(R*TT2))+ \
	Uc_alpha*((pot_chemical1betta[0] - pot_chemical1alfa[0])/(R*TT2))+ \
	Ucr_alpha*((pot_chemical1betta[2]-pot_chemical1alfa[2])/(R*TT2))-gtr/*-Y[4]*(Vm/M)/(R*TT2)*/;// iron

	
	
	F[2] = (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2) \
	- (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2)-mucr; //cr-fe
	
	
	std::vector <double> delta(2),omega(2),F_omega(2);
	
	
	omega[0] = (Uc_gamma-U0c)/(Uc_gamma-Uc_alpha);
	omega[1] = (Ucr_gamma-U0cr)/(Ucr_gamma-Ucr_alpha);

	
	//for (i =0; i < 2; i++)
	//{
	//	F_omega[i]=omega[i]/pow((1.-omega[i]),(1.+1./50.));
		//F_omega[i]=0.95*pow((omega[i]/(1.-omega[i])),(1.+1./50.));
		//F_omega[i]=-2.*omega[i]/log(omega[i]);
	//}
	
	for (i=0;i<2;i++)
	 {	
		 delta[i]=(1.-omega[i]);
	 }	  
	
	for (i = 0; i < 2; i++)
		if (delta[i]<1.e-16) delta[i]=1.e-15;
	

	
	F[3] = (Uc_gamma-Uc_alpha)*Y[4]-(Diffbetta[1]*correc_D/( Rzin*delta[0]))* \
	(Uc_gamma-U0c); //C (gamma->alpha)
	F[4] = (Ucr_gamma-Ucr_alpha)*Y[4]-(Diffbetta[1]/(Rzin*delta[1]))* \
	(Ucr_gamma-U0cr); //CR (gamma ->alpha)
	
}



//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
int i,j;
int n;
int i1 = Y.size();
n = i1;
std::vector <double> F(i1);
std::vector <double> Ynew(i1);
std::vector < std::vector<double> > Fnew;

Fnew.resize(i1);
Jacobian3D.resize(i1);

for (i = 0; i < Jacobian3D.size(); i++)
{
Fnew[i].resize(i1);
Jacobian3D[i].resize(i1);
}	
for( i=0 ; i < i1 ; i++)
{
Ynew[i] = Y[i];
}  	
for( i=0 ; i < Jacobian3D.size() ; i++)
	for(j=0 ; j < Jacobian3D[i].size() ; j++) 
{
Ynew[j] = Y[j] + kdelta;
KTT(n,Ynew,F,concentration);
Fnew[i][j] = F[i];
KTT(n,Y,F,concentration); 
Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
Ynew[j] = Y[j]; 
}
}



//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT_sphere(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_sphere(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_sphere(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}

//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT_prof(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}


//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT_prof1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof1(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof1(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}

//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT_prof3(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof3(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof3(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}


void interface::JKT_prof_sphere(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian3D.resize(i1);
	
	for (i = 0; i < Jacobian3D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian3D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KTT_prof_sphere(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT_prof_sphere(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}

std::vector <double> interface::Zenerplane(std::vector <double> &Y,std::vector <double> &concentration)
{
	std::vector<double> delta(2);
	std::vector<double> omega(2);

	
	omega[0] = ( ( (ph2->Thermo->m[1]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7]) ) - concentration[0] ) \
				/( ( (ph2->Thermo->m[1]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7])  ) \
				  -( (ph1->Thermo->m[1]*Y[2])/( ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[3]) ) )   ;//c
	
	
	omega[1] = ( ( (ph2->Thermo->m[0]*Y[5])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7]) ) - concentration[2] ) \
				/( ( (ph2->Thermo->m[0]*Y[5])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7]) ) \
				-  ( (ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[3]) ) );//cr
	
	/*for (int i = 0; i < 2; i++){
		if (omega[i]==1.)omega[i]=0.999999;
		if (omega[i]==0.)omega[i]=0.000001;}*/
	
	
	for (int i =0; i < 2; i++)
	 {
		 delta[i] = 2*(1. - omega[i])/omega[i] ;
	 }	 
		
	for (int i = 0; i < 2; i++)
		if (delta[i]<1.e-16) delta[i]=1.e-15;


	return delta;			
				
}	


std::vector <double> interface::Zenerplane_cin(std::vector <double> &Y,std::vector <double> &concentration)
{
	std::vector<double> delta(2);
	std::vector<double> omega(2);
	
	
	
	omega[0] = ( ( (ph2->Thermo->m[1]*Y[2])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*(1.-Y[2])) ) - concentration[0] ) \
	/( ( (ph2->Thermo->m[1]*Y[2])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*(1.-Y[2]))  ) \
	  -( (ph1->Thermo->m[1]*Y[0])/( ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*(1.-Y[0])) ) )   ;//c
	//omega[0] = (x_cFeCrC[11]-concentration[0])/(x_cFeCrC[11]-x_cFeCrC[1]);
	
	
	omega[1] = ( ( (ph2->Thermo->m[0]*Y[3])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*(1.-Y[2])) ) - concentration[2] ) \
	/( ( (ph2->Thermo->m[0]*Y[3])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*(1.-Y[2])) ) \
	  -  ( (ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*(1.-Y[0])) ) );//cr
	//omega[1] = (x_crFeCrC[11]-concentration[2])/(x_crFeCrC[11]-x_crFeCrC[1]); 
	
	/*for (int i = 0; i < 2; i++){
	 if (omega[i]==1.)omega[i]=0.999999;
	 if (omega[i]==0.)omega[i]=0.000001;}*/
	
	
	for (int i =0; i < 2; i++)
	{
		delta[i] = 2*(1. - omega[i])/omega[i] ;
	}	 
	
	for (int i = 0; i < 2; i++)
		if (delta[i]<1.e-16) delta[i]=1.e-15;
	
	
	return delta;			
	
}	


std::vector <double> interface::Zenerplane2d(std::vector <double> &Y,std::vector <double> &concentration)
{
	std::vector<double> delta(1);
	std::vector<double> omega(1);

	
	omega[0] = ( ( (ph2->Thermo->m[1]*Y[4])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y[4]) ) - concentration[0] ) \
	/( ( (ph2->Thermo->m[1]*Y[4])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y[4])  ) \
	  -( (ph1->Thermo->m[1]*Y[1])/( ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y[1]) ) )   ;//c

	
	/*for (int i = 0; i < 2; i++){
	 if (omega[i]==1.)omega[i]=0.999999;
	 if (omega[i]==0.)omega[i]=0.000001;}*/
	
	
	for (int i =0; i < 1; i++)
	{
		delta[i] = 2*(1. - omega[i])/omega[i] ;
	}	 
	
	for (int i = 0; i < 1; i++)
		if (delta[i]<1.e-16) delta[i]=1.e-15;
	
	
	return delta;


}

std::vector <double> interface::Zenerplane2d_prof(std::vector <double> &Y,std::vector <double> &concentration)
{
	std::vector<double> delta(1);
	std::vector<double> omega(1);
	
	
	omega[0] =( ( (ph2->Thermo->m[1]*Y[1])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y[1]) ) - concentration[0] ) \
				/(( (ph2->Thermo->m[1]*Y[1])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]*Y[1]) ) \
				  -( (ph1->Thermo->m[1]*Y[0])/( ph1->Thermo->m[0]+ph1->Thermo->m[1]*Y[0]) ) );

	
	//omega[0] = (x_c[10]-concentration[0])/(x_c[10]-x_c[0]);//c
	
	
	/*for (int i = 0; i < 2; i++){
	 if (omega[i]==1.)omega[i]=0.999999;
	 if (omega[i]==0.)omega[i]=0.000001;}*/
	
	
	for (int i =0; i < 1; i++)
	{
		delta[i] = 2*(1. - omega[i])/omega[i] ;
	}	 
	
	for (int i = 0; i < 1; i++)
		if (delta[i]<1.e-16) delta[i]=1.e-15;
	
	
	return delta;
	
	
}


//_________________________kinetic transition Fe-Cr-Mo-C____________________
void interface::KT4(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(3);//bcc_A2: FE,CR,MO
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(3);//fcc_A1: FE,CR,MO
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	Y1alfa[0] = Y[0];//fe
	Y1alfa[1] = Y[1];//cr
	Y1alfa[2] = Y[2];//mo
	Y2alfa[0] = Y[3];//c
	Y2alfa[1] = Y[4];//va
	
	Y1betta[0]= Y[5];//fe
	Y1betta[1]= Y[6];//cr
	Y1betta[2]= Y[7];//mo
	Y2betta[0]= Y[8];//c 
	Y2betta[1]= Y[9];//va
	
	std::vector <double> pot_chemical1alfa(4);//bcc_A2: 0 - C; 1 - FE; 2 - CR; 3 - MO
	std::vector <double> pot_chemical1betta(4);//fcc_A1: 0 - C; 1 - FE; 2 - CR; 3 - MO
	
	int dim = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(dim);	
	Diffbetta = calc_diffusion_betta(concentration);
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	
	F[0] = (pot_chemical1alfa[0] - pot_chemical1betta[0])/(R*TT2); //C_bcc - C_fcc
	F[1] = (pot_chemical1alfa[1] - pot_chemical1betta[1])/(R*TT2); //Fe_bcc - Fe_fcc
	F[2] = (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2) - (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2); //(Cr_bcc - Fe_bcc) - (Cr_fcc - Fe_fcc)
	F[3] = (pot_chemical1alfa[3] - pot_chemical1alfa[1])/(R*TT2) - (pot_chemical1betta[3] - pot_chemical1betta[1])/(R*TT2); //(Mo_bcc - Fe_bcc) - (Mo_fcc - Fe_fcc)
	F[4] = Y[0] + Y[1] + Y[2] - 1.; //Fe+Cr+Mo=1
	F[5] = Y[3] + Y[4] - 1.; //C+Va=1
	F[6] = Y[5] + Y[6] + Y[7] - 1.; //Fe+Cr+Mo=1
	F[7] = Y[8] + Y[9] - 1.; //C+Va=1
	
	F[8] = (-((ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]))*Y[10] \
	+((ph2->Thermo->m[0]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]))*(Y[10]-(Diffbetta[1]/taille_gradient)) \
	+(Diffbetta[1]/taille_gradient)*concentration[2]); //CR (gamma->alpha)

	
	F[9] = (-((ph1->Thermo->m[0]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]))*Y[10] \
	+((ph2->Thermo->m[0]*Y[7])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]))*(Y[10]-(Diffbetta[2]/taille_gradient)) \
	+(Diffbetta[2]/taille_gradient)*concentration[3]); //Mo (gamma->alpha)
	


	F[10] = (-((ph1->Thermo->m[1]*Y[3])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]))*Y[10] \
	+((ph2->Thermo->m[1]*Y[8])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]))*(Y[10]-(Diffbetta[0]/taille_gradient)) \
	+(Diffbetta[0]/taille_gradient)*concentration[0]); //C (gamma->alpha)

}

//_______________Jacobien cinetique Fe-Cr-Mo-C____________________
void interface::JK4(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian4D,std::vector <double> &concentration)
{
	int i,j;
	int n;
	int i1 = Y.size();
	n = i1;
	std::vector <double> F(i1);
	std::vector <double> Ynew(i1);
	std::vector < std::vector<double> > Fnew;
	
	Fnew.resize(i1);
	Jacobian4D.resize(i1);
	
	for (i = 0; i < Jacobian4D.size(); i++)
	{
		Fnew[i].resize(i1);
		Jacobian4D[i].resize(i1);
	}
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian4D.size() ; i++)
		for(j=0 ; j < Jacobian4D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			KT4(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KT4(n,Y,F,concentration); 
			Jacobian4D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
	
}


