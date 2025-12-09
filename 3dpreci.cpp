/*
 *  3d.cpp
 *  
 *
 *  Created by nat on 22/03/10.
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
#include <map>
#include <time.h>


extern double temps,dt;
extern double T,dT;
extern float taille_gradient;
extern double dtemps;
extern float taille_gradient_alpha;
extern float taille_gradient_betta;
extern segment seg[20];
extern int nbseg;
double correc_D;
float CIN;
extern int def;
extern double M;
extern double Vm;
extern double x_c0;
double x_cr0;
double x_i;
double x_j;

double compteur_c;
double compteur_cr;
//double compteur_fe;
extern double x_cFeCrC[11];
extern double x_crFeCrC[11];
double diffFeCrC;
double fricFeCrC;
extern double VIT_TEMPORAIRE;
extern double Rzin;
double coeff_n;
double coeff_k;
double dc , dcr;
double Lcoef,Bcoef;

extern double gtr;
extern double muc;
extern double mucr;

extern sortie so;
extern sortie so1;
extern sortie so2;
extern sortie so3;
extern sortie so4;
extern sortie so5;
extern sortie so6;
extern sortie so7;
extern sortie so8;

extern matrice mat;

extern distribution_moyenne<sphere> cem;
extern distribution_moyenne<sphere> eps;

void calc3dequilibre();
bool calc3cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc3cinetique2(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc3cinetique_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage);

extern void ecriture_growth_dissolution_FeCCr(std::vector <double> &t, \
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

extern void ecriture3D(double temps,std::vector <double> &Gtr, \
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
extern void ecriture_profileC_FeCCr(double temps,double Rzin,double SURSAT_c,double SURSAT_cr,int p, double *tx, \
							 double *ty, double xi, double xf);
extern void ecriture_profilePOT_FeCCr(double temps,double Rzin,int m,std::vector < std::vector<double> > &P, \
							   double xi, double xf);
extern void cryptograph(std::vector <double> &L,std::vector <double> &B, \
				 std::vector <double> &Fx_c_bcc,std::vector <double> &Fx_cr_bcc, \
				 std::vector <double> &Fx_c_fcc,std::vector <double> &Fx_cr_fcc);

clock_t t_debut_main , t_debut_qhull , t_fin_qhull , t_debut_NR , t_fin_NR , t_fin_main ;

void calc3dequilibre()
{	
	t_debut_main = clock();
	int dim=3;
	def=3;
	
	newt nt;
    qhull qh(dim);
	
	double prevT=0;
	bool conv = false;
	bool flage = false;

	std::vector<double> prevVector(9);		
	std::vector<double> prevVectorcin(9);
	for(int i=0;i<9;i++){prevVectorcin[i]=0;}		
	for(int i=0;i<9;i++){prevVector[i]=0;}	
	
	std::vector <double> CONCENTRATION_NOMINAL ( 3 );

	//CONCENTRATION_NOMINAL[0] = 9.2149836e-3;//9.2149836e-3;      2e-2        		            // carbon -C
	//CONCENTRATION_NOMINAL[2] = 2.6608097e-2;//2.6608097e-2;  	5.e-1                            // chrome -CR
	//CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0]-CONCENTRATION_NOMINAL[2];// iron -FE
	CONCENTRATION_NOMINAL[0] = 0.4e-2;//9.2149836e-3;      2e-2        		            // carbon -C
	CONCENTRATION_NOMINAL[2] = 2.e-2;//2.6608097e-2;  	5.e-1                            // chrome -CR
	CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0]-CONCENTRATION_NOMINAL[2];// iron -FE
	x_i=CONCENTRATION_NOMINAL[0];
	x_j=CONCENTRATION_NOMINAL[2];
	
	
	//double iii,jjj;
	//for(iii=1;iii<=50;iii+=10)
		//for(jjj=1;jjj<=1000;jjj+=10)
		{		
			//CONCENTRATION_NOMINAL[0] = iii*1e-4;              		        // carbon -C (fraction mole)
			//CONCENTRATION_NOMINAL[2] = jjj*1e-4;							// nickel -Ni (fraction mole)
			//CONCENTRATION_NOMINAL[1] = 1-CONCENTRATION_NOMINAL[2]-CONCENTRATION_NOMINAL[0];// nickel -Ni (fraction mole)	
			
			//x_i=CONCENTRATION_NOMINAL[0];
			//x_j=CONCENTRATION_NOMINAL[2];
			
	printf ( "Begin calcul the system FE-2.e-5CR-2.e-3C equilibre thermodynamique (3D): \n" );
	for ( int iterseg=1;iterseg<=nbseg;iterseg++ )
	{
		temps=seg[iterseg].date;
		dt= ( seg[iterseg+1].date-temps ) /seg[iterseg+1].pas;
		
		T=seg[iterseg].temperature;
		dT= dt * ( seg[iterseg+1].temperature-seg[iterseg].temperature )
		/ ( seg[iterseg+1].date-seg[iterseg].date );
		
		while ( temps<seg[iterseg+1].date )
		{
			bool prevExts=false;
			bool good;

			qh=qhull(dim);
			temps+=dt;
			T+=dT;
			
			/*if (( T>=1100. ) && ( T<=1771.)) continue;
			if (prevT==0)
			{prevT=T;}
			else if ((T-prevT)<=dT)
			{
			  prevExts=true;
			  prevT=T;
			}*/
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
			std::vector <double> Y1 ( nt.nn );

			nt.fvec.resize ( nt.nn );
			nt.fvec1.resize ( nt.nn ); //Jacobian
			printf ( "T = %lf \n", T );
			
			mat.Thermo->read ( T );
			cem.classe0->Thermo->read ( T );
			
			/*if ((prevExts) && (good))
			{
				Y=prevVector;
			}
			else*/
			{
				t_debut_qhull=clock();
				qh.discretization (T,55);
				qh.build_hull();
				vector<point> ext=qh.extremums();
				t_fin_qhull = clock();
				
				//if(ext.empty()) {printf("no extremums") continue};
				if ( ( !ext.empty() ) )
				{
					printf ( "extremums:" );
					for ( int i=0;i<ext.size();i++ )
						cout<<i+1<<":"<<ext[i]<<endl;
					printf("c = %e\n",CONCENTRATION_NOMINAL[0]/(1.-CONCENTRATION_NOMINAL[0]));
					printf("cr = %e\n",CONCENTRATION_NOMINAL[2]/(1.-CONCENTRATION_NOMINAL[0]));
				}
				else
				{
					printf ( "No extremums T = %lf\n",T );
					printf("c = %e\n",CONCENTRATION_NOMINAL[0]/(1.-CONCENTRATION_NOMINAL[0]));
					printf("cr = %e\n",CONCENTRATION_NOMINAL[2]/(1.-CONCENTRATION_NOMINAL[0]));

					continue;
				}
			
				cout<<"Current extremums:"<<ext[0]<<";"<<ext[1]<<endl;
				Y[1] = ( ext[0].func()?ext[1][1]:ext[0][1] );  //CR
				if ( Y[1]==0. )
				{
					Y[1]=1.e-4;
				}
				Y[0] = 1.-Y[1];       //FE
				Y[2] = ( ext[0].func()?ext[1][0]:ext[0][0] );  //C
				if ( Y[2]==0. )
				{
					Y[2]=1.e-4;
				}
				Y[3] = 1.-Y[2];       //VA
				Y[5] =( ext[0].func()?ext[0][1]:ext[1][1] );   //CR
				if ( Y[5]==0. )
				{
					Y[5]=1.e-4;
				}
				Y[4] = 1.0-Y[5];      //FE
				Y[6] = ( ext[0].func()?ext[0][0]:ext[1][0] );  //C
				if ( Y[6]==0.0 )
				{
					Y[6]=1.e-4;
				}
				Y[7] = 1.-Y[6];       //VA
				Y[8] = 0.45;				
			}

	printf ( "\n" );printf ( "Vector Y[i] before:" );
	for ( int k=0 ; k < 9 ; k++ )
	{printf ( "%g ",Y[k] );}printf ( "\n" );

	t_debut_NR = clock();			
	int check=0;
	good=nt.calc ( Y, &check,CONCENTRATION_NOMINAL,3,conv);
	t_fin_NR = clock() ;
	t_fin_main = clock();
	double ff= CLOCKS_PER_SEC;
			
			
	for ( int l = 0; l < 9; l++ )
	{
		if ( Y[l]<0.0 )
		{
			break;
		}
	}
			printf ( "\n" );printf ( "Vector Y[i] after:" );
			for ( int k=0 ; k < 9 ; k++ ){
			printf ( "%g ",Y[k] );}printf ( "\n" );
			
			//convergance FRACTION_SITE dans FRACTION_MOLAIRE
			std::vector <double> X(7);
			 double x_c_fcc,x_cr_fcc,x_c_bcc,x_cr_bcc;			
			 x_cr_bcc = Y[1] / (1. + 3. * Y[2]) ;
			 x_c_bcc  = 3. * Y[2] / (1. + 3. * Y[2]) ;
			 x_cr_fcc = Y[5] / (1. + Y[6]) ;
			 x_c_fcc  = Y[6] / (1. + Y[6]) ;			
			 
			 X[1] = x_cr_bcc;//cr			
			 X[2] = x_c_bcc;//c
			 X[0] = 1-X[1]-X[2];//fe
			 X[4] = x_cr_fcc;//cr			
			 X[5] = x_c_fcc;//c
			 X[3] = 1-X[4]-X[5]; //fe
			 X[6] = Y[8];
			
			/*printf("\n\nFIN\n\n");
			printf("\n\nTemps CPU : temps (s)	\
				   \n-QuickHull       : %e  %e s		\
				   \n-NR              : %e  %e s		\
				   \n-Total           : %e  %e s		\
				   \n",							\
				   (double) t_fin_qhull - (double) t_debut_qhull , ((double) t_fin_qhull - (double) t_debut_qhull) / ff ,	\
				   (double) t_fin_NR - (double) t_debut_NR , ((double) t_fin_NR - (double) t_debut_NR) / ff ,  \
				   (double) t_fin_main - (double) t_debut_main , ((double) t_fin_main - (double) t_debut_main) / ff ) ;*/
			
			
			
	if (good)
	{
		so.ecriture ( Y );
		so2.ecritureX ( X );
		//prevVector=Y;
		//prevT=T;
		if((T>=590.) && (T<=1100.))
			calc3cinetique1(CONCENTRATION_NOMINAL,Y,prevVectorcin,flage);

	}
	
			
			
			/***************************************Cinetique*************************************************************/

		
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
    so.fermeture();

}

bool calc3cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage)
{
	newt nt1;
	nt1=newt();
	
	printf("CALCUL CINETIQUE: \n");	
	int j = 1;
	while ( mat.produit[j]!=NULL )
	{
		nt1.func     =mat.produit[1]->frontiere[1];
		nt1.nrfuncv  =&interface::KTT;
		
		nt1.jacobian =mat.produit[1]->frontiere[1];
		nt1.jacobfunc=&interface::JKT;
		
		j++;
	}
	if ( nt1.func == NULL ) printf ( "newton.func = null\n" );
	if ( nt1.jacobian==NULL) printf ("jacobian.func=null\n");

	int cal,cal1;
	int n=3000000;
	std::vector <double> Rlong(n), t(n), Xca(n), Xfea(n),\
	Xcra(n),Xcb(n), Xfeb(n),Xcrb(n),XV(n), \
	Gtrf(n),Gfricf(n),G(n);
	
	float v_ini;
	int diff_step, check,compteur;
	bool conv = false;
	double raison,temps;
	
	
	FILE *fpp=fopen("coeff_cin.txt","r") ;
	fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
	fscanf(fpp , "%i  \n" , &diff_step) ;	
	fclose(fpp);

	
	check=0;
	nt1.fvec.clear();
	nt1.fvec1.clear();
	nt1.nn = 9;
	std::vector <double> Ycinetique ( nt1.nn );
	nt1.fvec.resize ( nt1.nn );
	nt1.fvec1.resize ( nt1.nn ); 

	//*************#######Parameters########*********
	v_ini=0.0;
	Rzin=1e-3;
	dc =  1.5e-5 * exp( -1.421e5 /R / T )*pow(1.e6,2.);
	dcr = 35.e-5 * exp( -2.860e5 /R/ T )*pow(1.e6,2.); 
	raison = pow ( dc / dcr , 1. / ( (double) diff_step ) ) ;
	compteur = 0;

	double dtempsFeCCr=0.1;
	double tempsFeCCr=0.0;
	cal=0;
	cal1=0;	
	
	double before=0.;
	double after=0.;	
//while(tempsFeCCr<=40000)	
 {	
	/*if (flage)
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
		 printf("before=%e \n",before);
		 after=Rzin; //printf("after=%e \n",after);
		 printf("rez=%e \n",(after-before)/before);
		 while((after-before)/before>1e-4)
		 {
		 dtempsFeCCr=dtempsFeCCr/10.;
		 Rzin=before+prevVectorcin[8]*dtempsFeCCr;
		 printf("dt=%e R=%e \n",dtempsFeCCr,Rzin);
		 after=Rzin;
		 printf("after-before/before=%e \n",(after-before)/after);
		 }
		 printf("dt= %e Rzin=%e \n",dtempsFeCCr,Rzin);	
		 tempsFeCCr=tempsFeCCr+dtempsFeCCr;		 		
	}
	else
	{	
		if(cal1==0){
		Ycinetique[0] = Y[0];			//FE
		Ycinetique[1] = 1.-Y[0];		//CR
		Ycinetique[2] = Y[2];			//C
		Ycinetique[3] = 1.-Y[2];        //VA
		Ycinetique[4] = Y[4];			//FE
		Ycinetique[5] =	1.0-Y[4];       //CR
		Ycinetique[6] = Y[6];			//C
		Ycinetique[7] = 1.-Y[6];        //VA
		Ycinetique[8] = v_ini;}
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
			//printf("before=%e \n",before);
			after=Rzin; //printf("after=%e \n",after);
			printf("rez=%e \n",(after-before)/before);
			while((after-before)/before>1e-4)
			{
				dtempsFeCCr=dtempsFeCCr/10.;
				Rzin=before+prevVectorcin[8]*dtempsFeCCr;
				printf("dt=%e R=%e \n",dtempsFeCCr,Rzin);
				after=Rzin;
				printf("after-before/before=%e \n",(after-before)/after);
			}
			printf("dt= %e Rzin=%e \n",dtempsFeCCr,Rzin);	
			tempsFeCCr=tempsFeCCr+dtempsFeCCr;	
		}

	}*/
	 Ycinetique[0] = Y[0];			//FE
	 Ycinetique[1] = 1.-Y[0];		//CR
	 Ycinetique[2] = Y[2];			//C
	 Ycinetique[3] = 1.-Y[2];        //VA
	 Ycinetique[4] = Y[4];			//FE
	 Ycinetique[5] =	1.0-Y[4];       //CR
	 Ycinetique[6] = Y[6];			//C
	 Ycinetique[7] = 1.-Y[6];        //VA
	 Ycinetique[8] = v_ini;
	 

	correc_D=0;
	if(flage)
	{	
		correc_D=dc/dcr;
		printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
		for ( int k=0 ; k < 9 ; k++ )
		{printf ( "%e ",Ycinetique[k] );}printf ( "\n" );
		
		nt1.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
		
		printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
		for ( int k=0 ; k < 9 ; k++ )
		{printf ( "%e ",Ycinetique[k] );}printf ( "\n" );
		
		printf("tempsFeCCr=%e \n",tempsFeCCr);

	}	
	else
	{	
	compteur=0;	
	while(compteur<= diff_step)
	{	
		
		correc_D = pow ( raison , compteur ) ;

			printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
			for ( int k=0 ; k < 9 ; k++ )
			{printf ( "%e ",Ycinetique[k] );}printf ( "\n" );
		
			nt1.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
		
			printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
			for ( int k=0 ; k < 9 ; k++ )
			{printf ( "%e ",Ycinetique[k] );}printf ( "\n" );
			
					
		printf("diff_step= %i \n",diff_step);
		printf("compt=%i \n",compteur);
		printf("correc_D=%e \n", correc_D);
		
		compteur ++;
	}}	


	//************########Register the resulats############*************
	if(conv)
	{	
		//convergance FRACTION_SITE dans FRACTION_MOLAIRE
		std::vector <double> Xcinetique(7);
		double x_c_fcc,x_cr_fcc,x_c_bcc,x_cr_bcc,x_fe_fcc,x_fe_bcc;			
		x_cr_bcc = Ycinetique[1] / (1. + 3. * Ycinetique[2]) ;
		x_c_bcc  = 3. * Ycinetique[2] / (1. + 3. * Ycinetique[2]) ;
		x_fe_bcc = Ycinetique[0]/(1.+3*Ycinetique[2]);
		x_cr_fcc = Ycinetique[5] / (1. + Ycinetique[6]) ;
		x_c_fcc  = Ycinetique[6] / (1. + Ycinetique[6]) ;		
		x_fe_fcc = Ycinetique[4]/(1.+Ycinetique[6]);
		
		
		double U_c_bcc,U_c_fcc,U_cr_bcc,U_cr_fcc;
		U_c_bcc=x_c_bcc/(1-x_c_bcc);U_c_fcc=x_c_fcc/(1-x_c_fcc);
		U_cr_bcc=x_cr_bcc/(1-x_c_bcc);U_cr_fcc=x_cr_fcc/(1-x_c_fcc);
		printf("Before NR : X(C,BCC) = %e  X(C,FCC)= %e \n",x_c_bcc,x_c_fcc);
		printf("            X(CR,BCC) = %e  X(CR,FCC)= %e \n\n",x_cr_bcc,x_cr_fcc);	
		printf("After NR: U(C,BCC) = %e	U(C,FCC) = %e \n",U_c_bcc,U_c_fcc);
		printf("	U(CR,BCC) = %e	U(CR,FCC) = %e \n",U_cr_bcc,U_cr_fcc);
		printf("	U(CR0) = %e	U(C0) = %e \n",CONCENTRATION_NOMINAL[2]/(1-CONCENTRATION_NOMINAL[0]), \
			   CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]));
		
		
		Xcinetique[1] = x_cr_bcc;//cr			
		Xcinetique[2] = x_c_bcc;//c
		Xcinetique[0] = 1-Xcinetique[1]-Xcinetique[2];//fe
		Xcinetique[4] = x_cr_fcc;//cr			
		Xcinetique[5] = x_c_fcc;//c
		Xcinetique[3] = 1-Xcinetique[4]-Xcinetique[5]; //fe
		Xcinetique[6] = Ycinetique[8];
		
		
		prevVectorcin[0]=Ycinetique[0];//fe
		prevVectorcin[1]=Ycinetique[1];//cr
		prevVectorcin[2]=Ycinetique[2];//c
		prevVectorcin[3]=Ycinetique[3];//va
		prevVectorcin[4]=Ycinetique[4];//fe
		prevVectorcin[5]=Ycinetique[5];//cr
		prevVectorcin[6]=Ycinetique[6];//c
		prevVectorcin[7]=Ycinetique[7];//va
		prevVectorcin[8]=Ycinetique[8];
		
		printf("\n*******Convergence\n");
		so1.ecriture (Ycinetique);
		so3.ecritureX (Xcinetique);

		 Rlong[cal]=Rzin; t[cal]=tempsFeCCr;
		 Xca[cal]=x_c_bcc/(1-x_c_bcc);Xfea[cal]=x_fe_bcc/(1-x_c_bcc);
		 Xcra[cal]=x_cr_bcc/(1-x_c_bcc);
		 Xcb[cal]=x_c_fcc/(1-x_c_fcc);Xfeb[cal]=x_fe_fcc/(1-x_c_fcc);
		 Xcrb[cal]=x_cr_fcc/(1-x_c_fcc);
		 XV[cal]=Xcinetique[6];
		 cal++;
		 printf("tempsFeCCr=%e dt=%e \n",tempsFeCCr,dtempsFeCCr);
		 printf("cal=%d \n",cal);
		
		flage = true;
	}
	else{
		printf("\n****Pas de convergence*******\n");printf("tempsFeCCr=%e \n",tempsFeCCr);
		flage = false;/*break;*/}
	 cal1++;
 }	 
	// t.resize(cal);Rlong.resize(cal);Xca.resize(cal);
	// Xfea.resize(cal);Xcra.resize(cal);
	// Xcb.resize(cal);Xfeb.resize(cal);Xcrb.resize(cal);
	// XV.resize(cal);Gtrf.resize(cal);Gfricf.resize(cal);
	// ecriture_growth_dissolution_FeCCr(t,Rlong,Xca,Xfea,Xcra,Xcb,Xfeb,Xcrb,XV,Gtrf,Gfricf,G);
	
	//###################################
	if((T>=800.) && (T<=1110.)){
		printf("              Calcul paraequilibre \n");
		bool flage_para = false;
		std::vector<double> prevVectorcin_para(nt1.nn);
		calc3cinetique_para(CONCENTRATION_NOMINAL,Y,prevVectorcin_para,flage_para);}
	
}	

bool calc3cinetique_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage)
{
	newt nt2;
	nt2=newt();
	integer itt2;
	itt2=integer();
	

	int f = 1;
	while ( mat.produit[f]!=NULL )
	{
		nt2.func     =mat.produit[1]->frontiere[1];
		nt2.nrfuncv  =&interface::KTT_prof;
		
		nt2.jacobian =mat.produit[1]->frontiere[1];
		nt2.jacobfunc=&interface::JKT_prof;
		
		itt2.funcin	 =mat.produit[1]->frontiere[1];
		itt2.nrfuncin =&interface::funcFeCrC;
		
		itt2.funcpot  =mat.produit[1]->frontiere[1];
		itt2.nrfuncpot=&interface::funconeFeCrC;
		
		itt2.funcderivFeCrC = mat.produit[1]->frontiere[1];
		itt2.nrfuncderivFeCrC = &interface::funcderivmuFeCrC;
		
		f++;
	}
	if ( nt2.func == NULL ) printf ( "newton.func = null\n" );
	if ( nt2.jacobian==NULL) printf ("jacobian.func=null\n");
	if ( itt2.funcin == NULL)printf("integer.funcin=null\n");

	int n=3000000;
	std::vector <double> Rlong(n), t(n), Xca(n), Xfea(n),\
	Xcra(n),Xcb(n), Xfeb(n),Xcrb(n),XV(n), \
	Gtrf(n),Gfricf(n),G(n);

	
	double dtempsFeCCr, tempsFeCCr;
	double x_c_fcc,x_c_bcc,x_fe_fcc,x_fe_bcc, \
			x_cr_bcc,x_cr_fcc;

	double U_c_fcc,U_c_bcc,U_fe_fcc,U_fe_bcc, \
	U_cr_bcc,U_cr_fcc;
	
	double OMEGA_c,OMEGA_cr,SURSAT_c,SURSAT_cr;
	vector <double> Gtr(n),Gmuc(n),Gmucr(n), \
	Gfric(n),Gx_c_bcc(n),Gx_cr_bcc(n),Gx_c_fcc(n),Gx_cr_fcc(n), \
	GtrV(n),taille_t(n),Gcomp(n),Gx_c_bcc_eq(n),Gx_cr_bcc_eq(n), \
	Gx_c_fcc_eq(n),Gx_cr_fcc_eq(n),GVeq(n), \
	G_c_nom(n),G_cr_nom(n);
	
	vector <double> L(n),B(n),Fx_c_bcc(n),Fx_cr_bcc(n), \
	Fx_c_fcc(n),Fx_cr_fcc(n);
	
	
	float v_ini;
	int diff_step,tr,tr1,opt,check,cal1,cal;
	int n1=100000;
	double Vopt[n1];
	double param_cin;
	bool conv_para = false;
	
	FILE *fpp=fopen("coeff_cin.txt","r") ;
	fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
	fscanf(fpp , "%i  \n" , &diff_step) ;	
	fclose(fpp);
	
	
	check=0;
	nt2.fvec.clear();
	nt2.fvec1.clear();
	nt2.nn = 5;
	std::vector <double> Ycinetique_para ( nt2.nn );
	nt2.fvec.resize ( nt2.nn );
	nt2.fvec1.resize ( nt2.nn ); 
	
	
	
	//=========================Profile===============//
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
	//======================================================//
	

	dc =  1.5e-5 * exp( -1.421e5 /R / T )*pow(1.e6,2.);
	dcr = 35.e-5 * exp( -2.860e5 /R/ T )*pow(1.e6,2.); 
	double raison = pow ( dc / dcr , 1. / ( (double) diff_step ) ) ;
	int compteur = 0;
	
	taille_gradient=1.e-3;
	Rzin=1.e-2;
	M=0.035*exp(-17700/T)*pow(1e6,4);
	//M=1e14;
	Vm=1e13;	
	std::vector <double> Ycinetique_para_profile ( 9 );
	
	//printf("dc=%e dcr=%e \n",dc,dcr);
	cal=0;
	cal1=0;
	double before=0.;
	double after=0.;
	
	dtempsFeCCr=0.1;
	tempsFeCCr=0.;
	
	
	
	//coeff pour coeff de diffusion à l'interface
	coeff_n=1e13; //(Dc)
	coeff_k=1e-3; //(Dcr)
	Lcoef=0.;Bcoef=0.;

	

	
	//int n2 = 20;
	//double step=20/n2;	
	//for(double ii=1;ii<=n2;ii++)
	//	for(double jj=1;jj<=n2;jj++)
	//	{
	//		coeff_n=ii*step;
	//		coeff_k=jj*step;
	
	//int n2 = 100;
	//double step=100/n2;	
	//for(double ii=-100;ii<=n2;ii+=50)
	//	for(double jj=-100;jj<=n2;jj+=50)
	//	{		
	//		Lcoef=ii*step;
	//		Bcoef=jj*step;

	//double n2 = 50;
	//double step=50/n2;
	//for(float ii=1e-3;ii<n2;ii+=1e-3)
	
	//while(tempsFeCCr<=60000)
	{
	
		if(flage)
		{
			dtempsFeCCr=0.1;			
			Ycinetique_para[0]=prevVectorcin[0];//c_alpha	
			Ycinetique_para[1]=prevVectorcin[1];//cr_alpha	
			Ycinetique_para[2]=prevVectorcin[2];//c_gamma
			Ycinetique_para[3]=prevVectorcin[3];//cr_gamma
			Ycinetique_para[4]=prevVectorcin[4];
			VIT_TEMPORAIRE	= Ycinetique_para[4];
			before=Rzin;
			Rzin=Rzin+prevVectorcin[4]*dtempsFeCCr;
			//printf("before=%e \n",before);
			after=Rzin; //printf("after=%e \n",after);
			printf("rez=%e \n",(after-before)/before);
			while((after-before)/before>1e-1)
			{
				dtempsFeCCr=dtempsFeCCr/10.;
				Rzin=before+prevVectorcin[4]*dtempsFeCCr;
				printf("dt=%e R=%e \n",dtempsFeCCr,Rzin);
				after=Rzin;
				printf("after-before/before=%e \n",(after-before)/after);
			}
			printf("dt= %e \n",dtempsFeCCr);	
			tempsFeCCr=tempsFeCCr+dtempsFeCCr;
			//Rzin=ii*step;
		}	
		else
		{
			Ycinetique_para[0]=Y[2];//c_alpha	
			Ycinetique_para[1]=Y[1];//cr_alpha	
			Ycinetique_para[2]=Y[6];//c_gamma
			Ycinetique_para[3]=Y[5];//cr_gamma
			Ycinetique_para[4]=1e-9;//Y[8];
			VIT_TEMPORAIRE	= Ycinetique_para[4];
		}
	
		tr=0;

		if(flage)
		{
			correc_D=dc/dcr;

			printf ( "\n" );printf ( "Vector Ycinetique_para[i] before NR:" );
			for ( int k=0 ; k < 5 ; k++ )
			{printf ( "%e ",Ycinetique_para[k] );}printf ( "\n" );
			x_c_bcc=3*Ycinetique_para[0]/(1.+3.*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[2]/(1.+Ycinetique_para[2]);
			x_cr_bcc=Ycinetique_para[1]/(1.+3*Ycinetique_para[0]);
			x_cr_fcc=Ycinetique_para[3]/(1.+Ycinetique_para[2]);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);U_c_fcc=x_c_fcc/(1-x_c_fcc);
			U_cr_bcc=x_cr_bcc/(1-x_c_bcc);U_cr_fcc=x_cr_fcc/(1-x_c_fcc);
			printf("Before NR : X(C,BCC) = %e  X(C,FCC)= %e \n",x_c_bcc,x_c_fcc);
			printf("            X(CR,BCC) = %e  X(CR,FCC)= %e \n\n",x_cr_bcc,x_cr_fcc);	
			printf("After NR: U(C,BCC) = %e	U(C,FCC) = %e \n",U_c_bcc,U_c_fcc);
			printf("	U(CR,BCC) = %e	U(CR,FCC) = %e \n",U_cr_bcc,U_cr_fcc);
			printf("	U(CR0) = %e	U(C0) = %e \n",CONCENTRATION_NOMINAL[2]/(1-CONCENTRATION_NOMINAL[0]), \
				   CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]));
			
			
			nt2.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
			
			printf ( "\n" );printf ( "Vector Ycinetique_para[i] after NR:" );
			for ( int k=0 ; k < 5 ; k++ )
				   {printf ( "%e ",Ycinetique_para[k] );}printf ( "\n" );
			x_c_bcc=3*Ycinetique_para[0]/(1.+3.*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[2]/(1.+Ycinetique_para[2]);
			x_cr_bcc=Ycinetique_para[1]/(1.+3*Ycinetique_para[0]);
			x_cr_fcc=Ycinetique_para[3]/(1.+Ycinetique_para[2]);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);U_c_fcc=x_c_fcc/(1-x_c_fcc);
			U_cr_bcc=x_cr_bcc/(1-x_c_bcc);U_cr_fcc=x_cr_fcc/(1-x_c_fcc);
			for ( int k=0 ; k < 5 ; k++ )
			{printf ( "%e ",Ycinetique_para[k] );}printf ( "\n" );
			printf("Before NR : X(C,BCC) = %e  X(C,FCC)= %e \n",x_c_bcc,x_c_fcc);
			printf("            X(CR,BCC) = %e  X(CR,FCC)= %e \n\n",x_cr_bcc,x_cr_fcc);	
			printf("After NR: U(C,BCC) = %e	U(C,FCC) = %e \n",U_c_bcc,U_c_fcc);
			printf("	U(CR,BCC) = %e	U(CR,FCC) = %e \n",U_cr_bcc,U_cr_fcc);
			printf("	U(CR0) = %e	U(C0) = %e \n",CONCENTRATION_NOMINAL[2]/(1-CONCENTRATION_NOMINAL[0]), \
			CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]));
			printf("tempsFeCCr = %e \n",tempsFeCCr);
			printf("L= %e; B=%e \n",Lcoef,Bcoef);


			
		}
		else{
		compteur=0;	
		while(compteur<= diff_step)
		{
			correc_D = pow ( raison , compteur ) ;
			
			printf ( "\n" );printf ( "Vector Ycinetique_para[i] before NR:" );
			for ( int k=0 ; k < 5 ; k++ )
			{printf ( "%e ",Ycinetique_para[k] );}printf ( "\n" );
			x_c_bcc=3*Ycinetique_para[0]/(1.+3.*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[2]/(1.+Ycinetique_para[2]);
			x_cr_bcc=Ycinetique_para[1]/(1.+3*Ycinetique_para[0]);
			x_cr_fcc=Ycinetique_para[3]/(1.+Ycinetique_para[2]);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);U_c_fcc=x_c_fcc/(1-x_c_fcc);
			U_cr_bcc=x_cr_bcc/(1-x_c_bcc);U_cr_fcc=x_cr_fcc/(1-x_c_fcc);
			printf("Before NR : X(C,BCC) = %e  X(C,FCC)= %e \n",x_c_bcc,x_c_fcc);
			printf("            X(CR,BCC) = %e  X(CR,FCC)= %e \n\n",x_cr_bcc,x_cr_fcc);	
			printf("After NR: U(C,BCC) = %e	U(C,FCC) = %e \n",U_c_bcc,U_c_fcc);
			printf("	U(CR,BCC) = %e	U(CR,FCC) = %e \n",U_cr_bcc,U_cr_fcc);
			printf("	U(CR0) = %e	U(C0) = %e \n",CONCENTRATION_NOMINAL[2]/(1-CONCENTRATION_NOMINAL[0]), \
			CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]));		
		
			nt2.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
		
			printf ( "\n" );printf ( "Vector Ycinetique_para[i] after NR:" );
			for ( int k=0 ; k < 5 ; k++ )
			{printf ( "%e ",Ycinetique_para[k] );}printf ( "\n" );
			x_c_bcc=3*Ycinetique_para[0]/(1.+3.*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[2]/(1.+Ycinetique_para[2]);
			x_cr_bcc=Ycinetique_para[1]/(1.+3*Ycinetique_para[0]);
			x_cr_fcc=Ycinetique_para[3]/(1.+Ycinetique_para[2]);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);U_c_fcc=x_c_fcc/(1-x_c_fcc);
			U_cr_bcc=x_cr_bcc/(1-x_c_bcc);U_cr_fcc=x_cr_fcc/(1-x_c_fcc);
			for ( int k=0 ; k < 5 ; k++ )
			{printf ( "%e ",Ycinetique_para[k] );}printf ( "\n" );
			printf("Before NR : X(C,BCC) = %e  X(C,FCC)= %e \n",x_c_bcc,x_c_fcc);
			printf("            X(CR,BCC) = %e  X(CR,FCC)= %e \n\n",x_cr_bcc,x_cr_fcc);	
			printf("After NR: U(C,BCC) = %e	U(C,FCC) = %e \n",U_c_bcc,U_c_fcc);
			printf("	U(CR,BCC) = %e	U(CR,FCC) = %e \n",U_cr_bcc,U_cr_fcc);
			printf("	U(CR0) = %e	U(C0) = %e \n",CONCENTRATION_NOMINAL[2]/(1-CONCENTRATION_NOMINAL[0]), \
			CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]));
			
			
			printf("\n\n");
			printf("diff_step= %i \n",diff_step);
			printf("compt=%i \n",compteur);
			printf("correc_D=%e \n", correc_D);
			printf("L= %e; B=%e \n",Lcoef,Bcoef);

			
		compteur ++;
		}}


		
		
	if(conv_para)
	{
		/*
		printf("\n      DIFFERENTIAL EQUATIONS WITH P VARIABLES OF ORDER 1\n");
		printf("            of type yi' = f(y1,y2,...,yn), i=1..n\n\n");

		VIT_TEMPORAIRE=Ycinetique_para[4];
		double xc=3. * Ycinetique_para[0] / (1. + 3. * Ycinetique_para[0]);
		double xcr=1. * Ycinetique_para[1] / (1. + 3. * Ycinetique_para[0]);
		//double xc=3. * Y[2] / (1. + 3. * Y[2]);
		//double xcr=1. * Y[1] / (1. + 3. * Y[2]);
		
		yi[0]=xc/(1.-xc);
		yi[1]=xcr/(1.-xc);
		x_c0	= yi[0];
		x_cr0	= yi[1];	
		
	
		printf("X(C,BCC) = %e  X(CR,BCC) = %e V = %e correct_D=%e taille= %e\n",xc,xcr,VIT_TEMPORAIRE,correc_D, \
			   taille_gradient);
		printf("U(C,BCC) = %e  U(CR,BCC) = %e \n",yi[0],yi[1]);
		itt2.Eqdifp(v1,v2,v3,v4,xi,xf,yi,ndata,p,fi,CONCENTRATION_NOMINAL,Ycinetique_para_profile,P,P1);
		printf ( "\n" );printf ("End Runga-Kutta\n");
		printf("\n*******Convergence growth diffusion\n");*/
		
		std::vector <double> Xcinetique(7);
		x_c_bcc  = 3. * Ycinetique_para[0] / (1. + 3. * Ycinetique_para[0]);
		x_cr_bcc = Ycinetique_para[1]/(1.+3.*Ycinetique_para[0]);
		x_fe_bcc = 1.-x_c_bcc-x_cr_bcc;
		x_c_fcc  = Ycinetique_para[2] / (1. + Ycinetique_para[2]);
		x_cr_fcc = Ycinetique_para[3]/(1.+Ycinetique_para[2]);
		x_fe_fcc = 1.-x_c_fcc-x_cr_fcc;
		Xcinetique[0] = x_c_bcc;
		Xcinetique[1] = x_cr_bcc;
		Xcinetique[2] = x_fe_bcc;
		Xcinetique[3] = x_c_fcc;
		Xcinetique[4] = x_cr_fcc;
		Xcinetique[5] = x_fe_fcc;
		Xcinetique[6] = Ycinetique_para[4];
		so5.ecritureX ( Xcinetique );	
		
		U_c_bcc=x_c_bcc/(1-x_c_bcc);
		U_fe_bcc=x_fe_bcc/(1-x_c_bcc);
		U_cr_bcc=x_cr_bcc/(1-x_c_bcc);
		U_c_fcc=x_c_fcc/(1-x_c_fcc);
		U_fe_fcc=x_fe_fcc/(1-x_c_fcc);
		U_cr_fcc=x_cr_fcc/(1-x_c_fcc);
		
		std::vector <double> Ucinetique(7);
		Ucinetique[0] = U_c_bcc;
		Ucinetique[1] = U_cr_bcc;
		Ucinetique[2] = U_fe_bcc;
		Ucinetique[3] = U_c_fcc;
		Ucinetique[4] = U_cr_fcc;
		Ucinetique[5] = U_fe_fcc;
		Ucinetique[6] = Ycinetique_para[4];
		so4.ecritureX ( Ucinetique );
		
		
		OMEGA_c=(U_c_fcc-(CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0])))/(U_c_fcc-x_c_bcc);
		SURSAT_c=2.*Rzin*(1.-OMEGA_c)/OMEGA_c;
		OMEGA_cr=(x_cr_fcc-(CONCENTRATION_NOMINAL[2]/(1-CONCENTRATION_NOMINAL[0])))/(x_cr_fcc-x_cr_bcc);
		SURSAT_cr=2.*Rzin*(1.-OMEGA_cr)/OMEGA_cr;
		
		printf("muc=%e mucr=%e gtr=%e\n",muc,mucr,gtr);
		Gtr[tr]=gtr;Gmuc[tr]=muc;Gmucr[tr]=mucr;
		Gfric[tr]=Ycinetique_para[2]*(Vm/M)/(R*T);
		
		Gx_c_bcc[tr]=U_c_bcc;Gx_cr_bcc[tr]=U_cr_bcc;
		Gx_c_fcc[tr]=U_c_fcc;Gx_cr_fcc[tr]=U_cr_fcc;
		GtrV[tr]=Ycinetique_para[4];	
		taille_t[tr]=taille_gradient;
		Gcomp[tr]=tr;
		
		x_c_bcc  = 3. * Y[2] / (1. + 3. * Y[2]);
		x_cr_bcc = Y[1]/(1.+3.*Y[2]);
		x_fe_bcc = 1.-x_c_bcc-x_cr_bcc;
		x_c_fcc  = Y[6] / (1. + Y[6]);
		x_cr_fcc = Y[5]/(1.+Y[6]);
		x_fe_fcc = 1.-x_c_fcc-x_cr_fcc;
		Gx_c_bcc_eq[tr]=x_c_bcc/(1-x_c_bcc);
		Gx_cr_bcc_eq[tr]=x_cr_bcc/(1-x_c_bcc);
		Gx_c_fcc_eq[tr]=x_c_fcc/(1-x_c_fcc);
		Gx_cr_fcc_eq[tr]=x_cr_fcc/(1-x_c_fcc);
		GVeq[tr]=Y[8];
		G_c_nom[tr]=CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]);
		G_cr_nom[tr]=CONCENTRATION_NOMINAL[2]/(1-CONCENTRATION_NOMINAL[0]);
		tr++;
		Gtr.resize(tr);Gmuc.resize(tr);Gmucr.resize(tr); 
		Gfric.resize(tr);Gx_c_bcc.resize(tr);Gx_cr_bcc.resize(tr);
		Gx_c_fcc.resize(tr);Gx_cr_fcc.resize(tr);
		GtrV.resize(tr);taille_t.resize(tr);Gcomp.resize(tr);
		Gx_c_bcc_eq.resize(tr);Gx_cr_bcc_eq.resize(tr); 
		Gx_c_fcc_eq.resize(tr);Gx_cr_fcc_eq.resize(tr);GVeq.resize(tr); 
		G_c_nom.resize(tr);G_cr_nom.resize(tr);
		
		/*
		if((cal1==0)||(cal1==11000)||(cal1==15000)||(cal1==20000)||(cal1==35000)||(cal1==60000)||(cal1==110000)||(cal1==150000)||(cal1==210000))
		{
		ecriture3D(tempsFeCCr,Gtr,Gmuc,Gmucr,Gfric,Gx_c_bcc,Gx_cr_bcc, \
				Gx_c_fcc,Gx_cr_fcc,GtrV,taille_t,Gcomp,Gx_c_bcc_eq, \
				Gx_cr_bcc_eq,Gx_c_fcc_eq,Gx_cr_fcc_eq,GVeq, \
				G_c_nom,G_cr_nom);
		ecriture_profileC_FeCCr(tempsFeCCr,Rzin,SURSAT_c,SURSAT_cr,ndata,v1,v2,xi,xf);
		ecriture_profilePOT_FeCCr(tempsFeCCr,Rzin,ndata,P,xi,xf);
		}*/
		prevVectorcin[0]=Ycinetique_para[0];
		prevVectorcin[1]=Ycinetique_para[1];
		prevVectorcin[2]=Ycinetique_para[2];
		prevVectorcin[3]=Ycinetique_para[3];
		prevVectorcin[4]=Ycinetique_para[4];
		
		
		 Rlong[cal]=Rzin; t[cal]=tempsFeCCr;
		 Xca[cal]=Ucinetique[0];Xfea[cal]=Ucinetique[2];
		 Xcra[cal]=Ucinetique[1];
		 Xcb[cal]=Ucinetique[3];Xfeb[cal]=Ucinetique[5];
		 Xcrb[cal]=Ucinetique[4];
		 XV[cal]=Ucinetique[6];
		 Gtrf[cal]=muc;Gfricf[cal]=mucr;G[cal]=gtr;
		

		
		/*L[cal] = dc*coeff_n;B[cal]=dcr*coeff_k;
		Fx_c_bcc[cal]=U_c_bcc;Fx_cr_bcc[cal]=U_cr_bcc;
		Fx_c_fcc[cal]=U_c_fcc;Fx_cr_fcc[cal]=U_cr_fcc;*/
		
		/*L[cal] = Lcoef;B[cal]=Bcoef;
		Fx_c_bcc[cal]=U_c_bcc;Fx_cr_bcc[cal]=U_cr_bcc;
		Fx_c_fcc[cal]=U_c_fcc;Fx_cr_fcc[cal]=U_cr_fcc;*/
		
		/*
		 L[cal] = Rzin;B[cal]=Ycinetique_para[4];
		 Fx_c_bcc[cal]=U_c_bcc;Fx_cr_bcc[cal]=U_cr_bcc;
		 Fx_c_fcc[cal]=U_c_fcc;Fx_cr_fcc[cal]=U_cr_fcc;*/
		
		cal++;
		flage=true;
	}
	else{printf("No convergence growth diffusion and the second try");/*break;*/}
		
		printf("cal1=%d \n",cal1);
		cal1++;
	}
/*
	L.resize(cal);B.resize(cal);Fx_c_bcc.resize(cal);Fx_cr_bcc.resize(cal);
	Fx_c_fcc.resize(cal);Fx_cr_fcc.resize(cal);
	cryptograph(L,B,Fx_c_bcc,Fx_cr_bcc,Fx_c_fcc,Fx_cr_fcc);*/
	
	
	 t.resize(cal);Rlong.resize(cal);Xca.resize(cal);
	 Xfea.resize(cal);Xcra.resize(cal);
	 Xcb.resize(cal);Xfeb.resize(cal);Xcrb.resize(cal);
	 XV.resize(cal);Gtrf.resize(cal);Gfricf.resize(cal);
	 ecriture_growth_dissolution_FeCCr(t,Rlong,Xca,Xfea,Xcra,Xcb,Xfeb,Xcrb,XV,Gtrf,Gfricf,G);

	
}


bool calc3cinetique2(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage)
{
	newt nt1;
	
	nt1=newt();
	int j = 1;
	while ( mat.produit[j]!=NULL )
	{
		nt1.func     =mat.produit[1]->frontiere[1];
		nt1.nrfuncv  =&interface::KTT1;
		
		nt1.jacobian =mat.produit[1]->frontiere[1];
		nt1.jacobfunc=&interface::JKT1;
		j++;
	}
	if ( nt1.func == NULL ) printf ( "newton.func = null\n" );
	if ( nt1.jacobian==NULL) printf ("jacobian.func=null\n");
	
	float v_ini;
	int diff_step;
	bool conv = false;
	
	FILE *fpp=fopen("coeff_cin.txt","r") ;
	fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
	fscanf(fpp , "%i  \n" , &diff_step) ;	
	fclose(fpp);
	int check=0;
	nt1.fvec.clear();
	nt1.fvec1.clear();
	nt1.nn = 9;
	std::vector <double> Ycinetique ( nt1.nn );
	nt1.fvec.resize ( nt1.nn );
	nt1.fvec1.resize ( nt1.nn ); 
	
	/*if (conv)
	{
		if(prevVectorcin[8]<=0.0){prevVectorcin[8]=1;}
		Ycinetique=prevVectorcin;
		//taille_gradient_alpha=taille_gradient_alpha+prevVectorcin[8]*dtemps;
		//taille_gradient_betta=taille_gradient_betta+prevVectorcin[8]*dtemps;
	}
	else*/
	{
		Ycinetique[0] = Y[0];  //FE
		Ycinetique[1] = 1.-Y[0];       //CR
		Ycinetique[2] = Y[2];  //C
		Ycinetique[3] = 1.-Y[2];       //VA
		Ycinetique[4] = Y[4];			//FE
		Ycinetique[5] = 1.0-Y[4];      //CR
		Ycinetique[6] = Y[6];  //C
		Ycinetique[7] = 1.-Y[6];       //VA
		Ycinetique[8] = v_ini;
	}	
	
	
	
	dc =  1.5e-5 * exp( -1.421e5 /R / T ) ;
	dcr = 35.e-5 * exp( -2.860e5 /R/ T ) ; 
	double raison = pow ( dc / dcr , 1. / ( (double) diff_step ) ) ;
	
	int compteur = 0;
	while(compteur <= diff_step)
	{
		correc_D = pow ( raison , compteur ) ;
		printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
		for ( int k=0 ; k < 9 ; k++ )
		{printf ( "%e ",Ycinetique[k] );}printf ( "\n" );
		nt1.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
		printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
		for ( int k=0 ; k < 9 ; k++ )
		{printf ( "%e ",Ycinetique[k] );}printf ( "\n" );
		printf("diff_step= %i \n",diff_step);
		printf("compt=%i \n",compteur);
		compteur ++;
	}	
	
	if(conv)
	{
		prevVectorcin=Ycinetique;
		printf("\n*******Convergence\n");
		so1.ecriture ( Ycinetique );
		flage = true;
		return flage;
	}
	
	printf("\n****Pas de convergence*******\n");
	flage = false;
	return flage;
}

