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
double M;
double Vm;
double diffFeC;
extern double x_c[10];
double x_c0;
double x_c01;
double fricFeC;
double VIT_TEMPORAIRE;
double Rzin;
double y_c_Svoboda;
extern double gtr;
extern double muc;


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


extern double correc_D;
extern float CIN;



void calc2dequilibre();
bool calc2cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc2cinetique2(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &falge);
bool calc2cinetique_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &falge);
extern void ecriture2D(double temps,std::vector <double> &diff_d, \
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
extern void ecriture2D_(double temps, double Rzin,std::vector <double> &diff_d, \
				 std::vector <double> &diff,std::vector <double> &del, \
				 std::vector <double> &omega, \
				 std::vector <double> &Grka, \
				 std::vector <double> &Gnra,std::vector <double> &Gnr, \
				 std::vector <double> &Grk, std::vector <double> &Gv, \
				 std::vector <double> &Gtr, \
				 std::vector <int> &Gcompt_);
extern void ecriture_growth_dissolution_FeC(std::vector <double> &t,std::vector <double> &Rlong, \
									 std::vector <double> &Xca, \
									 std::vector <double> &Xfea, \
									 std::vector <double> &Xcb, \
									 std::vector <double> &Xfeb, \
									 std::vector <double> &XV, \
									 std::vector <double> &Gtrf, \
									 std::vector <double> &Gfricf); 
extern void ecriture_profileC_FeC(double temps,double Rzin,double SURSAT,int m, double *tx, \
						   double *ty);
extern void ecriture_profilePOT_FeC(double temps,double Rzin,int m, double *tx, std::vector < std::vector<double> > &P);


double regulafalsi(double a,double b,double x_c_gamma,std::vector <double> CONCENTRATION_NOMINAL);
double f_vitess_x1(double V, double x_c_gamma,std::vector <double> CONCENTRATION_NOMINAL);
double f_vitess_x2(double V, double x_c_gamma,std::vector <double> CONCENTRATION_NOMINAL);
double regulafalsi(double a,double b,double V,std::vector <double> CONCENTRATION_NOMINAL);

extern clock_t t_debut_main , t_debut_qhull , t_fin_qhull , t_debut_NR , t_fin_NR , t_fin_main ;


void calc2dequilibre()
{

def=2;
int dim = 2;
	newt nt;
	qhull qh(dim);
	
	double prevT=0;
	bool conv = false;
	bool flage= false;
	
	vector<double> prevVector(6);
	std::vector<double> prevVectorcin(7);	
	for(int i=0;i<7;i++){prevVectorcin[i]=0;}	
	for(int i=0;i<6;i++){prevVector[i]=0;}

	
	std::vector <double> CONCENTRATION_NOMINAL ( 2 );
	CONCENTRATION_NOMINAL[0] = 5.e-3;// carbon -C
	//4.4859334e-2(1.e-2);//9.2319223e-3(2.e-3);//8.6666998e-2(2.e-2);0.53755378 (2.e-1);4.6327466e-3(1.e-3)
	//2.3205926E-3 (5.e-4); 2.2831635e-2(0.5e-2)
	CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0];// iron -FE
	
	fstream file_v;
	char filename[30];
	sprintf ( filename,"CPU2d.txt");
	file_v.open ( filename,fstream::out );
	if ( file_v.fail() )
		cout<<"Error opening file"<<endl;
	
	double xx = 0;
	int N;
	//for(N=100;N<200;N+=2)
	//for(xx=0; xx<0.69; xx+=0.01)
	{
		t_debut_main = clock();
		double counter=time(0);
		
		//double Mc = 12.011;
		//double Mfe= 55.847;
		//double summ = xx/Mc + (1.-xx)/Mfe;
		//CONCENTRATION_NOMINAL[0] = xx/Mc/summ;
		//CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0];
	
	printf ( "Begin calcul the system FE - 2.e-3C equilibre thermodynamique (2D): \n" );
	for ( int iterseg=1;iterseg<=nbseg;iterseg++ )
	{
		
		temps=seg[iterseg].date;
		dt= ( seg[iterseg+1].date-temps ) /seg[iterseg+1].pas;

		T=seg[iterseg].temperature;
		dT= dt * ( seg[iterseg+1].temperature-seg[iterseg].temperature )
		/ ( seg[iterseg+1].date-seg[iterseg].date );
		while ( temps<seg[iterseg+1].date )
		{
			bool good;
			bool prevExts=false;
			
			qh=qhull(dim);
			temps+=dt;
			T+=dT;
			
			
			if (prevT==0)prevT=T;
			else if ((T-prevT)<=dT)
			{
				prevExts=true;
				prevT=T;
			}
			
			nt=newt();
			int i = 1;
			while ( mat.produit[i]!=NULL )
			{
				nt.func     =mat.produit[1]->frontiere[1];
				nt.nrfuncv  =&interface::SB;
				
				nt.jacobian =mat.produit[1]->frontiere[1];
				nt.jacobfunc=&interface::JB;
				
				qh.func1    =mat.produit[1]->frontiere[1]; 
				qh.nrfuncv1 =&interface::phase1_binaire;   
				
				qh.func2    =mat.produit[1]->frontiere[1]; 
				qh.nrfuncv2 =&interface::phase2_binaire;   
				i++;
			}
			if ( nt.func == NULL ) printf ( "newton.func = null\n" );
			if ( qh.func1== NULL ) printf ( "qhull.func = null\n" );  //
			
			
			mat.Thermo->read (T);
			cem.classe0->Thermo->read (T);
			
						
			nt.fvec.clear();
			nt.fvec1.clear();
			nt.nn = 6;
			std::vector <double> Y ( nt.nn );
			nt.fvec.resize ( nt.nn );
			nt.fvec1.resize ( nt.nn ); 
			int check = 0;
			
			/*if ((prevExts) && (good))
			{
				Y=prevVector;
			}
			
			else*/
			{
				t_debut_qhull=clock();
				N=100;
				qh.discretization (T,N/*25*/);
				qh.build_hull();
				t_fin_qhull = clock();

				
				vector<point> ext=qh.extremums();
				if ( ( !ext.empty() ) )
				{
					cout<<ext[0]<<endl<<ext[1]<<endl;;
				}
				
				else
				{
					printf ( "No extremums T = %lf\n",T );
					continue;
					//return 0;
				}
				
				Y[0] = 1.; //fe
				Y[1] = ext[0].func()==0?ext[0][0]:ext[1][0];//c
				if ( Y[1]>=0.0 )
				{
					Y[1]=1.e-8;
				}
				Y[2] = 1.-Y[1];//va
				Y[3] = 1.;//fe
				Y[4] = ext[0].func()==1?ext[0][0]:ext[1][0];//c
				if ( Y[4]<=0.0 )
				{
					Y[4]=1.e-8;
				}
				Y[5] = 1.-Y[4];//va
			}


			
			 
			printf ( "\n" );printf ( "Vector Y[i] before:" );
			for ( int k=0 ; k < 6 ; k++ ){printf ( "%g ",Y[k] );}printf ( "\n" );
			t_debut_NR = clock();			
			good=nt.calc ( Y, &check,CONCENTRATION_NOMINAL,2,conv);
			t_fin_NR = clock() ;
			t_fin_main = clock();
			double ff= CLOCKS_PER_SEC;

			for ( int l = 0; l < 6; l++ )
			{if ( Y[l]<0.0 ){break;}}
			
			printf ( "\n" );printf ( "Vector Y[i] after:" );
			for ( int k=0 ; k < 6; k++ ){printf ( "%g ",Y[k] );}printf ( "\n" );

			//Convegance FRACTION_SITE dans FRACTION_MOLAIRE
			std::vector <double> X(4);
			double x_c_fcc,x_c_bcc,x_fe_fcc,x_fe_bcc;			
			x_c_bcc  = 3. * Y[1] / (1. + 3. * Y[1]);
			x_fe_bcc = Y[0] / (1. + 3. * Y[1]);
			x_c_fcc  = Y[4] / (1. + Y[4]);			
			x_fe_fcc = Y[3] / (1. + Y[4]);
			X[0] = x_fe_bcc;
			X[1] = x_c_bcc;
			X[2] = x_fe_fcc;
			X[3] = x_c_fcc;
	
			
			if (good)
			{
				so.ecriture ( Y );
				so2.ecritureX ( X );
				//prevVector=Y;
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
			double qhull_cpu = (double) t_fin_qhull - (double) t_debut_qhull;
			double qhull_cpus=((double) t_fin_qhull - (double) t_debut_qhull) / ff;
			double nr_cpu = (double) t_fin_NR - (double) t_debut_NR;
			double nr_cpus= ((double) t_fin_NR - (double) t_debut_NR) / ff;
			double fin_cpu= (double) t_fin_main - (double) t_debut_main;
			double fin_cpus= ((double) t_fin_main - (double) t_debut_main) / ff ;
			//file_v<<N<<" "<<T<<" "<<qhull_cpus<<" "<<nr_cpus<<" "<<fin_cpus<<endl;
			/***************************************Cinetique*************************************************************/
			//if(good)
			//{
				//if((T>=600.) && (T<=2000.))
					
					
					calc2cinetique1(CONCENTRATION_NOMINAL,Y,prevVectorcin,flage);
					
				//calc2cinetique2(CONCENTRATION_NOMINAL,Y,prevVectorcin, flage);
			//}
		}
	}
		counter=time(0)-counter;
		cout<<"time2d:"<<counter<<endl;
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
	file_v.close();

	
	

			
}

bool calc2cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage)
{
	newt nt1;
    nt1=newt();
	
	bool good1;
	bool conv=false;
	
	printf("CALCUL CINETIQUE: \n");				
	nt1=newt();
	int j = 1;
	while ( mat.produit[j]!=NULL )
	{
		nt1.func     =mat.produit[1]->frontiere[1];
		nt1.nrfuncv  =&interface::KTB;
		
		nt1.jacobian =mat.produit[1]->frontiere[1];
		nt1.jacobfunc=&interface::JKB;

		j++;
	}
	if ( nt1.func == NULL ) printf ( "newton.func = null\n" );
	if (nt1.jacobian==NULL) printf ("jacobian.func=null\n");

	int cal,cal1;
	int n = 3000000;
	std::vector <double> Rlong(n), t(n), Xca(n), Xfea(n),\
	Xcb(n), Xfeb(n), XV(n), Gtrf(n),Gfricf(n);
	
	
	float v_ini;
	int diff_step;
	FILE *fpp=fopen("coeff_cin.txt","r") ;
	fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
	fscanf(fpp , "%i  \n" , &diff_step) ;	
	fclose(fpp);

	nt1.fvec.clear();
	nt1.fvec1.clear();
	nt1.nn = 7;
	std::vector <double> Ycinetique ( nt1.nn );
	nt1.fvec.resize ( nt1.nn );
	nt1.fvec1.resize ( nt1.nn ); //Jacobian
	int check=0;
	
	//############# Parameters
	v_ini=0.0;//1.e-;//1.e-2;//100.e-9;
	Rzin=1.e-3;
	double dc;
	dc =  1.5e-5 * exp( -1.421e5 /R / T )*pow(1.e6,2.);
	diff_step=10;
	double raison = dc/diff_step;
	correc_D=raison;
	int compteur = 0;


	double dtempsFeC=0.1;
	double tempsFeC=0.0;
	cal=0;
	cal1=0;	

	double before=0.;
	double after=0.;	
 //while(tempsFeC<=0.3)
 {	

	if (flage)
	{
		//dtempsFeC=0.1;
		Ycinetique[0]=prevVectorcin[0];
		Ycinetique[1]=prevVectorcin[1];
		Ycinetique[2]=prevVectorcin[2];
		Ycinetique[3]=prevVectorcin[3];
		Ycinetique[4]=prevVectorcin[4];
		Ycinetique[5]=prevVectorcin[5];
		Ycinetique[6]=prevVectorcin[6];
		/*before=Rzin;
		Rzin=Rzin+prevVectorcin[6]*dtempsFeC;
		//printf("before=%e \n",before);
		after=Rzin; //printf("after=%e \n",after);
		printf("rez=%e \n",(after-before)/before);
		while((after-before)/before>1e-4)
		{
			dtempsFeC=dtempsFeC/10.;
			Rzin=before+prevVectorcin[6]*dtempsFeC;
			printf("dt=%e R=%e \n",dtempsFeC,Rzin);
			after=Rzin;
			printf("after-before/before=%e \n",(after-before)/after);
		}
		printf("dt= %e \n",dtempsFeC);	
		tempsFeC=tempsFeC+dtempsFeC;*/

	}
	else 
	{
	Ycinetique[0] = Y[0]; //fe
	Ycinetique[1] = Y[1]; //c
	Ycinetique[2] = 1.-Y[1];//va
	Ycinetique[3] = Y[3];//fe
	Ycinetique[4] = Y[4];//c
	Ycinetique[5] = 1.-Y[4];//va
	Ycinetique[6] = v_ini;
	}

	 /*if(cal1>=1)
	 {
		correc_D = dc;
		
		printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
		for ( int k=0 ; k < 7 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );
		
		good1=nt1.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
		
		printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
		for ( int k=0 ; k < 7 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );		
		printf("correct_D=%e \n",correc_D);
		printf("temps=%e  <R>=%e \n",tempsFeC,Rzin);
	 }	
	 else*/
	 { 
 	 while(compteur < diff_step)
	 {	
		printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
		for ( int k=0 ; k < 7 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );
	
		good1=nt1.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
		
		printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
		for ( int k=0 ; k < 7 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );		
		printf("diff_step= %i \n",diff_step);
		printf("compt=%i \n",compteur);
		printf("correct_D=%e \n",correc_D);
		printf("temps=%e  <R>=%e \n",tempsFeC,Rzin);

		
		compteur ++;
		correc_D += raison;
	 }
    }
	


	//************########Register the resulats############*************
	if(conv)
	{
		//Convegance FRACTION_SITE dans FRACTION_MOLAIRE
		std::vector <double> Xcinetique(5);
		double x_c_fcc,x_c_bcc,x_fe_fcc,x_fe_bcc;			
		x_c_bcc  = 3. * Ycinetique[1] / (1. + 3. * Ycinetique[1]);
		x_fe_bcc = Ycinetique[0] / (1. + 3. * Ycinetique[1]);
		x_c_fcc  = Ycinetique[4] / (1. + Ycinetique[4]);			
		x_fe_fcc = Ycinetique[3] / (1. + Ycinetique[4]);
		Xcinetique[0] = x_fe_bcc;
		Xcinetique[1] = x_c_bcc;
		Xcinetique[2] = x_fe_fcc;
		Xcinetique[3] = x_c_fcc;
		Xcinetique[4] = Ycinetique[6];
		
		printf("\n*******Convergence\n");
		so1.ecriture ( Ycinetique );
		so3.ecritureX ( Xcinetique );
		prevVectorcin[0]=Ycinetique[0];
		prevVectorcin[1]=Ycinetique[1];
		prevVectorcin[2]=Ycinetique[2];
		prevVectorcin[3]=Ycinetique[3];
		prevVectorcin[4]=Ycinetique[4];
		prevVectorcin[5]=Ycinetique[5];
		prevVectorcin[6]=Ycinetique[6];
		
		/*Rlong[cal]=Rzin; t[cal]=tempsFeC;
		Xca[cal]=Xcinetique[1];Xfea[cal]=Xcinetique[0];
		Xcb[cal]=Xcinetique[3];Xfeb[cal]=Xcinetique[2];
		XV[cal]=Xcinetique[4];
		cal++;
		printf("tempsFeC=%e dt=%e \n",tempsFeC,dtempsFeC);
		printf("cal=%d \n",cal);*/
		flage=true;
	}
	else{	printf("\n****Pas de convergence*******\n");
			flage=false;}
	
	 cal1++;
}	

	/*t.resize(cal);Rlong.resize(cal);Xca.resize(cal);
	Xfea.resize(cal);Xcb.resize(cal);Xfeb.resize(cal);XV.resize(cal);
	Gtrf.resize(cal);Gfricf.resize(cal);
	ecriture_growth_dissolution_FeC(t,Rlong,Xca,Xfea,Xcb,Xfeb,XV,Gtrf,Gfricf);*/

	
	//###################################
	
	//printf("              Calcul paraequilibre \n");
	//if(Ycinetique[6]>=0.){
		//bool flage_para = false;
		//std::vector<double> prevVectorcin_para(3);
		//calc2cinetique_para(CONCENTRATION_NOMINAL,Ycinetique,prevVectorcin_para,flage_para);}


}

bool calc2cinetique_para(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y, \
						 std::vector<double> &prevVectorcin,bool &flage)
{	

		newt nt2;
		nt2 = newt();
		integer itt2;
		itt2=integer();
		
		
		int f = 1;
		while(mat.produit[f]!=NULL)
		{
			nt2.func	= mat.produit[1]->frontiere[1];
			nt2.nrfuncv	= &interface::KTB_prof;
			
			nt2.jacobian =mat.produit[1]->frontiere[1];
			nt2.jacobfunc=&interface::JKB_prof;
			
			itt2.funcin2d	 =mat.produit[1]->frontiere[1];
			itt2.nrfuncin2d =&interface::funcFeC;
			
			itt2.funcpot2d  =mat.produit[1]->frontiere[1];
			itt2.nrfuncpot2d=&interface::funconeFeC;
			
			itt2.funcderiv = mat.produit[1]->frontiere[1];
			itt2.nrfuncderiv = &interface::funcderivmuFeC;
			
			f++;		
		}
	
	
	 int n=3000000;
	 std::vector <double> Rlong(n), t(n), Xca(n), Xfea(n),\
	 Xcb(n), Xfeb(n), XV(n), Gtrf(n), Gfricf(n);
	 std::vector <double> Gtr(n), GtrV(n);
	 std::vector <double> Cnr(n), Crk(n);
	 std::vector <int> Ccompt(n);
	 std::vector <double> Gnra(n),Grka(n),omega(n),del(n), \
	diff(n),diff_d(n),Gcompt_(n),Ccompt_(n);
	

	
	 double temporaire,dc,raison,tempsFeC,dtempsFeC,hdiff,sdiff,param_cin,OMEGA,SURSAT;
	 double x_c_fcc,x_c_bcc,x_fe_fcc,x_fe_bcc;
	 double U_c_fcc,U_c_bcc,U_fe_fcc,U_fe_bcc;
	 int cal, cal1, tr, compteur,tr1,opt;
	 double cal2;
	 double Vopt[100000],tdiff[10];

	 float v_ini;
	 int diff_step;
	 FILE *fpp=fopen("coeff_cin.txt","r") ;
	 fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
	 fscanf(fpp , "%i  \n" , &diff_step) ;	
	 fclose(fpp);
	 bool conv_para = false;

	
	//######  The parameters for Newton-Raphson ########
	int check=0;
	nt2.fvec.clear();
	nt2.fvec1.clear();
	nt2.nn = 3;
	nt2.fvec.resize ( nt2.nn );
	nt2.fvec1.resize ( nt2.nn );
	std::vector <double> Ycinetique_para ( nt2.nn );

	
	//############# Parameters
	Rzin=1.e-3;//mean size of precipitation
	taille_gradient=1.e-3; // the size of interface
	Vm=1e13;
	M=0.035*exp(-17700/T)*pow(1e6,4);
	std::vector <double> Ycinetique_para_profile ( 7 );
	dtempsFeC=0.1;
	tempsFeC=0.0;
	cal=0;
	cal1=0;
	

	dc =  1.5e-5 * exp( -1.421e5 /R / T )*pow(1e6,2.);
	diff_step=10;
	raison = dc/diff_step;
	correc_D=raison;
	compteur=0;

	
	//DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER 1 FOR PROFILE PARAEQUILIBRE\n");
	//(Prediction-correction method)\n");
	//of type y' = f(x,y)\n\n");	
	int fi,m,s;
	double xi,xf,yi;
	xi=0;//begin value x
	xf=1;//0.00001;//end value x
	m=10;//number of points
	fi=10;//finesse
	s=2;
	double tx[100],ty[100];
	std::vector < std::vector<double> > P(m+1);
	for(int i = 0; i < P.size(); i++)
	{
		P[i].resize(s);
	}
	std::vector < std::vector<double> > P1(m+1);
	for(int i = 0; i < P1.size(); i++)
	{
		P1[i].resize(3);
	}
	
	
	cal=0;
	cal1=0;
	double before=0.;
	double after=0.;
	//Recursion of time tnew=told+delta(t)
	//while(tempsFeC<=0.3)
	{	
		if (flage)
		{
			dtempsFeC=0.1;
			Ycinetique_para[0] = prevVectorcin[0]; //c_alpha
			Ycinetique_para[1] = prevVectorcin[1];
			Ycinetique_para[2] = prevVectorcin[2]; //c
			VIT_TEMPORAIRE = Ycinetique_para[2];
			before=Rzin;
			Rzin=Rzin+prevVectorcin[2]*dtempsFeC;
			//printf("before=%e \n",before);
			after=Rzin; //printf("after=%e \n",after);
			printf("rez=%e \n",(after-before)/before);
			while((after-before)/before>1e-4)
			{
				dtempsFeC=dtempsFeC/10.;
				Rzin=before+prevVectorcin[2]*dtempsFeC;
				printf("dt=%e R=%e \n",dtempsFeC,Rzin);
				after=Rzin;
				printf("after-before/before=%e \n",(after-before)/after);
			}
			printf("dt= %e \n",dtempsFeC);	
			tempsFeC=tempsFeC+dtempsFeC;
		}	
		else		
		{
			Ycinetique_para[0] = Y[1]; //c_alpha
			Ycinetique_para[1] = Y[4];
			Ycinetique_para[2] = Y[6]; //c
			VIT_TEMPORAIRE = Ycinetique_para[2];
		}
	
		tr=0;

	
			correc_D=dc;
			taille_gradient=1e-3;		
			printf("\n\n");
			x_c_bcc=3*Ycinetique_para[0]/(1.+3*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[1]/(1.+Ycinetique_para[1]);
			printf ( "\n" );printf ( "Vector Ycinetique_para[i] before NR:" );
			for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );
			printf("Before NR : X(C,BCC) = %e	X(C,FCC) = %e  <R>=%e	taille=%e \n",x_c_bcc, \
				   x_c_fcc,Rzin,taille_gradient);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);
			U_c_fcc=x_c_fcc/(1-x_c_fcc);	
			printf("After NR: U(C,BCC) = %e	U(C,FCC) = %e\n",U_c_bcc,U_c_fcc);

	
			flage=nt2.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
	
			x_c_bcc=3*Ycinetique_para[0]/(1.+3*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[1]/(1.+Ycinetique_para[1]);
			printf ( "\n" );printf ( "Vector Ycinetique_para[i] after:" );
			for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );
			printf("\n\n"); 	
			printf("After NR: X(C,BCC) = %e	X(C,FCC) = %e  <R>=%e   taille=%e \n",x_c_bcc, \
				   x_c_fcc,Rzin,taille_gradient);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);
			U_c_fcc=x_c_fcc/(1-x_c_fcc);													 
		    printf("After NR: U(C,BCC) = %e	U(C,FCC) = %e\n",U_c_bcc,U_c_fcc);
			
			printf("diff_step= %i \n",diff_step);
			printf("compt=%i \n",compteur);
			printf("correct_D=%e \n",correc_D);
			printf("temps=%e  <R>=%e \n",tempsFeC,Rzin);


	//If Newton-Raphson gives the solution we calculate the profile of concentration 
	if(conv_para)
	{
		printf("\n DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER 1 FOR PROFILE PARAEQUILIBRE\n");
		printf("         (Prediction-correction method)\n");
		printf("              of type y' = f(x,y)\n\n");
		VIT_TEMPORAIRE=Ycinetique_para[2];
		double x_ci=3. * Ycinetique_para[0] / (1. + 3. * Ycinetique_para[0]);

		yi=x_ci/(1.-x_ci);
		x_c0=yi;
		correc_D=dc;		
		printf ( "\n" );printf ( "Start Runga-Kutta:" );
		printf("X(C,BCC) = %e  V = %lf correct_D=%e taille=%e \n",x_ci,VIT_TEMPORAIRE,correc_D, \
					   taille_gradient);
		printf("U(C,BCC) = %e \n",yi);
		compteur_c = itt2.equadiff_pc(tx,ty,xi,xf,yi,m,s,fi,CONCENTRATION_NOMINAL,Ycinetique_para,P,P1);
				//calcul G diff
		hdiff =(xf-xi)/m;
		sdiff = P1[0][0]+P1[m][0];
		tdiff[10];
		for(int k=0; k <=m; k++)
		{
			sdiff=sdiff+2*P1[k][0];
			tdiff[k]=sdiff*hdiff/2; 
		}
		diffFeC = (1./VIT_TEMPORAIRE)*(tdiff[10]-tdiff[0])/(R*T);	
		fricFeC = VIT_TEMPORAIRE*(Vm/M)/(R*T);
		printf ( "\n" );printf ("End Runga-Kutta\n");
		printf("X(C,BCC) = %e  X(C,FCC)= %e V = %lf mu_c= %e Gtr=%e	 Gfric=%e \n\n",ty[0]/(1.+ty[0]),ty[m]/(1.+ty[m]), \
		VIT_TEMPORAIRE,compteur_c/(R*T),diffFeC,fricFeC);
		printf("U(C,BCC) = %e  U(C,FCC)= %e \n\n",ty[0],ty[m]);

				
		printf("\n*******Convergence growth diffusion\n");
		//***********Register the concentration of txo side of interface in site fraction and mole fraction
		std::vector <double> Xcinetique(5);
		x_c_bcc  = 3. * Ycinetique_para[0] / (1. + 3. * Ycinetique_para[0]);
		x_fe_bcc = 1.-x_c_bcc;
		x_c_fcc  = Ycinetique_para[1] / (1. + Ycinetique_para[1]);			
		x_fe_fcc = 1.-x_c_fcc;
		Xcinetique[0] = x_fe_bcc;
		Xcinetique[1] = x_c_bcc;
		Xcinetique[2] = x_fe_fcc;
		Xcinetique[3] = x_c_fcc;
		Xcinetique[4] = Ycinetique_para[2];
		so5.ecritureX ( Xcinetique );	
		
		std::vector <double> Ucinetique(5);
		U_c_bcc  = x_c_bcc/(1-x_c_bcc);
		U_fe_bcc = x_fe_bcc/(1-x_c_bcc);
		U_c_fcc  = x_c_fcc/(1-x_c_fcc);			
		U_fe_fcc = x_fe_fcc/(1-x_c_fcc);
		Ucinetique[0] = U_c_bcc;
		Ucinetique[1] = U_fe_bcc;
		Ucinetique[2] = U_c_fcc;
		Ucinetique[3] = U_fe_fcc;
		Ucinetique[4] = Ycinetique_para[2];
		so4.ecritureX ( Ucinetique );	
	
		
		//************########Register the resulats############*************
		OMEGA=(U_c_fcc-(CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0])))/(U_c_fcc-U_c_bcc);
		SURSAT=2.*Rzin*(1.-OMEGA)/OMEGA;
		//printf("position=%e \n",SURSAT);
		Gtr[tr]=gtr;
		GtrV[tr]=Ycinetique_para[2];	
		Cnr[tr]=U_c_bcc;
		Gnra[tr]=U_c_fcc;
		omega[tr]=muc;
		diff[tr]=Ycinetique_para[2]*(Vm/M)/(R*T);
		diff_d[tr]=taille_gradient;
		Ccompt[tr]=tr;	
		Crk[tr]=3.*Y[1]/(1.+3.*Y[1])/(1-3.*Y[1]/(1.+3.*Y[1]));
		Grka[tr]=Y[4]/(1.+Y[4])/(1-Y[4]/(1.+Y[4]));
		del[tr]=Y[6];
		Ccompt_[tr]=CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]);
		tr++;
		Cnr.resize(tr);Gnra.resize(tr);GtrV.resize(tr);Ccompt.resize(tr);
		diff_d.resize(tr);diff.resize(tr);omega.resize(tr);Gtr.resize(tr);
		Crk.resize(tr);Grka.resize(tr);del.resize(tr);Ccompt_.resize(tr);
		
		if((cal1==0)||(cal1==11000)||(cal1==15000)||(cal1==20000)||(cal1==35000)||(cal1==60000)||(cal1==110000)||(cal1==150000)||(cal1==210000))
		{	
			for(int k=0;k<=m;k++){ty[k]=ty[k]/(ty[k]+1);}	
		ecriture2D(tempsFeC,diff_d,diff,omega,Cnr,Gnra,GtrV,Gtr,Ccompt,Crk,Grka,del,Ccompt_);
		//***********All the profile of concentration of carbon as function of space (mikrometr)
		//**********for all temps and mean size of precipitation	
		ecriture_profileC_FeC(tempsFeC,Rzin,SURSAT,m,tx,ty);
		//**********All the profile of chemical potentiel of carbon and iron
		ecriture_profilePOT_FeC(tempsFeC,Rzin,m,tx,P);
		}
		//***********Recursion for tnew=told+delta(t)
		prevVectorcin[0]=Ycinetique_para[0];
		prevVectorcin[1]=Ycinetique_para[1];
		prevVectorcin[2]=Ycinetique_para[2];
		
		printf("temps = %lf R = %e  cal=%i \n",tempsFeC,Rzin,cal);
		Rlong[cal]=Rzin; t[cal]=tempsFeC;
		Xca[cal]=Ucinetique[1];Xfea[cal]=Ucinetique[0];
		Xcb[cal]=Ucinetique[3];Xfeb[cal]=Ucinetique[2];
		XV[cal]=Ucinetique[4];Gtrf[cal]=gtr;Gfricf[cal]=muc;
		cal++;
		flage=true;
	}
		
   else
	{	
		printf("No convergence growth diffusion and the second try");
		if (flage)
		{
			dtempsFeC=0.1;
			Ycinetique_para[0] = prevVectorcin[0]; //c_alpha
			Ycinetique_para[1] = prevVectorcin[1];
			Ycinetique_para[2] = prevVectorcin[2]; //c
			VIT_TEMPORAIRE = Ycinetique_para[2];
			 before=Rzin;
			 Rzin=Rzin+prevVectorcin[2]*dtempsFeC;
			 //printf("before=%e \n",before);
			 after=Rzin; //printf("after=%e \n",after);
			 printf("rez=%e \n",(after-before)/before);
			 while((after-before)/before>1e-4)
			 {
			 dtempsFeC=dtempsFeC/10.;
			 Rzin=before+prevVectorcin[2]*dtempsFeC;
			 printf("dt=%e R=%e \n",dtempsFeC,Rzin);
			 after=Rzin;
			 printf("after-before/before=%e \n",(after-before)/after);
			 }
			 printf("dt= %e \n",dtempsFeC);	
			 tempsFeC=tempsFeC+dtempsFeC;

			
		}	
		else
		{
			Ycinetique_para[0] = Y[1]; //c_alpha
			Ycinetique_para[1] = Y[4];
			Ycinetique_para[2] = Y[6]; //c
			VIT_TEMPORAIRE = Ycinetique_para[2];
		}
		
		
		double dtaille;
		double dtaille_tem;
		taille_gradient=1.e-6;
		tr=0;
		conv_para=true;
		while((taille_gradient<=1.e-1) && (conv_para))
		{		
			correc_D=dc;		
			if(taille_gradient<=1.e-4){dtaille=1.e-1/100000.;}
			else{dtaille=1.e-1/1000000.;}
			
			x_c_bcc=3*Ycinetique_para[0]/(1.+3*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[1]/(1.+Ycinetique_para[1]);
			printf("\n\n");
			printf ( "\n" );printf ( "Vector Ycinetique_para[i] before NR:" );
			for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );
			printf("Before NR : X(C,BCC) = %e	X(C,FCC) = %e  <R>=%e	taille=%e \n",x_c_bcc, \
				   x_c_fcc,Rzin,taille_gradient);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);
			U_c_fcc=x_c_fcc/(1-x_c_fcc);			
			printf("Before NR : U(C,BCC) = %e	U(C,FCC) = %e \n",U_c_bcc,U_c_fcc);
			
			flage=nt2.calc ( Ycinetique_para, &check,CONCENTRATION_NOMINAL,1,conv_para);
			
			x_c_bcc=3*Ycinetique_para[0]/(1.+3*Ycinetique_para[0]);
			x_c_fcc=Ycinetique_para[1]/(1.+Ycinetique_para[1]);
			printf ( "\n" );printf ( "Vector Ycinetique_para[i] after:" );
			for ( int k=0 ; k < 3 ; k++ ){printf ( "%g ",Ycinetique_para[k] );}printf ( "\n" );
			printf("\n\n"); 	
			printf("After NR: X(C,BCC) = %e	X(C,FCC) = %e  <R>=%e   taille=%e \n",x_c_bcc, \
				   x_c_fcc,Rzin,taille_gradient);
			U_c_bcc=x_c_bcc/(1-x_c_bcc);
			U_c_fcc=x_c_fcc/(1-x_c_fcc);			
			printf("Before NR : U(C,BCC) = %e	U(C,FCC) = %e \n",U_c_bcc,U_c_fcc);
			
			if(conv_para)
			{
				printf("convergence\n");		
				x_c_bcc  = 3*Ycinetique_para[0]/(1.+3.*Ycinetique_para[0]);
				x_c_fcc  = Ycinetique_para[1]/(1.+Ycinetique_para[1]);	
				U_c_bcc=x_c_bcc/(1-x_c_bcc);
				U_c_fcc=x_c_fcc/(1-x_c_fcc);
				OMEGA=(U_c_fcc-(CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0])))/(U_c_fcc-U_c_bcc);
				SURSAT=2.*Rzin*(1.-OMEGA)/OMEGA;	 
				//printf("position=%e \n",SURSAT);	 
				Gtr[tr]=gtr;
				GtrV[tr]=Ycinetique_para[2];	
				Cnr[tr]=U_c_bcc;
				Gnra[tr]=U_c_fcc;
				omega[tr]=muc;
				diff[tr]=Ycinetique_para[2]*(Vm/M)/(R*T);
				diff_d[tr]=taille_gradient;
				Ccompt[tr]=tr;	
				Crk[tr]=3.*Y[1]/(1.+3.*Y[1])/(1-3.*Y[1]/(1.+3.*Y[1]));
				Grka[tr]=Y[4]/(1.+Y[4])/(1-Y[4]/(1.+Y[4]));
				del[tr]=Y[6];
				Ccompt_[tr]=CONCENTRATION_NOMINAL[0]/(1-CONCENTRATION_NOMINAL[0]);
				tr++;
				
				prevVectorcin[0]=Ycinetique_para[0];
				prevVectorcin[1]=Ycinetique_para[1];
				prevVectorcin[2]=Ycinetique_para[2];
				dtaille_tem=taille_gradient;
				flage=true;	 
			}	 
			else{printf("The limite of size of interface\n");conv_para=false;flage=false;break;}			 
			taille_gradient+=dtaille;
		}
		Cnr.resize(tr);Gnra.resize(tr);GtrV.resize(tr);Ccompt.resize(tr);
		diff_d.resize(tr);diff.resize(tr);omega.resize(tr);Gtr.resize(tr);
		Crk.resize(tr);Grka.resize(tr);del.resize(tr);Ccompt_.resize(tr);
		
		
		printf("\n DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER 1 FOR PROFILE PARAEQUILIBRE\n");
		printf("         (Prediction-correction method)\n");
		printf("              of type y' = f(x,y)\n\n");
		double x_ci=3. * prevVectorcin[0] / (1. + 3. * prevVectorcin[0]);
		yi=x_ci/(1.-x_ci);
		x_c0=yi;
		taille_gradient=dtaille_tem;
		VIT_TEMPORAIRE=prevVectorcin[2];
		
		printf ( "\n" );printf ( "Start Runga-Kutta:" );
		printf("X(C,BCC) = %e  V = %lf correct_D=%e taille=%e \n",x_ci,VIT_TEMPORAIRE,correc_D, \
			   taille_gradient);
		printf("U(C,BCC) = %e \n",yi);
		compteur_c = itt2.equadiff_pc(tx,ty,xi,xf,yi,m,s,fi,CONCENTRATION_NOMINAL,prevVectorcin,P,P1);
		//calcul G diff
		hdiff =(xf-xi)/m;
		sdiff = P1[0][0]+P1[m][0];
		tdiff[10];
		for(int k=0; k <=m; k++)
		{
			sdiff=sdiff+2*P1[k][0];
			tdiff[k]=sdiff*hdiff/2; 
		}
		diffFeC = (1./VIT_TEMPORAIRE)*(tdiff[10]-tdiff[0])/(R*T);	
		fricFeC = VIT_TEMPORAIRE*(Vm/M)/(R*T);
		printf ( "\n" );printf ("End Runga-Kutta\n");
		printf("X(C,BCC) = %e  X(C,FCC)= %e V = %lf mu_c= %e Gtr=%e	 Gfric=%e \n\n",ty[0]/(1.+ty[0]),ty[m]/(1.+ty[m]), \
			   VIT_TEMPORAIRE,compteur_c/(R*T),diffFeC,fricFeC,fricFeC);
		printf("U(C,BCC) = %e  U(C,FCC)= %e\n\n",ty[0],ty[m]);
		
		if((cal1==0)||(cal1==11000)||(cal1==15000)||(cal1==20000)||(cal1==35000)||(cal1==60000)||(cal1==110000)||(cal1==150000)||(cal1==209999))
		{
			//for(int k=0;k<=m;k++){ty[k]=ty[k]/(ty[k]+1);}
			ecriture2D(tempsFeC,diff_d,diff,omega,Cnr,Gnra,GtrV,Gtr,Ccompt,Crk,Grka,del,Ccompt_);
			ecriture_profileC_FeC(tempsFeC,Rzin,SURSAT,m,tx,ty);
			//**********All the profile of chemical potentiel of carbon and iron
			ecriture_profilePOT_FeC(tempsFeC,Rzin,m,tx,P);
		}	
		
		 U_c_bcc  = 3. * prevVectorcin[0] / (1. + 3. * prevVectorcin[0])/(1-3. * prevVectorcin[0]/(1. + 3. * prevVectorcin[0]));
		 U_fe_bcc = 1;
		 U_c_fcc  = prevVectorcin[1] / (1. + prevVectorcin[1])/(1-prevVectorcin[1]/(1. + prevVectorcin[1]));			
		 U_fe_fcc = 1;
		 printf("temps = %lf R = %e  cal=%i \n",tempsFeC,Rzin,cal);
		 Rlong[cal]=Rzin; t[cal]=tempsFeC;
		 Xca[cal]=U_c_bcc;Xfea[cal]=U_fe_bcc;
		 Xcb[cal]=U_c_fcc;Xfeb[cal]=U_fe_fcc;
		 XV[cal]=prevVectorcin[2];Gtrf[cal]=gtr;Gfricf[cal]=muc;
		 cal++;
	}	
		
		
		printf("cal1=%d \n",cal1);
		cal1++;
		
	}
	
	 t.resize(cal);Rlong.resize(cal);Xca.resize(cal);
	 Xfea.resize(cal);Xcb.resize(cal);Xfeb.resize(cal);XV.resize(cal);
	 Gtrf.resize(cal);Gfricf.resize(cal);
	 ecriture_growth_dissolution_FeC(t,Rlong,Xca,Xfea,Xcb,Xfeb,XV,Gtrf,Gfricf);
		
}

double f_vitess_x1(double V, double x_c_gamma,std::vector <double> CONCENTRATION_NOMINAL)
{	
	double discriminant=pow((CONCENTRATION_NOMINAL[0]+x_c_gamma),2)- \
	4*x_c_gamma*CONCENTRATION_NOMINAL[0]+2.*(correc_D/(Rzin*V))* \
	pow((x_c_gamma-CONCENTRATION_NOMINAL[0]),2);
	
	double sqrt_discriminant=sqrt(discriminant);
	
	return (CONCENTRATION_NOMINAL[0]+x_c_gamma-sqrt_discriminant)/2.;
}

double f_vitess_x2(double V, double x_c_gamma,std::vector <double> CONCENTRATION_NOMINAL)
{
	double discriminant=pow((CONCENTRATION_NOMINAL[0]+x_c_gamma),2)- \
	4*x_c_gamma*CONCENTRATION_NOMINAL[0]+2.*(correc_D/(Rzin*V))* \
	pow((x_c_gamma-CONCENTRATION_NOMINAL[0]),2);
	
	double sqrt_discriminant=sqrt(discriminant);
	
	return (CONCENTRATION_NOMINAL[0]+x_c_gamma+sqrt_discriminant)/2.;
}

double regulafalsi(double a,double b,double V,std::vector <double> CONCENTRATION_NOMINAL)
{
	/*
	int n;
	double c;
	
	for(n=0;n<50;n++)
	{
		c=(a*f_vitess(b,x_c_gamma,CONCENTRATION_NOMINAL))/( \
		f_vitess(b,x_c_gamma,CONCENTRATION_NOMINAL)-f_vitess(a,x_c_gamma,CONCENTRATION_NOMINAL));
		printf("c=%e a=%e b=%e\n",c,a,b);
		if(fabs(f_vitess(c,x_c_gamma,CONCENTRATION_NOMINAL))<1e-10)break;
		if(f_vitess(a,x_c_gamma,CONCENTRATION_NOMINAL)*f_vitess(c,x_c_gamma,CONCENTRATION_NOMINAL)<0.0)b=c;
		else a=c;
	}
	return c;*/
	int j;
	float fl,f,dx,swap,xl,rts;
	
	fl=f_vitess_x1(V,a,CONCENTRATION_NOMINAL);
	f=f_vitess_x1(V,b,CONCENTRATION_NOMINAL);
	if (fabs(fl) < fabs(f)) {
		rts=a;
		xl=b;
		swap=fl;
		fl=f;
		f=swap;
	} else {
		xl=a;
		rts=b;
	}
	for (j=1;j<=100;j++) {
		dx=(xl-rts)*f/(f-fl);
		xl=rts;
		fl=f;
		rts += dx;
		f=f_vitess_x1(V,rts,CONCENTRATION_NOMINAL);
		if (fabs(dx) < 1.e-7 || f == 0.0) return rts;
	}
	printf("Maximum number of iterations exceeded in rtsec\n");
	return 0.0;	
}


bool calc2cinetique2(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage)
{
	newt nt1;

	
	printf("CALCUL CINETIQUE: \n");				
	nt1=newt();
	int j = 1;
	while ( mat.produit[j]!=NULL )
	{
		nt1.func     =mat.produit[1]->frontiere[1];
		nt1.nrfuncv  =&interface::KTB1;
		
		nt1.jacobian =mat.produit[1]->frontiere[1];
		nt1.jacobfunc=&interface::JKB1;
		j++;
	}
	if ( nt1.func == NULL ) printf ( "newton.func = null\n" );
	if (nt1.jacobian==NULL) printf ("jacobian.func=null\n");
	
	float v_ini;
	int diff_step;
	bool conv = false;
	
	FILE *fpp=fopen("coeff_cin.txt","r") ;
	fscanf(fpp , "%e \n %e \n" ,&CIN,&v_ini) ;
	fscanf(fpp , "%i  \n" , &diff_step) ;	
	fclose(fpp);
	
	nt1.fvec.clear();
	nt1.fvec1.clear();
	nt1.nn = 7;
	std::vector <double> Ycinetique ( nt1.nn );
	nt1.fvec.resize ( nt1.nn );
	nt1.fvec1.resize ( nt1.nn ); //Jacobian
	int check=0;
	if (flage)
	{
		Ycinetique=prevVectorcin;
		//taille_gradient_alpha=taille_gradient_alpha+1.e-9*prevVectorcin[6]*dtemps;
		//taille_gradient_betta=taille_gradient_betta+1.e-9*prevVectorcin[6]*dtemps;
	}
	else
	{	
		Ycinetique[0] = Y[0]; //fe
		Ycinetique[1] = Y[1]; //c
		Ycinetique[2] = 1.-Y[1];//va
		Ycinetique[3] = Y[3];//fe
		Ycinetique[4] = Y[4];//c
		Ycinetique[5] = 1.-Y[4];//va
		Ycinetique[6] = v_ini;
	}
	
	
	printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
	for ( int k=0 ; k < 7 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );
	
	nt1.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
	
	printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
	for ( int k=0 ; k < 7 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );
	
	if(conv)
	{
		so1.ecriture ( Ycinetique );
		prevVectorcin=Ycinetique;
		flage=true;
		return true;
	}
	else{flage=false; return false;}
}	


	
