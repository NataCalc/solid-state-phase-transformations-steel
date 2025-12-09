/*
 *  4dpreci.cpp
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

void calc4dequilibre();
bool calc4cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);

extern clock_t t_debut_main , t_debut_qhull , t_fin_qhull , t_debut_NR , t_fin_NR , t_fin_main ;


void calc4dequilibre()
{	
	t_debut_main = clock();
	int dim=4;
	def=4;
	
	newt nt;
    qhull qh(dim);

	double counter=time(0);
	double prevT=0;
	bool conv = false;
	bool flage = false;
	
	vector<double> prevVector(11);
	for(int i=0;i<11;i++)
	{prevVector[i]=0;}	
	std::vector<double> prevVectorcin(11);
	for(int i=0;i<11;i++){prevVectorcin[i]=0;}
	
	std::vector <double> CONCENTRATION_NOMINAL ( 4 );
	//X(C)=8.6991208E-2, X(CR)=2.5118552E-2, X(FE)=0.88026677, X(MO)=7.6234728E-3 (2.E-2)
	CONCENTRATION_NOMINAL[0] = 9.268719e-3;              		                     				              // carbon -C
	CONCENTRATION_NOMINAL[2] = 2.6763257e-2;  		                                     			              // chrome -CR
	CONCENTRATION_NOMINAL[3] = 8.1226401e-3;									                                      // Mo
	CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0]-CONCENTRATION_NOMINAL[2]-CONCENTRATION_NOMINAL[3];// iron -FE
	printf ("Begin calcul system FE-2.e-5CR-2.e-3C-1.4e-2MO equilibre thermodynamique (4D): \n");
	
	fstream file_v;
	char filename[30];
	sprintf ( filename,"CPU4d.txt");
	file_v.open ( filename,fstream::out );
	if ( file_v.fail() )
		cout<<"Error opening file"<<endl;
	
	int N ;
	
	//for(N=3;N<10;N++)
	{		
	for (int iterseg=1;iterseg<=nbseg;iterseg++)
	{
		temps=seg[iterseg].date;
		dt=(seg[iterseg+1].date-temps)/seg[iterseg+1].pas;
		
		T=seg[iterseg].temperature;
		dT= dt * (seg[iterseg+1].temperature-seg[iterseg].temperature)
		/ (seg[iterseg+1].date-seg[iterseg].date);
		
		while (temps<seg[iterseg+1].date)
		{
			bool good;
			bool prevExts=false;
			qh=qhull(dim);
			temps+=dt;
			T+=dT;
			
			if ( ( T>=1130.0 ) && ( T<=1730.0 ) ) continue;
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
				nt.func=mat.produit[1]->frontiere[1];
				nt.nrfuncv= &interface::SQ4;
				
				nt.jacobian=mat.produit[1]->frontiere[1];
				nt.jacobfunc = &interface::JQ4;
				
				qh.func1    =mat.produit[1]->frontiere[1];
				qh.nrfuncv1 =&interface::phase1_quaternaire;
				
				qh.func2    =mat.produit[1]->frontiere[1];
				qh.nrfuncv2 =&interface::phase2_quaternaire;
				
				qh.Jacob1=mat.produit[1]->frontiere[1];
				qh.nrJacob1=&interface::JacobianQhull4Dalfa;
				
				qh.Jacob2=mat.produit[1]->frontiere[1];
				qh.nrJacob2=&interface::JacobianQhull4Dbetta;
				i++;
			}
			if ( nt.func == NULL ) printf ( "newton.func = null\n" );
			if ( qh.func1== NULL ) printf ( "qhull.func1 = null\n" );
			
			mat.Thermo->read ( T );
			cem.classe0->Thermo->read ( T );
			
			nt.fvec.clear();
			nt.fvec1.clear();
			nt.nn = 11;
			std::vector <double> Y ( nt.nn );
			nt.fvec.resize ( nt.nn );
			nt.fvec1.resize ( nt.nn ); //Jacobian
			int check=0;
			
			printf("T= %lf \n",T);
			/*if ((prevExts) && (good))
			{
				Y=prevVector;
			}
			else*/
			{
				t_debut_qhull=clock();
				N=35;
                qh.discretization (T,N);
                qh.build_hull();
                vector<point> ext=qh.extremums();
				t_fin_qhull = clock();
                printf ( "T = %lf \n", T );
                //cout<<"Extremums:"<<ext.size()<<endl;
                //cout<<qh<<endl;
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

				
				
				Y[1] = ext[ ( ext[0].func()==0 ) ?0:1][1];//cr
				if ( Y[1]<=0 ) Y[1]=1.e-4;
				
				Y[2] = ext[ ( ext[0].func()==0 ) ?0:1][2];//mo
				if ( Y[2]<=0 ) Y[2]=1.e-4;
				
				Y[0] = 1-Y[1]-Y[2];//fe
				if ( Y[0]<=0 ) Y[0]=1.e-4;
				
				Y[3] = ext[ ( ext[0].func()==0 ) ?0:1][0];//c
				if ( Y[3]<=0 ) Y[3]=1.e-4;
				
				Y[4] = 1-Y[3]; //va
				if ( Y[4]<=0 ) Y[4]=1.e-4;
				
				Y[6] = ext[ ( ext[0].func()==0 ) ?1:0][1];//cr
				if ( Y[6]<=0 ) Y[6]=1.e-4;
				
				Y[7] = ext[ ( ext[0].func()==0 ) ?1:0][2];//mo
				if ( Y[7]<=0 ) Y[7]=1.e-4;
				
				Y[5] = 1-Y[6]-Y[7]; //fe
				if ( Y[5]<=0 ) Y[5]=1.e-4;
				
				Y[8] = ext[ ( ext[0].func()==0 ) ?1:0][0];//c
				if ( Y[8]<=0 ) Y[8]=1.e-4;
				
				Y[9] = 1-Y[8]; //va
				if ( Y[9]<=0 ) Y[9]=1.e-4;
				
				Y[10]= 0.45;
				
			}
			
			
			printf ( "\n" );printf ( "Vector Y[i] before:" );
			for ( int k=0 ; k < 11 ; k++ ){printf ( "%g ",Y[k] );}printf ( "\n" );

			t_debut_NR = clock();			
			//good=nt.calc ( Y, &check,CONCENTRATION_NOMINAL,4,conv);
			t_fin_NR = clock() ;
			t_fin_main = clock();
			double ff= CLOCKS_PER_SEC;
			
			printf ( "\n" );printf ( "Vector Y[i] after:" );
			for ( int k=0 ; k < 11 ; k++ ){printf ( "%g ",Y[k] );}printf ( "\n" );
			
			//convergance FRACTION_SITE dans FRACTION_MOLAIRE
			std::vector <double> X(9);
			double x_c_fcc,x_cr_fcc,x_mo_fcc,x_c_bcc,x_cr_bcc,x_mo_bcc;			
			x_cr_bcc = Y[1] / (1. + 3. * Y[3]) ;
			x_mo_bcc = Y[2] / (1. + 3. * Y[3]) ;
			x_c_bcc  = 3. * Y[3] / (1. + 3. * Y[3]) ;
			x_cr_fcc = Y[6] / (1. + Y[8]) ;
			x_mo_fcc = Y[7] / (1. + Y[8]) ;
			x_c_fcc  = Y[8] / (1. + Y[8]) ;			
			
			X[1] = x_cr_bcc;//cr			
			X[2] = x_mo_bcc;//mo
			X[3] = x_c_bcc;//c
			X[0] = 1-X[1]-X[2]-X[3];//fe			
			X[5] = x_cr_fcc;//cr			
			X[6] = x_mo_fcc;//mo
			X[7] = x_c_fcc;//c
			X[4] = 1-X[5]-X[6]-X[7]; //fe
			X[8] = Y[10];
			
			
			
			//if (good)
			{
				so.ecriture ( Y );
				so2.ecritureX ( X );
				prevVector=Y;
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


			
			//if((prevExts) && (good))
			 //{
			 //if(calc4cinetique1(CONCENTRATION_NOMINAL,Y,prevVectorcin,flage))
			 //so1.ecriture ( Y );
			 //}
			
						
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
file_v.close();

    counter=time(0)-counter;
    cout<<"time4d:"<<counter<<endl;
}


bool calc4cinetique1(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage)
{
	newt nt1;
	bool good1;
	bool conv = false;
	
	nt1=newt();	
	int j = 1;
	while ( mat.produit[j]!=NULL )
	{
		nt1.func=mat.produit[1]->frontiere[1];
		nt1.nrfuncv= &interface::KT4;
		
		nt1.jacobian=mat.produit[1]->frontiere[1];
		nt1.jacobfunc = &interface::JK4;
		j++;
	}
	if ( nt1.func == NULL ) printf ( "newton.func = null\n" );
	
	int check=0;
	nt1.fvec.clear();
	nt1.fvec1.clear();	
	nt1.nn = 11;
	std::vector <double> Ycinetique ( nt1.nn );
	nt1.fvec.resize ( nt1.nn );
	nt1.fvec1.resize ( nt1.nn ); //Jacobian
	
	if (flage)
	{
		if(prevVectorcin[10]<=0.0){prevVectorcin[10]=100;}
		Ycinetique=prevVectorcin;
		//taille_gradient=taille_gradient+prevVectorcin[10]*dtemps;
	}
	else
	{
	Ycinetique[1] = Y[1]; //cr
	Ycinetique[2] = Y[2]; //mo
	Ycinetique[0] = Y[0]; //fe
	Ycinetique[3] = Y[3]; //c 
	Ycinetique[4] = Y[4]; //va
	Ycinetique[6] = Y[6]; //cr
	Ycinetique[7] = Y[7]; //mo
	Ycinetique[5] = Y[5]; //fe
	Ycinetique[8] = Y[8]; //c				
	Ycinetique[9] = Y[9]; //va
	Ycinetique[10]= 100; //vitess
	}			
	
	printf ( "\n" );printf ( "Vector Ycinetique[i] before:" );
	for ( int k=0 ; k < 11 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );
	
	good1=nt1.calc ( Ycinetique, &check,CONCENTRATION_NOMINAL,1,conv);
	
	printf ( "\n" );printf ( "Vector Ycinetique[i] after:" );
	for ( int k=0 ; k < 11 ; k++ ){printf ( "%g ",Ycinetique[k] );}printf ( "\n" );
	
	if (good1)
	{		
		//convergance FRACTION_SITE dans FRACTION_MOLAIRE
		std::vector <double> Xcinetique(9);
		double x_c_fcc,x_cr_fcc,x_mo_fcc,x_c_bcc,x_cr_bcc,x_mo_bcc;			
		x_cr_bcc = Ycinetique[1] / (1. + 3. * Ycinetique[3]) ;
		x_mo_bcc = Ycinetique[2] / (1. + 3. * Ycinetique[3]) ;
		x_c_bcc  = 3. * Ycinetique[3] / (1. + 3. * Ycinetique[3]) ;
		x_cr_fcc = Ycinetique[6] / (1. + Ycinetique[8]) ;
		x_mo_fcc = Ycinetique[7] / (1. + Ycinetique[8]) ;
		x_c_fcc  = Ycinetique[8] / (1. + Ycinetique[8]) ;			
		
		Xcinetique[1] = x_cr_bcc;//cr			
		Xcinetique[2] = x_mo_bcc;//mo
		Xcinetique[3] = x_c_bcc;//c
		Xcinetique[0] = 1-Xcinetique[1]-Xcinetique[2]-Xcinetique[3];//fe			
		Xcinetique[5] = x_cr_fcc;//cr			
		Xcinetique[6] = x_mo_fcc;//mo
		Xcinetique[7] = x_c_fcc;//c
		Xcinetique[4] = 1-Xcinetique[5]-Xcinetique[6]-Xcinetique[7]; //fe
		Xcinetique[8] = Ycinetique[10];
		
		so1.ecriture ( Ycinetique );
		so3.ecritureX ( Xcinetique );
		prevVectorcin=Ycinetique;
		flage=true;
		return true;
	}
	else{flage=false; return false;}

}
