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
#include <map>

extern double temps,dt;
extern double T,dT;
extern float taille_gradient;
extern double dtemps;
extern float taille_gradient_alpha;
extern float taille_gradient_betta;
extern segment seg[20];
extern int nbseg;
extern double correc_D;
extern float CIN;
extern int def;


extern sortie so;
extern sortie so1;
extern sortie so2;
extern sortie so3;

extern matrice mat;

extern distribution_moyenne<sphere> cem;
extern distribution_moyenne<sphere> eps;

void calc3dequilibreFeMoC();
bool calc3cinetique1FeMoC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);
bool calc3cinetique2FeMoC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage);


void calc3dequilibreFeMoC()
{	
	int dim=3;
	def=5;
	
	newt nt;
    qhull qh(dim);
	
	double counter=time(0);
	double prevT=0;
	bool conv = false;
	bool flage = false;

	std::vector<double> prevVector(9);		
	std::vector<double> prevVectorcin(9);
	for(int i=0;i<9;i++){prevVectorcin[i]=0;}		
	for(int i=0;i<9;i++){prevVector[i]=0;}	
	
	std::vector <double> CONCENTRATION_NOMINAL ( 3 );


	CONCENTRATION_NOMINAL[0] = 9.2319989e-3;              		                    // carbon -C w(c)=2.e-3
	CONCENTRATION_NOMINAL[2] = 1.1557801e-5;										// nickel -Mo w(mo)=2.e-5
	CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0]-CONCENTRATION_NOMINAL[2];// iron -FE
	printf ( "Begin calcul the system FE-2.e-5Mo-2.e-3C (3D): \n" );
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
			if (( T>=1109. ) && ( T<=1771.)) continue;
			if (prevT==0)
			{prevT=T;}
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
			
			if ((prevExts) && (good))
			{
				Y=prevVector;
			}
			else
			{
				
				qh.discretization (T,25);
				qh.build_hull();
				vector<point> ext=qh.extremums(2.e-5);
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
				Y[1] = ( ext[0].func()?ext[1][1]:ext[0][1] );  //Ni
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
				Y[5] =( ext[0].func()?ext[0][1]:ext[1][1] );   //Ni
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

		

			/*
			Y[0]=0.99998816; //Fe
			Y[1]=1.183555e-5; //Mo
			Y[2]=2.9496377e-4; //C
			Y[3]=0.99970504; //Va
			
			Y[4]=0.99998886; //Fe
			Y[5]=1.1135827e-5; //Mo
			Y[6]=3.5584933e-2; //C
			Y[7]=0.96441507; //Va
			Y[8]=0.45;*/
			

	printf ( "\n" );printf ( "Vector Y[i] before:" );
	for ( int k=0 ; k < 9 ; k++ )
	{printf ( "%g ",Y[k] );}printf ( "\n" );
	
	int check=0;
	good=nt.calc ( Y, &check,CONCENTRATION_NOMINAL,3,conv);
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
			double x_c_fcc,x_mo_fcc,x_c_bcc,x_mo_bcc;			
			x_mo_bcc = Y[1] / (1. + 3. * Y[2]) ;
			x_c_bcc  = 3. * Y[2] / (1. + 3. * Y[2]) ;
			x_mo_fcc = Y[5] / (1. + Y[6]) ;
			x_c_fcc  = Y[6] / (1. + Y[6]) ;			
			
			X[1] = x_mo_bcc;//mo			
			X[2] = x_c_bcc;//c
			X[0] = 1-X[1]-X[2];//fe
			X[4] = x_mo_fcc;//mo			
			X[5] = x_c_fcc;//c
			X[3] = 1-X[4]-X[5]; //fe
			X[6] = Y[8];
			
			
			
	if (good)
	{
		so.ecriture ( Y );
		so2.ecritureX ( X );
		prevVector=Y;
		prevT=T;
	}
	

/***************************************Cinetique*************************************************************/
	if(good)
		{
			if((T>=600.) && (T<=1100.))
			calc3cinetique1FeMoC(CONCENTRATION_NOMINAL,Y,prevVectorcin,flage);
			//calc3cinetique2FeNiC(CONCENTRATION_NOMINAL,Y,prevVectorcin, flage);
		}		
			
	}

}
	
    so1.fermeture();
	so2.fermeture();
	so3.fermeture();
    so.fermeture();
    counter=time(0)-counter;
    cout<<"time3d:"<<counter<<endl;
}

bool calc3cinetique1FeMoC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin,bool &flage)
{
	newt nt1;
	nt1=newt();
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

	
	
/*****************************************************************************************/
	if (flage)
	{
		Ycinetique=prevVectorcin;
		taille_gradient=taille_gradient+1.e-9*prevVectorcin[8]*dtemps;
	}
	else
	{		
	Ycinetique[0] = Y[0];			//FE
	Ycinetique[1] = 1.-Y[0];		//Mo
	Ycinetique[2] = Y[2];			//C
	Ycinetique[3] = 1.-Y[2];        //VA
	Ycinetique[4] = Y[4];			//FE
	Ycinetique[5] = 1.0-Y[4];       //Mo
	Ycinetique[6] = Y[6];			//C
	Ycinetique[7] = 1.-Y[6];        //VA
	Ycinetique[8] = v_ini;
	}	
	
	
	
/***************======================================================*********************/
	double dc , dcr;
	
	dc =  1.5e-5 * exp( -1.421e5 /R / T ) ;
	dcr = 3.e-5 * exp( 3.e5 /R/ T ) ; 
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
		//convergance FRACTION_SITE dans FRACTION_MOLAIRE
		std::vector <double> Xcinetique(7);
		double x_c_fcc,x_mo_fcc,x_c_bcc,x_mo_bcc;			
		x_mo_bcc = Ycinetique[1] / (1. + 3. * Ycinetique[2]) ;
		x_c_bcc  = 3. * Ycinetique[2] / (1. + 3. * Ycinetique[2]) ;
		x_mo_fcc = Ycinetique[5] / (1. + Ycinetique[6]) ;
		x_c_fcc  = Ycinetique[6] / (1. + Ycinetique[6]) ;			
		
		Xcinetique[1] = x_mo_bcc;//mo			
		Xcinetique[2] = x_c_bcc;//c
		Xcinetique[0] = 1-Xcinetique[1]-Xcinetique[2];//fe
		Xcinetique[4] = x_mo_fcc;//mo		
		Xcinetique[5] = x_c_fcc;//c
		Xcinetique[3] = 1-Xcinetique[4]-Xcinetique[5]; //fe
		Xcinetique[6] = Ycinetique[8];

		
		prevVectorcin=Ycinetique;
		printf("\n*******Convergence\n");
		so1.ecriture ( Ycinetique );
		so3.ecritureX ( Xcinetique );
		flage = true;
		return flage;
	}
	printf("\n****Pas de convergence*******\n");
	flage = false;
	return flage;
}	


bool calc3cinetique2FeMoC(std::vector <double> &CONCENTRATION_NOMINAL,std::vector <double> &Y,std::vector<double> &prevVectorcin, bool &flage)
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
		Ycinetique[1] = 1.-Y[0];       //MO
		Ycinetique[2] = Y[2];  //C
		Ycinetique[3] = 1.-Y[2];       //VA
		Ycinetique[4] = Y[4];			//FE
		Ycinetique[5] = 1.0-Y[4];      //MO
		Ycinetique[6] = Y[6];  //C
		Ycinetique[7] = 1.-Y[6];       //VA
		Ycinetique[8] = v_ini;
	}	
	
	
	double dc , dcr;
	
	dc =  1.5e-5 * exp( -1.421e5 /R / T ) ;
	dcr = 3.e-5 * exp( 3.e5 /R/ T ) ; 
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

