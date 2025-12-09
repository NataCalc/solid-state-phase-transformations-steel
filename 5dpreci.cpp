/*
 *  5dpreci.cpp
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

extern double temps,dt;
extern double T,dT;
extern float taille_gradient;
extern double dtemps;
extern float taille_gradient_alpha;
extern float taille_gradient_betta;
extern segment seg[20];
extern int nbseg;

extern sortie so;
extern sortie so1;

extern matrice mat;

extern distribution_moyenne<sphere> cem;
extern distribution_moyenne<sphere> eps;

void calc5dequilibre();

void calc5dequilibre()
{

	int dim=5;
	newt nt;
    qhull qh(dim);
	bool conv = false;
	std::vector <double> CONCENTRATION_NOMINAL ( 5 );
	CONCENTRATION_NOMINAL[0] = 9.3173484e-3;              		                     				              // carbon -C
	CONCENTRATION_NOMINAL[2] = 2.6903673e-2;  		                                     			              // chrome -CR
	CONCENTRATION_NOMINAL[3] = 8.1652564e-3;									                                  // Mo
	CONCENTRATION_NOMINAL[4] = 9.4947374e-2;									                                  // Co
	CONCENTRATION_NOMINAL[1] = 1.-CONCENTRATION_NOMINAL[0]-CONCENTRATION_NOMINAL[2]-CONCENTRATION_NOMINAL[3]-CONCENTRATION_NOMINAL[4];// iron -FE
	//::::::Thermo FE-C-CR-MO (BCC_A2 & FCC_A1):::::::::::
	printf ("Begin calcul system FE-CR-C-MO-CO (5D): \n");
	
	for(int iterseg=1;iterseg<=nbseg;iterseg++)
	{
		temps=seg[iterseg].date;
		dt=(seg[iterseg+1].date-temps)/seg[iterseg+1].pas;
		
		T=seg[iterseg].temperature;
		dT= dt * (seg[iterseg+1].temperature-seg[iterseg].temperature)
		/ (seg[iterseg+1].date-seg[iterseg].date);
		
		while(temps<seg[iterseg+1].date)
		{
			temps+=dt;
			T+=dT;
        	
			qh=qhull(dim);
			int j = 1;
			while ( mat.produit[j]!=NULL )
			{
				nt.func=mat.produit[1]->frontiere[1];
				nt.nrfuncv= &interface::SQ5;
				
				nt.jacobian=mat.produit[1]->frontiere[1];
				nt.jacobfunc = &interface::JQ5;
				
				qh.func1    =mat.produit[1]->frontiere[1];
				qh.nrfuncv1 =&interface::phase1_quinary;
				
				qh.func2    =mat.produit[1]->frontiere[1];
				qh.nrfuncv2 =&interface::phase2_quinary;
				
				j++;
			}
			if ( nt.func == NULL ) printf ( "newton.func = null\n" );
			if ( qh.func1== NULL ) printf ( "qhull.func1 = null\n" );
			
			//::::Data reading the files ferrite1.txt and cementite1.txt(Begin):::::
			mat.Thermo->read ( T );
			cem.classe0->Thermo->read ( T );
			//::::Data reading the files ferrite1.txt and cementite1.txt(End):::::
			
			qh.discretization (T,17);
			qh.build_hull();
			
			vector<point> ext=qh.extremums();
			printf ( "T = %lf \n", T );
			//cout<<"Extremums:"<<ext.size()<<endl;
			//cout<<qh<<endl;
			if ( ( !ext.empty() ) )
			{
				printf ( "extremums:" );
				//for ( int i=0;i<ext.size();i++ )
				//cout<<i+1<<":"<<ext[i]<<endl;
			}
			else
			{
				printf ( "No extremums T = %lf\n",T );
				continue;
			}
			
			nt.nn = 13;
			std::vector <double> Y ( nt.nn );
			nt.fvec.resize ( nt.nn );
			nt.fvec1.resize ( nt.nn ); //Jacobian
			
			Y[0] = 0.88014401;//1.-Y[1]-Y[2]-Y[3];//0.87915264;   //fe -
			Y[1] = 2.2694631E-2;//0.0211688;//2.3468693E-2;  //cr - y
			Y[2] = 6.5101503E-5;//0.00835335;//3.879751E-4; //mo - z
			Y[3] = 9.7096257E-2;//0.0971345;//9.6990691E-2;//co-d
			Y[4] = 2.6999152E-10;//1.E-10;//1.7474716E-6; //c -x
			Y[5] = 1.;//1.;//0.99999825;    //va
			
			Y[6] = 9.687299E-5;//1.-Y[7]-Y[8]-Y[9];//3.3265268E-3;   //fe
			Y[7] = 0.36766457;//0.470588;//0.33429753;   //cr - y
			Y[8] = 0.63223833;//1.e-10;//0.66233691; //mo - z
			Y[9] = 2.2493395E-7;//1.e-10;//3.9030294E-5;//co-d
			Y[10] =0.72711408;//0.705882;//0.79222387;   //c  - x
			Y[11] =0.27288592;//1.-Y[10];//0.20777613;   //va
			Y[12]= 0.45;
			
			/*
			 Y[1] = ext[ ( ext[0].function==0 ) ?0:1][1]; //cr - y
			 if ( Y[1]<=0 ) Y[1]=1.e-10;
			 Y[2] = ext[ ( ext[0].function==0 ) ?0:1][2]; //mo - z
			 if ( Y[2]<=0 ) Y[1]=1.e-10;
			 Y[3] = ext[ ( ext[0].function==0 ) ?0:1][3];//co-d
			 if ( Y[3]<=0 ) Y[1]=1.e-10;
			 Y[0] = 1.-Y[1]-Y[2]-Y[3]; //fe -
			 if ( Y[0]<=0 ) Y[1]=1.e-10;
			 Y[4] = ext[ ( ext[0].function==0 ) ?0:1][0]; ///c -x
			 if ( Y[4]<=0 ) Y[1]=1.e-10;
			 Y[5] = 1. - Y[4]; //va
			 
			 Y[7] = ext[ ( ext[0].function==0 ) ?1:0][1]; //cr - y
			 if ( Y[7]<=0 ) Y[1]=1.e-10;
			 Y[8] = ext[ ( ext[0].function==0 ) ?1:0][2]; //mo - z
			 if ( Y[8]<=0 ) Y[1]=1.e-10;
			 Y[9] = ext[ ( ext[0].function==0 ) ?1:0][3];//co-d
			 if ( Y[9]<=0 ) Y[1]=1.e-10;
			 Y[6] = 1.-Y[7]-Y[8]-Y[9]; //fe
			 if ( Y[6]<=0 ) Y[1]=1.e-10;
			 if ( Y[10]<=0 ) Y[1]=1.e-10;
			 Y[11] =1.-Y[10]; //va
			 Y[12]= 0.45;
			 */
			int check=0;
			nt.calc ( Y, &check,CONCENTRATION_NOMINAL,0,conv);
			
			printf ( "\n" );
			printf ( "Vector Y[i]:" );
			for ( int k=0 ; k < 13 ; k++ )
			{
				printf ( "%g ",Y[k] );
			}
			printf ( "\n" );
			so.ecriture ( Y );
        }
		
	}
	so.fermeture();
}