/*
 *  equilibre.cpp
 *  
 *
 *  Created by nat on 12/03/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "phase.hpp"
#include "donnees.h"

//::::::::::::::::::::::FE-C::::::::::::::::::::::::
void interface::SB(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(1); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(1);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	Y1alfa[0]  =  Y[0]; //bcc_A2: FE
	Y2alfa[0]  =  Y[1]; //bcc_A2: C
	Y2alfa[1]  =  Y[2]; //bcc_A2: VA
	Y1betta[0] =  Y[3]; //fcc_A1: FE
	Y2betta[0] =  Y[4]; //fcc_A1: C
	Y2betta[1] =  Y[5]; //fcc_A1: VA
	
	std::vector <double> pot_chemicalalfa(2); //bcc_A2 : 0 - "C"; 1 - "FE"
	std::vector <double> pot_chemicalbetta(2);//fcc_A1 : 0 - "C"; 1 - "FE"
	
	
	pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
	pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
	F[0] = (pot_chemicalbetta[0]-pot_chemicalalfa[0])/(R*TT2); //C_bcc - C_fcc
	F[1] = (pot_chemicalbetta[1]-pot_chemicalalfa[1])/(R*TT2); //Fe_bcc - Fe_fcc
	F[2] = Y[0] - 1.;        //Fe-1.
	F[3] = Y[3] - 1.;        //Fe-1.
	F[4] = Y[1] + Y[2] - 1.; //C+VA-1.
	F[5] = Y[4] + Y[5] - 1.; //C+VA-1.
	
	/*for (int i=0; i < 6; i++)
	 {printf("%e ", F[i]);
	 }printf("\n");*/	
}

//::::::::::::::::::::::FE-Ni::::::::::::::::::::::::
void interface::SB_FeNi(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2); //bcc_A2: FE
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: FE
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	Y1alfa[0]  =  Y[0]; //bcc_A2: FE
	Y1alfa[1]  =  Y[1]; //bcc_A2: Ni
	Y2alfa[0]  =  0.; //bcc_A2: C
	Y2alfa[1]  =  1.; //bcc_A2: Va

	Y1betta[0] =  Y[2]; //fcc_A1: FE
	Y1betta[1] =  Y[3]; //fcc_A1: Ni
	Y2betta[0] =  0.; //fcc_A1: C 
	Y2betta[1] =  1.; //fcc_A1: Va 

	std::vector <double> pot_chemicalalfa(3); //bcc_A2 : 0 - "C"; 1 - "FE"
	std::vector <double> pot_chemicalbetta(3);//fcc_A1 : 0 - "C"; 1 - "FE"
	

	pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
	pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
	//printf("pot_chemicalbetta[2] = %e pot_chemicalalpha[2]=%e \n",pot_chemicalbetta[2],pot_chemicalalfa[2]);
	//printf("pot_chemicalbetta[1] = %e pot_chemicalalpha[1]=%e \n",pot_chemicalbetta[1],pot_chemicalalfa[1]);

	F[0] = (pot_chemicalbetta[2]-pot_chemicalalfa[2])/(R*TT2); //C_bcc - C_fcc
	F[1] = (pot_chemicalbetta[1] - pot_chemicalalfa[1])/(R*TT2); //(Cr_bcc - Fe_bcc)-(Cr_fcc - Fe_fcc)
	//F[1] = (pot_chemicalbetta[1])/(R*TT2) - (pot_chemicalalfa[1])/(R*TT2); //(Cr_bcc - Fe_bcc)-(Cr_fcc - Fe_fcc)
	F[2] = Y[0] + Y[1] - 1.;        //Fe+Ni-1.
	F[3] = Y[2] + Y[3] - 1.;        //Fe+Ni-1.

	/*for (int i=0; i < 4; i++)
	 {printf("%e ", F[i]);
	 }printf("\n");*/	
}


//_________________________Jacobian FE-C____________________
void interface::JB(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
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
			SB(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			SB(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
}


//_________________________Jacobian FE-C____________________
void interface::JB_FeNi(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
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
			SB_FeNi(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			SB_FeNi(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
}


//_______________________system Fe-Cr-C______________________________
void interface::ST(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2); //bcc_A2: FE,CR
	std::vector <double> Y2alfa(2); //bcc_A2: C,VA
	std::vector <double> Y1betta(2);//fcc_A1: FE,CR
	std::vector <double> Y2betta(2);//fcc_A1: C,VA
	
	Y1alfa[0]  = Y[0];//bcc_A2: FE
	Y1alfa[1]  = Y[1];//bcc_A2: CR
	Y2alfa[0]  = Y[2];//bcc_A2: C
	Y2alfa[1]  = Y[3];//bcc_A2: VA
	Y1betta[0] = Y[4];//fcc_A1: FE
	Y1betta[1] = Y[5];//fcc_A1: CR
	Y2betta[0] = Y[6];//fcc_A1: C
	Y2betta[1] = Y[7];//fcc_A1: VA
	
	std::vector <double> pot_chemical1alfa(3);//bcc_A2: 0 - C; 1 - FE; 2 - CR
	std::vector <double> pot_chemical1betta(3);//fcc_A1: 0 - C; 1 - FE; 2 - CR
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	F[0] = (pot_chemical1alfa[0] - pot_chemical1betta[0])/(R*TT2); // C_bcc - C_fcc
	F[1] = (pot_chemical1alfa[1] - pot_chemical1betta[1])/(R*TT2); // Fe_bcc - F_fcc
	F[2] = (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2) - (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2); //(Cr_bcc - Fe_bcc)-(Cr_fcc - Fe_fcc)
	F[3] = Y[0] + Y[1] - 1.; //Fe+Cr=1
	F[4] = Y[4] + Y[5] - 1.; //Fe+Cr=1
	F[5] = Y[2] + Y[3] - 1.; //C+Va=1
	F[6] = Y[6] + Y[7] - 1.; //C+Va=1
	F[7] = (ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[3])*Y[8] \
	+(ph2->Thermo->m[0]*Y[5])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7])*(1-Y[8])-concentration[2]; //CR
	
	F[8] = (ph1->Thermo->m[1]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[3])*Y[8] \
	+(ph2->Thermo->m[1]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7])*(1-Y[8])-concentration[0]; //C
}

//________________________Jacobien Fe-Cr-C_______________________________
void interface::JT(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
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
		if(Y[i]==0.0) Y[i]=0.0000001;
	} 
	
	for( i=0 ; i < i1 ; i++)
	{
		Ynew[i] = Y[i];
	}  
	
	for( i=0 ; i < Jacobian3D.size() ; i++)
		for(j=0 ; j < Jacobian3D[i].size() ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			ST(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			ST(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}

void interface::Jacobian3D(std::vector <double> &Y, std::vector <double> &F, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(2);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(2);
	std::vector <double> Y2betta(2);
	int size = Y.size();
	Jacobian2D.resize(size);
	
	for (int i = 0; i < Jacobian2D.size(); i++)
	{
		Jacobian2D[i].resize(size);
	}
	
	Y1alfa[0] = Y[0];//fe
	Y1alfa[1] = Y[1];//cr
	Y2alfa[0] = Y[2];//c
	Y2alfa[1] = Y[3];//va
	
	Y1betta[0]= Y[4];//fe
	Y1betta[1]= Y[5];//cr
	Y2betta[0]= Y[6];//c 
	Y2betta[1]= Y[7];//va
	
	std::vector < std::vector<double> > Aalfa(3), Abetta(3);
	for(int i = 0 ; i < Aalfa.size() ; i++)
	{
		Aalfa[i].resize(4);
	}  
	for(int i = 0 ; i < Abetta.size() ; i++)
	{
		Abetta[i].resize(4);
	}  
	
	Aalfa =ph1->Thermo->chemical_potential2(Y1alfa,Y2alfa);
	Abetta=ph2->Thermo->chemical_potential2(Y1betta,Y2betta);
	
	
	Jacobian2D[0][0]=Aalfa[0][0];Jacobian2D[0][1]=Aalfa[0][1];Jacobian2D[0][2]=Aalfa[0][2];
	Jacobian2D[0][3]=Aalfa[0][3];Jacobian2D[0][4]=Aalfa[0][4];Jacobian2D[0][5]=-Abetta[0][0];
	Jacobian2D[0][6]=-Abetta[0][1];Jacobian2D[0][7]=-Abetta[0][2];Jacobian2D[0][8]=-Abetta[0][3];
	Jacobian2D[0][9]=-Abetta[0][4];Jacobian2D[0][10]=0.0;Jacobian2D[1][0]=Aalfa[1][0];
	Jacobian2D[1][1]=Aalfa[1][1];Jacobian2D[1][2]=Aalfa[1][2];Jacobian2D[1][3]=Aalfa[1][3];
	Jacobian2D[1][4]=Aalfa[1][4];Jacobian2D[1][5]=-Abetta[1][0];Jacobian2D[1][6]=-Abetta[1][1];
	Jacobian2D[1][7]=-Abetta[1][2];Jacobian2D[1][8]=-Abetta[1][3];Jacobian2D[1][9]=-Abetta[1][4];
	Jacobian2D[1][10]=0.0;Jacobian2D[2][0]=Aalfa[2][0]-Aalfa[1][0];Jacobian2D[2][1]=Aalfa[2][1]-Aalfa[1][1];
	Jacobian2D[2][2]=Aalfa[2][2]-Aalfa[1][2];Jacobian2D[2][3]=Aalfa[2][3]-Aalfa[1][3];
	Jacobian2D[2][4]=Aalfa[2][4]-Aalfa[1][4];Jacobian2D[2][5]=-Abetta[2][0]+Abetta[1][0];
	Jacobian2D[2][6]=-Abetta[2][1]+Abetta[1][1];Jacobian2D[2][7]=-Abetta[2][2]+Abetta[1][2];
	Jacobian2D[2][8]=-Abetta[2][3]+Abetta[1][3];Jacobian2D[2][9]=-Abetta[2][4]+Abetta[1][4];Jacobian2D[2][10]=0.0;
	Jacobian2D[3][0]=Aalfa[3][0]-Aalfa[1][0];Jacobian2D[3][1]=Aalfa[3][1]-Aalfa[1][1];
	Jacobian2D[3][2]=Aalfa[3][2]-Aalfa[1][2];Jacobian2D[3][3]=Aalfa[3][3]-Aalfa[1][3];
	Jacobian2D[3][4]=Aalfa[3][4]-Aalfa[1][4];Jacobian2D[3][5]=-Abetta[3][0]+Abetta[1][0];
	Jacobian2D[3][6]=-Abetta[3][1]+Abetta[1][1];Jacobian2D[3][7]=-Abetta[3][2]+Abetta[1][2];
	Jacobian2D[3][8]=-Abetta[3][3]+Abetta[1][3];Jacobian2D[3][9]=-Abetta[3][4]+Abetta[1][4];
	Jacobian2D[3][10]=0.0;Jacobian2D[4][0]=1.0;Jacobian2D[4][1]=1.0;Jacobian2D[4][2]=1.0;
	Jacobian2D[4][3]=0.0;Jacobian2D[4][4]=0.0;Jacobian2D[4][5]=0.0;Jacobian2D[4][6]=0.0;
	Jacobian2D[4][7]=0.0;Jacobian2D[4][8]=0.0;Jacobian2D[4][9]=0.0;Jacobian2D[4][10]=0.0;
	Jacobian2D[5][0]=0.0;Jacobian2D[5][1]=0.0;Jacobian2D[5][2]=0.0;Jacobian2D[5][3]=1.0;
	Jacobian2D[5][4]=1.0;Jacobian2D[5][5]=0.0;Jacobian2D[5][6]=0.0;Jacobian2D[5][7]=0.0;
	Jacobian2D[5][8]=0.0;Jacobian2D[5][9]=0.0;Jacobian2D[5][10]=0.0;
	Jacobian2D[6][0]=0.0;Jacobian2D[6][1]=0.0;Jacobian2D[6][2]=0.0;Jacobian2D[6][3]=0.0;
	Jacobian2D[6][4]=0.0;Jacobian2D[6][5]=1.0;Jacobian2D[6][6]=1.0;Jacobian2D[6][7]=1.0;
	Jacobian2D[6][8]=0.0;Jacobian2D[6][9]=0.0;Jacobian2D[6][10]=0.0;
	Jacobian2D[7][0]=0.0;Jacobian2D[7][1]=0.0;Jacobian2D[7][2]=0.0;Jacobian2D[7][3]=0.0;
	Jacobian2D[7][4]=0.0;Jacobian2D[7][5]=0.0;Jacobian2D[7][6]=0.0;Jacobian2D[7][7]=0.0;
	Jacobian2D[7][8]=1.0;Jacobian2D[7][9]=1.0;Jacobian2D[7][10]=0.0;
	Jacobian2D[8][0]=0.0;Jacobian2D[8][1]=ph1->Thermo->m[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10];
	Jacobian2D[8][2]=0.0;Jacobian2D[8][3]=0.0;
	Jacobian2D[8][4]=(ph1->Thermo->m[0]*ph1->Thermo->m[1]*Y[1])/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]),2.)*Y[10];
	Jacobian2D[8][5]=0.0;Jacobian2D[8][6]=ph2->Thermo->m[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10]);
	Jacobian2D[8][7]=0.0;Jacobian2D[8][8]=0.0;
	Jacobian2D[8][9]=(ph2->Thermo->m[0]*ph2->Thermo->m[1]*Y[6])/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]),2.)*(1-Y[10]);
	Jacobian2D[8][10]=(ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]) \
	-(ph2->Thermo->m[0]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]);
	Jacobian2D[9][0]=0.0;Jacobian2D[9][1]=0.0;
	Jacobian2D[9][2]=ph1->Thermo->m[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10];
	Jacobian2D[9][3]=0.0;
	Jacobian2D[9][4]=(ph1->Thermo->m[0]*ph1->Thermo->m[1]*Y[2])/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]),2.)*Y[10];
	Jacobian2D[9][5]=0.0;Jacobian2D[9][6]=0.0;Jacobian2D[9][7]=ph2->Thermo->m[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10]);
	Jacobian2D[9][8]=0.0;
	Jacobian2D[9][9]=(ph2->Thermo->m[0]*ph2->Thermo->m[1]*Y[7])/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]),2.)*(1-Y[10]);
	Jacobian2D[9][10]=(ph1->Thermo->m[0]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]) \
	-(ph2->Thermo->m[0]*Y[7])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]);
	Jacobian2D[10][0]=0.0;Jacobian2D[10][1]=0.0;Jacobian2D[10][2]=0.0;
	Jacobian2D[10][3]=ph1->Thermo->m[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10];
	Jacobian2D[10][4]=(ph1->Thermo->m[1]*ph1->Thermo->m[1]*Y[3])/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]),2.)*Y[10];
	Jacobian2D[10][5]=0.0;Jacobian2D[10][6]=0.0;Jacobian2D[10][7]=0.0;
	Jacobian2D[10][8]=ph2->Thermo->m[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10]);
	Jacobian2D[10][9]=(ph2->Thermo->m[1]*ph2->Thermo->m[1]*Y[8])/pow((ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]),2.)*(1-Y[10]);
	Jacobian2D[10][10]= (ph1->Thermo->m[1]*Y[3])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]) \
	-(ph2->Thermo->m[1]*Y[8])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]);
	
	for (int i = 0; i < 10; i++)
	{
		Jacobian2D[0][i] =Jacobian2D[0][i]/(R*TT2);
		Jacobian2D[1][i] =Jacobian2D[1][i]/(R*TT2);
		Jacobian2D[2][i] =Jacobian2D[2][i]/(R*TT2);
		Jacobian2D[3][i] =Jacobian2D[3][i]/(R*TT2);
	}
	
}


//_________________________system Fe-Cr-Mo-C________________________
void interface::SQ4(int,std::vector <double> &Y,std::vector <double> &F,std::vector <double> &concentration)
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
	F[8] = (ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10] \
	+(ph2->Thermo->m[0]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10])-concentration[2]; //CR
	
	F[9] = (ph1->Thermo->m[0]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10] \
	+(ph2->Thermo->m[0]*Y[7])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10])-concentration[3]; //Mo
	
	F[10] = (ph1->Thermo->m[1]*Y[3])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10] \
	+(ph2->Thermo->m[1]*Y[8])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10])-concentration[0]; //C
	
}


//_________________________Jacobian Fe-Cr-Mo-C__________________________
void interface::JQ4(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian4D,std::vector <double> &concentration)
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
			SQ4(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			SQ4(n,Y,F,concentration); 
			Jacobian4D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}


void interface::Jacobian4D(std::vector <double> &Y, std::vector <double> &F, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(3);
	std::vector <double> Y2alfa(2);
	std::vector <double> Y1betta(3);
	std::vector <double> Y2betta(2);
	
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
	
	std::vector < std::vector<double> > Aalfa(4), Abetta(4);
	for(int i = 0 ; i < Aalfa.size() ; i++)
	{
		Aalfa[i].resize(5);
	}  
	for(int i = 0 ; i < Abetta.size() ; i++)
	{
		Abetta[i].resize(5);
	} 
	
	Aalfa =ph1->Thermo->chemical_potential2(Y1alfa,Y2alfa);
	Abetta=ph2->Thermo->chemical_potential2(Y1betta,Y2betta);
	
	
	int size = Y.size();
	Jacobian2D.resize(size);
	
	for (int i = 0; i < Jacobian2D.size(); i++)
	{
		Jacobian2D[i].resize(size);
	}
	
	
	Jacobian2D[0][0]=Aalfa[0][0];Jacobian2D[0][1]=Aalfa[0][1];Jacobian2D[0][2]=Aalfa[0][2];
	Jacobian2D[0][3]=Aalfa[0][3];Jacobian2D[0][4]=Aalfa[0][4];Jacobian2D[0][5]=-Abetta[0][0];
	Jacobian2D[0][6]=-Abetta[0][1];Jacobian2D[0][7]=-Abetta[0][2];Jacobian2D[0][8]=-Abetta[0][3];
	Jacobian2D[0][9]=-Abetta[0][4];Jacobian2D[0][10]=0.0;Jacobian2D[1][0]=Aalfa[1][0];
	Jacobian2D[1][1]=Aalfa[1][1];Jacobian2D[1][2]=Aalfa[1][2];Jacobian2D[1][3]=Aalfa[1][3];
	Jacobian2D[1][4]=Aalfa[1][4];Jacobian2D[1][5]=-Abetta[1][0];Jacobian2D[1][6]=-Abetta[1][1];
	Jacobian2D[1][7]=-Abetta[1][2];Jacobian2D[1][8]=-Abetta[1][3];Jacobian2D[1][9]=-Abetta[1][4];
	Jacobian2D[1][10]=0.0;Jacobian2D[2][0]=Aalfa[2][0]-Aalfa[1][0];Jacobian2D[2][1]=Aalfa[2][1]-Aalfa[1][1];
	Jacobian2D[2][2]=Aalfa[2][2]-Aalfa[1][2];Jacobian2D[2][3]=Aalfa[2][3]-Aalfa[1][3];
	Jacobian2D[2][4]=Aalfa[2][4]-Aalfa[1][4];Jacobian2D[2][5]=-Abetta[2][0]+Abetta[1][0]; 
	Jacobian2D[2][6]=-Abetta[2][1]+Abetta[1][1];Jacobian2D[2][7]=-Abetta[2][2]+Abetta[1][2];
	Jacobian2D[2][8]=-Abetta[2][3]+Abetta[1][3];Jacobian2D[2][9]=-Abetta[2][4]+Abetta[1][4];
	Jacobian2D[2][10]=0.0;Jacobian2D[3][0]=Aalfa[3][0]-Aalfa[1][0];Jacobian2D[3][1]=Aalfa[3][1]-Aalfa[1][1];
	Jacobian2D[3][2]=Aalfa[3][2]-Aalfa[1][2];Jacobian2D[3][3]=Aalfa[3][3]-Aalfa[1][3];
	Jacobian2D[3][4]=Aalfa[3][4]-Aalfa[1][4];Jacobian2D[3][5]=-Abetta[3][0]+Abetta[1][0];
	Jacobian2D[3][6]=-Abetta[3][1]+Abetta[1][1];Jacobian2D[3][7]=-Abetta[3][2]+Abetta[1][2];
	Jacobian2D[3][8]=-Abetta[3][3]+Abetta[1][3];Jacobian2D[3][9]=-Abetta[3][4]+Abetta[1][4];
	Jacobian2D[3][10]=0.0;Jacobian2D[4][0]=1.0;Jacobian2D[4][1]=1.0;Jacobian2D[4][2]=1.0;
	Jacobian2D[4][3]=0.0;Jacobian2D[4][4]=0.0;Jacobian2D[4][5]=0.0;Jacobian2D[4][6]=0.0;
	Jacobian2D[4][7]=0.0;Jacobian2D[4][8]=0.0;Jacobian2D[4][9]=0.0;Jacobian2D[4][10]=0.0;
	Jacobian2D[5][0]=0.0;Jacobian2D[5][1]=0.0;Jacobian2D[5][2]=0.0;Jacobian2D[5][3]=1.0;
	Jacobian2D[5][4]=1.0;Jacobian2D[5][5]=0.0;Jacobian2D[5][6]=0.0;Jacobian2D[5][7]=0.0;
	Jacobian2D[5][8]=0.0;Jacobian2D[5][9]=0.0;Jacobian2D[5][10]=0.0;
	Jacobian2D[6][0]=0.0;Jacobian2D[6][1]=0.0;Jacobian2D[6][2]=0.0;Jacobian2D[6][3]=0.0;
	Jacobian2D[6][4]=0.0;Jacobian2D[6][5]=1.0;Jacobian2D[6][6]=1.0;Jacobian2D[6][7]=1.0;
	Jacobian2D[6][8]=0.0;Jacobian2D[6][9]=0.0;Jacobian2D[6][10]=0.0;
	Jacobian2D[7][0]=0.0;Jacobian2D[7][1]=0.0;Jacobian2D[7][2]=0.0;Jacobian2D[7][3]=0.0;
	Jacobian2D[7][4]=0.0;Jacobian2D[7][5]=0.0;Jacobian2D[7][6]=0.0;Jacobian2D[7][7]=0.0;
	Jacobian2D[7][8]=1.0;Jacobian2D[7][9]=1.0;Jacobian2D[7][10]=0.0;
	Jacobian2D[8][0]=0.0;Jacobian2D[8][1]=ph1->Thermo->m[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10];
	Jacobian2D[8][2]=0.0;Jacobian2D[8][3]=0.0;Jacobian2D[8][4]=(ph1->Thermo->m[0]*ph1->Thermo->m[1]*Y[1])/pow((ph1->Thermo->m[0] \
																											   +ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]),2.)*Y[10];
	Jacobian2D[8][5]=0.0;Jacobian2D[8][6]=ph2->Thermo->m[0]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10]);
	Jacobian2D[8][7]=0.0;Jacobian2D[8][8]=0.0;Jacobian2D[8][9]=(ph2->Thermo->m[0]*ph2->Thermo->m[1]*Y[6])/pow((ph2->Thermo->m[0] \
																											   +ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]),2.)*(1-Y[10]);
	Jacobian2D[8][10]=(ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]) \
	-(ph2->Thermo->m[0]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]);
	Jacobian2D[9][0]=0.0;Jacobian2D[9][1]=0.0;Jacobian2D[9][2]=ph1->Thermo->m[0]/(ph1->Thermo->m[0]+ph1->Thermo->m[1] \
																				  -ph1->Thermo->m[1]*Y[4])*Y[10];Jacobian2D[9][3]=0.0;
	Jacobian2D[9][4]=(ph1->Thermo->m[0]*ph1->Thermo->m[1]*Y[2])/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]),2.)*Y[10];
	Jacobian2D[9][5]=0.0;Jacobian2D[9][6]=0.0;Jacobian2D[9][7]=ph2->Thermo->m[0]/(ph2->Thermo->m[0] \
																				  +ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10]);
	Jacobian2D[9][8]=0.0;Jacobian2D[9][9]=(ph2->Thermo->m[0]*ph2->Thermo->m[1]*Y[7])/pow((ph2->Thermo->m[0] \
																						  +ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]),2.)*(1-Y[10]);
	Jacobian2D[9][10]=(ph1->Thermo->m[0]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]) \
	-(ph2->Thermo->m[0]*Y[7])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]);
	Jacobian2D[10][0]=0.0;Jacobian2D[10][1]=0.0;Jacobian2D[10][2]=0.0;
	Jacobian2D[10][3]=ph1->Thermo->m[1]/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4])*Y[10];
	Jacobian2D[10][4]=(ph1->Thermo->m[1]*ph1->Thermo->m[1]*Y[3])/pow((ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]),2.)*Y[10];
	Jacobian2D[10][5]=0.0;Jacobian2D[10][6]=0.0;Jacobian2D[10][7]=0.0;
	Jacobian2D[10][8]=ph2->Thermo->m[1]/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9])*(1-Y[10]);
	Jacobian2D[10][9]=(ph2->Thermo->m[1]*ph2->Thermo->m[1]*Y[8])/pow((ph2->Thermo->m[0] \
																	  +ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]),2.)*(1-Y[10]);
	Jacobian2D[10][10]= (ph1->Thermo->m[1]*Y[3])/(ph1->Thermo->m[0]+ph1->Thermo->m[1] \
												  -ph1->Thermo->m[1]*Y[4])-(ph2->Thermo->m[1]*Y[8])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]);
	
	for (int i = 0; i < 10; i++)
	{
		Jacobian2D[0][i] =Jacobian2D[0][i]/(R*TT2);
		Jacobian2D[1][i] =Jacobian2D[1][i]/(R*TT2);
		Jacobian2D[2][i] =Jacobian2D[2][i]/(R*TT2);
		Jacobian2D[3][i] =Jacobian2D[3][i]/(R*TT2);
	}
	
}


//______________________________system Fe-Cr-Mo-Co-C-Va________________________________
void interface::SQ5(int,std::vector <double> &Y,std::vector <double> &F,std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	std::vector <double> Y1alfa(4);//bcc_A2: FE,CR,MO,CO
	std::vector <double> Y2alfa(2);//bcc_A2: C,VA
	std::vector <double> Y1betta(4);//fcc_A1:FE,CR,MO,CO
	std::vector <double> Y2betta(2);//fcc_A1:C,VA
	
	Y1alfa[0] = Y[0];//fe
	Y1alfa[1] = Y[1];//cr
	Y1alfa[2] = Y[2];//mo
	Y1alfa[3] = Y[3];//co
	Y2alfa[0] = Y[4];//c
	Y2alfa[1] = Y[5];//va
	
	Y1betta[0]= Y[6];//fe
	Y1betta[1]= Y[7];//cr
	Y1betta[2]= Y[8];//mo
	Y1betta[3]= Y[9];//co
	Y2betta[0]= Y[10];//c 
	Y2betta[1]= Y[11];//va
	
	std::vector <double> pot_chemical1alfa(5);//bcc_A2: 0 - C; 1 - FE, 2 - CR; 3 - MO; 4 - CO
	std::vector <double> pot_chemical1betta(5);//fcc_A1: 0 - C; 1 - FE; 2 - CR; 3 - MO; 4 - CO
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);//bcc_A2
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);//fcc_A1
	
	
	F[0] = (pot_chemical1alfa[0] - pot_chemical1betta[0])/(R*TT2); // C_bcc - C_fcc
	F[1] = (pot_chemical1alfa[1] - pot_chemical1betta[1])/(R*TT2); // Fe_bcc - Fe_fcc
	F[2] = (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2) - (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2); //(Cr_bcc-Fe_bcc)-(Cr_fcc-Fe_fcc)
	F[3] = (pot_chemical1alfa[3] - pot_chemical1alfa[1])/(R*TT2) - (pot_chemical1betta[3] - pot_chemical1betta[1])/(R*TT2); //(Mo_bcc-Fe_bcc)-(Mo_fcc-Fe_fcc)
	F[4] = (pot_chemical1alfa[4] - pot_chemical1alfa[1])/(R*TT2) - (pot_chemical1betta[4] - pot_chemical1betta[1])/(R*TT2); //(Co_bcc-Fe_bcc)-(Co_fcc-Fe_fcc)
	F[5] = Y[0] + Y[1] + Y[2] + Y[3]- 1.; //Fe+Cr+Mo+Co=1
	F[6] = Y[4] + Y[5] - 1.; //C+Va=1
	F[7] = Y[6] + Y[7] + Y[8] + Y[9] - 1.; //Fe+Cr+Mo+Co=1
	F[8] = Y[10] + Y[11] - 1.; //C+Va=1
	F[9] = (ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[5])*Y[12] \
	+(ph2->Thermo->m[0]*Y[7])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[11])*(1-Y[12])-concentration[2]; //CR
	
	F[10] = (ph1->Thermo->m[0]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[5])*Y[12] \
	+(ph2->Thermo->m[0]*Y[8])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[11])*(1-Y[12])-concentration[3]; //Mo	
	
	F[11] = (ph1->Thermo->m[0]*Y[3])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[5])*Y[12] \
	+(ph2->Thermo->m[0]*Y[9])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[11])*(1-Y[12])-concentration[4]; //Co
	
	F[12] = (ph1->Thermo->m[1]*Y[4])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[5])*Y[12] \
	+(ph2->Thermo->m[1]*Y[10])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[11])*(1-Y[12])-concentration[0]; //C
	
}


//___________________________Jacobien Fe-Cr-Mo-Co-C-Va______________________________
void interface::JQ5(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
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
			SQ5(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			SQ5(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}
