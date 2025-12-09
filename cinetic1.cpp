/*
 *  cinetic1.cpp
 *  
 *
 *  Di(T)=Di^0 * exp (-Qi/RT)
 *  delta = delta_old + Vitess * dtemps
 *  
 */

#include "phase.hpp"
#include "donnees.h"

extern float Diff[10];
extern float taille_gradient;
extern float taille_gradient_alpha;
extern float taille_gradient_betta;
extern double correc_D;
extern float CIN;


std::vector <double> interface::calc_diffusion_alpha(std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	
	int n = ph1->Thermo->Dex.size();
	std::vector <double> Diffalpha(n);
	
	for(int i=0;i<n;i++)
	{
		Diffalpha[i]=ph1->Thermo->Dex[i]*exp(-ph1->Thermo->Qex[i]/R/TT2);
	}
	return Diffalpha;
}

std::vector <double> interface::calc_diffusion_betta(std::vector <double> &concentration)
{
	double TT2 = 0.;
	TT2 = ph2->Thermo->T1;
	
	int n = ph2->Thermo->Dex.size();
	std::vector <double> Diffbetta(n);

	
	for(int i=0;i<n;i++){
		Diffbetta[i]=ph2->Thermo->Dex[i]*exp(-ph2->Thermo->Qex[i]/R/TT2);
	}
	
	return Diffbetta;
}

//___________________kinetic transition Fe-C_____________________
void interface::KTB1(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
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
	
	int dim = ph1->Thermo->Dex.size();
	std::vector <double> Diffalpha(dim);
	std::vector <double> Diffbetta(dim);

	Diffalpha = calc_diffusion_alpha(concentration);
	Diffbetta = calc_diffusion_betta(concentration);
	
	
	pot_chemicalalfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);// ph1->Thermo : bcc_A2
	pot_chemicalbetta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);// ph2->Thermo : fcc_A1
	
	F[0] = (pot_chemicalalfa[0]-pot_chemicalbetta[0])/(R*TT2); //C_bcc - C_fcc
	F[1] = (pot_chemicalalfa[1]-pot_chemicalbetta[1])/(R*TT2); //Fe_bcc - Fe_fcc
	F[2] = Y[0] - 1.;        //Fe-1.
	F[3] = Y[3] - 1.;        //Fe-1.
	F[4] = Y[1] + Y[2] - 1.; //C+VA-1.
	F[5] = Y[4] + Y[5] - 1.; //C+VA-1.	
	

	F[6]= (((ph2->Thermo->m[1]*Y[4])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[5]))*(Y[6]-(Diffbetta[0]/taille_gradient_betta)) \
		-((ph1->Thermo->m[1]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[2]))*(Y[6]+(Diffalpha[0]/taille_gradient_alpha)) \
		+concentration[0]*(Diffbetta[0]/taille_gradient_betta + Diffalpha[0]/taille_gradient_alpha)); //C 
		
	//for (int i=0; i < 7; i++)
	 //{printf("%e ", F[i]);
	 //}printf("\n");
}

//____________________________Jacobien kinetic transition Fe-C_____________________
void interface::JKB1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian2D,std::vector <double> &concentration)
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
			KTB1(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTB1(n,Y,F,concentration); 
			Jacobian2D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}		
}	

//_________________________kinetic transition Fe-Cr-C____________________
void interface::KTT1(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
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
	
	int dim = ph1->Thermo->Dex.size();
	std::vector <double> Diffalpha(dim);
	std::vector <double> Diffbetta(dim);
	
	Diffalpha = calc_diffusion_alpha(concentration);
	Diffbetta = calc_diffusion_betta(concentration);
	//printf("Diffalpha_c =%e Diffalpha_cr=%e \n",Diffalpha[0],Diffalpha[1]);
	//printf("Diffbetta_c =%e Diffbetta_cr=%e \n",Diffbetta[0],Diffbetta[1]);
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(Y1alfa,Y2alfa);
	pot_chemical1betta =  ph2->Thermo->chemical_potential(Y1betta,Y2betta);
	
	F[0] = (pot_chemical1alfa[0] - pot_chemical1betta[0])/(R*TT2); // carbone
	F[1] = (pot_chemical1alfa[1] - pot_chemical1betta[1])/(R*TT2); // iron
	F[2] = (pot_chemical1alfa[2] - pot_chemical1alfa[1])/(R*TT2) - (pot_chemical1betta[2] - pot_chemical1betta[1])/(R*TT2); //cr-fe
	F[3] = Y[0] + Y[1] - 1; //printf("%i\n",F[3]);
	F[4] = Y[4] + Y[5] - 1; //printf("%i\n",F[4]);
	F[5] = Y[2] + Y[3] - 1; //printf("%i\n",F[5]);
	F[6] = Y[6] + Y[7] - 1; //printf("%i\n",F[6]);
	
	std::vector <double> delta(2);
	delta = Zenerplane(Y,concentration);

	F[7]=  ( (ph2->Thermo->m[1]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7]) )*( Y[8]-(Diffbetta[1]*correc_D/taille_gradient_betta/delta[0])/CIN ) \
	-( (ph1->Thermo->m[1]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[3]) )*( Y[8]+(Diffalpha[1]*correc_D/taille_gradient_alpha/delta[0])/CIN ) \
	+concentration[0]*( (Diffbetta[1]*correc_D/taille_gradient_betta/delta[0])/CIN + (Diffalpha[1]*correc_D/taille_gradient_alpha/delta[0])/CIN ); //C 
	
	F[8] = ( (ph2->Thermo->m[1]*Y[5])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[7]) )*( Y[8]-(Diffbetta[1]/taille_gradient_betta/delta[1])/CIN ) \
	-( (ph1->Thermo->m[1]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[3]) )*( Y[8]+(Diffalpha[1]/taille_gradient_alpha/delta[1])/CIN ) \
	+concentration[2]*( (Diffbetta[1]/taille_gradient_betta/delta[1])/CIN + (Diffalpha[1]/taille_gradient_alpha/delta[1])/CIN )/*Diffbetta[1]/(Diffbetta[1]*correc_D)*/; //CR
	

}

//_______________Jacobien cinetique Fe-Cr-C____________________
void interface::JKT1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian3D,std::vector <double> &concentration)
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
			KTT1(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KTT1(n,Y,F,concentration); 
			Jacobian3D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
}

//_________________________kinetic transition Fe-Cr-Mo-C____________________
void interface::KT4_1(int n, std::vector <double> &Y, std::vector <double> &F,std::vector <double> &concentration)
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
	
	int dim = ph1->Thermo->Dex.size();
	std::vector <double> Diffalpha(dim);
	std::vector <double> Diffbetta(dim);
	
	Diffalpha = calc_diffusion_alpha(concentration);
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
	
	F[8] = -((ph1->Thermo->m[0]*Y[1])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]))*Y[10] \
	+((ph2->Thermo->m[0]*Y[6])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]))*(Y[10]-(Diff[0]/taille_gradient)) \
	+(Diff[0]/taille_gradient)*concentration[2]/**(Diffbetta[0]/Diffbetta[1])*/; //CR
	
	F[9] = -((ph1->Thermo->m[0]*Y[2])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]))*Y[10] \
	+((ph2->Thermo->m[0]*Y[7])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]))*(Y[10]-(Diff[0]/taille_gradient)) \
	+(Diff[0]/taille_gradient)*concentration[3]/**(Diffbetta[0]/Diffbetta[2])*/; //Mo
	
	F[10] = -((ph1->Thermo->m[1]*Y[3])/(ph1->Thermo->m[0]+ph1->Thermo->m[1]-ph1->Thermo->m[1]*Y[4]))*Y[10] \
	+((ph2->Thermo->m[1]*Y[8])/(ph2->Thermo->m[0]+ph2->Thermo->m[1]-ph2->Thermo->m[1]*Y[9]))*(Y[10]-(Diff[0]/taille_gradient)) \
	+(Diff[0]/taille_gradient)*concentration[0]; //C
}

//_______________Jacobien cinetique Fe-Cr-Mo-C____________________
void interface::JK4_1(std::vector <double> &Y,std::vector <double> &fvec, std::vector < std::vector<double> > &Jacobian4D,std::vector <double> &concentration)
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
			KT4_1(n,Ynew,F,concentration);
			Fnew[i][j] = F[i];
			KT4_1(n,Y,F,concentration); 
			Jacobian4D[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
	
}
