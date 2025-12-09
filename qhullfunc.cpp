/*
 *  qhullfunc.cpp
 *  
 *
 *  Created by nat on 12/03/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "phase.hpp"
#include "donnees.h"

//_____L'énergie de Gibbs system FE-C_______________________
double interface::phase1_binaire(double c,double fe,double cr, double mo, double co, double al)
{
	//For QhuickHull	BCC
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	int i;
	
	std::vector <double> y1alfa(1);//bcc_A2: FE
	std::vector <double> y2alfa(2);//bcc_A2: C,VA
	
	
	
	y1alfa[0]=1.;
	y2alfa[0]=c;
	y2alfa[1]=1.-y2alfa[0];
	
	G = ph1->Thermo->functionGibbs(TT2,y1alfa,y2alfa);
	return G;
}


//_____L'énergie de Gibbs system FE-C_______________________
double interface::phase1_binaireFeNi(double c,double fe,double cr, double mo, double co, double al)
{
	//For QhuickHull	BCC
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	int i;
	
	std::vector <double> y1alfa(2);//bcc_A2: FE
	std::vector <double> y2alfa(2);//bcc_A2: C,VA
	
	
	
	y1alfa[1]=c;
	y1alfa[0]=1.-c;
	y2alfa[0]=1.;
	y2alfa[1]=0.;
	
	G = ph1->Thermo->functionGibbs(TT2,y1alfa,y2alfa);
	return G;
}


double interface::phase2_binaire(double c,double fe,double cr, double mo, double co, double al)
{
	//For QuickHull FCC	
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph2->Thermo->T1;
	int i;
	
	std::vector <double> y1betta(1);//fcc_A1: FE
	std::vector <double> y2betta(2);//fcc_A1: C,VA
	
	
	
	y1betta[0]=1.;
	y2betta[0]=c;
	y2betta[1]=1.-y2betta[0];
	G = ph2->Thermo->functionGibbs(TT2,y1betta,y2betta);
	return G;
}


double interface::phase2_binaireFeNi(double c,double fe,double cr, double mo, double co, double al)
{
	//For QuickHull FCC	
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph2->Thermo->T1;
	int i;
	
	std::vector <double> y1betta(2);//bcc_A2: FE
	std::vector <double> y2betta(2);//bcc_A2: C,VA
	
	y1betta[1]=c;
	y1betta[0]=1.-c;
	y2betta[0]=1.;
	y2betta[1]=0.;
	
	G = ph2->Thermo->functionGibbs(TT2,y1betta,y2betta);
	return G;
}

//___________L'énergie de Gibbs Fe-C-Cr________________________
double interface::phase1_ternaire(double c,double fe,double cr, double mo, double co, double al)
{
	//For QuickHull BCC	
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	int i;
	
	std::vector <double> y1alfa(2);//bcc_A2: FE,CR
	std::vector <double> y2alfa(2);//bcc_A2: C,VA
	
	y1alfa[0]=1.0-cr;
	y1alfa[1]=cr;
	y2alfa[0]=c;
	y2alfa[1]=1.0-y2alfa[0];
	
	G = ph1->Thermo->functionGibbs(TT2,y1alfa,y2alfa);
	return G;
}

double interface::phase2_ternaire(double c,double fe,double cr, double mo, double co, double al)
{
	//For QuickHull FCC	
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph2->Thermo->T1;
	int i;
	
	std::vector <double> y1betta(2);//fcc_A1: FE,CR
	std::vector <double> y2betta(2);//fcc_A1: C,VA
	
	y1betta[0]=1.0-cr;
	y1betta[1]=cr;
	y2betta[0]=c;
	y2betta[1]=1.0-y2betta[0];
	
	G = ph2->Thermo->functionGibbs(TT2,y1betta,y2betta);
	
	return G;
}


//__________L'énergie de Gibbs de system Fe-Cr-Mo-C______________________
double interface::phase1_quaternaire(double c,double fe,double cr, double mo, double co, double al)
{
	//For QuickHull BCC
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	int i;
	
	std::vector <double> y1alfa(3);
	std::vector <double> y2alfa(2);
	
	y1alfa[0]=1.0-cr-mo;
	y1alfa[1]=cr;
	y1alfa[2]=mo;
	y2alfa[0]=c;
	y2alfa[1]=1.-y2alfa[0];
	
	G = ph1->Thermo->functionGibbs(TT2,y1alfa,y2alfa);
	return G;
}

double interface::phase2_quaternaire(double c,double fe,double cr, double mo, double co, double al)
{
	//For QuickHull FCC
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph2->Thermo->T1;
	int i;
	
	std::vector <double> y1betta(3);
	std::vector <double> y2betta(2);
	
	y1betta[0]=1.0-cr-mo;
	y1betta[1]=cr;
	y1betta[2]=mo;
	y2betta[0]=c;
	y2betta[1]=1.-y2betta[0];
	
	G = ph2->Thermo->functionGibbs(TT2,y1betta,y2betta);
	
	return G;
}


//___________L'énergie de Gibbs Fe-Cr-Mo-Co-C___________________
double interface::phase1_quinary(double c,double fe,double cr, double mo, double co, double al)
{
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph1->Thermo->T1;
	int i;
	
	std::vector <double> y1alfa(4);
	std::vector <double> y2alfa(2);
	
	y1alfa[0]=1.0-cr-mo-co;
	y1alfa[1]=cr;
	y1alfa[2]=mo;
	y1alfa[3]=co;
	y2alfa[0]=c;
	y2alfa[1]=1.-y2alfa[0];
	
	G = ph1->Thermo->functionGibbs(TT2,y1alfa,y2alfa);
	return G;
}


double interface::phase2_quinary(double c,double fe,double cr, double mo, double co, double al)
{
	double G = 0.0;
	double TT2 = 0.;
	TT2 = ph2->Thermo->T1;
	int i;
	
	std::vector <double> y1betta(4);
	std::vector <double> y2betta(2);
	
	y1betta[0]=1.0-cr-mo-co;
	y1betta[1]=cr;
	y1betta[2]=mo;
	y1betta[3]=co;
	y2betta[0]=c;
	y2betta[1]=1.-y2betta[0];
	
	G = ph2->Thermo->functionGibbs(TT2,y1betta,y2betta);
	
	return G;
} 


//_______La fonction de stabilite Fe-Cr-C_________________________
std::vector < std::vector<double> > interface::JacobianQhull3Dalfa(double c,double fe,double cr, double mo, double co, double al)
{	
	//mu_C - chemical potential "C";
	//mu_Va- chemical potential "VA";
	//mu_Cr- chemical potential "Cr";
	//mu_Fe- chemical potential "Fe";
	
	//Jacobien_bcc = yc*yva*ycr*yfe*[d(mu_C - mu_Va)/d(yc) d(mu_C - mu_Va)/d(ycr) ]
	//	                            [d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr)]
	//mu_VA=0.0;
	//Jacobien_bcc = yc*yva*ycr*yfe*[d(mu_C)/d(yc)         d(mu_C)/d(ycr)         ]
	//	                            [d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr)]	
	
	int n = 2;
	double TT2 = 0.;
	TT2= ph1->Thermo->T1;
	std::vector <double> F(2);
	std::vector <double> Ynew(2);
	std::vector <double> Y(2);
	std::vector <double> pot_chemical1alfa(3);
	for (int i = 0; i < 2; i++)
	{
		F[i]= 0.0;
		Ynew[i]=0.0;
		Y[i]=0.0;
	}
	
	for (int i = 0; i < 3; i++)
	{
		pot_chemical1alfa[i]= 0.0;
	}
	
	std::vector <double> y1alfa(2);
	std::vector <double> y2alfa(2);
	
	for (int i= 0; i < 2; i++)
	{
		y1alfa[i]= 0.0;
		y2alfa[i]= 0.0;
		
	}
	
	y1alfa[0]=1.0-cr;
	y1alfa[1]=cr;
	y2alfa[0]=c;
	y2alfa[1]=1.0-y2alfa[0];
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(y1alfa,y2alfa);
	F[0]=pot_chemical1alfa[0]/(R*TT2); //CARBON;
	F[1]=(pot_chemical1alfa[2]-pot_chemical1alfa[1])/(R*TT2); //CHROME
	
	
	//_________
	//Jacobian
	std::vector < std::vector<double> > Fnew, Jacobian;	
	Fnew.resize(2);
	Jacobian.resize(2);
	
	for (int i = 0; i < Jacobian.size(); i++)
	{
		Fnew[i].resize(2);
		Jacobian[i].resize(2);
	}
	
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
		{
			Fnew[i][j] = 0.0;
			Jacobian[i][j]= 0.0;
		}
	
	Y[0]=y2alfa[0];
	Y[1]=y1alfa[1];
	
	for(int i = 0 ; i < 2 ; i++)
	{
		Ynew[i] = Y[i];
	} 	
	
	for( int i = 0 ; i < 2 ; i++)
		for(int j = 0 ; j < 2 ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			y2alfa[0]=Ynew[0];
			y2alfa[1]=1.0-y2alfa[0];
			y1alfa[1]=Ynew[1];
			y1alfa[0]=1.0-y1alfa[1];
			
			pot_chemical1alfa  =  ph1->Thermo->chemical_potential(y1alfa,y2alfa);
			F[0]=pot_chemical1alfa[0]/(R*TT2); //C
			F[1]=(pot_chemical1alfa[2]-pot_chemical1alfa[1])/(R*TT2); //CR-FE
			
			Fnew[i][j] = F[i];
			y2alfa[0]=Y[0];
			y2alfa[1]=1.0-y2alfa[0];
			y1alfa[1]=Y[1];
			y1alfa[0]=1.0-y1alfa[1];
			
			pot_chemical1alfa  =  ph1->Thermo->chemical_potential(y1alfa,y2alfa);
			F[0]=pot_chemical1alfa[0]/(R*TT2); //C
			F[1]=(pot_chemical1alfa[2]-pot_chemical1alfa[1])/(R*TT2); //CR-FE
			
			Jacobian[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
	Jacobian[0][0]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[0][0]; Jacobian[0][1]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[0][1];
	Jacobian[1][0]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[1][0]; Jacobian[1][1]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[1][1];
	
	return Jacobian;
	
}

std::vector < std::vector<double> > interface::JacobianQhull3Dbetta(double c,double fe,double cr, double mo, double co, double al)
{
	//mu_C - chemical potential "C";
	//mu_Va- chemical potential "VA";
	//mu_Cr- chemical potential "Cr";
	//mu_Fe- chemical potential "Fe";
	
	//Jacobien_fcc = yc*yva*ycr*yfe*[d(mu_C - mu_Va)/d(yc) d(mu_C - mu_Va)/d(ycr) ]
	//	                            [d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr)]
	//mu_VA=0.0;
	//Jacobien_fcc = yc*yva*ycr*yfe*[d(mu_C)/d(yc)         d(mu_C)/d(ycr)         ]
	//	                            [d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr)]	
	int n = 2;
	double TT2 = 0.;
	TT2= ph2->Thermo->T1;
	std::vector <double> F(2);
	std::vector <double> Ynew(2);
	std::vector <double> Y(2);
	std::vector <double> pot_chemical1betta(3);
	for (int i = 0; i < 2; i++)
	{
		F[i]= 0.0;
		Ynew[i]=0.0;
		Y[i]=0.0;
	}
	
	for (int i = 0; i < 3; i++)
	{
		pot_chemical1betta[i]= 0.0;
	}
	
	std::vector <double> y1betta(2);
	std::vector <double> y2betta(2);
	
	for (int i= 0; i < 2; i++)
	{
		y1betta[i]= 0.0;
		y2betta[i]= 0.0;
		
	}
	
	y1betta[0]=1.0-cr;
	y1betta[1]=cr;
	y2betta[0]=c;
	y2betta[1]=1.0-y2betta[0];
	pot_chemical1betta  =  ph2->Thermo->chemical_potential(y1betta,y2betta);
	F[0]=pot_chemical1betta[0]/(R*TT2); //C
	F[1]=(pot_chemical1betta[2]-pot_chemical1betta[1])/(R*TT2); //CR-FE
	
	
	//_________
	//Jacobian
	std::vector < std::vector<double> > Fnew, Jacobian;	
	Fnew.resize(2);
	Jacobian.resize(2);
	
	for (int i = 0; i < Jacobian.size(); i++)
	{
		Fnew[i].resize(2);
		Jacobian[i].resize(2);
	}
	
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
		{
			Fnew[i][j] = 0.0;
			Jacobian[i][j]= 0.0;
		}
	
	Y[0]=y2betta[0];
	Y[1]=y1betta[1];
	
	for(int i = 0 ; i < 2 ; i++)
	{
		Ynew[i] = Y[i];
	} 	
	
	for( int i = 0 ; i < 2 ; i++)
		for(int j = 0 ; j < 2 ; j++) 
		{
			
			Ynew[j] = Y[j] + kdelta;
			y2betta[0]=Ynew[0];
			y2betta[1]=1.0-y2betta[0];
			y1betta[1]=Ynew[1];
			y1betta[0]=1.0-y1betta[1];
			
			pot_chemical1betta  =  ph2->Thermo->chemical_potential(y1betta,y2betta);
			F[0]=pot_chemical1betta[0]/(R*TT2); //C
			F[1]=(pot_chemical1betta[2]-pot_chemical1betta[1])/(R*TT2); //Cr-Fe
			
			Fnew[i][j] = F[i];
			y2betta[0]=Y[0];
			y2betta[1]=1.0-y2betta[0];
			y1betta[1]=Y[1];
			y1betta[0]=1.0-y1betta[1];
			
			pot_chemical1betta  =  ph2->Thermo->chemical_potential(y1betta,y2betta);
			F[0]=pot_chemical1betta[0]/(R*TT2); //C
			F[1]=(pot_chemical1betta[2]-pot_chemical1betta[1])/(R*TT2); //Cr-Fe
			
			Jacobian[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
			
		}
	
	
	Jacobian[0][0]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[0][0]; Jacobian[0][1]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[0][1];
	Jacobian[1][0]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[1][0]; Jacobian[1][1]=c*cr*(1.0-c)*(1.0-cr)*Jacobian[1][1];
	
	return Jacobian;	
	
}


//_________________La fonction de stabilite Fe-Cr-Mo-C__________________________
std::vector < std::vector<double> > interface::JacobianQhull4Dalfa(double c,double fe,double cr, double mo, double co, double al)
{
	//mu_C - chemical potential "C";
	//mu_Va- chemical potential "VA";
	//mu_Cr- chemical potential "Cr";
	//mu_Fe- chemical potential "Fe";
	//mu_Mo- chemical potential "Mo";
	
	//Jacobien_bcc = yc*yva*ycr*yfe*ymo*[d(mu_C - mu_Va)/d(yc) d(mu_C - mu_Va)/d(ycr)  d(mu_C - mu_Va)/d(ymo) ]
	//									[d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr) d(mu_Cr - mu_Fe)/f(ymo)]
	//									[d(mu_Mo-mu_Fe)/d(yc)  d(mu_Mo - mu_Fe)/f(ycr) d(mu_Mo - mu_Fe)/f(ymo)]
	//mu_VA=0.0;
	//Jacobien_bcc = yc*yva*ycr*yfe*ymo*[d(mu_C)/d(yc)		   d(mu_C)/d(ycr)		   d(mu_C)/d(ymo)		  ]
	//									[d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr) d(mu_Cr - mu_Fe)/f(ymo)]
	//									[d(mu_Mo-mu_Fe)/d(yc)  d(mu_Mo - mu_Fe)/f(ycr) d(mu_Mo - mu_Fe)/f(ymo)]
	int n = 3;
	double TT2 = 0.;
	TT2= ph1->Thermo->T1;
	std::vector <double> F(3);
	std::vector <double> Ynew(3);
	std::vector <double> Y(3);
	std::vector <double> pot_chemical1alfa(4);
	for (int i = 0; i < 3; i++)
	{
		F[i]= 0.0;
		Ynew[i]=0.0;
		Y[i]=0.0;
	}
	
	for (int i = 0; i < 4; i++)
	{
		pot_chemical1alfa[i]= 0.0;
	}
	
	std::vector <double> y1alfa(3);
	std::vector <double> y2alfa(2);
	
	for (int i= 0; i < 3; i++)
	{
		y1alfa[i]= 0.0;
		
	}
	
	for (int i=0; i < 2; i++)
	{
		y2alfa[i]=0.0;
	}
	
	y1alfa[0]=1.0-cr-mo;
	y1alfa[1]=cr;
	y1alfa[2]=mo;
	y2alfa[0]=c;
	y2alfa[1]=1.0-y2alfa[0];
	
	
	
	pot_chemical1alfa  =  ph1->Thermo->chemical_potential(y1alfa,y2alfa);
	F[0]=pot_chemical1alfa[0]/(R*TT2); //CARBON;
	F[1]=(pot_chemical1alfa[2]-pot_chemical1alfa[1])/(R*TT2); //Cr-Fe
	F[2]=(pot_chemical1alfa[3]-pot_chemical1alfa[1])/(R*TT2); //Mo-Fe
	
	
	
	//_________
	//Jacobian
	std::vector < std::vector<double> > Fnew, Jacobian;	
	Fnew.resize(3);
	Jacobian.resize(3);
	
	for (int i = 0; i < Jacobian.size(); i++)
	{
		Fnew[i].resize(3);
		Jacobian[i].resize(3);
	}
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			Fnew[i][j] = 0.0;
			Jacobian[i][j]= 0.0;
		}
	
	Y[0]=y2alfa[0];//C
	Y[1]=y1alfa[1];//CR
	Y[2]=y1alfa[2];//MO
	
	for(int i = 0 ; i < 3 ; i++)
	{
		Ynew[i] = Y[i];
	} 	
	
	for( int i = 0 ; i < 3 ; i++)
		for(int j = 0 ; j < 3 ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			y2alfa[0]=Ynew[0]; //C
			y2alfa[1]=1.0-y2alfa[0];//VA
			y1alfa[1]=Ynew[1];//CR
			y1alfa[2]=Ynew[2];//MO
			y1alfa[0]=1.0-y1alfa[1]-y1alfa[2];//FE
			pot_chemical1alfa  =  ph1->Thermo->chemical_potential(y1alfa,y2alfa);
			F[0]=pot_chemical1alfa[0]/(R*TT2); //C
			F[1]=(pot_chemical1alfa[2]-pot_chemical1alfa[1])/(R*TT2); //CR
			F[2]=(pot_chemical1alfa[3]-pot_chemical1alfa[1])/(R*TT2); //MO-FE
			Fnew[i][j] = F[i];
			y2alfa[0]=Y[0]; //C
			y2alfa[1]=1.0-y2alfa[0]; //VA
			y1alfa[1]=Y[1]; //CR
			y1alfa[2]=Y[2]; //MO
			y1alfa[0]=1.0-y1alfa[1]-y1alfa[2]; //FE
			pot_chemical1alfa  =  ph1->Thermo->chemical_potential(y1alfa,y2alfa);
			F[0]=pot_chemical1alfa[0]/(R*TT2); //C
			F[1]=(pot_chemical1alfa[2]-pot_chemical1alfa[1])/(R*TT2); //CR-FE
			F[2]=(pot_chemical1alfa[3]-pot_chemical1alfa[1])/(R*TT2); //MO-FE
			Jacobian[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
	Jacobian[0][0]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[0][0]; 
	Jacobian[0][1]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[0][1];
	Jacobian[0][2]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[0][2];
	Jacobian[1][0]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[1][0];
	Jacobian[1][1]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[1][1];
	Jacobian[1][2]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[1][2];
	Jacobian[2][0]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[2][0];
	Jacobian[2][1]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[2][1];
	Jacobian[2][2]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[2][2];
	
	
	return Jacobian;	
}



std::vector < std::vector<double> > interface::JacobianQhull4Dbetta(double c,double fe,double cr, double mo, double co, double al)
{
	//mu_C - chemical potential "C";
	//mu_Va- chemical potential "VA";
	//mu_Cr- chemical potential "Cr";
	//mu_Fe- chemical potential "Fe";
	//mu_Mo- chemical potential "Mo";
	
	//Jacobien_fcc = yc*yva*ycr*yfe*ymo*[d(mu_C - mu_Va)/d(yc) d(mu_C - mu_Va)/d(ycr)  d(mu_C - mu_Va)/d(ymo) ]
	//									[d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr) d(mu_Cr - mu_Fe)/f(ymo)]
	//									[d(mu_Mo-mu_Fe)/d(yc)  d(mu_Mo - mu_Fe)/f(ycr) d(mu_Mo - mu_Fe)/f(ymo)]
	//mu_VA=0.0;
	//Jacobien_fcc = yc*yva*ycr*yfe*ymo*[d(mu_C)/d(yc)		   d(mu_C)/d(ycr)		   d(mu_C)/d(ymo)		  ]
	//									[d(mu_Cr-mu_Fe)/d(yc)  d(mu_Cr - mu_Fe)/f(ycr) d(mu_Cr - mu_Fe)/f(ymo)]
	//									[d(mu_Mo-mu_Fe)/d(yc)  d(mu_Mo - mu_Fe)/f(ycr) d(mu_Mo - mu_Fe)/f(ymo)]
	int n = 3;
	double TT2 = 0.;
	TT2= ph2->Thermo->T1;
	std::vector <double> F(3);
	std::vector <double> Ynew(3);
	std::vector <double> Y(3);
	std::vector <double> pot_chemical1betta(4);
	for (int i = 0; i < 3; i++)
	{
		F[i]= 0.0;
		Ynew[i]=0.0;
		Y[i]=0.0;
	}
	
	for (int i = 0; i < 4; i++)
	{
		pot_chemical1betta[i]= 0.0;
	}
	
	std::vector <double> y1betta(3);
	std::vector <double> y2betta(2);
	
	for (int i= 0; i < 3; i++)
	{
		y1betta[i]= 0.0;
		
	}
	for (int i=0; i < 2; i++)
	{
		y2betta[i]=0.0;
	}
	y1betta[0]=1.0-cr-mo;
	y1betta[1]=cr;
	y1betta[2]=mo;
	y2betta[0]=c;
	y2betta[1]=1.0-y2betta[0];
	
	
	
	pot_chemical1betta  =  ph2->Thermo->chemical_potential(y1betta,y2betta);
	F[0]=pot_chemical1betta[0]/(R*TT2); //CARBON;
	F[1]=(pot_chemical1betta[2]-pot_chemical1betta[1])/(R*TT2); //Cr-Fe
	F[2]=(pot_chemical1betta[3]-pot_chemical1betta[1])/(R*TT2); //Mo-Fe
	
	//cout<<"2:"<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
	//cout<<"FE:"<<y1betta[0]<<"CR:"<<y1betta[1]<<"MO:"<<y1betta[2]<<endl;
	//cout<<"C:"<<y2betta[0]<<"VA:"<<y2betta[1]<<endl;
	
	//_________
	//Jacobian
	std::vector < std::vector<double> > Fnew, Jacobian;	
	Fnew.resize(3);
	Jacobian.resize(3);
	
	for (int i = 0; i < Jacobian.size(); i++)
	{
		Fnew[i].resize(3);
		Jacobian[i].resize(3);
	}
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			Fnew[i][j] = 0.0;
			Jacobian[i][j]= 0.0;
		}
	
	Y[0]=y2betta[0];//C
	Y[1]=y1betta[1];//CR
	Y[2]=y1betta[2];//MO
	
	for(int i = 0 ; i < 3 ; i++)
	{
		Ynew[i] = Y[i];
	} 	
	
	for( int i = 0 ; i < 3 ; i++)
		for(int j = 0 ; j < 3 ; j++) 
		{
			Ynew[j] = Y[j] + kdelta;
			y2betta[0]=Ynew[0]; //C
			y2betta[1]=1.0-y2betta[0];//VA
			y1betta[1]=Ynew[1];//CR
			y1betta[2]=Ynew[2];//MO
			y1betta[0]=1.0-y1betta[1]-y1betta[2];//FE
			pot_chemical1betta  =  ph2->Thermo->chemical_potential(y1betta,y2betta);
			F[0]=pot_chemical1betta[0]/(R*TT2); //C
			F[1]=(pot_chemical1betta[2]-pot_chemical1betta[1])/(R*TT2); //CR
			F[2]=(pot_chemical1betta[3]-pot_chemical1betta[1])/(R*TT2); //MO-FE
			Fnew[i][j] = F[i];
			y2betta[0]=Y[0]; //C
			y2betta[1]=1.0-y2betta[0]; //VA
			y1betta[1]=Y[1]; //CR
			y1betta[2]=Y[2]; //MO
			y1betta[0]=1.0-y1betta[1]-y1betta[2]; //FE
			pot_chemical1betta  =  ph2->Thermo->chemical_potential(y1betta,y2betta);
			F[0]=pot_chemical1betta[0]/(R*TT2); //C
			F[1]=(pot_chemical1betta[2]-pot_chemical1betta[1])/(R*TT2); //CR-FE
			F[2]=(pot_chemical1betta[3]-pot_chemical1betta[1])/(R*TT2); //MO-FE
			Jacobian[i][j] = (Fnew[i][j] - F[i])/kdelta;
			Ynew[j] = Y[j]; 
		}
	
	Jacobian[0][0]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[0][0]; 
	Jacobian[0][1]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[0][1];
	Jacobian[0][2]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[0][2];
	Jacobian[1][0]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[1][0];
	Jacobian[1][1]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[1][1];
	Jacobian[1][2]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[1][2];
	Jacobian[2][0]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[2][0];
	Jacobian[2][1]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[2][1];
	Jacobian[2][2]=c*cr*mo*(1.0-c)*(1.0-cr-mo)*Jacobian[2][2];
	
	
	return Jacobian;	
}

