/*
 *  thermo.cpp
 *  
 *
 *  Created by nat on 15/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "thermo.h"
#include "donnees.h"
#include "string.h"

//====================================Class Gibbs=================================================
Gibbs::Gibbs() {
	Tc=0;
	Bmag=0;
	p=0;
	A=0;
	n=0;
	fp=0;
	T1=0;
	datafile=0;

}
//Constructeur
Gibbs::Gibbs(const char *nomfichier) : datafile(nomfichier) {
	Tc=0;
	Bmag=0;
	p=0;
	A=0;
	n=0;
	fp=0;
	T1=0;
} 

//Polynom
double Gibbs::polynom(double T1, Array <double,1> &a, Array <double,1> &Gpower) 
{
	int size=a.extent(0);
	double dummy = 0.;
	
	dummy = a(0) * T1 * log(T1);
	
	for (int k=1; k<size; k++) 
		dummy += a(k) * pow(T1,Gpower(k-1));
	return dummy;
}

//To consider quantity of blanks in a line for a text file
int Gibbs::calc_spaces(char *str) 
{
	char *tmp;
	int ans = 0;
	bool sp = false;
	while (*str == ' ')
		++str;
	for (tmp = str; *tmp; ++tmp)
	{
		if (isdigit(*tmp) && sp)
		{
			++ans;
			sp = false;
		}
		if (*tmp == ' ')
			sp = true;
	}
	return ans;
}

//To choose an interval of temperatures and to write down a matrix with the fixed index on temperature GHSER
void Gibbs::rewriteGHSERp(double T1, Array <double,1> &Tr, Array <double,4> &GHSERp, Array <double,3> &GHSERint)
{
	int v = Tr.extent(0); 
	int t , tGHSERp;
	
	for (t=0; t<v-1; t++)
	{
		if (T1==Tr(t))
		{
			tGHSERp = t; GHSERint=GHSERp(tGHSERp,Range::all(),Range::all(),Range::all());
		}
		if (T1==Tr(t+1))
		{
			tGHSERp = t + 1; GHSERint=GHSERp(tGHSERp,Range::all(),Range::all(),Range::all());
		} 
		else if ((T1 - Tr(t))>0 && (T1 - Tr(t+1))<0)
		{
			tGHSERp = t; GHSERint=GHSERp(tGHSERp,Range::all(),Range::all(),Range::all());
		}
	}
}

//To choose an interval of temperatures and to write down a matrix with the fixed index on temperature for interaction parameter Li,j:k
void Gibbs::rewriteL(double T1, Array <double,1> &Tr, Array <double,6> &Aini, Array <double,5> &Aint)
{
	int v = Tr.extent(0);
	int t, tini;
	
	for (t=0; t<v-1; t++)
	{
		if (T1==Tr(t))
		{
			tini = t; Aint=Aini(tini,Range::all(),Range::all(),Range::all(),Range::all(),Range::all());
		}
		if (T1==Tr(t+1))
		{
			tini = t + 1; Aint=Aini(tini,Range::all(),Range::all(),Range::all(),Range::all(),Range::all());
		}
		else if ((T1 - Tr(t))>0 && (T1 - Tr(t+1))<0)
		{
			tini = t; Aint=Aini(tini,Range::all(),Range::all(),Range::all(),Range::all(),Range::all());
		}
	}
}

//To choose an interval of temperatures and to write down a matrix with the fixed index on temperature for interaction parameter Li,j:k,l
void Gibbs::rewriteLtetra(double T1, Array <double,1> &Tr, Array <double,7> &A0, Array <double,6> &A0_int)
{
	int v = Tr.extent(0);
	int t, tini;
	
	for (t=0; t<v-1; t++)
	{
		if (T1==Tr(t))
		{
			tini = t; A0_int=A0(tini,Range::all(),Range::all(),Range::all(),Range::all(),Range::all(),Range::all());
		}
		if (T1==Tr(t+1))
		{
			tini = t + 1; A0_int=A0(tini,Range::all(),Range::all(),Range::all(),Range::all(),Range::all(),Range::all());
		} 
		else if ((T1 - Tr(t))>0 && (T1 - Tr(t+1))<0)
		{
    		tini = t;A0_int=A0(tini,Range::all(),Range::all(),Range::all(),Range::all(),Range::all(),Range::all());
		}
	}
}

//To write down a matrix GHSER with the fixed index on temperature in a polynom
void Gibbs::rewriteGHSERfinal(double T1, Array <double,3> &Aint, Array <double,1> &Gpower,Array <double,2> &Bfinal)
{
	int size = Gpower.extent(0);
	Array <double,1> a(size+1);
	a = 0.;
	double dummy = 0.;
	int i1 = Aint.extent(0);
	int j1 = Aint.extent(1);
	int k1 = Aint.extent(2);
	int i, j, k;
	
	for (i = 0; i < i1; i++)
		for (j = 0; j < j1; j++)
			for (k = 0; k < k1; k++)
			{
				if (Aint(i,j,k)!=0)
				{ 
					a = Aint(i,j,Range::all());
					dummy = polynom(T1,a,Gpower);
					Bfinal(i,j) = dummy;
				}
			} 
}

//To write down a matrix Li,j,:k with the fixed index on temperature in a polynom
void Gibbs::rewriteLfinal(double T1, Array <double,5> &Aint,Array <double,1> &Gpower, Array <double,4> &Bfinal)
{
	int size = Gpower.extent(0);
	Array <double,1> a(size+1);
	a = 0.;
	double dummy = 0.;
	int i1 = Aint.extent(0);
	int j1 = Aint.extent(1);
	int k1 = Aint.extent(2);
	int m1 = Aint.extent(3);
	int n1 = Aint.extent(4);
	int i, j, k, m, n;
	
	for (i = 0; i < i1; i++)
		for (j = 0; j < j1; j++)
			for (k = 0; k < k1; k++)
				for (m = 0; m < m1; m++)
				{
					for (n = 0; n < n1; n++)
					{
						if (Aint(i,j,k,m,n)!=0)
						{
							a = Aint(i,j,k,m,Range::all());
							dummy = polynom(T1,a,Gpower); 
							Bfinal(i,j,k,m) = dummy;
						}
					}
				}
}

//To write down a matrix Li,j:k,l with the fixed index on temperature in a polynom
void Gibbs::rewriteL_tetra_final(double T1, Array <double,6> &A0_int,Array <double,1> &Gpower, Array <double,5> &A0_final)
{
	int size = Gpower.extent(0);
	Array <double,1> a(size+1); 
	a = 0.;
	double dummy = 0.;
	int i1 = A0_int.extent(0);
	int j1 = A0_int.extent(1);
	int k1 = A0_int.extent(2);
	int m1 = A0_int.extent(3);
	int n1 = A0_int.extent(4);
	int l1 = A0_int.extent(5);
	int i, j, k, m, n, l;
	
	for (i = 0; i < i1; i++)
		for (j = 0; j < j1; j++)
			for (k = 0; k < k1; k++)
				for (m = 0; m < m1; m++)
					for (n = 0; n < n1; n++)
						for (l = 0; l < l1; l++) 
						{
							if (A0_int(i,j,k,m,n,l)!=0)
							{
								a = A0_int(i,j,k,m,n,Range::all()); 
								dummy = polynom(T1,a,Gpower);
								A0_final(i,j,k,m,n) = dummy;
							}
						}
}

//Fixes all elements in a matrix Li,j:k, not equal to zero and writes down them in a vector. Indexes of these elements register in a vector, accordingly
void Gibbs::scanmatrixL(Array <double,4> &Afinal,Array <double,1> &A, Array <int,1> &iA, Array <int,1> &jA, Array <int,1> &kA, Array <int,1> &lA)
{
	int i1 = Afinal.extent(0);
	int j1 = Afinal.extent(1);
	int f1 = Afinal.extent(2);
	int l1 = Afinal.extent(3);
	int i, j , f, l;
	int vel = 0;
    
	for(i = 0; i < i1; i++)
		for(l = 0; l < l1; l++)
			for (j = 0; j < j1; j++)
				for (f = 0; f < f1; f++)
					if (Afinal(i,j,f,l)!=0)
					{
						iA(vel)=i;
						jA(vel)=j;       
						kA(vel)=f;
						lA(vel)=l;   
						A(vel++)=Afinal(i,j,f,l); 
					}
	
	A.resizeAndPreserve(vel);
	iA.resizeAndPreserve(vel);
	jA.resizeAndPreserve(vel);
	kA.resizeAndPreserve(vel);
	lA.resizeAndPreserve(vel);
	
}

//Fixes all elements in a matrix Li,j:k,l, not equal to zero and writes down them in a vector. Indexes of these elements register in a vector, accordingly
void Gibbs::scanmatrixLtetra(Array <double,5> &Afinal,Array <double,1> &A,Array <int,1> &iA,Array <int,1> &jA,Array <int,1> &kA,Array <int,1> &lA,Array <int,1> &mA)
{
	int i1 = Afinal.extent(0);
	int j1 = Afinal.extent(1);
	int k1 = Afinal.extent(2);
	int l1 = Afinal.extent(3);
	int m1 = Afinal.extent(4);
	int i, j , k, l, m;
	int vel = 0;
	
	for(i = 0; i < i1; i++)
		for (j = 0; j < j1; j++)
			for (k = 0; k < k1; k++)
				for(l = 0; l < l1; l++)
					for(m = 0; m < m1; m++)
						if (Afinal(i,j,k,l,m)!=0)
						{
							iA(vel)=i;
							jA(vel)=j;       
							kA(vel)=k;
							lA(vel)=l;
							mA(vel)=m;   
							A(vel++)=Afinal(i,j,k,l,m); 
						}
	
	A.resizeAndPreserve(vel);
	A.resizeAndPreserve(vel);
	iA.resizeAndPreserve(vel);
	jA.resizeAndPreserve(vel);
	kA.resizeAndPreserve(vel);
	lA.resizeAndPreserve(vel);
	mA.resizeAndPreserve(vel);
}

//Function for Curie temperature (TC) and the average Bohr magneton number per atom (BMAG) - magnetic model
double Gibbs::functionTC(double T,std::vector <double> &y1, std::vector <double> &y2)
{
	double dummy = 0.;
	int i, j;
	int ig = TCfinal.extent(0); int iex = TC1.extent(0);
	int jg = TCfinal.extent(1);
	
	double dummy1 = 0., dummy2 = 0.;
	
			//reference
	for (i = 0; i < ig; i++)
		for (j = 0; j < jg; j++)
		{
			dummy1 += TCfinal(i,j) * y1[i] * y2[j];
		}	 
			//ex
	//cout << lTC1 <<endl;
	for (i = 0; i < iex; i++)
    {
		if ((jTC1(i) < 2) && (def!=5))
        {
			dummy2 += TC1(i) * y1[iTC1(i)] * y1[jTC1(i)] * y2[kTC1(i)] * pow((y1[jTC1(i)] - y1[iTC1(i)]),lTC1(i));
        }
		else
        {
 			dummy2 += TC1(i) * y1[iTC1(i)] * y1[jTC1(i)] * y2[kTC1(i)] * pow((y1[iTC1(i)] - y1[jTC1(i)]),lTC1(i));
        }
    }
	//cout<<dummy2<<endl;
	//cout<<TCfinal<<endl;
	//cout << dummy2<<endl;
	dummy = dummy1 + dummy2;
	return dummy;
}

double Gibbs::functionBMAG(double T,std::vector <double> &y1, std::vector <double> &y2)
{
	double dummy = 0.;
	int i, j;
	int ig = BMAGfinal.extent(0); int iex = BMAG1.extent(0);
	int jg = BMAGfinal.extent(1);
	double dummy1 = 0., dummy2 = 0.;
	
			//reference
	for (i = 0; i < ig; i++)
		for (j = 0; j < jg; j++)
		{
			dummy1 += BMAGfinal(i,j) * y1[i] * y2[j];
		} 
			//ex
	//cout << lBMAG1 <<endl;

	for (i = 0; i < iex; i++)
	{
		if ((jBMAG1(i) < 2) && (def!=5))
		{
			dummy2 += BMAG1(i) * y1[iBMAG1(i)] * y1[jBMAG1(i)] * y2[kBMAG1(i)] * pow((y1[jBMAG1(i)] - y1[iBMAG1(i)]),lBMAG1(i));  
		}
		else
		{
			dummy2 += BMAG1(i) * y1[iBMAG1(i)] * y1[jBMAG1(i)] * y2[kBMAG1(i)] * pow((y1[iBMAG1(i)] - y1[jBMAG1(i)]),lBMAG1(i));  
		} 
		
	}
	
	//Gm
	dummy = dummy1 + dummy2;
	return dummy;
	
}

double Gibbs::fun_tc(double Tc)
{
	if (Tc>=0.0)
	{
		return Tc;
	}
	
	if (Tc<0.0)
	{
		return -Tc/n;
	}
}

double Gibbs::fun_bmag(double Bmag)
{
	if (Bmag>=0.0)
	{
		return Bmag;
	}	
	
	if (Bmag<0.0)
	{
		return -Bmag/n;
	}
}	

double Gibbs::magnetic(double T,double TC,double BMAG)
{
	A = 0.;
	A = 518./1125. + 11692./15975. * (1/p - 1);
	double tau = 0.0;
	tau = T/fun_tc(TC);
	
	if (tau <= 1.)
	{
		return (1.-1./A*(79./140./p*pow(tau,-1.)+474./497.*(1./p-1.)*(pow(tau,3.)/6.+pow(tau,9.)/135.+pow(tau,15.)/600.)))*log(fun_bmag(BMAG)+1.)*R*T; 
	}
	if (tau > 1.)
	{
		return -1./A*(1./10.*pow(tau,-5.)+1./315.*pow(tau,-15.)+1./1500.*pow(tau,-25.))*log(fun_bmag(BMAG)+1.)*R*T;
	}
	
}

vector <double> Gibbs::magnetic_derivation(double T,double TC,double BMAG,Array <double,1> &dTC,Array <double,1> &dBMAG)
{
	int m = dTC.extent(0);
	int i;
	std::vector <double> MAG(m);
	double tau = 0.0;
	tau = T/fun_tc(TC);	
	
	if (TC < 0)
	{
		if (tau <= 1.)
		{
			for (i = 0; i < m; i++)
			{
				MAG[i]= -1./A*dTC(i)*log(fun_bmag(BMAG)+1.)*R*T*(79./140./p*(pow(T,-1.)/(-1./n)) \
					    +474./497.*(1./p-1.)*(-1./2.*pow(T,3.)/(pow(-TC,4.)*pow((-1./n),3.)) \
					    -1./15.*pow(T,9.)/(pow(-TC,10.)*pow((-1./n),9.)) \
					    -1./40.*pow(T,15.)/(pow(-TC,16.)*pow((-1./n),15)))) \
				        +(1.-1./A*(79./140./p*pow(tau,-1.)+474./497.*(1./p-1.) \
						*(pow(tau,3.)/6.+pow(tau,9.)/135.+pow(tau,15.)/600.)))*R*T*dBMAG(i)*1./(-n)*1./(fun_bmag(BMAG)+1.);
			}
		}
		
		if (tau > 1.)
		{
			for (i = 0; i < m; i++)
			{
				MAG[i]= -1./A*dTC(i)*R*T*log(fun_bmag(BMAG)+1.)*(1./2.*pow(T,-5.)/(pow((-1./n),-5.)*pow(-TC,-4.)) \
					    +1./21.*pow(T,-15.)/(pow((-1./n),-15)*pow(-TC,-14.))+1./60.*pow(T,-25.)/(pow((-1./n),-25)*pow(-TC,-24.))) \
				        -1./A*(1./10.*pow(tau,-5.)+1./315.*pow(tau,-15.)+1./1500.*pow(tau,-25.))*R*T*dBMAG(i)*(-1./n)*1./(fun_bmag(BMAG)+1.); 
			}        
		}
		
	} 
	
	else if (TC > 0)
	{
		if (tau <= 1.)
		{ 
            for (i = 0; i < m; i++)
			{
				MAG[i]= -1./A*dTC(i)*(79./140./p/T+474./497.*(1./p-1.)*(-1./2.*(pow(T,3.)/pow(TC,4.))-1./15.*(pow(T,9.)/pow(TC,10.)) \
					    -1./40.*(pow(T,15.)/pow(TC,16.))))*log(fun_bmag(BMAG)+1.)*R*T+(1.-1./A*(79./140./p*pow(tau,-1.) \
						+474./497.*(1./p-1.)*(pow(tau,3.)/6.+pow(tau,9.)/135.+pow(tau,15.)/600.)))*R*T*dBMAG(i)*1./(fun_bmag(BMAG)+1.);
			}
		}
		if (tau > 1.)
		{
			for (i = 0; i < m; i++)
			{
				MAG[i]= -1./A*dTC(i)*(1./2.*pow(TC,4.)/pow(T,5.)+1./21.*pow(TC,14.)/pow(T,15.)+1./60.*pow(TC,24.)/pow(T,25.))*R*T*log(fun_bmag(BMAG)+1.) \
						-1./A*(1./10.*pow(tau,-5.)+1./315.*pow(tau,-15.)+1./1500.*pow(tau,-25.))*R*T*dBMAG(i)*1./(fun_bmag(BMAG)+1.); 
			}
		}
		
	}
	return MAG;
}

vector < std::vector<double> > Gibbs::magnetic_derivation2(double T,double TC,double BMAG,Array <double,1> &dTC,Array <double,1> &dBMAG,	\
														   Array <double,1> &dTC1,Array <double,1> &dBMAG1,	\
														   Array <double,2> &dTC2,Array <double,2> &dBMAG2, int kk)
{
	//cout<<"TC="<<TC<<endl;
	int m = dTC.extent(0);
	int k = dTC2.extent(0);
	int l = dTC2.extent(1);
	double B,C,D;

	double tau = 0.0;
	tau = T/fun_tc(TC);
	
	int i,j,z;
	std::vector < std::vector<double> >MAG(k);
	for( i = 0 ; i < MAG.size() ; i++)
	{
		MAG[i].resize(l);
	}

	
	if (TC < 0)
	{
		if(tau<=1.)
		{
			//printf("Tc= %lf, T= %lf \n",TC,T);
			//printf("tau= %lf \n",tau);
			//printf("TC<0; TAU<=1 kk =%i \n",kk);
			B = 79./140./p*(pow(T,-1.)/(-1./n))													\
					  +474./497.*(1./p-1.)*( -1./2.*pow(T,3.)/(pow(-TC,4.)*pow((-1./n),3.))		\
											 -1./15.*pow(T,9.)/(pow(-TC,10.)*pow((-1./n),9.))	\
											 -1./40.*pow(T,15.)/(pow(-TC,16.)*pow((-1./n),15)));
			
			//D = 474./497.*(1./p-1.)*( 2.*pow(T,3.)/(pow((-1./n),3.)*pow(-TC,5.))				\
										  +2./3.*pow(T,9.)/(pow((-1./n),9.)*pow(-TC,11.))		\
										  +2./5.*pow(T,15.)/(pow((-1./n),15.)*pow(-TC,17.))		\
										  );
			
			D = 474./497.*(1./p-1.)*( -2.*pow(T,3.)/(pow((-1./n),3.)*pow(-TC,5.))				\
									 -2./3.*pow(T,9.)/(pow((-1./n),9.)*pow(-TC,11.))		\
									 -2./5.*pow(T,15.)/(pow((-1./n),15.)*pow(-TC,17.))		\
									 );
			C = 1.-1./A*(79./140./p*pow(tau,-1.)+474./497.*(1./p-1.)							\
								*(pow(tau,3.)/6.+pow(tau,9.)/135.+pow(tau,15.)/600.));
			
			if(kk==1){
			for(j = 0; j < k; j++)
				for (z = 0; z < l; z++)
			{
				MAG[j][z] = -R*T/A*dTC2(j,z)*log(fun_bmag(BMAG)+1.)*B							\
							-R*T/A*dTC1(z)*dTC(j)*log(fun_bmag(BMAG)+1.)*D						\
							+R*T*dBMAG2(j,z)*(-1./n)/(fun_bmag(BMAG)+1.)*C						\
							-R*T*dBMAG1(z)*dBMAG(j)*(pow((-1./n),2.)/pow((fun_bmag(BMAG)+1.),2))*C \
							-R*T/A*dBMAG1(z)*dTC(j)*B*(-1./n)/(fun_bmag(BMAG)+1.)					\
							-R*T/A*dBMAG(j)*dTC1(z)*B*(-1./n)/(fun_bmag(BMAG)+1.);

			}};
			
			if(kk==2){
				for(j = 0; j < k; j++)
					for (z = 0; z < l; z++)
					{						
						MAG[j][z] = -R*T/A*dTC2(j,z)*log(fun_bmag(BMAG)+1.)*B							\
						-R*T/A*dTC(z)*dTC(j)*log(fun_bmag(BMAG)+1.)*D						\
						+R*T*dBMAG2(j,z)*(-1./n)/(fun_bmag(BMAG)+1.)*C						\
						-R*T*dBMAG(z)*dBMAG(j)*(pow((-1./n),2.)/pow((fun_bmag(BMAG)+1.),2))*C \
						-R*T/A*dBMAG(z)*dTC(j)*B*(-1./n)/(fun_bmag(BMAG)+1.)					\
						-R*T/A*dBMAG(j)*dTC(z)*B*(-1./n)/(fun_bmag(BMAG)+1.);
						
					}};
			
			
		}
		
		if (tau > 1.)
		{	
			//printf("Tc= %lf, T= %lf \n",TC,T);
			//printf("tau= %lf \n",tau);
			//printf("TC<0; TAU>1 kk=%i \n",kk);


			C = 1./10.*pow(tau,-5.)+1./315.*pow(tau,-15.)+1./1500.*pow(tau,-25.);

			B = 1./2.*pow((-1./n),5.)*pow(TC,4.)/pow(T,5.) \
			+1./21.*pow((-1./n),15.)*pow(TC,14.)/pow(T,15.) \
			+1./60.*pow((-1./n),25.)*pow(TC,24.)/pow(T,25.);
			

			//D = 2.*pow((-1./n),5.)*pow(-TC,3.)/pow(T,5.) \
			+2./3.*pow((-1./n),15.)*pow(-TC,13.)/pow(T,15.) \
			+2./5.*pow((-1./n),25.)*pow(-TC,23.)/pow(T,25.);
			
			D = -2.*pow((-1./n),5.)*pow(-TC,3.)/pow(T,5.) \
			-2./3.*pow((-1./n),15.)*pow(-TC,13.)/pow(T,15.) \
			-2./5.*pow((-1./n),25.)*pow(-TC,23.)/pow(T,25.);
			


		  if (kk==2){
			for(j = 0; j < k; j++)
				for (z = 0; z < l; z++)
				{
					MAG[j][z]= -R*T/A*dTC2(j,z)*log(fun_bmag(BMAG)+1.)*B						\
							   -R*T/A*dTC(j)*dTC(z)*log(fun_bmag(BMAG)+1.)*D					\
							   -R*T/A*C*dBMAG2(j,z)*(-1./n)/(fun_bmag(BMAG)+1.)					\
							   +R*T/A*C*dBMAG(j)*dBMAG(z)*(pow((-1./n),2.)/pow((fun_bmag(BMAG)+1.),2)) \
							   -R*T/A*B*((-1./n)/(fun_bmag(BMAG)+1.))*dTC(j)*dBMAG(z)						\
							   -R*T/A*B*((-1./n)/(fun_bmag(BMAG)+1.))*dTC(z)*dBMAG(j);	
				}};

			
			if (kk==1){
				for(j = 0; j < k; j++)
					for (z = 0; z < l; z++)
					{
						MAG[j][z]= -R*T/A*dTC2(j,z)*log(fun_bmag(BMAG)+1.)*B						\
						-R*T/A*dTC1(z)*dTC(j)*log(fun_bmag(BMAG)+1.)*D					\
						-R*T/A*C*dBMAG2(j,z)*(-1./n)/(fun_bmag(BMAG)+1.)					\
						+R*T/A*C*dBMAG1(z)*dBMAG(j)*(pow((-1./n),2.)/pow((fun_bmag(BMAG)+1.),2)) \
						-R*T/A*B*((-1./n)/(fun_bmag(BMAG)+1.))*dTC(j)*dBMAG1(z)						\
						-R*T/A*B*((-1./n)/(fun_bmag(BMAG)+1.))*dTC1(z)*dBMAG(j);	
					}};

			
		}	
	}
	
	else if (TC > 0)
	{
		if (tau <= 1.)
		{ 
			//printf("Here\n");

			//printf("tau= %lf \n",tau);
			//printf("TC>0; TAU<=1 \n");
			B = -1./A*(79./140./p/T+474./497.*(1./p-1.)											\
				*(-1./2.*(pow(T,3.)/pow(TC,4.))													\
				-1./15.*(pow(T,9.)/pow(TC,10.))													\
				-1./40.*(pow(T,15.)/pow(TC,16.))));
			
			C = 1.-1./A*(79./140./p*pow(tau,-1.)												\
						 +474./497.*(1./p-1.)*(pow(tau,3.)/6.							\
													  +pow(tau,9.)/135.+pow(tau,15.)/600.));
			
			D = -1./A*474./497.*(1./p-1.)*(2.*pow(T,3.)/pow(TC,5.)+2./3.*pow(T,9.)/pow(TC,11.)	\
										   +2./5.*pow(T,15.)/pow(TC,17.));
			
			//if(BMAG<0.)BMAG=fun_bmag(BMAG);
			
			//cout<<"B="<<B<<"C="<<C<<"D="<<D<<endl;
			if (kk==2){
			for(j = 0; j < k; j++)
				for (z = 0; z < l; z++)
				{
					MAG[j][z]= R*T*B*dTC2(j,z)*log(BMAG+1.) \
							   +R*T*dTC(j)*dTC(z)*log(BMAG+1.)*D \
							   +R*T*dBMAG2(j,z)*C/(BMAG+1.) \
							   -R*T*dBMAG(j)*dBMAG(z)*C/pow((BMAG+1.),2.) \
							   +R*T*dBMAG(z)*B*dTC(j)/(BMAG+1.)	\
							   +R*T*dBMAG(j)*B*dTC(z)/(BMAG+1.);	
				}};
			/*cout<<"dTC2="<<dTC2<<endl;
			cout<<"log(BMAG+1.)="<<log(BMAG+1.)<<endl;
			cout<<"BMAG+1.="<<(BMAG+1.)<<endl;
			cout<<"dTC="<<dTC<<endl;
			cout<<"dBMAG2="<<dBMAG2<<endl;
			cout<<"1./(BMAG+1.)="<<1./(BMAG+1.)<<endl;
			cout<<"dBMAG="<<dBMAG<<endl;
			cout<<"dTC="<<dTC<<endl;
			cout<<" kk=2   "<<MAG[0][0]<<" "<<MAG[0][1]<<" "<<endl;
			cout<<" kk=2   "<<MAG[1][0]<<" "<<MAG[1][1]<<" "<<endl;
			*/
			if (kk==1){
				for(j = 0; j < k; j++)
					for (z = 0; z < l; z++)
					{
						MAG[j][z]= R*T*B*dTC2(j,z)*log(BMAG+1.) \
						+R*T*dTC1(z)*dTC(j)*log(BMAG+1.)*D \
						+R*T*dBMAG2(j,z)*C/(BMAG+1.) \
						-R*T*dBMAG1(z)*dBMAG(j)*C/pow((BMAG+1.),2.) \
						+R*T*dBMAG1(z)*B*dTC(j)/(BMAG+1.)	\
						+R*T*dBMAG(j)*B*dTC1(z)/(BMAG+1.);	
					}};
			//cout<<" kk=1   "<<MAG[0][0]<<" "<<MAG[0][1]<<" "<<endl;
			//cout<<" kk=1   "<<MAG[1][0]<<" "<<MAG[1][1]<<" "<<endl;
			
			/*cout<<"1:="<<R*T*pow(dTC(0),2.)*log(BMAG+1.)*D<<endl;
			cout<<"D="<<D<<endl;
			cout<<"Tc[0]="<<dTC(0)<<endl;
			cout<<"ln(b+1)="<<log(BMAG+1)<<endl;
			cout<<"R="<<R<<endl;
			cout<<"T="<<T<<endl;
			cout<<"A="<<A<<endl;
			cout<<"p="<<p<<endl;
			cout<<"B="<<B<<endl;
			cout<<"2:"<<2.*R*T*dBMAG(j)*B*dTC(j)/(BMAG+1.)<<endl;
			cout<<"Bmag="<<dBMAG(0)<<endl;*/
			
		}
		
		if (tau > 1.)
		{
			//printf("tau= %lf \n",tau);
			//printf("TC>0; TAU>1 \n");
			B = 1./2.*pow(TC,4.)/pow(T,5.)														\
				+1./21.*pow(TC,14.)/pow(T,15.)													\
				+1./60.*pow(TC,24.)/pow(T,25.);
			
			C = 1./10.*pow(tau,-5.)																\
				+1./315.*pow(tau,-15.)															\
				+1./1500.*pow(tau,-25.);
			
			D = 2.*pow(TC,3.)/pow(T,5.)															\
				+2./3.*pow(TC,13.)/pow(T,15.)													\
			+2./5.*pow(TC,23.)/pow(T,25.);	
			
			//C:C
			if(kk==2){
			for(j = 0; j < k; j++)
				for (z = 0; z < l; z++)
				{
					MAG[j][z] = -R*T/A*dTC2(j,z)*log(BMAG+1.)*B									\
					-R*T/A*dTC(j)*dTC(z)*log(BMAG+1.)*D										\
					-R*T/A*C*dBMAG2(j,z)/(BMAG+1.)												\
					+R*T/A*C*dBMAG(j)*dBMAG(z)/pow((BMAG+1.),2.)									\
					-R*T/A*dTC(j)*dBMAG(z)/(BMAG+1.)*B											\
					-R*T/A*dTC(z)*dBMAG(j)/(BMAG+1.)*B;
				}};
			
			//Fe:C
			if(kk==1){
				for(j = 0; j < k; j++)
					for (z = 0; z < l; z++)
					{
						MAG[j][z] = -R*T/A*dTC2(j,z)*log(BMAG+1.)*B									\
						-R*T/A*dTC1(z)*dTC(j)*log(BMAG+1.)*D										\
						-R*T/A*C*dBMAG2(j,z)/(BMAG+1.)												\
						+R*T/A*C*dBMAG1(z)*dBMAG(j)/pow((BMAG+1.),2.)									\
						-R*T/A*dTC1(z)*dBMAG(j)/(BMAG+1.)*B											\
						-R*T/A*dTC(j)*dBMAG1(z)/(BMAG+1.)*B;
					}};
		}
	}
	
	/*cout<<"B="<<B<<endl;
	cout<<"C="<<C<<endl;
	cout<<"D="<<D<<endl;
	cout<<"log(fun_bmag(BMAG)+1.)="<<log(fun_bmag(BMAG)+1.)<<endl;
	cout<<"(fun_bmag(BMAG)+1.)="<<fun_bmag(BMAG)+1.<<endl;*/
	//cout<<"TC="<<TC<<endl;
	//cout<<"tau="<<tau<<endl;
	//cout<<"    "<<MAG[0][0]<<" "<<MAG[0][1]<<" "<<endl;
	//cout<<"    "<<MAG[1][0]<<" "<<MAG[1][1]<<" "<<endl;
	return MAG;
}



double Gibbs::ylogy(double y)
{
	if (y <= 0.0)
		return 0.0;
	else return y*log(y);
}

double Gibbs::logy(double y)
{
	if(y<=0.0)
		return -1.e33;
	else return log(y);
}

//::::::The modeling quantity is the molar Gibbs energy which is divided into four terms::::::
double Gibbs::functionGibbs(double T,std::vector <double> &y1, std::vector <double> &y2)
{
	double dummyGm = 0.;
	int i, j;
	int ig = GHSERfinal.extent(0);   int i1 = y1.size(); int iex = L1.extent(0); int itetra = Ltetra.extent(0);
	int jg = GHSERfinal.extent(1);   int j1 = y2.size(); int jex = L2.extent(0); int kex = L3.extent(0);
	double dummy1 = 0., dummy2 = 0., dummy2a = 0., dummy2total = 0.,dummy3 = 0., dummy3a = 0., dummy3b = 0., dummy3total = 0.;
	double dummy = 0.,dummy3c = 0.;
	
	//=====reference=====
	for (i = 0; i < ig; i++)
		for (j = 0; j < jg; j++)
		{
			
			dummy1 += GHSERfinal(i,j)* y1[i] * y2[j];
			
		}
	//cout << GHSERfinal<<endl;
	//cout<<"Gref ="<<dummy1<<endl;
	//printf("yfe=%e yni=%e yc=%e yva=%e \n",y1[0],y1[1],y2[0],y2[1]);
	
	//=====mix.ideal=====
	for (i = 0; i < i1; i++)
	{
		
		dummy2 += m[0] * R * T * ylogy(y1[i]);
		//printf("y1[%i]=%e \n",i,y1[i]);
		
	}
	
	for (j = 0; j < j1; j++)
	{
		
		dummy2a += m[1] * R * T * ylogy(y2[j]);
		//printf("y2[%i]=%e \n",i,y2[i]);

	}
	
	dummy2total=dummy2+dummy2a+dummy1;
	//cout<<"(Gmix+ref)="<<dummy2total<<endl;
	//======ex===========
	for (i = 0; i < iex; i++)
	{
		if ((jL1(i) < 2) && (def!=5))
		{
			dummy3 += L1(i) * y1[iL1(i)] * y1[jL1(i)] * y2[kL1(i)] * pow((y1[jL1(i)] - y1[iL1(i)]),lL1(i));
			
		}
		else
		{
			dummy3 += L1(i) * y1[iL1(i)] * y1[jL1(i)] * y2[kL1(i)] * pow((y1[iL1(i)] - y1[jL1(i)]),lL1(i));
			
		}
		
	}
	//cout << L1 <<endl;
	//cout << dummy3<<endl;

	//cout<<"Gex1 ="<<dummy3<<endl;
	
	for (j = 0; j < jex; j++)
	{
		if ((jL2(j) < 2) && (def!=5))
		{           
            dummy3a += L2(j) * y2[iL2(j)] * y2[jL2(j)] * y1[kL2(j)] * pow((y2[jL2(j)] - y2[iL2(j)]),lL2(j));
			
		}
		else
		{
            dummy3a += L2(j) * y2[iL2(j)] * y2[jL2(j)] * y1[kL2(j)] * pow((y2[iL2(j)] - y2[jL2(j)]),lL2(j));
			
		}   
		
	}
	//cout<<"Gex2="<<dummy3a<<endl;
	for (i = 0; i < itetra; i++)
	{
		
		dummy3b += Ltetra(i) * y1[iLtetra(i)] * y1[jLtetra(i)] * y2[kLtetra(i)] * y2[lLtetra(i)];
		
	}
	//cout<<"Gex3 ="<<dummy3b<<endl;
	
	
	for (i = 0; i < kex; i++)
	{
		if (mL3(i)>=2)
		{
            if (mL3(i)==0)
			{
				dummy3c += pow(y1[iL3(i)],2.) * y1[jL3(i)] * y1[kL3(i)] * y2[lL3(i)]*L3(i);
			}
            if (mL3(i)==1)
			{
				dummy3c += y1[iL3(i)] * pow(y1[jL3(i)],2.) * y1[kL3(i)] * y2[lL3(i)]*L3(i);
			}
            if (mL3(i)==2)
			{
				dummy3c += y1[iL3(i)] * y1[jL3(i)] * pow(y1[kL3(i)],2) * y2[lL3(i)]*L3(i);
			}  
			
		}       
		else
		{
			dummy3c += L3(i) * y1[iL3(i)] * y1[jL3(i)] * y1[kL3(i)] * y2[lL3(i)];
		}
	}
	//cout<<"Gex4="<<dummy3c<<endl;
	
	dummy3total = dummy3 + dummy3a + dummy3b + dummy3c;
	//cout<<"Gex="<<dummy3total<<endl;
	//cout<<L1<<" "<<L2<<" "<<L3<<" "<<endl;
	//======mag=========
	Tc = 0.;Bmag = 0.;
	Tc = functionTC(T,y1,y2);
	Bmag = functionBMAG(T,y1,y2);
	
	//printf ("Tc = %lf Bmag = %lf p= %lf n =%lf T= %e\n",Tc,Bmag,p,n,T);
	dummy = magnetic(T,Tc,Bmag);
	//cout<<"Gmag ="<<dummy<<endl;
	
	
	//=======Gm==========
	dummyGm = 1./(m[0] + m[1]*y2[0])*(dummy2total + dummy3total + dummy);
	
	//cout<<"G="<<dummyGm<<endl;
	
	
	return dummyGm;
}

//The modeling is the derivation for molar Gibbs energy 
void Gibbs::function_derivation_Gibbs(double T,std::vector <double> &y1, std::vector <double> &y2)
{  
	int i, j, a;
	int ig =  GHSERfinal.extent(0);    
	int jg =  GHSERfinal.extent(1);
	int iex = L1.extent(0);
	int jex = L2.extent(0);
	int itetra = Ltetra.extent(0);
	int kex = L3.extent(0);
	
	
	Array <double,1> dGmint1(ig), dGmint2(jg), dLint1(ig), dLint2(jg), dLtetra1(ig), dLtetra2(jg),dL3_1(ig),dL3_2(jg);
	
	dGmint1  = 0.;
	dGmint2  = 0.;
	dLint1   = 0.;
	dLint2   = 0.;
	dLtetra1 = 0.;
	dLtetra2 = 0.; 
	dL3_1    = 0.;
	dL3_2    = 0.; 
	dGm1.resize(ig);
	dGm2.resize(jg);
	
	//for (i = 0; i < ig; i++) dGm1[i]=0.0;
	//for (j = 0; j < jg; j++) dGm2[j]=0.0;
	
	
	//::::::::::::::::::Gref et Gmix:::::::::::::::
	//========dG1[i:j]/dy1[i] or dG2[i:j]/dy2[j]
	Array <double,1> dummy(ig), dummy1(ig);
	dummy  = 0.;
	dummy1 = 0.;
	
	for (i = 0; i < ig; i++)
	{
		
        dummy(i) +=m[0] * R * T * (logy(y1[i]) + 1);
		
	}
	
	
	for (i = 0; i < ig; i++)
		for (j = 0; j < jg; j++)
		{
			
			dummy1(i) += GHSERfinal(i,j) * y2[j];
			
		}

	
	for (i = 0; i < ig; i++)
	{
		
        dGmint1(i) = dummy(i) + dummy1(i); 
		
	}
	
	//cout<<"GHSERfinal(i,j)="<<" "<<GHSERfinal<<endl;
	//cout<<"dGmint1(i)="<<" "<<dGmint1<<endl;
	
	Array <double,1> dummy2(jg), dummy3(jg); 
	dummy2 = 0.;
	dummy3 = 0.;
	
	for (j = 0; j < jg; j++)
	{
		
        dummy2(j) += m[1] * R * T * (logy(y2[j]) + 1);
		
	}

	
	for (i = 0; i < ig; i++)
		for (j = 0; j < jg; j++)
		{
			
			dummy3(j) += GHSERfinal(i,j) * y1[i];
			
		}
	
	for (j = 0; j < jg; j++)
	{
		
        dGmint2(j) = dummy2(j) + dummy3(j);
		
	}
	
	//cout<<"dGmint2(i)="<<" "<<dGmint2<<endl;
	
	//:::::::::::::::::::::Lex(fe-cr-mo-c-va)::::::::::::::::::::::::::::::::::::: 
	//ex dL1i,j:k,l/dy1[i,j] or dL2i,j:k,l/dy2[k,l]
	for (a = 0; a < itetra; a++)
	{
		
        dLtetra1(iLtetra(a)) += Ltetra(a) * y2[kLtetra(a)] * y2[lLtetra(a)] * y1[jLtetra(a)];//========dGm/dy_Fe
        dLtetra1(jLtetra(a)) += Ltetra(a) * y2[kLtetra(a)] * y2[lLtetra(a)] * y1[iLtetra(a)];//========dGm/dy_Cr
        dLtetra2(kLtetra(a)) += Ltetra(a) * y2[lLtetra(a)] * y1[iLtetra(a)] * y1[jLtetra(a)];//========dGm/dy_C
        dLtetra2(lLtetra(a)) += Ltetra(a) * y2[kLtetra(a)] * y1[iLtetra(a)] * y1[jLtetra(a)];//========dGm/dy_Va
		
	}
	
	//cout<<"dLtetra1(i)="<<" "<<dLtetra1<<endl;
	//cout<<"dLtetra2(i)="<<" "<<dLtetra2<<endl;
	//:::::::::::::::::::::::::::::::::::
	
	for (a = 0; a < kex; a++)
    {
		if (mL3(a)>=2)
        {
            if (mL3(i)==0)
			{
				dL3_1(iL3(a)) += 2. * L3(a) * y1[iL3(i)] * y1[jL3(a)] * y1[kL3(a)] * y2[lL3(a)]; //Fe
				dL3_1(jL3(a)) += L3(a) * pow(y1[iL3(a)],2.) * y1[kL3(a)] * y2[lL3(a)]; //Cr
				dL3_1(kL3(a)) += L3(a) * pow(y1[iL3(a)],2.) * y1[jL3(a)] * y2[lL3(a)]; //Mo
				dL3_2(lL3(a)) += L3(a) * pow(y1[iL3(a)],2.) * y1[jL3(a)] * y1[kL3(a)]; //C,Va
			}
            if (mL3(i)==1)
			{
				dL3_1(iL3(a)) += L3(a) * pow(y1[jL3(a)],2.) * y1[kL3(a)] * y2[lL3(a)]; //Fe
				dL3_1(jL3(a)) += 2.* L3(a) * y1[iL3(a)] * y1[jL3(a)] * y1[kL3(a)] * y2[lL3(a)]; //Cr
				dL3_1(kL3(a)) += L3(a) * y1[iL3(a)] * pow(y1[jL3(a)],2.) * y2[lL3(a)]; //Mo
				dL3_2(lL3(a)) += L3(a) * y1[iL3(a)] * pow(y1[jL3(a)],2.) * y1[kL3(a)]; //C,Va
			}
            if (mL3(i)==2)
			{
				dL3_1(iL3(a)) += L3(a) * y1[jL3(a)] * pow(y1[kL3(a)],2.) * y2[lL3(a)]; //Fe
				dL3_1(jL3(a)) += L3(a) * y1[iL3(a)] * pow(y1[kL3(a)],2.) * y2[lL3(a)]; //Cr
				dL3_1(kL3(a)) += 2. * L3(a) * y1[iL3(a)] * y1[jL3(a)] * y1[kL3(a)] * y2[lL3(a)]; //Mo
				dL3_2(lL3(a)) += L3(a) * y1[iL3(a)] * y1[jL3(a)] * pow(y1[kL3(a)],2.); //C,Va
			}  
        } 
		else
        {
			dL3_1(iL3(a)) += L3(a) * y1[jL3(a)] * y1[kL3(a)] * y2[lL3(a)]; //Fe
			dL3_1(jL3(a)) += L3(a) * y1[iL3(a)] * y1[kL3(a)] * y2[lL3(a)]; //Cr
			dL3_1(kL3(a)) += L3(a) * y1[iL3(a)] * y1[jL3(a)] * y2[lL3(a)]; //Mo
			dL3_2(lL3(a)) += L3(a) * y1[iL3(a)] * y1[jL3(a)] * y1[kL3(a)]; //C,Va
        }  
		
    } 
	
	/*cout<<"dL3_1 ="<<dL3_1<<endl;
	cout<<"dL3_2="<<dL3_2<<endl;
	cout<<"L3 ="<<L3<<endl;cout<<"iL3 ="<<iL3<<endl;cout<<"jL3 ="<<jL3<<endl;cout<<"kL3 ="<<kL3<<endl;cout<<"lL3="<<lL3<<endl;
	cout<<"L2 ="<<L2<<endl;cout<<"iL2 ="<<iL2<<endl;cout<<"jL2 ="<<jL2<<endl;cout<<"kL2 ="<<kL2<<endl;cout<<"lL2 ="<<lL2<<endl;
	
	//::::::::::::::::::::::Lex(fe/cr/mo/c-va) et (c-va/fe-cr-mo)::::::::::::::::::
	cout<<"L1 ="<<L1<<endl;cout<<"iL1 ="<<iL1<<endl;cout<<"jL1 ="<<jL1<<endl;cout<<"kL1 ="<<kL1<<endl;cout<<"lL1="<<lL1<<endl;
	cout<<"L2 ="<<L2<<endl;cout<<"iL2 ="<<iL2<<endl;cout<<"jL2 ="<<jL2<<endl;cout<<"kL2 ="<<kL2<<endl;cout<<"lL2 ="<<lL2<<endl;*/
	//dL1i,j:k/dy1[i,j]
	for (a = 0; a < iex; a++)
	{
        double check = 0.0;
        if ((jL1(a) < 2) && (def!=5))
		{
			if ((y1[jL1(a)] - y1[iL1(a)])==0.0) {check=1.0;}
			else {check = (y1[jL1(a)] - y1[iL1(a)]);} 
			//===dGm/dy_Fe
			dLint1(iL1(a)) +=L1(a) * y2[kL1(a)] * y1[jL1(a)] * pow((y1[jL1(a)] - y1[iL1(a)]),lL1(a)) \
							-lL1(a) * L1(a) * y2[kL1(a)] * y1[jL1(a)] * y1[iL1(a)] * pow(check,(lL1(a)-1.));
			
			//===dGm/dy_Cr
			dLint1(jL1(a)) +=L1(a) * y2[kL1(a)] * y1[iL1(a)] * pow((y1[jL1(a)] - y1[iL1(a)]),lL1(a)) \
							+lL1(a) * L1(a) * y2[kL1(a)] * y1[iL1(a)] * y1[jL1(a)] * pow(check,(lL1(a)-1.));
			//
			dLint2(kL1(a)) += L1(a) * y1[iL1(a)] * y1[jL1(a)] * pow((y1[jL1(a)] - y1[iL1(a)]),lL1(a));
			
		}
		else
		{
			if ((y1[iL1(a)] - y1[jL1(a)])==0.0) {check=1.0;}
			else {check = (y1[iL1(a)] - y1[jL1(a)]);} 
			//===dGm/dy_Fe
			dLint1(iL1(a)) +=L1(a) * y2[kL1(a)] * y1[jL1(a)] * pow((y1[iL1(a)] - y1[jL1(a)]),lL1(a)) \
							+lL1(a) * L1(a) * y2[kL1(a)] * y1[jL1(a)] * y1[iL1(a)] * pow(check,(lL1(a)-1.));
			
			//===dGm/dy_Cr
			dLint1(jL1(a)) +=L1(a) * y2[kL1(a)] * y1[iL1(a)] * pow((y1[iL1(a)] - y1[jL1(a)]),lL1(a)) \
							-lL1(a) * L1(a) * y2[kL1(a)] * y1[iL1(a)] * y1[jL1(a)] * pow(check,(lL1(a)-1.));
			//
			dLint2(kL1(a)) += L1(a) * y1[iL1(a)] * y1[jL1(a)] * pow((y1[iL1(a)] - y1[jL1(a)]),lL1(a));
			
		}
		
	}
	
    //dL2k,l:i/dy2[k,l]
	for (a = 0; a < jex; a++)
	{
        double check = 0.0;
        if ((y2[jL2(a)] - y2[iL2(a)])==0.0) {check = 1.0;}
        else {check = (y2[jL2(a)] - y2[iL2(a)]);}
		//===dGm/dy_C
        dLint2(iL2(a)) += L2(a) * y2[jL2(a)] * y1[kL2(a)] * pow((y2[jL2(a)] - y2[iL2(a)]),lL2(a)) \
						 -lL2(a) * L2(a) * y2[iL2(a)] * y2[jL2(a)] * y1[kL2(a)] * pow(check,(lL2(a) - 1.));
		//===dGm/dy_Va
        dLint2(jL2(a)) +=L2(a) * y2[iL2(a)] * y1[kL2(a)] * pow((y2[jL2(a)] - y2[iL2(a)]),lL2(a)) \
						+lL2(a) * L2(a) * y2[iL2(a)] * y2[jL2(a)] * y1[kL2(a)] * pow(check,(lL2(a) - 1.));
		//
        dLint1(kL2(a)) += L2(a) * y2[jL2(a)] * y2[iL2(a)] * pow((y2[jL2(a)] - y2[iL2(a)]),lL2(a));
		
	}
	
	//double privet=L1(0)*y2[0]*y1[1]+L1(1)*y2[1]*y1[1]+L1(2)*y2[0]*y1[1]*(y1[0]-y1[1])+L1(2)*y2[0]*y1[1]*y1[0] \
	+L1(3)*y2[1]*y1[1]*(y1[0]-y1[1])+L1(3)*y2[1]*y1[1]*y1[1]+L2(0)*y2[1]*y2[0];
	//cout<<"privet="<<privet<<endl;
	
	//cout<<"dLint1="<<dLint1<<endl;
	//cout<<"dLint2="<<dLint2<<endl;
	//===================================TC===================
	int ig1 = TCfinal.extent(0);  
	int jg1 = TCfinal.extent(1);
	int iex1 = TC1.extent(0); 
	Array <double,1> dTCint1(ig1),dTCint2(jg1),dTCLint1(ig1),dTCLint2(jg1);
	derivTC1.resize(ig1);derivTC2.resize(jg1);
	dTCint1  = 0.;
	dTCLint1 = 0.;
	dTCint2  = 0.;
	dTCLint2 = 0.;
	derivTC1 = 0.;
	derivTC2 = 0.;
	
	for (i = 0; i < ig1; i++)
		for (j = 0; j < jg1; j++)
		{
			dTCint1(i) += TCfinal(i,j) * y2[j];
		}
	
	
	for (i = 0; i < ig1; i++)
		for (j = 0; j < jg1; j++)
		{
			dTCint2(j) += TCfinal(i,j) * y1[i];
		} 
	//cout<<"dTCint2="<<dTCint2<<endl;
	
	
	//TCL1
	for (a = 0; a < iex1; a++)
	{
        double check = 0.0; 
		if ((jTC1(a) < 2) && (def!=5))
		{

			if ((y1[jTC1(a)] - y1[iTC1(a)])==0.0) {check = 1.0;}
			else {check = (y1[jTC1(a)] - y1[iTC1(a)]);}
			
			dTCLint1(iTC1(a)) +=TC1(a) * y2[kTC1(a)] * y1[jTC1(a)] * pow((y1[jTC1(a)] - y1[iTC1(a)]),lTC1(a)) \
								- lTC1(a) * TC1(a) * y2[kTC1(a)] * y1[jTC1(a)] * y1[iTC1(a)] * pow(check,(lTC1(a)-1.));
			
			dTCLint1(jTC1(a)) +=TC1(a) * y2[kTC1(a)] * y1[iTC1(a)] * pow((y1[jTC1(a)] - y1[iTC1(a)]),lTC1(a)) \
								+ lTC1(a) * TC1(a) * y2[kTC1(a)] * y1[iTC1(a)] * y1[jTC1(a)] * pow(check,(lTC1(a)-1.));
			
			dTCLint2(kTC1(a)) += TC1(a) * y1[iTC1(a)] * y1[jTC1(a)] * pow((y1[jTC1(a)] - y1[iTC1(a)]),lTC1(a));
			
		}
		else
		{

			if ((y1[iTC1(a)] - y1[jTC1(a)])==0.0) {check = 1.0;}
			else {check = (y1[iTC1(a)] - y1[jTC1(a)]);}
			
			dTCLint1(iTC1(a)) +=TC1(a) * y2[kTC1(a)] * y1[jTC1(a)] * pow((y1[iTC1(a)] - y1[jTC1(a)]),lTC1(a)) \
								+ lTC1(a) * TC1(a) * y2[kTC1(a)] * y1[jTC1(a)] * y1[iTC1(a)] * pow(check,(lTC1(a)-1.));
			
			dTCLint1(jTC1(a)) +=TC1(a) * y2[kTC1(a)] * y1[iTC1(a)] * pow((y1[iTC1(a)] - y1[jTC1(a)]),lTC1(a)) \
								- lTC1(a) * TC1(a) * y2[kTC1(a)] * y1[iTC1(a)] * y1[jTC1(a)] * pow(check,(lTC1(a)-1.));
			
			dTCLint2(kTC1(a)) += TC1(a) * y1[iTC1(a)] * y1[jTC1(a)] * pow((y1[iTC1(a)] - y1[jTC1(a)]),lTC1(a));
			
		}
	}
	
	derivTC1 = dTCint1 + dTCLint1;
	derivTC2 = dTCint2 + dTCLint2;
	//cout<<"derivTC1="<<derivTC1<<endl;
	//cout<<"derivTC2="<<derivTC2<<endl;
	
	//====================================BMAG=========================
	int ig2 = BMAGfinal.extent(0); 
	int jg2 = BMAGfinal.extent(1); 
	int iex2 = BMAG1.extent(0);
	Array <double,1> dBMAGint1(ig2), dBMAGLint1(ig2), dBMAGint2(jg2),dBMAGLint2(jg2);
	derivBMAG1.resize(ig2); derivBMAG2.resize(jg2);
	dBMAGint1  = 0.;
	dBMAGLint1 = 0.;
	dBMAGint2  = 0.;
	dBMAGLint2 = 0.;
	derivBMAG1 = 0.;
	derivBMAG2 = 0.;
	
	
	for (i = 0; i < ig2; i++)
		for (j = 0; j < jg2; j++)
		{
			dBMAGint1(i) += BMAGfinal(i,j) * y2[j];
		} 
	
	for (i = 0; i < ig2; i++)
		for (j = 0; j < jg2; j++)
		{
			dBMAGint2(j) += BMAGfinal(i,j) * y1[i];
		} 
	
	//BMAGL1
	for (a = 0; a < iex2; a++)
	{
		double check = 0.0;
		if ((jBMAG1(a) < 2) && (def!=5))
		{

			if ((y1[jBMAG1(a)] - y1[iBMAG1(a)])==0.0) {check = 1.0;}
			else {check = (y1[jBMAG1(a)] - y1[iBMAG1(a)]);}
			
			dBMAGLint1(iBMAG1(a)) +=BMAG1(a) * y2[kBMAG1(a)] * y1[jBMAG1(a)] * pow((y1[jBMAG1(a)] - y1[iBMAG1(a)]),lBMAG1(a)) \
									-lBMAG1(a) * BMAG1(a) * y2[kBMAG1(a)] * y1[jBMAG1(a)] * y1[iBMAG1(a)] * pow(check,(lBMAG1(a)-1.));
			
			dBMAGLint1(jBMAG1(a)) +=BMAG1(a) * y2[kBMAG1(a)] * y1[iBMAG1(a)] * pow((y1[jBMAG1(a)] - y1[iBMAG1(a)]),lBMAG1(a)) \
									+lBMAG1(a) * BMAG1(a) * y2[kBMAG1(a)] * y1[iBMAG1(a)] * y1[jBMAG1(a)] * pow(check,(lBMAG1(a)-1.));
			
			dBMAGLint2(kBMAG1(a)) += BMAG1(a) * y1[iBMAG1(a)] * y1[jBMAG1(a)] * pow((y1[jBMAG1(a)] - y1[iBMAG1(a)]),lBMAG1(a));
			
		}  
        else
		{

			if ((y1[iBMAG1(a)] - y1[jBMAG1(a)])==0.0) {check = 1.0;}
			else {check = (y1[iBMAG1(a)] - y1[jBMAG1(a)]);}
			
			dBMAGLint1(iBMAG1(a)) +=BMAG1(a) * y2[kBMAG1(a)] * y1[jBMAG1(a)] * pow((y1[iBMAG1(a)] - y1[jBMAG1(a)]),lBMAG1(a)) \
									+lBMAG1(a) * BMAG1(a) * y2[kBMAG1(a)] * y1[jBMAG1(a)] * y1[iBMAG1(a)] * pow(check,(lBMAG1(a)-1.));
			
			dBMAGLint1(jBMAG1(a)) +=BMAG1(a) * y2[kBMAG1(a)] * y1[iBMAG1(a)] * pow((y1[iBMAG1(a)] - y1[jBMAG1(a)]),lBMAG1(a)) \
									-lBMAG1(a) * BMAG1(a) * y2[kBMAG1(a)] * y1[iBMAG1(a)] * y1[jBMAG1(a)] * pow(check,(lBMAG1(a)-1.));
			
			dBMAGLint2(kBMAG1(a)) += BMAG1(a) * y1[iBMAG1(a)] * y1[jBMAG1(a)] * pow((y1[iBMAG1(a)] - y1[jBMAG1(a)]),lBMAG1(a));
			
		}
	}
	
	derivBMAG1 = dBMAGint1 + dBMAGLint1;
	derivBMAG2 = dBMAGint2 + dBMAGLint2;
 
	/*cout<<"dBMAGLint1="<<dBMAGLint1<<endl;
	cout<<"dBMAGLint2="<<dBMAGLint2<<endl;
	cout<<"derivBMAG1="<<derivBMAG1<<endl;
	cout<<"derivBMAG2="<<derivBMAG2<<endl;
	cout<<"TC="<<Tc<<endl;
	cout<<"n="<<n<<"p="<<p<<endl;*/
	
	std::vector <double> MAG1(ig1);
	std::vector <double> MAG2(jg1);
	MAG1 = magnetic_derivation(T,Tc,Bmag,derivTC1,derivBMAG1);
	MAG2 = magnetic_derivation(T,Tc,Bmag,derivTC2,derivBMAG2);
	

	//for (i = 0; i < ig1; i++)
	//printf ("MAG1 =%g \n",MAG1[i]);
	
	//for (i = 0; i < jg1; i++)
	//printf ("MAG2 =%lf \n",MAG2[i]); 
	//===Total=
	
	for (i = 0; i < ig; i++)
    {
        
		dGm1[i] = dGmint1(i) + dLint1(i)+ dLtetra1(i)+MAG1[i]+dL3_1(i);
		
    }
	
	//if(T==702.){
	//for (int i = 0; i < ig; i++){
	//printf("dGm1[%i]=%g ",i,dGm1[i]);}
    //printf ("\n");
	
	for (i = 0; i < jg; i++)
    {
		
		dGm2[i] = dGmint2(i) + dLint2(i) + dLtetra2(i)+MAG2[i]+dL3_2(i);
		
    }
	
	//if(T==702.){
	//for (int i = 0; i < ig; i++)
	//printf("dGm2[%i]=%lf \n",i,dGm2[i]);
	
}







//::::::::::::::::Jacobian::::::::::::::::::::::::
void Gibbs::function_derivation_Gibbs2(double T,std::vector <double> &y1, std::vector <double> &y2)
{
	int i, j, a;
	int ig     = GHSERfinal.extent(0);    
	int jg     = GHSERfinal.extent(1);    
	int iex    = L1.extent(0);
	int jex    = L2.extent(0);
	int kex    = L3.extent(0);
	
	//:::::::::::::::::::::::::::
	Array <double,2> Gmix1(ig,ig), Gref1(ig,jg);
	Gmix1  = 0.;
	Gref1  = 0.;
	for (i = 0; i < ig; i++)
	{
		
        Gmix1(i,i) +=m[0] * R * T * 1./y1[i]/* + R * T * m[0]*/;
        if (y1[i]==0.0) {/*printf ("y1=0.0, error nan\n");*/y1[i]=1.e-9;}
		
	}
	
	//cout<<"Gmix1="<<Gmix1<<endl;
	
	for (i = 0; i < ig; i++)
		for (j = 0; j < jg; j++)
		{
			
			Gref1(i,j) += GHSERfinal(i,j);
			
		}
	
	//cout<<"Gref1="<<Gref1<<endl;
	//::::::::::::::::::::::::::::
	Array <double,2> Gmix2(jg,jg), Gref2(jg,ig); 
	Gmix2 = 0.;
	Gref2 = 0.;
	
	for (j = 0; j < jg; j++)
	{
		
        Gmix2(j,j) += m[1] * R * T * 1./y2[j]/* + R * T* m[1]*/;
        if (y2[j]==0.0) {/*printf ("y2=0.0, error nan\n");*/y2[j]=1.e-9;}
		
	}
	//cout<<"Gmix2="<<Gmix2<<endl;
	
	
	for (i = 0; i < ig; i++)
		for (j = 0; j < jg; j++)
		{
			
			Gref2(j,i) += GHSERfinal(i,j);
			
		}
	//cout<<"Gref2="<<Gref2<<endl;
	
	//:::::::::::::::::::::::::::::::::::::
	Array <double,2> dL3_1_2(ig,ig),dL3_2_2(ig,jg),dL3_3_2(jg,ig);
	dL3_1_2 = 0.;
	dL3_2_2 = 0.;
	dL3_3_2 = 0.; 
	
	for (a = 0; a < kex; a++)
    {
		
		dL3_1_2(iL3(a),jL3(a)) += L3(a)  * y1[kL3(a)] * y2[lL3(a)]; //Fe_CR
		dL3_1_2(iL3(a),kL3(a)) += L3(a)  * y1[jL3(a)] * y2[lL3(a)]; //Fe_MO
		dL3_2_2(iL3(a),lL3(a)) += L3(a)  * y1[jL3(a)] * y1[kL3(a)]; //Fe_C/VA    
		
		dL3_1_2(jL3(a),iL3(a)) += L3(a)  * y1[kL3(a)] * y2[lL3(a)]; //Cr_Fe
		dL3_1_2(jL3(a),kL3(a)) += L3(a)  * y1[iL3(a)] * y2[lL3(a)]; //Cr_Mo
		dL3_2_2(jL3(a),lL3(a)) += L3(a)  * y1[iL3(a)] * y1[kL3(a)]; //Cr_C/Va
		
		dL3_1_2(kL3(a),iL3(a)) += L3(a)  * y1[jL3(a)] * y2[lL3(a)]; //Mo_Fe
		dL3_1_2(kL3(a),jL3(a)) += L3(a)  * y1[iL3(a)] * y2[lL3(a)]; //Mo_Cr
		dL3_2_2(kL3(a),lL3(a)) += L3(a)  * y1[iL3(a)] * y1[jL3(a)]; //Mo_C/VA
		
		dL3_3_2(lL3(a),iL3(a)) += L3(a)  * y1[jL3(a)] * y1[kL3(a)]; //C,Va_Fe
		dL3_3_2(lL3(a),jL3(a)) += L3(a)  * y1[iL3(a)] * y1[kL3(a)]; //C,Va_Cr
		dL3_3_2(lL3(a),kL3(a)) += L3(a)  * y1[iL3(a)] * y1[jL3(a)]; //C,Va_Mo
		
    }
	
	//:::::::::::::::::::::::::::::::::::::
	Array <double,2> dLint1_2(ig,ig), dLint1_3(ig,jg), dLint2_4(jg,jg), dLint2_5(jg,ig);
	dLint1_2 = 0.;
	dLint1_3 = 0.;
	dLint2_4 = 0.;
	dLint2_5 = 0.;
	
	//cout<<"dL3_1_2="<<dL3_1_2<<endl;
	//cout<<"dL3_2_2="<<dL3_2_2<<endl;

	//cout<<"L1 ="<<L1<<endl;cout<<"iL1 ="<<iL1<<endl;cout<<"jL1 ="<<jL1<<endl;cout<<"kL1 ="<<kL1<<endl;cout<<"lL1="<<lL1<<endl; 
	for (a = 0; a < iex; a++)
	{
		double check = 0.0;
		if (jL1(a)<2 && (def!=5))
		{
			if ((y1[jL1(a)]-y1[iL1(a)])==0.0) {check = 1.0;}
			else {check = (y1[jL1(a)]-y1[iL1(a)]);}
			
			dLint1_2(iL1(a),iL1(a)) += -lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										+lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2)) \
										-lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1));//Fe_Fe    
			
			
			dLint1_2(jL1(a),jL1(a)) += lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										+lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2)) \
										+lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1));//Cr_Cr   
			
			dLint1_2(iL1(a),jL1(a)) += L1(a)*y2[kL1(a)]*pow((y1[jL1(a)]-y1[iL1(a)]),lL1(a)) \
										-lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										+lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										-lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2));//Fe_Cr/Mo
			
			dLint1_2(jL1(a),iL1(a)) += L1(a)*y2[kL1(a)]*pow((y1[jL1(a)]-y1[iL1(a)]),lL1(a)) \
										-lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										+lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										-lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2));//Cr/Mo_Fe   
			
			dLint1_3(iL1(a),kL1(a)) += L1(a)*y1[jL1(a)]*pow((y1[jL1(a)]-y1[iL1(a)]),lL1(a)) \
									   -lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1)); //Fe/Cr:C/Va
			
			dLint1_3(jL1(a),kL1(a)) += L1(a)*y1[iL1(a)]*pow((y1[jL1(a)]-y1[iL1(a)]),lL1(a)) \
										+lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1)); //Fe/Cr:C/Va
			
			dLint2_5(kL1(a),iL1(a)) += L1(a)*y1[jL1(a)]*pow((y1[jL1(a)]-y1[iL1(a)]),lL1(a)) \
										-lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1));//c/va:fe/cr      
			
			dLint2_5(kL1(a),jL1(a)) += L1(a)*y1[iL1(a)]*pow((y1[jL1(a)]-y1[iL1(a)]),lL1(a)) \
										+lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1));//c/va:cr/mo
			
		}
		else
		{
			if ((y1[iL1(a)]-y1[jL1(a)])==0.0) {check = 1.0;}
			else {check = (y1[iL1(a)]-y1[jL1(a)]);}
			dLint1_2(iL1(a),iL1(a)) += lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										+lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2)) \
										+lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1));//Fe_Fe     
			
			
			dLint1_2(jL1(a),jL1(a)) += -lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										+lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2)) \
										-lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1));//Cr_Cr 
			
			
			dLint1_2(iL1(a),jL1(a)) += L1(a)*y2[kL1(a)]*pow((y1[iL1(a)]-y1[jL1(a)]),lL1(a)) \
										+lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										-lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										-lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2));//Fe_Cr/Mo
			
			
			dLint1_2(jL1(a),iL1(a)) += L1(a)*y2[kL1(a)]*pow((y1[iL1(a)]-y1[jL1(a)]),lL1(a)) \
										+lL1(a)*L1(a)*y1[iL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										-lL1(a)*L1(a)*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-1)) \
										-lL1(a)*(lL1(a)-1)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*y2[kL1(a)]*pow(check,(lL1(a)-2));//Cr/Mo_Fe   
			
			
			dLint1_3(iL1(a),kL1(a)) += L1(a)*y1[jL1(a)]*pow((y1[iL1(a)]-y1[jL1(a)]),lL1(a)) \
										+lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1)); //Fe/Cr:C/Va
			
			dLint1_3(jL1(a),kL1(a)) += L1(a)*y1[iL1(a)]*pow((y1[iL1(a)]-y1[jL1(a)]),lL1(a)) \
										-lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1)); //Fe/Cr:C/Va
			
			dLint2_5(kL1(a),iL1(a)) += L1(a)*y1[jL1(a)]*pow((y1[iL1(a)]-y1[jL1(a)]),lL1(a)) \
										+lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1));
			
			dLint2_5(kL1(a),jL1(a)) += L1(a)*y1[iL1(a)]*pow((y1[iL1(a)]-y1[jL1(a)]),lL1(a)) \
										-lL1(a)*L1(a)*y1[iL1(a)]*y1[jL1(a)]*pow(check,(lL1(a)-1));
			
		}
		
	}
	//cout<<"dLint1_2="<<dLint1_2<<endl;
	//cout<<"dLint1_3="<<dLint1_3<<endl;
	//cout<<"dLint2_5="<<dLint2_5<<endl;

	//cout<<"L2 ="<<L2<<endl;cout<<"iL2 ="<<iL2<<endl;cout<<"jL2 ="<<jL2<<endl;cout<<"kL2 ="<<kL2<<endl;cout<<"lL2="<<lL2<<endl; 
	for (a = 0; a < jex; a++)
	{
		double check = 0.0; 
		if ((y2[jL2(a)]-y2[iL2(a)])==0.0) {check = 1.0;}
		else {check = (y2[jL2(a)]-y2[iL2(a)]);}
		
		dLint1_3(kL2(a),iL2(a)) += L2(a)*y2[jL2(a)]*pow((y2[jL2(a)]-y2[iL2(a)]),lL2(a)) \
									-lL2(a)*L2(a)*y2[iL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1));//fe,cr,mo:c
		
		dLint1_3(kL2(a),jL2(a)) += L2(a)*y2[iL2(a)]*pow((y2[jL2(a)]-y2[iL2(a)]),lL2(a)) \
									+lL2(a)*L2(a)*y2[iL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1));//fe,cr,mo:va
		
		dLint2_4(iL2(a),iL2(a)) += -lL2(a)*L2(a)*y1[kL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1)) \
									+lL2(a)*(lL2(a)-1)*L2(a)*y1[kL2(a)]*y2[iL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-2)) \
									-lL2(a)*L2(a)*y1[kL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1));//c:c 
		
		dLint2_4(jL2(a),jL2(a)) +=  lL2(a)*L2(a)*y1[kL2(a)]*y2[iL2(a)]*pow(check,(lL2(a)-1)) \
									+lL2(a)*(lL2(a)-1)*L2(a)*y1[kL2(a)]*y2[iL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-2)) \
									+lL2(a)*L2(a)*y1[kL2(a)]*y2[iL2(a)]*pow(check,(lL2(a)-1));//va:va
		
		dLint2_4(iL2(a),jL2(a)) +=  L2(a)*y1[kL2(a)]*pow((y2[jL2(a)]-y2[iL2(a)]),lL2(a)) \
									-lL2(a)*L2(a)*y1[kL2(a)]*y2[iL2(a)]*pow(check,(lL2(a)-1)) \
									+lL2(a)*L2(a)*y1[kL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1)) \
									-lL2(a)*(lL2(a)-1)*L2(a)*y1[kL2(a)]*y2[jL2(a)]*y2[iL2(a)]*pow(check,(lL2(a)-2));//c:va
		
		dLint2_4(jL2(a),iL2(a)) +=  L2(a)*y1[kL2(a)]*pow((y2[jL2(a)]-y2[iL2(a)]),lL2(a)) \
									+lL2(a)*L2(a)*y1[kL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1)) \
									-lL2(a)*L2(a)*y1[kL2(a)]*y2[iL2(a)]*pow(check,(lL2(a)-1)) \
									-lL2(a)*(lL2(a)-2)*L2(a)*y1[kL2(a)]*y2[jL2(a)]*y2[iL2(a)]*pow(check,(lL2(a)-2));//va:c  
		
		dLint2_5(iL2(a),kL2(a)) += L2(a)*y2[jL2(a)]*pow((y2[jL2(a)]-y2[iL2(a)]),lL2(a)) \
									-lL2(a)*L2(a)*y2[iL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1));
		
		dLint2_5(jL2(a),kL2(a)) += L2(a)*y2[iL2(a)]*pow((y2[jL2(a)]-y2[iL2(a)]),lL2(a)) \
									+lL2(a)*L2(a)*y2[iL2(a)]*y2[jL2(a)]*pow(check,(lL2(a)-1));
	}
	
	//cout<<"dLint1_3="<<dLint1_3<<endl;
	//cout<<"dLint2_4="<<dLint2_4<<endl;
	//cout<<"dLint2_5="<<dLint2_5<<endl;

	
	//:::::::::::::::::::::::::::::::::::::::::::::::
	//===================================TC(Begin)===================
	int ig1  = TCfinal.extent(0); 
	//cout<<"ig1="<<ig1<<endl;
	int jg1  = TCfinal.extent(1);
	//cout<<"jg1="<<jg1<<endl;
	int iex1 = TC1.extent(0); 
	
	Array <double,2> dTcref1(ig1,jg1), dTcref2(jg1,ig1), dTc1(ig1,ig1), dTc2(ig1,jg1), dTc3(jg1,ig1);
	dTcref1 = 0.;
	dTcref2 = 0.;
	dTc1    = 0.;
	dTc2    = 0.;
	dTc3    = 0.;
	
	//cout<<TCfinal<<endl;
	
	for (i = 0; i < ig1; i++)
		for (j = 0; j < jg1; j++)
		{
			//if (p==0.28)
			//{
			//	dTcref1(i,j) +=  -TCfinal(i,j)/n;
			//}
			//else
			//{
				dTcref1(i,j) +=  TCfinal(i,j);
			//}	
		}
	
	for (i = 0; i < ig1; i++)
		for (j = 0; j < jg1; j++)
		{
			//if (p==0.28)
			//{
			//	dTcref2(j,i) += -TCfinal(i,j)/n;
			//}
			//else
			//{
				dTcref2(j,i) += TCfinal(i,j);
			//}	
		}
	
	//cout<<"dTcref1="<<dTcref1<<endl;
	//cout<<"dTcref2="<<dTcref2<<endl;
	//cout<<"TC1 ="<<TC1<<endl;cout<<"iTC1 ="<<iTC1<<endl;cout<<"jTC1 ="<<jTC1<<endl;cout<<"kTC1 ="<<kTC1<<endl;cout<<"lTC1="<<lTC1<<endl;
	
	
	//::::::::::::::::::::::::::::::::::::::
	for (a = 0; a < iex1; a++)
    {
		double check = 0.0; 
		/*if (p==0.28)
		{
			if ((jTC1(a)<2) && (def!=5))
			{
				if ((y1[jTC1(a)]-y1[iTC1(a)])==0.0) {check = 1.0;}
				else {check = (y1[jTC1(a)]-y1[iTC1(a)]);}
				
				dTc1(iTC1(a),iTC1(a)) += -lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										-lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Fe_Fe 
				
				dTc1(jTC1(a),jTC1(a)) += lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										+lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Cr_Cr   
				
				dTc1(iTC1(a),jTC1(a)) += (-TC1(a)/n)*y2[kTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Fe_Cr/Mo
				
				dTc1(jTC1(a),iTC1(a)) += (-TC1(a)/n)*y2[kTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Cr/Mo_Fe
				
				dTc2(iTC1(a),kTC1(a)) += (-TC1(a)/n)*y1[jTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc2(jTC1(a),kTC1(a)) += (-TC1(a)/n)*y1[iTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										+lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc3(kTC1(a),iTC1(a)) += (-TC1(a)/n)*y1[jTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1));//c/va:fe/cr      
				
				dTc3(kTC1(a),jTC1(a)) += (-TC1(a)/n)*y1[iTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										+lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1));//c/va:cr/mo
				
			} 
			
			else
			{
				if ((y1[iTC1(a)]-y1[jTC1(a)])==0.0) {check = 1.0;}
				else {check = (y1[iTC1(a)]-y1[jTC1(a)]);}
				
				dTc1(iTC1(a),iTC1(a)) += lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										+lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Fe_Fe     
				
				dTc1(jTC1(a),jTC1(a)) += -lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										-lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Cr_Cr 
				
				dTc1(iTC1(a),jTC1(a)) += (-TC1(a)/n)*y2[kTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Fe_Cr/Mo
				
				dTc1(jTC1(a),iTC1(a)) += (-TC1(a)/n)*y2[kTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(-TC1(a)/n)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Cr/Mo_Fe
				
				dTc2(iTC1(a),kTC1(a)) += (-TC1(a)/n)*y1[jTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc2(jTC1(a),kTC1(a)) += (-TC1(a)/n)*y1[iTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										-lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc3(kTC1(a),iTC1(a)) += (-TC1(a)/n)*y1[jTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //c/va:fe/cr
				
				dTc3(kTC1(a),jTC1(a)) += (-TC1(a)/n)*y1[iTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										-lTC1(a)*(-TC1(a)/n)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //c/va:cr/mo
			}
			
		}*/
		
		//else
		//{	 
			if ((jTC1(a)<2) && (def!=5))
			{
				if ((y1[jTC1(a)]-y1[iTC1(a)])==0.0) {check = 1.0;}
				else {check = (y1[jTC1(a)]-y1[iTC1(a)]);}
				
				dTc1(iTC1(a),iTC1(a)) += -lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										-lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Fe_Fe 
				
				dTc1(jTC1(a),jTC1(a)) += lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										+lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Cr_Cr   
				
				dTc1(iTC1(a),jTC1(a)) += TC1(a)*y2[kTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Fe_Cr/Mo
				
				dTc1(jTC1(a),iTC1(a)) += TC1(a)*y2[kTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Cr/Mo_Fe
				
				dTc2(iTC1(a),kTC1(a)) += TC1(a)*y1[jTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc2(jTC1(a),kTC1(a)) += TC1(a)*y1[iTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										+lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc3(kTC1(a),iTC1(a)) += TC1(a)*y1[jTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										-lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1));//c/va:fe/cr      
				
				dTc3(kTC1(a),jTC1(a)) += TC1(a)*y1[iTC1(a)]*pow((y1[jTC1(a)]-y1[iTC1(a)]),lTC1(a)) \
										+lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1));//c/va:cr/mo
				
			} 
			
			else
			{
				if ((y1[iTC1(a)]-y1[jTC1(a)])==0.0) {check = 1.0;}
				else {check = (y1[iTC1(a)]-y1[jTC1(a)]);}
				
				dTc1(iTC1(a),iTC1(a)) += lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										+lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Fe_Fe     
				
				dTc1(jTC1(a),jTC1(a)) += -lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										+lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2)) \
										-lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1));//Cr_Cr 
				
				dTc1(iTC1(a),jTC1(a)) += TC1(a)*y2[kTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Fe_Cr/Mo
				
				dTc1(jTC1(a),iTC1(a)) += TC1(a)*y2[kTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*TC1(a)*y1[iTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*TC1(a)*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-1)) \
										-lTC1(a)*(lTC1(a)-1)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*y2[kTC1(a)]*pow(check,(lTC1(a)-2));//Cr/Mo_Fe
				
				dTc2(iTC1(a),kTC1(a)) += TC1(a)*y1[jTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc2(jTC1(a),kTC1(a)) += TC1(a)*y1[iTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										-lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //Fe/Cr:C/Va
				
				dTc3(kTC1(a),iTC1(a)) += TC1(a)*y1[jTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										+lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //c/va:fe/cr
				
				dTc3(kTC1(a),jTC1(a)) += TC1(a)*y1[iTC1(a)]*pow((y1[iTC1(a)]-y1[jTC1(a)]),lTC1(a)) \
										-lTC1(a)*TC1(a)*y1[iTC1(a)]*y1[jTC1(a)]*pow(check,(lTC1(a)-1)); //c/va:cr/mo
				
			}
		//}
    }
	//cout<<"dTc3="<<dTc3<<endl;
	//cout<<"dTc1="<<dTc1<<endl;
	//cout<<"dTc2="<<dTc2<<endl;
	//::::::End:::::::::
	
	//::::::::::::::BMAG:::::::::::::
	//====================================BMAG=========================
	int ig2  = BMAGfinal.extent(0); 
	int jg2  = BMAGfinal.extent(1); 
	int iex2 = BMAG1.extent(0);
	
	Array <double,2> dBref1(ig2,jg2), dBref2(jg2,ig2), dB1(ig2,ig2), dB2(ig2,jg2), dB3(jg2,ig2);
	dBref1 = 0.;
	dBref2 = 0.;
	dB1    = 0.;
	dB2    = 0.;
	dB3    = 0.;
	
	for (i = 0; i < ig2; i++)
		for (j = 0; j < jg2; j++)
		{
			//if (p==0.28)
			//{	
			//	dBref1(i,j) +=  -BMAGfinal(i,j)/n;
			//}
			//else
			//{
				dBref1(i,j) +=  BMAGfinal(i,j);
			//}		 
		}
	
	for (i = 0; i < ig2; i++)
		for (j = 0; j < jg2; j++)
		{
			//if (p==0.28)
			//{ 	
			//	dBref2(j,i) += -BMAGfinal(i,j)/n;
			//}
			//else
			//{
				dBref2(j,i) += BMAGfinal(i,j);
			//}	
		}
	
	//cout<<"dBref1="<<dBref1<<endl;
	//cout<<"dBref2="<<dBref2<<endl;
	//cout<<"BMAG1 ="<<BMAG1<<endl;cout<<"iBMAG1 ="<<iBMAG1<<endl;cout<<"jBMAG1 ="<<jBMAG1<<endl;cout<<"kBMAG1 ="<<kBMAG1<<endl;cout<<"lBMAG1="<<lBMAG1<<endl;
	
	for(a = 0; a < iex2; a++)
	{
		double check = 0.0;
		/*if (p==0.28)
		{
			if ((jBMAG1(a)<2) && (def!=5))
			{
				if ((y1[jBMAG1(a)]-y1[iBMAG1(a)])==0.0) {check = 1.0;}
				else {check = (y1[jBMAG1(a)]-y1[iBMAG1(a)]);}	 
				
				dB1(iBMAG1(a),iBMAG1(a)) += -lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Fe_Fe 
				
				dB1(jBMAG1(a),jBMAG1(a)) += lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Cr_Cr   
				
				dB1(iBMAG1(a),jBMAG1(a)) += (-BMAG1(a)/n)*y2[kBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Fe_Cr/Mo
				
				dB1(jBMAG1(a),iBMAG1(a)) += (-BMAG1(a)/n)*y2[kBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Cr/Mo_Fe
				
				dB2(iBMAG1(a),kBMAG1(a)) += (-BMAG1(a)/n)*y1[jBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB2(jBMAG1(a),kBMAG1(a)) += (-BMAG1(a)/n)*y1[iBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB3(kBMAG1(a),iBMAG1(a)) += (-BMAG1(a)/n)*y1[jBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1));//c/va:fe/cr      
				
				dB3(kBMAG1(a),jBMAG1(a)) += (-BMAG1(a)/n)*y1[iBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1));//c/va:cr/mo
				
			}  
			else
			{
				if ((y1[iBMAG1(a)]-y1[jBMAG1(a)])==0.0) {check = 1.0;}
				else {check = (y1[iBMAG1(a)]-y1[jBMAG1(a)]);}	   
				
				dB1(iBMAG1(a),iBMAG1(a)) += lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Fe_Fe     
				
				dB1(jBMAG1(a),jBMAG1(a)) += -lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Cr_Cr 
				
				dB1(iBMAG1(a),jBMAG1(a)) += (-BMAG1(a)/n)*y2[kBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Fe_Cr/Mo
				
				dB1(jBMAG1(a),iBMAG1(a)) += (-BMAG1(a)/n)*y2[kBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Cr/Mo_Fe
				
				dB2(iBMAG1(a),kBMAG1(a)) += (-BMAG1(a)/n)*y1[jBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB2(jBMAG1(a),kBMAG1(a)) += (-BMAG1(a)/n)*y1[iBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB3(kBMAG1(a),iBMAG1(a)) += (-BMAG1(a)/n)*y1[jBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //c/va:fe/cr
				
				dB3(kBMAG1(a),jBMAG1(a)) += (-BMAG1(a)/n)*y1[iBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*(-BMAG1(a)/n)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //c/va:cr/mo
			}
		}*/
		
		//else
		//{			
			if ((jBMAG1(a)<2) && (def!=5))
			{
				if ((y1[jBMAG1(a)]-y1[iBMAG1(a)])==0.0) {check = 1.0;}
				else {check = (y1[jBMAG1(a)]-y1[iBMAG1(a)]);}	 
				
				dB1(iBMAG1(a),iBMAG1(a)) += -lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											-lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Fe_Fe 
				
				dB1(jBMAG1(a),jBMAG1(a)) += lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											+lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Cr_Cr   
				
				dB1(iBMAG1(a),jBMAG1(a)) += BMAG1(a)*y2[kBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Fe_Cr/Mo
				
				dB1(jBMAG1(a),iBMAG1(a)) += BMAG1(a)*y2[kBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Cr/Mo_Fe
				
				dB2(iBMAG1(a),kBMAG1(a)) += BMAG1(a)*y1[jBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB2(jBMAG1(a),kBMAG1(a)) += BMAG1(a)*y1[iBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB3(kBMAG1(a),iBMAG1(a)) += BMAG1(a)*y1[jBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1));//c/va:fe/cr      
				
				dB3(kBMAG1(a),jBMAG1(a)) += BMAG1(a)*y1[iBMAG1(a)]*pow((y1[jBMAG1(a)]-y1[iBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1));//c/va:cr/mo
				
			}  
			else
			{
				if ((y1[iBMAG1(a)]-y1[jBMAG1(a)])==0.0) {check = 1.0;}
				else {check = (y1[iBMAG1(a)]-y1[jBMAG1(a)]);}	   
				
				dB1(iBMAG1(a),iBMAG1(a)) += lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											+lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Fe_Fe     
				
				dB1(jBMAG1(a),jBMAG1(a)) += -lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											+lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2)) \
											-lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1));//Cr_Cr 
				
				dB1(iBMAG1(a),jBMAG1(a)) += BMAG1(a)*y2[kBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Fe_Cr/Mo
				
				dB1(jBMAG1(a),iBMAG1(a)) += BMAG1(a)*y2[kBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*BMAG1(a)*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-1)) \
											-lBMAG1(a)*(lBMAG1(a)-1)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*y2[kBMAG1(a)]*pow(check,(lBMAG1(a)-2));//Cr/Mo_Fe
				
				dB2(iBMAG1(a),kBMAG1(a)) += BMAG1(a)*y1[jBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB2(jBMAG1(a),kBMAG1(a)) += BMAG1(a)*y1[iBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //Fe/Cr:C/Va
				
				dB3(kBMAG1(a),iBMAG1(a)) += BMAG1(a)*y1[jBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											+lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //c/va:fe/cr
				
				dB3(kBMAG1(a),jBMAG1(a)) += BMAG1(a)*y1[iBMAG1(a)]*pow((y1[iBMAG1(a)]-y1[jBMAG1(a)]),lBMAG1(a)) \
											-lBMAG1(a)*BMAG1(a)*y1[iBMAG1(a)]*y1[jBMAG1(a)]*pow(check,(lBMAG1(a)-1)); //c/va:cr/mo
			}
		//} 
	}
	//cout<<"dB2="<<dB2<<endl;
	//cout<<"dB1="<<dB1<<endl;
	//cout<<"dB3="<<dB3<<endl;
	
	Array <double,2> dTC_1(ig1,jg1), dTC_2(jg1,ig1), dBmag_1(ig2,jg2), dBmag_2(jg2,ig2);
	dTC_1 = 0.;
	dTC_2 = 0.;
	dBmag_1= 0.;
	dBmag_2= 0.;
	
	//:::::::::::::
	//[c:c][c:va]
	//[va:c][va:va]
	Array <double,2> dTC_3(jg,jg),dBmag_3(jg,jg);
	dTC_3=0.;
	dBmag_3=0.;
	
	
	//:::::::::::::
	//[fe:c][fe:va]
	//[cr:c][cr:va]
	//[mo:c][mo:va]
	dTC_1 = dTcref1 + dTc2;
	//cout<<"dTC_1="<<dTC_1<<endl;
	
	//:::::::::::::
	//[c:fe][c:cr][c:mo][...]
	//[va:fe][va:cr][va:mo][...]
	dTC_2 = dTcref2 + dTc3;
	//cout<<"dTC_2="<<dTC_2<<endl;
	dBmag_1= dBref1 + dB2;
	//cout<<"dBmag_1="<<dBmag_1<<endl;
	dBmag_2= dBref2 + dB3;
	//cout<<"dBmag_2="<<dBmag_2<<endl;

	
	
	//:::::::
	//[fe:fe][fe:cr][fe:mo]
	//[cr:fe][cr:cr][cr:mo]
	//[mo:fe][mo:cr][mo:mo] 
	std::vector < std::vector<double> >MAG1(ig1);
	for(int i = 0 ; i < MAG1.size() ; i++)
	{
		MAG1[i].resize(ig1);
	}
	
	/*cout<<"derivTC1="<<derivTC1<<endl;
	cout<<"derivTC2="<<derivTC2<<endl;	
	cout<<"derivBMAG1="<<derivBMAG1<<endl;
	cout<<"derivBMAG2="<<derivBMAG2<<endl;
	cout<<"dB1="<<dB1<<endl;
	cout<<"dTc1="<<dTc1<<endl;
	cout<<"dB2="<<dB2<<endl;
	cout<<"dTc2="<<dTc2<<endl;

	
	cout<<"dTC_2="<<dTC_2<<endl;
	cout<<"dTC_1="<<dTC_1<<endl;
	
	cout<<"derivBMAG1="<<derivBMAG1<<endl;
	cout<<"Tc="<<Tc<<endl;
	cout<<"Bmag="<<Bmag<<endl;
	cout<<"derivTC1="<<derivTC1<<endl;
	cout<<"dTc1="<<dTc1<<endl;
	cout<<"dB1="<<dB1<<endl;*/
	
	MAG1 = magnetic_derivation2(T,Tc,Bmag,derivTC1,derivBMAG1,derivTC2,derivBMAG2,dTc1,dB1,2);
	/*cout<<"Tc="<<Tc<<endl;
	cout<<"T/Tc="<<T/Tc<<endl;
	cout<<"MAG1="<<"[fe:fe]"<<" "<<"[fe:cr]"<<endl;
	cout<<"     "<<"[cr:fe]"<<" "<<"[cr:cr]"<<endl;	
	cout<<"    "<<MAG1[0][0]<<" "<<MAG1[0][1]<<" "<<endl;
	cout<<"    "<<MAG1[1][0]<<" "<<MAG1[1][1]<<" "<<endl;
	printf ("\n");*/
	
	//:::::::
	//[c:fe] [c:cr] [c:mo]
	//[va:fe][va:cr][va:mo]
	std::vector < std::vector<double> >MAG2(jg1);
	for(int i = 0 ; i < MAG2.size() ; i++)
	{
		MAG2[i].resize(ig1);
	}
	
	MAG2 = magnetic_derivation2(T,Tc,Bmag,derivTC2,derivBMAG2,derivTC1,derivBMAG1,dTC_2,dBmag_2,1);
	
	
	/*cout<<"MAG2="<<"[c:fe]"<<" "<<"[c:cr]"<<endl;
	cout<<"     "<<"[va:fe]"<<" "<<"[va:cr]"<<endl;

	
	cout<<"    "<<MAG2[0][0]<<" "<<MAG2[0][1]<<endl;
	cout<<"    "<<MAG2[1][0]<<" "<<MAG2[1][1]<<endl;
	printf ("\n");*/
	
	
	 /*cout<<"derivTC1="<<derivTC1<<endl;
	 cout<<"derivBMAG1="<<derivBMAG1<<endl;
	cout<<"derivTC2="<<derivTC2<<endl;
	cout<<"derivBMAG2="<<derivBMAG2<<endl;
	 cout<<"TC="<<Tc<<endl;
	 cout<<"Bmag="<<Bmag<<endl;
	 cout<<"dTC_1="<<dTC_1<<endl;
	 cout<<"dBmag_1="<<dBmag_1<<endl;
	cout<<"T/TC="<<T/Tc<<endl;*/

	
	//:::::::
	//[fe:c][fe:va]
	//[cr:c][cr:va]
	//[mo:c][mo:va]
	std::vector < std::vector<double> >MAG3(ig1);
	for(int i = 0 ; i < MAG3.size() ; i++)
	{
		MAG3[i].resize(jg1);
	}
	
	MAG3 = magnetic_derivation2(T,Tc,Bmag,derivTC1,derivBMAG1,derivTC2,derivBMAG2,dTC_1,dBmag_1,1);
	
	/*cout<<"MAG3="<<"[fe:c]"<<" "<<"[fe:va]"<<endl;
	cout<<"     "<<"[cr:c]"<<" "<<"[cr:va]"<<endl;
	cout<<"     "<<"[mo:c]"<<" "<<"[mo:va]"<<endl;
	cout<<"     "<<MAG3[0][0]<<" "<<MAG3[0][1]<<endl;
	cout<<"     "<<MAG3[1][0]<<" "<<MAG3[1][1]<<endl;
	//cout<<"     "<<MAG3[2][0]<<" "<<MAG3[2][1]<<endl; 
	printf ("\n");*/
	//cout<<"MAG3="<<MAG3<<endl;	
	
	//:::::::
	//[c:c][c:va]
	//[va:c][va:va]	
    std::vector < std::vector<double> >MAG4(jg);
	for(int i = 0 ; i < MAG4.size() ; i++)
	{
		MAG4[i].resize(jg);
	}
	
	MAG4 = magnetic_derivation2(T,Tc,Bmag,derivTC2,derivBMAG2,derivTC1,derivBMAG1,dTC_3,dBmag_3,2);
	
	/*cout<<"derivTC2="<<derivTC2<<endl;
	cout<<"derivBMAG2="<<derivBMAG2<<endl;
	cout<<"TC="<<Tc<<endl;
	cout<<"Bmag="<<Bmag<<endl;*/
	/*cout<<"MAG4="<<"[c:c]"<<" "<<"[c:va]"<<endl;
	cout<<"     "<<"[va:c]"<<" "<<"[va:va]"<<endl;
	cout<<"     "<<MAG4[0][0]<<" "<<MAG4[0][1]<<endl;
	cout<<"     "<<MAG4[1][0]<<" "<<MAG4[1][1]<<endl;
	printf ("\n");*/
	//cout<<"MAG4="<<MAG4<<endl;	

	
	
	dG_1.resize(ig,ig), dG_2.resize(ig,jg), dG_3.resize(jg,jg), dG_4.resize(jg,ig);
	dG_1 = 0.0;
	dG_2 = 0.0;
	dG_3 = 0.0;
	dG_4 = 0.0;
	
	//::::
	//cout<<"[fe:fe]"<<" "<<"[fe:cr]"<<endl;// [fe:mo]
	//cout<<"[cr:fe]"<<" "<<"[cr:cr]"<<endl;// [cr:mo]
	//  [mo:fe] [mo:cr] [mo:mo] 
	for(int i = 0 ; i < ig ; i++)
	{
		for (int j = 0; j < ig; j++)
		{
			dG_1(i,j) = Gmix1(i,j) + dL3_1_2(i,j) + dLint1_2(i,j)+ MAG1[i][j];
		}
	} 
	//cout<<"dG_1= "<<dG_1<<endl;
	//cout<<"Gmix1="<<Gmix1<<endl;
	//cout<<"dL3_1_2="<<dL3_1_2<<endl;
	//cout<<"dLint1_2="<<dLint1_2<<endl;
	
	//::::
	//cout<<"[fe:c]"<<" "<<"[fe:va]"<<endl;
	//cout<<"[cr:c]"<<" "<<"[cr:va]"<<endl;
	//  [mo:c][mo:va]
	for(int i = 0 ; i < ig ; i++)
	{
		for (int j = 0; j < jg; j++)
		{
			dG_2(i,j) = Gref1(i,j) + dL3_2_2(i,j) + dLint1_3(i,j) + MAG3[i][j];
		}
    }
	//cout<<"dG_2="<<dG_2<<endl;
	//cout<<"dL3_2_2="<<dL3_2_2<<endl;
	//cout<<"dLint1_3="<<dLint1_3<<endl;
	//cout<<"Gref1="<<Gref1<<endl;
	//::::
	//cout<<"[c:c]"<<" "<<"[c:va]"<<endl;
	//cout<<"[va:c]"<<" "<<"[va:va]"<<endl;
	for(int i = 0 ; i < jg ; i++)
	{
		for (int j = 0; j < jg; j++)
		{
    dG_3(i,j) = Gmix2(i,j) + dLint2_4(i,j) + MAG4[i][j];
		}
	}	
	//cout<<"dG_3="<<dG_3<<endl;
	
	
	//for(int i = 0 ; i < jg ; i++)
	//{
	//	printf("\n");
	//	for (int j = 0; j < jg; j++)
	//	{
	//		printf("%e ",dG_3(i,j));
	//	}
	//}printf("\n");	
	
	
	//cout<<"dLint2_4="<<dLint2_4<<endl;
	//::::
	//cout<<"[c:fe]"<<" "<<"[c:cr]"<<endl;// [c:mo]
	//cout<<"[va:fe]"<<" "<<"[va:cr]"<<endl;//[va:mo]
	for(int i = 0 ; i < jg ; i++)
	{
		for (int j = 0; j < ig; j++)
		{
			dG_4(i,j) = Gref2(i,j) + dL3_3_2(i,j) + dLint2_5(i,j) + MAG2[i][j];
		}
	}
	//cout<<"dG_4="<<dG_4<<endl;
}

vector < std::vector<double> > Gibbs::chemical_potential2(std::vector <double> &y1,std::vector <double> &y2)
{
	
	int i, j, k;
	int size1 = y1.size();
	int size2 = y2.size();
	double dummy1 = 0.0, dummy2 = 0.0, dummytotal = 0.0;
	double dummy11=0.0, dummy21=0.0;
	int n = size1+size2;
	int m1 = size1 + size2 -1;
	std::vector < std::vector<double> > potentiel(m1);
	for(i = 0 ; i < potentiel.size() ; i++)
	{
		potentiel[i].resize(n);
	}  
	
	
	double Gm = 0;
	Gm = functionGibbs(T1,y1,y2);
	function_derivation_Gibbs(T1,y1,y2); 
	function_derivation_Gibbs2(T1,y1,y2);
	
	std::vector <double> d1(size1),d2(size2);
	//:::fe,cr,mo
	for (i = 0; i < size1; i++)
		for (j = 0; j < size1; j++)
		{
			d1[i] += y1[j]*dG_1(j,i);
		}
	for (i = 0; i < size1; i++)
	{
		d1[i] =d1[i]+dGm1[i];
	}
	
	for (j = 0; j < size1; j++)
		for (i = 0; i < size2; i++)
		{
			d1[j] += y2[i]*dG_4(i,j);
		}
	for (j = 0; j < size1; j++)
		dummy1+=d1[j];
	//::::c,va
	for (i = 0; i < size2; i++)
		for (j = 0; j < size2; j++)
		{
			d2[i]  += y2[j]*dG_3(j,i); 
		} 
	
	for (i = 0; i < size2; i++)  
	{
		d2[i] = dGm2[i] + d2[i];
	} 
	
	for (i = 0; i < size2; i++)
		for (j = 0; j < size1; j++)
		{
			d2[i] += y1[j]*dG_2(j,i);
		}
	
	//:::::carbon:::::
	//cout<<"dG_3="<<dG_3<<endl;
	std::vector <double> dummy(size2);
	for (j = 0; j < size2; j++)
	{
		dummy[j] =(dG_3(0,j) - dG_3(1,j))/m[1];
	}
    potentiel[0][n-2] = dummy[0];
    potentiel[0][n-1] = dummy[1];
	
	
	//:::::C, Va::::::::::::::::
	for (i = 0; i < size1; i++)
	{
		potentiel[0][i] = (dG_4(0,i) - dG_4(1,i))/m[1]; 
	}
	
	//printf ("\n");
	//for (i = 0; i < n; i++)
	//{
    //printf ("%g ",potentiel[0][i]);
	//}
	//printf ("\n");
	
	//::::Fe, Cr, Mo:::
	//cout<<"dG_1="<<dG_1<<endl;
	std::vector < std::vector<double> > pt(size1);
	for(i = 0 ; i < pt.size() ; i++)
	{
		pt[i].resize(size1);
	}  
	
	
	for (i = 0; i < size1; i++)
		for (j = 0; j < size1; j++)
		{
			potentiel[i+1][j] =/*(m[0]+m[1]*y2[0])*/dGm1[j]+(dG_1(i,j)+dG_4(1,j)-d1[j])*1./m[0];
		}
	
	for (i = 0; i < size1; i++)
	{
        potentiel[i+1][n-2]= /*(m[0] + m[1]*y2[0])*dGm2[0]+m[1]*Gm+*/dGm2[0]+1./m[0]*(dG_2(i,0)+dG_3(1,0)-d2[0]);
        potentiel[i+1][n-1]= /*(m[0] + m[1]*y2[0])*dGm2[1]+*/dGm2[1]+1./m[0]*(dG_2(i,1)+dG_3(1,1)-d2[1]);
	}
	
	
	printf ("\n");
	for (i = 0; i < m1; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf ("%lf ",potentiel[i][j]);
		} 
		printf ("\n");
	}
	printf ("\n");
	
	return potentiel;
}

//::::::::::::::::Potential chimique::::::::::::::
vector <double> Gibbs::chemical_potential(std::vector <double> &y1,std::vector <double> &y2)
{
	int i;
	int i1 = y1.size();
	int i2 = y2.size();
	double dummy1 = 0., dummy2 = 0., dummy = 0.;
	int n = i1 + 1;
	std::vector <double> pot_chemicalphase(n);
	
	//for(i=0; i< i1; i++)printf("y1[%i]=%e \n",i,y1[i]);
	//for(i=0; i< i2; i++)printf("y2[%i]=%e \n",i,y2[i]);
	//printf("\n");
	
	double Gm = 0;
	Gm = functionGibbs(T1,y1,y2); //cout<<"Gm ="<<Gm<<endl;
	
	function_derivation_Gibbs(T1,y1,y2); 
	
	//cout << dGm1[0]<<" "<<dGm1[1]<<endl;
	//cout << dGm2[0]<<" "<<dGm2[1]<<endl;
	for (i = 0; i < i1; i++)
	{
		
		dummy1 += y1[i] * dGm1[i];
		
	}
	
	for (i = 0; i < i2; i++)
	{
		
		dummy2 += y2[i] * dGm2[i];
		
	} 
	
	dummy = dummy1 + dummy2;
	//dummy = dummy1;

	
	pot_chemicalphase[0] = 1./m[1]*(dGm2[0] - dGm2[1]); //CARBONE
	
	for (i = 0; i < i1; i++)
	{
		
		pot_chemicalphase[i+1] = (m[0] + m[1]*y2[0])*Gm + 1/m[0]*(dGm1[i] + dGm2[1]- dummy);
		//pot_chemicalphase[i+1] = (m[0] + m[1]*y2[0])*Gm + 1/m[0]*(dGm1[i]- dummy);

		
	}
	
	//function_derivation_Gibbs2(T1,y1,y2);//ubrat'
    //printf ("Gm =%g  carbon=%g fer=%g chrome =%g\n",(m[0] + m[1]*y2[0])*Gm,pot_chemicalphase[0],pot_chemicalphase[1],pot_chemicalphase[2]);
	/*cout<<"dGm2[0]="<<dGm2[0]<<endl;
	cout<<"dGm2[1]="<<dGm2[1]<<endl;
	cout<<"dGm1[0]="<<dGm1[0]<<endl;
	cout<<"dGm1[1]="<<dGm1[1]<<endl;
	cout<<"Gm="<<Gm<<endl;*/
	//cout<<"mu_C-mu_FE="<<pot_chemicalphase[0]-pot_chemicalphase[1]<<endl;
	
	return pot_chemicalphase;
	
}

#include <cerrno>

//::::::La function de lire:::::::::::::
void Gibbs::read(double TT1)
{ 
	
	T1 = TT1;
	//cout<<"I'm here\n="<<T1;
	//========================================================================================================
	Array <double,1> Tr; Tr = 0;          
	Array <double,1> Gpower;Gpower = 0;           
	Array <double,1> L1power; L1power = 0;         
	Array <double,1> L2power; L2power = 0; 
	Array <double,1> L3power; L3power = 0;  //     
	Array <double,1> L0power; L0power = 0;        
	Array <double,1> TCpower; TCpower = 0;        
	Array <double,1> TC1power;TC1power = 0;       
	Array <double,1> BMAGpower; BMAGpower = 0;    
	Array <double,1> BMAG1power; BMAG1power = 0;
	Array <double,4> GHSERp, TC, BMAG;
	Array <double,6> L1ini, L2ini, TC1ini, BMAG1ini;
	Array <double,7> L0ini,L3ini; 
	Array <double,3> GHSERint, TCint, BMAGint; 
	Array <double,5> L1int, L2int, TC1int, BMAG1int;
	Array <double,6> L0int,L3int;
	Array <double,4> L1final, L2final, TC1final, BMAG1final;
	Array <double,5> L0final,L3final;
	
	int size2  = 0; // m
	int size3  = 0, size4  = 0, size5  = 0, size6  = 0; //pour GHSERp[temperature][FE/CR][C/VA][L'ORDRE]
	int size7  = 0, size8  = 0, size9  = 0, size10 = 0, size11 = 0, size12 = 0;//pour L1[temperature][FE][CR][C/VA][L'OR][L'OR]
	int size13 = 0, size14 = 0, size15 = 0, size16 = 0, size17 = 0, size18 = 0; //pour L2
	int size19 = 0, size20 = 0, size21 = 0, size22 = 0, size23 = 0, size24 = 0, size25 = 0; //pour L0
	int size26 = 0, size27 = 0, size28 = 0, size29 = 0; //pour TC
	int size30 = 0, size31 = 0, size32 = 0, size33 = 0, size34 = 0, size35 = 0; //pour TC1
	int size36 = 0, size37 = 0, size38 = 0, size39 = 0; //pour BMAG
	int size40 = 0, size41 = 0, size42 = 0, size43 = 0, size44 = 0, size45 = 0; //pour BMAG1ini
	int size46 = 0, size47 = 0, size48 = 0, size49 = 0, size50 = 0, size51 = 0, size52 = 0; //pour L3ini
	int sizeD  = 0, sizeQ  =0;
	
	//=========================================================================================================
	int nm; //Number of a matrix which it is filled
	fp=fopen(datafile,"rt");
	if(!fp)
	{
		printf("Error open file:%d",errno);
		exit(0);
	}
	char str[256];
	int nb_lignes = 0;
	nm = 0;
	bool wasNumbers = false;
	while (!feof(fp))
	{
		str[0] = 0;
		char *tmp = fgets(str, sizeof(str), fp);
		if (!tmp || str[0] == 0) //If have found nothing,
			break; //we leave
		
		if (str[0] == '/') //Comment
		{
			if (wasNumbers)//If we considered numbers from a matrix
			{
				wasNumbers = false;
				++nm; //Means we pass to the second matrix
			}
			continue;
		}
		
		int sp = calc_spaces(str); //We consider quantity of blanks
		
		
		if (sp == 4 && nm == 4) nm = 5;
		if (sp == 6 && nm == 7) nm = 8;
		if (sp == 7 && nm == 13) nm = 14;
		
		if (nm == 0)
		{
			sscanf(str, "%d", &size2);	
			m.resize(size2);
			wasNumbers = true;
		}
		if (nm == 1)
		{
			int i ;
			double v = 0.;
			sscanf(str, "%d %lf", &i, &v);
			m[i] = v;
			wasNumbers = true;
		}
		if (nm == 2)
		{
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				Tr.resize(sp+1);
				Tr(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}	 
			continue;
		}
		if (nm == 3)
		{
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				Gpower.resize(sp+1);
				Gpower(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}	
			continue;
		}
		if (nm == 4)
		{
			sscanf(str, "%d %d %d %d", &size3, &size4, &size5, &size6);
			GHSERp.resizeAndPreserve(size3,size4,size5,size6);GHSERp = 0;
			wasNumbers = true;
		}
		if (nm == 5)
		{
			int t , i , j , l;
			double v = 0.;
			sscanf(str, "%d %d %d %d %lf", &t, &i, &j, &l, &v);
			GHSERp(t,i,j,l) = v;
			wasNumbers = true;
		}
		if (nm == 6)
		{
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				L1power.resize(sp+1); 
				L1power(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}	

			continue;
		}
		if (nm == 7)
		{
			sscanf(str, "%d %d %d %d %d %d", &size7, &size8, &size9, &size10, &size11, &size12);
			L1ini.resizeAndPreserve(size7,size8,size9,size10,size11,size12);L1ini = 0;
			wasNumbers = true;
		}
		if (nm == 8)
		{
			int t , i, j, k, l, m;
			double v = 0.;
			sscanf(str, "%d %d %d %d %d %d %lf", &t, &i, &j, &k, &l, &m, &v);
			L1ini(t,i,j,k,l,m) = v;
			wasNumbers = true;	 
		}
		if (nm == 9)
		{
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				L2power.resize(sp+1);
				L2power(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}	 
			continue;
		}
		if (nm == 10)
		{
			sscanf(str, "%d %d %d %d %d %d", &size13, &size14, &size15, &size16, &size17, &size18);
			L2ini.resizeAndPreserve(size13,size14,size15,size16,size17,size18); L2ini = 0;
			wasNumbers = true;
		}
		if (nm == 11)
		{
			int t, i, j, k, l, m;
			double v = 0.;
			sscanf(str, "%d %d %d %d %d %d %lf", &t, &i, &j, &k, &l, &m, &v);
			L2ini(t,i,j,k,l,m) = v;
			wasNumbers = true;	
		}
		if (nm ==12)
		{
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				L3power.resize(sp+1); 
				L3power(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}
			continue;
		} 
		if (nm==13)
		{
			sscanf(str, "%d %d %d %d %d %d %d", &size46, &size47, &size48, &size49, &size50, &size51, &size52);
			L3ini.resizeAndPreserve(size46,size47,size48,size49,size50,size51,size52); L3ini = 0;
			wasNumbers = true;
		} 
		if (nm==14)
		{
			int t, i, j , k, l, m, n ;
			double v = 0.;
			sscanf(str, "%d %d %d %d %d %d %d %lf", &t, &i, &j, &k, &l, &m, &n, &v);
			L3ini(t,i,j,k,l,m,n) = v;
			wasNumbers = true;
		}
		if (nm == 15)
		{
			nm = 16;
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				L0power.resize(sp+1); 
				L0power(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}
			continue;
		}
		if (nm == 16)
		{
			sscanf(str, "%d %d %d %d %d %d %d", &size19, &size20, &size21, &size22, &size23, &size24, &size25);
			L0ini.resizeAndPreserve(size19,size20,size21,size22,size23,size24,size25); L0ini = 0;
			wasNumbers = true;	 
		}
		if (nm == 17)
		{  
			int t, i, j , k, l, m, n ;
			double v = 0.;
			sscanf(str, "%d %d %d %d %d %d %d %lf", &t, &i, &j, &k, &l, &m, &n, &v);
			L0ini(t,i,j,k,l,m,n) = v;
			wasNumbers = true;
		}
		if (nm == 18)
		{
			nm = 19;
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				TCpower.resize(sp+1); 
				TCpower(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}	 
		 continue;
		}
		if (nm == 19)
		{
			sscanf(str, "%d %d %d %d", &size26, &size27, &size28, &size29);
			TC.resizeAndPreserve(size26,size27,size28,size29); TC = 0;
			wasNumbers = true;
		}
		if (nm == 20)
		{
			int t , i , j , l ;
			double v = 0.;
			sscanf(str, "%d %d %d %d %lf", &t, &i, &j, &l, &v);
			TC(t,i,j,l) = v;
			wasNumbers = true;
		}
		if (nm == 21)
		{
			nm = 22;
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				TC1power.resize(sp+1); 
				TC1power(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}
			continue;
		}
		if (nm == 22)
		{
			sscanf(str, "%d %d %d %d %d %d", &size30, &size31, &size32, &size33, &size34, &size35);
			TC1ini.resizeAndPreserve(size30,size31,size32,size33,size34,size35); TC1ini = 0;
			wasNumbers = true;
		}
		if (nm == 23)
		{
			int t , i, j, k, l, m;
			double v = 0.;
			sscanf(str, "%d %d %d %d %d %d %lf", &t, &i, &j, &k, &l, &m, &v);
			TC1ini(t,i,j,k,l,m) = v;
			wasNumbers = true;
		}
		if (nm == 24)
		{
			nm = 25;
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				BMAGpower.resize(sp+1); 
				BMAGpower(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}
			continue;
		}
		if (nm == 25)
		{
			sscanf(str, "%d %d %d %d", &size36, &size37, &size38, &size39);
			BMAG.resizeAndPreserve(size36,size37,size38,size39); BMAG = 0;
			wasNumbers = true;
		}
		if (nm == 26)
		{
			int t , i , j , l ;
			double v = 0.;
			sscanf(str, "%d %d %d %d %lf", &t, &i, &j, &l, &v);
			BMAG(t,i,j,l) = v;
			wasNumbers = true;
		}
		if (nm == 27)
		{
			nm = 28;
			char *tmp = str;
			for (int i = 0; i<=sp; i++){
				double a = 0.;
				sscanf(tmp,"%lf",&a);
				BMAG1power.resize(sp+1); 
				BMAG1power(i) = a;
				tmp = strchr(tmp, ' ');
				if(!tmp)
					break;
				++tmp;
				wasNumbers = true;
			}
			continue;
		}
		if (nm == 28)
		{
			sscanf(str, "%d %d %d %d %d %d", &size40, &size41, &size42, &size43, &size44, &size45);
			BMAG1ini.resizeAndPreserve(size40,size41,size42,size43,size44,size45); BMAG1ini = 0;
			wasNumbers = true;
		}
		if (nm == 29)
		{
			int t , i, j, k, l, m;
			double v = 0.;
			sscanf(str, "%d %d %d %d %d %d %lf", &t, &i, &j, &k, &l, &m, &v);
			BMAG1ini(t,i,j,k,l,m) = v;
			wasNumbers = true;
		} 
		if (nm == 30)
		{
			sscanf(str, "%lf %lf",&p,&n);
			wasNumbers = true;
		}
		if (nm == 31)
		{
			sscanf(str, "%d", &sizeD);
			Dex.resize(sizeD);
			wasNumbers = true;
		}
		if (nm == 32)
		{
			int i ;
			float v = 0.;
			sscanf(str, "%d %e", &i, &v);
			Dex[i] = v;
			wasNumbers = true;
		}
		
		if (nm == 33)
		{
			sscanf(str, "%d", &sizeQ);
			Qex.resize(sizeQ);
			wasNumbers = true;
		}
		if (nm == 34)
		{
			int i ;
			float v = 0.;
			sscanf(str, "%d %e", &i, &v);
			Qex[i] = v;
			wasNumbers = true;
		}
	}
	fclose(fp);

	/* for (int i = 0; i<2;i++)
	 printf("%lf", m[i]);
	 cout<<"Tr="<<Tr<<endl;
	 cout<<"Gpower ="<<Gpower<<endl;
	 cout<<"GHSERp="<<GHSERp<<endl;
	 cout<<"L1power ="<<L1power<<endl;
	 cout<<"L1ini ="<<L1ini<<endl;
	 cout<<"L2power ="<<L2power<<endl;
	 cout<<"L2ini ="<<L2ini<<endl;
	 cout<<"L3power="<<L3power<<endl;
	 cout<<"L3ini="<<L3ini<<endl; 
	 cout<<"L0power ="<<L0power<<endl;
	 cout<<"L0ini ="<<L0ini<<endl;
	 cout<<"TCpower ="<<TCpower<<endl;
	 cout<<"TC ="<<TC<<endl;
	 cout<<"TC1ini"<<TC1ini<<endl;
	 cout<<"BMAG ="<<BMAG<<endl;
	 cout<<"BMAG1ini"<<BMAG1ini<<endl;
	 cout<<"TC1power"<<TC1power<<endl;
	 cout<<"BMAGpower="<<BMAGpower<<endl;
	 cout<<"BMAG1power ="<<BMAG1power<<endl;
	 cout<<"BMAG1ini="<<BMAG1ini<<endl;
	 printf ("%lf %lf \n",p,n);
	 
	for (int i = 0; i<3;i++)
	{printf("%e ", Dex[i]);} printf("\n");
	for (int i= 0; i < 3; i++)
	{printf("%e ", Qex[i]);}printf("\n");*/

	 
	
	GHSERint.resizeAndPreserve(size4,size5,size6); GHSERint = 0.;
	L1int.resizeAndPreserve(size8,size9,size10,size11,size12); L1int = 0.;
	L2int.resizeAndPreserve(size14,size15,size16,size17,size18); L2int = 0.;
	L3int.resizeAndPreserve(size47,size48,size49,size50,size51,size52); L3int = 0.;//
	L0int.resizeAndPreserve(size20,size21,size22,size23,size24,size25); L0int = 0.;
	TCint.resizeAndPreserve(size27,size28,size29); TCint = 0.;
	TC1int.resizeAndPreserve(size31,size32,size33,size34,size35); TC1int = 0.; 
	BMAGint.resizeAndPreserve(size37,size38,size39); BMAGint = 0.;
	BMAG1int.resizeAndPreserve(size41,size42,size43,size44,size45); BMAG1int = 0.;
	GHSERfinal.resizeAndPreserve(size4,size5); GHSERfinal = 0.;
	L1final.resizeAndPreserve(size8,size9,size10,size11); L1final = 0.;
	L2final.resizeAndPreserve(size14,size15,size16,size17); L2final = 0.;
	L3final.resizeAndPreserve(size47,size48,size49,size50,size51); L3final = 0.;//
	L0final.resizeAndPreserve(size20,size21,size22,size23,size24); L0final = 0.;
	TCfinal.resizeAndPreserve(size27,size28); TCfinal = 0.;
	TC1final.resizeAndPreserve(size31,size32,size33,size34); TC1final = 0.;
	BMAGfinal.resizeAndPreserve(size37,size38); BMAGfinal = 0.;
	BMAG1final.resizeAndPreserve(size41,size42,size43,size44); BMAG1final = 0.;
	L1.resizeAndPreserve(size8*size9*size10*size11);  L1 = 0.; 
	iL1.resizeAndPreserve(size8*size9*size10*size11); iL1 = 0;                   
	jL1.resizeAndPreserve(size8*size9*size10*size11); jL1 = 0;                   
	kL1.resizeAndPreserve(size8*size9*size10*size11); kL1 = 0;                   
	lL1.resizeAndPreserve(size8*size9*size10*size11); lL1 = 0; 
	L2.resizeAndPreserve(size14*size15*size16*size17); L2 = 0.;
	iL2.resizeAndPreserve(size14*size15*size16*size17);  iL2 = 0;                  
	jL2.resizeAndPreserve(size14*size15*size16*size17);  jL2 = 0;                  
	kL2.resizeAndPreserve(size14*size15*size16*size17);  kL2 = 0;                 
	lL2.resizeAndPreserve(size14*size15*size16*size17);  lL2 = 0; 
	L3.resizeAndPreserve(size47*size48*size49*size50*size51); L3 = 0.;//
	iL3.resizeAndPreserve(size47*size48*size49*size50*size51); iL3 = 0;//
	jL3.resizeAndPreserve(size47*size48*size49*size50*size51); jL3 = 0;//
	kL3.resizeAndPreserve(size47*size48*size49*size50*size51); kL3 = 0;//
	lL3.resizeAndPreserve(size47*size48*size49*size50*size51); lL3 = 0;//
	mL3.resizeAndPreserve(size47*size48*size49*size50*size51); mL3 = 0;//
	Ltetra.resizeAndPreserve(size20*size21*size22*size23*size24); Ltetra = 0.;
	iLtetra.resizeAndPreserve(size20*size21*size22*size23*size24); iLtetra = 0;
	jLtetra.resizeAndPreserve(size20*size21*size22*size23*size24); jLtetra = 0;
	kLtetra.resizeAndPreserve(size20*size21*size22*size23*size24); kLtetra = 0;
	lLtetra.resizeAndPreserve(size20*size21*size22*size23*size24); lLtetra = 0;
	mLtetra.resizeAndPreserve(size20*size21*size22*size23*size24); mLtetra = 0;
	TC1.resizeAndPreserve(size31*size32*size33*size34);  TC1 = 0.;
	iTC1.resizeAndPreserve(size31*size32*size33*size34); iTC1 = 0;                   
	jTC1.resizeAndPreserve(size31*size32*size33*size34); jTC1 = 0;                   
	kTC1.resizeAndPreserve(size31*size32*size33*size34); kTC1 = 0;                   
	lTC1.resizeAndPreserve(size31*size32*size33*size34); lTC1 = 0;
	BMAG1.resizeAndPreserve(size41*size42*size43*size44); BMAG1 = 0.;
	iBMAG1.resizeAndPreserve(size41*size42*size43*size44); iBMAG1 = 0;                   
	jBMAG1.resizeAndPreserve(size41*size42*size43*size44); jBMAG1 = 0;                   
	kBMAG1.resizeAndPreserve(size41*size42*size43*size44); kBMAG1 = 0;                   
	lBMAG1.resizeAndPreserve(size41*size42*size43*size44); lBMAG1 = 0;
	
	//There are calculations while the temperature varies
	//===================1===========
	rewriteGHSERp(T1,Tr,GHSERp,GHSERint);//cout<<"GHSERint="<<GHSERint<<endl;
	rewriteL(T1,Tr,L1ini,L1int);//cout<<"L1int"<<L1int<<endl;
	rewriteL(T1,Tr,L2ini,L2int);//cout<<"L2int"<<L2int<<endl;
	rewriteLtetra(T1,Tr,L0ini,L0int);//cout<<"L0int="<<L0int<<endl;
	rewriteLtetra(T1,Tr,L3ini,L3int);//cout<<"L3int="<<L3int<<endl;
	
	rewriteGHSERp(T1,Tr,TC,TCint);//cout<<"TCint="<<TCint<<endl;
	rewriteL(T1,Tr,TC1ini,TC1int);//cout<<"TC1int="<<TC1int<<endl;
	rewriteGHSERp(T1,Tr,BMAG,BMAGint);//cout<<"BMAGint="<<BMAGint<<endl;
	rewriteL(T1,Tr,BMAG1ini,BMAG1int);//cout<<"BAMG1int="<<BMAG1int<<endl;
	
	//===================2===========
	rewriteGHSERfinal(T1,GHSERint,Gpower,GHSERfinal);//cout<<"GHSERfinal="<<GHSERfinal<<endl;
	rewriteLfinal(T1,L1int,L1power,L1final);//cout<<"L1final ="<<L1final<<endl;
	rewriteLfinal(T1,L2int,L2power,L2final);//cout<<"L2final="<<L2final<<endl;
	rewriteL_tetra_final(T1,L3int,L3power,L3final);//cout<<"L3final="<<L3final<<endl;
	rewriteL_tetra_final(T1,L0int,L0power,L0final);//cout<<"L0final"<<L0final<<endl;
	rewriteGHSERfinal(T1,TCint,TCpower,TCfinal);//cout<<"TCfinal="<<TCfinal<<endl;
	rewriteLfinal(T1,TC1int,TC1power,TC1final);//cout<<"TC1final ="<<TC1final<<endl;
	rewriteGHSERfinal(T1,BMAGint,BMAGpower,BMAGfinal);//cout<<"BMAGfinal ="<<BMAGfinal<<endl;
	rewriteLfinal(T1,BMAG1int,BMAG1power,BMAG1final);//cout<<"BMAG1final ="<<BMAG1final<<endl;
	
	//===================3=============
	scanmatrixL(L1final,L1,iL1,jL1,kL1,lL1);
	//cout<<"L1 ="<<L1<<endl;cout<<"iL1 ="<<iL1<<endl;cout<<"jL1 ="<<jL1<<endl;cout<<"kL1 ="<<kL1<<endl;cout<<"lL1 ="<<lL1<<endl;
	scanmatrixL(L2final,L2,iL2,jL2,kL2,lL2);
	//cout<<"L2 ="<<L2<<endl;cout<<"iL2 ="<<iL2<<endl;cout<<"jL2 ="<<jL2<<endl;cout<<"kL2 ="<<kL2<<endl;cout<<"lL2 ="<<lL2<<endl; 
	scanmatrixLtetra(L3final,L3,iL3,jL3,kL3,lL3,mL3);
	//cout<<"L3="<<L3<<endl;cout<<"iL3="<<iL3<<endl;cout<<"jL3"<<jL3<<endl;cout<<"kL3="<<kL3<<endl;cout<<"lL3"<<lL3<<endl;
	scanmatrixLtetra(L0final,Ltetra,iLtetra,jLtetra,kLtetra,lLtetra,mLtetra);
	//cout<<"Ltetra="<<Ltetra<<endl;cout<<"iLtetra"<<iLtetra<<endl;cout<<"jLtetra="<<jLtetra<<endl;cout<<"kLtetra="<<kLtetra<<endl;
	//cout<<"lLtetra="<<lLtetra<<endl;
	scanmatrixL(TC1final,TC1,iTC1,jTC1,kTC1,lTC1);
	//cout<<"TC1 ="<<TC1<<endl;cout<<"iTC1 ="<<iTC1<<endl;cout<<"jTC1 ="<<jTC1<<endl;cout<<"kTC1 ="<<kTC1<<endl;cout<<"lTC1 ="<<lTC1<<endl;
	scanmatrixL(BMAG1final,BMAG1,iBMAG1,jBMAG1,kBMAG1,lBMAG1);
	//cout<<"BMAG1 ="<<BMAG1<<endl;cout<<"iBMAG1 ="<<iBMAG1<<endl;cout<<"jBMAG1 ="<<jBMAG1<<endl;cout<<"kBMAG1 ="<<kBMAG1<<endl;cout<<"lBMAG1 ="<<lBMAG1<<endl;
	
}
