#include "newt.h"
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>

//extern fstream file_v;

using namespace std;
//=================================================
void newt::echange(double &a,double &b) const
{
   double tmp = a;
   a = b;
   b = tmp;
}

//========================Norm of a vector=========================================================
double newt::norme2(std::vector <double> &V)
{
	double norme = 0.0;
	for(unsigned int i=0; i<nn; i++)
		norme += V[i]*V[i];

	return sqrt(norme);
}

//==============Transpose of a matrice (Lines are interchanged the position i1 and i2)======================================
void newt::SwapRows(int i1, int i2,std::vector < std::vector<double> > &A)
{
  assert((i1 < nn) && (i2 < nn));
  for (int j = 0; j < nn; ++j)
    echange(A[i1][j], A[i2][j]);

}
//==============Return of an identical matrix===========================================================
void newt::IdentityMatrix(int n, std::vector < std::vector<double> > &A)
{
  std::vector < std::vector<double> > res(n);
  for(int i=0 ; i < res.size() ; i++)
  {
   res[i].resize(n);
  }
  for(int i=0 ; i < res.size() ; i++)
  res[i][i] = 1;
  A = res;
}

//==============LUP-decomposition (C= L + U - E)=======================================================
bool newt::LUP(std::vector < std::vector<double> > &A,std::vector < std::vector<double> > &C, std::vector < std::vector<double> > &P)
{
  IdentityMatrix(nn, P);

  for(int i=0 ; i < C.size() ; i++)
   for(int j=0 ; j < C[i].size() ; j++)
    {
      C[i][j] = A[i][j];
    }

  for( int i = 0; i < nn; i++ )
  {
  double pivotValue = 0;
  int pivot = -1;
  for( int row = i; row < nn; row++ )
    {
    if( fabs(C[row][i]) > pivotValue )
       {
        pivotValue = fabs(C[row][i]);
        pivot = row;
       }
    }
    if( pivotValue == 0 ) 
     {
    cout<<"matrix is singular"<<endl;
    printf ("i =%i pivot = %i \n",i,pivot);
    return false;
     }

    SwapRows(pivot, i, P);
    SwapRows(pivot, i, C);
     for( int j = i+1; j < nn; j++ ) 
      {
         C[j][i] /= C[i][i];
        for( int k = i+1; k < nn; k++ )
           C[j][k] -= C[j][i] * C[i][k];
      }
   }
	return true;	
}

//==============Inverse matrix a method LUP-decomposition=================================================
bool newt::Inverse(std::vector < std::vector<double> > &A, std::vector < std::vector<double> > &K)
{
  std::vector < std::vector<double> > P(nn);
  std::vector < std::vector<double> > C(nn);
  std::vector < std::vector<double> > X(nn);
  for(int i = 0; i < P.size(); i++)
   {
    P[i].resize(nn);
   }
  for(int i = 0; i < C.size(); i++)
   {
    C[i].resize(nn);
   }
 
  for(int i = 0 ; i < X.size() ; i++)
   {
    X[i].resize(nn);
   } 

	if(!LUP(A, C, P)) return false;
   for(int k = nn-1; k >= 0; k--) 
      {
        X[k][k] = 1;
        for( int j = nn-1; j > k; j--) X[k][k] -= C[k][j] * X[j][k];
        X[k][k] /= C[k][k];
        for( int i = k-1; i >= 0; i-- )
           {
            for( int j = nn-1; j > i; j-- ) 
              {
                X[i][k] -= C[i][j] * X[j][k];
                X[k][i] -= C[j][i] * X[k][j];
              }
            X[i][k] /= C[i][i];
           }
      }

 for(int i = 0 ; i < nn ; i++)
  for(int j = 0 ; j < nn ; j++)
   for (int k = 0; k < nn; k++)  
     {
       K[i][j] += X[i][k] * P[k][j];
     }
	return true;
}

double newt::fmin(std::vector <double> &x,std::vector <double> &concentration)
{
 int i;
 double sum;

 if (func == NULL)
 {
   printf("Chert, func is NULL!!!\n");
   return 0;
 }
 if (nrfuncv == NULL)
 {
   printf("Chert, nrfuncv is NULL!!!\n");
   return 0;
 }

  (func->*nrfuncv)(GetN(), x, fvec,concentration);

 for (sum = 0.0, i = 0; i < nn; ++i)
   sum += fvec[i]*fvec[i];
  return 0.5*sum;
}

//                                                                                           check
bool newt::calc(std::vector <double> &x, int *check,std::vector <double> &concentration, int s,bool &conv)
{
	conv=false;
	
	double counter=time(0);
	std::vector <double> g(nn), xold(nn);
	double test,den,temp;

  double fold,F= fmin(x,concentration); //fvec est calcule par fmin par appel de *nrfuncv
  double sum = 0;
  for (int i = 0; i < nn; ++i)
  sum += SQRnewt(x[i]);

	double stpmax = 0.0;
	
	if((s!=1) && (s!=3)){
		stpmax = FMAXnewt(sqrt(sum), (double)nn);
		    }
	
	else if(s==3){stpmax = 1.0e-1*FMAXnewt(sqrt(sum), (double)nn);}
	else if(s==1){stpmax = correc_STPMXnewt*STPMXnewt*FMAXnewt(sqrt(sum), (double)nn);}
	else if(s==5){stpmax = 1e-10*correc_STPMXnewt*STPMXnewt*FMAXnewt(sqrt(sum), (double)nn);}

    int iter = 0;
	
	if(s!=1) printf("\nNewton equilibre \n");
	if(s==1) printf("\nNewton cinetique \n");
	
    while (true)
    {
		    iter++;
			printf ("Iteration = %i \n",iter);
			printf ("F = %e \n",F);
		
			

			double x_c_fcc,x_cr_fcc,x_c_bcc,x_cr_bcc;
			x_cr_fcc = x[5] / (1. + x[6]) ;
			x_c_fcc  = x[6] / (1. + x[6]) ;
			x_cr_bcc = x[1] / (1. + 3. * x[2]) ;
			x_c_bcc  = 3. * x[2] / (1. + 3. * x[2]) ;
			
			printf("\n") ;
			printf("\nX(FCC , CR) = %g , X(FCC , C) =%g , X(BCC , CR) = %g , X(BCC , C) = %g\n",x_cr_fcc , x_c_fcc , x_cr_bcc , x_c_bcc) ;
			printf("\nVecteur Y : %e %e %e %e %e %e %e %e %e \n ", x[0] , x[1] , x[2] , x[3] , x[4] , x[5] , x[6] , x[7] , x[8]) ;	
		
		//file_v<<x_cr_fcc<<" "<<x_c_fcc<<" "<<x_cr_bcc<<" "<<x_c_bcc<<endl;

/*		double x_c_fcc,x_cr_fcc,x_c_bcc,x_cr_bcc;
		x_c_fcc  = x[1] / (1. + x[1]) ;
		x_c_bcc  = 3. * x[0] / (1. + 3. * x[0]) ;
		
		printf("\n") ;
		printf("\n X(FCC , C) =%g , X(BCC , C) = %g\n",x_c_fcc , x_c_bcc) ;
		printf("\nVecteur Y : %e %e %e \n ", x[0] , x[1] , x[2]) ;	*/
		

		/*if(s==1){
		double omega_c=(x_c_fcc-concentration[0])/(x_c_fcc-x_c_bcc);
			double omega_cr=(x_cr_fcc-concentration[2])/(x_cr_fcc-x_cr_bcc);
			printf("omega_c=%lf omega_cr=%lf \n",omega_c,omega_cr);
		
		
		double delta_c = 2 * (1. - omega_c)/omega_c ;
		double delta_cr = 2 * (1. - omega_cr)/omega_cr ;
		printf("delta_c=%lf delta_cr=%lf \n",delta_c,delta_cr);}*/
		
		
		(func->*nrfuncv)(nn, x, fvec,concentration);
		
		//for(int i=0;i<nn;i++) printf("\nx[%i] = %e",i,x[i]) ; printf("\n") ;  
		//for(int i=0;i<nn;i++) printf("\nFvec[%i] = %e",i,fvec[i]);
		
		
		
		/*TEST*/
		test = 0;
        for (int i = 0; i < nn; ++i)
			if (fabs(fvec[i]) > test)
				test = fabs(fvec[i]);
		if (test < 0.01*TOLXnewt)
		{
			printf("fabs(fvec[i])<0.01TOLXnewt)");
			conv=true;
			break;
		}
		
		
		std::vector < std::vector<double> > J(nn);
		for(int i = 0 ; i < J.size() ; i++)
		{
			J[i].resize(nn);
		} 
		
		std::vector < std::vector<double> > K(nn);
		for(int i = 0; i < K.size(); i++)
		{
			K[i].resize(nn);
		}
		(jacobian->*jacobfunc)(x,fvec1,J,concentration);

		
					/*printf("Matrice:\n");
						for(int i=0 ; i < nn ; i++)
						{for(int j=0 ; j < nn ; j++)
							{printf("%g ",J[i][j]);} printf("\n");}*/
		
		
		for (int i = 0; i < nn; ++i)
		{
			sum = 0;
			for (int j=0; j < nn; ++j)
                sum +=J[j][i]*fvec[j];
			g[i] = sum;
		}
		
        if(!Inverse(J,K))return false;
					/*printf("Matrice inverse:\n");
						for(int i=0 ; i < nn ; i++){
						for(int j=0 ; j < nn ; j++){
							printf("%g ",K[i][j]);}printf("\n");}*/
		
		
		std::vector <double> dx(nn);
		double indice;
        for (int k = 0; k < nn; ++k)
		{
			indice = 0.0;
			for (int l = 0; l < nn; ++l)
				indice += (-1.) * K[k][l] * fvec[l];
			dx[k] = indice; 
		}
		
		fold = F;
		for (int i = 0; i < nn; ++i)
			xold[i] = x[i];
		
		
					//printf("\n");printf ("Vector dx[i] avant lnsrch :");
					//for(int i=0 ; i < nn ; i++){printf ("\ndx[%i] = %e ",i,dx[i]);}printf("\n");
		lnsrch(nn, xold, fold, g, dx, x, &F, stpmax, check,concentration,s);
		
		
					//printf("\n");printf ("Vector dx après lnsrch :");
						//for(int i=0 ; i < nn ; i++){
							//printf ("\ndx[%i] = %e ",i,dx[i]);}printf("\n");
		
		
		/************************************************************************************************************/
						//printf("\n");printf ("Vector dx after party[i]:");
						//for(int i=0 ; i < nn ; i++){printf ("%g ",dx[i]);}printf("\n");
		
						//printf ("Vector x[i]:");for(int i=0 ; i < nn ; i++){
							//printf ("%g ",x[i]);}printf("\n");
		
		/*TEST*/
		double test = 0;
		for (int i = 0; i < nn; ++i)
			if (fabs(fvec[i]) > test)
				test = fabs(fvec[i]);
		if (test < epsilonnewt)
		{
			printf ("sortir1 : F < epsilonnewt\n");
			conv=true;
			break;
		}
		
		if (*check) 
		{ 
			test = 0.0;
			double den = FMAXnewt(F, 0.5*nn);
			for (int i = 0; i < nn; i++) 
			{
				temp =fabs(g[i])*FMAXnewt(fabs(x[i]),1.0)/den;
				if (temp > test) 
				{
					test=temp;
				}
			}
			*check = (test < TOLMIN ? 1 : 0);
			printf("Gradient de F nul : *check = (test < TOLMIN)\n");
			return true;
		}
		
		test = 0.0; 
        for (int i= 0; i < nn; i++) 
		{
			temp = (fabs(x[i]-xold[i]))/FMAXnewt(fabs(x[i]),1.0);
	        if (temp > test) test=temp;
		}
		if (test < TOLXnewt* 1.e-6 ){ printf ("sortir2 : Max ( dx / x ) < TOLXnewt \n");return false;}
		
		
		double norm = norme2(dx); //We consider norm
		//printf ("norm = %g \n",norm);
		
		//if (iter == 1)
		//break;
		
		if (norm < epsilonnewt) //If < eps 
		{
			printf ("sortir3\n");  
			conv=true;
			break;
		} //We leave
		
		if (s==1){
			if (iter>=300)
				return false;}
		if(s!=1){if(iter>=2000.)
			return false;}
		if (s==5){
			if (iter>=300)
				return false;}
		//break;
		
    }
	
	//fclose(fp);
	
	counter=time(0)-counter;
	cout<<"timenewt:"<<counter<<endl;
	
  *check = 0;
  return true;
}


//                                                                                                              dx                                  F                                                                    cin(1) ou eq(0)
void newt::lnsrch(int n, std::vector <double> &xold, double fold, std::vector <double> &g, std::vector <double> &p, std::vector <double> &x,double *f, double stpmax, int *check,std::vector <double> &concentration,int s)
{
									//printf("\nLancement de lnsrch\n");
	
  	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
	test,tmplam;
	
	*check=0;
	for (sum=0.0, i=0; i<n; i++)
	{  
		sum += p[i]*p[i];
	}
	sum=sqrt(sum);
	if (sum > stpmax)
	{
									//printf("\n remise à l'echelle de dx\n");
		for (i=0; i<n; i++)
		{
			p[i] *= stpmax/sum;		//Cette correction permet de rester dans des domaines realistes de X.
		}
	}
	for (slope=0.0, i=0; i<n; i++)
	{
		slope += g[i] * p[i];
	} 
	if (slope>=0.0) printf("Roundoff problem in lnsrch = %lf \n",slope);

	test=0.0;
	for (i=0; i<n; i++)
	{
		temp = fabs(p[i])/FMAXnewt(fabs(xold[i]),1.0);
		if (temp > test)
		{
			test = temp;
		}
	}
	alamin = TOLXnewt/test*1.e-8;
 	alam = 1.0;
	
	

	
	for (;;) 
	{
		for (i=0; i<n; i++)
		{  
			x[i] = xold[i] + alam * p[i];
		}



//Test de validite du domaine, equilibre

		if((s!=1) && (s!=3) && (s!=5))
		{


			bool flage = false;
			for (i = 0; i < n; i++)
			{
				if ((x[i]<0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]<0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((1.e-1),k);
						}
						k++;
					}
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((1.e-1),k);
						}
						k++; 
					}
				}
			}}


		if(s==3)
		{
			bool flage = false;
			for (i = 0; i < n; i++)
			{
				if ((x[i]<0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]<0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++;
					}
				}
			}
			if (flage)
			{
				printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++; 
					}
				}
			}}
		
		// Test de validite du domaine, cinetique
		if(s==1)
		{


			bool flage = false;
			for (i = 0; i < n-1; i++)
			{
				if ((x[i]<0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]<0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);//1e-6
						}
						k++;
					}
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++; 
					}
				}
			}
			
			double x_cr_bcc , x_cr_fcc; //alpha = bcc ; betta = fcc
			x_cr_bcc = x[1] / (1. + 3. * x[2] );
			x_cr_fcc = x[5] / (1. + x[6] ) ;

		}
		if(s==1)
		{

			if (x[n-1] < 0.)
			{
				printf("\nVitesse negative : %e",x[n-1]) ;
				int k = 1 ;
				while (x[n-1] < 0.)
				{
					for (int j = 0; j < n; j++)
					{
						x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
					}
					k++; 
				}
				printf("\nVitesse negative. Correction %i fois , facteur = %e ",k,pow((9.e-1),k)) ;	
			}
		}
		
		
		
				// Test de validite du domaine, vitesse
			
		//Test de validite du domaine, surstaturation positive et inferieure a 1
		if(s==1)
		{

			double v;
			double x_c_fcc,x_cr_fcc,x_c_bcc,x_cr_bcc;
			double xc0,xx0;
			double yca,yxa,ycg,yxg;
			double uca,uxa,ucg,uxg;
			double uc0,ux0;
			double omc,omx;
			
			yca=x[0];yxa=x[1];ycg=x[2];yxg=x[3];v=x[4];
			
			xc0 = concentration[0] ;
			xx0 = concentration[2] ;
			
			x_cr_fcc = yxg / (1. + ycg) ;
			x_c_fcc  = ycg / (1. + ycg) ;
			x_cr_bcc = yxa / (1. + 3. * yca) ;
			x_c_bcc  = 3. * yca / (1. + 3. * yca) ;
			
			uca = x_c_bcc / (1. - x_c_bcc) ;
			uxa = x_cr_bcc / (1. - x_c_bcc) ;
			ucg = x_c_fcc / (1. - x_c_fcc) ;
			uxg = x_cr_fcc / (1. - x_c_fcc) ;

			uc0 = xc0 / (1. - xc0) ;
			ux0 = xx0 / (1. - xc0) ;
			
			omc = (ucg - uc0) / (ucg-uca) ;
			omx = (uxg - ux0) / (uxg-uxa) ;
			
			int k=1 ;
//			printf("\nuxg=%e , ux0=%e , ucg=%e , uc0=%e , omx=%e , omc=%e\n",uxg,ux0,ucg,uc0,omx,omc) ;
			
			if(( uxg <= ux0 || ucg <= uc0 || omx >=1. || omx <=0. || omc >=1. || omc <=0. ) && uxg < 0.5  )
			{
			printf("\nerreur : uxg=%e < ux0=%e OU ucg=%e < uc0=%e OU omx = %e OU omc = %e\n",uxg,ux0,ucg,uc0,omx,omc) ;
				while(uxg <= ux0 || ucg <= uc0 || omx >=1. || omx <=0. || omc >=1. || omc <=0.)
				{
					for(int j = 0;j < n; j++)
					{
						x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
					}

					yca=x[0];yxa=x[1];ycg=x[2];yxg=x[3];v=x[4];
					x_cr_fcc = yxg / (1. + ycg) ;
					x_c_fcc  = ycg / (1. + ycg) ;
					x_cr_bcc = yxa / (1. + 3. * yca) ;
					x_c_bcc  = 3. * yca / (1. + 3. * yca) ;
					
					uca = x_c_bcc / (1. - x_c_bcc) ;
					uxa = x_cr_bcc / (1. - x_c_bcc) ;
					ucg = x_c_fcc / (1. - x_c_fcc) ;
					uxg = x_cr_fcc / (1. - x_c_fcc) ;
					
					omc = (ucg - uc0) / (ucg-uca) ;
					omx = (uxg - ux0) / (uxg-uxa) ;
					
					k++;
				}
				printf("\ncorrection %i fois, facteur= %e\n",k,pow(9.e-1,k)) ;
			}
		}
		
	
		
		// Test de validite du domaine, cinetique
		if(s==5)
		{
			bool flage = false;
			for (i = 0; i < n-1; i++)
			{
				if ((x[i]<=0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]<=0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);//1e-6
						}
						k++;
					}
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++; 
					}
				}
			}

		}
		
		*f = fmin(x,concentration);
		if (alam < alamin)
		{
			for (i = 0; i<n; i++)
			{
				x[i]=xold[i] + alam * p[i];
			}

			//printf("alam (%g) < alamin (%g)\n" , alam , alamin );
			*check=1;
			return;
		}
		else if (*f <= fold+ALFnewt*alam*slope/*1e-2*/)
		{
			//printf("diminution de F suffisante, alam = %e\n",alam);
			return;
		}

		
		else {                                                                          
			if (alam == 1.0)
			{
				//printf("\nalam==1 \n");
				if (s==1){tmplam = -slope/(2.0*(*f - fold - slope));return;}
				if (s!=1){tmplam = -slope/(2.0*(*f - fold - slope));return;} // Gradient de F nul
			}
			else  
			{		if (s==1) {printf("\nune seule iteration lnsrch, alam = %g",alam) ; return; }
					//printf ("alam!=1!, alam = %g\n" , alam );
				rhs1 = *f - fold - alam * slope;
				rhs2=f2 - /*fold2*/fold - alam2 * slope;
				a=(rhs1/(alam * alam) - rhs2/(alam2 * alam2))/(alam - alam2);
				b=(-alam2 * rhs1/(alam * alam) + alam * rhs2/(alam2 * alam2))/(alam - alam2);
				if (a == 0.0)
				{
					tmplam = -slope/(2.0 * b);//printf ("alam=0\n");
				}
				else 
				{
						//printf ("koren'\n");
					disc=b * b - 3.0 * a * slope;
					if (disc<0.0)
					{
						tmplam = 0.5 * alam;
						//printf("net kornei\n");
					}
					else if (b<=0.0)
					{
						tmplam=(-b + sqrt(disc))/(3.0 * a);
						//printf ("KU>0!\n");
					}  
					else
					{
						tmplam= -slope/(b + sqrt(disc));
						//printf ("KU<0!\n");
					}
				}
				if (tmplam > (0.5 * alam))
				{
					tmplam=0.5 * alam;
					//printf ("Net KU!\n"); 
				}
			}
		}
		alam2 = alam;
		f2 = *f;
		fold2 = fold;
		alam = FMAXnewt(tmplam,0.1*alam);
	}
}



//                                                                                                              dx                                  F                                                                    cin(1) ou eq(0)
void newt::lnsrch_old(int n, std::vector <double> &xold, double fold, std::vector <double> &g, std::vector <double> &p, std::vector <double> &x,double *f, double stpmax, int *check,std::vector <double> &concentration,int s)
{
									//printf("\nLancement de lnsrch\n");
	
  	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
	test,tmplam;
	
	*check=0;
	for (sum=0.0, i=0; i<n; i++)
	{  
		sum += p[i]*p[i];
	}
	sum=sqrt(sum);
	if (sum > stpmax)
	{
									//printf("\n remise à l'echelle de dx\n");
		for (i=0; i<n; i++)
		{
			p[i] *= stpmax/sum;		//Cette correction permet de rester dans des domaines realistes de X.
		}
	}
	for (slope=0.0, i=0; i<n; i++)
	{
		slope += g[i] * p[i];
	} 
	if (slope>=0.0) printf("Roundoff problem in lnsrch = %lf \n",slope);
	
	test=0.0;
	for (i=0; i<n; i++)
	{
		temp = fabs(p[i])/FMAXnewt(fabs(xold[i]),1.0);
		if (temp > test)
		{
			test = temp;
		}
	}
	alamin = TOLXnewt/test*1.e-8;
 	alam = 1.0;
	
	
	
	
	
	for (;;) 
	{
		for (i=0; i<n; i++)
		{  
			x[i] = xold[i] + alam * p[i];
		}



//Test de validite du domaine, equilibre

		if((s!=1) && (s!=3) && (s!=5))
		{
			bool flage = false;
			for (i = 0; i < n; i++)
			{
				if ((x[i]<0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]<0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((1.e-1),k);
						}
						k++;
					}
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((1.e-1),k);
						}
						k++; 
					}
				}
			}}

		if(s==3)
		{
			bool flage = false;
			for (i = 0; i < n; i++)
			{
				if ((x[i]<0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]<0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++;
					}
				}
			}
			if (flage)
			{
				printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++; 
					}
				}
			}}
		
		// Test de validite du domaine, cinetique
		if(s==1)
		{
			bool flage = false;
			for (i = 0; i < n-1; i++)
			{
				if ((x[i]<0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]<0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);//1e-6
						}
						k++;
					}
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++; 
					}
				}
			}
			
			double x_cr_bcc , x_cr_fcc; //alpha = bcc ; betta = fcc
			x_cr_bcc = x[1] / (1. + 3. * x[2] );
			x_cr_fcc = x[5] / (1. + x[6] ) ;

		}
		
		
		// Test de validite du domaine, cinetique
		if(s==5)
		{
			bool flage = false;
			for (i = 0; i < n-1; i++)
			{
				if ((x[i]<=0.0)||(x[i]>1.0))
				{
					flage = true;
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]<=0.0)
					{
						for(int j = 0;j < n; j++)
						{
							x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);//1e-6
						}
						k++;
					}
				}
			}
			if (flage)
			{
				//printf("\nCorrection : sortie du domaine de definition"); 
				int k = 0;
				for (i = 0; i < n-1; i++)
				{
					while(x[i]>1)
					{   
						for (int j = 0; j < n; j++)
						{
                            x[j] = xold[j] + alam * p[j]*pow((9.e-1),k);
						}
						k++; 
					}
				}
			}

		}
		
		*f = fmin(x,concentration);
		if (alam < alamin)
		{
			for (i = 0; i<n; i++)
			{
				x[i]=xold[i] + alam * p[i];
			}

			//printf("alam (%g) < alamin (%g)\n" , alam , alamin );
			*check=1;
			return;
		}
		else if (*f <= fold+ALFnewt*alam*slope/*1e-2*/)
		{
			//printf("diminution de F suffisante, alam = %e\n",alam);
			return;
		}

		
		else {                                                                          
			if (alam == 1.0)
			{
				//printf("\nalam==1 \n");
				if (s==1){tmplam = -slope/(2.0*(*f - fold - slope));return;}
				if (s!=1){tmplam = -slope/(2.0*(*f - fold - slope));return;} // Gradient de F nul
			}
			else  
			{		//if (s==1) {printf("\nune seule iteration lnsrch, alam = %g",alam) ; return; }
					//printf ("alam!=1!, alam = %g\n" , alam );
				rhs1 = *f - fold - alam * slope;
				rhs2=f2 - /*fold2*/fold - alam2 * slope;
				a=(rhs1/(alam * alam) - rhs2/(alam2 * alam2))/(alam - alam2);
				b=(-alam2 * rhs1/(alam * alam) + alam * rhs2/(alam2 * alam2))/(alam - alam2);
				if (a == 0.0)
				{
					tmplam = -slope/(2.0 * b);//printf ("alam=0\n");
				}
				else 
				{
						//printf ("koren'\n");
					disc=b * b - 3.0 * a * slope;
					if (disc<0.0)
					{
						tmplam = 0.5 * alam;
						//printf("net kornei\n");
					}
					else if (b<=0.0)
					{
						tmplam=(-b + sqrt(disc))/(3.0 * a);
						//printf ("KU>0!\n");
					}  
					else
					{
						tmplam= -slope/(b + sqrt(disc));
						//printf ("KU<0!\n");
					}
				}
				if (tmplam > (0.5 * alam))
				{
					tmplam=0.5 * alam;
					//printf ("Net KU!\n"); 
				}
			}
		}
		alam2 = alam;
		f2 = *f;
		fold2 = fold;
		alam = FMAXnewt(tmplam,0.1*alam);
	}
}
