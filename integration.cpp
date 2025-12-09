/*
 *  integration.cpp
 *  
 *
 *  Created by nat on 06/07/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "integration.h"





double Wang(double r)
{
	//return r*r*(3.-2.*r);
	return r;
}

double dWang(double r)
{
	//return 6.*r-6.*r*r;
	return 1;
}


/**********************************************************
 *        SOLVING DIFFERENTIAL EQUATIONS OF ORDER 1        *
 *             (Prediction-correction method)              *
 * ------------------------------------------------------- */
/*********************************************************
 *             Prediction-correction method               *
 * ------------------------------------------------------ *
 * INPUTS:                                                *
 *	      xi      begin x value                        *
 *	      xf      end x value                          *
 *           y1      begin y value (at x=xi)              *
 *           m	      number of points to calculate        *
 *           fi      finesse (number of intermediate      *
 *                   points (for example 10)              *
 * OUTPUTS:                                               *
 *           tx      table of x values (m values)         *
 *           ty      table of y values (m values)         * 
 *                                                        *
 * DESCRIPTION:                                           *
 * The algorithm has the following steps:                 *
 *     1. calculate y2, y3, y4 using a classical Runge-   *
 *        Kutta method of order 4.                        *
 *     2. from point (x4, y4), first estimate y(n+1) by   *
 *        formula:                                        * 
 *        y(n+1)=y(n) + H/24(55y'(n)-59y'(n-1)+37y'(n-2)  *
 *               -9y'(n-3)                                *
 *        then continue with formula:                     *
 *        y(n+1)=y(n) + H/24(9y'(n+1)+19y'(n)-5y'(n-1)    *
 *               +y'(n-2),                                *
 *        noting that y'(n+1)=f(x(n+1),y(n+1)) with the   *
 *        estimated value of y(n+1), until convergence is *
 *        obtained.                                       *
 *********************************************************/

//classical Runge-Kutta method of order 4
void runge_kutta(double *h,double *x, double *y,vector <double> &concentration, \
				 vector <double> &Ysol,interface *iface,double (interface::*fp)(double , double ,vector <double> &, vector <double> &)) 
{
    double a,b,c,d;
	a=*h*(iface->*fp)(*x,*y,concentration,Ysol);
	b=*h*(iface->*fp)(*x+*h/2,*y+a/2,concentration,Ysol);
	c=*h*(iface->*fp)(*x+*h/2,*y+b/2,concentration,Ysol);
	*x+=*h;
	d=*h*(iface->*fp)(*x,*y+c,concentration,Ysol);
	*y =*y + (a+b+b+c+c+d)/6;
}


double equadiff_pc(double *tx,double *ty,double xi,double xf,double yi,int m,int n,int fi,vector <double> &concentration, \
				 vector <double> &Ysol,interface *iface,double (interface::*fp)(double , double ,vector <double> &, vector <double> &), \
				 std::vector < std::vector<double> > &P,void (interface::*fp1)(double *tx, double *ty,vector <double> &, \
																			   vector <double> &,int,int, \
																			   std::vector < std::vector<double> > &))
{
	//void runge_kutta(double h, double x, double *y,double (*fp)(double, double)); 
    double z=yi,y,w,p[4],x,h;
    int i,j,k; long ni;
    if ((m>100) || (fi<1)) return 0.;
    h = (xf-xi)/fi/m;
	//cout<<"xi="<<xi<<endl;
	//cout<<"yi="<<yi<<endl;
	p[3]=(iface->*fp)(xi,yi,concentration,Ysol);
	//cout<<"p[3]="<<p[3]<<endl;
    tx[0]=xi; ty[0]=yi;
    for (k=0, i=1; i<=m; i++) {
		ni=(long) (i-1)*fi-1;
		for (j=1; j<=fi; j++) 
			{
			x=xi+h*(ni+j);
			if (++k<4) {
				runge_kutta(&h,&x,&z,concentration,Ysol,iface,fp);
				p[3-k]=(iface->*fp)(x,z,concentration,Ysol);
			}
			else 
			{
				x+=h;
				w=z+h/24*(55*p[0]-59*p[1]+37*p[2]-9*p[3]);
				do 
				 {
					y=w;
					w=z+h/24*(9*(iface->*fp)(x,y,concentration,Ysol)+19*p[0]-5*p[1]+p[2]); 
					//printf("fabs(y-w)=%e \n",fabs(y-w));
				 } while (fabs(y-w) > 1e-10);
				z=w; p[3]=p[2]; p[2]=p[1];
				p[1]=p[0]; p[0]=(iface->*fp)(x,z,concentration,Ysol);
			 }
		} // j loop
		tx[i]=x; ty[i]=z;
    } // i loop
	Affiche_pc(tx,ty,m);
	(iface->*fp1)(tx,ty,concentration,Ysol,m,n,P);
	Affiche_pot(tx,m,n,P);
	return P[m][0]-P[0][0];
}

void Affiche_pc(double *tx,double *ty,int n) 
{ 
    printf("\n");
    printf("        X               Y     \n");
    printf(" -----------------------------\n");
    for (int i=0; i<=n; i++) 
		printf("%12.6f     %12.6f\n",tx[i],ty[i]); 
    printf(" -----------------------------\n");
}

void Affiche_pot(double *tx,int m, int n, std::vector < std::vector<double> > &P)
{
	int i,j;
	printf("\n");
    printf("		X			C               FE     \n");
    printf(" -----------------------------\n");
    for (int i=0; i<=m; i++) 
		printf("%12.6f		%12.6f     %12.6f\n",tx[i],P[i][0],P[i][1]); 
    printf(" -----------------------------\n");
}


void Affichefinal_1(int n, vector <double> &Y, \
					int m, vector <double> &Yder, \
					double G1, double G2, \
					double dGm2aC, double dGm2aVA, double dGm2bC,double dGm2bVA, \
					double ddGm2aCC, double ddGm2aCVA, double ddGm2aVAVA, double ddGm2aVAC, \
					double ddGm2bCC, double ddGm2bCVA, double ddGm2bVAVA, double ddGm2bVAC, \
					double x)
{
	int i,j;
	printf("Position  x :   %12.6f \n",x);
	printf("\n");
	printf("Fraction site :Fe_alfa,C_alfa,VAalfa,Fe_gamma,C_gamma,VAgamma \n");
	for(i=0; i<n; i++){printf("Y[%i]=%e ",i,Y[i]);}printf("\n");
	printf("\n");
	printf("1-ere derivation fraction site:(dYc/dc)alfa,(dYva/dc)alfa,(dYc/dc)gamma,(dYva/dc)gamma \n");
	for(i=0; i<m/2; i++){printf("Y[%i]=%e ",i,Yder[i]);}printf("\n");
	printf("\n");
	printf("2-eme derivation fraction site:(d(dYc/dc)dc)dc)alfa,(d(dYva/dc)dc)alfa,(d(dYc/dc)dc)gamma,(d(dYva/dc)dc)gamma \n");
	for(i=m/2; i<m; i++){printf("Y[%i]=%e ",i,Yder[i]);}printf("\n");
	printf("\n  Energy Gibbs : Galfa   Ggamma  \n");
	printf("%e	  %e\n",G1,G2);
	printf("\n 1-ere derivation l'energie de Gibbs : d(Galfa)/dYc   d(Galfa)/dYva  \n");
	printf("%e	  %e\n",dGm2aC,dGm2aVA);
	printf("\n 1-ere derivation l'energie de Gibbs : d(Ggamma)/dYc   d(Ggamma)/dYva  \n");
	printf("%e	  %e\n",dGm2bC,dGm2bVA);
	printf("\n 2-eme derivation l'energie de Gibbs : d(d(Galfa)/dYc)dYc d(d(Galfa)/dYc)dYva \n");
	printf("%e	  %e\n",ddGm2aCC,ddGm2aCVA);
	printf("\n 2-eme derivation l'energie de Gibbs : d(d(Galfa)/dYva)dYc d(d(Galfa)/dYva)dYva \n");
	printf("%e	  %e\n",ddGm2aVAC,ddGm2aVAVA);
	printf("\n 2-eme derivation l'energie de Gibbs : d(d(Ggamma)/dYc)dYc d(d(Ggamma)/dYc)dYva \n");
	printf("%e	  %e\n",ddGm2bCC,ddGm2bCVA);
	printf("\n 2-eme derivation l'energie de Gibbs : d(d(Ggamma)/dYva)dYc d(d(Ggamma)/dYva)dYva \n");
	printf("%e	  %e\n",ddGm2bVAC,ddGm2bVAVA);
	printf("\n");
}





/***************************************************************************
 *         SOLVING DIFFERENTIAL SYSTEMS WITH P VARIABLES OF ORDER 1         *
 *                 of type yi' = f(y1,y2,...,yn), i=1..n                    *
 * ------------------------------------------------------------------------ *
 *  INPUTS:                                                                 *
 *    m         number of points to display                                 *
 *    xi, xf    begin, end values of variable x                             *
 *    Yi        table of begin values of functions at xi                    *
 *    p         number of independant variables                             *
 *    fi        finesse (number of intermediary points)                     *
 *                                                                          *
 *  OUTPUTS:                                                                *
 *    t1,t2     real vectors storing the results for first two functions,   *
 *              y1 and y2.                                                  *
 * ------------------------------------------------------------------------ *      
 *  EXAMPLE:    y1'=y2+y3-3y1, y2'=y1+y3-3y2, y3'=y1+y2-3y3                 *
 *              Exact solution :  y1 = 1/3 (exp(-4x)  + 2 exp(-x))          *
 *                                y2 = 1/3 (4exp(-4x) + 2 exp(-x))          *
 *                                y3 = 1/3 (-5exp(-4x)+ 2 exp(-x))          *
 ***************************************************************************/
//Example: y1'=y2+y3-3y1, y2'=y1+y3-3y2, y3'=y1+y2-3y3
/*double fp2(int k, double x, double *y) 
{
	switch(k) {
		case 0: return y[1]+y[2]-3*y[0];
		case 1: return y[0]+y[2]-3*y[1];
		case 2: return y[0]+y[1]-3*y[2];
		default: return 0;
	}
}*/

/* Example #2: integrate system of equations from x=0 to PI*Sqrt(2) */
/* with initial conditions: y1(0)=3, y2(0)=0, y3(0)=4, y4(0)=0      */
/*double fp3(int k, double x, double *y) 
{
	switch(k) { 
		case 0: return y[1];
		case 1: return (-4*y[0]-3*y[2]);
		case 2: return y[3];
		case 3: return (-8*y[0]-2*y[2]);
		default: return 0;
	}
}*/

// print results (p>1)
void Display(double *t1,double *t2,double *t3,double *t4,int m,int p,double xi,double xf)
{
	int    i;
	double h,x;
	h=(xf-xi)/(m-1); 
	x=xi-h;
	printf("\n      X");
	for (i=1; i<=p; i++) printf("         Y%d", i);
	printf("\n"); 
	printf("--------------------------------------------------------\n");
	for (i=1; i<m+1; i++) {
		x += h;
		if (p==2)
			printf(" %9.6f  %e  %e \n",x,t1[i],t2[i]);
		else if (p==3)
			printf(" %9.6f  %9.6f  %9.6f  %9.6f\n",x,t1[i],t2[i],t3[i]);
		else if (p>3)	     
			printf(" %9.6f  %9.6f  %9.6f  %9.6f  %9.6f\n",x,t1[i],t2[i],t3[i],t4[i]);
	}
	printf("--------------------------------------------------------\n");
}

void Affiche_pot3d(int m,double xi,double xf,std::vector < std::vector<double> > &P)
{
	int i,j;
	double h,x;
	h=(xf-xi)/(m-1);
	x=xi-h;

	printf("\n");
    printf("		X		C		FE		CR		CR-FE\n");
    printf(" -----------------------------\n");
    for (int i=1; i<m+1; i++)
	{
		x += h;
		printf("%9.6f		%e     %e	%e	%e \n",x,P[i][0],P[i][1],P[i][2],P[i][3]); 
	}	
    printf(" -----------------------------\n");
}
/***************************************************************************
 *         SOLVING DIFFERENTIAL SYSTEMS WITH P VARIABLES OF ORDER 1         *
 *                 of type yi' = f(y1,y2,...,yn), i=1..n                    *
 * ------------------------------------------------------------------------ *
 *  INPUTS:                                                                 *
 *    m         number of points to display                                 *
 *    xi, xf    begin, end values of variable x                             *
 *    Yi        table of begin values of functions at xi                    *
 *    p         number of independant variables                             *
 *    fi        finesse (number of intermediary points)                     *
 *                                                                          *
 *  OUTPUTS:                                                                *
 *    t1,t2     real vectors storing the results for first two functions,   *
 *              y1 and y2.                                                  *
 * ------------------------------------------------------------------------ *      
 *  EXAMPLE:    y1'=y2+y3-3y1, y2'=y1+y3-3y2, y3'=y1+y2-3y3                 *
 *              Exact solution :  y1 = 1/3 (exp(-4x)  + 2 exp(-x))          *
 *                                y2 = 1/3 (4exp(-4x) + 2 exp(-x))          *
 *                                y3 = 1/3 (-5exp(-4x)+ 2 exp(-x))          *
 ***************************************************************************/

void Eqdifp(double *t1,double *t2,double *t3,double *t4,double xi,double xf,double *Yi,
            int m,int p,int fi,vector <double> &concentration, \
			vector <double> &Ysol,interface *iface,double (interface::*fp)(int k, double x, double *y,vector <double> &, vector <double> &), \
			std::vector < std::vector<double> > &P,void (interface::*fp1)(int, double, double,double *,double *, \
																		  std::vector < std::vector<double> > &P))
{

	double h,x;
	double ta[10],tb[10],tc[10],td[10],y[10],z[10];
	int    i,j,k,ni;
	if (fi<1) return;
	h = (xf - xi) / fi / (m-1);
	p--;
	t1[1]=Yi[0];
	t2[1]=Yi[1];
	t3[1]=Yi[2];
	t4[1]=Yi[3];
	for (k=0; k<p+1; k++) {
		y[k]=Yi[k]; z[k]=Yi[k];
	}
	for (i=1; i<m+1; i++) {
		ni=(i-1)*fi-1;
		for (j=1; j<fi+1; j++) {
			x=xi+h*(ni+j);
			for (k=0; k<p+1; k++)  y[k]=z[k];
			for (k=0; k<p+1; k++)  ta[k]=h*(iface->*fp)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)  y[k]=z[k]+ta[k]/2;
			x=x+h/2;
			for (k=0; k<p+1; k++)  tb[k]=h*(iface->*fp)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)  y[k]=z[k]+tb[k]/2;
			for (k=0; k<p+1; k++)  tc[k]=h*(iface->*fp)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)  y[k]=z[k]+tc[k];
			x=x+h/2;
			for (k=0; k<p+1; k++)  td[k]=h*(iface->*fp)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)
				z[k]=z[k]+(ta[k]+2*tb[k]+2*tc[k]+td[k])/6;
		}
		t1[i+1]=z[0];
		t2[i+1]=z[1];
		t3[i+1]=z[2];
		t4[i+1]=z[3];
	}
	Display(t1,t2,t3,t4,m,p+1,xi,xf);
	(iface->*fp1)(m,xi,xf,t1,t2,P);
	Affiche_pot3d(m,xi,xf,P);
}

void RK4n(int n, double x, double *y, double h, double *y1, \
		  vector <double> &concentration,vector <double> &Ysol, \
		  interface *iface,void (interface::*fp)(double x, double *,double *yp,vector <double> &, vector <double> &)) 
{
	
	
	double c1[10],c2[10],c3[10],c4[10],yy[10],h2; 
	int i;
	(iface->*fp)(x,y,c1,concentration,Ysol);
	h2=h/2.0;
	for (i=0; i<n; i++) yy[i]=y[i]+h2*c1[i];
	(iface->*fp)(x+h2,yy,c2,concentration,Ysol);
	for (i=0; i<n; i++) yy[i]=y[i]+h2*c2[i];
	(iface->*fp)(x+h2,yy,c3,concentration,Ysol);
	for (i=0; i<n; i++) yy[i]=y[i]+h*c3[i];
	(iface->*fp)(x+h,yy,c4,concentration,Ysol);
	for (i=0; i<n; i++)
		y1[i]=y[i]+/*h*c2[i];/*/h*(c1[i]+2.0*c2[i]+2.0*c3[i]+c4[i])/6.0;
		
}





 
