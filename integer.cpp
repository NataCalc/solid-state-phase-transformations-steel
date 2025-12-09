#include "integer.h"
#include <iostream>

using namespace std;

double x_c[10];
double x_cFeCrC[11];
double x_crFeCrC[11];


void integer::Eqdifp(double *t1,double *t2,double *t3,double *t4,double xi,double xf,double *Yi,
			int m,int p,int fi,vector <double> &concentration, \
			vector <double> &Ysol,std::vector < std::vector<double> > &P,std::vector < std::vector<double> > &P1)
{
	if (funcin == NULL)
	{
		printf("Chert, func is NULL!!!\n");
		return;
	}
	if (nrfuncin == NULL)
	{
		printf("Chert, nrfunc is NULL!!!\n");
		return;
	}
	
	
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
	for (i=1; i<m+1; i++) 
	{
		ni=(i-1)*fi-1;
		for (j=1; j<fi+1; j++) {
			x=xi+h*(ni+j);
			for (k=0; k<p+1; k++)  y[k]=z[k];
			for (k=0; k<p+1; k++)  ta[k]=h*(funcin->*nrfuncin)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)  y[k]=z[k]+ta[k]/2;
			x=x+h/2;
			for (k=0; k<p+1; k++)  tb[k]=h*(funcin->*nrfuncin)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)  y[k]=z[k]+tb[k]/2;
			for (k=0; k<p+1; k++)  tc[k]=h*(funcin->*nrfuncin)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)  y[k]=z[k]+tc[k];
			x=x+h/2;
			for (k=0; k<p+1; k++)  td[k]=h*(funcin->*nrfuncin)(k,x,y,concentration,Ysol);
			for (k=0; k<p+1; k++)
				z[k]=z[k]+(ta[k]+2*tb[k]+2*tc[k]+td[k])/6;
		}
		t1[i+1]=z[0];
		t2[i+1]=z[1];
		t3[i+1]=z[2];
		t4[i+1]=z[3];
	}

	Display(t1,t2,t3,t4,m,p+1,xi,xf);	
	//for (i=1; i<=m+1; i++){x_cFeCrC[i]=t1[i];x_crFeCrC[i]=t2[i];}
	(funcpot->*nrfuncpot)(m,xi,xf,t1,t2,P);
	Affiche(m,xi,xf,P);
	(funcderivFeCrC->*nrfuncderivFeCrC)(m,xi,xf,t1,t2,concentration,Ysol,P1);
	AfficheFeCrC(m,xi,xf,P1);
	
}

void integer::Display(double *t1,double *t2,double *t3,double *t4,int m,int p,double xi,double xf)
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


void integer::Affiche(int m,double xi,double xf,std::vector < std::vector<double> > &P)
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

void integer::AfficheFeCrC(int m,double xi,double xf,std::vector < std::vector<double> > &P)
{
	int i,j;
	double h,x;
	h=(xf-xi)/(m-1);
	x=xi-h;
	
	printf("\n");
    printf("		X		dmuC/dx		d(muCr-muFe)/dx		Lc		Lcrcr\n");
    printf(" -----------------------------\n");
    for (int i=1; i<m+1; i++)
	{
		x += h;
		printf("%9.6f		%e     %e	 %e		%e\n",x,P[i][0],P[i][1],P[i][2],P[i][3]); 
	}	
    printf(" -----------------------------\n");
}


void integer::runge_kutta(double *h,double *x, double *y,vector <double> &concentration, \
				 vector <double> &Ysol) 
{

    double a,b,c,d;
	a=*h*(funcin2d->*nrfuncin2d)(*x,*y,concentration,Ysol);
	b=*h*(funcin2d->*nrfuncin2d)(*x+*h/2,*y+a/2,concentration,Ysol);
	c=*h*(funcin2d->*nrfuncin2d)(*x+*h/2,*y+b/2,concentration,Ysol);
	*x+=*h;
	d=*h*(funcin2d->*nrfuncin2d)(*x,*y+c,concentration,Ysol);
	*y =*y + (a+b+b+c+c+d)/6;
}


double integer::equadiff_pc(double *tx,double *ty,double xi,double xf,double yi,int m,int n,int fi,vector <double> &concentration, \
				   vector <double> &Ysol, std::vector < std::vector<double> > &P,std::vector < std::vector<double> > &P1)
{
	//int vv = Ysol.size();
	//vector <double> Ysol1(vv),Ysol2(vv);
	//Ysol1=Ysol; Ysol2=Ysol;
	//VIT_TEMPORAIRE = Ysol[vv-1];
    double z=yi,y,w,p[4],x,h;
    int i,j,k; long ni;
    if ((m>100) || (fi<1)) return 0.;
    h = (xf-xi)/fi/m; //!!!!!!!
	//cout<<"xi="<<xi<<endl;
	//cout<<"yi="<<yi<<endl;
	p[3]=(funcin2d->*nrfuncin2d)(xi,yi,concentration,Ysol); //!!!!!!!
	//cout<<"p[3]="<<p[3]<<endl;
    tx[0]=xi; ty[0]=yi; //!!!!!!!!!!
    for (k=0, i=1; i<=m; i++) {
		ni=(long) (i-1)*fi-1;
		for (j=1; j<=fi; j++) 
		{
			x=xi+h*(ni+j);//!!!!!!!!!!!!!!!!
			if (++k<4) {
				runge_kutta(&h,&x,&z,concentration,Ysol);
				p[3-k]=(funcin2d->*nrfuncin2d)(x,z,concentration,Ysol);
			}
			else 
			{
				x+=h;
				w=z+h/24*(55*p[0]-59*p[1]+37*p[2]-9*p[3]);
				do 
				{
					y=w;
					w=z+h/24*(9*(funcin2d->*nrfuncin2d)(x,y,concentration,Ysol)+19*p[0]-5*p[1]+p[2]); 
					//printf("fabs(y-w)=%e \n",fabs(y-w));
				} while (fabs(y-w) > 1e-10);
				z=w; p[3]=p[2]; p[2]=p[1];
				p[1]=p[0]; p[0]=(funcin2d->*nrfuncin2d)(x,z,concentration,Ysol);
			}
		} // j loop
		tx[i]=x; ty[i]=z;
    } // i loop

	//double ty1[10];
	//for(i=0;i<=m;i++)
	//{ty1[i]=ty[i]/(1.+ty[i]);}	
	Affiche_pc(tx,ty,m);	
	(funcpot2d->*nrfuncpot2d)(tx,ty,concentration,Ysol,m,n,P);
	//for (int i=0; i<=m; i++){x_c[i]=ty[i];}
	Affiche_pot(tx,m,n,P);
	(funcderiv->*nrfuncderiv)(tx,ty,concentration,Ysol,m,n,P1);
	Affiche_potderiv(tx,m,n,P1);
	
	return -P[0][0]+P[m][0];
	
}

void integer::Affiche_pc(double *tx,double *ty,int n) 
{ 
    printf("\n");
    printf("        X               Y     \n");
    printf(" -----------------------------\n");
    for (int i=0; i<=n; i++){ 
		printf("%12.6f     %e\n",tx[i],ty[i]);} 
    printf(" -----------------------------\n");
}

void integer::Affiche_pot(double *tx,int m, int n, std::vector < std::vector<double> > &P)
{
	int i,j;
	printf("\n");
    printf("		X		C           FE     \n");
    printf(" -----------------------------\n");
    for (int i=0; i<=m; i++) 
		printf("%12.6f		%12.6f     %12.6f\n",tx[i],P[i][0],P[i][1]); 
    printf(" -----------------------------\n");
}

void integer::Affiche_potderiv(double *tx,int m, int n, std::vector < std::vector<double> > &P)
{
	int i,j;
	printf("\n");
    printf("		X	   dMU/dx	L	dMU/dc	\n");
    printf(" -----------------------------\n");
    for (int i=0; i<=m; i++) 
		printf("%12.6f		%e		%e		%e \n",tx[i],P[i][0],P[i][1],P[i][2]); 
    printf(" -----------------------------\n");
}


void integer::psevdo_FeC(double y, vector <double> &Ysol, std::vector < std::vector<double> > &P)
{
	(funcFeC->*nrfuncFeC)(y,Ysol,P);
}


double integer::Euler_Romberg(int NC, double H, double ER,double X0, double Y0, \
							double *X, double *Y, int *NL, \
							vector <double> &concentration,vector <double> &Ysol, \
							std::vector < std::vector<double> > &P, \
							std::vector < std::vector<double> > &P1)
{
	int N;
	double ET, XC,YC;
	double T[20];
	int J,K,L,LM,M,MM;


	
	//Initial conditions
	X[0]=X0; Y[0]=Y0;
	
	//main integration loop
	for (N=0; N<=NC; N++)
	{
		XC=X[N]; YC=Y[N];
		T[1]=Y[N] + H*(funcin2d->*nrfuncin2d)(XC,YC,concentration,Ysol);
		L=1; LM=2;
		
		do
		{
			XC=X[N]; YC=Y[N];
			for(J=1; J<=LM; J++)
			{
				XC=XC+H/LM; YC=YC+H/LM*(funcin2d->*nrfuncin2d)(XC,YC,concentration,Ysol);
			}
			T[L+1]=YC; M=1; K=L; MM=2; ET=1.0;
			if(K>1)
				do
				{
					T[K]=(MM*T[K+1]-T[K])/(MM-1);
					ET=fabs(T[K]-T[K-1]);
					M++; K--; MM *= 2;
				}while(ET>=ER && K>1);
			if(K==1)
			{
				L++; LM *= 2;
			}
		}while(L<10 && ET>=ER);
		X[N+1]=X[N]+H; Y[N+1]=T[K];
	}
	Affiche_pc(X,Y,NC);
	int n;
	(funcpot2d->*nrfuncpot2d)(X,Y,concentration,Ysol,NC,n,P);
	//for (int i=0; i<=NC; i++){x_c[i]=Y[i];}
	Affiche_pot(X,NC,n,P);
	(funcderiv->*nrfuncderiv)(X,Y,concentration,Ysol,NC,n,P1);
	Affiche_potderiv(X,NC,n,P1);
	
	return P[0][0]-P[NC][0];
}

double integer::Equadiff_Rung(double *t,double xi,double xf, double yi,int m,int fi, \
							  vector <double> &concentration,vector <double> &Ysol, \
							  std::vector < std::vector<double> > &P, \
							  std::vector < std::vector<double> > &P1)
{
	int i,j,ni;
	double a,b,c,d,h,x,y;
	
	if(fi<1) return 0;
	h = (xf-xi)/fi/(m-1);
	y=yi;
	t[1]=yi;

	
	for(i=1;i<m+1;i++){
		ni=(i-1)*fi-1;
		for(j=1;j<fi+1;j++){
			x=xi+h*(ni+j);
			a=h*(funcin2d->*nrfuncin2d)(x,y,concentration,Ysol);
			b=h*(funcin2d->*nrfuncin2d)(x+h/2,y+a/2,concentration,Ysol);
			c=h*(funcin2d->*nrfuncin2d)(x+h/2,y+b/2,concentration,Ysol);
			x=x+h;
			d=h*(funcin2d->*nrfuncin2d)(x,y+c,concentration,Ysol);
			y=y+(a+b+b+c+c+d)/6;
		}
		t[i+1]=y;
	}
	h = (xf-xi)/(m-1);
	x=xi-h;

	printf("\n      X           Y     \n");
	printf("------------------------\n");
	for (i=1; i<m+1; i++) {
		x+=h;
		printf(" %9.6f   %9.6f\n", x, t[i]);
	}
	
	double xx[m-1];
	double ty[m-1];
	x=xi-h;
	for(i=0;i<=m-1;i++){
		x+=h;
		xx[i]=x;}
	
	i=0;
	for(int j=1; j<m+1;j++){		
		t[i]=t[j];i++;}
	
	//for(i=0; i<=m-1; i++ )
	//{printf("x=%lf \n",xx[i]);}
	
	//for(i=0; i<=m-1; i++ )
	//{printf("x=%lf \n",t[i]);}
	int n;
	(funcpot2d->*nrfuncpot2d)(xx,t,concentration,Ysol,m-1,n,P);
	Affiche_pot(xx,m-1,n,P);
	(funcderiv->*nrfuncderiv)(xx,t,concentration,Ysol,m-1,n,P1);
	Affiche_potderiv(xx,m-1,n,P1);
	
	return P[0][0]-P[m-1][0];
}


