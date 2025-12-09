#include "tools.h"
#include <cmath>
#include <iostream>

using namespace std;

double det ( double **p, int n )
{
    int i,h,x,y;
    if ( n==1 )
        return p[0][0];
    double d=0;
    //for ( i=0;i<n;i++ )
        //d+=p[i][i];
    //if ( d>-1e-10&&d<1e-10 )
        //return 0;
    d=0;
    for ( i=0; i<n; ++i )
    {
        double** add = new double*[n-1];
        for ( h=0; h<n-1; h++ )
            add[h] = new double[n-1];
        for ( y=0; y<n-1; y++ )
            for ( x=0; x<n; x++ )
            {
                if ( x==i ) continue;
                if ( x<i )
                {
                    if ( y<n )
                        add[x][y] = p[x][y+1];
                }
                else
                {
                    if ( y<n&&x<n )
                        add[x-1][y] = p[x][y+1];
                }
            }
        d += ( i%2?-p[i][0]:p[i][0] ) *det ( add, n-1 );
        for ( h=0; h<n-1; h++ )
            delete[] add[h];
        delete[] add;
    }
    return d;
}

double det ( vector<vector<double> >& p, int n )
{
    int i,h,x,y;
    if ( n==1 )
        return p[0][0];
    double d=0;
    for ( i=0;i<n;i++ )
        d+=p[i][i];
    if ( d>-1e-10&&d<1e-10 )
        return 0;
    d=0;
    for ( i=0; i<n; ++i )
    {
        double** add = new double*[n-1];
        for ( h=0; h<n-1; h++ )
            add[h] = new double[n-1];
        for ( y=0; y<n-1; y++ )
            for ( x=0; x<n; x++ )
            {
                if ( x==i ) continue;
                if ( x<i )
                {
                    if ( y<n )
                        add[x][y] = p[x][y+1];
                }
                else
                {
                    if ( y<n&&x<n )
                        add[x-1][y] = p[x][y+1];
                }
            }
        d += ( i%2?-p[i][0]:p[i][0] ) *det ( add, n-1 );
        for ( h=0; h<n-1; h++ )
            delete[] add[h];
        delete[] add;
    }
    return d;
}

void solve ( const double ** M,double const* b,double *x,int n )
{
    //________________________________________________________________
    // resheniye sistemi lineinix uravnenii metodom Kramera
    // opisanie metoda tut: http://en.wikipedia.org/wiki/Cramer's_rule
    // M - matrica (levaya chast' uravneniya
    // b - pravaya chast' uravneniya
    // x - sudi zapishutsya otveti
    // n - razmernost'


    //________________________________________________________________________
    //La décision du système des équations linéaires par la méthode de Kramera
    //La description de la méthode ici :http://en.wikipedia.org/wiki/Cramer's_rule
    //М est la matrice (une gauche partie de l'équation)
    //b est la partie droite de l'équation
    //dans x s'écrit ici la décision
    //n est la dimension

    //________________________________________________________________________
    //The decision of system of the linear equations Kramer's method
    //The method description here:http://en.wikipedia.org/wiki/Cramer's_rule
    //M is a matrix (the left member of equation)
    //b is the right member of equation
    //the decision is written in x
    //n is  dimension
    double** M1=new double*[n];
    for ( int i=0; i<n; i++ )
        M1[i] = new double[n];
    
  /*cout<<"before:"<<endl;
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
      cout<<M[i][j]<<" ";
    cout<<b[i]<<endl;
  }*/
    //_____________________________
    //vichisleniye opredelitelya
    //Le calcul du déterminant
    //Determinant calculation
    double dets=det ( ( double** ) M,n );
    double sum_x=0;

    //____________________________________________________________________________________________
    //kazhdii stolbec matritsi zamenyaem na vector b i vichislyaem opredelitel'
    //Chaque colonne de la matrice est remplacée par le vecteur b et nous calculons le déterminant
    //Each column of a matrix is replaced with a vector b and we calculate a determinant
    for ( int i=0;i<n;i++ )
    {

        for ( int j=0;j<n;j++ )
            for ( int k=0;k<n;k++ )
                M1[j][k]= ( ( k==i ) ? ( b[j] ) : ( M[j][k] ) );
        //__________________________________________________________________
        //delenie matritsi teckuchei iteratsii na glavnuu
        //La division de la matrice de l'itération en cours en la principale
        //Division of a matrix of current iteration into the main
        x[i]= ( det ( M1,n ) /dets );
    }

    //_______________________________
    //ochistka pamyati
    //Le nettoyage de la mémoire
    //Memory clearing
    for ( int i=0; i<n; i++ )
        delete[] M1[i];
    delete[] M1;
}

bool solve(const double** matr,double* b, int n)
{
  cout<<"before:"<<endl;
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
      cout<<matr[i][j]<<" ";
    cout<<b[i]<<endl;
  }
   bool res=true;
   double eps=.000001;
   double max;
   int max_i;
   double lead, a_div_lead;
   double** a=new double*[n];
   for(int i=0;i<n;i++)
   {
     a[i]=new double[n];
     for(int j=0;j<n;j++)
       a[i][j]=matr[i][j];
   }
   
   for(int k=0; k<n; k++)
   {

      max=0;
      max_i=-1;

      for(int i=k; i<n; i++)
      {
	 if(fabs(a[i][k])>max)
	 {
	    max=fabs(a[i][k]);
	    max_i=i;
	 }
      }

      if(max_i==-1 || fabs(a[max_i][k])<eps)
      {
	for(int i=0;i<n;i++)
	  delete[] a[i];
	delete[] a;
	cout<<"no squares"<<endl;
	return false;
        res=false;
	break;
      }

      lead=a[k][k];

      for(int j=k; j<n; j++)
        a[k][j]/=lead;
      b[k]/=lead;

      for(int i=0; i<n; i++)
      {
        a_div_lead=a[i][k]/a[k][k];

        if(i!=k)
        {
          for(int j=k; j<n; j++)
            a[i][j]-=a[k][j]*a_div_lead;
          b[i]-=b[k]*a_div_lead;
        }
      }
  }
  cout<<"after:"<<endl;
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
      cout<<a[i][j]<<" ";
    cout<<b[i]<<endl;
  }
  for(int i=0;i<n;i++)
    delete[] a[i];
  delete[] a;
  return res;
}