#include "simplex.h"
#include <vector>
#include <cmath>

using namespace std;

//====================================================================================================
void simplex::amoeba(int nn,std::vector < std::vector<double> > &simplex,std::vector <double> &concentration, std::vector <double> &Y)
{
 int ihi, ilo, inhi, j;
 double fxr = 0.0, fxe = 0.0, fxs = 0.0, a = 0.0, b= 0.0;
 std::vector <double> midpoint(nn), next(nn), next_e(nn), next_s(nn);
 std::vector <double> fx(nn+1);
 evaluate_simplex(nn,simplex,concentration,fx);

 while (true)
  {
    simplex_extremes(nn,fx,&ihi,&ilo,&inhi);
    simplex_bearings(nn,simplex,midpoint,ihi);
    if(check_tol(fx[ihi],fx[ilo],tol)) {printf ("zdes\n");break;}
    fxr = update_simplex(nn,ihi,simplex,midpoint,next,1.0,concentration);

      //4
      if (fxr < fx[ilo])
       {
        fxe = update_stretch(nn,next,midpoint,next_e,2.0,concentration);
         
         //4a 
          if (fxe < fx[ilo])
           {
             for (j = 0 ; j < nn; j++)
              {
                simplex[ihi][j] = next_e[j];
              }
                evaluate_simplex(nn,simplex,concentration,fx);
                if (control(nn,fx))break;             
           }
             
             if (fxe > fx[ilo])
              {
                for (j = 0; j < nn; j++)
                  {
                    simplex[ihi][j] = next[j];
                  }
                 evaluate_simplex(nn,simplex,concentration,fx);
                 if (control(nn,fx))break;
              }
         }
              //4b
               if ((fxr < fx[inhi])&&(fxr > fx[ilo]))
                {
                  for (j = 0; j < nn; j++)
                   {
                     simplex[ihi][j] = next[j];
                   }
                 evaluate_simplex(nn,simplex,concentration,fx);
                 if (control(nn,fx))break;
                 }  

                //4c
                  if ((fxr < fx[ihi])&&(fxr > fx[inhi]))
                    {
                       for (j = 0; j < nn; j++)
                         {
                           a = next[j];
                           b = simplex[ihi][j];
                           next[j]  = b;
                           simplex[ihi][j] = a;
                         }      
                           a = fxr;
                           b = fx[ihi];
                           fxr     = b;
                           fx[ihi] = a; 

                           fxs = update_simplex(nn,ihi,simplex,midpoint,next_s,-0.5,concentration);
                    }
        
                   //4d
                      if (fxr > fx[ihi])
                       {
                           fxs = update_simplex(nn,ihi,simplex,midpoint,next_s,-0.5,concentration);
                       }
            //end 4
 
      //6  
        if (fxs < fx[ihi])
          {
              for (j = 0; j < nn; j++)
                {
                  simplex[ihi][j] = next_s[j];
                }
              evaluate_simplex(nn,simplex,concentration,fx);
              if (control(nn,fx))break;
          }
          
         //7
           if (fxs > fx[ihi])
             {
               contract_simplex(nn,ilo,simplex,fx,concentration);
             } 
    }
 
  //std::vector <double> Y(nn);
  for (j = 0; j < nn; j++)
   { 
     Y[j] = simplex[ilo][j];
     printf ("%e ",Y[j]);
   }printf ("\n");
        

}

//====================================================================================================
void simplex::evaluate_simplex(int nn,std::vector < std::vector<double> > &simplex,std::vector <double> &concentration, std::vector <double> &fx)
{
 int i;
 std::vector <double> Y(nn);
 for (i = 0; i < nn+1; i++)
  { 
    Y[0]=simplex[i][0];Y[1]=simplex[i][1];Y[2]=simplex[i][2];Y[3]=simplex[i][3];Y[4]=simplex[i][4];Y[5]=simplex[i][5];Y[6]=simplex[i][6];Y[7]=simplex[i][7];
    fx[i] = (func->*nrfuncv)(nn,Y,concentration);
  }
}

//===================================================================================================
void simplex::simplex_extremes(int nn,std::vector <double> &fx,int *ihi, int *ilo, int *inhi)
{
 int i;

  if (fx[0] > fx[1])
   {
     *ihi = 0; *ilo = *inhi = 1;
   }
   else
    {
      *ihi = 1; *ilo = *inhi = 0;
    }
 
   for (i = 2; i < nn+1; i++)
    {
      if (fx[i]<=fx[*ilo])
       { 
         *ilo = i;
       }
        else if (fx[i] > fx[*ihi])
         {
           *inhi = *ihi;
           *ihi  = i;
         }
          else if (fx[i] > fx[*inhi])
            {
              *inhi = i;
            } 
    }    
}

//=====================================================================================================
void simplex::simplex_bearings(int nn,std::vector < std::vector<double> > &simplex,std::vector <double> &midpoint,int ihi)
{
  int i, j;
  int alfa = 1;
  std::vector <double>  dummy(nn); 
  for (j = 0; j < nn; j++)
    {
      midpoint[j] = 0.0;
    }

 for (j = 0; j < nn; j++)
  {
    for (i = 0; i < nn+1; i++)
     {
      if (i != ihi)
       {
          dummy[j]     = simplex[i][j];
          midpoint[j] += dummy[j];
       }
     }
  }

  for (j = 0; j < nn; j++)
   {
     midpoint[j] /=nn;
   }
}

//======================================================================================================
int simplex::check_tol(double fmax, double fmin, double ftol)
{
 double delta = fabs(fmax-fmin);
 double accuracy = (fabs(fmax)+fabs(fmin))*ftol;
 return (delta < (accuracy + ZEPS));
}

//======================================================================================================
double simplex::update_simplex(int nn,int ihi,std::vector < std::vector<double> > &simplex,std::vector <double> &midpoint,std::vector <double> &next,double scale, std::vector <double> &concentration)
{
 int i;

 double fx = 0.0;

 for (i = 0; i < nn; i++)
  {
   next[i] = (1+scale)*midpoint[i] - scale*simplex[ihi][i]; 
  }

   fx   = (func->*nrfuncv)(nn,next,concentration);
  
 return fx;
}

//==========================================================================================================================
double simplex::update_stretch(int nn,std::vector <double> &next,std::vector <double> &midpoint,std::vector <double> &next_e,double scale, std::vector <double> &concentration)
{
  int i;
  double fxe = 0.0;
 
  for (i = 0; i < nn; i++)
   {
      next_e[i] = (1-scale)*midpoint[i]+scale*next[i];
   }


   fxe = (func->*nrfuncv)(nn,next_e,concentration);

 return fxe;
}
//===========================================================================================================================

void simplex::contract_simplex(int nn, int ilo,std::vector < std::vector<double> > &simplex,std::vector <double> &fx,std::vector <double> &concentration)
{
  int i, j;
  std::vector <double> Y(nn);
  double dummy = 0.0;  

  for (int i = 0; i < nn+1; i++)
    {
      if (i!=ilo)
       {
         for (j = 0; j < nn; j++)
          {
            dummy = (simplex[i][j]-simplex[ilo][j])*0.5;
            simplex[i][j] = simplex[ilo][j]+ dummy;
            Y[0]=simplex[i][0];Y[1]=simplex[i][1];Y[2]=simplex[i][2];Y[3]=simplex[i][3];Y[4]=simplex[i][4];Y[5]=simplex[i][5];Y[6]=simplex[i][6];Y[7]=simplex[i][7];
              fx[i] = (func->*nrfuncv)(nn,Y,concentration);
          }
       }
    }
}

//===========================================================================================================================
int simplex::control(int nn, std::vector <double> fx)
{
 int i;
 double f = 0.0, s = 0.0;

  for (i = 0; i < nn+1; i++)
   {
     f += fx[i]/(nn+1);
   }

    for (i = 0; i < nn+1; i++)
     {
       s +=pow((fx[i]-f),2.)/(nn+1);
     }

  return (sqrt(s)<tol);
}



