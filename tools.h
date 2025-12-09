#include <vector>

using namespace std;

double det ( double **p, int n );
double det ( vector<vector<double> >& p, int n );
void solve ( const double **const M,double const* b,double *x,int n );
bool solve( const double** matr,double* b, int n);
