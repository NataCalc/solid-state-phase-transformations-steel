#include "point.h"
#include <iostream>

using namespace std;
const double delta=1e-10;

/***************************************************************
*    default constructor (set all members of point to null)    *
****************************************************************/
point::point()
{
    coords=0;
    next=0;
    prev=0;
    function=0;
}

/***************************************************************
*    constructor build point from array of coordinates         *
*    coord - array of coordinates                              *
*    dim - dimension and size of array coordinate              *
*    func - function                                           *
****************************************************************/
point::point ( double* coord,int dim,int func )
{
    /*if ( dim<2&&dim>6 )
    {
        print_error ( "point::point(double*,int,int) : ERROR : source dimension too long. Not correct point?\n" );
        abort();
    }
    if ( !coord )
    {
        print_error ( "point::point(double*,int,int) : ERROR : source coordinates not allocated\n" );
        abort();
    }*/
	
	
    this->dim=dim;				//or (*this).dim
    function=func;  
    coords=new double[dim];     //Allocation memory under coordinates
    for ( int i=0;i<dim;i++ )
        coords[i]=coord[i];		//Filling of a vector by the coordinates
    next=0;						//Zeroing of the previous and following points
    prev=0;						
}
/***************************************************************
*                      copy constructor                        *
****************************************************************/
point::point ( point const& pt )
{
    /*if ( pt.dim<2&&pt.dim>6 )
    {
        print_error ( "point::point(point const&) : ERROR : source dimension too long. Not correct point?\n" );
        abort();
    }
    if ( !pt.coords )
    {
        print_error ( "point::point(point const&) : ERROR : source coordinates not allocated\n" );
        abort();
    }*/
    dim=pt.dim;
    function=pt.function;
    coords=new double[dim];
    for ( int i=0;i<dim;i++ )
        coords[i]=pt.coords[i];
    next=0;
    prev=0;
}

/***************************************************************
*    constructor build point from coordinates                  *
*    x,y - coordinates                                         *
*    func - function                                           *
****************************************************************/
point::point ( double x,double y,int func )
{
    dim=2;
    function=func;
    coords=new double[dim];
    coords[0]=x;
    coords[1]=y;
    next=0;
    prev=0;
}

/***************************************************************
*    constructor build point from coordinates                  *
*    x,y,z - coordinates                                       *
*    func - function                                           *
****************************************************************/
point::point ( double x,double y,double z,int func )
{
    dim=3;
    function=func;
    coords=new double[dim];
    coords[0]=x;
    coords[1]=y;
    coords[2]=z;
    next=0;
    prev=0;
}

/***************************************************************
*    constructor build point from coordinates                  *
*    x,y,z,t - coordinates                                     *
*    func - function                                           *
****************************************************************/
point::point ( double x,double y,double z,double t,int func )
{
    dim=4;
    function=func;
    coords=new double[dim];
    coords[0]=x;
    coords[1]=y;
    coords[2]=z;
    coords[3]=t;
    next=0;
    prev=0;
}

/***************************************************************
*    constructor build point from coordinates                  *
*    x,y,z,t,w - coordinates                                   *
*    func - function                                           *
****************************************************************/
point::point ( double x,double y,double z,double t,double w,int func )
{
    dim=5;
    function=func;
    coords=new double[dim];
    coords[0]=x;
    coords[1]=y;
    coords[2]=z;
    coords[3]=t;
    coords[4]=w;
    next=0;
    prev=0;
}

/***************************************************************
*    constructor build point from coordinates                  *
*    x,y,z,t,w,u - coordinates                                 *
*    func - function                                           *
****************************************************************/
point::point ( double x,double y,double z,double t,double w,double u,int func )
{
    dim=6;
    function=func;
    coords=new double[dim];
    coords[0]=x;
    coords[1]=y;
    coords[2]=z;
    coords[3]=t;
    coords[4]=w;
    coords[5]=u;
    next=0;
    prev=0;
}

/***************************************************************
*    destructor - free memory of point                         *
****************************************************************/
point::~point()
{
    if ( coords )
        delete[] coords;
}

/***************************************************************
*    iterator - return coordinate of point                     *
****************************************************************/
double point::operator[] ( int i ) const
{
    /*if ( !coords )
    {
        print_error ( "point::operator[]: ERROR: getting coordinate from not allocated point coords\n" );
        abort();
    }
    if ( i>=dim )
    {
        print_error ( "point::operator[]: ERROR: getting coordinate greater dimension\n" );
        abort();
    }*/
    return coords[i];
}

/***************************************************************
*    equal operator - compare two points                       *
****************************************************************/
bool point::operator== ( const point& pt )
{
    /*if ( pt.dim<2&&pt.dim>6 )
    {
        print_error ( "point::operator== : ERROR : source dimension too long. Not correct point?\n" );
        abort();
    }
    if ( dim!=pt.dim )
    {
        print_error ( "point::operetor== : ERROR : dimensions not equal\n" );
        abort();
    }
    if ( !pt.coords )
    {
        print_error ( "point::operator== : ERROR : source coordinates not allocated\n" );
        abort();
    }*/
    if ( &pt==this )
        return true;
    for ( int i=0;i<dim;i++ )
        //if(( (pt[i]<(coords[i]-delta)) || (pt[i]>(coords[i]+delta) ) ) )
        if ( pt[i]!=coords[i] )
            return false;
    return true;
}

/***************************************************************
*    not equal operator - use equal operator                   *
****************************************************************/
bool point::operator!= ( const point& pt )
{
    /*if ( pt.dim<2&&pt.dim>6 )
    {
        print_error ( "point::operator!= : ERROR : source dimension too long. Not correct point?\n" );
        abort();
    }
    if ( dim!=pt.dim )
    {
        print_error ( "point::operetor!=: ERROR : dimensions not equal\n" );
        abort();
    }
    if ( !pt.coords )
    {
        print_error ( "point::operator!= : ERROR : source coordinates not allocated\n" );
        abort();
    }*/
    if ( &pt==this )
        return false;
    if ( *this==pt )
        return false;
    return true;
}

/***************************************************************
*    assignment operator                                       *
****************************************************************/
point& point::operator= ( point const& pt )
{
    /*if ( pt.dim<2&&pt.dim>6 )
    {
        print_error ( "point::operator=(point const&) : ERROR : source dimension too long. Not correct point?\n" );
        abort();
    }
    if ( !pt.coords )
    {
        print_error ( "point::operator=(point const&) : ERROR : source coordinates not allocated\n" );
        abort();
    }
    if ( coords!=0 )
    {
        delete[] coords;
        coords=0;
    }*/
    dim=pt.dim;
    coords=new double[dim];
    for ( int i=0;i<dim;i++ )
        coords[i]=pt[i];
    function=pt.function;
    return *this;
}

/***************************************************************
*    substract operator - calculate lenght of vector           *
****************************************************************/
point point::operator- ( point& pt )
{
    double* coord=new double[dim];
    for ( int i=0;i<dim;i++ )
        coord[i]=coords[i]-pt[i];
    return point ( coord,pt.dim,0 );
}

/***************************************************************
*    convert point to CONCENTRATION_MOLE format                *
****************************************************************/
point& point::to_conc_mole_format()
{
    double* new_coords=new double[dim-1];
    double m[2];
    m[0]=1.0;
    m[1]= ( function==0 ) ?3.0:1.0;
    new_coords[0]=coords[0]*m[1]/ ( m[0]+m[1]*coords[0] );
    for ( int i = 1; i<dim-1; i++ )
        new_coords[i]=coords[i]*m[0]/ ( m[0]+m[1]*coords[0] );
    delete[] coords;
    coords=new_coords;
    dim-=1;
    //in_conc_mole_format=true;
    return *this;
}

/***************************************************************
*    convert point from CONCENTRATION_MOLE format              *
****************************************************************/
point& point::from_conc_mole_format()
{
    //if(!in_conc_mole_format)
    //  return *this;
    double* new_coords=new double[dim];
    double m[2];
    m[0]=1;
    m[1]= ( function==0 ) ?3.0:1.0;
    new_coords[0]= m[0]/m[1]*coords[0]/ ( 1-coords[0] ); // - dlya koordinati x, dlya pervoi koordinati
    for ( int i=1;i<dim;i++ )
        new_coords[i]= ( coords[i]* ( 1+m[1]/m[0]*new_coords[0] ) );  //dlya vsex ostal'nix
    delete[] coords;
    //in_conc_mole_format=false;
    coords=new_coords;
    return *this;
}

/***************************************************************
*    output operator for point objects                         *
****************************************************************/
ostream& operator<< ( std::ostream &os, const point &pt )
{
    os<<"(";
    for ( int i=0;i<pt.dimensions();i++ )
    {
        os<< ( double ) pt[i];
        if ( i!=pt.dimensions()-1 )
            os<<";";//";";
    }
    os<<"("<<pt.func() <<")";
    os<<")";
    return os;
}
