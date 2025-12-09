#include "hyperplane.h"
#include "tools.h"
#include "cstdio"
#include "iostream"
#include <cmath>



const double delta=1e-10;

/***************************************************************
*  default constructor (set all members of hyperplane to null) *
****************************************************************/
hyperplane::hyperplane()
{
    points=0;
    edges=0;
    equation=0;
    prev=next=0;
    dim=0;
    visited=old=false;
    border=false;
    no_points=false;
}

/***************************************************************
*                   copy constructor                           *
****************************************************************/
hyperplane::hyperplane ( hyperplane const& plane )
{
    /*if ( plane.dim<2&&plane.dim>6 )
    {
        print_error ( "hyperplane::hyperplane(hyperplane const&) : ERROR : source dimension too long. Not correct hyperplane?\n" );
        abort();
    }
    if ( !plane.points )
    {
        print_error ( "hypeplane::hyperplane(hyperplane&) : ERROR : dimensions not equal\n" );
        abort();
    }*/
    dim=plane.dim;
    points=new point[dim];		// allocte memory
    for ( int i=0;i<dim;i++ )	// copy all points
        points[i]=plane[i];
    init_plane ( false );		// initialize edges and hyperplane equation
    for ( int i=0;i<dim;i++ )	// copy pointers to neighbors
        edges[i].neighbor=plane.edges[i].neighbor; //rus:kopirovanie sosedei isxodnoi
												   //ploskosti v sozdavaemuu (gran' soseda= gran' soseda ploskosti)
    outside=plane.outside;
    visited=old=false;
    next=prev=0;
    border=false;
    no_points=false;
}

/***************************************************************
*    constructor build hyperplane from points                  *
*    pt1,pt2 - points                                          *
*    no_equation - need calculate equation or no               *
****************************************************************/
hyperplane::hyperplane ( point& pt1,point& pt2,bool no_equation )
{
    dim=2;
    points=new point[dim];
    points[0]=pt1;
    points[1]=pt2;
    init_plane ( no_equation );
    visited=old=false;
    border=false;
    no_points=false;
}

/***************************************************************
*    constructor build hyperplane from points                  *
*    pt1,pt2,pt3 - points                                      *
*    no_equation - need calculate equation or no               *
****************************************************************/
hyperplane::hyperplane ( point& pt1,point& pt2,point& pt3,bool no_equation )
{
    dim=3;
    points=new point[dim];
    points[0]=pt1;
    points[1]=pt2;
    points[2]=pt3;
    init_plane ( no_equation );
    visited=old=false;
    border=false;
    no_points=false;
}

/***************************************************************
*    constructor build hyperplane from points                  *
*    pt1,pt2,pt3,pt4 - points                                  *
*    no_equation - need calculate equation or no               *
****************************************************************/
hyperplane::hyperplane ( point& pt1,point& pt2,point& pt3,point& pt4,bool no_equation )
{
    dim=4;
    points=new point[dim];
    points[0]=pt1;
    points[1]=pt2;
    points[2]=pt3;
    points[3]=pt4;
    init_plane ( no_equation );
    visited=old=false;
    border=false;
    no_points=false;
}

/***************************************************************
*    constructor build hyperplane from points                  *
*    pt1,pt2,pt3,pt4,pt5 - points                              *
*    no_equation - need calculate equation or no               *
****************************************************************/
hyperplane::hyperplane ( point& pt1,point& pt2,point& pt3,point& pt4,point& pt5,bool no_equation )
{
    dim=5;
    points=new point[dim];
    points[0]=pt1;
    points[1]=pt2;
    points[2]=pt3;
    points[3]=pt4;
    points[4]=pt5;
    init_plane ( no_equation );
    visited=old=false;
    border=false;
    no_points=false;
}

/***************************************************************
*    constructor build hyperplane from points                  *
*    pt1,pt2,pt3,pt4,pt5,pt6 - points                          *
*    no_equation - need calculate equation or no               *
****************************************************************/
hyperplane::hyperplane ( point& pt1,point& pt2,point& pt3,point& pt4,point& pt5,point& pt6,bool no_equation )
{
    dim=6;
    points=new point[dim];
    points[0]=pt1;
    points[1]=pt2;
    points[2]=pt3;
    points[3]=pt4;
    points[4]=pt5;
    points[5]=pt6;
    init_plane ( no_equation );
    visited=old=false;
    border=false;
    no_points=false;
}

/***************************************************************
*  constructor build hyperplane from point and array of points *
*    pt - point                                                *
*    pts - array of points (array size dim-1)                  *
*    dim - dimension                                           *
****************************************************************/
hyperplane::hyperplane ( point& pt,point** pts,int dim )
{
    int pt_dim=pt.dimensions();
    /*if ( pt_dim<2&&pt_dim>6 )
    {
        print_error ( "hyperplane::hyperplane(point&,point*,int) : ERROR : source point dimension too long. Not correct point?\n" );
        abort();
    }*/
    this->dim=dim;
    points=new point[dim];
    points[dim-1]=pt;		// copy point
    for ( int i=0;i<dim-1;i++ ) // copy points from array
        points[i]=* ( pts[i] );
    init_plane ( false );	// init hyperplane
    visited=old=false;
    border=false;
    no_points=false;
}

/***************************************************************
*                 destructor - free memory                     *
****************************************************************/
hyperplane::~hyperplane()
{
    if ( points )
        delete[] points;
    if ( equation )
        delete[] equation;
    if ( !edges )
        return;
    for ( int i=0;i<dim;i++ )
    {
        delete[] edges[i].points;
    }
    delete[] edges;
}

/***************************************************************
*           iterator  - return point with i index              *
****************************************************************/
point& hyperplane::operator[] ( int i ) const
{
    /*if ( i>=dim||i<0 )
    {
        print_error ( "hyperplane::operator[] : ERROR : index too many long\n" );
        abort();
    }*/
    return points[i];
}

/***************************************************************
*         equal operator - compare two hyperplanes             *
****************************************************************/
bool hyperplane::operator== ( const hyperplane& plane ) const
{
    /*if ( plane.dim<2&&plane.dim>6 )
    {
        print_error ( "hyperplane::operator== : ERROR : source dimension too long. Not correct hyperplane?\n" );
        abort();
    }
    if ( !plane.points )
    {
        print_error ( "hypeplane::operator== : ERROR : dimensions not equal\n" );
        abort();
    }*/
    int cur_pt=0;	// current point 
    int i=0;
    while ( i<dim )
    {
        if ( cur_pt==dim ) return true;		//if cur_pt==dim, hyperplanes is equal
        if ( points[i]==plane[i] )	//if points is equal
        {
            cur_pt++;			// to next point
            i=0;			// to first point
        }
    }
    return false;
}

/***************************************************************
* not equal operator - compare hyperplanes (use equal operator)*
****************************************************************/
bool hyperplane::operator!= ( const hyperplane& plane ) const
{
    /*if ( plane.dim<2&&plane.dim>6 )
    {
        print_error ( "hyperplane::operator!= : ERROR : source dimension too long. Not correct hyperplane?\n" );
        abort();
    }
    if ( !plane.points )
    {
        print_error ( "hypeplane::operator!= : ERROR : dimensions not equal\n" );
        abort();
    }*/
    if ( *this==plane )
        return false;
    return true;
}

/***************************************************************
*  assignment operator                                         * 
****************************************************************/
hyperplane& hyperplane::operator= ( const hyperplane& plane )
{
    /*if ( plane.dim<2&&plane.dim>6 )
    {
        print_error ( "hyperplane::operator= : ERROR : source dimension too long. Not correct hyperplane?\n" );
        abort();
    }
    if ( !plane.points )
    {
        print_error ( "hypeplane::operator= : ERROR : dimensions not equal\n" );
        abort();
    }*/
    dim=plane.dim;
    if ( !points )
        points=new point[dim];
    for ( int i=0;i<dim;i++ )
        points[i]=plane[i];
    if ( !equation )
        equation=new double[dim+1];
    for ( int i=0;i<dim+1;i++ )
        equation[i]=plane.equation[i];
    outside=plane.outside;
    visited=false;

    if ( !edges )
    {
        edges=new edge[dim];
        for ( int i=0;i<dim;i++ )
            edges[i].points=new point*[dim-1];
    }
    outside=plane.outside;
    for ( int i=0;i<dim;i++ )
    {
        edges[i].neighbor=plane.edges[i].neighbor;
        edges[i].deleted=false;
        for ( int j=0;j<dim;j++ )
            if ( i!=j )
            {
                edges[i].points[j<i?j:j-1]=& ( points[j] );
            }
    }
    return *this;
}

/***************************************************************
*           calculate hyperplane equation coefficients         *
****************************************************************/
void hyperplane::calc_equation()
{
    double** matr = new double*[dim];	// allocate memory
    for ( int i=0; i<dim; i++ )
        matr[i] = new double[dim];

    for ( int i=0; i<dim+1; i++ )
    {
        for ( int y=0;y<dim;y++ )
            for ( int x=0;x<dim;x++ )
            {
                if ( y==i )
                    matr[x][y] = 1;
                else
                    matr[x][y] = points[x][y];
            }
        equation[i]=det ( matr,dim );
    }
    equation[dim]=-equation[dim];
    for ( int i=0; i<dim; i++ )
        delete[] matr[i];
    delete[] matr;
}

/***************************************************************
*           initialize plane                                   *
****************************************************************/
void hyperplane::init_plane ( bool no_equation )
{
    next=prev=0;
    equation=new double[dim+1];	// allocate memory for coefficients
    if ( !no_equation )
        calc_equation();	// calculate equation coefficients
    edges=new edge[dim];
    for ( int i=0;i<dim;i++ )	// make edges
    {
        edges[i].points=new point*[dim-1];	// allocate memory
        edges[i].neighbor=0;			// neighbor is null
        edges[i].deleted=false;
        for ( int j=0;j<dim;j++ )	// add pointers to points to edge
            if ( i!=j )
            {
                edges[i].points[j<i?j:j-1]=& ( points[j] );
            }
    }
}

/***************************************************************
*  distance - calculate distance from hyperplane to point pt   *
****************************************************************/
double hyperplane::distance ( const point& pt ) const
{
    int pt_dim=pt.dimensions();
    /*if ( pt_dim<2&&pt_dim>6 )
    {
        print_error ( "hyperplane::distance : ERROR : source point dimension too long. Not correct point?\n" );
        abort();
    }*/
    double nominator=0;
    double denominator=0;
    for ( int i=0;i<dim;i++ )
    {
        nominator+=pt[i]*equation[i];
        denominator+=equation[i]*equation[i];
    }
    nominator+=equation[dim];
    return nominator/denominator;
	
}

/***************************************************************
*             if point in hyperplane - return true             *
****************************************************************/
bool hyperplane::is_point_in_hyperplane ( const point& pt ) const
{
    int pt_dim=pt.dimensions();
    /*if ( pt_dim<2&&pt_dim>6 )
    {
        print_error ( "hyperplane::is_point_in_hyperplane : ERROR : source point dimension too long. Not correct point?\n" );
        abort();
    }*/
    double res=0;
    for ( int i=0;i<dim;i++ )
    {
        res+=pt[i]*equation[i];
    }
    res+=equation[dim];
    if ( ( res>-delta ) && ( res<delta ) )
        return true;
    //double dist=distance(pt);
    //if ( ( dist>-delta ) && ( dist<delta ) )
        //return true;
    return false;
}

/***************************************************************
*             calculate outside point of hyperplane            *
*    plane - hyperplane of qhull ( for example: first plane )  *
****************************************************************/
void hyperplane::calc_outside ( const hyperplane& plane )
{
    /*if ( plane.dim<2&&plane.dim>6 )
    {
        print_error ( "hyperplane::calc_outside : ERROR : source dimension too long. Not correct hyperplane?\n" );
        abort();
    }
    if ( !plane.points )
    {
        print_error ( "hypeplane::calc_outside : ERROR : dimensions not equal\n" );
        abort();
    }*/
    for ( int i=0;i<dim;i++ )
    {
        if ( !is_point_in_hyperplane ( plane[i] ) ) 
        {
            outside=get_side ( plane[i] ) ?false:true; 
													   	   
            return;
        }
    }
}

/***************************************************************
*             calculate side of point pt                       *
****************************************************************/
bool hyperplane::get_side ( const point& pt ) const
{
    int pt_dim=pt.dimensions();
    /*if ( pt_dim<2&&pt_dim>6 )
    {
        print_error ( "hyperplane::get_side : ERROR : source point dimension too long. Not correct point?\n" );
        abort();
    }*/
    if ( equation[dim]/*<-delta||equation[dim]>delta*/!=0 ) // if plane not contain start of hyperspace
    {
        double dist=distance ( pt );	// get distance
        return dist>0?true:false;	// return side
    }
    // calculate side this last coordinate
    double zn=0;
    for ( int i=0;i<dim-1;i++ )
        zn+=equation[i]*pt[i];
    zn+=equation[dim];
    zn=zn/ ( -equation[dim-1] );
    return ( pt[dim-1]>zn ) ?true:false;	// return side
}

/***************************************************************
*  return true if point pt is outside point of hyperplane      *
****************************************************************/
bool hyperplane::is_outside_point ( point& pt )
{
    int pt_dim=pt.dimensions();
    /*if ( pt_dim<2&&pt_dim>6 )
    {
        print_error ( "hyperplane::is_outside_point : ERROR : source point dimension too long. Not correct point?\n" );
        abort();
    }*/
    if ( is_point_in_hyperplane ( pt ) )
        return true;
    if ( get_side ( pt ) ==outside )
        return true;
    return false;
}

/***************************************************************
*        return poiner to array of edges of plane              *
****************************************************************/
const edge* hyperplane::get_edges() const
{
    /*if ( !edges )
    {
        print_error ( "hyperplane::get_edges : ERROR : getting edges from not initialized hyperplane\n" );
        abort();
    }*/
    return edges;
}

/***************************************************************
*     add plane to neigbors if planes have common edges        *
****************************************************************/
bool hyperplane::add_neighbor ( const hyperplane* plane )
{
    /*if ( !plane )
    {
        print_error ( "hyperplane::add_neightbor : ERROR : neightbor is null?\n" );
        abort();
    }
    if ( plane==this )
    {
        print_error ( "hyperplane::add_neightbor : ERROR : neightbor plane equal self plane?\n" );
        abort();
    }
    if ( plane->dim<2&&plane->dim>6 )
    {
        print_error ( "hyperplane::add_neightbor : ERROR : source dimension too long. Not correct hyperplane?\n" );
        abort();
    }

    if ( plane->points==0 )
    {
        print_error ( "hypeplane::add_neightbor : ERROR : dimensions not equal\n" );
        abort();
    }*/
    //cout<<"add_neighbor:"<<endl;
    //this->print(false);
    //plane->print(false);
    for ( int i=0;i<dim;i++ )	// for each edge
    {
        int cur_pt=0;		// index of current point
        int j=0;
        while ( j<dim )		// for each point in plane
        {
            if ( cur_pt==dim-1 )// if all points of edge have in plane
            {
                if ( edges[i].neighbor )				//if neigbor have
                    edges[i].neighbor->del_neighbor ( this );		// delete self from neigbor neighbors
                edges[i].neighbor=const_cast<hyperplane*> ( plane );	// add neighbor
                //cout<<"Plane added"<<endl;
                return true;						// return true
            }
            if ( * ( edges[i].points[cur_pt] ) == ( *plane ) [j] ) 	// if points is equal
            {
                cur_pt++;	// next point
                j=0;		// to first point of plane
                continue;
            }
            ++j;
        }
    }
    return false;
}


/***************************************************************
*     delete plane from neigbors if it have                    *
****************************************************************/
bool hyperplane::del_neighbor ( const hyperplane* plane )
{
    for ( int i=0;i<dim;i++ )		// for each edges
        if ( plane==edges[i].neighbor )	// if pointers is equal
        {
            edges[i].neighbor=0;	// set this neighbor to null
            return true;		
        }
    return false;
}

/***************************************************************
*     check  membership point of the simplex                   *
****************************************************************/
bool is_point_in_simplex ( point& pt,hyperplane& plane )
{
    //_________________________
    //razmernost' prostranstva
    //La dimension de l'espace
    //Dimension of space
    int n=plane.dimensions();
    //cout<<"plane:"<<plane<<endl;
    //cout<<"pt:"<<pt<<endl;
    //_______
    //vektori
    //vecteur
    //vector
    point* vects=new point[n];
    //______________________________________
    //vichislenieye vectora ,vershina-tochka
    //Le calcul du vecteur ,le sommet-pointчч
    //Calculation of a vector,top-point
    point vect=pt-plane[0];
    for ( int i=1;i<n;i++ )
    {
        vects[i-1]=plane[i]-plane[0];
    }

    //__________________________________
    //matrica - sistema uravnenii
    //La matrice (le système l'équation)
    //Matrix (system the equation)
    double** M=new double*[n-1];
    for ( int i=0; i<n-1; i++ )

        //________________________________
        //vidileniye pamyati
        //La mise en relief de la mémoire
        //Memory allocation
        M[i] = new double[n-1];

    //________________________
    //zapolneniye matrici
    //Le remplissage de la matrice
    //Matrix filling
    for ( int i=0;i<n-1;i++ )
    {
        for ( int j=0;j<n-1;j++ )
            M[i][j]=vects[j][i];
    }
    //__________________________________________________________________________________
    //vidileniye pamyati dlya pravoi chasti uravneniya i korney
    //La mise en relief de la mémoire pour la partie droite de l'équation et les racines
    //Allocation of memory for the right member of equation and roots
    double* x=new double[n-1];
    double* b=new double[n-1];

    //________________________________________________
    //zapolneniye pravoi chasti uravneniya
    //Le remplissage de la partie droite de l'équation
    //Filling of the right member of equation
    for ( int i=0;i<n-1;i++ )
        b[i]=vect[i];

    //_______________________
    //reshenie sistemi
    //La décision du système
    //The system decision
    solve ( ( const double** ) M,b,x,n-1 );
    bool result=true;
    //if(!solve ( ( const double** ) M,b,n-1 ))
      //result=false;
    //for(int i=0;i<n-1;i++)
      //x[i]=b[i];
    double xsum=0;

    //_______________________
    //ochistka pamyati
    //Le nettoyage de la mémoire
    //Memory clearing
    for ( int i=0; i<n-1; i++ )
        delete[] M[i];
    delete[] M;
    delete[] b;

    //_________________________________________
    //predpolagaem, chto tochka v bazise
    //Suppose que le point dans la base
    //Assumes that a point in basis
    

    //_________________________________________
    //slozhenie vsex mnozhitelei vectorov
    //L'addition de tous les multiplicateurs des vecteurs
    //Addition of all multipliers of vectors
    for ( int i=0;i<n-1;i++ )
    { 
        //_________________________________________________
        //esli odin iz nix <0, sostoyanie v false
        //Si un d'eux <0, l'état à false
        //If one of them <0, a condition in false
        if ( x[i]<0 ) result=false;
        if(  isnan(x[i])) result=false; //{ cout<<"x is nan"<<endl; result=false;}else cout<<"x["<<i<<"]="<<x[i]<<endl;
	
        //____________
        //slozheniye
        //addition
        xsum+=x[i];
    }

    //___________________________________________
    //esli summa bol'she odnogo, sostoyanie v false
    //Si la somme est plus grande qu'un, l'état à false
    //If the sum more than one, a condition in false
    if ( xsum>1/*-0.0001*/ ) result=false;

    //______________
    //ochistka pamyati
    //Le nettoyage de la mémoire
    //Memory clearing
    delete[] x;
    delete[] vects;
    return result;
}


/***************************************************************
*     convert all points to concentration_mole format          *
****************************************************************/
hyperplane& hyperplane::to_conc_mole_format()
{
    for ( int i=0;i<dim;i++ )		// for each point
        points[i].to_conc_mole_format();// convert point
}

ostream& operator<< ( std::ostream &os, const hyperplane &hp )
{
    os<<"{";
    for ( int i=0;i<hp.dimensions();i++ )
    {
        os<<hp[i];
        if ( i!=hp.dimensions()-1 )
            os<<",";
    }
    os<<"("<< ( hp.get_outside() ?"under":"after" ) <<")";
    os<<"}";
    return os;
}
