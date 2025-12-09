#include "hyperplane.h"
#include "phase.hpp"

#include <vector>

using namespace std;

class qhull;
ostream& operator<< ( ostream &os, const qhull &qh );

class qhull
{
    int dim;			// dimension
    //int discritization;
    int num_points;		// number of points in hull
    int num_planes;		// number of planes in hull
    point_list source_points;	// list of sourse points
    point_list hull_points;	// list of hull points
    hyperplane* firstplane;	// firstplane
    hyperplane* start;		// pointer to start of list of hyperplanes
    hyperplane* end;		// pointer to end of list of hyperplanes
    hyperplane* add_plane ( hyperplane* plane );			// add plane to end of hull
    hyperplane* add_plane ( hyperplane* plane,hyperplane* prev );	// end plane after prev 
    hyperplane* del_plane ( hyperplane* plane );			// delete plane from hull
    point* next_point ( hyperplane& plane );				// get next point
    void set_neighbors ();		// set neighbors to new heperplanes
    hyperplane* add_point ( hyperplane* plane,point* point );		// add point to hull
    hyperplane* find_visiblies ( hyperplane* plane,point* point );	// find all visible planes
    void remove_inner_points();


public:
    qhull ( int d );				// constructor
    ~qhull();					// destructor
    int size();					// get size of hull
    void build_hull();				// build convex hull
    bool is_correct();				// if hull is correct - return true
    bool is_contain_all_points();		// if hull contain all points - return trie
    //	discretozation
    void discretization ( double T,int n );
    void discretization2d ( double T,int n );
    void discretization3d ( double T,int n );
    void discretization4d ( double T,int n );
    void discretization5d ( double T,int n );
    void discretization6d ( double T,int n );
    // find extremums
    vector<point> extremums(double cr=2.5e-3,double mo=1.4e-2);
    vector<point> extremums2d();
    vector<point> extremums3d(double cr);
    vector<point> extremums4d(double cr,double mo);
    vector<point> extremums5d();
    vector<point> extremums6d();
    void print() const;
    /*void print_list()
    {
        printf("start=%x\n",start);
        hyperplane* cur=start;
        while (cur)
        {
            cur->print_ptrs();
            cur->print(true);
            cur=cur->next;

        }
        printf("\nend=%x\n",end);
    }*/
    void write_tecplot_file();
    interface *func1;
    interface *func2;

    interface *Jacob1;
    interface *Jacob2;

    double ( interface::*nrfuncv1 ) ( double c,double fe,double cr, double mo, double co, double al );
    double ( interface::*nrfuncv2 ) ( double c,double fe,double cr, double mo, double co, double al );

    std::vector < std::vector<double> > ( interface::*nrJacob1 ) ( double c,double fe,double cr, double mo, double co, double al );
    std::vector < std::vector<double> > ( interface::*nrJacob2 ) ( double c,double fe,double cr, double mo, double co, double al );

};
