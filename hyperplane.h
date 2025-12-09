#include "point.h"
#include <vector>

using namespace std;
class hyperplane;

bool is_point_in_simplex ( point& pt,hyperplane& plane );

// edge
//Rus:gran'
struct edge
{
    hyperplane* neighbor;	// neighbor for this edge
							//Rus:smezhnaya giperploskost' po etoi grani
    point** points;			// pointers to points
							//Rus:ssilka na tochki
							//Example:Int mas[10];
							//dinamique:int* mas=new int[10];
							//int** mas=new int*[10]; - kazhdii element massiva ukazatel', a ne chislo.
							//Ukazatel'(int*) na chislo tipa int. Eta massiv ukazatelei
							//Pervaya zvezdochka, potomu chto v massive xranyatsya ukazateli
							//Vtaraya *, potomu chto beretsya ne sam massiv , a ukazatel' na nego.
	
    bool deleted;			//flag. set true if this is horyzont edge
							//Rus: esli eto liniya gorizonta, znachit istina
};

class hyperplane
{
    friend class qhull;
private:
    hyperplane* next;	// pointer to next plane in list
						// Rus: ukazatel' na sledushuu ploskost' v liste
    hyperplane* prev;	// pointer to prevous plane in list
						// Rus: ukazatel' ne predidushuyu ploskost' v liste
    int dim;			// dimension
    bool outside;		// outside
						// Rus:vnechnya storona	
    bool visited;		// if plane visited, this flag is true.
						/* Rus:Poisk vidimix ploskostei proisxodit po sosedyam
						Popadanie na odnu i tu zhe ploskost' cherez sosedei
						razlichnimi putyami privedet k zatciklivaniu.
						chtobi ne prosmatrivat' uzhe prosmotrennie ploskosti 
						proveryaetsya etot flag. vse giperploskoti, u kotorix
						imeetsya etot flag, yavlyautsya tak zhe vidimimi i udalyautsya
						pri dobavlenii tochki */
    bool old;			// if it is plane is new - false. (for optimization)
						/* Rus:etot flag raven false dlya giperploskostei, dobavlennix
						 na etape dobavleniya tekushei tochki. Dlya ostal'nix giperploskostei true.
						 Ustanavlivaetsya etot flag v qhull::set_neighbors, nuzhen dlya togo, chtoby
						 iskat' sosedei tol'ko dlya poslednix dobavlennix giperploskostei sredi nix zhe.*/
    bool border;		// set true if this plane is horyzont plane
						// Rus:istina, esli ploskost' vxodit v ensemble ploskostei, formiruushix liniu gorizonta
    bool no_points;		/* Rus:eto tozhe dlya optimizatcii. Esli dlya giperploskosti net vneshnix tochek, to
						 ustanavlivaetsya etot flag. Eto, chtobi propuskat' takie giperploskosti pri postroenii
						 */
    point* points;		// pointer to points (array of points)
						/* Rus:zdes' massiv tochek, obrashenie k kotoromu proisxodit cherez ukazatel'. Zdes'
						 elementom massiva yavlyaetsya tochka, a vishe elementom massiva yavlyaetsya ukazatel'
						 na tochku, kotoraya naxoditsya v etom massive.
						 */
    double* equation;	// pointer to equation (coefficients of hyperplane equation)
						// Rus: ssilka na koefficienti uravneniya ploskosti	
    edge* edges;		// pointer to edges of plane
						/* Rus:ssilka na grani ploskosti. Eto massiv granei (obektov edge). Kol-vo
						   granei ravno izmereniu prostranstva*/	
    void calc_equation();	// calculate equation coefficients
							// raschet koeff-ov uravneniya ploskosti
    void init_plane ( bool no_equation );	//initialize plane
											//nachal'naya ploskost' esli no_equation=true
    hyperplane();		// default constructor
public:
    hyperplane ( hyperplane const& plane );	//copy constructor
    hyperplane ( point& pt1,point& pt2,bool no_equation=false );											//2d
																			//no_equation=true, uravnenie ne rasschitivaetsya
																			//v ekstremumax takoe vozmozhno
    hyperplane ( point& pt1,point& pt2,point& pt3,bool no_equation=false );									//3d
    hyperplane ( point& pt1,point& pt2,point& pt3,point& pt4,bool no_equation=false );						//4d
    hyperplane ( point& pt1,point& pt2,point& pt3,point& pt4,point& pt5,bool no_equation=false );			//5d
    hyperplane ( point& pt1,point& pt2,point& pt3,point& pt4,point& pt5,point& pt6,bool no_equation=false );//6d
    hyperplane ( point& pt1,point** pts,int dim );				// build plane from point and array of dim-1 points
																// Rus:sozdanie ploskosti s razmernost'u dim-1 dlya poiska ekstremumov			
    ~hyperplane();	//destructor
    bool operator== ( const hyperplane& plane ) const;	//equal operator
    bool operator!= ( const hyperplane& plane ) const;	//not equal operator
    hyperplane& operator= ( const hyperplane& plane );	//assignment operator
    double distance ( const point& pt ) const;		// calculate distance to point
													// Rus:raschet distancii mezhdu tochkoi i ploskost'u
    void calc_outside ( const hyperplane& plane );	// calculate outside of plane
													// Rus:functciya ustanavlivaet vneshnuu storonu po svoemu parametru
    void set_outside ( bool side )	// set outside of plane
    {								// Rus:functciya ustanavlivaet vneshnuu storonu po svoemu parametru
        outside=side;
    }
    bool get_outside() const		// get outside of plane
    {									
        return outside;
    }
    bool is_point_in_hyperplane ( const point& pt ) const;	// if point in hyperplane - return true
															// Rus: esli tochka v ploskosti - return true				
    point& operator[] ( int i ) const;						// iterator
    bool get_side ( const point& pt ) const;				// return side of plane
															// Rus: vozvrat giperploskosti	
    const edge* get_edges() const;							// return pointer to array of edges
															// Rus:vozvrat ukazatelya na gran'
    bool add_neighbor ( const hyperplane* plane );			// add neighbor to plane
															// Rus:Functciya vipolnyaet sravnenie tochek i ishet
															// est' li u giperploskosti gran', po kotoroi ona smezhna
															// s ploskost'u. Esli takoi grani net,to nikuda ne dobavlyaetsya sosed
    bool del_neighbor ( const hyperplane* plane );			// del neighbor from plane
															// Rus: udalyaem odnogo soseda, ukazatel' na kotorii peredaetsya suda.	
    hyperplane& to_conc_mole_format();						// convert points to CONCENTRATION_MOLE format
    bool is_outside_point ( point& pt );					// return true if pt is outside point
															
    bool is_visited()			// is neighour hyperplane visited
    {								
        return visited;
    }
    void set_visited()			// set neighbour hyperplane as visited
    {								
        visited=true;
    }
    bool is_old()				// is hyperplane in the convex envelope is old
    {							
        return old;
    }
    bool set_old()				// set hyperplane in the convexe envelope is old
    {							
        old=true;
    }
    bool is_deleted_neighbor ( int i )	// set neighbor i as horyzont
    {									// Rus: ustanavlivaut i-yu gran kak gran' gorizonta
        return edges[i].deleted;
    }
    bool set_deleted_neighbors ( hyperplane* plane )	// set hyperplane plane as non visible neighbor
    {													// Rus: sredi sosedei ishetsya sosed ploskosti.
														// idet sravnenie ukazateleli, sravneniya tochek zdes' net
        for ( int i=0;i<dim;i++ )						// for each edges
														// Rus: dlya kazhdoi grani
            if ( edges[i].neighbor==plane )				// if neighbor for this edge equal plane
            {											// Rus:esli u tekushei grani sosed plane, to ustanavlivaetsya
														// to ustanavlivaem etu gran' kak gorizont	
                edges[i].deleted=true;					// set this edge as horyzont
                return true;							
            }
        return false;
    }
    bool is_border()		// if plane is border - return true
    {						
        return border;
    }
    void set_border ( int i )	// set this plane as horyzont
    {							
        border=true;
        edges[i].deleted=true;	//Rus: udalenie vsex vidimix ploskostei	
    }
    
    bool is_no_points(){return no_points;}
    void set_no_points(){no_points=true;}
    int dimensions() const	// return dimension of hyperplane
    {
        return dim;
    }
    const double* get_coeffs()	// return pointer to array of hyperplane eqution coefficients
    {
        return equation;
    }
    /*void print_ptrs()
    {
      printf("(%x <= %x => %x)\n",prev,this,next);
      for(int i=0;i<dim;i++)
    printf(" %x",edges[i].neighbor);
      printf("\n");

    }
    /*
    void print(bool all)const
    {
      std::printf("(");
      if(points)
      {
    for(int i=0;i<dim;i++)
    {
      points[i].print();
      if(i!=dim-1) std::printf(",");
    }
    std::printf(" (%d) )\n",dim);
      }
      if(!all) return;
      printf("equation:");
      for(int i=0;i<dim+1;i++)
      {
    printf("%f",equation[i]);
    if(i!=dim) std::printf(",");
      }
      //printf("\n");
      if(edges)
    for(int i=0;i<dim;i++)
    {
      printf("edges[%d]:(",i);
      for(int j=0;j<dim-1;j++)
      {
        edges[i].points[j]->print();
        if(j!=dim-2) std::printf(",");
      }
      printf(") =>");
      if(edges[i].neighbor)
        edges[i].neighbor->print(false);
      else printf("none");
        //printf("\n");
    }
      if(next)
      {
    printf("next:");
    next->print(false);
      }
      if(prev)
      {
    printf("\nprev:");
    prev->print(false);
    printf("\n");
      }
    }
    */
};

ostream& operator<< ( ostream &os, const hyperplane &hp );
