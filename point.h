#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

#define print_error printf

//The other kind creates an unnamed type. 
//This is used when you want names for 
//constants but don't plan to use the type
//to declare variables, function arguments, etc. For example, you can write 
enum
{
    x=0,y,z,t,w,u
};

//Class point
class point
{
    friend class qhull;
    friend class point_list;
private:
	
	//The previous and next values (point) are established
	//in correct positions at addition and removal of point.
    point* next;		// The next point (for list)
    point* prev;		// The previous point (for list)


	int function;		// Defines an accessory of a point of one of two function
    double* coords;		// Vector with coordinates
    int dim;			// dimension of space
public:

    point();		    //default constructor
	
	//It is used at creation of a point from another points
	//example:
	//point pt=point(x,y,func); // creation of point
	//point pt1=point(pt);      // constructor of the copy

    point ( point const& pt );					// copy constructor
    point ( double* coord,int dim,int func );	// build point from array
    point ( double x,double y,int func );										//2d
    point ( double x,double y,double z,int func );								//3d
    point ( double x,double y,double z,double t,int func );						//4d
    point ( double x,double y,double z,double t,double w,int func );			//5d
    point ( double x,double y,double z,double t,double w,double u,int func );	//6d
    ~point();								// destructor
    double operator[] ( int i ) const;		// iterator
    bool operator== ( const point& pt );	// equal operator
    bool operator!= ( const point& pt );	// not equal operator
    point& operator= ( point const& pt );	// assignment operator
    point operator- ( point& pt );			// substract operator
    int dimensions() const
    {
        return dim;							// return dimension of point
    }
    int func() const
    {
        return function;					// return function of point
    }
    point& to_conc_mole_format();			// convert point to CONCENTRATION_MOLE format
    point& from_conc_mole_format();			// convert point from CONCENTRATION_MOLE format
};

// class container for list of points. 
//Points are organised in the list with use point* next and point* prev.
//On the basis of this class point_list in a class qhull there are two fields variables a class qhull:
//							-- source_points - Includes points after descritization which be deleted 
//											   at addition in a convex qhull. Points from the list 
//											   which lie in a convex qhull also be deleted.
//							-- - hull_points -Contains the points which have entered into 
//											   a convex qhull. They are added here in process 
//											   of addition of points in a convex qhull.
class point_list
{
    friend class qhull;
    point* begin;			// first element pointer -
							//The pointer to value of the first element with which
							//have filled page (list) of points
    point* end;				// end element pointer -
							//The pointer to value of last element with which have 
							//filled sheet of points	
public:
    point_list()			// constructor
    {
        begin=end=0;
    }
    ~point_list()			// destructor (free memory)
    {
        point* cur=begin;	//We establish the pointer on the beginning of the list of points
        while ( cur!=0 )	//While the pointer on a point not is a zero we do following actions
        {
            point* tmp=cur; //We save value of the pointer in the temporary variable
            cur=cur->next;	//We fill the pointer of a current point by the following
							//point in the list. Here, if the pointer on the following
							//point is equal to zero (cur==0) there is the exit from a cycle.
            delete tmp;		//Clearing of the memory in which the pointer of the "tmp" pointe
							//out (a current point of the list).
        }
    }
    void add ( point& pt )					// add point to list
    {
        point* new_point=new point ( pt );  //Allocation of the memory under the new point
        if ( begin==0 )						//If begin == 0 that the list is empty.											
        {									//I.e. the added element will be the first and last (unique) in the list
            *new_point=pt;					//We assign to a point with the address new_point point pt
            new_point->next=0;				//Pointer on following and previous points are established in a zero
            new_point->prev=0;
            begin=end=new_point;
        }
        else								//If the list is not empty, i.e. in it already there are points
        {
            end->next=new_point;			//At last element of the list the pointer on the following 
											//element it is equal to zero. We establish it in an added point.
											//Rus:u poslednego elementa spiska ukazatel' na sledushii element
											//raven nulu.Ustanavlivaem ukazatel' na dobavlyaemuu tochku.

			
            new_point->prev=end;			//We establish the pointer on last point of a list (last point 
											//still has not changed)
											//Rus:ustanavlivaem ukazatel' na predidushii element dobavlyaemoi
											//tochki, na poslednii element obolochki (poslednii element eshe 
											//ne izmenilsya)

            end=new_point;
        }
    }
	
	//The point delete on which the pointer is transferred in function delete.
	//The point can be the first, last or in the middle of list with points
    point* del ( point* pt )				// delete point from list
    {
        point* ret=pt->next;				//The pointer on the following point after a deleted point 
											//is brought in ret
        if ( !pt )							//If in function it is transferred the pointer zero
        {
            print_error ( "point_list::del: ERROR : pointer in deleted point equal zero" );
            abort();
        }
        if ( pt==begin&&pt==end )			//If this point is unique in list (the first and last)
        {
            begin=end=0;					//We null list
        }
        else
        {
            if ( begin==pt )				//If this point is the first point in list
            {
                begin=pt->next;				//The pointer of the beginning of the list is point
											//out to a point, following the deleted point
											//Rus: stavitsya ukazatel' nachala spiska na tochku, 
											//sledushuyu za udalyaemoi
                pt->next->prev=0;			//The pointer of this following point is point
											//out to the previous point (it was the deleted point)
											//and point out in zero.
											//Rus:ukazatel' etoi tochki na predidushii element
											//(im bila udalyaemaya tochka) ustanavlivaetsya v nol'
            }
            else if ( end==pt )				//If a point in the end of list
            {
                pt->prev->next=0;			//The pointer on the following point from point1 is
											//established in a zero. Point1 is previous a deleted point.
											//Rus:ukazatel' na sledushuyu tochku ot tochki; kotorya 
											//yavlyaetsya predidushej udalyaemoi.Ukazatel' ustanavlivaem 
											//v pozitciu nol'
                end=pt->prev;				//The pointer of the end of the list is established on a 
											//point is previous by the deleted.
											//Rus:ukazatel' kontca spiska ustanavlivaetsya na tochku,
											// predidushej udalyaemoi            
			}
            else
            {
                pt->prev->next=pt->next;	//The pointer from the previous point of a deleted point
											//is established on the following point after the deleted. 
											//[1,...,2] - > [1->2]	
                pt->next->prev=pt->prev;	//The pointer from the following point after the deleted 
											//point is established on the previous point before the
											//deleted point [1,...,2] -> [1<-2]
											//Rus: ukazatel' ot sleduyushhei tochki posle udalennoi
											//ustanavlivaetsya na predidushuyu tochku pered udalennoj
											//i obratnyj protcess
            }
        }
        delete pt; //Delete a point after definition her place in list
        return ret;
    }
    bool empty()	// return true if list is empty
    {
        return begin?false:true;
    }
    int size()				// return size of list
    {
        int size=0;
        point* cur=begin;
        while ( cur )		//While the list is not empty
        {
            size++;			//Increase the size
            cur=cur->next;	//Transition to the following point
        }
        return size;
    }
};

// output operator for point
ostream& operator<< ( ostream &os, const point &pt );
