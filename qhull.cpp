#include "qhull.h"
#include "tools.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;



extern std::vector <double> CONCENTRATION_MASStoCONCENTRATION_MOLE ( std::vector <double> &C,std::vector <double> &M );
extern double x_i;
extern double x_j;
/***************************************************************
*  constructor (set dim and all members of hyperplane to null) *
****************************************************************/
qhull::qhull(int dim)
{
    this->dim=dim;
    start=end=0;
    firstplane=0;
    num_points=dim;
    num_planes=0;
    func1=func2=0;
    Jacob1=Jacob2=0;
}

/***************************************************************
*                destructor - free memory                      *
****************************************************************/
qhull::~qhull()
{
    hyperplane* cur=start;
    while ( cur )
    {
        hyperplane* tmp=cur;
        cur=cur->next;
        delete tmp;
    }
    if (firstplane)
        delete firstplane;
    
}

/***************************************************************
*                 add plane to end of hull                     *
****************************************************************/
hyperplane* qhull::add_plane ( hyperplane* plane )
{
    hyperplane* new_plane=new hyperplane();	// create new plane
    *new_plane=*plane;				
    if ( !start )		// if hull is empty
    {
        start=new_plane;	// set start to new plane
        new_plane->next=0;	
        new_plane->prev=0;
    }
    else
    {
        end->next=new_plane;	// set old end next plane to new plane
        new_plane->prev=end;
        new_plane->next=0;
    }
    end=new_plane;		// set end to new plane
    num_planes++;		// increase number of planes
    return new_plane;
}

/***************************************************************
*                   add plane after prev                       *
****************************************************************/
hyperplane* qhull::add_plane ( hyperplane* plane,hyperplane* prev )
{
    hyperplane* new_plane=new hyperplane(); //create new plane
    *new_plane=*plane;
    if ( !start )		// hull is empty
    {
        start=new_plane;	// 
        end=new_plane;
        new_plane->next=0;
        new_plane->prev=0;
    }
    else
    {
        new_plane->next=prev->next;
        new_plane->prev=prev;
        prev->next=new_plane;
        if ( new_plane->next )
            new_plane->next->prev=new_plane;

        if ( prev==end )
            end=new_plane;
    }
    num_planes++;
    return new_plane;
}

/***************************************************************
*                   delete plane from hull                     *
****************************************************************/
hyperplane* qhull::del_plane ( hyperplane* plane )
{
    hyperplane* ret=plane->next;
    if ( plane==start&&plane==end )
    {
        start=end=0;
    }
    else
    {
        if ( plane==start )
        {
            start=plane->next;
            plane->next->prev=0;
        }
        else if ( plane==end )
        {
            end=plane->prev;
            plane->prev->next=0;
        }
        else
        {
            plane->prev->next=plane->next;
            plane->next->prev=plane->prev;
        }
    }
    const edge* edges=plane->get_edges();
    for ( int i=0;i<dim;i++ )
        if ( edges[i].neighbor )
        {
            edges[i].neighbor->del_neighbor ( plane );
        }
    delete plane;
    num_planes--;
    return ret;
}

/***************************************************************
*            find next point to add to hull                    *
****************************************************************/
point* qhull::next_point ( hyperplane& plane )
{
    double max_dist=0;
    point* max_pt=0;
    for ( point* cur=source_points.begin;cur!=0;cur=cur->next )	// for each points
    {
        if ( plane.is_outside_point ( *cur ) )		 	// if point is outside
	{
	    double dist=fabs ( plane.distance ( *cur ) );	//calc distance
	    if ( dist>max_dist )				// if distance > max distance
	    {
		max_dist=dist;					// update max daistance
		max_pt=cur;					// save pointer to point with max distance
	    }
	}
    }	
    return max_pt;	// return pointer to point
}

/***************************************************************
*            find next point to add to hull                    *
****************************************************************/
void qhull::build_hull()
{
	
    int no_points=0;
    hyperplane* cur=start;
    while ( true )	// for ever
    {
	if(cur->is_no_points())
	{
	  no_points++;		// 
            if(no_points==num_planes)	//
	      break;
            cur=cur->next;		// to next plane
            if ( !cur )			// if it is end of hull
	    {
	      no_points=0;		//
	      cur=start;		// set current hyperplane to first
	    }
            continue;
	}
        point* pt=next_point ( *cur );	// get point
        if ( !pt )			// if no point
        {
	    cur->set_no_points();
            no_points++;		// 
            if(no_points==num_planes)	//
	      break;
            cur=cur->next;		// to next plane
            if ( !cur )			// if it is end of hull
	    {
	      no_points=0;		//
	      cur=start;		// set current hyperplane to first
	    }
            continue;
        }
        //cout<<"Add point:"<<*pt<<" for plane:"<<*cur<<endl;
        cur=add_point ( cur,pt );	// add point to hull
	//write_tecplot_file();
	if(num_planes<100)
	  remove_inner_points();
        if ( !cur )			// if it is end of hull
	{
	    no_points=0;		//
            cur=start;			// set current hyperplane to first
	}
    }
//write_tecplot_file();
    //cout<<"Planes:"<<num_planes<<endl;
    //cout<<"Points:"<<num_points<<endl;
    /*if(!is_contain_all_points()) 		// if hull not contain all source points
    {
	cout<<"Hull not contain all points"<<endl;
        cout<<"Number of points:"<<num_points<<endl;
        cout<<"Number of planes:"<<num_planes<<endl;
	write_tecplot_file();
        abort();				// abort
    }*/
}

/***************************************************************
*                      add point to to hull                    *
*      plane - plane, to each add hull                         *
*      pt - point to add                                       *
****************************************************************/
hyperplane* qhull::add_point ( hyperplane* plane, point* pt )
{
    hyperplane* from;
    hyperplane* ret=plane;
    find_visiblies ( plane,pt);	// find all visible planes
    hyperplane* cur=start;	// set current hyperplane to start
    while ( cur )		// for each planes
    {
	if ( cur->is_visited() ) // if plane visited( visible)
	{
	  if ( cur->is_border() )//if plane contain horyzont edge
	  {
            const edge* edges=cur->get_edges(); // get edges of plane
            for ( int i=0;i<dim;i++ )	// for each edge
            {
                if ( edges[i].deleted ) // if this edge is horyzont
                {
                    hyperplane new_plane=hyperplane ( *pt,edges[i].points,dim );	// make new plane
                    new_plane.calc_outside ( *firstplane );				// calculate outside for new plane
                    ret=add_plane ( &new_plane,ret->next?ret->next:ret );		// add plane to hull
		    if(edges[i].neighbor)	// set horyzont neighbors
		    {
		      ret->add_neighbor(edges[i].neighbor);	// add neighbor for new plane
		      edges[i].neighbor->add_neighbor(ret);	// add neighbor to neighbor neighbor
		    }
		    //cout<<"Add plane:"<<*ret<<endl;
                }
            }
	  }
	  //cout<<"Del plane:"<<*cur<<endl;
	  //cur=del_plane ( cur );
	}

        /*else*/ cur=cur->next;
    }
    cur=start;
    while(cur)	// foreach planes
      if(cur->is_visited())	// if plane is visited (visible)
      {
	//cout<<"Del plane:"<<*cur<<endl;
	cur=del_plane(cur);	// delete plane
      }
      else cur=cur->next;	//else - to next plane
    num_points++;		// increase number of points
    //cout<<"points:"<<num_points<<endl;
    //cout<<"planes:"<<num_planes<<endl;
    set_neighbors ();	// correct neighbors
    //print_list();
    //print();
    hull_points.add ( *pt );	// add point to list of hull points
    source_points.del ( pt );	// delete point from list of source points
    //cout<<*this<<endl;
    /*if ( !is_correct() )	// if hull not correct
    {
        cout<<"Hull not correct"<<endl;
        cout<<"Number of points:"<<num_points<<endl;
        cout<<"Number of planes:"<<num_planes<<endl;
	write_tecplot_file();
        abort();		// abort
    }*/
    return ret->next;		// return next plane for add next point
}

/***************************************************************
*                      add point to to hull                    *
*      plane - plane, to each add hull                         *
*      pt - point to add                                       *
****************************************************************/
hyperplane* qhull::find_visiblies ( hyperplane* plane, point* pt)
{
    plane->set_visited();			// set plane visited
    const edge* const edges=plane->get_edges();	// get edges
    for ( int i=0;i<dim;i++ )		// for each edges
    {
        if ( edges[i].neighbor )	// if have neighbor for this edge
        {
            if ( edges[i].neighbor->is_visited() )	// if neighbor alteady visited
                continue;				//continue
            if ( edges[i].neighbor->is_outside_point ( *pt ) )	// if point is outside point for neighbor
                find_visiblies ( edges[i].neighbor,pt);		//recurse call for neighbor
            else plane->set_border ( i );			// else - this is horyzont
        }
        else plane->set_border ( i );	// else this is zoryzont
    }
    return 0;
}

/***************************************************************
*                      correct neighbors                       *
****************************************************************/
void qhull::set_neighbors ()
{
    hyperplane* begin=start;
    for ( hyperplane* cur=begin;cur;cur=cur->next ) // for each planes
    {
	if(cur->is_old())	// if plane is old
            continue;		// continue
        int counter=0;		
        for ( hyperplane* cur1=cur->next;cur1;cur1=cur1->next )	// foreach planes
        {
            if(counter==dim)	// if counter == dim, all neighbors is added
	      break;
            if(cur1->is_old())  // if plane is old
		continue;	// continue
            if ( cur->add_neighbor ( cur1 ) )	// add neighbor to plane
            {
                counter++;
                if ( !cur1->add_neighbor ( cur ) )// add neighbor to neigbor
                    cout<<"qhull::set_neighbors : ERROR : Not add neighbor to neighbor"<<endl;
            }
        }
        cur->set_old();	// set plane as old
    }
}

/***************************************************************
*                    return size of plane                      *
****************************************************************/
inline int qhull::size()
{
    return num_planes;
}

/***************************************************************
*                         print hull                           *
****************************************************************/
void qhull::print() const
{
    for ( hyperplane* cur=start;cur!=0;cur=cur->next )
        cout<<*cur<<endl;
}

/***************************************************************
*          write hull to file for tecplot                      *
****************************************************************/
void qhull::write_tecplot_file()
{
    fstream file;
    char filename[30];
    //memset ( filename,0,30 );
    int num_points=hull_points.size();
    sprintf ( filename,"tecplot/%d.txt",num_points );
    //if(!file9.is_open())
    file.open ( filename,fstream::out );
    if ( file.fail() )
        cout<<"Error opening file"<<endl;
    file<<"TITLE = \"Convex hull 3D for "<<num_points-dim<<" points\""<<endl;
    file<<"VARIABLES = \"X\", \"Y\", \"Z\""<<endl;
    file<<"ZONE NODES="<<num_points <<", ELEMENTS="<<size() <<", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE"<<endl<<endl;
	//FETRIANGLE  FEQUADRILATERAL
	
	for ( int i=0;i<dim;i++ )
	{
		point* cur=hull_points.begin;
		while ( cur )
		{
			file<< ( *cur ) [i]<<" ";
			cur=cur->next;
		}
		file<<endl;
	}//
	point* cur=hull_points.begin;
	while ( cur )
	{
		file<< ( *cur ).func()<<" ";
		cur=cur->next;
	}
	file<<endl;//
    hyperplane* cur_plane=start;
    while ( cur_plane )
    {
        for ( int j=0;j<dim;j++ )
        {
            point* cur_point=hull_points.begin;
            for ( int k=0;cur_point;k++ )
            {
                if ( ( *cur_point ) == ( *cur_plane ) [j] )
                    file<<k+1<<" ";
                cur_point=cur_point->next;
            }
        }
        file<<endl;
        cur_plane=cur_plane->next;
    }
    file.close();
}

/***************************************************************
*                    check hull for correct                    *
****************************************************************/
bool qhull::is_correct()
{
    hyperplane* cur=start;
    int no_correct=0;
    while ( cur )
    {
        point* cur_pt=hull_points.begin;
        while ( cur_pt )
        {
            if ( cur->is_outside_point ( *cur_pt ) && !cur->is_point_in_hyperplane(*cur_pt))
            {//return false;
                cout<<"hull not convexed"<<endl;
                cout<<*cur<<endl;
                cout<<*cur_pt<<endl;
		no_correct++;
            }
            cur_pt=cur_pt->next;
        }

        const edge* edges=cur->get_edges();
        int neighbors=0;
        for ( int i=0;i<dim;i++ )
        {
            if ( edges[i].neighbor )
                neighbors++;
            else
            {
                int equ=0;
                for ( int j=0;j<dim;j++ )
                {
                    if ( equ==dim-1 )
                    {
                        neighbors++;
                        break;
                    }
                    if ( * ( edges[i].points[equ] ) == ( *firstplane ) [j] )
                    {
                        equ++;
                        j=0;
                    }
                }
            }
        }
        if ( neighbors!=dim )
        {
            cout<<"neighbor not correct"<<endl;
            cout<<*cur<<endl;
            no_correct++;
        }
        cur=cur->next;
    }
    if ( no_correct )
    {
        //print_list();
        return false;
    }
    return true;
}

/***************************************************************
*       check hull for contain all source points               *
****************************************************************/
bool qhull::is_contain_all_points()
{
    hyperplane* cur=start;
    int no_correct=0;
    while ( cur )
    {
        point* cur_pt=source_points.begin;
        while ( cur_pt )
        {
            if ( cur->is_outside_point ( *cur_pt ) && !cur->is_point_in_hyperplane(*cur_pt))
            {
                cout<<"hull not contain point"<<endl;
                cout<<*cur<<endl;
                cout<<*cur_pt<<endl;
		no_correct++;
            }
            cur_pt=cur_pt->next;
        }
        cur=cur->next;
    }
    if ( no_correct )
    {
        return false;
    }
    return true;
}

void qhull::remove_inner_points()
{
  point* cur_pt=source_points.begin;
  while(cur_pt)
  {
    hyperplane* cur=start;
    int cnt=0;
    while(cur)
    {
      if(!cur->is_outside_point(*cur_pt))
	cnt++;
      cur=cur->next;
    }
    if(cnt==num_planes)
      cur_pt=source_points.del(cur_pt);
    else cur_pt=cur_pt->next;
  }
}

void qhull::discretization ( double T, int n )
{
    switch (dim)
    {
    case 2:
        discretization2d(T,n);
        break;
    case 3:
        discretization3d(T,n);
        break;
    case 4:
        discretization4d(T,n);
        break;
    case 5:
        discretization5d(T,n);
        break;
    case 6:
        discretization6d(T,n);
        break;
    }
    add_plane(firstplane);
    point* cur=source_points.begin;
    fstream ofile;
    char filename[30];
    //memset ( filename,0,30 );
    int num_points=hull_points.size();
    sprintf ( filename,"points-%d.txt",dim );
    //if(!file9.is_open())
    ofile.open ( filename,fstream::out );
    if ( ofile.fail() )
        cout<<"Error opening file"<<endl;
    else
    {
      cur=source_points.begin;
      while(cur)
      {
	for(int i=0;i<dim;i++)
	  ofile<<(*cur)[i]<<" ";
	ofile<<endl;
	cur=cur->next;
      }
      ofile.close();
    }
    fstream func1("function1.txt",fstream::out);
    fstream func2("function2.txt",fstream::out);
    cur=source_points.begin;
    while(cur)
    {
      if((*cur).func()==0)
      {
	for(int i=0;i<dim;i++)
	  func1<<(*cur)[i]<<" ";
	func1<<endl;
      }
      else
      {
	for(int i=0;i<dim;i++)
	  func2<<(*cur)[i]<<" ";
	func2<<endl;
      }
      if (!firstplane->is_outside_point ( *cur ) && !firstplane->is_point_in_hyperplane(*cur))
	cur=source_points.del(cur);
      else cur=cur->next;
    }
    func1.close();
    func2.close();
}
void qhull::discretization2d (double T,int n)
{
    double step=1./n;
    point pt1,pt2;
    for (int i=0;i<=n;i++)
    {
		
		/*double G[2];
		G[0]= ( func1->*nrfuncv1 ) ( i*step,0,0,0,0,0 );
		G[1]= ( func2->*nrfuncv2 ) ( i*step,0,0,0,0,0 );
		double min=G[0]>G[1]?G[1]:G[0];
		int func=G[0]>G[1]?1:0;
		point pt=point(i*step,min,func);
		pt.to_conc_mole_format();
		pt=point(pt[0],pt[1],min,func);
		printf("%e	%lf \n",min,pt[0]);
		if (i==0) pt1=pt;
        else if (i==n) pt2=pt;
        else source_points.add ( pt );*/
		
		
		
        double G[2];
        G[0]= ( func1->*nrfuncv1 ) ( i*step,0,0,0,0,0 );		
		//printf("%e	%e	\n",G[0],i*step);
		//printf("%e	%e	%lf %d \n",G[0],G[1],i*step,func);
        G[1]= ( func2->*nrfuncv2 ) ( i*step,0,0,0,0,0 );
		//printf("%e	%e	\n",G[1],i*step);

        double min=G[0]>G[1]?G[1]:G[0];
        int func=G[0]>G[1]?1:0;
		point pt=point(i*step,min,func);
        if (i==0) pt1=pt;
        else if (i==n) pt2=pt;
        else source_points.add ( pt );
		 
    }
    hull_points.add(pt1);
    hull_points.add(pt2);
    firstplane=new hyperplane(pt1,pt2);
    firstplane->set_outside(true);
}
void qhull::discretization3d (double T,int n )
{
    /*std::vector <double> M ( dim );
    std::vector <double> CONCENTRATION_MOLE ( dim );
    M[0]=12.011;//"c"
    M[1]=51.996;//"cr"
    M[2]=55.847;//"fe"
    std::vector <double> CONCENTRATION ( dim );
    CONCENTRATION[0]=2.e-3;   //"c"
    CONCENTRATION[1]=2.5e-2;  //"cr"
    CONCENTRATION[2]=1-CONCENTRATION[0]-CONCENTRATION[1];//"fe"
    CONCENTRATION_MOLE=CONCENTRATION_MASStoCONCENTRATION_MOLE ( CONCENTRATION,M );//*/
    //point conc_mole(0.5442 , 0.82531,1);
    
	double step=0.7/n;
	double stepy=1./n;//travaille en x-discretization
	//double step=1./n;
    point pt1,pt2,pt3;
    double minx=0;
    double maxx=1;
    double miny=0;
    double maxy=1;
    std::vector < std::vector<double> > Jacobian(dim-1);
    for ( int i = 0; i < ( dim-1 ); i++ )
        Jacobian[i].resize ( dim-1 );
    for (int l = 0; l < (dim-1); l++)
	for(int k = 0; k < (dim-1); k++)
		Jacobian[l][k]=0.0;
    for (int i=0;i<=n;i++)
        for (int j=0;j<=n;j++)
        {
			
			/*
			double G[2];
            G[0]= ( func1->*nrfuncv1 ) ( i*step,0,j*step,0,0,0 );
            G[1]= ( func2->*nrfuncv2 ) ( i*step,0,j*step,0,0,0 );
            double min=G[0]>G[1]?G[1]:G[0];
            int func=G[0]>G[1]?1:0;
            point pt=point(i*step,j*step,min,func);
            if (i==0&&j==0) pt1=pt;
            else if (i==n&&j==0) pt2=pt;
            else if (i==0&&j==n) pt3=pt;*/
			
			
			/*double G[2];
			G[0]= ( func1->*nrfuncv1 ) ( i*step,0,j*step,0,0,0 );
			G[1]= ( func2->*nrfuncv2 ) ( i*step,0,j*step,0,0,0 );
			double min=G[0]>G[1]?G[1]:G[0];
			int func=G[0]>G[1]?1:0;
			point pt=point(i*step,j*step,min,func);
			pt.to_conc_mole_format();
			pt=point(pt[0],pt[1],min,func);

			if (i==0&&j==0) pt1=pt;
			else if (i==n&&j==0) pt2=pt;
			else if (i==0&&j==n) pt3=pt;
			

            else
            {
				
                double dets=0.0;
                if ( pt.func()==0 )
                    Jacobian = ( Jacob1->*nrJacob1 ) ( i*step,0,j*step,0,0,0 );
                else
                    Jacobian = ( Jacob2->*nrJacob2 ) ( i*step,0,j*step,0,0,0 );
                dets=det ( Jacobian,dim-1 );
                if ( dets<0.0 )
                    continue;
				
                if ( T<=1090. )
                    source_points.add(pt);
                else  if ( pt[0]+pt[1]<1.0 )
					source_points.add ( pt );
				//printf("%lf	%lf	%lf %d\n",pt[0],pt[1],min,func);

			}*/
			
			
			double G[2];
			G[0]= ( func1->*nrfuncv1 ) ( (i*i*step*step)/(3-3*i*i*step*step),0,j*stepy/(1-i*i*step*step),0,0,0 );
			G[1]= ( func2->*nrfuncv2 ) ( (i*i*step*step)/(1-i*i*step*step),0,j*stepy/(1-i*i*step*step),0,0,0 );
			double min=G[0]>G[1]?G[1]:G[0];
			int func=G[0]>G[1]?1:0;
						
			point pt=point(i*i*step*step,j*stepy,min,func);


			
			if (i==0&&j==0) pt1=pt;
			else if (i==n&&j==0) pt2=pt;
			else if (i==0&&j==n) pt3=pt;

            else
            {
				
                double dets=0.0;
                if ( pt.func()==0 )
                    Jacobian = ( Jacob1->*nrJacob1 ) ((i*i*step*step)/(3-3*i*i*step*step),0,j*stepy/(1-i*i*step*step),0,0,0 );
                else
                    Jacobian = ( Jacob2->*nrJacob2 ) ((i*i*step*step)/(1-i*i*step*step),0,j*stepy/(1-i*i*step*step),0,0,0 );
                dets=det ( Jacobian,dim-1 );
                if ( dets<0.0 )
                    continue;
				 else if ( pt[0]+pt[1]<1.0 ){
					source_points.add ( pt );
				//printf("%lf	%lf	%lf %d\n",pt[0],pt[1],min,func);	 
					 }
				

            }
			
			/*double G[2];
			 G[0]= ( func1->*nrfuncv1 ) ( i*i*step*step,0,j*step,0,0,0 );
			 G[1]= ( func2->*nrfuncv2 ) ( i*i*step*step,0,j*step,0,0,0 );
			 double min=G[0]>G[1]?G[1]:G[0];
			 int func=G[0]>G[1]?1:0;
			 
			point pt=point(i*i*step*step,j*step,min,func);
			pt.to_conc_mole_format();
			pt=point(pt[0],pt[1],min,func);
			 
			 
			 
			 if (i==0&&j==0) pt1=pt;
			 else if (i==n&&j==0) pt2=pt;
			 else if (i==0&&j==n) pt3=pt;
			 
			 else
			 {
			 
			 double dets=0.0;
			 if ( pt.func()==0 )
			 Jacobian = ( Jacob1->*nrJacob1 ) (i*i*step*step,0,j*step,0,0,0 );
			 else
			 Jacobian = ( Jacob2->*nrJacob2 ) (i*i*step*step,0,j*step,0,0,0 );
			 dets=det ( Jacobian,dim-1 );
			 if ( dets<0.0 )
			 continue;
			 
			 source_points.add ( pt );
			 //printf("%lf	%lf	%lf %d\n",pt[0],pt[1],min,func);
			 }*/
			 
	
	}		
        /*cout<<"maxx="<<maxx<<";minx="<<minx<<endl;
	cout<<"maxy="<<maxy<<";miny="<<miny<<endl;
	double stepx=(maxx-minx)/n;
	double stepy=(maxy-miny)/n;
	for(int i=0;i<=n;i++)
	  for(int j=0;j<=n;j++)
	  {
	    double G[2];
            G[0]= ( func1->*nrfuncv1 ) (minx+i*stepx,0,miny+j*stepy,0,0,0 );
            G[1]= ( func2->*nrfuncv2 ) (minx+ i*stepx,0,miny+j*stepy,0,0,0 );
            double min=G[0]>G[1]?G[1]:G[0];
            int func=G[0]>G[1]?1:0;
            point pt=point(minx+i*stepx,miny+j*stepy,min,func);
	    source_points.add(pt);
	    cout<<"Add point:"<<pt<<endl;
	  }*/
	hull_points.add(pt1);
	hull_points.add(pt2);
	hull_points.add(pt3);
    firstplane=new hyperplane(pt1,pt2,pt3);
    firstplane->set_outside(false);
    const double* coeffs=firstplane->get_coeffs();
    //for(int i=0;i<=dim;i++)
      //cout<<"coeffs["<<i<<"]="<<coeffs[i]<<endl;
	//for(int i=0;i<=n-1;i++){
	//	stepx=(1./(n*n))*(2*i+1);
	//	printf("%d	%lf	%lf \n",i,i*i*pow(stepx,2),i*stepx);}

}
void qhull::discretization4d (double T,int n )
{
    double step=1./n;
    point pt1,pt2,pt3,pt4;
    std::vector < std::vector<double> > Jacobian(dim-1);
    for ( int i = 0; i < ( dim-1 ); i++ )
        Jacobian[i].resize ( dim-1 );
    for (int l = 0; l < (dim-1); l++)
	for(int k = 0; k < (dim-1); k++)
		Jacobian[l][k]=0.0;
    for (int i=0;i<=n;i++)
        for (int j=0;j<=n;j++)
            for (int k=0;k<=n;k++)
            {
				
				/*
				double G[2];
				G[0]= ( func1->*nrfuncv1 ) ( i*step,0,j*step,k*step,0,0 );
				G[1]= ( func2->*nrfuncv2 ) ( i*step,0,j*step,k*step,0,0 );
				double min=G[0]>G[1]?G[1]:G[0];
				int func=G[0]>G[1]?1:0;
				point pt=point(i*step,j*step,k*step,min,func);
				pt.to_conc_mole_format();
				pt=point(pt[0],pt[1],pt[2],min,func);
				 if (i==0&&j==0&&k==0) pt1=pt;
				 else if (i==n&&j==0&&k==0) pt2=pt;
				 else if (i==0&&j==n&&k==0) pt3=pt;
				 else if (i==0&&j==0&&k==n) pt4=pt;*/
				
				
                double G[2];
                G[0]= ( func1->*nrfuncv1 ) ( i*step,0,j*step,k*step,0,0 );
                G[1]= ( func2->*nrfuncv2 ) ( i*step,0,j*step,k*step,0,0 );
                double min=G[0]>G[1]?G[1]:G[0];
                int func=G[0]>G[1]?1:0;
                point pt=point(i*step,j*step,k*step,min,func);
                if (i==0&&j==0&&k==0) pt1=pt;
                else if (i==n&&j==0&&k==0) pt2=pt;
                else if (i==0&&j==n&&k==0) pt3=pt;
                else if (i==0&&j==0&&k==n) pt4=pt;

                else
                {
                    double dets=0.0;
                    if ( pt.func()==0 )
                        Jacobian = ( Jacob1->*nrJacob1 ) ( i*step,0,j*step,k*step,0,0 );
                    else
                        Jacobian = ( Jacob2->*nrJacob2 ) ( i*step,0,j*step,k*step,0,0 );
                    dets=det ( Jacobian,dim-1 );
                    if ( dets<0.0 )
                        continue;
                    if ( T<=1130.0 )
		    {
                        source_points.add(pt);
		    }
                    else if ( pt[0]+pt[1]+pt[2]<=1.0 )
			  source_points.add ( pt );
                }
            }
	    hull_points.add(pt1);
	    hull_points.add(pt2);
	    hull_points.add(pt3);
	    hull_points.add(pt4);
    firstplane=new hyperplane(pt1,pt2,pt3,pt4);
    firstplane->set_outside(true);
    cout<<"firstplane:"<<*firstplane<<endl;
    const double* coeffs=firstplane->get_coeffs();
    for(int i=0;i<=dim;i++)
      cout<<"coeff["<<i<<"]="<<coeffs[i]<<endl;
}
void qhull::discretization5d (double T,int n )
{
    double step=1./n;
    point pt1,pt2,pt3,pt4,pt5;
    std::vector < std::vector<double> > Jacobian(dim-1);
    for ( int i = 0; i < ( dim-1 ); i++ )
        Jacobian[i].resize ( dim-1 );
    for (int i=0;i<=n;i++)
        for (int j=0;j<=n;j++)
            for (int k=0;k<=n;k++)
                for (int l=0;l<=n;l++)
                {
                    double G[2];
                    G[0]= ( func1->*nrfuncv1 ) ( i*step,0,j*step,k*step,l*step,0 );
                    G[1]= ( func2->*nrfuncv2 ) ( i*step,0,j*step,k*step,l*step,0 );
                    double min=G[0]>G[1]?G[1]:G[0];
                    int func=G[0]>G[1]?1:0;
                    point pt=point(i*step,j*step,k*step,min,func);
                    if (i==0&&j==0&&k==0&&l==0) pt1=pt;
                    else if (i==n&&j==0&&k==0&&l==0) pt2=pt;
                    else if (i==0&&j==n&&k==0&&l==0) pt3=pt;
                    else if (i==0&&j==0&&k==n&&l==0) pt4=pt;
                    else if (i==0&&j==0&&k==0&&l==1) pt5=pt;
                    else
                    {
                        double dets=0.0;
                        if ( pt.func()==0 )
                            Jacobian = ( Jacob1->*nrJacob1 ) ( i*step,0,j*step,k*step,0,0 );
                        else
                            Jacobian = ( Jacob2->*nrJacob2 ) ( i*step,0,j*step,k*step,0,0 );
                        dets=det ( Jacobian,dim-1 );
                        if ( dets<0.0 )
                            continue;
                        if ( T<=1118.0 )
                            source_points.add(pt);
                        else if ( pt[0]+pt[1]+pt[2]<=1.0 )
			      source_points.add ( pt );
                    }
                }
		hull_points.add(pt1);
	    hull_points.add(pt2);
	    hull_points.add(pt3);
	    hull_points.add(pt4);
	    hull_points.add(pt5);
    firstplane=new hyperplane(pt1,pt2,pt3,pt4,pt5);
    firstplane->set_outside(false);
}
void qhull::discretization6d (double T,int n )
{
    double step=1./n;
    point pt1,pt2,pt3,pt4,pt5,pt6;
    std::vector < std::vector<double> > Jacobian(dim-1);
    for ( int i = 0; i < ( dim-1 ); i++ )
        Jacobian[i].resize ( dim-1 );
    for (int i=0;i<=n;i++)
        for (int j=0;j<=n;j++)
            for (int k=0;k<=n;k++)
                for (int l=0;l<=n;l++)
                    for (int m=0;m<=n;m++)
                    {
                        double G[2];
                        G[0]= ( func1->*nrfuncv1 ) ( i*step,0,j*step,k*step,l*step,m*step );
                        G[1]= ( func2->*nrfuncv2 ) ( i*step,0,j*step,k*step,l*step,m*step );
                        double min=G[0]>G[1]?G[1]:G[0];
                        int func=G[0]>G[1]?1:0;
                        point pt=point(i*step,j*step,k*step,min,func);
                        if (i==0&&j==0&&k==0&&l==0&&m==0) pt1=pt;
                        else if (i==n&&j==0&&k==0&&l==0&&m==0) pt2=pt;
                        else if (i==0&&j==n&&k==0&&l==0&&m==0) pt3=pt;
                        else if (i==0&&j==0&&k==n&&l==0&&m==0) pt4=pt;
                        else if (i==0&&j==0&&k==0&&l==n&&m==0) pt5=pt;
                        else if (i==0&&j==0&&k==0&&l==0&&m==n) pt6=pt;
                        else
                        {
                            double dets=0.0;
                            if ( pt.func()==0 )
                                Jacobian = ( Jacob1->*nrJacob1 ) ( i*step,0,j*step,k*step,l*step,m*step );
                            else
                                Jacobian = ( Jacob2->*nrJacob2 ) ( i*step,0,j*step,k*step,l*step,m*step);
                            dets=det ( Jacobian,dim-1 );
                            if ( dets<0.0 )
                                continue;
                            if ( T<=1118.0 )
                                source_points.add(pt);
                            else if ( ( T>=1118.0 ) && ( T<=1730.0 ) ) continue;
                            else if ( T>=1822.0 ) continue;
                            else if ( pt[0]+pt[1]+pt[2]<=1.0 )
				  source_points.add ( pt );
                        }
                    }
    firstplane=new hyperplane(pt1,pt2,pt3,pt4,pt5,pt6);
    firstplane->set_outside(false);
}

vector<point> qhull::extremums(double cr,double mo)
{
    hyperplane* cur=start;
    while (cur)
    {
        int func=(*cur)[0].func();
        int one_func=1;
        for (int i=1;i<dim;i++)
        {
            if ((*cur)[i].func()==func)
                one_func++;
        }
        if (one_func==dim)
            cur=del_plane(cur);
        else cur=cur->next;
    }
    vector<point> exts;
    if (!start)
        return exts;
    switch (dim)
    {
    case 2:
        exts.push_back((*start)[0]);
        exts.push_back((*start)[1]);
        return exts;			
		//exts.push_back(((*start)[0]).from_conc_mole_format());
	    //exts.push_back((*start)[1]).from_conc_mole_format());	
    case 3:
	return extremums3d(cr);
    case 4:
	return extremums4d(cr,mo);
    case 5:
        return extremums5d();
    case 6:
        return extremums6d();
    }
    /*
        if ( el==6 )
        {
            std::vector <double> M ( el );
            M[0]=12.011;//"c"
            M[1]=51.996;//"cr"
            M[2]=95.940;//"mo"
            M[3]=58.933;//"co"
            M[4]=26.982;//"al"
            M[5]=55.847;//fe
            std::vector <double> CONCENTRATION ( el );
            CONCENTRATION[0]=2.e-3;   //"c"
            CONCENTRATION[1]=2.5e-2;  //"cr"
            CONCENTRATION[2]=1.4e-2;  //"mo"
            CONCENTRATION[4]=10e-2;   //"co"
            CONCENTRATION[5]=0.9e-2;  //"al"
            CONCENTRATION[6]=1.-CONCENTRATION[0]-CONCENTRATION[1]-CONCENTRATION[2]-CONCENTRATION[3]-CONCENTRATION[4]-CONCENTRATION[5];//"fe"

            CONCENTRATION_MOLE=CONCENTRATION_MASStoCONCENTRATION_MOLE ( CONCENTRATION,M );//
        }
        if ( el==7 )
        {
            std::vector <double> M ( el );
            M[0]=12.011;//"c"
            M[1]=51.996;//"cr"
            M[2]=95.940;//"mo"
            M[3]=58.933;//"co"
            M[4]=26.982;//"al"
            M[5]=58.690;//"ni"
            M[6]=55.847;//"fe"
            std::vector <double> CONCENTRATION ( el );
            CONCENTRATION[0]=2.e-3;   //"c"
            CONCENTRATION[1]=2.5e-2;  //"cr"
            CONCENTRATION[2]=1.4e-2;  //"mo"
            CONCENTRATION[3]=10e-2;   //"co"
            CONCENTRATION[4]=0.9e-2;  //"al"
            CONCENTRATION[5]=14e-2;   //"ni"
            CONCENTRATION[6]=1.-CONCENTRATION[0]-CONCENTRATION[1]-CONCENTRATION[2]-CONCENTRATION[3]-CONCENTRATION[4]-CONCENTRATION[5];//"fe"

            CONCENTRATION_MOLE=CONCENTRATION_MASStoCONCENTRATION_MOLE ( CONCENTRATION,M );//
        }
    */

    /*
    	    if(el==5)
    	    {
    	      if ( line.size() ==2 )
                    {
                        cout<<"line==2"<<endl;
                        vector<vector<double> > coeffs;
                        cout<<"pts1:"<<pts1<<endl;
                        coeffs.push_back ( hyperplane ( pts1 ).getCoeffs() );
                        int is_null=0;
                        for ( int j=0;j<coeffs[0].size();j++ )
                            if ( coeffs[0][j]==0 )
                                is_null++;
                        if ( is_null==coeffs[0].size() )
                            continue;
                        line.push_back ( pts1[0] );
    		    line.push_back ( pts1[1] );
                        coeffs.push_back ( hyperplane ( line ).getCoeffs() );
                        line.pop_back();
    		    line.pop_back();
                        line.push_back ( pts1[1] );
    		    line.push_back ( pts1[2] );
                        coeffs.push_back ( hyperplane ( line ).getCoeffs() );
                        line.pop_back();
    		    line.pop_back();
    		    line.push_back ( pts1[3] );
    		    line.push_back ( pts1[4] );
                        coeffs.push_back ( hyperplane ( line ).getCoeffs() );
                        line.pop_back();
    		    line.pop_back();
                        double** M=new double*[el-1];
                        for ( int i=0; i<el-1; i++ )
                            M[i] = new double[el-1];
                        for ( int i=0;i<el-1;i++ )
                        {
                            for ( int j=0;j<el-1;j++ )
                            {
                                M[i][j]=coeffs[i][j];
                                cout<<M[i][j]<<" ";
                            }
                            cout<<" | "<<coeffs[i][el-1];
                            cout<<endl;
                        }
                        double* b=new double[el-1];
                        for ( int i=0; i<el-1; i++ )
                            b[i] = -coeffs[i][el-1];
                        double* x=new double[el-1];
                        solve ( M,b,x,el-1 );
                        for ( int i=0; i<el-1; i++ )
                            delete[] M[i];
                        delete[] M;
                        delete[] b;
                        vector<double> newpt1;
                        for ( int i=0; i<el-1; i++ )
                            newpt1.push_back ( x[i] );
                        delete[] x;
                        pts.push_back ( point ( newpt1,pts1[0].function ) );
                    }
                    else
                    {
    		  //pts.pop_back();
                        pts.clear();
    		    vector<point> cur;
                        vector<vector<double> > coeffs;
    		    line.pop_back();
    		    cur=line;

    		    cur.push_back(pts1[0]);
    		    cur.push_back(pts1[1]);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
    		    cur.clear();

    		    cur=pts1;
    		    cur.push_back(pt);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
                        cur.clear();

                        cur.push_back ( line[0] );
    		    cur.push_back ( pt );
    		    cur.push_back(pts1[0]);
    		    cur.push_back(pts1[1]);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
                        cur.clear();

    		    cur.push_back ( line[1] );
    		    cur.push_back ( pt );
    		    cur.push_back(pts1[1]);
    		    cur.push_back(pts1[2]);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
                        cur.clear();

                        double** M=new double*[el-1];
                        for ( int i=0; i<el-1; i++ )
                            M[i] = new double[el-1];
                        for ( int i=0;i<el-1;i++ )
                        {
                            for ( int j=0;j<el-1;j++ )
                                M[i][j]=coeffs[i][j];
                        }
                        double* b=new double[el-1];
                        for ( int i=0; i<el-1; i++ )
                            b[i] = -coeffs[i][el-1];
                        double* x=new double[el-1];
                        solve ( M,b,x,el-1 );
                        vector<double> newpt1;
                        for ( int i=0; i<el-1; i++ )
                            newpt1.push_back ( x[i] );
                        pts.push_back ( point ( newpt1,pts1[0].function ) );


                        coeffs.clear();
    		    cur=line;
                        cur.push_back(pts1[0]);
    		    cur.push_back(pts1[1]);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
    		    cur.clear();

    		    cur=line;
    		    cur.push_back(pt);
    		    cur.push_back(pts1[1]);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
                        cur.clear();

                        cur=line;
    		    cur.push_back ( pt );
    		    cur.push_back(pts1[2]);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
                        cur.clear();

    		    cur=line;
    		    cur.push_back ( pt );
    		    cur.push_back(pts1[0]);
                        coeffs.push_back ( hyperplane ( cur ).getCoeffs() );
                        cur.clear();
                        //pts1.pop_back();
                        for ( int i=0;i<el-1;i++ )
                        {
                            for ( int j=0;j<el-1;j++ )
                                M[i][j]=coeffs[i][j];
                        }
                        for ( int i=0; i<el-1; i++ )
                            b[i] = -coeffs[i][el-1];
                        solve ( M,b,x,el-1 );
                        for ( int i=0; i<el-1; i++ )
                            delete[] M[i];
                        delete[] M;
                        delete[] b;
                        newpt1.clear();
                        for ( int i=0; i<el-1; i++ )
                            newpt1.push_back ( x[i] );
                        delete[] x;
                        pts.push_back ( point ( newpt1,line[0].function ) );
                    }
         */
}

vector<point> qhull::extremums3d(double cr)
{
    std::vector <double> M ( dim );
    std::vector <double> CONCENTRATION_MOLE ( dim );
    vector<point> exts;
    M[0]=12.011;//"c"
    //M[1]=51.996;//"cr"
	M[1]=58.69;//"ni" 
    M[2]=55.847;//"fe"
    std::vector <double> CONCENTRATION ( dim );
    CONCENTRATION[0]=2.e-3;   //"c"
    CONCENTRATION[1]=cr;  //"cr"
    CONCENTRATION[2]=1-CONCENTRATION[0]-CONCENTRATION[1];//"fe"
    CONCENTRATION_MOLE=CONCENTRATION_MASStoCONCENTRATION_MOLE ( CONCENTRATION,M );//
	CONCENTRATION_MOLE[0]=x_i;
	CONCENTRATION_MOLE[1]=x_j;
    point conc_mole(CONCENTRATION_MOLE[0],CONCENTRATION_MOLE[1],0);
    hyperplane* cur=start;
    hyperplane plane;
    while (cur)
    {

        plane=*cur;
		
		//plane.to_conc_mole_format();

        if (!is_point_in_simplex(conc_mole,plane))
        {
            cur=cur->next;
            continue;
        }
        else break;
    }
    if (!cur)
        return exts;
	//cout<<"plane:"<<plane<<endl;
    cout<<"plane:"<<*cur<<endl;
    int first_function=0;
    for (int i=0;i<dim;i++)
    {
        if (plane[i].func()==0)
            first_function++;
    }
    int one=0;
    point pt;
    point pts[2];
    int pts_i=0;
    for (int i=0;i<dim;i++)
    {
        if ((first_function==1&&plane[i].func()==0)||(first_function==2&&plane[i].func()==1))
        {
            one=i;
            pt=plane[i];
        }
        else {pts[pts_i]=plane[i];pts_i++;}
    }

    hyperplane line1=hyperplane(pt,conc_mole);
    hyperplane line2=hyperplane(pts[0],pts[1]);
    const double* coeffs1=line1.get_coeffs();
    const double* coeffs2=line2.get_coeffs();
    vector<double> newpt1;
    //double x= ( ( coeffs1[1]*coeffs2[2]-coeffs1[2]*coeffs2[1] ) / ( coeffs1[0]*coeffs2[1]-coeffs1[1]*coeffs2[0] ) );
    //double y= ( ( coeffs1[2]*coeffs2[0]-coeffs1[0]*coeffs2[2] ) / ( coeffs1[0]*coeffs2[1]-coeffs1[1]*coeffs2[0] ) );
    double** matr=new double*[2];
    for(int i=0;i<2;i++)
      matr[i]=new double[2];
    matr[0][0]=coeffs1[0];
    matr[0][1]=coeffs1[1];
    matr[1][0]=coeffs2[0];
    matr[1][1]=coeffs2[1];
    double dets=det(matr,2);
    matr[0][0]=-coeffs1[2];
    matr[0][1]=coeffs1[1];
    matr[1][0]=-coeffs2[2];
    matr[1][1]=coeffs2[1];
    double dets1=det(matr,2);
    matr[0][0]=coeffs1[0];
    matr[0][1]=-coeffs1[2];
    matr[1][0]=coeffs2[0];
    matr[1][1]=-coeffs2[2];
    double dets2=det(matr,2);
    double x=dets1/dets;
    double y=dets2/dets;
    for(int i=0;i<2;i++)
      delete[] matr[i];
    delete[] matr;
    //pt.from_conc_mole_format();
    exts.push_back(pt.from_conc_mole_format());
    exts.push_back(point(x,y,pts[0].func()).from_conc_mole_format());
    return exts;
}

vector<point> qhull::extremums4d(double cr,double mo)
{
    std::vector <double> M ( dim );
    vector<point> exts;
    std::vector <double> CONCENTRATION_MOLE ( dim );
    M[0]=12.011;//"c"
    M[1]=51.996;//"cr"
    M[2]=95.940;//"mo"
    M[3]=55.847;//"fe"

    std::vector <double> CONCENTRATION ( dim );
    CONCENTRATION[0]=2.e-3;   //"c"
    CONCENTRATION[1]=cr;//2.5e-2;  //"cr"
    CONCENTRATION[2]=mo;//1.4e-2;   //"mo"
    CONCENTRATION[3]=1-CONCENTRATION[0]-CONCENTRATION[1]-CONCENTRATION[2];

    CONCENTRATION_MOLE=CONCENTRATION_MASStoCONCENTRATION_MOLE ( CONCENTRATION,M );//
    point conc_mole(CONCENTRATION_MOLE[0],CONCENTRATION_MOLE[1],CONCENTRATION_MOLE[2],0);
    hyperplane* cur=start;
    hyperplane plane;
    while (cur)
    {
        plane=*cur;
		plane.to_conc_mole_format();

        if (!is_point_in_simplex(conc_mole,plane))
        {
            cur=cur->next;
            continue;
        }
        else break;//{cout<<"conc_mole in plane"<<plane<<endl;cur=cur->next;}
    }
    if (!cur)
        return exts;
    cout<<"plane:"<<plane<<endl;
    cout<<"plane:"<<*cur<<endl;
    int first_function=0;
    for (int i=0;i<dim;i++)
    {
        if (plane[i].func()==0)
            first_function++;
    }
    cout<<plane<<endl;
    point tri1[3];
    point tri2[3];
    point tri3[3];
    point pt;
    int cur_i=0;
    if(first_function==1||first_function==3)
    {
      for (int i=0;i<dim;i++)
        if ((first_function==1&&plane[i].func()==0)||(first_function==3&&plane[i].func()==1))
        {
            pt=plane[i];
        } else { tri1[cur_i]=plane[i]; cur_i++;}
        point ext1=pt;
        /*tri2[0]=tri3[0]=pt;
        tri2[1]=tri3[1]=conc_mole;
        tri2[2]=tri1[0];
        tri3[2]=tri1[1];*/
	cout<<tri1[0]<<endl<<tri1[1]<<endl<<tri1[2]<<endl;
        hyperplane plane1=hyperplane(tri1[0],tri1[1],tri1[2]);
        hyperplane plane2=hyperplane(pt,conc_mole,tri1[0]);
        hyperplane plane3=hyperplane(pt,conc_mole,tri1[1]);
        const double* coeffs1=plane1.get_coeffs();
        const double* coeffs2=plane2.get_coeffs();
        const double* coeffs3=plane3.get_coeffs();
        const double** M=new const double*[dim-1];
        M[0]=coeffs1;
        M[1]=coeffs2;
        M[2]=coeffs3;
        double* b=new double[dim-1];
        for ( int i=0; i<dim-1; i++ )
            b[i] = -M[i][dim-1];
        double* x=new double[dim-1];
        solve ( M,b,x,dim-1 );
        
        point ext2=point(x[0],x[1],x[2],plane1[0].func());
	delete[] b;
        delete[] M;
        delete[] x;
        exts.push_back(ext1.from_conc_mole_format());
        exts.push_back(ext2.from_conc_mole_format());
    }
    else
    {
      int cur_i1=0;
      int cur_i2=0;
        for (int i=0;i<dim;i++)
            if (plane[i].func()==0)
	    {
	      tri1[cur_i1]=plane[i];
	      cur_i1++;
	    }
            else
	    {
	      tri2[cur_i2]=plane[i];
	      cur_i2++;
	    }
	    cout<<tri1[0]<<endl<<tri1[1]<<endl<<tri2[0]<<endl<<tri2[1]<<endl;
        hyperplane plane1=hyperplane(tri2[0],tri2[1],conc_mole);
        hyperplane plane2=hyperplane(tri1[0],tri1[1],tri2[0]);
        hyperplane plane3=hyperplane(tri1[0],tri1[1],tri2[1]);

        const double* coeffs1=plane1.get_coeffs();
        const double* coeffs2=plane2.get_coeffs();
        const double* coeffs3=plane3.get_coeffs();
        const double** M=new const double*[dim-1];
        M[0]=coeffs1;
        M[1]=coeffs2;
        M[2]=coeffs3;
        double* b=new double[dim-1];
	for(int i=0;i<3;i++)
	{
	  for(int j=0;j<=3;j++)
	  {
	    cout<<M[i][j]<<" ";
	  }
	  cout<<endl;
	}
        for ( int i=0; i<dim-1; i++ )
            b[i] = -M[i][dim-1];
        double* x=new double[dim-1];
	
        solve ( M,b,x,dim-1 );
	cout<<"x:";
	for(int j=0;j<=3;j++)
	    cout<<b[j]<<" ";
	cout<<endl;
        point ext1=point(x[0],x[1],x[2],plane2[0].func());

        plane1=hyperplane(tri1[0],tri1[1],conc_mole);
        plane2=hyperplane(tri2[0],tri2[1],tri1[0]);
        plane3=hyperplane(tri2[0],tri2[1],tri1[1]);

        coeffs1=plane1.get_coeffs();
        coeffs2=plane2.get_coeffs();
        coeffs3=plane3.get_coeffs();
        M[0]=coeffs1;
        M[1]=coeffs2;
        M[2]=coeffs3;
        for ( int i=0; i<dim-1; i++ )
            b[i] = -M[i][dim-1];
	for(int i=0;i<3;i++)
	{
	  for(int j=0;j<=3;j++)
	  {
	    cout<<M[i][j]<<" ";
	  }
	  cout<<endl;
	}
        solve ( M,b,x,dim-1 );
	cout<<"x:";
	for(int j=0;j<=3;j++)
	    cout<<b[j]<<" ";
	cout<<endl;
        point ext2=point(x[0],x[1],x[2],plane2[0].func());
        delete[] b;
        delete[] M;
        delete[] x;
        exts.push_back(ext1.from_conc_mole_format());
        exts.push_back(ext2.from_conc_mole_format());
    }
    return exts;
}

vector<point> qhull::extremums5d()
{
    std::vector <double> CONCENTRATION_MOLE ( dim );
    std::vector <double> M ( dim );
    vector<point> exts;
    M[0]=12.011;//"c"
    M[1]=51.996;//"cr"
    M[2]=95.940;//"mo"
    M[3]=58.933;//"co"
    M[4]=55.847;//"fe"

    std::vector <double> CONCENTRATION ( dim );
    CONCENTRATION[0]=2.e-3;   //"c"
    CONCENTRATION[1]=2.5e-2;  //"cr"
    CONCENTRATION[2]=1.4e-2;  //"mo"
    CONCENTRATION[3]=10e-2;  //"co"
    CONCENTRATION[4]=1-CONCENTRATION[0]-CONCENTRATION[1]-CONCENTRATION[2]-CONCENTRATION[3];

    CONCENTRATION_MOLE=CONCENTRATION_MASStoCONCENTRATION_MOLE ( CONCENTRATION,M );//
    point conc_mole(CONCENTRATION_MOLE[0],CONCENTRATION_MOLE[1],CONCENTRATION_MOLE[2],CONCENTRATION_MOLE[3],0);
    hyperplane* cur=start;
    hyperplane plane;
    while (cur)
    {
        point pt1((*cur)[0][0],(*cur)[0][1],(*cur)[0].func());
        point pt2((*cur)[1][0],(*cur)[1][1],(*cur)[1].func());
        point pt3((*cur)[2][0],(*cur)[2][1],(*cur)[2].func());
        pt1.to_conc_mole_format();
        pt2.to_conc_mole_format();
        pt3.to_conc_mole_format();
        plane=hyperplane(pt1,pt2,pt3);

        if (!is_point_in_simplex(conc_mole,plane))
        {
            cur=cur->next;
            continue;
        }
        else break;
    }
    if (!cur)
        return exts;
    int first_function=0;
    for (int i=0;i<dim;i++)
    {
        if (plane[i].func()==0)
            first_function++;
    }

    point tri1[3];
    point tri2[3];
    point tri3[3];
    point pt;
    int one=-1;
    for (int i=0;i<dim;i++)
        if ((first_function==1&&plane[i].func()==0)||(first_function==3&&plane[i].func()==1))
        {
            one=i;
            pt=plane[i];
        } else tri1[i<3?i:i-1]=plane[i];
    if (one!=-1)
    {
        point ext1=pt;
        tri2[0]=tri3[0]=pt;
        tri2[1]=tri3[1]=conc_mole;
        tri2[2]=tri1[0];
        tri3[2]=tri1[1];
        hyperplane plane1=hyperplane(tri1[0],tri1[1],tri1[2]);
        hyperplane plane2=hyperplane(tri2[0],tri2[1],tri2[2]);
        hyperplane plane3=hyperplane(tri3[0],tri3[1],tri3[2]);
        const double* coeffs1=plane1.get_coeffs();
        const double* coeffs2=plane2.get_coeffs();
        const double* coeffs3=plane3.get_coeffs();
        const double** M=new const double*[dim-1];
        M[0]=coeffs1;
        M[1]=coeffs2;
        M[2]=coeffs3;
        double* b=new double[dim-1];
        for ( int i=0; i<dim-1; i++ )
            b[i] = -M[i][dim-1];
        double* x=new double[dim-1];
        solve ( M,b,x,dim-1 );
        delete[] b;
        delete[] M;
        point ext2=point(x[0],x[1],x[2],plane1[0].func());
        delete[] x;
    }
    else
    {
        for (int i=0;i<dim;i++)
            if (plane[i].func()==0)
                tri2[i<2?i:i-2]=plane[i];
            else tri2[i<2?i:i-2]=plane[i];
        hyperplane plane1=hyperplane(tri1[0],tri1[1],conc_mole);
        hyperplane plane2=hyperplane(tri1[0],tri1[1],tri2[0]);
        hyperplane plane3=hyperplane(tri1[0],tri1[1],tri2[1]);

        const double* coeffs1=plane1.get_coeffs();
        const double* coeffs2=plane2.get_coeffs();
        const double* coeffs3=plane3.get_coeffs();
        const double** M=new const double*[dim-1];
        M[0]=coeffs1;
        M[1]=coeffs2;
        M[2]=coeffs3;
        double* b=new double[dim-1];
        for ( int i=0; i<dim-1; i++ )
            b[i] = -M[i][dim-1];
        double* x=new double[dim-1];
        solve ( M,b,x,dim-1 );
        point ext1=point(x[0],x[1],x[2],plane1[0].func());

        plane1=hyperplane(tri2[0],tri2[1],conc_mole);
        plane2=hyperplane(tri2[0],tri2[1],tri1[0]);
        plane3=hyperplane(tri2[0],tri2[1],tri1[1]);

        coeffs1=plane1.get_coeffs();
        coeffs2=plane2.get_coeffs();
        coeffs3=plane3.get_coeffs();
        M[0]=coeffs1;
        M[1]=coeffs2;
        M[2]=coeffs3;
        for ( int i=0; i<dim-1; i++ )
            b[i] = -M[i][dim-1];
        solve ( M,b,x,dim-1 );

        point ext2=point(x[0],x[1],x[2],plane1[0].func());
        delete[] b;
        delete[] M;
        delete[] x;
    }
}

vector<point> qhull::extremums6d()
{
    std::vector <double> M ( dim );
    vector<point> exts;
    std::vector <double> CONCENTRATION_MOLE ( dim );
    return exts;
}

ostream& operator<< ( std::ostream &os, const qhull &qh )
{
    qh.print();
    return os;
}
