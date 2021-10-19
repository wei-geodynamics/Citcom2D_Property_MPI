C0755 6774 element_definitions.h
/*	Header file containing the appropriate definitions for the elements being used 
in the problem.   Include in all of the files linking in to the element program.  */

#define ELN 8  /* for this element, the max nodes/vpoints/ppoints, used in indexing */
#define ELV 8
#define ELP 1
#define EL1N 4
#define EL1V 4
#define EL1P 1
#define GNVI ELV*ELN
#define GNPI ELP*ELN
#define GN1VI EL1V*EL1N
#define GN1PI EL1P*EL1N

#define GNVINDEX(n,v) ((ELV*((n)-1)) + ((v)-1))
#define GNPINDEX(n,p) ((ELP*((n)-1)) + ((p)-1))

#define GNVXINDEX(d,n,v) ((ELV*ELN*(d))+(ELV*((n)-1)) + ((v)-1))
#define GNPXINDEX(d,n,p) ((ELP*ELN*(d))+(ELP*((n)-1)) + ((p)-1))
#define GNVXSHORT(d,i) ((ELV*ELN*(d))+(i))
#define GNPXSHORT(d,i) ((ELP*ELN*(d))+(i))

/*
#define GNVINDEX(n,v) ((ELN*((v)-1)) + (n-1))
#define GNPINDEX(n,p) ((ELN*((p)-1)) + (n-1))
#define GNVXINDEX(d,n,v) ((ELV*ELN*(d))+(ELN*((v)-1)) + (n-1))
#define GNPXINDEX(d,n,p) ((ELP*ELN*(d))+(ELN*((p)-1)) + (n-1))
*/


#define GMVINDEX(n,v) ((EL1V*((n)-1)) + (v-1))
#define GMPINDEX(n,p) ((EL1P*((n)-1)) + (p-1))
#define GMVXINDEX(d,n,v) ((EL1V*EL1N*(d))+(EL1V*((n)-1)) + (v-1))
#define GMPXINDEX(d,n,p) ((EL1P*EL1N*(d))+(EL1P*((n)-1)) + (p-1))

#define GMVGAMMA(i,n) (4*(i) + (n))


/* Element definitions */

static const int enodes[] = {0,2,4,8};
static const int pnodes[] = {0,1,1,1};
static const int vpoints[] = {0,2,4,8};
static const int ppoints[] = {0,1,1,1};
static const int spoints[] = {0,1,2,4};
static const int onedvpoints[] = {0,0,2,4};
static const int onedppoints[] = {0,0,1,1};
static const int additional_dof[] = {0,0,1,2};

/* More cumbersome, but more optimizable versions ! */

#define ENODES3D 8
#define ENODES2D 4
#define VPOINTS3D 8
#define VPOINTS2D 4
#define PPOINTS3D 1
#define PPOINTS2D 1
#define SPOINTS3D 4
#define SPOINTS2D 2
#define V1DPOINTS3D 4
#define V1DPOINTS2D 2
#define P1DPOINTS3D 1
#define P1DPOINTS2D 1

/*  As we use arrays starting at index 1  */
static const int sfnarraysize[] = {0,3,5,9};
static const int onedsfnarsize[] = {0,0,3,5};
static const int gptarraysize[] = {0,3,5,9};
static const int pptarraysize[] = {0,2,2,2};
static const int node_array_size[] = {0,3,5,9};
static const int pnode_array_size[] = {0,2,2,2};
static const int onedgptarsize[] = {0,0,3,5};
static const int onedpptarsize[] = {0,0,2,2};
/*  As we use arrays starting at index 0  */
static const int loc_mat_size[] = {0,4,8,24};
static const int stored_mat_size[] = {0,10,36,300};
static const int node_mat_size[] = {0,3,18,81};
static const int ploc_mat_size[] = {0,1,1,1};
/*  Inter-node information */
static const int max_eqn_interaction[]={0,3,21,81};
static const int max_els_per_node[]={0,2,4,8};

#define MAX_EQN_INT_3D 81
#define MAX_EQN_INT_2D 21

static const int seq[]={0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210,231,253,276,300};



#define BETA  0.57735026918962576451      

/*  4-pt co-ords  */

/*   Declare everything to be static: no changes are made to this data ... it is 
easiest to redeclare it all each time and have just one file of definitions   */

static const struct One_d_int_points {
	double x[2];
	double weight[3]; 	 } 
			g_1d[5] = 
                        {	{{0.0,0.0},{0.0,0.0,0.0}},
				{{-BETA,-BETA},{0.0,1.0,1.0}},
				{{ BETA,-BETA},{0.0,1.0,1.0}}, 
				{{ BETA, BETA},{0.0,1.0,1.0}}, 
				{{-BETA, BETA},{0.0,1.0,1.0}}   } ;

static const struct One_d_int_points 
			l_1d[5] = 
                        {	{{0.0,0.0},{0.0,0.0,0.0}},
				{{-1.0,-1.0},{0.0,1.0,1.0}},
				{{ 1.0,-1.0},{0.0,1.0,1.0}}, 
				{{ 1.0, 1.0},{0.0,1.0,1.0}}, 
				{{-1.0, 1.0},{0.0,1.0,1.0}}   } ;

static const struct One_d_int_points 
			p_1d[2] = 
			{	{{0.0,0.0},{0.0,0.0,0.0}},
				{{0.0,0.0},{0.0,2.0,4.0}} };
	                                                  	   
static const struct Int_points {
	double x[3]; 
	double weight[3];    }
                        g_point[9] = 
			{	{{0.0,0.0,0.0},{0.0,0.0,0.0}},
				{{-BETA,-BETA,-BETA},{0.0,1.0,1.0}},
				{{