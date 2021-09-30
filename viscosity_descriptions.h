/* in this file define the contents of the VISC_OPT data structure
   which is used to store information used to create predefined 
   viscosity fields, those determined from prior input, those
   related to temperature/pressure/stress/anything else. */


struct VISC_OPT {
    void (* update_viscosity)();
  
    int update_allowed;		/* determines whether visc field can evolve */
    int EQUIVDD;			/* Whatever the structure, average in the end */
    int equivddopt;
    int proflocx;			/* use depth dependence from given x,y location */
    int proflocy;
    int SMOOTH;
    int smooth_cycles;
  

    char STRUCTURE[20];		/* which option to determine viscosity field, one of .... */
    int FROM_SYSTEM;
    int FROM_FILE;
    int FROM_SPECS;
  
				/* System ... */
    int RHEOL;			/* 1,2 */
    int rheol_layers;
    int num_mat;

    int nlm;
    int n410;
    int nlith;
    int ncrust1;
    int ncrust2;
    int ndd;
    double zlm;
    double z410;
    double zlith;
    double zcrust1;
    double zcrust2;
    double zdd;
    double z1000;
    double z300;    
    double zbasalt;

    int FREEZE;
    double freeze_thresh;
    double freeze_value;

    int MAX;
    double max_value;
    int MIN;
    double min_value;

    int SDEPV;
    int lower_diff;
    double sdepv_misfit;
    int sdepv_normalize;
    double sdepv_expt[40];
    double sdepv_trns[40];

    int TDEPV;
    int TDEPV_AVE;
    double N0[40];
    double E[40],T0[40];
    double T[40],Z[40];

    int PDEPV;
     
    int weak_blobs;
    double weak_blobx[40];
    double weak_bloby[40];
    double weak_blobz[40];
    double weak_blobwidth[40];
    double weak_blobmag[40];
   
    int weak_zones;
    double weak_zonex1[40];
    double weak_zoney1[40];
    double weak_zonez1[40];
    double weak_zonex2[40];
    double weak_zoney2[40];
    double weak_zonez2[40];
  
    double weak_zonewidth[40];
    double weak_zonemag[40];
  
    int guess;
    char old_file[100];
				/* Specification info */
  
				/* Prespecified viscosity parameters */
    char VISC_OPT[20];

    int layers;			/* number of layers with properties .... */
    double layer_depth[40];
    double layer_visc[40];

    int SLABLVZ;			/* slab structure imposed on top of 3 layer structure */
    int slvzd1,slvzd2,slvzd3;	        /* layer thicknesses (nodes) */
    int slvzD1,slvzD2;		        /* slab posn & length */
    double slvzn1,slvzn2,slvzn3,slvzN;   /* viscosities */

    int COSX;
    double cosx_epsilon;
    double cosx_k;
    int cosx_exp;
 
    int EXPX;
    double expx_epsilon;
 
    /* MODULE BASED VISCOSITY VARIATIONS */

    int RESDEPV;
    double RESeta0[40];

    int CHEMDEPV;
    double CH0[40];
    double CHEMeta0[40];

    int visc_platebond,visc_platebond_selfadapt;
    double z_weakzone_platebond;
    double width_weakzone_platebond;

    double left_weakzone_platebond;
    double left_weakzone_platebond_ratio;
    double right_weakzone_platebond;
    double visc_reduce_platebond;
    int visc_platebond_const;
    double weakzone_ratio;
  
} viscosity;
