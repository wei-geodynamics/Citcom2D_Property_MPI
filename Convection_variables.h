struct CONVECTION { /* information controlling convection problems */   
    char old_T_file[100];
   
    double temp_blob_x[40];
    double temp_blob_y[40];
    double temp_blob_z[40];
    double temp_blob_radius[40]; /* +/- */
    double temp_blob_T[40];
    double temp_blob_bg[40];     /* Reference level if sticky */
    int temp_blob_sticky[40];
    int temp_blobs;
 
    double temp_zonex1[40];
    double temp_zonex2[40];
    double temp_zonez1[40];
    double temp_zonez2[40];
    double temp_zoney1[40];
    double temp_zoney2[40];
    double temp_zonehw[40];
    double temp_zonemag[40];
    int temp_zone_sticky[40];
    int temp_zones;

    double half_space_age;
    int half_space_cooling;

    int number_of_perturbations;
    double perturb_mag[33];
    double perturb_k[33];
    double perturb_ll[33];
    double perturb_mm[33];

    struct SOURCES {
	    int number;
	    double t_offset;
	    double Q[10];
	    double lambda[10];
	}  heat_sources;

    double elasticity1;
  
} convection;


