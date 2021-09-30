/* This file contains the definitions of variables which are passed as arguments */
/* to functions across the whole filespace of CITCOM. #include this file everywhere !*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__osf__)
void *Malloc1();
#endif

#define Malloc0(a) Malloc1((a), __FILE__, __LINE__)

/* #define Malloc0 malloc */

#define LIDN 0x1
#define VBX 0x2
#define VBZ 0x4
#define VBY 0x8
#define TBX 0x10
#define TBZ 0x20
#define TBY 0x40
#define TZEDGE 0x80
#define TXEDGE 0x100
#define TYEDGE 0x200
#define VXEDGE 0x400
#define VZEDGE 0x800
#define VYEDGE 0x1000
#define INTX 0x2000
#define INTZ 0x4000
#define INTY 0x8000
#define SBX 0x10000
#define SBZ 0x20000
#define SBY 0x40000
#define FBX 0x80000
#define FBZ 0x100000
#define FBY 0x200000

#define CBX 0x400000
#define CBZ 0x800000
#define CBY 0x1000000
#define HBX 0x2000000
#define HBZ 0x4000000
#define HBY 0x8000000

#define OFFSIDE 0x10000000
#define SIDEE 0x800000

#define SKIP 0x1000000
#define SKIPID 0x1
#define ZEROID 0x2

#define REFINE1 0x1
#define REFINE2 0x2
#define GREFINE1 0x4
#define GREFINE2 0x8

#define LIDE 1

#ifndef COMPRESS_BINARY
#define COMPRESS_BINARY "/usr/bin/compress"
#endif

#define MAX_LEVELS 12
#define MAX_F 10
#define MAX_S 30

/* Macros */

#define max(A, B) (((A) > (B)) ? (A) : (B))
#define min(A, B) (((A) < (B)) ? (A) : (B))
#define SWAP(a, b)  \
    {               \
        temp = (a); \
        (a) = (b);  \
        (b) = temp; \
    }

typedef double higher_precision;  /* matrix coeffs etc */
typedef double higher_precision1; /* intermediate calculations for finding above coeffs */

/* Common structures */

struct Rect
{
    int numb;
    char overlay[40];
    double x1[40];
    double x2[40];
    double z1[40];
    double z2[40];
    double y1[40];
    double y2[40];
    double halfw[40];
    double mag[40];
};

struct Circ
{
    int numb;
    char overlay[40];
    double x[40];
    double z[40];
    double y[40];
    double rad[40];
    double mag[40];
    double halfw[40];
};

struct Harm
{
    int numb;
    int harms;
    char overlay[40];
    double off[40];
    double x1[40];
    double x2[40];
    double z1[40];
    double z2[40];
    double y1[40];
    double y2[40];
    double kx[20][40];
    double kz[20][40];
    double ky[20][40];
    double ka[20][40];
    double phx[20][40];
    double phz[20][40];
    double phy[20][40];
};

struct Erfc
{
    int numb;
};

struct RectBc
{
    int numb;
    char norm[40];
    double intercept[40];
    double x1[40];
    double x2[40];
    double z1[40];
    double z2[40];
    double halfw[40];
    double mag[40];
};

struct CircBc
{
    int numb;
    char norm[40];
    double intercept[40];
    double x[40];
    double z[40];
    double rad[40];
    double mag[40];
    double halfw[40];
};

struct PolyBc
{
    int numb;
    int order;
    char norm[40];
    double intercept[40];
    double x1[40];
    double x2[40];
    double z1[40];
    double z2[40];
    double ax[20][40];
    double az[20][40];
};

struct HarmBc
{
    int numb;
    int harms;
    char norm[40];
    double off[40];
    double intercept[40];
    double x1[40];
    double x2[40];
    double z1[40];
    double z2[40];
    double kx[20][40];
    double kz[20][40];
    double ka[20][40];
    double phx[20][40];
    double phz[20][40];
};

struct Shape_function_dA
{
    double vpt[8];
    double spt[4];
    double ppt[1];
};

struct Shape_function1_dA
{
    double vpt[6 * 4];
    double ppt[6 * 1];
};

struct Shape_function1
{
    double vpt[4 * 4]; /* node & gauss pt */
    double ppt[4 * 1];
};

struct Shape_function
{
    double vpt[8 * 8]; /* node & gauss pt */
    double spt[8 * 4]; /* node & gauss pt */
    double ppt[8 * 1];
};

struct Shape_function_dx
{
    double vpt[3 * 8 * 8]; /* dirn & node & gauss pt */
    double spt[3 * 8 * 4]; /* dirn & node & gauss pt */
    double ppt[3 * 8 * 1];
};

struct Shape_function1_dx
{
    double vpt[2 * 4 * 4]; /* dirn & node & gauss pt */
    double ppt[2 * 4 * 1];
};

struct EG
{
    higher_precision g[24][1];
};

struct EK2
{
    double k[8 * 8];
};

struct EK
{
    double k[24 * 24];
};

struct MEK
{
    double nint[9];
};

struct NK
{
    higher_precision *k;
    int *map;
};

struct COORD
{
    double centre[4];
    double size[4];
    double recip_size[4];
    double area;
};

struct SUBEL
{
    int sub[9];
};

struct ID
{
    int doff[6];
}; /* can  be 1 or 2 or 3 */
struct IEN
{
    int node[9];
};
struct FNODE
{
    double node[9];
};
struct SIEN
{
    int node[5];
};
struct LM
{
    struct
    {
        int doff[4];
    } node[9];
};

struct NEI
{
    int *nels;
    int *lnode;
    int *element;
};

struct Segment
{

    int zlayers;
    int nzlayer[40];
    double zzlayer[40];
    int xlayers;
    int nxlayer[40];
    double xxlayer[40];
    int ylayers;
    int nylayer[40];
    double yylayer[40];
};

struct IBM_DX
{
    double *x1;
    double *x2;
    int nox;
    int noz;
};

struct BOUND
{
    int bound[4][3];
};

struct PASS
{
    int pass[4][3];
};

struct Parallel
{
    char machinename[160];
    int me;
    int nproc;
    int nprocx;
    int nprocz;
    int nprocy;
    int nprocxy;
    int nproczy;
    int nprocxz;
    int automa;
    int idb;
    int me_loc[4];
    int num_b;

    int *IDD[MAX_LEVELS];
    int *ELE_ORDER[MAX_LEVELS];
    int *NODE_ORDER[MAX_LEVELS];
    struct BOUND *NODE[MAX_LEVELS];
    struct BOUND NUM_NNO[MAX_LEVELS];
    struct BOUND NUM_PASS[MAX_LEVELS];
    struct BOUND NUM_ELE[MAX_LEVELS];
    struct BOUND *IDPASS[MAX_LEVELS];

    int TNUM_PASS[MAX_LEVELS];
    int START_ELE[MAX_LEVELS];
    int START_NODE[MAX_LEVELS];
    struct PASS NUM_NEQ[MAX_LEVELS];
    struct PASS NUM_NODE[MAX_LEVELS];
    struct PASS PROCESSOR[MAX_LEVELS];
    struct PASS *EXCHANGE_ID[MAX_LEVELS];
    struct PASS *EXCHANGE_NODE[MAX_LEVELS];

    int me_sph;
    int no_neighbors;
    int neighbors[27];
    int *neighbors_rev;
    int traces_receive_number[27];
    int traces_transfer_number[27];
    int *traces_transfer_index[27];
};

struct MESH_DATA
{            /* general information concerning the fe mesh */
    int nsd; /* Spatial extent 1,2,3d*/
    int dof; /* degrees of freedom per node */
    int levmax;
    int levmin;
    int levels;
    int mgunitx;
    int mgunitz;
    int mgunity;
    int NEQ[MAX_LEVELS]; /* All other values refer to the biggest mesh (& lid)  */
    int NNO[MAX_LEVELS];
    int NNOV[MAX_LEVELS];
    int NLNO[MAX_LEVELS];
    int NPNO[MAX_LEVELS];
    int NEL[MAX_LEVELS];
    int NOX[MAX_LEVELS];
    int NOZ[MAX_LEVELS];
    int NOY[MAX_LEVELS];
    int NMX[MAX_LEVELS];
    int NXS[MAX_LEVELS];
	int NZS[MAX_LEVELS];
	int NYS[MAX_LEVELS];
    int NNX[MAX_LEVELS][4];
    int ELX[MAX_LEVELS];
    int ELZ[MAX_LEVELS];
    int ELY[MAX_LEVELS];
    int LNDS[MAX_LEVELS];
    int LELS[MAX_LEVELS];
    int neqd;
    int neq;
    int nno;
    int nnov;
    int nlno;
    int npno;
    int nel;
    int snel;
    int nex[4];
    int elz;
    int ely;
    int nnx[4]; /* general form of ... */
    int nox;
    int elx;
    int noz;
    int noy;
    int *exs;
    int ezs;
    int eys;
    int nxs;
    int nzs;
    int nys;
    int nmx;
    int nsf; /* nodes for surface observables */
    int toptbc, topcbc;
    int bottbc, botcbc;
    int topvbc;
    int botvbc;
    int sidevbc;

    char topvbc_file[100];
    char botvbc_file[100];
    char sidevbc_file[100];
    char gridfile[4][100];

    int periodic_x;
    int periodic_y;
    double layer[4]; /* dimensionless dimensions */
    double lidz;
    double bl1width[4], bl2width[4], bl1mag[4], bl2mag[4];
    double hwidth[4], magnitude[4], offset[4], width[4]; /* grid compression information */
    int fnodal_malloc_size;
    int dnodal_malloc_size;
    int feqn_malloc_size;
    int deqn_malloc_size;
    int bandwidth;
    int null_source;
    int null_sink;
    int matrix_size[MAX_LEVELS];
};

struct HAVE
{ /* horizontal averages */
    double *T;
    double *Vi;
    double *Rho;
    double *f;
    double *F;
    double *vrms;
    double *V[4];

    double *T_phase;
    double *P_phase;
    double *density_phase;
    double *Vp_phase;
    double *Vs_phase;
};

struct SLICE
{ /* horizontally sliced data, including topography */
    double *tpg;
    double *tpgb;
    double *grv;
    double *geo;
    double *geok;
    double *grvk;
    double *grvb;
    double *geob;
    double *geobk;
    double *grvbk;
    double *tpgk;
    double *tpgbk;
    double *shflux;
    double *bhflux;
    double *cen_mflux;
    double *vxsurf[3]; /* surface velocity vectors */
    double *vline;     /* for kernels, velocity at force term */
    double *vlinek;
    double *tpglong;
    double *tpgblong;
    double *grvlong;
    double *geolong;
    double *grvblong;
    double *geoblong;

    int minlong;
    int maxlong;
};

struct BAVE
{
    double T;
    double Vi;
    double V[4];
};

struct TOTAL
{
    double melt_prod;
};

struct MONITOR
{
    char node_output[100][6]; /* recording the format of the output data */
    char sobs_output[100][6]; /* recording the format of the output data */
    int node_output_cols;
    int sobs_output_cols;

    int solution_cycles;

    double time_scale;
    double length_scale;
    double viscosity_scale;
    double geoscale;
    double tpgscale;
    double grvscale;

    double delta_v_last_soln;
    double elapsed_time;
    double elapsed_time_vsoln;
    double elapsed_time_vsoln1;
    double reference_stress;
    double incompressibility;
    double vdotv;
    double nond_av_heat_fl;
    double nond_av_adv_hfl;
    double cpu_time_elapsed;
    double cpu_time_on_vp_it;
    double cpu_time_on_forces;
    double cpu_time_on_mg_maps;
    double tpgkmag;
    double grvkmag;

    double Nusselt;
    double Vmax;
    double Vsrms;
    double Vrms;
    double Vrms_surface;
    double Vrms_base;
    double F_surface;
    double F_base;
    double Frat_surface;
    double Frat_base;
    double T_interior;
    double T_maxvaried;
    double Sigma_max;
    double Sigma_interior;
    double Vi_average;
};

struct CONTROL
{
    int PID;

    char output_written_external_command[500]; /* a unix command to run when output files have been created */

    int ORTHO, ORTHOZ;             /* indicates levels of mesh symmetry */
    char B_is_good[MAX_LEVELS];    /* general information controlling program flow */
    char Ahat_is_good[MAX_LEVELS]; /* general information controlling program flow */
    char old_P_file[100];
    char data_file[100];
    char data_file1[100];

    char which_data_files[1000];
    char which_horiz_averages[1000];
    char which_running_data[1000];
    char which_observable_data[1000];

    char PROBLEM_TYPE[20]; /* one of ... */
    int KERNEL;
    int stokes;
    int CONVECTION;
    char comp_adv_method[20];
    int SLAB;
    int PLUME;
    char GEOMETRY[20]; /* one of ... */
    int CART2D;
    int CART2pt5D;
    int CART3D;
    int AXI;
    char SOLVER_TYPE[20]; /* one of ... */
    int DIRECT;
    int CONJ_GRAD;
    int NMULTIGRID;
    int EMULTIGRID;
    int DIRECTII;
    char NODE_SPACING[20]; /* turns into ... */
    int GRID_TYPE;
    int COMPRESS;
    int DX;
    int CONMAN;
    int visc_heating;
    int adi_heating;
    int latent_heating;

    int composition;
    double comp_diff;
    double z_comp;
    int composition_phasechange;

    int initialTOption;
    int ocean_lith;
    int platemodel;
    double inter_temp;
    double age_left;
    double age_right;
    double age_loc_left;
    double age_loc_right;
    double depth_lith;
    double depth_ocean_lith;
    int ocean_lith_margin;
    double ocean_lith_margin_curve;
    double dip_center_x;
    double dip_center_z;
    double dip_margin, dip_margin_left;
    double depth_lith_margin;
    int phasefile_Complete;
    int llmax;

    int imposevelo;
    int age_total;
    char velo_file_pre[255];
    char velo_file_post[255];
    double timescale;
    double age_total_double;

    int continent_lith;
    double depth_continent_lith;
    double depth_continent_crust;
    double continent_loc_left;
    double continent_loc_right;
    int continent_platemodel;
    double continent_age;

    double velo_surf_loc_mid;
    double velo_surf_mag_right;
    double velo_surf_width_right;
    int trechmigrate;
    double velo_surf_mag_left;
    double velo_surf_width_left;
    double velo_surf_loc_mid_rate;
    double velo_surf_width_mid;
    double velo_surf_corner_right;
    double velo_surf_loc_left_overshoot;

    double temp_slabvisc;
    double viscincreaseslab;
    int slab_visc;
    double slab_visc_depth;
    int visc_const_cor;
    int visc_leftcor;
    int visc_rightcor;
    int visc_mid;
    double visc_weakzone;
    double x_weakzone_leftcor;
    double x_weakzone_rightcor;
    double x_weakzone_mid_left;
    double x_weakzone_mid_right;
    double z_weakzone_leftcor;
    double z_weakzone_rightcor;
    double z_weakzone_mid;
    double dip_weakzone_mid;
    int visc_mid_dip;

    int CBF;
    int phasevisc;
    double phaseTop;
    double phaseBot;
    double phaseffactorLa;
    double phaseffactorSm;
    int phasevisc_slab;
    double ViscReduce;
    double phasevisc_slab_T;
    /* Wei add 2020 0806 */
    int phasevisc_C, phasevisc_d, phasevisc_dmax;
    double phasevisc_C_time;
    double phasevisc_d0, phasevisc_dp, phasevisc_dm, phasevisc_E, phasevisc_d_phasereduce, phasevisc_tratio, phasevisc_dmax_value;
    int lesstimeinte;
    double lesstimeinte_number;
    int Visc_C;
    double ViscReduce_C;
    /* */
    int dfact;
    double penalty;
    int augmented_Lagr;
    double augmented;
    int macroele;
    int faults;
    int NASSEMBLE;
    int comparison;
    int crust;
    int restart;
    int restart_frame;
    double plate_vel;

    double tole_comp;

    double sob_tolerance;
    int tracer_correction;
    int meshReadIn;
    char fileMeshX[200], fileMeshZ[200];
    double Ra_temp, Ra_comp, Ra_comp_a;
    double Ra_670, clapeyron670, transT670, width670;
    double Ra_410, clapeyron410, transT410, width410;
    double Ra_670_basalt, clapeyron670_basalt, transT670_basalt, width670_basalt;
    int temperature_perturbation;
    int phasefile, phasefile_C, phasefile_buoyancy, phasefile_C_num_nno, phasefile_C_num_marker, phasefile_C_num_element, phasefile_buoyancy_correction, phasefile_C_flavor;
    double phasefile_buoyancy_depth, phasefile_buoyancy_continent;
    double phasefile_buoyancy_crust, phasefile_buoyancy_crust_depth;
    char phasefile_basa[200], phasefile_pyro[200], phasefile_harz[200], phasefile_PREM[200];
    int phasefile_noT, phasefile_noP, phasefile_noPREM;
    int phasefile_C_mat_mineral[200];
    int crust_generation, harz_generation;
    double depth_harz;
    double *phase_T, *phase_P;
    double *phase_basa_density, *phase_basa_Vp, *phase_basa_Vs;
    double *phase_pyro_density, *phase_pyro_Vp, *phase_pyro_Vs;
    double *phase_harz_density, *phase_harz_Vp, *phase_harz_Vs;
    double *phase_PREM_depth, *phase_PREM_density, *phase_PREM_Vp, *phase_PREM_Vs, *phase_PREM_P;
    double *phase_density, *phase_Vp, *phase_Vs;
    double adi_um, adi_lm;
    double Ts;
    double VBXtopval;
    double VBXbotval;
    double VBYtopval;
    double VBYbotval;

    double TBCtopval, CBCtopval;
    double TBCbotval, CBCbotval;

    double Q0;

    int precondition;
    int vprecondition;
    int keep_going;
    int v_steps_low;
    int v_steps_high;
    int v_steps_upper;
    int max_vel_iterations, min_vel_iterations;
    int p_iterations;
    int max_same_visc;
    double max_res_red_each_p_mg;
    double sub_stepping_factor;
    int mg_cycle;
    int true_vcycle;
    int down_heavy;
    int up_heavy;
    int depth_dominated;
    int eqn_viscosity;
    int eqn_zigzag;
    int verbose;
    double accuracy;
    double vaccuracy;

    int total_iteration_cycles;
    int total_v_solver_calls;

    int record_every;
    int record_all_until;

    int print_convergence;
    int sdepv_print_convergence;

    /* modules */
    int MELTING_MODULE;
    int CHEMISTRY_MODULE;
};

struct Radioactive_Heat
{
    int num;
    double total;
    double concen_u;
    double percent[10];
    double heat_g[10];
    double decay_t[10];
    double concen[10];
};

struct DATA
{
    double layer_km;
    double layer_meter;
    double grav_acc;
    double therm_exp;
    double therm_exp_factor;
    double visc_factor;
    double Cp;
    double disptn_number;
    double therm_diff;
    double therm_cond;
    double density;
    double res_density;
    double res_density_X;
    double melt_density;
    double density_above;
    double density_below;

    double gas_const;
    double surf_heat_flux;
    double ref_viscosity;
    double melt_viscosity;
    double permeability;
    double grav_const;
    double surf_temp;
    double youngs_mod;
    double Te;
    double ref_temperature;
    double Tsurf;
    double T_sol0;
    double delta_S;
    double dTsol_dz;
    double dTsol_dF;
    double dT_dz;
};

struct All_variables
{
#include "Convection_variables.h"
#include "viscosity_descriptions.h"
#include "temperature_descriptions.h"
#include "advection.h"

    FILE *fp;
    FILE *filed[20];
    struct HAVE Have;
    struct BAVE Bulkave;
    struct TOTAL Total;
    struct MESH_DATA mesh;
    struct MESH_DATA lmesh;
    struct CONTROL control;
    struct MONITOR monitor;
    struct DATA data;
    struct SLICE slice;
    struct COORD *eco;
    struct IBM_DX ibm_dx;
    struct IEN *ien; /* global */
    struct SIEN *fault_ien;
    struct SIEN *sien;
    struct ID *id;
    struct COORD *ECO[MAX_LEVELS];
    struct IEN *IEN[MAX_LEVELS];   /* global at each level */
    struct FNODE *TWW[MAX_LEVELS]; /* */
    struct ID *ID[MAX_LEVELS];
    struct NEI NEI[MAX_LEVELS];
    struct SUBEL *EL[MAX_LEVELS];
    struct EG *elt_del[MAX_LEVELS];
    struct EK *elt_k[MAX_LEVELS];
    struct Segment segment;
    struct Radioactive_Heat rad_heat;
    struct Parallel parallel;
    higher_precision *Eqn_k[MAX_LEVELS];
    int *Node_map[MAX_LEVELS];
    int *Node_eqn[MAX_LEVELS];
    int *Node_k_id[MAX_LEVELS];

    float *RVV[27], *PVV[27];
    double *RXX[27], *PXX[27];
    int *RINS[27], *PINS[27];

    double *BI[MAX_LEVELS]; /* inv of  diagonal elements of K matrix */
    double *BPI[MAX_LEVELS];
    double *V[4];    /* velocity X[dirn][node] can save memory */
    double *V1[4];   /* velocity X[dirn][node] can save memory */
    double *Vest[4]; /* velocity X[dirn][node] can save memory */

    double *VO[4], *XMC[4], *XMCpred[4];
    double *Vpred[4];
    int *C12, *CElement;
    int *traces_leave;
    int *traces_leave_index;

    double *P, *F, *H, *S, *U;
    double *Psi;
    double *NP;
    double *edot;                    /* strain rate invariant */
    double *MASS[MAX_LEVELS], *Mass; /* lumped mass matrix (diagonal elements) for p-g solver etc. */
    double *tw;
    double *stress;
    double *XP[4], XG1[4], XG2[4], *X[4],*SX[4], *XX[MAX_LEVELS][4], *Interp[MAX_LEVELS][4];
    double *ZZ;
    double *T, *C, *CE, *IHeat, *buoyancy, *T_old;

    double *T_phase, *P_phase, *density_phase, *Vp_phase, *Vs_phase;

    double *Tdot, *Cdot;
    double *Cphasedot, *Cphase, *Cphase_marker, *CphaseE, *Cphase_marker_old, *Cphase_old, *CphaseE_old, *Cphasedotnum, *Cphase_node;
    double *Tphase_node, *Pphase_node, *Tphase_marker, *Pphase_marker;
    double *d, *dE, *d_marker, *d_marker_old, *d_node, *d_dot, *d_dotnum;
    int C_phasefile_num_nno, C_phasefile_num_element, C_phasefile_num_marker;
    int *C_phasefile_marker_int[100];
    double *C_phasefile_nno[100], *C_phasefile_element[100], *C_phasefile_marker_double[100];

    double *heating_visc, *heating_adi, *heating_latent;
    double *Fas670, *Fas410, *Fas670_basalt;
    double *Fas670_b, *Fas410_b, *Fas670_basalt_b;
    double *Fas670_all_b;
    double *Vi, *EVi;
    double *VI[MAX_LEVELS];        /* viscosity has to soak down to all levels */
    double *EVI[MAX_LEVELS];       /* element viscosity has to soak down to all levels */
    double *VB[4], *TB[4], *CB[4]; /* boundary conditions for V,T defined everywhere */
    double *TW[MAX_LEVELS];        /* nodal weightings */

    double *harm_geoid_from_bncy[2], *harm_geoid_from_bncy_botm[2], *harm_tpgt[2], *harm_tpgb[2];
    double *harm_geoid_from_tpgt[2], *harm_geoid_from_tpgb[2];
    double *harm_geoid[2];

    int num_zero_resid[MAX_LEVELS];
    int *zero_resid[MAX_LEVELS];
    int *surf_element;
    int *surf_node;
    int *mat;           /* properties of mat */
    unsigned int *node; /* properties of node */
    unsigned int *NODE[MAX_LEVELS];
    unsigned int *Element;
    unsigned int *eqn;
    unsigned int *EQN[MAX_LEVELS];

    double **global_K; /* direct solver stuff */
    double **factor_K;
    double *global_F;
    struct LM *lmd;
    struct LM *lm;
    struct LM *LMD[MAX_LEVELS];

    struct Shape_function1 M; /* master-element shape funtions */
    struct Shape_function1_dx Mx;
    struct Shape_function N;
    struct Shape_function_dx Nx;
    struct Shape_function1 L; /* master-element shape funtions */
    struct Shape_function1_dx Lx;

    void (*build_forcing_term)();
    void (*iterative_solver)();
    void (*next_buoyancy_field)();
    void (*obtain_gravity)();
    void (*problem_settings)();
    void (*problem_derived_values)();
    void (*problem_allocate_vars)();
    void (*problem_boundary_conds)();
    void (*problem_node_positions)();
    void (*problem_update_node_positions)();
    void (*problem_initial_fields)();
    void (*problem_update_bcs)();
    void (*special_process_new_velocity)();
    void (*special_process_new_buoyancy)();
    void (*solve_stokes_problem)();
    void (*solver_allocate_vars)();
    void (*transform)();

    double (*node_space_function[3])();
};
