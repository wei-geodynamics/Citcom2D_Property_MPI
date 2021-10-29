/* Set up the finite element problem to suit: returns with all memory */
/* allocated, temperature, viscosity, node locations and how to use */
/* them all established. 8.29.92 or 29.8.92 depending on your nationality*/
#include <signal.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"

int Emergency_stop;

void read_instructions(E, argc, argv) struct All_variables *E;
int argc;
char **argv;
{
  int get_process_identifier();

  void allocate_common_vars();
  void common_initial_fields();
  void read_initial_settings();
  void global_default_values();
  void global_derived_values();
  void construct_ien();
  void construct_masks();
  void construct_shape_functions();
  void construct_id();
  void construct_lm();
  void construct_sub_element();
  void mass_matrix();
  void construct_node_ks();
  void construct_node_maps();
  void construct_mat_group();
  void interuption();
  void set_up_nonmg_aliases();
  void check_bc_consistency();
  void node_locations();
  void allocate_velocity_vars();
  void parallel_domain_decomp1();
  void parallel_shuffle_ele_and_id();
  void parallel_communication_routs();
  void parallel_process_termination();
  void setup_parser();

  double start_time, CPU_time0(), vmag;
  double vdot();

  int *temp, i;
  void initial_phasefile();

  /* =====================================================
       Global interuption handling routine defined once here
       =====================================================  */
  if (E->parallel.me == 0)
    start_time = CPU_time0();
  Emergency_stop = 0;
  signal(SIGINT, interuption);
  signal(SIGTERM, interuption);

  E->control.PID = get_process_identifier();

  /* ==================================================
       Initialize from the command line 
       from startup files. (See Parsing.c).
       ==================================================  */
  if (E->parallel.me == 0)
    fprintf(stderr, "ok1\n");
  setup_parser(E, argc, argv);

  global_default_values(E);

  if (E->parallel.me == 0)
    fprintf(stderr, "ok2\n");
  read_initial_settings(E);

  if (E->parallel.me == 0)
    fprintf(stderr, "ok3\n");
  (E->problem_derived_values)(E); /* call this before global_derived_  */
  global_derived_values(E);

  if (E->control.phasefile)
  {
    fprintf(stderr, "start phasefile\n");
    initial_phasefile(E);
  }

  if (E->parallel.me == 0)
    fprintf(stderr, "ok4\n");
  parallel_domain_decomp1(E);

  if (E->parallel.me == 0)
    fprintf(stderr, "ok5\n");
  allocate_common_vars(E);
  if (E->parallel.me == 0)
    fprintf(stderr, "ok6\n");
  (E->problem_allocate_vars)(E);
  (E->solver_allocate_vars)(E);

  if (E->parallel.me == 0)
    fprintf(stderr, "ok6a\n");
  /* logical domain */
  construct_ien(E);
  if (E->parallel.me == 0)
    fprintf(stderr, "ok9\n");
  construct_sub_element(E);

  /* physical domain */
  node_locations(E);
  //allocate_velocity_vars(E);
  (E->problem_boundary_conds)(E);

  construct_masks(E); /* order is important here */

  if (E->parallel.me == 0)
    fprintf(stderr, "ok10\n");
  construct_id(E);
  construct_lm(E);
  if (E->parallel.me == 0)
    fprintf(stderr, "ok11\n");
  construct_mat_group(E);

  check_bc_consistency(E);

  if (E->parallel.me == 0)
    fprintf(stderr, "ok13\n");

  parallel_shuffle_ele_and_id(E);

  if (E->parallel.me == 0)
    fprintf(stderr, "ok14\n");
  parallel_communication_routs(E);
  fprintf(stderr, "ok15\n");
  construct_shape_functions(E);
  mass_matrix(E);
  fprintf(stderr, "ok16\n");
  (E->problem_initial_fields)(E); /* temperature/chemistry/melting etc */
  common_initial_fields(E);       /* velocity/pressure/viscosity (viscosity must be done LAST) */
  fprintf(stderr, "ok17\n");
  shutdown_parser(E);
  /*
*/
  fprintf(stderr, "ok18\n");
  return;
}

/* ===================================
   Functions which set up details 
   common to all problems follow ...
   ===================================  */

void allocate_common_vars(E) struct All_variables *E;

{
  double **dmatrix();
  double **fmatrix();
  void set_up_nonmg_aliases();
  void allocate_velocity_vars();
  int i, j, l, nno_l, npno_l, nozl, nnov_l, nxyz, nox, noy, noz;

  E->mesh.fnodal_malloc_size = (E->lmesh.nno + 2) * sizeof(double);
  E->mesh.dnodal_malloc_size = (E->lmesh.nno + 2) * sizeof(double);
  E->mesh.feqn_malloc_size = (E->mesh.nsd * E->lmesh.nno + 2) * sizeof(double);
  E->mesh.deqn_malloc_size = (E->mesh.nsd * E->lmesh.nno + 2) * sizeof(double);

  E->S = (double *)malloc((E->lmesh.npno + 1) * sizeof(double));

  E->Cphase = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->CphaseE = (double *)malloc((E->lmesh.nel + 1) * sizeof(double));
  E->Cphase_old = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->CphaseE_old = (double *)malloc((E->lmesh.nel + 1) * sizeof(double));
  //    E->d        = (double *) malloc((E->mesh.nno+1)*sizeof(double));
  //    E->dE       = (double *) malloc((E->mesh.nel+1)*sizeof(double));
  E->T_old = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->P = (double *)malloc((E->lmesh.npno + 1) * sizeof(double));

  E->F = (double *)malloc((E->mesh.nsd * E->lmesh.nnov + 1) * sizeof(double));
  E->U = (double *)malloc((E->lmesh.neq + 2) * sizeof(double));
  E->T = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->C = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->CE = (double *)malloc((E->lmesh.nel + 1) * sizeof(double));
  E->buoyancy = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->NP = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->heatflux = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->heatflux_adv = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->edot = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));

  E->Fas670 = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->Fas670_b = (double *)malloc((E->lmesh.nsf + 1) * sizeof(double));
  E->Fas410 = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->Fas410_b = (double *)malloc((E->lmesh.nsf + 1) * sizeof(double));
  E->Fas670_basalt = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->Fas670_basalt_b = (double *)malloc((E->lmesh.nsf + 1) * sizeof(double));
  E->Fas670_all_b = (double *)malloc((E->lmesh.nsf + 1) * sizeof(double));
  E->heating_adi = (double *)malloc((E->lmesh.nel + 1) * sizeof(double));
  E->heating_visc = (double *)malloc((E->lmesh.nel + 1) * sizeof(double));
  E->heating_latent = (double *)malloc((E->lmesh.nel + 1) * sizeof(double));

  E->diffusivity = (double *)malloc((E->lmesh.noz + 1) * sizeof(double));
  E->expansivity = (double *)malloc((E->lmesh.noz + 1) * sizeof(double));
  for (i = 1; i <= E->mesh.nsd; i++)
  {
    E->TB[i] = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
    E->V[i] = (double *)malloc((E->lmesh.nnov + 1) * sizeof(double));
    E->VB[i] = (double *)malloc((E->lmesh.nnov + 1) * sizeof(double));
    E->CB[i] = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  }

  E->Have.f = (double *)malloc((E->mesh.noz + 1) * sizeof(double));
  E->Have.F = (double *)malloc((E->mesh.noz + 1) * sizeof(double));
  E->Have.T = (double *)malloc((E->mesh.noz + 1) * sizeof(double));
  E->Have.C = (double *)malloc((E->mesh.noz + 1) * sizeof(double));
  E->Have.Vi = (double *)malloc((E->mesh.noz + 1) * sizeof(double));
  E->Have.Rho = (double *)malloc((E->mesh.noz + 1) * sizeof(double));
  E->Have.vrms = (double *)malloc((E->mesh.noz + 1) * sizeof(double));
  E->Have.Tadi = (double *)malloc((E->mesh.noz + 1) * sizeof(double));

  E->segment.Tz = (double *)malloc((E->lmesh.noz + 1) * sizeof(double));
  E->segment.Vx = (double *)malloc((E->lmesh.noz + 1) * sizeof(double));

  E->stress = (double *)malloc((E->lmesh.nsf * 12 + 12) * sizeof(double));
  E->slice.tpg = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.tpgb = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.vline = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.vlinek = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.shflux = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.bhflux = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.cen_mflux = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.vxsurf[1] = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));
  E->slice.vxsurf[2] = (double *)malloc((E->lmesh.nsf + 2) * sizeof(double));

  E->mat = (int *)malloc((E->lmesh.nel + 2) * sizeof(int));

  E->XP[1] = (double *)malloc((E->lmesh.nox + 1) * sizeof(double));
  E->XP[2] = (double *)malloc((E->lmesh.noz + 1) * sizeof(double));
  E->lmesh.rnoz = 50 * E->lmesh.elz + 1;

  /* set up memory for different grids  */
  for (i = E->mesh.levmin; i <= E->mesh.levmax; i++)
  {
    for (j = 1; j <= E->mesh.nsd; j++)
      E->XX[i][j] = (double *)malloc((E->lmesh.NNO[i] + 1) * sizeof(double));

    E->MASS[i] = (double *)malloc((E->lmesh.NNO[i] + 1) * sizeof(double));

    E->ECO[i] = (struct COORD *)malloc((E->lmesh.NNO[i] + 2) * sizeof(struct COORD));
    E->IEN[i] = (struct IEN *)malloc((E->lmesh.NNO[i] + 2) * sizeof(struct IEN));
    E->ID[i] = (struct ID *)malloc((E->lmesh.NNO[i] + 2) * sizeof(struct ID));
    E->GNX[i] = (struct Shape_function_dx *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct Shape_function_dx));
    E->GDA[i] = (struct Shape_function_dA *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct Shape_function_dA));
    E->EL[i] = (struct SUBEL *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct SUBEL));
    E->LMD[i] = (struct LM *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct LM));

    E->EVI[i] = (double *)malloc((E->lmesh.NEL[i] + 2) * vpoints[E->mesh.nsd] * sizeof(double));

    E->TW[i] = (double *)malloc((E->lmesh.NNO[i] + 2) * sizeof(double));
    E->VI[i] = (double *)malloc((E->lmesh.NNO[i] + 2) * sizeof(double));
    E->NODE[i] = (unsigned int *)malloc((E->lmesh.NNO[i] + 2) * sizeof(unsigned int));
    E->TWW[i] = (struct FNODE *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct FNODE));

    E->NEI[i].nels = (int *)malloc((E->lmesh.NNO[i] + 2) * sizeof(int));
    E->NEI[i].lnode = (int *)malloc((E->lmesh.NNO[i] + 2) * enodes[E->mesh.nsd] * sizeof(int));
    E->NEI[i].element = (int *)malloc((E->lmesh.NNO[i] + 2) * enodes[E->mesh.nsd] * sizeof(int));

    E->elt_del[i] = (struct EG *)malloc((E->lmesh.NEL[i] + 1) * sizeof(struct EG));

    E->BI[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
    E->BPI[i] = (double *)malloc((E->lmesh.NPNO[i] + 1) * sizeof(double));
    E->control.B_is_good[i] = 0;
  }

  E->temp = (double *)malloc((E->lmesh.NEQ[E->mesh.levmax] + 2) * sizeof(double));
  E->Element = (unsigned int *)malloc((E->lmesh.nel + 2) * sizeof(unsigned int));

  for (i = E->mesh.levmin; i <= E->mesh.levmax; i++)
  {
    nox = E->lmesh.NOX[i];
    noy = E->lmesh.NOY[i];
    noz = E->lmesh.NOZ[i];
    if (E->mesh.nsd == 2)
    {
      nxyz = max(nox, noz);
    }
    else if (E->mesh.nsd == 3)
    {
      nxyz = max(nox * noz, nox * noy);
      nxyz = max(nxyz, noz * noy);
    }

    E->parallel.IDD[i] = (int *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(int));
    E->parallel.NODE[i] = (struct BOUND *)malloc((nxyz + 2) * sizeof(struct BOUND));
    E->parallel.IDPASS[i] = (struct BOUND *)malloc((10) * sizeof(struct BOUND));
    E->parallel.EXCHANGE_NODE[i] = (struct PASS *)malloc((nxyz + 2) * sizeof(struct PASS));
    E->parallel.EXCHANGE_ID[i] = (struct PASS *)malloc((nxyz * E->mesh.nsd + 3) * sizeof(struct PASS));

    if (i == E->mesh.levmax)
    {
      E->sien = (struct SIEN *)malloc((nxyz + 2) * sizeof(struct SIEN));
      E->surf_element = (int *)malloc((nxyz + 2) * sizeof(int));
      E->surf_node = (int *)malloc((E->lmesh.nsf + 2) * sizeof(int));
    }
  }

  for (l = E->mesh.levmin; l <= E->mesh.levmax; l++)
  {
    for (i = 1; i <= E->lmesh.NNO[l]; i++)
    {
      E->NODE[l][i] = (INTX | INTY | INTZ); /* and any others ... */
      E->VI[l][i] = 1.0;
      E->TW[l][i] = 0.0;
    }

    for (i = 0; i < E->lmesh.NEQ[l]; i++)
    {
      E->BI[l][i] = 0.0;
      E->parallel.IDD[l][i] = 0;
    }
  }

  for (i = 1; i <= E->lmesh.nnov; i++)
    for (j = 1; j <= E->mesh.nsd; j++)
      E->V[j][i] = E->VB[j][i] = 0.0;

  for (i = 1; i <= E->lmesh.nno; i++)
    for (j = 1; j <= E->mesh.nsd; j++)
      E->TB[j][i] = 0.0;

  for (i = 1; i <= E->lmesh.npno; i++)
    E->P[i] = 0.0;

  for (i = 1; i <= E->lmesh.nel; i++)
  {
    E->mat[i] = 1;
    E->heating_visc[i] = 0;
    E->heating_adi[i] = 0;
    E->heating_latent[i] = 1;
  }

  for (i = 1; i <= E->lmesh.noz; i++)
  {
    E->diffusivity[i] = E->expansivity[i] = 1.0;
    E->Have.Tadi[i] = .0;
  }

  for (i = 1; i <= E->lmesh.nno; i++)
    E->T[i] = E->buoyancy[i] = 0.0;
  for (i = 0; i <= E->lmesh.neq + 1; i++)
    E->U[i] = 0.0;
  E->harm_geoid_from_bncy[0] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid_from_bncy[1] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid_from_bncy_botm[0] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid_from_bncy_botm[1] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_tpgt[0] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_tpgt[1] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_tpgb[0] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_tpgb[1] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid_from_tpgt[0] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid_from_tpgt[1] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid_from_tpgb[0] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid_from_tpgb[1] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid[0] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  E->harm_geoid[1] = (double *)malloc((E->control.llmax + 1) * sizeof(double));
  if (E->control.phasefile)
  {
    E->T_phase = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
    E->P_phase = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
    E->density_phase = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
    E->Vp_phase = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
    E->Vs_phase = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  }

  set_up_nonmg_aliases(E);
  //allocate_velocity_vars(E);
  return;
}

/*  =========================================================  */

void allocate_velocity_vars(E) struct All_variables *E;

{
  int i, j, l;

  E->lmesh.nnov = E->lmesh.nno;
  E->lmesh.NEQ[E->mesh.levmax] = E->lmesh.nnov * E->mesh.nsd;

  E->F = (double *)malloc((E->mesh.nsd * E->lmesh.nnov + 1) * sizeof(double));
  E->U = (double *)malloc((E->mesh.nsd * E->lmesh.nnov + 1) * sizeof(double));

  for (i = 1; i <= E->mesh.nsd; i++)
  {
    E->V[i] = (double *)malloc((E->lmesh.nnov + 1) * sizeof(double));
    E->VB[i] = (double *)malloc((E->lmesh.nnov + 1) * sizeof(double));
  }

  for (i = E->mesh.levmin; i <= E->mesh.levmax; i++)
  {
    E->BI[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));

    E->EQN[i] = (unsigned int *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(unsigned int));
  }

  E->temp = (double *)malloc((E->lmesh.NEQ[E->mesh.levmax] + 2) * sizeof(double));
  E->Element = (unsigned int *)malloc((E->lmesh.nel + 2) * sizeof(unsigned int));

  for (l = E->mesh.levmin; l <= E->mesh.levmax; l++)
    for (i = 0; i < E->lmesh.NEQ[l]; i++)
    {
      E->BI[l][i] = 0.0;
      E->EQN[l][i] = 0;
      E->parallel.IDD[l][i] = 0;
    }

  for (i = 0; i < E->lmesh.NEQ[E->mesh.levmax]; i++)
    E->U[i] = 0.0;

  for (i = 1; i <= E->lmesh.nnov; i++)
    for (j = 1; j <= E->mesh.nsd; j++)
      E->V[j][i] =
          E->VB[j][i] = 0.0;

  return;
}

/*  =========================================================  */

void interuption()

{
  if (Emergency_stop++)
    exit(0);
  fprintf(stderr, "Cleaning up before exit\n");
  return;
}

void global_default_values(E) struct All_variables *E;
{
  FILE *fp;

  /* FIRST: values which are not changed routinely by the user */

  E->control.v_steps_low = 10;
  E->control.v_steps_upper = 1;
  E->control.max_res_red_each_p_mg = 1.0e-3;
  E->control.accuracy = 1.0e-6;
  E->control.vaccuracy = 1.0e-8;
  E->control.true_vcycle = 0;
  E->control.depth_dominated = 0;
  E->control.eqn_zigzag = 0;
  E->control.verbose = 0; /* debugging/profiles */

  /* SECOND: values for which an obvious default setting is useful */

  E->control.ORTHO = 1;  /* for orthogonal meshes by default */
  E->control.ORTHOZ = 1; /* for orthogonal meshes by default */

  E->control.KERNEL = 0;
  E->control.stokes = 0;
  E->control.CONVECTION = 0;
  E->control.CART2D = 0;
  E->control.CART3D = 0;
  E->control.CART2pt5D = 0;
  E->control.AXI = 0;
  E->control.CONJ_GRAD = 0;
  E->control.NMULTIGRID = 0;
  E->control.EMULTIGRID = 0;
  E->control.COMPRESS = 1;
  E->control.augmented_Lagr = 0;
  E->control.augmented = 0.0;

  /* Default: all optional modules set to `off' */
  E->control.MELTING_MODULE = 0;
  E->control.CHEMISTRY_MODULE = 0;

  E->control.composition = 0;
  E->control.comp_diff = 0.0;
  E->control.composition_phasechange = 0;

  E->control.GRID_TYPE = 1;
  E->mesh.hwidth[1] = E->mesh.hwidth[2] = E->mesh.hwidth[3] = 1.0; /* divide by this one ! */
  E->mesh.magnitude[1] = E->mesh.magnitude[2] = E->mesh.magnitude[3] = 0.0;
  E->mesh.offset[1] = E->mesh.offset[2] = E->mesh.offset[3] = 0.0;

  E->parallel.automa = 0;
  E->parallel.nprocx = 1;
  E->parallel.nprocz = 1;
  E->parallel.nprocy = 1;

  E->mesh.levmax = 0;
  E->mesh.levmin = 0;
  E->mesh.nox = 1;
  E->mesh.nxs = 1;
  E->lmesh.nox = 1;
  E->lmesh.nxs = 1;
  E->mesh.noz = 1;
  E->mesh.nzs = 1;
  E->lmesh.noz = 1;
  E->lmesh.nzs = 1;
  E->mesh.noy = 1;
  E->mesh.nys = 1;
  E->lmesh.noy = 1;
  E->lmesh.nys = 1;

  E->monitor.T_interior = 1.0;

  E->viscosity.guess = 0;
  sprintf(E->viscosity.old_file, "initialize");

  E->control.precondition = 0; /* for larger visc contrasts turn this back on  */
  E->control.vprecondition = 1;

  E->mesh.toptbc = 1; /* fixed t */
  E->mesh.bottbc = 1;
  E->mesh.topvbc = 0; /* stress */
  E->mesh.botvbc = 0;
  E->mesh.sidevbc = 0;
  E->mesh.periodic_x = 0; /* reflection is default*/
  E->mesh.periodic_y = 0;
  E->control.VBXtopval = 0.0;
  E->control.VBYtopval = 0.0;
  E->control.VBXbotval = 0.0;
  E->control.VBYbotval = 0.0;

  E->data.layer_km = 2800.0; /* Earth, whole mantle defaults */
  E->data.grav_acc = 9.81;
  E->data.therm_exp = 3.28e-5;
  E->data.therm_exp_factor = 1.0;
  E->data.visc_factor = 1.0;
  E->data.Cp = 1200.0;
  E->data.therm_diff = 8.0e-7;
  E->data.therm_cond = 3.168;
  E->data.density = 3340.0;
  E->data.res_density = 3295.0; /* density when X = ... */
  E->data.res_density_X = 0.3;
  E->data.melt_density = 2800.0;
  E->data.permeability = 3.0e-10;
  E->data.density_above = 1030.0; /* sea water */
  E->data.gas_const = 8.3;
  E->data.surf_heat_flux = 4.4e-2;
  E->data.grav_const = 6.673e-11;
  E->data.surf_temp = 0.0;
  E->data.youngs_mod = 1.0e11;
  E->data.Te = 0.0;
  E->data.T_sol0 = 1373.0; /* Dave's values 1991 (for the earth) */
  E->data.Tsurf = 273.0;
  E->data.dTsol_dz = 3.4e-3;
  E->data.dTsol_dF = 440.0;
  E->data.dT_dz = 0.48e-3;
  E->data.delta_S = 250.0;
  E->data.ref_temperature = 2 * 1350.0; /* fixed temperature ... delta T */

  /* THIRD: you forgot and then went home, let's see if we can help out */

  sprintf(E->control.data_file, "citcom.tmp.%d", getpid());

  E->control.NASSEMBLE = 0;

  E->mesh.layer[1] = E->mesh.layer[2] = E->mesh.layer[3] = 1.0;
  E->monitor.elapsed_time = 0.0;

  return;
}

void global_derived_values(E) struct All_variables *E;

{
  int d, lx, lz, ly, i, nox, noz, noy;

  if (E->control.NMULTIGRID || E->control.EMULTIGRID)
  {
    E->mesh.levmax = E->mesh.levels - 1;
    E->mesh.nox = E->mesh.mgunitx * (int)pow(2.0, ((double)E->mesh.levmax)) + 1;
    E->mesh.noz = E->mesh.mgunitz * (int)pow(2.0, ((double)E->mesh.levmax)) + 1;
    if (E->mesh.nsd == 3)
      E->mesh.noy = E->mesh.mgunity * (int)pow(2.0, ((double)E->mesh.levmax)) + 1;
  }

  if (E->mesh.nsd != 3)
    E->mesh.noy = 1;

  E->mesh.nnx[1] = E->mesh.nox;
  E->mesh.nnx[2] = E->mesh.noz;
  E->mesh.nnx[3] = E->mesh.noy;
  E->mesh.elx = E->mesh.nox - 1;
  E->mesh.elz = E->mesh.noz - 1;
  E->mesh.ely = max(E->mesh.noy - 1, 1);

  E->mesh.nel = E->mesh.elx * E->mesh.ely * E->mesh.elz;
  E->mesh.nno = E->mesh.nox * E->mesh.noy * E->mesh.noz;
  E->mesh.nnov = E->mesh.nno;
  E->mesh.neq = E->mesh.nnov * E->mesh.nsd;

  E->mesh.npno = E->mesh.nel;
  E->mesh.nsf = E->mesh.nox * E->mesh.noy;
  for (i = E->mesh.levmax; i >= E->mesh.levmin; i--) /* set up dimensions for different grids  */
  {
    if (E->control.NMULTIGRID || E->control.EMULTIGRID)
    {
      nox = E->mesh.mgunitx * (int)pow(2.0, (double)i) + 1;
      noz = E->mesh.mgunitz * (int)pow(2.0, (double)i) + 1;
      if (E->mesh.nsd == 3)
        noy = E->mesh.mgunity * (int)pow(2.0, (double)i) + 1;
      else
        noy = 1;
    }
    else
    {
      nox = E->mesh.nox;
      noz = E->mesh.noz;
      noy = E->mesh.noy;
    }

    E->mesh.ELX[i] = nox - 1;
    E->mesh.ELZ[i] = noz - 1;
    E->mesh.ELY[i] = max(noy - 1, 1);
    E->mesh.NNO[i] = nox * noz * noy;
    E->mesh.NEL[i] = (nox - 1) * (noz - 1) * max((noy - 1), 1);
    ;
    E->mesh.NPNO[i] = E->mesh.NEL[i];
    E->mesh.NOX[i] = nox;
    E->mesh.NOZ[i] = noz;
    E->mesh.NOY[i] = noy;
    E->mesh.NNX[i][1] = nox;
    E->mesh.NNX[i][2] = noz;
    E->mesh.NNX[i][3] = noy;

    E->mesh.NNOV[i] = E->mesh.NNO[i];
    E->mesh.NEQ[i] = E->mesh.nsd * E->mesh.NNOV[i];
  }

  if (E->control.print_convergence)
    fprintf(stderr, "Problem has %d x %d x %d nodes\n", E->mesh.nox, E->mesh.noz, E->mesh.noy);

  return;
}

void read_initial_settings(E) struct All_variables *E;

{
  void set_convection_defaults();
  void set_2dc_defaults();
  void set_3dc_defaults();
  void set_cg_defaults();
  void set_mg_defaults();
  char logfile[100];
  FILE *fp;
  int m;

  m = E->parallel.me;

  /* first the problem type (defines subsequent behaviour) */

  input_string("Problem", E->control.PROBLEM_TYPE, NULL);
  if (strcmp(E->control.PROBLEM_TYPE, "convection") == 0)
  {
    E->control.CONVECTION = 1;
    set_convection_defaults(E);
  }

  else if (strcmp(E->control.PROBLEM_TYPE, "convection-chemical") == 0)
  {
    E->control.CONVECTION = 1;
    E->control.CHEMISTRY_MODULE = 1;
    set_convection_defaults(E);
  }

  else
  {
    if (E->parallel.me == 0)
      fprintf(E->fp, "Unable to determine problem type, assuming convection ... \n");
    E->control.CONVECTION = 1;
    set_convection_defaults(E);
  }

  input_string("Geometry", E->control.GEOMETRY, NULL);
  if (strcmp(E->control.GEOMETRY, "cart2d") == 0)
  {
    E->control.CART2D = 1;
    set_2dc_defaults(E);
  }
  else if (strcmp(E->control.GEOMETRY, "axi") == 0)
  {
    E->control.AXI = 1;
  }
  else if (strcmp(E->control.GEOMETRY, "cart2pt5d") == 0)
  {
    E->control.CART2pt5D = 1;
    set_2pt5dc_defaults(E);
  }
  else if (strcmp(E->control.GEOMETRY, "cart3d") == 0)
  {
    E->control.CART3D = 1;
    set_3dc_defaults(E);
  }
  else
  {
    fprintf(E->fp, "Unable to determine geometry, assuming cartesian 2d ... \n");
    E->control.CART2D = 1;
    set_2dc_defaults(E);
  }

  input_string("Solver", E->control.SOLVER_TYPE, NULL);
  if (strcmp(E->control.SOLVER_TYPE, "cgrad") == 0)
  {
    E->control.CONJ_GRAD = 1;
    set_cg_defaults(E);
  }
  else if (strcmp(E->control.SOLVER_TYPE, "multigrid") == 0)
  {
    E->control.NMULTIGRID = 1;
    set_mg_defaults(E);
  }
  else if (strcmp(E->control.SOLVER_TYPE, "multigrid-el") == 0)
  {
    E->control.EMULTIGRID = 1;
    set_mg_defaults(E);
  }
  else
  {
    fprintf(stderr, "Unable to determine how to solve, specify Solver=VALID_OPTION \n");
    exit(0);
  }

  /* admin */

  input_string("Spacing", E->control.NODE_SPACING, "regular");
  if (strcmp(E->control.NODE_SPACING, "regular") == 0)
    E->control.GRID_TYPE = 1;
  else if (strcmp(E->control.NODE_SPACING, "bound_lyr") == 0)
    E->control.GRID_TYPE = 2;
  else if (strcmp(E->control.NODE_SPACING, "region") == 0)
    E->control.GRID_TYPE = 3;
  else if (strcmp(E->control.NODE_SPACING, "ortho_files") == 0)
    E->control.GRID_TYPE = 4;
  else
  {
    E->control.GRID_TYPE = 1;
  }

  /* Information on which files to print, which variables of the flow to calculate and print.
       Default is no information recorded (apart from special things for given applications.
     */

  input_string("datatypes", E->control.which_data_files, "");
  input_string("averages", E->control.which_horiz_averages, "");
  input_string("timelog", E->control.which_running_data, "");
  input_string("observables", E->control.which_observable_data, "");

  input_string("datafile", E->control.data_file, "initialize");
  input_string("restart_datafile", E->control.data_file1, "initialize");
  input_string("process_command", E->control.output_written_external_command, "");
  input_boolean("CONMAN", &(E->control.CONMAN), "off");

  input_boolean("imposevelo", &(E->control.imposevelo), "off");
  input_int("age_total", &(E->control.age_total), "130");
  input_string("velo_file_pre", E->control.velo_file_pre, "SurfVelo130_");
  input_string("velo_file_post", E->control.velo_file_post, ".dat");
  input_double("timescale", &(E->control.timescale), "261100.0");
  input_double("age_total_double", &(E->control.age_total_double), "130.0");

  input_boolean("ocean_lith", &(E->control.ocean_lith), "off");
  input_boolean("platemodel", &(E->control.platemodel), "off");
  input_double("inter_temp", &(E->control.inter_temp), "0.5196");
  input_double("age_left", &(E->control.age_left), "100.0");
  input_double("age_right", &(E->control.age_right), "0.0");
  input_double("age_loc_left", &(E->control.age_loc_left), "1.0");
  input_double("age_loc_right", &(E->control.age_loc_right), "2.0");
  input_double("depth_lith", &(E->control.depth_lith), "0.95644599");
  input_boolean("ocean_lith_margin", &(E->control.ocean_lith_margin), "off");
  input_double("ocean_lith_margin_curve", &(E->control.ocean_lith_margin_curve), "1.0");
  input_double("dip_center_x", &(E->control.dip_center_x), "1.0");
  input_double("dip_center_z", &(E->control.dip_center_z), "1.0");
  input_int("initialTOption", &(E->control.initialTOption), "2");
  input_int("initialCOption", &(E->control.initialCOption), "1");

  input_double("depth_lith_margin", &(E->control.depth_lith_margin), "0.9303");
  input_double("dip_margin", &(E->control.dip_margin), "45.0");
  input_double("dip_margin_left", &(E->control.dip_margin_left), "45.0");
  input_double("depth_ocean_lith", &(E->control.depth_ocean_lith), "0.95644599");

  input_boolean("continent_lith", &(E->control.continent_lith), "off");
  input_boolean("continent_platemodel", &(E->control.continent_platemodel), "off");
  input_double("continent_age", &(E->control.continent_age), "20");
  input_double("depth_continent_lith", &(E->control.depth_continent_lith), "0.95644599");
  input_double("depth_continent_crust", &(E->control.depth_continent_crust), "0.95644599");

  input_double("continent_loc_left", &(E->control.continent_loc_left), "0.0");
  input_double("continent_loc_right", &(E->control.continent_loc_right), "1.0");

  input_boolean("slab_visc", &(E->control.slab_visc), "off");
  input_double("slab_visc_depth", &(E->control.slab_visc_depth), "0.76655");

  input_double("viscincreaseslab", &(E->control.viscincreaseslab), "10.0");
  input_double("temp_slabvisc", &(E->control.temp_slabvisc), "0.05");

  input_boolean("visc_const_cor", &(E->control.visc_const_cor), "off");
  input_boolean("visc_leftcor", &(E->control.visc_leftcor), "off");
  input_boolean("visc_rightcor", &(E->control.visc_rightcor), "off");
  input_boolean("visc_mid", &(E->control.visc_mid), "off");
  input_boolean("visc_mid_dip", &(E->control.visc_mid_dip), "off");

  input_double("visc_weakzone", &(E->control.visc_weakzone), "0.001");
  input_double("x_weakzone_leftcor", &(E->control.x_weakzone_leftcor), "0.01");
  input_double("x_weakzone_rightcor", &(E->control.x_weakzone_rightcor), "0.99");
  input_double("x_weakzone_mid_left", &(E->control.x_weakzone_mid_left), "0.45");
  input_double("x_weakzone_mid_right", &(E->control.x_weakzone_mid_right), "0.55");
  input_double("dip_weakzone_mid", &(E->control.dip_weakzone_mid), "30.0");

  input_double("z_weakzone_leftcor", &(E->control.z_weakzone_leftcor), "0.97");
  input_double("z_weakzone_rightcor", &(E->control.z_weakzone_rightcor), "0.97");
  input_double("z_weakzone_mid", &(E->control.z_weakzone_mid), "0.97");

  input_boolean("CBF", &(E->control.CBF), "off");
  input_boolean("phasevisc", &(E->control.phasevisc), "off");
  input_double("phaseTop", &(E->control.phaseTop), "0.0871");
  input_double("phaseBot", &(E->control.phaseBot), "0.0871");
  input_double("phaseffactorLa", &(E->control.phaseffactorLa), "0.95");
  input_double("phaseffactorSm", &(E->control.phaseffactorSm), "0.05");
  input_boolean("phasevisc_slab", &(E->control.phasevisc_slab), "off");
  input_double("ViscReduce", &(E->control.ViscReduce), "0.01");
  input_double("phasevisc_slab_T", &(E->control.phasevisc_slab_T), "-0.01");

  input_boolean("phasevisc_C", &(E->control.phasevisc_C), "off");
  input_double("phasevisc_C_time", &(E->control.phasevisc_C_time), "1e6");

  input_boolean("phasevisc_d", &(E->control.phasevisc_d), "off");
  input_double("phasevisc_d0", &(E->control.phasevisc_d0), "1.0");
  input_double("phasevisc_dp", &(E->control.phasevisc_dp), "2.0");
  input_double("phasevisc_dm", &(E->control.phasevisc_dm), "3.0");
  input_double("phasevisc_d_phasereduce", &(E->control.phasevisc_d_phasereduce), "0.03333");
  input_double("phasevisc_E", &(E->control.phasevisc_E), "9.21");
  input_double("phasevisc_tratio", &(E->control.phasevisc_tratio), "1.0000");
  input_double("phasevisc_dmax_value", &(E->control.phasevisc_dmax_value), "1.0");
  input_boolean("phasevisc_dmax", &(E->control.phasevisc_dmax), "off");
  input_boolean("lesstimeinte", &(E->control.lesstimeinte), "off");
  input_double("lesstimeinte_number", &(E->control.lesstimeinte_number), "essential,1.0000");

  input_boolean("Visc_C", &(E->control.Visc_C), "off");
  input_double("ViscReduce_C", &(E->control.ViscReduce_C), "0.01");

  input_double("velo_surf_loc_mid", &(E->control.velo_surf_loc_mid), "1.0");
  input_double("velo_surf_mag_right", &(E->control.velo_surf_mag_right), "-1000.0");
  input_double("velo_surf_width_right", &(E->control.velo_surf_width_right), "50.0");
  input_boolean("trechmigrate", &(E->control.trechmigrate), "off");
  input_double("velo_surf_mag_left", &(E->control.velo_surf_mag_left), "-1000.0");
  input_double("velo_surf_width_left", &(E->control.velo_surf_width_left), "50.0");
  input_double("velo_surf_width_mid", &(E->control.velo_surf_width_mid), "50.0");
  input_double("velo_surf_corner_right", &(E->control.velo_surf_corner_right), "50.0");
  input_double("velo_surf_loc_mid_rate", &(E->control.velo_surf_loc_mid_rate), "0.003484");
  input_double("velo_surf_loc_left_overshoot", &(E->control.velo_surf_loc_left_overshoot), "0.02998955");

  /* As early as possible, set up the log file to 
       record information about the progress of the 
       program as it runs 
       */

  sprintf(logfile, "%s.log", E->control.data_file);

  if ((fp = fopen(logfile, "w")) == NULL)
    E->fp = stdout;
  else
    E->fp = fp;

  if (E->control.NMULTIGRID || E->control.EMULTIGRID)
  {
    input_int("mgunitx", &(E->mesh.mgunitx), "1");
    input_int("mgunitz", &(E->mesh.mgunitz), "1");
    input_int("mgunity", &(E->mesh.mgunity), "1");
    input_int("levels", &(E->mesh.levels), "0");
  }

  input_boolean("node_assemble", &(E->control.NASSEMBLE), "off");
  /* general mesh structure */

  input_boolean("parallel_auto", &(E->parallel.automa), "off", m);
  if (E->parallel.automa == 0)
  {
    input_int("nprocx", &(E->parallel.nprocx), "1", m);
    input_int("nprocz", &(E->parallel.nprocz), "1", m);
    input_int("nprocy", &(E->parallel.nprocy), "1", m);
  }

  input_boolean("verbose", &(E->control.verbose), "off");
  input_boolean("see_convergence", &(E->control.print_convergence), "off");
  input_boolean("COMPRESS", &(E->control.COMPRESS), "on");
  input_double("sobtol", &(E->control.sob_tolerance), "0.0001");

  input_int("obs_maxlongk", &(E->slice.maxlong), "100,1");
  input_int("obs_minlongk", &(E->slice.minlong), "1,1");

  input_int("stokes_flow_only", &(E->control.stokes), "0");

  input_int("slab_nz", &(E->control.SLAB), "1");
  input_int("plume_nz", &(E->control.PLUME), "1");

  /* for phase change    */

  input_double("Ra_670", &(E->control.Ra_670), "0.0");
  input_double("clapeyron670", &(E->control.clapeyron670), "0.0");
  input_double("transT670", &(E->control.transT670), "0.0");
  input_double("width670", &(E->control.width670), "0.0");

  input_double("Ra_670_basalt", &(E->control.Ra_670_basalt), "0.0");
  input_double("clapeyron670_basalt", &(E->control.clapeyron670_basalt), "0.0");
  input_double("transT670_basalt", &(E->control.transT670_basalt), "0.0");
  input_double("width670_basalt", &(E->control.width670_basalt), "0.0");

  input_double("Ra_410", &(E->control.Ra_410), "0.0");
  input_double("clapeyron410", &(E->control.clapeyron410), "0.0");
  input_double("transT410", &(E->control.transT410), "0.0");
  input_double("width410", &(E->control.width410), "0.0");

  input_boolean("temperature_perturbation", &(E->control.temperature_perturbation), "off");

  input_boolean("meshReadIn", &(E->control.meshReadIn), "off");
  if (E->control.meshReadIn)
  {
    input_string("fileMeshX", E->control.fileMeshX, "MeshX.txt");
    input_string("fileMeshZ", E->control.fileMeshZ, "MeshZ.txt");
  }

  input_boolean("tracer_correction", &(E->control.tracer_correction), "off");
  input_boolean("phasefile", &(E->control.phasefile), "off");
  input_boolean("phasefile_Complete", &(E->control.phasefile_Complete), "off");
  input_boolean("phasefile_buoyancy", &(E->control.phasefile_buoyancy), "off");
  input_boolean("phasefile_buoyancy_correction", &(E->control.phasefile_buoyancy_correction), "off");
  input_double("phasefile_buoyancy_depth", &(E->control.phasefile_buoyancy_depth), "0.0");
  input_double("phasefile_buoyancy_continent", &(E->control.phasefile_buoyancy_continent), "0.0");
  input_double("phasefile_buoyancy_crust", &(E->control.phasefile_buoyancy_crust), "0.0");
  input_double("phasefile_buoyancy_crust_depth", &(E->control.phasefile_buoyancy_crust_depth), "300.0");
  input_double("adi_um", &(E->control.adi_um), "0.5");
  input_double("adi_lm", &(E->control.adi_lm), "0.3");

  input_string("phasefile_basa", E->control.phasefile_basa, "basalt_less_all.dat");
  input_string("phasefile_pyro", E->control.phasefile_pyro, "pyrolite_less_all.dat");
  input_string("phasefile_harz", E->control.phasefile_harz, "harzburgite_less_all.dat");
  input_int("phasefile_noP", &(E->control.phasefile_noP), "313");
  input_int("phasefile_noT", &(E->control.phasefile_noT), "313");
  input_int("phasefile_noPREM", &(E->control.phasefile_noPREM), "199");
  input_string("phasefile_PREM", E->control.phasefile_PREM, "PREM.txt");
  input_double("depth_harz", &(E->control.depth_harz), "0.010453");
  input_boolean("phasefile_C", &(E->control.phasefile_C), "off");
  /*input_int("phasefile_C_num_nno", &(E->control.phasefile_C_num_nno), "3");
  input_int("phasefile_C_num_element", &(E->control.phasefile_C_num_element), "1");
  input_int("phasefile_C_num_marker", &(E->control.phasefile_C_num_marker), "1");*/
  input_int("phasefile_C_flavor", &(E->control.phasefile_C_flavor), "3");
  E->control.phasefile_C_num_nno = E->control.phasefile_C_flavor;
  E->control.phasefile_C_num_element = E->control.phasefile_C_flavor;
  E->control.phasefile_C_num_marker = E->control.phasefile_C_flavor;
  input_int_vector("phasefile_C_mat_mineral", E->control.phasefile_C_flavor, (E->control.phasefile_C_mat_mineral));

  input_boolean("crust_generation", &(E->control.crust_generation), "off");
  input_boolean("harz_generation", &(E->control.harz_generation), "off");

  input_int("restart", &(E->control.restart), "0");

  input_int("topvbc", &(E->mesh.topvbc), "0");
  input_int("botvbc", &(E->mesh.botvbc), "0");
  input_int("sidevbc", &(E->mesh.sidevbc), "0");

  input_boolean("periodicx", &(E->mesh.periodic_x), "off");
  input_boolean("periodicy", &(E->mesh.periodic_y), "off");
  input_boolean("depthdominated", &(E->control.depth_dominated), "off");
  input_boolean("eqnzigzag", &(E->control.eqn_zigzag), "off");
  input_boolean("eqnviscosity", &(E->control.eqn_viscosity), "off");

  input_double("topvbxval", &(E->control.VBXtopval), "0.0");
  input_double("botvbxval", &(E->control.VBXbotval), "0.0");
  input_double("topvbyval", &(E->control.VBYtopval), "0.0");
  input_double("botvbyval", &(E->control.VBYbotval), "0.0");

  input_int("toptbc", &(E->mesh.toptbc), "1");
  input_int("bottbc", &(E->mesh.bottbc), "1");
  input_double("toptbcval", &(E->control.TBCtopval), "0.0");
  input_double("bottbcval", &(E->control.TBCbotval), "1.0");

  input_double("blyr_hwx1", &(E->mesh.bl1width[1]), "nodefault");
  input_double("blyr_hwz1", &(E->mesh.bl1width[2]), "nodefault");
  input_double("blyr_hwy1", &(E->mesh.bl1width[3]), "nodefault");
  input_double("blyr_hwx2", &(E->mesh.bl2width[1]), "nodefault");
  input_double("blyr_hwz2", &(E->mesh.bl2width[2]), "nodefault");
  input_double("blyr_hwy2", &(E->mesh.bl2width[3]), "nodefault");
  input_double("blyr_mgx1", &(E->mesh.bl1mag[1]), "nodefault");
  input_double("blyr_mgz1", &(E->mesh.bl1mag[2]), "nodefault");
  input_double("blyr_mgy1", &(E->mesh.bl1mag[3]), "nodefault");
  input_double("blyr_mgx2", &(E->mesh.bl2mag[1]), "nodefault");
  input_double("blyr_mgz2", &(E->mesh.bl2mag[2]), "nodefault");
  input_double("blyr_mgy2", &(E->mesh.bl2mag[3]), "nodefault");

  input_double("region_wdx", &(E->mesh.width[1]), "nodefault");
  input_double("region_wdz", &(E->mesh.width[2]), "nodefault");
  input_double("region_wdy", &(E->mesh.width[3]), "nodefault");
  input_double("region_hwx", &(E->mesh.hwidth[1]), "nodefault");
  input_double("region_hwz", &(E->mesh.hwidth[2]), "nodefault");
  input_double("region_hwy", &(E->mesh.hwidth[3]), "nodefault");
  input_double("region_mgx", &(E->mesh.magnitude[1]), "nodefault");
  input_double("region_mgz", &(E->mesh.magnitude[2]), "nodefault");
  input_double("region_mgy", &(E->mesh.magnitude[3]), "nodefault");
  input_double("region_ofx", &(E->mesh.offset[1]), "nodefault");
  input_double("region_ofz", &(E->mesh.offset[2]), "nodefault");
  input_double("region_ofy", &(E->mesh.offset[3]), "nodefault");

  input_string("gridxfile", E->mesh.gridfile[1], " ");
  input_string("gridzfile", E->mesh.gridfile[2], " ");
  input_string("gridyfile", E->mesh.gridfile[3], " ");

  input_double("layerd", &(E->data.layer_km), "2800.0");

  E->data.layer_meter = E->data.layer_km * 1000;
  /* for layers    */
  E->viscosity.zlm = 0.76655;
  E->viscosity.zlith = 0.95644599;
  input_int("nz_lmantle", &(E->viscosity.nlm), "1");
  input_int("nz_410", &(E->viscosity.n410), "1");
  input_int("nz_lith", &(E->viscosity.nlith), "1");
  input_int("nz_mid_moho", &(E->viscosity.ncrust2), "1");
  input_double("z_lmantle", &(E->viscosity.zlm), "0.76655");
  input_double("z_410", &(E->viscosity.z410), "0.85714");
  input_double("z_lith", &(E->viscosity.zlith), "0.95644599");
  input_double("z_mid_moho", &(E->viscosity.zcrust2), "0.0");
  input_double("z_300", &(E->viscosity.z300), "0.89547");
  input_double("z_1000", &(E->viscosity.z1000), "0.6515679");
  input_double("z_basalt", &(E->viscosity.zbasalt), "0.76655");

  input_double("z_comp", &(E->viscosity.zcrust1), "0.0");

  input_double("dimenx", &(E->mesh.layer[1]), "nodefault");
  input_double("dimenz", &(E->mesh.layer[2]), "nodefault");

  input_int("nodex", &(E->mesh.nox), "nodefault,1,nomax");
  input_int("nodez", &(E->mesh.noz), "nodefault,1,nomax");
  input_int("nodey", &(E->mesh.noy), "1,1,nomax");
  input_boolean("aug_lagr", &(E->control.augmented_Lagr), "off");
  input_double("aug_number", &(E->control.augmented), "0.0");

  input_double("tole_compressibility", &(E->control.tole_comp), "0.0");
  input_boolean("orthogonal", &(E->control.ORTHO), "on");

  input_int("storage_spacing", &(E->control.record_every), "10");
  input_int("storage_always_before", &(E->control.record_all_until), "5");

  input_boolean("precond", &(E->control.precondition), "off");
  input_boolean("vprecond", &(E->control.vprecondition), "on");
  input_int("mg_cycle", &(E->control.mg_cycle), "2,0,nomax");
  input_int("down_heavy", &(E->control.down_heavy), "1,0,nomax");
  input_int("up_heavy", &(E->control.up_heavy), "1,0,nomax");
  input_double("accuracy", &(E->control.accuracy), "1.0e-4,0.0,1.0");
  if (E->control.NMULTIGRID || E->control.EMULTIGRID)
  { /* multigrid: set previous default */
    E->control.relaccCG = 1e-3;
    input_double("relaccCG", &(E->control.relaccCG), "1.0e-3,0.0,1.0");
  }
  else
  { /* conj.grad.: set previous default */
    E->control.relaccCG = 1e0;
    input_double("relaccCG", &(E->control.relaccCG), "1.0e0,0.0,1.0");
  }

  E->control.relaccMG = 1e0;
  input_double("relaccMG", &(E->control.relaccMG), "1.0e0,0.0,1.0");
  fprintf(stderr, "relaccMG = %g\n", E->control.relaccMG);

  input_int("viterations", &(E->control.max_vel_iterations), "250,0,nomax");
  input_int("viterations_min", &(E->control.min_vel_iterations), "2,0,nomax");

  input_int("vhighstep", &(E->control.v_steps_high), "1,0,nomax");
  input_int("vlowstep", &(E->control.v_steps_low), "250,0,nomax");
  input_int("vupperstep", &(E->control.v_steps_upper), "1,0,nomax");
  input_int("piterations", &(E->control.p_iterations), "100,0,nomax");
  input_int("maxsamevisc", &(E->control.max_same_visc), "25,0,nomax");

  input_int("llmax", &(E->control.llmax), "50,0,nomax");

  /* data section */

  input_double("Ts", &(E->control.Ts), "0.0");
  input_double("ReferenceT", &(E->data.ref_temperature), "2600.0");

  E->data.ref_temperature = E->data.ref_temperature - E->control.Ts;
  E->control.Ts = E->control.Ts / E->data.ref_temperature;

  E->rad_heat.num = 0;

  input_int("num_radioactives", &(E->rad_heat.num), "0");
  fprintf(E->fp, "n_rad %d\n", E->rad_heat.num);

  if (E->rad_heat.num != 0)
  {
    input_double("concen_u", &(E->rad_heat.concen_u), "0.0");
    input_double_vector("percent", E->rad_heat.num, (E->rad_heat.percent));
    input_double_vector("heat_g", E->rad_heat.num, (E->rad_heat.heat_g));
    input_double_vector("decay_time", E->rad_heat.num, (E->rad_heat.decay_t));
    E->rad_heat.concen[0] = E->rad_heat.concen_u;
    E->rad_heat.concen[1] = E->rad_heat.concen_u;
    E->rad_heat.concen[2] = E->rad_heat.concen_u * 4;
    E->rad_heat.concen[3] = E->rad_heat.concen_u * 10000;
    fprintf(E->fp, "Rad_heat %.4e %.4e %.4e %.4e\n", E->rad_heat.percent[0], E->rad_heat.heat_g[0], E->rad_heat.decay_t[0], E->rad_heat.concen[0]);
    fprintf(E->fp, "Rad_heat %.4e %.4e %.4e %.4e\n", E->rad_heat.percent[1], E->rad_heat.heat_g[1], E->rad_heat.decay_t[1], E->rad_heat.concen[1]);
    fprintf(E->fp, "Rad_heat %.4e %.4e %.4e %.4e\n", E->rad_heat.percent[2], E->rad_heat.heat_g[2], E->rad_heat.decay_t[2], E->rad_heat.concen[2]);
    fprintf(E->fp, "Rad_heat %.4e %.4e %.4e %.4e\n", E->rad_heat.percent[3], E->rad_heat.heat_g[3], E->rad_heat.decay_t[3], E->rad_heat.concen[3]);
    fflush(E->fp);
  }
  else
    input_double("Q0", &(E->control.Q0), "0.0");

  input_double("gravacc", &(E->data.grav_acc), "9.81");
  input_double("thermexp", &(E->data.therm_exp), "3.28e-5");
  input_double("thermexp_factor", &(E->data.therm_exp_factor), "1");
  input_double("visc_factor", &(E->data.visc_factor), "3.28e-5");
  input_double("cp", &(E->data.Cp), "1200.0");
  input_double("thermdiff", &(E->data.therm_diff), "8.0e-7");
  input_double("thermcond", &(E->data.therm_cond), "3.168");
  input_double("density", &(E->data.density), "3340.0");
  input_double("mdensity", &(E->data.melt_density), "2800.0");
  input_double("density_above", &(E->data.density_above), "1030.0");
  input_double("density_below", &(E->data.density_below), "5400.0");
  input_double("rdensity", &(E->data.res_density), "3295.0");
  input_double("heatflux", &(E->data.surf_heat_flux), "4.4e-2");
  input_double("refvisc", &(E->data.ref_viscosity), "nodefault");
  input_double("meltvisc", &(E->data.melt_viscosity), "nodefault");
  input_double("surftemp", &(E->data.surf_temp), "0.0");
  input_double("youngs", &(E->data.youngs_mod), "1.0e11");
  input_double("Te", &(E->data.Te), "0.0");
  input_double("Tsol0", &(E->data.T_sol0), "1373.0");
  input_double("dTsoldz", &(E->data.dTsol_dz), "3.4e-3");
  input_double("dTsoldF", &(E->data.dTsol_dF), "440.0");
  input_double("dTdz", &(E->data.dT_dz), "0.48e-3");
  input_double("deltaS", &(E->data.delta_S), "250.0");
  input_double("gasconst", &(E->data.gas_const), "8.3"); /* not much cause to change these ! */
  input_double("gravconst", &(E->data.grav_const), "6.673e-11");
  input_double("permeability", &(E->data.permeability), "3.0e-10");

  E->monitor.time_scale = E->data.layer_meter * E->data.layer_meter / (E->data.therm_diff * 3600.0 * 24.0 * 365.25); /* years*/

  (E->problem_settings)(E);

  return;
}

void check_bc_consistency(E) struct All_variables *E;

{
  int i, lev;

  for (i = 1; i <= E->lmesh.nno; i++)
  {
    if ((E->node[i] & VBX) && (E->node[i] & SBX))
      printf("Inconsistent x velocity bc at %d\n", i);
    if ((E->node[i] & VBZ) && (E->node[i] & SBZ))
      printf("Inconsistent z velocity bc at %d\n", i);
    if ((E->node[i] & VBY) && (E->node[i] & SBY))
      printf("Inconsistent y velocity bc at %d\n", i);
    if ((E->node[i] & TBX) && (E->node[i] & FBX))
      printf("Inconsistent x temperature bc at %d\n", i);
    if ((E->node[i] & TBZ) && (E->node[i] & FBZ))
      printf("Inconsistent z temperature bc at %d\n", i);
    if ((E->node[i] & TBY) && (E->node[i] & FBY))
      printf("Inconsistent y temperature bc at %d\n", i);
  }

  for (lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
    for (i = 1; i <= E->lmesh.NNO[lev]; i++)
    {
      if ((E->NODE[lev][i] & VBX) && (E->NODE[lev][i] & SBX))
        printf("Inconsistent x velocity bc at %d,%d\n", lev, i);
      if ((E->NODE[lev][i] & VBZ) && (E->NODE[lev][i] & SBZ))
        printf("Inconsistent z velocity bc at %d,%d\n", lev, i);
      if ((E->NODE[lev][i] & VBY) && (E->NODE[lev][i] & SBY))
        printf("Inconsistent y velocity bc at %d,%d\n", lev, i);
  /* Tbc's not applicable below top level */ }

  return;
}

void set_up_nonmg_aliases(E) struct All_variables *E;

{ /* Aliases for functions only interested in the highest mg level */
  int i;

  E->eco = E->ECO[E->mesh.levmax];
  E->ien = E->IEN[E->mesh.levmax];
  E->id = E->ID[E->mesh.levmax];
  E->lm = E->LMD[E->mesh.levmax];
  E->Vi = E->VI[E->mesh.levmax];
  E->EVi = E->EVI[E->mesh.levmax];
  E->node = E->NODE[E->mesh.levmax];
  E->tw = E->TW[E->mesh.levmax];
  E->Mass = E->MASS[E->mesh.levmax];
  E->gDA = E->GDA[E->mesh.levmax];
  E->gNX = E->GNX[E->mesh.levmax];
  for (i = 1; i <= E->mesh.nsd; i++)
    E->X[i] = E->XX[E->mesh.levmax][i];

  return;
}

report(E, string) struct All_variables *E;
char *string;
{
  if (E->control.verbose && E->parallel.me == 0)
  {
    fprintf(stderr, "%s\n", string);
    fflush(stderr);
  }
  return;
}

record(E, string) struct All_variables *E;
char *string;
{
  if (E->control.verbose)
  {
    fprintf(E->fp, "%s\n", string);
    fflush(E->fp);
  }

  return;
}

/* =============================================================
   Initialize values which are not problem dependent.
   NOTE: viscosity may be a function of all previous
   input fields (temperature, pressure, velocity, chemistry) and 
   so is always to be done last.
   ============================================================= */
void common_initial_fields(E) struct All_variables *E;
{
  void initial_pressure();
  void initial_velocity();
  void read_viscosity_option();

  report(E, "Initialize pressure field");
  initial_pressure(E);
  report(E, "Initialize velocity field");
  initial_velocity(E);
  report(E, "Initialize viscosity field");
  get_viscosity_option(E);
  return;
}

/* ========================================== */

void initial_pressure(E) struct All_variables *E;
{
  int i, node, ii;

  if (!E->control.restart)
    for (i = 1; i <= E->lmesh.npno; i++)
      E->P[i] = 0.0;

  return;
}

void initial_velocity(E) struct All_variables *E;
{
  int i, node, ii;

  for (i = 1; i <= E->lmesh.nnov; i++)
  {
    E->V[1][i] = 0.0;
    E->V[2][i] = 0.0;
  }

  if (E->mesh.dof == 3)
  {
    for (i = 1; i <= E->lmesh.nnov; i++)
      E->V[3][i] = 0.0;
  }

  return;
}

void initial_phasefile(E) struct All_variables *E;
{
  fprintf(stderr, "input phasefile\n");
  int m, n, num;
  char tempstring[255], input_s[255];
  FILE *fp_read;
  int keep_going;
  double temp_P, temp_T, temp_den, temp_Vp, temp_Vs, temp1, temp2, temp3, temp4, temp5, temp6;
  E->control.phase_basa_density = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));
  E->control.phase_basa_Vp = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));
  E->control.phase_basa_Vs = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));

  E->control.phase_pyro_density = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));
  E->control.phase_pyro_Vp = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));
  E->control.phase_pyro_Vs = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));

  E->control.phase_harz_density = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));
  E->control.phase_harz_Vp = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));
  E->control.phase_harz_Vs = (double *)malloc((E->control.phasefile_noP * E->control.phasefile_noT + 1) * sizeof(double));

  E->control.phase_P = (double *)malloc((E->control.phasefile_noP + 1) * sizeof(double));
  E->control.phase_T = (double *)malloc((E->control.phasefile_noT + 1) * sizeof(double));

  E->control.phase_PREM_depth = (double *)malloc((E->control.phasefile_noPREM + 1) * sizeof(double));
  E->control.phase_PREM_density = (double *)malloc((E->control.phasefile_noPREM + 1) * sizeof(double));
  E->control.phase_PREM_Vp = (double *)malloc((E->control.phasefile_noPREM + 1) * sizeof(double));
  E->control.phase_PREM_Vs = (double *)malloc((E->control.phasefile_noPREM + 1) * sizeof(double));
  E->control.phase_PREM_P = (double *)malloc((E->control.phasefile_noPREM + 1) * sizeof(double));

  fprintf(stderr, "phasefile malloc\n");
  sprintf(tempstring, "%s", E->control.phasefile_basa);
  fp_read = fopen(tempstring, "r");
  if (fp_read == NULL)
    fprintf(stderr, "cannot open phase file %s\n", tempstring);
  for (m = 1; m <= E->control.phasefile_noT; m++)
  {
    for (n = 1; n <= E->control.phasefile_noP; n++)
    {
      num = n + E->control.phasefile_noP * (m - 1);
      fgets(input_s, 200, fp_read);
      sscanf(input_s, "%lf %lf %lf %lf %lf", &temp_P, &temp_T, &E->control.phase_basa_density[num], &E->control.phase_basa_Vp[num], &E->control.phase_basa_Vs[num]);
      if (m == 1)
      {
        E->control.phase_P[n] = temp_P / 10000;
      }
      if (n == 1)
      {
        E->control.phase_T[m] = temp_T - 273.15;
      }
    }
  }
  fclose(fp_read);
  fprintf(stderr, "basa finished\n");
  sprintf(tempstring, "%s", E->control.phasefile_pyro);
  fprintf(stderr, "start open phase file %s\n", tempstring);
  fp_read = fopen(tempstring, "r");
  fprintf(stderr, "start open phase file %s\n", tempstring);
  if (fp_read == NULL)
    fprintf(stderr, "cannot open phase file %s\n", tempstring);
  for (m = 1; m <= E->control.phasefile_noT; m++)
  {
    for (n = 1; n <= E->control.phasefile_noP; n++)
    {
      num = n + E->control.phasefile_noP * (m - 1);
      fgets(input_s, 200, fp_read);
      sscanf(input_s, "%lf %lf %lf %lf %lf", &temp_P, &temp_T, &E->control.phase_pyro_density[num], &E->control.phase_pyro_Vp[num], &E->control.phase_pyro_Vs[num]);
    }
  }
  fclose(fp_read);

  sprintf(tempstring, "%s", E->control.phasefile_harz);
  fp_read = fopen(tempstring, "r");
  if (fp_read == NULL)
    fprintf(stderr, "cannot open phase file %s\n", tempstring);
  for (m = 1; m <= E->control.phasefile_noT; m++)
  {
    for (n = 1; n <= E->control.phasefile_noP; n++)
    {
      num = n + E->control.phasefile_noP * (m - 1);
      fgets(input_s, 200, fp_read);
      sscanf(input_s, "%lf %lf %lf %lf %lf", &temp_P, &temp_T, &E->control.phase_harz_density[num], &E->control.phase_harz_Vp[num], &E->control.phase_harz_Vs[num]);
    }
  }
  fclose(fp_read);

  sprintf(tempstring, "%s", E->control.phasefile_PREM);
  fp_read = fopen(tempstring, "r");
  if (fp_read == NULL)
    fprintf(stderr, "cannot open phase file %s\n", tempstring);
  for (m = 1; m <= E->control.phasefile_noPREM; m++)
  {
    fgets(input_s, 200, fp_read);
    sscanf(input_s, "%lf %lf %lf %lf", &E->control.phase_PREM_depth[m], &E->control.phase_PREM_density[m], &E->control.phase_PREM_Vp[m], &E->control.phase_PREM_Vs[m]);
  }
  fclose(fp_read);

  E->control.phase_PREM_P[1] = 1 / 1000; // in GPa
  for (m = 2; m <= E->control.phasefile_noPREM; m++)
  {
    E->control.phase_PREM_P[m] = E->control.phase_PREM_P[m - 1] + (E->control.phase_PREM_density[m] + E->control.phase_PREM_density[m - 1]) / 2.0 * 1e3 * 9.8 * (E->control.phase_PREM_depth[m] - E->control.phase_PREM_depth[m - 1]) * 1e3 / 1e9;
    fprintf(stderr, "PREM_P %lf %lf\n", E->control.phase_PREM_depth[m], E->control.phase_PREM_P[m]);
  }

  return;
}
