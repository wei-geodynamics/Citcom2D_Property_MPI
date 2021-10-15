/*  Here are the routines which process the results of each velocity solution, and call
    the relevant output routines. At this point, the velocity and pressure fields have
    been calculated and stored at the nodes. The only properties of the velocity field
    which are already known are those required to check convergence of the iterative
    scheme and so on. */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <stdlib.h> /* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"

void process_new_velocity(E, ii) struct All_variables *E;
int ii;
{
  void output_velo_related();
  void get_STD_topo();
  void get_CBF_topo();
  void averages();
  static FILE *fp;
  char output_file[255];
  int i;
  static double *temp, vrms, cc1;
  double return_bulk_value();
  static int been_here = 0;

  if (been_here == 0)
  {
    E->monitor.length_scale = E->data.layer_km / E->mesh.layer[2]; /* km */
    E->monitor.time_scale = pow(E->data.layer_km * 1000.0, 2.0) /  /* Million years */
                            (E->data.therm_diff * 3600.0 * 24.0 * 365.25 * 1.0e6);
    sprintf(output_file, "%s/vrms_time", E->control.data_file);
    fp = fopen(output_file, "w");
    been_here++;
  }

  if (E->control.stokes || ((ii % E->control.record_every) == 0))
  {

    /*
      get_CBF_topo(E,E->slice.tpg,E->slice.tpgb);  
*/

    get_STD_topo(E, E->slice.tpg, E->slice.tpgb, ii);

    averages(E);

    output_velo_related(E, ii); /* also topo */
  }
  if ((ii % 1) == 0)
  {
    cc1 = return_bulk_value(E, E->C, 0.2, 1);
    fprintf(fp, "%.4e %.4e %.4e\n", E->monitor.elapsed_time, E->monitor.Vrms, cc1);
  }
  fflush(fp);

  return;
}

/* ===============================================   */

void get_surface_velo(E, SV) struct All_variables *E;
double *SV;
{

  int el, els, i, m, node, lev;
  char output_file[255];
  FILE *fp;

  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];
  const int nno = E->lmesh.nno;

  lev = E->mesh.levmax;

  m = 0;

  for (node = 1; node <= nno; node++)
    if ((node - 1) % E->lmesh.noz == 0)
    {
      i = (node - 1) / E->lmesh.noz + 1;
      SV[(i - 1) * 2 + 1] = E->V[1][node];
      SV[(i - 1) * 2 + 2] = E->V[3][node];
    }

  return;
}

/* ===============================================   */

void get_ele_visc(E, EV) struct All_variables *E;
double *EV;
{

  int el, j, lev;

  const int nel = E->lmesh.nel;
  const int vpts = vpoints[E->mesh.nsd];

  lev = E->mesh.levmax;

  for (el = 1; el <= nel; el++)
  {
    EV[el] = 0.0;
    for (j = 1; j <= vpts; j++)
      EV[el] += E->EVI[lev][(el - 1) * vpts + j];

    EV[el] /= vpts;
  }

  return;
}

void get_surf_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY) struct All_variables *E;
double *SXX, *SYY, *SZZ, *SXY, *SXZ, *SZY;
{
  int i, node, stride;

  stride = E->lmesh.nsf * 6;

  for (node = 1; node <= E->lmesh.nno; node++)
    if (((node - 1) % E->lmesh.noz) == 0)
    {
      i = (node - 1) / E->lmesh.noz + 1;
      E->stress[(i - 1) * 6 + 1] = SXX[node];
      E->stress[(i - 1) * 6 + 2] = SZZ[node];
      E->stress[(i - 1) * 6 + 3] = SYY[node];
      E->stress[(i - 1) * 6 + 4] = SXY[node];
      E->stress[(i - 1) * 6 + 5] = SXZ[node];
      E->stress[(i - 1) * 6 + 6] = SZY[node];
    }
    else if (((node - 2) % E->lmesh.noz) == 0)
    {
      i = (node - 2) / E->lmesh.noz + 1;
      E->stress[stride + (i - 1) * 6 + 1] = SXX[node];
      E->stress[stride + (i - 1) * 6 + 2] = SZZ[node];
      E->stress[stride + (i - 1) * 6 + 3] = SYY[node];
      E->stress[stride + (i - 1) * 6 + 4] = SXY[node];
      E->stress[stride + (i - 1) * 6 + 5] = SXZ[node];
      E->stress[stride + (i - 1) * 6 + 6] = SZY[node];
    }

  return;
}

void averages(E) struct All_variables *E;
{

  int lev, i, j, el;
  double *temp, vrmssqr;
  void return_horiz_ave();
  void plume_buoyancy_flux();
  double return_bulk_value();

  lev = E->mesh.levmax;

  temp = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));

  return_horiz_ave(E, E->T, E->Have.T);
  visc_from_gint_to_nodes(E, E->EVI[lev], temp, lev);

  return_horiz_ave(E, temp, E->Have.Vi);

  if (E->mesh.nsd == 2)
    for (i = 1; i <= E->lmesh.nno; i++)
    {
      temp[i] = E->V[1][i] * E->V[1][i];
    }
  else
    for (i = 1; i <= E->lmesh.nno; i++)
    {
      temp[i] = E->V[1][i] * E->V[1][i] + E->V[3][i] * E->V[3][i];
    }

  return_horiz_ave(E, temp, E->Have.vrms);

  for (i = 1; i <= E->lmesh.noz; i++)
    E->Have.vrms[i] = sqrt(E->Have.vrms[i]);

  //plume_buoyancy_flux(E);

  if (E->mesh.nsd == 2)
    for (i = 1; i <= E->lmesh.nno; i++)
      temp[i] = E->V[1][i] * E->V[1][i] + E->V[2][i] * E->V[2][i];
  else if (E->mesh.nsd == 3)
    for (i = 1; i <= E->lmesh.nno; i++)
      temp[i] = E->V[1][i] * E->V[1][i] + E->V[2][i] * E->V[2][i];
  +E->V[3][i] * E->V[3][i];

  vrmssqr = return_bulk_value(E, temp, -0.1, 1);
  E->monitor.Vrms = sqrt(vrmssqr);

  free((void *)temp);

  return;
}
