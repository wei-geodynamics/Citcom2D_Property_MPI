/*  Here are the routines which process the results of each buoyancy solution, and call
    any relevant output routines. Much of the information has probably been output along
    with the velocity field. (So the velocity vectors and other data are fully in sync).
    However, heat fluxes and temperature averages are calculated here (even when they
    get output the next time around the velocity solver);
    */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <stdlib.h> /* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"

void process_temp_field(E, ii) struct All_variables *E;
int ii;
{
  void heat_flux();
  void output_temp();

  if (((ii % E->control.record_every) == 0))
  {
    heat_flux(E);
    /*
      output_temp(E,ii);
*/
  }

  return;
}
/* ===================
    Surface heat flux  
   =================== */

void heat_flux(E) struct All_variables *E;
{
  int e, i, j, node, lnode;
  double *mass, *flux, *SU, *RU, *inp, *outp;
  double VZ[9], u[9], T[9], dTdz[9], area, uT;
  double tempb, tempt, hfb, hft, areab, areat;
  void return_horiz_sum();

  struct Shape_function GN;
  struct Shape_function_dA dOmega;
  struct Shape_function_dx GNx;

  const int dims = E->mesh.nsd, dofs = E->mesh.dof;
  const int vpts = vpoints[dims];
  const int ppts = ppoints[dims];
  const int ends = enodes[dims];
  const int nno = E->lmesh.nno;
  const int lev = E->mesh.levmax;

  mass = (double *)malloc((1 + nno) * sizeof(double));
  flux = (double *)malloc((1 + nno) * sizeof(double));
  RU = (double *)malloc((1 + E->lmesh.nsf) * sizeof(double));
  SU = (double *)malloc((1 + E->lmesh.nsf) * sizeof(double));
  inp = (double *)malloc((6) * sizeof(double));
  outp = (double *)malloc((6) * sizeof(double));

  for (i = 1; i <= nno; i++)
  {
    mass[i] = 0.0;
    flux[i] = 0.0;
  }

  for (e = 1; e <= E->lmesh.nel; e++)
  {

    for (j = 1; j <= ends; j++)
      VZ[j] = E->V[2][E->ien[e].node[j]];

    for (i = 1; i <= ppts; i++)
    {
      u[i] = 0.0;
      T[i] = 0.0;
      dTdz[i] = 0.0;
      for (j = 1; j <= ends; j++)
      {
        u[i] += VZ[j] * E->N.ppt[GNPINDEX(j, i)];
        T[i] += E->T[E->ien[e].node[j]] * E->N.ppt[GNPINDEX(j, i)];
        dTdz[i] += -E->T[E->ien[e].node[j]] * E->gNX[e].ppt[GNPXINDEX(1, j, i)];
      }
    }

    uT = 0.0;
    area = 0.0;
    for (i = 1; i <= ppts; i++)
    {
      uT += u[i] * T[i] * E->gDA[e].ppt[i] + dTdz[i] * E->gDA[e].ppt[i];
      area += E->gDA[e].ppt[i];
    }

    uT /= area;
    for (j = 1; j <= ends; j++)
    {
      flux[E->ien[e].node[j]] += uT * E->gDA[e].ppt[1];
      mass[E->ien[e].node[j]] += E->gDA[e].ppt[1];
    }
  } /* end of e */

  for (i = 1; i <= E->lmesh.nsf; i++)
  {
    RU[i] = flux[E->surf_node[i]];
    SU[i] = mass[E->surf_node[i]];
    flux[E->surf_node[i]] = RU[i];
    mass[E->surf_node[i]] = SU[i];
    RU[i] = flux[E->surf_node[i] + 1];
    SU[i] = mass[E->surf_node[i] + 1];
    flux[E->surf_node[i] + 1] = RU[i];
    mass[E->surf_node[i] + 1] = SU[i];
  }
  for (i = 1; i <= E->lmesh.nsf; i++)
    E->slice.shflux[i] = -(2 * flux[E->surf_node[i]] / mass[E->surf_node[i]] - flux[E->surf_node[i] + 1] / mass[E->surf_node[i] + 1]);

  for (i = 1; i <= E->lmesh.nsf; i++)
  {
    RU[i] = flux[E->surf_node[i] + E->lmesh.noz - 1];
    SU[i] = mass[E->surf_node[i] + E->lmesh.noz - 1];
    flux[E->surf_node[i] + E->lmesh.noz - 1] = RU[i];
    mass[E->surf_node[i] + E->lmesh.noz - 1] = SU[i];
    RU[i] = flux[E->surf_node[i] + E->lmesh.noz - 2];
    SU[i] = mass[E->surf_node[i] + E->lmesh.noz - 2];
    flux[E->surf_node[i] + E->lmesh.noz - 2] = RU[i];
    mass[E->surf_node[i] + E->lmesh.noz - 2] = SU[i];
  }
  for (i = 1; i <= E->lmesh.nsf; i++)
    E->slice.bhflux[i] = -(2 * flux[E->surf_node[i] + E->lmesh.noz - 1] /
                               mass[E->surf_node[i] + E->lmesh.noz - 1] -
                           flux[E->surf_node[i] + E->lmesh.noz - 2] /
                               mass[E->surf_node[i] + E->lmesh.noz - 2]);

  areat = areab = hft = hfb = 0.0;

  for (i = 1; i <= E->lmesh.snel; i++)
  {
    tempb = tempt = 0.0;
    for (j = 1; j <= enodes[dims - 1]; j++)
    {
      tempb += E->slice.bhflux[E->sien[i].node[j]];
      tempt += E->slice.shflux[E->sien[i].node[j]];
    }
    e = i * E->lmesh.elz;
    hfb += tempb * E->eco[e].area;
    areab += E->eco[e].area;
    e = (i - 1) * E->lmesh.elz + 1;
    hft += tempt * E->eco[e].area;
    areat += E->eco[e].area;
  }

  inp[0] = hfb;
  inp[1] = areab;
  inp[2] = hft;
  inp[3] = areat;

  return_horiz_sum(E, inp, outp, 4);

  E->slice.Nub = outp[0] / (outp[1] * enodes[dims - 1]);
  E->slice.Nut = outp[2] / (outp[3] * enodes[dims - 1]);

  free((void *)flux);
  free((void *)mass);
  free((void *)RU);
  free((void *)SU);

  return;
}
