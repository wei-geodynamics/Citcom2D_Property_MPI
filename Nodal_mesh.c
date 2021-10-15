/* Functions relating to the building and use of mesh locations ... */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

/* =================================================
   Standard node positions including mesh refinement 

   =================================================  */

void node_locations(E) struct All_variables *E;
{
  int lev, nodel, i, j, k, ii, d, node, ix, iy, iz;
  double *XG[4], dx[4], dx1, dx2;
  double *XX[4], dxx[40], x00;
  double PI = 3.141592635;
  int n00, nox, noz, noy, fn;

  const int dims = E->mesh.nsd;

  void inject_scalar();
  void inject_node_fvector();

  /* Option to deform mesh to test non-orthogonality, see
   * tracer-impl.schrift, p23), see also around line 166: 
   * 0=no deform, 1=deform */
  E->control.testweirdmesh = 0;

  input_int("z_grid_layers", &(E->segment.zlayers), "1");
  input_double_vector("zz", E->segment.zlayers, (E->segment.zzlayer));
  input_int_vector("nz", E->segment.zlayers, (E->segment.nzlayer));

  input_int("x_grid_layers", &(E->segment.xlayers), "1");
  input_double_vector("xx", E->segment.xlayers, (E->segment.xxlayer));
  input_int_vector("nx", E->segment.xlayers, (E->segment.nxlayer));

  if (dims == 3)
  {
    input_int("y_grid_layers", &(E->segment.ylayers), "1");
    input_double_vector("yy", E->segment.ylayers, (E->segment.yylayer));
    input_int_vector("ny", E->segment.ylayers, (E->segment.nylayer));
  }

  nox = E->mesh.nox;
  noz = E->mesh.noz;
  noy = E->mesh.noy;

  for (d = 1; d <= E->mesh.nsd; d++)
  {
    XX[d] = (double *)malloc((2 + E->mesh.nnx[d]) * sizeof(double));
    XG[d] = (double *)malloc((2 + 1) * sizeof(double));
    E->mesh.globalXX[d] = (double *)malloc((2 + E->mesh.nnx[d]) * sizeof(double)); /*VM 20130429*/
  }

  for (d = 1; d <= E->mesh.nsd; d++)
  { /* for even space */
    dx[d] = E->mesh.layer[d] / (E->mesh.nnx[d] - 1);
    XX[d][1] = 0.0;
    XX[d][E->mesh.nnx[d]] = E->mesh.layer[d];
    for (i = 2; i < E->mesh.nnx[d]; i++)
      XX[d][i] = XX[d][i - 1] + dx[d];
  }

  for (j = 1; j < E->segment.xlayers; j++)
    dxx[j] = (E->segment.xxlayer[j] - E->segment.xxlayer[j - 1]) / (E->segment.nxlayer[j] - E->segment.nxlayer[j - 1]);
  j = 1;
  for (i = 2; i < E->mesh.nnx[1]; i++)
  {
    if (i <= E->segment.nxlayer[j])
      XX[1][i] = XX[1][i - 1] + dxx[j];
    if (i == E->segment.nxlayer[j])
      j++;
  }

  for (j = 1; j < E->segment.zlayers; j++)
    dxx[j] = (E->segment.zzlayer[j] - E->segment.zzlayer[j - 1]) / (E->segment.nzlayer[j] - E->segment.nzlayer[j - 1]);
  j = 1;
  for (i = 2; i < E->mesh.nnx[2]; i++)
  {
    if (i <= E->segment.nzlayer[j])
      XX[2][i] = XX[2][i - 1] + dxx[j];
    if (i == E->segment.nzlayer[j])
      j++;
  }

  /* VM - Global array with nodes coordinates */
  for (i = 1; i <= E->mesh.nnx[1]; i++)
    E->mesh.globalXX[1][i] = (double)XX[1][i];
  for (i = 1; i <= E->mesh.nnx[2]; i++)
    E->mesh.globalXX[2][i] = (double)XX[2][i];
  if (E->mesh.nsd == 3)
  {
    for (i = 1; i <= E->mesh.nnx[3]; i++)
      E->mesh.globalXX[3][i] = (double)XX[3][i];
  }

  lev = E->mesh.levmax;
  nox = E->lmesh.NOX[lev];
  noy = E->lmesh.NOY[lev];
  noz = E->lmesh.NOZ[lev];

  /* Fill E->XL array, a per-processor version of XX */
  E->XL[1] = (double *)malloc((nox + 1) * sizeof(double));
  E->XL[2] = (double *)malloc((noz + 1) * sizeof(double));
  if (E->mesh.nsd == 3)
    E->XL[3] = (double *)malloc((noy + 1) * sizeof(double));
  for (ix = 1; ix <= nox; ix++)
    E->XL[1][ix] = XX[1][ix - 1 + E->lmesh.NXS[lev]];
  for (iz = 1; iz <= noz; iz++)
    E->XL[2][iz] = XX[2][iz - 1 + E->lmesh.NZS[lev]];
  if (E->mesh.nsd == 3)
    for (iy = 1; iy <= noy; iy++)
      E->XL[3][iy] = XX[3][iy - 1 + E->lmesh.NYS[lev]];
  if (E->mesh.nsd == 3)
    for (k = 1; k <= noy; k++)
      for (i = 1; i <= nox; i++)
        for (j = 1; j <= noz; j++)
        {
          nodel = j + (i - 1) * noz + (k - 1) * noz * nox;
          E->XX[lev][1][nodel] = XX[1][i - 1 + E->lmesh.NXS[lev]];
          E->XX[lev][2][nodel] = XX[2][j - 1 + E->lmesh.NZS[lev]];
          E->XX[lev][3][nodel] = XX[3][k - 1 + E->lmesh.NYS[lev]];
        }
  else if (E->mesh.nsd == 2)
  {
    for (i = 1; i <= nox; i++)
      for (j = 1; j <= noz; j++)
      {
        nodel = j + (i - 1) * noz;
        E->XX[lev][1][nodel] = XX[1][i - 1 + E->lmesh.NXS[lev]];
        E->XX[lev][2][nodel] = XX[2][j - 1 + E->lmesh.NZS[lev]];
      }
    /* Testing:
    */
    for (ix = 1; ix <= nox; ix++)
      fprintf(E->fp, "XL[1][%d]=%.9e, E->XX=%.9e\n", ix,
              E->XL[1][ix], E->XX[lev][1][1 + (ix - 1) * noz]);
    for (iz = 1; iz <= noz; iz++)
      fprintf(E->fp, "XL[2][%d]=%.9e\n", iz, E->XL[2][iz]);
    if (E->mesh.nsd == 3)
      for (iy = 1; iy <= noy; iy++)
        fprintf(E->fp, "XL[3][%d]=%.9e\n", ix, E->XL[3][iy]);

    for (j = 1; j <= E->lmesh.noz; j++)
      E->XP[2][j] = E->XX[E->mesh.levmax][2][j];
    for (i = 1; i <= E->lmesh.nox; i++)
      E->XP[1][i] = E->XX[E->mesh.levmax][1][1 + (i - 1) * E->lmesh.noz];
    E->XG1[1] = 1.0e-6;
    E->XG1[2] = 1.0e-6;

    E->XG2[1] = E->mesh.layer[1] - 1.0e-6;
    E->XG2[2] = E->mesh.layer[2] - 1.0e-6;

    if (E->control.testweirdmesh == 1)
    {
      for (i = 1; i <= nox; i++)
        for (j = 1; j <= noz; j++)
        {
          nodel = j + (i - 1) * noz;
          E->XX[lev][1][nodel] +=
              0.25 * sin(E->XX[lev][1][nodel] * PI) *
              (E->XX[lev][2][nodel] - 0.5);
          fprintf(stderr, "Adjustment x-coor = %g, z=%g\n",
                  0.25 * sin(E->XX[lev][1][nodel] * PI) *
                      (E->XX[lev][2][nodel] - 0.5),
                  E->XX[lev][2][nodel]);
        }
    }
  }
  for (lev = E->mesh.levmax; lev > E->mesh.levmin; lev--)
  {
    inject_node_fvector(E, lev, E->XX[lev], E->XX[lev - 1]);
  }

  if (E->control.verbose)
  {
    for (lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
    {
      fprintf(E->fp, "output_coordinates \n");
      if (dims == 2)
        for (i = 1; i <= E->lmesh.NNO[lev]; i++)
          fprintf(E->fp, "%d %g %g\n", i, E->XX[lev][1][i], E->XX[lev][2][i]);
      else if (dims == 3)
        for (i = 1; i <= E->lmesh.NNO[lev]; i++)
          fprintf(E->fp, "%d %g %g %g\n", i, E->XX[lev][1][i], E->XX[lev][2][i], E->XX[lev][3][i]);
    }
  }

  for (d = 1; d <= E->mesh.nsd; d++)
  {
    free((void *)XX[d]);
    free((void *)XG[d]);
  }

  return;
}

void dlogical_mesh_to_real(E, data, level) struct All_variables *E;
double *data;
int level;

{
  int i, j, n1, n2;

  if (E->mesh.periodic_x)
    for (i = 1; i <= E->mesh.NOZ[level]; i++)
      for (j = 1; j <= E->mesh.NOY[level]; j++)
      {
        n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
        n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

        data[n2] = data[n1];
      }

  if (E->mesh.periodic_y)
    for (i = 1; i <= E->mesh.NOZ[level]; i++)
      for (j = 1; j <= E->mesh.NOX[level]; j++)
      {
        n1 = i + (j - 1) * E->mesh.NOZ[level];
        n2 = n1 + (E->mesh.NOY[level] - 1) * E->mesh.NOZ[level] * E->mesh.NOX[level];

        data[n2] = data[n1];
      }

  if (E->mesh.periodic_y && E->mesh.periodic_x) /* then need to do the 1st one again */
    for (i = 1; i <= E->mesh.NOZ[level]; i++)
      for (j = 1; j <= E->mesh.NOY[level]; j++)
      {
        n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
        n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

        data[n2] = data[n1];
      }

  return;
}

void flogical_mesh_to_real(E, data, level) struct All_variables *E;
double *data;
int level;

{
  int i, j, n1, n2;

  if (E->mesh.periodic_x)
    for (i = 1; i <= E->mesh.NOZ[level]; i++)
      for (j = 1; j <= E->mesh.NOY[level]; j++)
      {
        n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
        n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

        data[n2] = data[n1];
      }

  if (E->mesh.periodic_y)
    for (i = 1; i <= E->mesh.NOZ[level]; i++)
      for (j = 1; j <= E->mesh.NOX[level]; j++)
      {
        n1 = i + (j - 1) * E->mesh.NOZ[level];
        n2 = n1 + (E->mesh.NOY[level] - 1) * E->mesh.NOZ[level] * E->mesh.NOX[level];

        data[n2] = data[n1];
      }

  if (E->mesh.periodic_y && E->mesh.periodic_x) /* then need to do the 1st one again */
    for (i = 1; i <= E->mesh.NOZ[level]; i++)
      for (j = 1; j <= E->mesh.NOY[level]; j++)
      {
        n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
        n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

        data[n2] = data[n1];
      }

  return;
}

void p_to_nodes(E, P, PN, lev)
    /* NB The use of E->TW is incorrect in combi with mesh refinement */
    struct All_variables *E;
double *P;
double *PN;
int lev;

{
  int e, element, node, j;
  void e_exchange_node_fc();
  void exchange_node_f20();

  for (node = 1; node <= E->lmesh.NNO[lev]; node++)
    PN[node] = 0.0;

  for (element = 1; element <= E->lmesh.NEL[lev]; element++)
  {

    for (j = 1; j <= enodes[E->mesh.nsd]; j++)
    {
      node = E->IEN[lev][element].node[j];
      PN[node] += P[element] * E->TW[lev][node];
    }
  }

  exchange_node_f20(E, PN, lev);

  return;
}

void edot_to_nodes(E, P, PN, lev) /* As p_to_nodes, but input in single 
                                  precision double */
                                  /* NB The use of E->TW is incorrect in combi with mesh refinement */
    struct All_variables *E;
double *P;
double *PN;
int lev;

{
  int e, element, node, j;
  void e_exchange_node_fc();
  void exchange_node_f20();

  for (node = 1; node <= E->lmesh.NNO[lev]; node++)
    PN[node] = 0.0;

  for (element = 1; element <= E->lmesh.NEL[lev]; element++)
  {

    for (j = 1; j <= enodes[E->mesh.nsd]; j++)
    {
      node = E->IEN[lev][element].node[j];
      PN[node] += P[element] * E->TW[lev][node];
    }
  }

  exchange_node_f20(E, PN, lev);

  return;
}

void p_to_centres(E, PN, P, lev) struct All_variables *E;
double *PN;
double *P;
int lev;

{
  int p, element, node, j;
  double weight;

  for (p = 1; p <= E->lmesh.NEL[lev]; p++)
    P[p] = 0.0;

  weight = 1.0 / ((double)enodes[E->mesh.nsd]);

  for (p = 1; p <= E->lmesh.NEL[lev]; p++)
    for (j = 1; j <= enodes[E->mesh.nsd]; j++)
      P[p] += PN[E->IEN[lev][p].node[j]] * weight;

  return;
}

void v_to_intpts(E, VN, VE, lev) struct All_variables *E;
double *VN, *VE;
int lev;
{

  int e, i, j, k;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
    for (i = 1; i <= vpts; i++)
    {
      VE[(e - 1) * vpts + i] = 0.0;
      for (j = 1; j <= ends; j++)
        VE[(e - 1) * vpts + i] += VN[E->IEN[lev][e].node[j]] * E->N.vpt[GNVINDEX(j, i)];
    }
}

void v_to_nodes(E, VE, VN, lev) struct All_variables *E;
double *VE, *VN;
int lev;
{
  int e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];
  for (i = 1; i <= E->lmesh.NNO[lev]; i++)
    VN[i] = 0.0;

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
    for (j = 1; j <= ends; j++)
    {
      n = E->IEN[lev][e].node[j];
      for (i = 1; i <= vpts; i++)
        VN[n] += E->N.vpt[GNVINDEX(j, i)] * E->TW[lev][n] * VE[(e - 1) * vpts + i];
    }
  flogical_mesh_to_real(E, VN, E->mesh.levmax);
  return;
}

void visc_to_intpts(E, VN, VE, lev) struct All_variables *E;
double *VN, *VE;
int lev;
{

  int e, i, j, k;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
    for (i = 1; i <= vpts; i++)
    {
      VE[(e - 1) * vpts + i] = 0.0;
      for (j = 1; j <= ends; j++)
        VE[(e - 1) * vpts + i] += log(VN[E->IEN[lev][e].node[j]]) * E->N.vpt[GNVINDEX(j, i)];
      VE[(e - 1) * vpts + i] = exp(VE[(e - 1) * vpts + i]);
    }
}

void visc_to_nodes(E, VE, VN, lev) struct All_variables *E;
double *VE, *VN;
int lev;
{
  int e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];
  double temp_visc;

  for (i = 1; i <= E->lmesh.NNO[lev]; i++)
    VN[i] = 0.0;

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
    for (j = 1; j <= ends; j++)
    {
      n = E->IEN[lev][e].node[j];
      temp_visc = 0.0;
      for (i = 1; i <= vpts; i++)
        temp_visc += E->TW[lev][n] * log(E->N.vpt[GNVINDEX(j, i)] * VE[(e - 1) * vpts + i]);
      VN[n] += exp(temp_visc);
    }
  return;
}

void visc_from_ele_to_gint(E, VN, VE, lev) struct All_variables *E;
double *VE, *VN;
int lev;
{

  int m, e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
    for (i = 1; i <= vpts; i++)
    {
      VE[(e - 1) * vpts + i] = VN[e];
    }

  return;
}

void visc_from_gint_to_ele(E, VE, VN, lev) struct All_variables *E;
double *VE, *VN;
int lev;
{
  int m, e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];
  double temp_visc;

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
  {
    temp_visc = 0.0;
    for (i = 1; i <= vpts; i++)
      temp_visc += VE[(e - 1) * vpts + i];
    temp_visc = temp_visc / vpts;

    VN[e] = temp_visc;
  }

  return;
}

void visc_from_gint_to_elel(E, VE, VN, lev) struct All_variables *E;
double *VE, *VN;
int lev;
{
  int m, e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];
  double temp_visc;

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
  {
    temp_visc = 1.0;
    for (i = 1; i <= vpts; i++)
      temp_visc *= VE[(e - 1) * vpts + i];
    temp_visc = pow(temp_visc, 1. / vpts);

    VN[e] = temp_visc;
  }

  return;
}

void visc_from_gint_to_nodes(E, VE, VN, lev) struct All_variables *E;
double *VE, *VN;
int lev;
{
  int m, e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];
  double temp_visc;
  void exchange_node_f20();

  for (i = 1; i <= E->lmesh.NNO[lev]; i++)
    VN[i] = 0.0;

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
  {
    temp_visc = 0.0;
    for (i = 1; i <= vpts; i++)
      temp_visc += VE[(e - 1) * vpts + i];
    temp_visc = temp_visc / vpts;

    for (j = 1; j <= ends; j++)
    {
      n = E->IEN[lev][e].node[j];
      VN[n] += E->TWW[lev][e].node[j] * temp_visc;
    }
  }

  exchange_node_f20(E, VN, lev);

  for (n = 1; n <= E->lmesh.NNO[lev]; n++)
    VN[n] *= E->MASS[lev][n];

  return;
}

void fld_from_elem_to_nodes(E, efld, nfld, lev) struct All_variables *E;
double *efld, *nfld;
int lev;
{
  int m, e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];
  void exchange_node_f20();

  for (i = 1; i <= E->lmesh.NNO[lev]; i++)
    nfld[i] = 0.0;

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
  {
    for (j = 1; j <= ends; j++)
    {
      n = E->IEN[lev][e].node[j];
      nfld[n] += E->TWW[lev][e].node[j] * efld[e];
    }
  }

  exchange_node_f20(E, nfld, lev);

  for (n = 1; n <= E->lmesh.NNO[lev]; n++)
    nfld[n] *= E->MASS[lev][n];

  return;
}

void visc_from_nodes_to_gint(E, VN, VE, lev) struct All_variables *E;
double *VE, *VN;
int lev;
{

  int m, e, i, j, k, n;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  const int ends = enodes[nsd];
  double temp_visc;
  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
    for (i = 1; i <= vpts; i++)
      VE[(e - 1) * vpts + i] = 0.0;

  for (e = 1; e <= E->lmesh.NEL[lev]; e++)
    for (i = 1; i <= vpts; i++)
    {
      temp_visc = 0.0;
      for (j = 1; j <= ends; j++)
        temp_visc += E->N.vpt[GNVINDEX(j, i)] * VN[E->IEN[lev][e].node[j]];
      VE[(e - 1) * vpts + i] = temp_visc;
    }

  return;
}
