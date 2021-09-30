#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>

/* ========================================== */

void velocity_boundary_conditions(E) struct All_variables *E;
{
  void velocity_refl_vert_bc();
  void velocity_imp_vert_bc();
  void horizontal_bc();
  void velocity_apply_periodic_bcs();
  int lv;
  int node, d;
  int age_run, newnum, tempint;
  double *velosurf;
  char tempstring[255], input_s[255];
  FILE *fp_read;
  int keep_going, velonum;
  double loc_mid;
  for (lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
  {
    if (E->mesh.botvbc == 0)
    {

      horizontal_bc(E, E->VB, 1, 1, 0.0, VBX, 0, lv);
      horizontal_bc(E, E->VB, 1, 2, 0.0, VBZ, 1, lv);
      horizontal_bc(E, E->VB, 1, 1, E->control.VBXbotval, SBX, 1, lv);
      horizontal_bc(E, E->VB, 1, 2, 0.0, SBZ, 0, lv);
      if (E->mesh.nsd == 3)
      {
        horizontal_bc(E, E->VB, 1, 3, E->control.VBYbotval, SBY, 1, lv);
        horizontal_bc(E, E->VB, 1, 3, 0.0, VBY, 0, lv);
      }
    }
    if (E->mesh.topvbc == 0)
    {
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, 0.0, VBX, 0, lv);
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, 0.0, VBZ, 1, lv);
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, E->control.VBXtopval, SBX, 1, lv);
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, 0.0, SBZ, 0, lv);
      if (E->mesh.nsd == 3)
      {
        horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, E->control.VBYtopval, SBY, 1, lv);
        horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, 0.0, VBY, 0, lv);
      }
    }
  }

  velocity_refl_vert_bc(E); /* default */

  if (E->mesh.periodic_x || E->mesh.periodic_y)
    velocity_apply_periodic_bcs(E);

  for (lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
  {
    if (E->mesh.botvbc >= 1)
    {
      horizontal_bc(E, E->VB, 1, 1, E->control.VBXbotval, VBX, 1, lv);
      horizontal_bc(E, E->VB, 1, 2, 0.0, VBZ, 1, lv);
      horizontal_bc(E, E->VB, 1, 1, 0.0, SBX, 0, lv);
      horizontal_bc(E, E->VB, 1, 2, 0.0, SBZ, 0, lv);
      if (E->mesh.nsd == 3)
      {
        horizontal_bc(E, E->VB, 1, 3, E->control.VBYbotval, VBY, 1, lv);
        horizontal_bc(E, E->VB, 1, 3, 0.0, SBY, 0, lv);
      }
    }
    if (E->mesh.topvbc >= 1)
    {
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, E->control.VBXtopval, VBX, E->mesh.topvbc, lv);
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, 0.0, VBZ, 1, lv);
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, 0.0, SBX, 0, lv);
      horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, 0.0, SBZ, 0, lv);
      if (E->mesh.nsd == 3)
      {
        horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, E->control.VBYtopval, VBY, 1, lv);
        horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, 0.0, SBY, 0, lv);
      }
    }
  }

  if (E->control.trechmigrate)
  {
    if (E->parallel.me_loc[2] == E->parallel.nprocz - 1)
    {
      loc_mid = E->control.velo_surf_loc_mid;
      loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
      printf("loc_mid=%lf E->control.velo_surf_loc_mid_rate*E->monitor.elapsed_time*E->control.timescale=%lf\n", loc_mid, E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale);
      for (newnum = 1; newnum <= E->lmesh.nox; newnum++)
      {
        node = newnum * E->lmesh.noz;
        /*first left part */
        if (E->X[1][node] <= loc_mid)
        {
          E->VB[1][node] = E->control.velo_surf_mag_left * (tanh((E->X[1][node]) * E->control.velo_surf_width_left));
          /*overshoot*/
          if (E->X[1][node] >= loc_mid - E->control.velo_surf_loc_left_overshoot)
          {
            E->VB[1][node] += E->control.velo_surf_mag_left * sin((E->X[1][node] - loc_mid + E->control.velo_surf_loc_left_overshoot) / E->control.velo_surf_loc_left_overshoot * 3.141592653);
          }
        }
        /* then right part*/
        else
        {
          /*first right corner*/
          if (E->X[1][node] >= E->mesh.layer[1] - 0.2)
          {
            E->VB[1][node] = E->control.velo_surf_mag_right * (tanh((E->mesh.layer[1] - E->X[1][node]) * E->control.velo_surf_width_right));
          }
          /*then normal part*/
          else
          {
            E->VB[1][node] = (E->control.velo_surf_mag_right - E->control.velo_surf_mag_left) * tanh((E->X[1][node] - loc_mid) * E->control.velo_surf_width_mid) + E->control.velo_surf_mag_left;
          }
        } /*end of judge x location*/
      }   /*end of newnum loop*/
    }

  } /*end of trechmigrate*/

  if (E->control.imposevelo)
  {
    velosurf = (double *)malloc((E->lmesh.nox + 3) * sizeof(double));

    age_run = (int)(E->monitor.elapsed_time * E->control.timescale);

    sprintf(tempstring, "%s%d%s", E->control.velo_file_pre, E->control.age_total - age_run, E->control.velo_file_post);

    fp_read = fopen(tempstring, "r");
    if (fp_read == NULL)
      fprintf(stderr, "cannot open plate boundary file %s\n", tempstring);
    keep_going = 1;
    velonum = 0;

    while (keep_going)
    {
      if (fgets(input_s, 200, fp_read) == NULL)
      {
        keep_going = 0;
      }
      else
      {
        velonum++;
        sscanf(input_s, "%d %lf", &tempint, &velosurf[velonum]);
        E->VB[1][velonum * E->lmesh.noz] = velosurf[velonum];
        if (velonum <= 5)
          E->VB[1][velonum * E->lmesh.noz] = tanh((E->X[1][velonum * E->lmesh.noz])) * velosurf[velonum];
      }
    }
    for (newnum = velonum + 1; newnum <= E->lmesh.nox; newnum++)
    {
      E->VB[1][E->lmesh.noz * newnum] = velosurf[velonum];
      if (newnum >= E->lmesh.nox - 5)
        E->VB[1][newnum * E->lmesh.noz] = tanh((E->mesh.layer[1] - E->X[1][E->lmesh.noz * newnum])) * velosurf[velonum];
    }

  } /*end of reading velofile*/

  /*    for (node=1;node<=E->mesh.nnov;node++)
      fprintf(E->fp,"VB== %d %lf %lf \n",node,E->VB[1][node],E->VB[2][node]);
*/
  /*
if (E->control.verbose)     {
  if(E->mesh.nsd==3)
    for (node=1;node<=E->mesh.nnov;node++)
      fprintf(E->fp,"VB== %d %lf %lf %lf\n",node,E->VB[1][node],E->VB[2][node],E->VB[3][node]);
  else
    for (node=1;node<=E->mesh.nnov;node++)
      fprintf(E->fp,"VB== %d %lf %lf \n",node,E->VB[1][node],E->VB[2][node]);

  for(lv=E->mesh.levmax;lv>=E->mesh.levmin;lv--)    {
    fprintf(E->fp,"VBB level=%d %d\n",lv,E->mesh.NNO[lv]);
    for (node=1;node<=E->mesh.NNO[lv];node++)
      fprintf(E->fp,"VB== %d %u %u %u\n",node,E->NODE[lv][node]&VBX,E->NODE[lv][node]&VBZ,E->NODE[lv][node]&VBY);
    }
  }
*/
  if (E->control.imposevelo)
  {
    free(velosurf);
  }

  return;
}
/* ========================================== */

void composition_boundary_conditions(E) struct All_variables *E;
{
  void composition_refl_vert_bc();
  void compositions_conform_bcs();
  void horizontal_bc();

  if (E->mesh.botcbc == 1)
  {
    horizontal_bc(E, E->CB, 1, 2, E->control.CBCbotval, CBZ, 1, E->mesh.levmax);
    horizontal_bc(E, E->CB, 1, 2, E->control.CBCbotval, HBZ, 0, E->mesh.levmax);
  }
  else
  {
    horizontal_bc(E, E->CB, 1, 2, E->control.CBCbotval, CBZ, 0, E->mesh.levmax);
    horizontal_bc(E, E->CB, 1, 2, E->control.CBCbotval, HBZ, 1, E->mesh.levmax);
  }

  if (E->mesh.topcbc == 1)
  {
    horizontal_bc(E, E->CB, E->mesh.noz, 2, E->control.CBCtopval, CBZ, 1, E->mesh.levmax);
    horizontal_bc(E, E->CB, E->mesh.noz, 2, E->control.CBCtopval, HBZ, 0, E->mesh.levmax);
  }
  else
  {
    horizontal_bc(E, E->CB, E->mesh.noz, 2, E->control.CBCtopval, CBZ, 0, E->mesh.levmax);
    horizontal_bc(E, E->CB, E->mesh.noz, 2, E->control.CBCtopval, HBZ, 1, E->mesh.levmax);
  }

  composition_refl_vert_bc(E); /* default */

  /* zero flux is used for all the compositional field */

  /*
  compositions_conform_bcs(E);  
*/

  return;
}

/* ========================================== */

void temperature_boundary_conditions(E) struct All_variables *E;
{
  void temperature_refl_vert_bc();
  void temperatures_conform_bcs();
  void horizontal_bc();
  void temperature_apply_periodic_bcs();
  void temperature_imposed_vert_bcs();

  if (E->mesh.bottbc == 1)
  {
    horizontal_bc(E, E->TB, 1, 2, E->control.TBCbotval, TBZ, 1, E->mesh.levmax);
    horizontal_bc(E, E->TB, 1, 2, E->control.TBCbotval, FBZ, 0, E->mesh.levmax);
  }
  else
  {
    horizontal_bc(E, E->TB, 1, 2, E->control.TBCbotval, TBZ, 0, E->mesh.levmax);
    horizontal_bc(E, E->TB, 1, 2, E->control.TBCbotval, FBZ, 1, E->mesh.levmax);
  }

  if (E->mesh.toptbc == 1)
  {
    horizontal_bc(E, E->TB, E->mesh.noz, 2, E->control.TBCtopval, TBZ, 1, E->mesh.levmax);
    horizontal_bc(E, E->TB, E->mesh.noz, 2, E->control.TBCtopval, FBZ, 0, E->mesh.levmax);
  }
  else
  {
    horizontal_bc(E, E->TB, E->mesh.noz, 2, E->control.TBCtopval, TBZ, 0, E->mesh.levmax);
    horizontal_bc(E, E->TB, E->mesh.noz, 2, E->control.TBCtopval, FBZ, 1, E->mesh.levmax);
  }

  temperature_refl_vert_bc(E); /* default */

  if (E->mesh.periodic_x || E->mesh.periodic_y)
    temperature_apply_periodic_bcs(E);

  temperatures_conform_bcs(E);

  return;
}

/* ========================================== */

void velocity_refl_vert_bc(E) struct All_variables *E;
{
  int i, j, ii, jj;
  int node1, node2;
  int level, nox, noy, noz;
  const int dims = E->mesh.nsd;

  /* for two YOZ planes if 3-D, or two OZ side walls for 2-D */
  if (E->parallel.me_loc[1] == 0 || E->parallel.me_loc[1] == E->parallel.nprocx - 1)
    for (j = 1; j <= E->lmesh.noy; j++)
      for (i = 1; i <= E->lmesh.noz; i++)
      {
        node1 = i + (j - 1) * E->lmesh.noz * E->lmesh.nox;
        node2 = node1 + (E->lmesh.nox - 1) * E->lmesh.noz;

        E->VB[1][node1] = 0.0;
        E->VB[1][node2] = 0.0;
        if ((i != 1) && (i != E->lmesh.noz))
        {
          E->VB[2][node1] = 0.0;
          E->VB[2][node2] = 0.0;
        }
      } /* end loop for i and j */

  /* for two XOZ planes if 3-D */

  if (E->mesh.nsd == 3)
  {
    if (E->parallel.me_loc[2] == 0 || E->parallel.me_loc[2] == E->parallel.nprocy - 1)
      for (j = 1; j <= E->lmesh.nox; j++)
        for (i = 1; i <= E->lmesh.noz; i++)
        {
          node1 = i + (j - 1) * E->lmesh.noz;
          node2 = node1 + (E->lmesh.noy - 1) * E->lmesh.noz * E->lmesh.nox;
          if (E->mesh.nsd == 3)
          {
            E->VB[3][node2] = 0.0;
            E->VB[3][node1] = 0.0;
          }

          if ((i != 1) && (i != E->lmesh.noz))
          {
            E->VB[2][node2] = 0.0;
            E->VB[2][node1] = 0.0;
          }

        } /* end of loop i & j */

  } /* end of if */

  /* all vbc's apply at all levels  */
  for (level = E->mesh.levmax; level >= E->mesh.levmin; level--)
  {
    noz = E->lmesh.NOZ[level];
    noy = E->lmesh.NOY[level];
    nox = E->lmesh.NOX[level];
    if (E->parallel.me_loc[1] == 0 || E->parallel.me_loc[1] == E->parallel.nprocx - 1)
      for (j = 1; j <= noy; j++)
        for (i = 1; i <= noz; i++)
        {
          node1 = i + (j - 1) * noz * nox;
          node2 = node1 + (nox - 1) * noz;
          E->NODE[level][node1] = E->NODE[level][node1] | VBX;
          E->NODE[level][node1] = E->NODE[level][node1] & (~SBX);
          if ((i != 1) && (i != noz))
          {
            E->NODE[level][node1] = E->NODE[level][node1] & (~VBY);
            E->NODE[level][node1] = E->NODE[level][node1] | SBY;
            E->NODE[level][node1] = E->NODE[level][node1] & (~VBZ);
            E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
          }
          E->NODE[level][node2] = E->NODE[level][node2] | VBX;
          E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
          if ((i != 1) && (i != noz))
          {
            E->NODE[level][node2] = E->NODE[level][node2] & (~VBY);
            E->NODE[level][node2] = E->NODE[level][node2] | SBY;
            E->NODE[level][node2] = E->NODE[level][node2] & (~VBZ);
            E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
          }
        } /* end for loop i & j */

    if (E->mesh.nsd == 3)
    {
      if (E->parallel.me_loc[2] == 0 || E->parallel.me_loc[2] == E->parallel.nprocy - 1)
        for (j = 1; j <= nox; j++)
          for (i = 1; i <= noz; i++)
          {
            node1 = i + (j - 1) * noz;
            node2 = node1 + (noy - 1) * noz * nox;
            ii = i + E->lmesh.NZS[level] - 1;
            jj = j + E->lmesh.NXS[level] - 1;
            if (E->parallel.me_loc[2] == 0)
            {
              E->NODE[level][node1] = E->NODE[level][node1] | VBY;
              E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
              if ((ii != 1) && (ii != E->mesh.NOZ[level]))
              {
                E->NODE[level][node1] = E->NODE[level][node1] & (~VBZ);
                E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
              }
              if ((jj != 1) && (jj != E->mesh.NOX[level]) && (ii != 1) && (ii != E->mesh.NOZ[level]))
              {
                E->NODE[level][node1] = E->NODE[level][node1] & (~VBX);
                E->NODE[level][node1] = E->NODE[level][node1] | SBX;
              }
            }
            if (E->parallel.me_loc[2] == E->parallel.nprocy - 1)
            {
              E->NODE[level][node2] = E->NODE[level][node2] | VBY;
              E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
              if ((ii != 1) && (ii != E->mesh.NOZ[level]))
              {
                E->NODE[level][node2] = E->NODE[level][node2] & (~VBZ);
                E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
              }
              if ((jj != 1) && (jj != E->mesh.NOX[level]) && (ii != 1) && (ii != E->mesh.NOZ[level]))
              {
                E->NODE[level][node2] = E->NODE[level][node2] & (~VBX);
                E->NODE[level][node2] = E->NODE[level][node2] | SBX;
              }
            }
          } /* end for loop i & j  */

    } /* end for if dims=3 */
  }   /* end for loop level */

  return;
}

/* =========================================== */
void temperature_refl_vert_bc(E) struct All_variables *E;
{
  int i, j;
  int node1, node2;
  const int dims = E->mesh.nsd;

  /* Temps and bc-values  at top level only */
  if (E->parallel.me_loc[1] == 0 || E->parallel.me_loc[1] == E->parallel.nprocx - 1)
    for (j = 1; j <= E->lmesh.noy; j++)
      for (i = 1; i <= E->lmesh.noz; i++)
      {
        node1 = i + (j - 1) * E->mesh.noz * E->mesh.nox;
        node2 = node1 + (E->mesh.nox - 1) * E->mesh.noz;
        E->node[node1] = E->node[node1] & (~TBX);
        E->node[node1] = E->node[node1] | FBX;
        E->TB[1][node1] = 0.0;
        E->node[node2] = E->node[node2] & (~TBX);
        E->node[node2] = E->node[node2] | FBX;
        E->TB[1][node2] = 0.0;
      } /* end for loop i & j */

  if (E->mesh.nsd == 3)
  {
    if (E->parallel.me_loc[2] == 0 || E->parallel.me_loc[2] == E->parallel.nprocy - 1)
      for (j = 1; j <= E->lmesh.nox; j++)
        for (i = 1; i <= E->lmesh.noz; i++)
        {
          node1 = i + (j - 1) * E->lmesh.noz;
          E->node[node1] = E->node[node1] & (~TBY);
          E->node[node1] = E->node[node1] | FBY;
          E->TB[3][node1] = 0.0;
        }
    for (j = 1; j <= E->lmesh.nox; j++)
      for (i = 1; i <= E->lmesh.noz; i++)
      {
        node2 = i + (j - 1) * E->lmesh.noz + (E->lmesh.noy - 1) * E->lmesh.noz * E->lmesh.nox;
        E->node[node2] = E->node[node2] & (~TBY);
        E->node[node2] = E->node[node2] | FBY;
        E->TB[3][node2] = 0.0;
      } /* end loop for i and j */
  }     /* end for if ==3    */

  return;
}

/* =========================================== */
void composition_refl_vert_bc(E) struct All_variables *E;
{
  int i, j;
  int node1, node2;
  const int dims = E->mesh.nsd;

  /* Temps and bc-values  at top level only */

  for (j = 1; j <= E->mesh.noy; j++)
    for (i = 1; i <= E->mesh.noz; i++)
    {
      node1 = i + (j - 1) * E->lmesh.noz * E->lmesh.nox;
      node2 = node1 + (E->lmesh.nox - 1) * E->lmesh.noz;
      E->node[node1] = E->node[node1] & (~CBX);
      E->node[node1] = E->node[node1] | HBX;
      E->CB[1][node1] = 0.0;
      E->node[node2] = E->node[node2] & (~CBX);
      E->node[node2] = E->node[node2] | HBX;
      E->CB[1][node2] = 0.0;
    } /* end for loop i & j */

  if (E->mesh.nsd == 3)
  {
    for (j = 1; j <= E->lmesh.nox; j++)
      for (i = 1; i <= E->lmesh.noz; i++)
      {
        node1 = i + (j - 1) * E->mesh.noz;
        E->node[node1] = E->node[node1] & (~CBY);
        E->node[node1] = E->node[node1] | HBY;
        E->CB[3][node1] = 0.0;
      }
    for (j = 1; j <= E->lmesh.nox; j++)
      for (i = 1; i <= E->lmesh.noz; i++)
      {
        node2 = i + (j - 1) * E->lmesh.noz + (E->lmesh.noy - 1) * E->lmesh.noz * E->lmesh.nox;
        E->node[node2] = E->node[node2] & (~CBY);
        E->node[node2] = E->node[node2] | HBY;
        E->CB[3][node2] = 0.0;
      } /* end loop for i and j */
  }     /* end for if ==3    */

  return;
}

/*  =========================================================  */

void horizontal_bc(E, BC, ROW, dirn, value, mask, onoff, level) struct All_variables *E;
double *BC[];
int ROW;
int dirn;
double value;
unsigned int mask;
int onoff;
int level;

{
  int i, j, node, rowl;
  const int dims = E->mesh.nsd;

  /* safety feature */
  if (dirn > E->mesh.nsd)
    return;

  if (ROW == 1)
    rowl = 1;
  else
    rowl = E->lmesh.NOZ[level];

  if (ROW == 1 && E->parallel.me_loc[3] == 0 || ROW == E->mesh.NOZ[level] && E->parallel.me_loc[3] == E->parallel.nprocz - 1)
  {

    /* turn bc marker to zero */
    if (onoff == 0)
    {
      for (j = 1; j <= E->lmesh.NOY[level]; j++)
        for (i = 1; i <= E->lmesh.NOX[level]; i++)
        {
          node = rowl + (i - 1) * E->lmesh.NOZ[level] + (j - 1) * E->lmesh.NOX[level] * E->lmesh.NOZ[level];
          E->NODE[level][node] = E->NODE[level][node] & (~mask);
        } /* end for loop i & j */
    }

    /* turn bc marker to one */
    else
    {
      for (j = 1; j <= E->lmesh.NOY[level]; j++)
        for (i = 1; i <= E->lmesh.NOX[level]; i++)
        {
          node = rowl + (i - 1) * E->lmesh.NOZ[level] + (j - 1) * E->lmesh.NOX[level] * E->lmesh.NOZ[level];
          E->NODE[level][node] = E->NODE[level][node] | (mask);
          if (level == E->mesh.levmax)
          { /* NB */
            BC[dirn][node] = value;
          }
          /**/
        } /* end for loop i & j */
    }

  } /* end for if ROW */

  return;
}

void velocity_apply_periodic_bcs(E) struct All_variables *E;
{
  int n1, n2, level;
  int i, j, ii, jj;
  const int dims = E->mesh.nsd;

  fprintf(E->fp, "Periodic boundary conditions\n");

  return;
}

void temperature_apply_periodic_bcs(E) struct All_variables *E;
{
  int n1, n2, e1, level;
  int i, j, ii, jj;
  const int dims = E->mesh.nsd;

  fprintf(E->fp, "Periodic temperature boundary conditions\n");

  return;
}

void strip_bcs_from_residual(E, Res, level) struct All_variables *E;
double *Res;
int level;
{
  int i;
  const int dims = E->mesh.nsd;

  for (i = 1; i <= E->mesh.NNO[level]; i++)
  {
    if (E->NODE[level][i] & OFFSIDE)
      continue;
    if ((E->NODE[level][i] & VBX) != 0)
      Res[E->ID[level][i].doff[1]] = 0.0;
    if ((E->NODE[level][i] & VBZ) != 0)
      Res[E->ID[level][i].doff[2]] = 0.0;
    if (3 == dims && ((E->NODE[level][i] & VBY) != 0))
      Res[E->ID[level][i].doff[3]] = 0.0;
  }

  return;
}

void temperatures_conform_bcs(E) struct All_variables *E;
{
  int node;
  unsigned int type;

  for (node = 1; node <= E->lmesh.nno; node++)
  {
    if (E->node[node] & OFFSIDE)
      continue;

    type = (E->node[node] & (TBX | TBZ | TBY));

    switch (type)
    {
    case 0: /* no match, next node */
      break;
    case TBX:
      E->T[node] = E->TB[1][node];
      break;
    case TBZ:
      E->T[node] = E->TB[2][node];
      break;
    case TBY:
      E->T[node] = E->TB[3][node];
      break;
    case (TBX | TBZ): /* clashes ! */
      E->T[node] = 0.5 * (E->TB[1][node] + E->TB[2][node]);
      break;
    case (TBX | TBY): /* clashes ! */
      E->T[node] = 0.5 * (E->TB[1][node] + E->TB[3][node]);
      break;
    case (TBZ | TBY): /* clashes ! */
      E->T[node] = 0.5 * (E->TB[2][node] + E->TB[3][node]);
      break;
    case (TBZ | TBY | TBX): /* clashes ! */
      E->T[node] = 0.3333333 * (E->TB[1][node] + E->TB[2][node] + E->TB[3][node]);
      break;
    }

    /* next node */
  }

  return;
}

void velocities_conform_bcs(E, U) struct All_variables *E;
double *U;
{
  int node, d;

  const unsigned int typex = VBX;
  const unsigned int typez = VBZ;
  const unsigned int typey = VBY;
  const int addi_dof = additional_dof[E->mesh.nsd];

  const int dofs = E->mesh.dof;
  const int nno = E->lmesh.nno;
  int age_run, newnum, tempint;
  double *velosurf;
  char tempstring[255], input_s[255];
  FILE *fp_read;
  int keep_going, velonum;
  double loc_mid;

  if (E->control.imposevelo)
  {
    velosurf = (double *)malloc((E->mesh.nox + 3) * sizeof(double));

    age_run = (int)(E->monitor.elapsed_time * E->control.timescale);

    sprintf(tempstring, "%s%d%s", E->control.velo_file_pre, E->control.age_total - age_run, E->control.velo_file_post);

    fp_read = fopen(tempstring, "r");
    if (fp_read == NULL)
      fprintf(stderr, "cannot open plate boundary file %s\n", tempstring);
    keep_going = 1;
    velonum = 0;
    while (keep_going)
    {
      if (fgets(input_s, 200, fp_read) == NULL)
      {
        keep_going = 0;
      }
      else
      {
        velonum++;
        sscanf(input_s, "%d %lf", &tempint, &velosurf[velonum]);
        E->VB[1][velonum * E->lmesh.noz] = velosurf[velonum];
        if (velonum <= 5)
          E->VB[1][velonum * E->lmesh.noz] = tanh((E->X[1][velonum * E->lmesh.noz])) * velosurf[velonum];
      }
    }
    for (newnum = velonum + 1; newnum <= E->lmesh.nox; newnum++)
    {
      E->VB[1][E->lmesh.noz * newnum] = velosurf[velonum];
      if (newnum >= E->lmesh.nox - 5)
        E->VB[1][newnum * E->lmesh.noz] = tanh((E->mesh.layer[1] - E->X[1][E->lmesh.noz * newnum])) * velosurf[velonum];
    }

  } /*end of reading velofile*/

  if (E->control.trechmigrate)
  {
    if (E->parallel.me_loc[2] == E->parallel.nprocz - 1)
    {
      loc_mid = E->control.velo_surf_loc_mid;
      loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
      printf("loc_mid=%lf E->control.velo_surf_loc_mid_rate*E->monitor.elapsed_time*E->control.timescale=%lf\n", loc_mid, E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale);
      for (newnum = 1; newnum <= E->lmesh.nox; newnum++)
      {
        node = newnum * E->lmesh.noz;
        /*first left part */
        if (E->X[1][node] <= loc_mid)
        {
          E->VB[1][node] = E->control.velo_surf_mag_left * (tanh((E->X[1][node]) * E->control.velo_surf_width_left));
          /*overshoot*/
          if (E->X[1][node] >= loc_mid - E->control.velo_surf_loc_left_overshoot)
          {
            E->VB[1][node] += E->control.velo_surf_mag_left * sin((E->X[1][node] - loc_mid + E->control.velo_surf_loc_left_overshoot) / E->control.velo_surf_loc_left_overshoot * 3.141592653);
          }
        }
        /* then right part*/
        else
        {
          /*first right corner*/
          if (E->X[1][node] >= E->mesh.layer[1] - 0.2)
          {
            E->VB[1][node] = E->control.velo_surf_mag_right * (tanh((E->mesh.layer[1] - E->X[1][node]) * E->control.velo_surf_width_right));
          }
          /*then normal part*/
          else
          {
            E->VB[1][node] = (E->control.velo_surf_mag_right - E->control.velo_surf_mag_left) * tanh((E->X[1][node] - loc_mid) * E->control.velo_surf_width_mid) + E->control.velo_surf_mag_left;
          }

        } /*end of judge x location*/
      }   /*end of newnum loop*/
    }

  } /*end of trechmigrate*/

  for (node = 1; node <= nno; node++)
  {
    if (E->node[node] & OFFSIDE)
      continue;

    if (E->node[node] & typex)
      U[E->id[node].doff[1]] = E->VB[1][node];
    if (E->node[node] & typez)
      U[E->id[node].doff[2]] = E->VB[2][node];
    if (3 == dofs && E->node[node] & typey)
      U[E->id[node].doff[3]] = E->VB[3][node];
  }

  if (E->control.imposevelo)
  {
    free(velosurf);
  }

  return;
}
