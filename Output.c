/* Routine to process the output of the finite element cycles 
   and to turn them into a coherent suite  files  */

#include <fcntl.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h> /* for "system" command */
#ifndef __sunos__   /* string manipulations */
#include <strings.h>
#else
#include <string.h>
#endif
#include <mpi.h>
#include "element_definitions.h"
#include "global_defs.h"
void output_velo_related(E, file_number) struct All_variables *E;
int file_number;
{
  FILE *fp1, *fp0;
  double surf, botm;
  int el, els, i, j, k, ii, jj, kk, m, node, fd;
  int nox, noz, noy, nfx, nfz, nfy1, nfy2, size1, size2;
  char output_file[255];
  static double *SV, *EV;
  static int been_here = 0;
  double *eedot, *nedot;
  double *eedot11, *nedot11, *eedot22, *nedot22, *eedot12, *nedot12;
  double *edot, *edot11, *edot22, *edot12;
  int n, e;
  void strain_rate_2_inv_moreout();
  double temp, temp11, temp22, temp12;
  void get_surface_velo();
  void get_ele_visc();
  void coordinates_dx();
  void frames_for_dx();
  void return_horiz_ave();
  void output_phasefile();
  void visc_from_gint_to_nodes();
  const int nno = E->lmesh.nno;
  int num1, num2;
  if (been_here == 0)
    SV = (double *)malloc((nno + 1) * sizeof(double));

  if (been_here == 0 && E->control.restart == 0)
  {
    been_here++;

    sprintf(output_file, "%s/coord.%d", E->control.data_file, E->parallel.me);
    fp0 = fopen(output_file, "w");
    fprintf(fp0, "%6d %6d %.5e\n", E->mesh.nno, E->advection.timesteps, E->monitor.elapsed_time);
    for (i = 1; i <= E->lmesh.nno; i++)
      fprintf(fp0, "%6d %.3e %.3e %.5e %.5e %.5e %.5e %.4e\n", i, E->X[1][i], E->X[2][i], E->V[1][i], E->V[2][i], E->T[i], E->C[i], E->Vi[i]);
    fclose(fp0);
  }

  if ((E->advection.timesteps % (E->control.record_every)) == 0)
  {

    sprintf(output_file, "%s/temp_comp.%d.%d", E->control.data_file, E->parallel.me, file_number);
    fp0 = fopen(output_file, "w");

    fprintf(fp0, "%6d %6d %.5e\n", E->lmesh.nno, E->advection.timesteps, E->monitor.elapsed_time);
    if (E->control.composition)
    {
      if ((E->advection.timesteps % (1 * E->control.record_every)) == 0)
      {
        for (i = 1; i <= E->lmesh.nno; i++)
        {
          if (E->control.phasevisc_C)
          {
            if (E->control.phasevisc_d)
            {
              fprintf(fp0, "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", E->T[i], E->C[i], E->V[1][i], E->V[2][i], E->Vi[i], E->Cdot[i], E->Cphasedotnum[i], E->Cphasedot[i], E->Cphase_node[i], E->d_node[i]);
            }
          }
          else if (E->control.phasefile_C || E->control.phasefile_Complete)
          {
            fprintf(fp0, "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", E->T[i], E->C[i], E->V[1][i], E->V[2][i], E->Vi[i], E->C_phasefile_nno[0][i], E->C_phasefile_nno[1][i], E->C_phasefile_nno[2][i]);
          }
          else
          {
            fprintf(fp0, "%.4e %.4e %.4e %.4e %.4e\n", E->T[i], E->C[i], E->V[1][i], E->V[2][i], E->Vi[i]);
          }
        }
        for (i = 1; i <= E->lmesh.nel; i++)
          fprintf(fp0, "%.4e\n", E->P[i]);
      }
      else
        for (i = 1; i <= E->lmesh.nno; i++)
          fprintf(fp0, "%.4e %.4e %.4e\n", E->T[i], E->C[i], E->Vi[i]);
    }
    else
    {
      if ((E->advection.timesteps % (E->control.record_every)) == 0)
      {
        for (i = 1; i <= E->lmesh.nno; i++)
          fprintf(fp0, "%.4e %.4e %.4e %.4e\n", E->T[i], E->V[1][i], E->V[2][i], E->Vi[i]);
        for (i = 1; i <= E->lmesh.nel; i++)
          fprintf(fp0, "%.4e\n", E->P[i]);
      }
      else
        for (i = 1; i <= E->lmesh.nno; i++)
          fprintf(fp0, "%.4e %.4e\n", E->T[i], E->Vi[i]);
    }

    fclose(fp0);
  }

  if (E->control.composition && (E->advection.timesteps % (1 * E->control.record_every)) == 0)
  {
    sprintf(output_file, "%s/traces.%d.%d", E->control.data_file, E->parallel.me, file_number);
    fp0 = fopen(output_file, "w");
    num1 = E->advection.markers;
    MPI_Allreduce(&num1, &num2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    fprintf(fp0, "%d %d %d %.5e\n", E->advection.markers, num2, E->advection.timesteps, E->monitor.elapsed_time);
    for (i = 1; i <= E->advection.markers; i++)
    {
      // fprintf(fp0,"%.5e %.5e %d %d\n",E->XMC[1][i],E->XMC[2][i],E->CElement[i],E->C12[i]);
      if (E->control.phasevisc_C)
      {
        if (E->control.phasevisc_d)
        {
          fprintf(fp0, "%.5e %.5e %d %d %.5e %.5e %.5e %.5e %.5e %.5e\n", E->XMC[1][i], E->XMC[2][i], E->CElement[i], E->C12[i], E->Cphase_marker[i], E->Cphase_marker_old[i], E->d_marker[i], E->d_marker_old[i], E->d_dotnum[i], E->d_dot[i]);
        }
      }
      else if (E->control.phasefile_C || E->control.phasefile_Complete)
      {
        fprintf(fp0, "%.5e %.5e %d %d %d\n", E->XMC[1][i], E->XMC[2][i], E->CElement[i], E->C12[i], E->C_phasefile_marker_int[0][i]);
      }
      else
        fprintf(fp0, "%.5e %.5e %d %d\n", E->XMC[1][i], E->XMC[2][i], E->CElement[i], E->C12[i]);
    }
    for (i = 1; i <= E->lmesh.nel; i++)
    {
      //      fprintf(fp0,"%.4e\n",E->CE[i]);
      if (E->control.phasevisc_C)
      {
        fprintf(fp0, "%.4e %.4e\n", E->CE[i], E->CphaseE[i]);
      }
      else
        fprintf(fp0, "%.4e\n", E->CE[i]);
    }
    fclose(fp0);
    if (E->control.phasefile)
    {
      output_phasefile(E, file_number);
    }
  }

  if ((E->advection.timesteps % (E->control.record_every)) == 0)
  {

    surf = 0.0;
    botm = 0.0;
    for (i = 1; i <= E->lmesh.nox; i++)
    {
      j = i * E->lmesh.noz;
      if (i > 1)
      {
        surf += (E->slice.shflux[i] + E->slice.shflux[i - 1]) * 0.5 *
                (E->X[1][j] - E->X[1][j - E->lmesh.noz]);
        botm += (E->slice.bhflux[i] + E->slice.bhflux[i - 1]) * 0.5 *
                (E->X[1][j] - E->X[1][j - E->lmesh.noz]);
      }
    }
    surf = surf / E->X[1][E->lmesh.nno];
    botm = botm / E->X[1][E->lmesh.nno];

    for (i = 1; i <= E->lmesh.nno; i++)
      SV[i] = sqrt(E->V[1][i] * E->V[1][i] + E->V[2][i] * E->V[2][i]);

    return_horiz_ave(E, SV, E->Have.vrms);
    return_horiz_ave(E, E->Vi, E->Have.Vi);
    return_horiz_ave(E, E->T, E->Have.T);
    if (E->control.composition)
      return_horiz_ave(E, E->C, E->Have.Rho);

    sprintf(output_file, "%s/geoid.%d", E->control.data_file, file_number);
    fp1 = fopen(output_file, "w");
    fprintf(fp1, "%d %.5e %d\n", E->advection.timesteps, E->monitor.elapsed_time, E->control.llmax);
    for (i = 0; i < E->control.llmax; i++)
    {
      fprintf(fp1, "%d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", i + 1, E->harm_geoid[0][i], E->harm_geoid[1][i], E->harm_geoid_from_bncy[0][i], E->harm_geoid_from_bncy[1][i], E->harm_geoid_from_tpgt[0][i], E->harm_geoid_from_tpgt[1][i], E->harm_geoid_from_tpgb[0][i], E->harm_geoid_from_tpgb[1][i], E->harm_tpgt[0][i], E->harm_tpgt[1][i], E->harm_tpgb[0][i], E->harm_tpgb[1][i]);
    }
    fclose(fp1);

    sprintf(output_file, "%s/topo_hf.%d.%d", E->control.data_file, E->parallel.me, file_number);
    fp1 = fopen(output_file, "w");

    fprintf(fp1, "%6d %6d %.5e %.5e %.5e %.5e\n", E->lmesh.nno, E->advection.timesteps, E->monitor.elapsed_time, surf, botm, E->rad_heat.total);
    for (i = 1; i <= E->lmesh.noz; i++)
      fprintf(fp1, "%.4e %.5e %.5e %.5e %.5e\n", E->X[2][i], E->Have.T[i], E->Have.Rho[i], E->Have.Vi[i], E->Have.vrms[i]);
    for (i = 1; i <= E->lmesh.nox; i++)
    {
      j = i * E->lmesh.noz;
      fprintf(fp1, "%.4e %.5e %.5e %.5e %.5e\n", E->X[1][j], E->slice.tpg[i], E->slice.tpgb[i], E->slice.shflux[i], E->slice.bhflux[i]);
    }

    fclose(fp1);
    /*  ouput strainrate */
    eedot = (double *)malloc((2 + E->lmesh.nel) * sizeof(double));
    nedot = (double *)malloc((2 + nno) * sizeof(double));
    eedot11 = (double *)malloc((2 + E->lmesh.nel) * sizeof(double));
    nedot11 = (double *)malloc((2 + nno) * sizeof(double));
    eedot22 = (double *)malloc((2 + E->lmesh.nel) * sizeof(double));
    nedot22 = (double *)malloc((2 + nno) * sizeof(double));
    eedot12 = (double *)malloc((2 + E->lmesh.nel) * sizeof(double));
    nedot12 = (double *)malloc((2 + nno) * sizeof(double));

    edot = (double *)malloc((E->lmesh.NEL[E->mesh.levmax] + 2) * vpoints[E->mesh.nsd] * sizeof(double));
    edot11 = (double *)malloc((E->lmesh.NEL[E->mesh.levmax] + 2) * vpoints[E->mesh.nsd] * sizeof(double));
    edot22 = (double *)malloc((E->lmesh.NEL[E->mesh.levmax] + 2) * vpoints[E->mesh.nsd] * sizeof(double));
    edot12 = (double *)malloc((E->lmesh.NEL[E->mesh.levmax] + 2) * vpoints[E->mesh.nsd] * sizeof(double));

    strain_rate_2_inv_moreout(E, eedot, eedot11, eedot22, eedot12, 1);

    for (i = 1; i <= E->lmesh.nel; i++)
    {
      temp = eedot[i];
      temp11 = eedot11[i];
      temp22 = eedot22[i];
      temp12 = eedot12[i];
      for (jj = 1; jj <= vpoints[E->mesh.nsd]; jj++)
      {
        /*for(kk=1;kk<=ends;kk++) {
                temp += eedot

            }*/
        edot[(i - 1) * vpoints[E->mesh.nsd] + jj] = temp;
        edot11[(i - 1) * vpoints[E->mesh.nsd] + jj] = temp11;
        edot22[(i - 1) * vpoints[E->mesh.nsd] + jj] = temp22;
        edot12[(i - 1) * vpoints[E->mesh.nsd] + jj] = temp12;
      }
    }
    visc_from_gint_to_nodes(E, edot, nedot, E->mesh.levmax);
    visc_from_gint_to_nodes(E, edot11, nedot11, E->mesh.levmax);
    visc_from_gint_to_nodes(E, edot22, nedot22, E->mesh.levmax);
    visc_from_gint_to_nodes(E, edot12, nedot12, E->mesh.levmax);

    sprintf(output_file, "%s/strainrate.%d.%d", E->control.data_file, E->parallel.me, file_number);
    fp1 = fopen(output_file, "w");

    for (n = 1; n <= nno; n++)
      fprintf(fp1, "%.4e %.4e %.4e %.4e\n", nedot[n], nedot11[n], nedot22[n], nedot12[n]);
    fclose(fp1);
    free((void *)eedot);
    free((void *)nedot);
    free((void *)eedot11);
    free((void *)nedot11);
    free((void *)eedot22);
    free((void *)nedot22);
    free((void *)eedot12);
    free((void *)nedot12);
    free((void *)edot);
    free((void *)edot11);
    free((void *)edot22);
    free((void *)edot12);

    /* ouput strainrate end */
  }

  return;
}

void output_plumes(E, file_number) struct All_variables *E;
int file_number;
{
  static FILE *fp0;
  int i, n1, n2;
  char output_file[255];
  static int been_here = 0;

  const int nno = E->lmesh.nno;

  if (been_here == 0)
  {
    sprintf(output_file, "%s/plume", E->control.data_file);
    fp0 = fopen(output_file, "w");
    been_here++;
  }

  fprintf(fp0, "%6d %6d %.5e\n", E->lmesh.nox, E->advection.timesteps, E->monitor.elapsed_time);
  for (i = 1; i <= E->lmesh.nox; i++)
  {
    n1 = (i - 1) * E->lmesh.noz + E->control.PLUME;
    n2 = (i - 1) * E->lmesh.noz + E->control.SLAB;
    fprintf(fp0, "%.4e %.4e %.4e %.4e\n", E->T[n1], E->T[n2], E->V[1][n1], E->V[1][n2]);
  }
  fflush(fp0);

  return;
}

void output_phasefile(E, file_number) struct All_variables *E;
int file_number;
{
  static FILE *fp0;
  int i, n1, n2;
  char output_file[255];
  static int been_here = 0;
  int m, n, l, num, number1, number2, number3, number4, count_mineral, count_mat, temp_num_mineral;
  int left, right, mid;
  const int nno = E->lmesh.nno;
  int mark_P, mark_T, mark_P0_PREM = 160, mark_P_PREM;
  double coord_x, coord_z, temp, comp, coord_x_km, coord_z_km, T_C, P_GPa;
  double delta, delta1, delta2, delta_P, delta_T, delta_1, delta_2, delta_3, delta_4;
  double den_basa, Vp_basa, Vs_basa, den_pyro, Vp_pyro, Vs_pyro, den_harz, Vp_harz, Vs_harz;
  double r1, t1, loc_mid, tempdist;
  /*note here mark_p0=158 because CMB pressure  */
  static int mark_P_store[1001], mark_P_PREM_store[1001];
  double den_mineral[100], Vp_mineral[100], Vs_mineral[100];
  if (!been_here)
  {
    for (n = 1; n <= E->lmesh.noz; n++)
    {
      mark_P_store[n] = E->control.phasefile_noP;
      mark_P_PREM_store[n] = mark_P0_PREM;
    }
    mark_P = E->control.phasefile_noP;
    mark_P_PREM = mark_P0_PREM;

    for (n = 1; n <= E->lmesh.noz; n++)
    {
      num = n;
      coord_x = E->X[1][num];
      coord_z = E->X[2][num];
      temp = E->T[num];
      comp = E->C[num];
      coord_x_km = coord_x * E->data.layer_km;
      coord_z_km = (1.0 - coord_z) * E->data.layer_km;
      for (l = mark_P_PREM; l >= 2; l--)
      {
        if (coord_z_km <= E->control.phase_PREM_depth[l] && coord_z_km >= E->control.phase_PREM_depth[l - 1])
        {
          mark_P_PREM_store[n] = l;
          break;
        }
      }
      mark_P_PREM = mark_P_PREM_store[n];
      delta = E->control.phase_PREM_depth[mark_P_PREM] - E->control.phase_PREM_depth[mark_P_PREM - 1];
      delta1 = E->control.phase_PREM_depth[mark_P_PREM] - coord_z_km;
      delta2 = coord_z_km - E->control.phase_PREM_depth[mark_P_PREM - 1];
      P_GPa = E->control.phase_PREM_P[mark_P_PREM] * delta2 / delta + E->control.phase_PREM_P[mark_P_PREM - 1] * delta1 / delta;
      for (l = mark_P; l >= 2; l--)
      {
        if (P_GPa <= E->control.phase_P[l] && P_GPa >= E->control.phase_P[l - 1])
        {
          mark_P_store[n] = l;
          break;
        }
      }
      if (n == E->lmesh.noz && E->parallel.me_loc[2] == E->parallel.nprocz - 1)
      {
        mark_P_PREM_store[n] = 2;
        mark_P_store[n] = 2;
      }

      mark_P = mark_P_store[n];
//      fprintf(stderr, "%d %d %d\n", n, mark_P_PREM_store[n], mark_P_store[n]);
    }
    been_here++;
  }

  //  if (been_here==0)  {
  sprintf(output_file, "%s/phasefile_all.%d.%d", E->control.data_file, file_number, E->parallel.me);
  fp0 = fopen(output_file, "w");
  //    been_here++;
  //   }
  fprintf(fp0, "%6d %6d %.5e\n", E->lmesh.nno, E->advection.timesteps, E->monitor.elapsed_time);
  for (m = 1; m <= E->lmesh.nox; m++)
  {
    mark_P = E->control.phasefile_noP;
    mark_T = E->control.phasefile_noT;
    mark_P_PREM = mark_P0_PREM;
    for (n = 1; n <= E->lmesh.noz; n++)
    {
      num = n + E->lmesh.noz * (m - 1);
      mark_P = mark_P_store[n];
      mark_P_PREM = mark_P_PREM_store[n];

      coord_x = E->X[1][num];
      coord_z = E->X[2][num];
      temp = E->T[num];
      comp = E->C[num];
      coord_x_km = coord_x * E->data.layer_km;
      coord_z_km = (1.0 - coord_z) * E->data.layer_km;
      if (coord_z_km < 0.0)
        coord_z_km = 0.0;
      if (coord_z_km < (1.0 - E->viscosity.zlith) * E->data.layer_km)
      {
        T_C = temp * E->data.ref_temperature;
      }
      else if (coord_z_km < (1 - E->viscosity.zlm) * E->data.layer_km)
      {
        /* adiabatic temperature gradient is 0.5 K/km in the upper mantle, 0.3 K/km in the lower mantle */
        T_C = temp * E->data.ref_temperature + (coord_z_km - (1.0 - E->viscosity.zlith) * E->data.layer_km) * E->control.adi_um;
      }
      else
      {
        T_C = temp * E->data.ref_temperature + (coord_z_km - (1.0 - E->viscosity.zlm) * E->data.layer_km) * E->control.adi_lm + ((1.0 - E->viscosity.zlm) * E->data.layer_km - (1.0 - E->viscosity.zlith) * E->data.layer_km) * E->control.adi_um;
      }

      /* trace from CMB to surface and get pressure from PREM density*/
      /*            for(l=mark_P_PREM;l>=2;l--) {
                if(coord_z_km <= E->control.phase_PREM_depth[l] && coord_z_km >= E->control.phase_PREM_depth[l-1]) {
                    mark_P_PREM = l;
                    break; 
                }
            }       */
      //fprintf(stderr,"l=%d coord_z_km=%.4e  E->control.phase_PREM_depth[l] = %.4e\n",l,coord_z_km, E->control.phase_PREM_depth[l]);
      delta = E->control.phase_PREM_depth[mark_P_PREM] - E->control.phase_PREM_depth[mark_P_PREM - 1];
      delta1 = E->control.phase_PREM_depth[mark_P_PREM] - coord_z_km;
      delta2 = coord_z_km - E->control.phase_PREM_depth[mark_P_PREM - 1];
      P_GPa = E->control.phase_PREM_P[mark_P_PREM] * delta2 / delta + E->control.phase_PREM_P[mark_P_PREM - 1] * delta1 / delta;
      /* after finish calculating T in C and P in GPa, get location in the precalculated phase diagram */
      /*            for(l=mark_P;l>=2;l--) {
                if(P_GPa<=E->control.phase_P[l] && P_GPa>=E->control.phase_P[l-1]) {
                    mark_P = l;
                    break;
                }
            } */
      /* note T is not ordered from max to min, so should use O(logn) algorithm instead of O(1) */
      left = 2;
      right = E->control.phasefile_noT;
      mid = (E->control.phasefile_noT + 1) / 2;
      while (left < right - 1)
      {
        if (T_C < E->control.phase_T[mid])
        {
          right = mid;
        }
        else
        {
          left = mid;
        }
        mid = (left + right) / 2;
      }
      if (T_C < E->control.phase_T[mid])
      {

        mark_T = mid;
      }
      else
      {
        mark_T = mid + 1;
      }

      //fprintf(stderr,"mark_P_PREm=%d mark_P=%d mark_T=%d\n",mark_P_PREM,mark_P,mark_T);
      //fprintf(stderr,"coord_z_km=%.4e T_C=%.4e P_GPa=%.4e\n",coord_z_km,T_C,P_GPa);

      delta_P = E->control.phase_P[mark_P] - E->control.phase_P[mark_P - 1];
      delta_1 = E->control.phase_P[mark_P] - P_GPa;
      delta_2 = P_GPa - E->control.phase_P[mark_P - 1];
      delta_T = E->control.phase_T[mark_T] - E->control.phase_T[mark_T - 1];
      delta_3 = E->control.phase_T[mark_T] - T_C;
      delta_4 = T_C - E->control.phase_T[mark_T - 1];
      number1 = mark_P - 1 + E->control.phasefile_noP * (mark_T - 1 - 1);
      number2 = mark_P + E->control.phasefile_noP * (mark_T - 1 - 1);
      number3 = mark_P - 1 + E->control.phasefile_noP * (mark_T - 1);
      number4 = mark_P + E->control.phasefile_noP * (mark_T - 1);
      if (delta_2 < 0)
      {
        delta_1 = 1.0;
        delta_2 = 0.0;
        delta_P = 1.0;
      }
      if (delta_4 < 0)
      {
        delta_3 = 1.0;
        delta_4 = 0.0;
        delta_T = 1.0;
      }

      /* now calculate density from phase diagram */
      // 1 basalt, 0 pyrolite, 2 harzburgite
      den_mineral[1] = E->control.phase_basa_density[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_basa_density[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_basa_density[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_basa_density[number4] * delta_2 / delta_P * delta_4 / delta_T;
      Vp_mineral[1] = E->control.phase_basa_Vp[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_basa_Vp[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_basa_Vp[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_basa_Vp[number4] * delta_2 / delta_P * delta_4 / delta_T;
      Vs_mineral[1] = E->control.phase_basa_Vs[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_basa_Vs[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_basa_Vs[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_basa_Vs[number4] * delta_2 / delta_P * delta_4 / delta_T;
      den_mineral[0] = E->control.phase_pyro_density[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_pyro_density[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_pyro_density[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_pyro_density[number4] * delta_2 / delta_P * delta_4 / delta_T;
      Vp_mineral[0] = E->control.phase_pyro_Vp[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_pyro_Vp[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_pyro_Vp[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_pyro_Vp[number4] * delta_2 / delta_P * delta_4 / delta_T;
      Vs_mineral[0] = E->control.phase_pyro_Vs[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_pyro_Vs[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_pyro_Vs[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_pyro_Vs[number4] * delta_2 / delta_P * delta_4 / delta_T;
      den_mineral[2] = E->control.phase_harz_density[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_harz_density[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_harz_density[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_harz_density[number4] * delta_2 / delta_P * delta_4 / delta_T;
      Vp_mineral[2] = E->control.phase_harz_Vp[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_harz_Vp[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_harz_Vp[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_harz_Vp[number4] * delta_2 / delta_P * delta_4 / delta_T;
      Vs_mineral[2] = E->control.phase_harz_Vs[number1] * delta_1 / delta_P * delta_3 / delta_T + E->control.phase_harz_Vs[number2] * delta_2 / delta_P * delta_3 / delta_T + E->control.phase_harz_Vs[number3] * delta_1 / delta_P * delta_4 / delta_T + E->control.phase_harz_Vs[number4] * delta_2 / delta_P * delta_4 / delta_T;

      E->density_phase[num] = den_mineral[1] * E->C[num] + den_mineral[2] * (1.0 - E->C[num]);
      E->Vp_phase[num] = Vp_mineral[1] * E->C[num] + Vp_mineral[2] * (1.0 - E->C[num]);
      E->Vs_phase[num] = Vs_mineral[1] * E->C[num] + Vs_mineral[2] * (1.0 - E->C[num]);

      if (E->control.phasefile_C || E->control.phasefile_Complete)
      {
        for (count_mat = 0; count_mat < E->control.phasefile_C_flavor; count_mat++)
        {
          if (count_mat == 0)
          {
            E->density_phase[num] = 0.0;
            E->Vp_phase[num] = 0.0;
            E->Vs_phase[num] = 0.0;
          }
          temp_num_mineral = E->control.phasefile_C_mat_mineral[count_mat];
          E->density_phase[num] += den_mineral[temp_num_mineral] * E->C_phasefile_nno[count_mat][num];
          E->Vp_phase[num] += Vp_mineral[temp_num_mineral] * E->C_phasefile_nno[count_mat][num];
          E->Vs_phase[num] += Vs_mineral[temp_num_mineral] * E->C_phasefile_nno[count_mat][num];
        }
      }

      /* continent correction */
      if (E->control.phasefile_buoyancy_correction)
      {

        if (E->control.phasefile_buoyancy_continent > 0.0)
        {
          if (E->control.trechmigrate)
          {
            loc_mid = E->control.velo_surf_loc_mid;
            loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
            t1 = E->X[1][num];
            r1 = E->X[2][num];
            tempdist = (t1 - loc_mid) * tan(E->control.dip_margin * 3.14159265 / 180.0) + 1.0 - r1;
            if (tempdist <= (E->viscosity.right_weakzone_platebond) * tan(E->control.dip_margin * 3.14159265 / 180.0) && r1 >= E->viscosity.zlith)
            {
              E->density_phase[num] -= E->control.phasefile_buoyancy_continent;
            }
          }
        } /* end continent correction */
        /* crust correction */
        if (E->control.phasefile_buoyancy_crust > 0.0)
        {
          if (coord_z_km <= E->control.phasefile_buoyancy_crust_depth)
          {
            if (E->C[num] > 0.0)
            {
              E->density_phase[num] -= E->control.phasefile_buoyancy_crust;
            }
          }
        }
      }

      fprintf(fp0, "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", E->density_phase[num], E->Vp_phase[num], E->Vs_phase[num], den_mineral[0], den_mineral[1], Vp_mineral[0], Vp_mineral[1], Vs_mineral[0], Vs_mineral[1]);
    } // end of noz
  }   // end of nox
  fflush(fp0);
  return;
}

/*      ----------------------------------- */
void coordinates_dx(E) struct All_variables *E;
{

  int i;

  E->ibm_dx.x1 = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));
  E->ibm_dx.x2 = (double *)malloc((E->lmesh.nno + 1) * sizeof(double));

  E->ibm_dx.nox = E->lmesh.nox;
  E->ibm_dx.noz = E->lmesh.noz;

  for (i = 1; i <= E->lmesh.nno; i++)
  {
    E->ibm_dx.x1[i] = E->X[2][i] * sin(E->X[1][i]);
    E->ibm_dx.x2[i] = E->X[2][i] * cos(E->X[1][i]);
  }

  return;
}

/*      ----------------------------------- */
void frames_for_dx(E, file_number) struct All_variables *E;
int file_number;
{

  int i;
  static int nframe = 0;
  FILE *fp;
  char output_file[255];
  const double offset1 = -1.2;
  const double offset2 = 0.2;

  nframe++;

  sprintf(output_file, "%s/mv.%03d.dx", E->control.data_file, nframe);
  fp = fopen(output_file, "w");
  for (i = 1; i <= E->lmesh.nno; i++)
    fprintf(fp, "%.4e %.4e %.4e\n", E->ibm_dx.x1[i] + offset1, E->ibm_dx.x2[i], E->T[i]);
  fclose(fp);

  sprintf(output_file, "%s/nv.%03d.dx", E->control.data_file, nframe);
  fp = fopen(output_file, "w");
  for (i = 1; i <= E->lmesh.nno; i++)
    fprintf(fp, "%.4e %.4e %.4e\n", E->ibm_dx.x1[i] + offset2, E->ibm_dx.x2[i], E->C[i]);
  fclose(fp);

  sprintf(output_file, "%s/mv.%03d.general", E->control.data_file, nframe);
  fp = fopen(output_file, "w");
  fprintf(fp, "file = mv.%03d.dx\n", nframe);
  fprintf(fp, "grid = %2d x %2d\n", E->ibm_dx.nox, E->ibm_dx.noz);
  fprintf(fp, "format = ascii\n");
  fprintf(fp, "interleaving = field\n");
  fprintf(fp, "majority = row\n");
  fprintf(fp, "field = locations, field0\n");
  fprintf(fp, "structure = 2-vector, scalar\n");
  fprintf(fp, "type = double, double\n");
  fprintf(fp, "\n");
  fprintf(fp, "end\n");
  fclose(fp);

  sprintf(output_file, "%s/nv.%03d.general", E->control.data_file, nframe);
  fp = fopen(output_file, "w");
  fprintf(fp, "file = nv.%03d.dx\n", nframe);
  fprintf(fp, "grid = %2d x %2d\n", E->ibm_dx.nox, E->ibm_dx.noz);
  fprintf(fp, "format = ascii\n");
  fprintf(fp, "interleaving = field\n");
  fprintf(fp, "majority = row\n");
  fprintf(fp, "field = locations, field0\n");
  fprintf(fp, "structure = 2-vector, scalar\n");
  fprintf(fp, "type = double, double\n");
  fprintf(fp, "\n");
  fprintf(fp, "end\n");
  fclose(fp);

  return;
}

void output_velo_related_binary(E, file_number) struct All_variables *E;
int file_number;
{
  int el, els, i, j, k, ii, m, node, fd;
  int nox, noz, noy, nfx, nfz, nfy1, nfy2, size1, size2;
  char output_file[255];
  static double *SV, *EV;
  static int been_here = 0;
  void output_phasefile();
  void get_surface_velo();
  void get_ele_visc();
  const int nno = E->lmesh.nno;
  /*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.velo",E->control.data_file);
    E->filed[10]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_t",E->control.data_file);
    E->filed[11]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_b",E->control.data_file);
    E->filed[12]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.visc",E->control.data_file);
    E->filed[13]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.fas670",E->control.data_file);
    E->filed[14]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.stress",E->control.data_file);
    E->filed[9]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0)  {
    ii = E->mesh.nsf;
    SV = (double *) malloc ((2*ii+2)*sizeof(double));

    size2 = (E->mesh.nel+1)*sizeof(double);
    EV = (double *) malloc (size2);
    been_here++;
    }

  ii = E->mesh.nsf;
  size2 = 2*(ii+2)*sizeof(double);
    get_surface_velo (E,SV);
  write(E->filed[10],SV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(double);
  write(E->filed[11],E->slice.tpg,size2);
  write(E->filed[12],E->slice.tpgb,size2);

  size2 = (E->mesh.nel+1)*sizeof(double);
    get_ele_visc (E,EV);
  write(E->filed[13],EV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(double);
  write(E->filed[14],E->Fas670_b,size2);

  size2 = (2*E->mesh.nsf+1)*sizeof(double);
  write(E->filed[9],E->stress,size2);
*/
  return;
}

/* ====================================================================== */

void output_temp(E, file_number) struct All_variables *E;
int file_number;
{
  int nno, i, j, fd;
  static int *temp1;
  static int been_here = 0;
  static int size2, size1;
  char output_file[255];
  /*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.temp",E->control.data_file);
    E->filed[5]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0) {
    temp1 = (int *) malloc ((E->mesh.noy*6)*sizeof(int));

    sprintf(output_file,"%s.mesh",E->control.data_file);
    E->filed[1]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.x",E->control.data_file);
    E->filed[2]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.z",E->control.data_file);
    E->filed[3]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.y",E->control.data_file);
    E->filed[4]=open(output_file,O_RDWR | O_CREAT, 0644);

    size1 = (E->mesh.noy*6)*sizeof(int);
    size2= (E->mesh.nno+1)*sizeof(double);

    temp1[1] = E->mesh.nno;
    temp1[3] = size2;
    temp1[5] = E->mesh.nsf;
    temp1[6] = E->mesh.nel;

        write(E->filed[1],temp1,size1);
        write(E->filed[2],E->X[1],size2);
        write(E->filed[3],E->X[2],size2);
        write(E->filed[4],E->X[3],size2);

    close(E->filed[1]);
    close(E->filed[2]);
    close(E->filed[3]);
    close(E->filed[4]);

    been_here++;
    }


    write(E->filed[5],E->T,size2);
*/
  return;
}

void output_strainrate(E, file_number) struct All_variables *E;
int file_number;
{
  int i, j, fd, n;
  static int *temp1;
  static int been_here = 0;
  static int size2, size1;
  char output_file[255];
  double one, two, scale, stress_magnitude, depth, exponent1;
  double *eedot, *nedot, SS[5];
  double *eedot11, *nedot11, *eedot22, *nedot22, *eedot12, *nedot12;
  FILE *fp1;
  void strain_rate_2_inv_moreout();
  int e, l, z, jj, kk;
  double temp0;
  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int nno = E->lmesh.nno;
  const int ends = enodes[E->lmesh.nsd];
  void visc_from_gint_to_nodes();
  eedot = (double *)malloc((2 + nel) * sizeof(double));
  nedot = (double *)malloc((2 + nno) * sizeof(double));
  eedot11 = (double *)malloc((2 + nel) * sizeof(double));
  nedot11 = (double *)malloc((2 + nno) * sizeof(double));
  eedot22 = (double *)malloc((2 + nel) * sizeof(double));
  nedot22 = (double *)malloc((2 + nno) * sizeof(double));
  eedot12 = (double *)malloc((2 + nel) * sizeof(double));
  nedot12 = (double *)malloc((2 + nno) * sizeof(double));

  one = 1.0;
  two = 2.0;
  strain_rate_2_inv_moreout(E, eedot, eedot11, eedot22, eedot12, 1);
  sprintf(output_file, "%s/strainrate.%d.%d", E->control.data_file, file_number, E->parallel.me);
  fp1 = fopen(output_file, "w");
  fprintf(fp1, "%6d %6d %.5e\n", E->lmesh.nno, E->advection.timesteps, E->monitor.elapsed_time);
  visc_from_gint_to_nodes(E, eedot, nedot, E->mesh.levmax);
  visc_from_gint_to_nodes(E, eedot11, nedot11, E->mesh.levmax);
  visc_from_gint_to_nodes(E, eedot22, nedot22, E->mesh.levmax);
  visc_from_gint_to_nodes(E, eedot12, nedot12, E->mesh.levmax);
  for (n = 1; n <= nno; n++)
    fprintf(fp1, "%.4e %.4e %.4e %.4e\n", nedot[n], nedot11[n], nedot22[n], nedot12[n]);
  fclose(fp1);
  free((void *)eedot);
  free((void *)nedot);
  free((void *)eedot11);
  free((void *)nedot11);
  free((void *)eedot22);
  free((void *)nedot22);
  free((void *)eedot12);
  free((void *)nedot12);

  return;
}
