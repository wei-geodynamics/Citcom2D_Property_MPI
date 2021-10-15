#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void phase_change(E, Bb, Bb_b, Bt, Bt_b) struct All_variables *E;
double *Bb, *Bb_b, *Bt, *Bt_b;
{
  static int been_here = 0;

  FILE *fp1, *fp2;
  char output_file[255];

  int i, j, k, n, ns;
  double e_pressure, temp;

  static double pt5 = 0.5;
  static double one = 1.0;

  if (been_here++ == 0)
  {

    E->control.Ra_670 = E->control.Ra_670 * E->control.Ra_temp / (E->data.density * E->data.therm_exp * E->data.ref_temperature);
    E->control.Ra_410 = E->control.Ra_410 * E->control.Ra_temp / (E->data.density * E->data.therm_exp * E->data.ref_temperature);

    E->control.clapeyron670 = E->control.clapeyron670 * E->data.ref_temperature /
                              (E->data.density * E->data.grav_acc * E->data.layer_meter);
    E->control.clapeyron410 = E->control.clapeyron410 * E->data.ref_temperature /
                              (E->data.density * E->data.grav_acc * E->data.layer_meter);

    E->control.width670 = E->data.layer_meter / E->control.width670;
    E->control.width410 = E->data.layer_meter / E->control.width410;

    E->control.transT670 = E->control.transT670 / E->data.ref_temperature;
    E->control.transT410 = E->control.transT410 / E->data.ref_temperature;

    fprintf(E->fp, "Rab410 670=%lf %lf Clap410 670=%lf %lf Di=%lf %lf \n", E->control.Ra_410, E->control.Ra_670, E->control.clapeyron410, E->control.clapeyron670, E->data.disptn_number, E->data.ref_viscosity);
  }

  for (i = 1; i <= E->lmesh.nno; i++)
  {
    e_pressure = E->viscosity.zlm - E->X[2][i] -
                 E->control.clapeyron670 * (E->T[i] - E->control.transT670);
    Bb[i] = pt5 * (one + tanh(E->control.width670 * e_pressure));
  }

  for (i = 1; i <= E->lmesh.nno; i++)
  {
    e_pressure = E->viscosity.z410 - E->X[2][i] -
                 E->control.clapeyron410 * (E->T[i] - E->control.transT410);
    Bt[i] = pt5 * (one + tanh(E->control.width410 * e_pressure));
  }

  if (E->advection.timesteps % E->control.record_every == 0)
  {

    for (j = 1; j <= E->lmesh.nox; j++)
    {
      Bb_b[j] = 0.0;
      for (i = 1; i < E->lmesh.noz; i++)
      {
        n = (j - 1) * E->lmesh.noz + i;
        if (Bb[n] >= pt5 && Bb[n + 1] <= pt5)
        {
          Bb_b[j] = (E->X[2][n + 1] - E->X[2][n]) * (pt5 - Bb[n]) / (Bb[n + 1] - Bb[n]) + E->X[2][n];
          break;
        }
      }
    }

    for (j = 1; j <= E->lmesh.nox; j++)
    {
      Bt_b[j] = 0.0;
      for (i = 1; i < E->lmesh.noz; i++)
      {
        n = (j - 1) * E->lmesh.noz + i;
        if (Bt[n] >= pt5 && Bt[n + 1] <= pt5)
        {
          Bt_b[j] = (E->X[2][n + 1] - E->X[2][n]) * (pt5 - Bt[n]) / (Bt[n + 1] - Bt[n]) + E->X[2][n];
          break;
        }
      }
    }

    sprintf(output_file, "%s/fas.%d", E->control.data_file, E->advection.timesteps);
    fp1 = fopen(output_file, "w");
    /*  for (j=1;j<=E->mesh.nox;j++)
      fprintf(fp1,"%.4e %.5e %.5e %.5e %.5e\n",E->X[1][j*E->mesh.noz],Bt_b[j],Bb_b[j],Bt[i],Bb[i]); */
    fclose(fp1);
  }

  return;
}

void phase_change_basalt(E, Bb, Bb_b, Bt, Bt_b, Bb_new, Bb_new_b, Bb_all_b) struct All_variables *E;
double *Bb, *Bb_b, *Bt, *Bt_b, *Bb_new, *Bb_new_b, *Bb_all_b;
{
  static int been_here = 0;

  FILE *fp1, *fp2;
  char output_file[255];

  int i, j, k, n, ns;
  double e_pressure, temp;

  static double pt5 = 0.5;
  static double one = 1.0;

  if (been_here++ == 0)
  {

    E->control.Ra_670 = E->control.Ra_670 * E->control.Ra_temp / (E->data.density * E->data.therm_exp * E->data.ref_temperature);
    E->control.Ra_410 = E->control.Ra_410 * E->control.Ra_temp / (E->data.density * E->data.therm_exp * E->data.ref_temperature);
    E->control.Ra_670_basalt = E->control.Ra_670_basalt * E->control.Ra_temp / (E->data.density * E->data.therm_exp * E->data.ref_temperature);

    E->control.clapeyron670 = E->control.clapeyron670 * E->data.ref_temperature /
                              (E->data.density * E->data.grav_acc * E->data.layer_meter);
    E->control.clapeyron410 = E->control.clapeyron410 * E->data.ref_temperature /
                              (E->data.density * E->data.grav_acc * E->data.layer_meter);
    E->control.clapeyron670_basalt = E->control.clapeyron670_basalt * E->data.ref_temperature /
                                     (E->data.density * E->data.grav_acc * E->data.layer_meter);

    E->control.width670 = E->data.layer_meter / E->control.width670;
    E->control.width410 = E->data.layer_meter / E->control.width410;

    E->control.transT670 = E->control.transT670 / E->data.ref_temperature;
    E->control.transT410 = E->control.transT410 / E->data.ref_temperature;

    E->control.width670_basalt = E->data.layer_meter / E->control.width670_basalt;

    E->control.transT670_basalt = E->control.transT670_basalt / E->data.ref_temperature;

    fprintf(E->fp, "Rab410 670=%lf %lf Clap410 670=%lf %lf Di=%lf %lf \n", E->control.Ra_410, E->control.Ra_670, E->control.clapeyron410, E->control.clapeyron670, E->data.disptn_number, E->data.ref_viscosity);
  }

  for (i = 1; i <= E->lmesh.nno; i++)
  {
    e_pressure = E->viscosity.zlm - E->X[2][i] -
                 E->control.clapeyron670 * (E->T[i] - E->control.transT670);
    Bb[i] = pt5 * (one + tanh(E->control.width670 * e_pressure));
  }

  for (i = 1; i <= E->lmesh.nno; i++)
  {
    e_pressure = E->viscosity.z410 - E->X[2][i] -
                 E->control.clapeyron410 * (E->T[i] - E->control.transT410);
    Bt[i] = pt5 * (one + tanh(E->control.width410 * e_pressure));
  }

  for (i = 1; i <= E->lmesh.nno; i++)
  {
    e_pressure = E->viscosity.zbasalt - E->X[2][i] -
                 E->control.clapeyron670_basalt * (E->T[i] - E->control.transT670_basalt);
    Bb_new[i] = pt5 * (one + tanh(E->control.width670_basalt * e_pressure));
  }

  if (E->advection.timesteps % E->control.record_every == 0)
  {

    for (j = 1; j <= E->lmesh.nox; j++)
    {
      Bb_b[j] = 0.0;
      for (i = 1; i < E->lmesh.noz; i++)
      {
        n = (j - 1) * E->lmesh.noz + i;
        if (Bb[n] >= pt5 && Bb[n + 1] <= pt5)
        {
          Bb_b[j] = (E->X[2][n + 1] - E->X[2][n]) * (pt5 - Bb[n]) / (Bb[n + 1] - Bb[n]) + E->X[2][n];
          break;
        }
      }
    }

    for (j = 1; j <= E->lmesh.nox; j++)
    {
      Bt_b[j] = 0.0;
      for (i = 1; i < E->lmesh.noz; i++)
      {
        n = (j - 1) * E->lmesh.noz + i;
        if (Bt[n] >= pt5 && Bt[n + 1] <= pt5)
        {
          Bt_b[j] = (E->X[2][n + 1] - E->X[2][n]) * (pt5 - Bt[n]) / (Bt[n + 1] - Bt[n]) + E->X[2][n];
          break;
        }
      }
    }

    for (j = 1; j <= E->lmesh.nox; j++)
    {
      Bb_new_b[j] = 0.0;
      for (i = 1; i < E->lmesh.noz; i++)
      {
        n = (j - 1) * E->lmesh.noz + i;
        if (Bb_new[n] >= pt5 && Bb_new[n + 1] <= pt5)
        {
          Bb_new_b[j] = (E->X[2][n + 1] - E->X[2][n]) * (pt5 - Bb_new[n]) / (Bb_new[n + 1] - Bb_new[n]) + E->X[2][n];
          break;
        }
      }
    }

    for (j = 1; j <= E->lmesh.nox; j++)
    {
      Bb_all_b[j] = Bb_b[j];
      for (i = 1; i < E->lmesh.noz; i++)
      {
        n = (j - 1) * E->lmesh.noz + i;
        if (E->C[n] >= 0.5 && E->X[2][n] - Bb_new_b[j] >= -5.0 / 2870.0 && E->X[2][n] - Bb_new_b[j] <= 5.0 / 2870.0)
        {
          Bb_all_b[j] = Bb_new_b[j];
          break;
        }
      }
    }

    sprintf(output_file, "%s/fas.%d", E->control.data_file, E->advection.timesteps);
    fp1 = fopen(output_file, "w");
    for (j = 1; j <= E->lmesh.nox; j++)
      //    fprintf(fp1,"%.4e %.5e %.5e\n",E->X[1][j*E->mesh.noz],Bt_b[j],Bb_b[j]);
      fprintf(fp1, "%.4e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", E->X[1][j * E->mesh.noz], Bt_b[j], Bb_b[j], Bb_new_b[j], Bb_all_b[j], Bt[j], Bb[j], Bb_new[j]);
    fclose(fp1);

    /*
  sprintf(output_file,"%s/fas_no.%d",E->control.data_file,E->advection.timesteps);
  fp1=fopen(output_file,"w");
  for (j=1;j<=E->mesh.nno;j++)
      fprintf(fp1,"%.4e %.5e %.5e %.5e\n",E->X[2][j],Bt[j],Bb[j],Bb_new[j]);
  fclose(fp1);
*/
  }

  return;
}
