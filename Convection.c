/* Assumes parameter list is opened and reads the things it needs. 
   Variables are initialized etc, default values are set */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h> /* for "system" command */
#include <strings.h>

void set_convection_defaults(E) struct All_variables *E;
{
  void PG_timestep_with_melting();
  void PG_timestep();
  void PG_timestep_particle();
  void read_convection_settings();
  void convection_derived_values();
  void convection_allocate_memory();
  void convection_boundary_conditions();
  void node_locations();
  void convection_initial_fields();
  void twiddle_thumbs();

  input_int("composition", &(E->control.composition), "0");
  input_double("comp_diffusivity", &(E->control.comp_diff), "0");
  input_string("comp_adv_method", E->control.comp_adv_method, NULL);
  input_int("composition_phasechange", &(E->control.composition_phasechange), "0");

  if (E->control.composition && strcmp(E->control.comp_adv_method, "field") == 0)
    E->next_buoyancy_field = PG_timestep;
  else if (E->control.composition && strcmp(E->control.comp_adv_method, "particle") == 0)
  {
    E->next_buoyancy_field = PG_timestep_particle;
  }
  else
  {
    E->next_buoyancy_field = PG_timestep;
  }

  E->special_process_new_buoyancy = twiddle_thumbs;
  E->problem_settings = read_convection_settings;
  E->problem_derived_values = convection_derived_values;
  E->problem_allocate_vars = convection_allocate_memory;
  E->problem_boundary_conds = convection_boundary_conditions;
  E->problem_initial_fields = convection_initial_fields;
  E->problem_node_positions = node_locations;
  E->problem_update_node_positions = twiddle_thumbs;
  E->problem_update_bcs = twiddle_thumbs;

  sprintf(E->control.which_data_files, "Temp,Strf,Pres");
  sprintf(E->control.which_horiz_averages, "Temp,Visc,Vrms");
  sprintf(E->control.which_running_data, "Step,Time,");
  sprintf(E->control.which_observable_data, "Shfl");

  return;
}

void read_convection_settings(E) struct All_variables *E;

{
  void advection_diffusion_parameters();
  double density_diff;

  /* parameters */

  input_double("rayleigh", &(E->control.Ra_temp), "essential");

  E->data.ref_viscosity = E->data.grav_acc * E->data.density * E->data.therm_exp * E->data.ref_temperature * E->data.layer_meter * E->data.layer_meter * E->data.layer_meter / (E->control.Ra_temp * E->data.therm_diff);

  input_double("rayleigh_comp", &(E->control.Ra_comp), "essential");

  density_diff = E->control.Ra_comp * E->data.ref_viscosity * E->data.therm_diff / (E->data.grav_acc * E->data.layer_meter * E->data.layer_meter * E->data.layer_meter);

  fprintf(E->fp, "Ra_temp=%.5e Ra_comp=%.5e %.5e %.5e\n", E->control.Ra_temp, E->control.Ra_comp, E->data.ref_viscosity, density_diff);

  input_boolean("halfspace", &(E->convection.half_space_cooling), "off");
  input_double("halfspage", &(E->convection.half_space_age), "nodefault");

  input_int("temperature_blobs", &(E->convection.temp_blobs), "0");
  input_double_vector("temperature_blobx", E->convection.temp_blobs, E->convection.temp_blob_x);
  input_double_vector("temperature_bloby", E->convection.temp_blobs, E->convection.temp_blob_y);
  input_double_vector("temperature_blobz", E->convection.temp_blobs, E->convection.temp_blob_z);
  input_double_vector("temperature_blobsize", E->convection.temp_blobs, E->convection.temp_blob_radius);
  input_double_vector("temperature_blobDT", E->convection.temp_blobs, E->convection.temp_blob_T);
  input_double_vector("temperature_blobbg", E->convection.temp_blobs, E->convection.temp_blob_bg);
  input_int_vector("temperature_blobsticky", E->convection.temp_blobs, E->convection.temp_blob_sticky);

  input_int("temperature_zones", &(E->convection.temp_zones), "0");
  input_double_vector("temperature_zonex1", E->convection.temp_zones, E->convection.temp_zonex1);
  input_double_vector("temperature_zonex2", E->convection.temp_zones, E->convection.temp_zonex2);
  input_double_vector("temperature_zonez1", E->convection.temp_zones, E->convection.temp_zonez1);
  input_double_vector("temperature_zonez2", E->convection.temp_zones, E->convection.temp_zonez2);
  input_double_vector("temperature_zoney1", E->convection.temp_zones, E->convection.temp_zoney1);
  input_double_vector("temperature_zoney2", E->convection.temp_zones, E->convection.temp_zoney2);
  input_double_vector("temperature_zoney2", E->convection.temp_zones, E->convection.temp_zoney2);
  input_double_vector("temperature_zoney2", E->convection.temp_zones, E->convection.temp_zoney2);
  input_double_vector("temperature_zonehw", E->convection.temp_zones, E->convection.temp_zonehw);
  input_double_vector("temperature_zonemag", E->convection.temp_zones, E->convection.temp_zonemag);
  input_int_vector("temperature_zonesticky", E->convection.temp_zones, E->convection.temp_zone_sticky);

  input_int("num_perturbations", &(E->convection.number_of_perturbations), "0,0,32");
  input_double_vector("perturbmag", E->convection.number_of_perturbations, E->convection.perturb_mag);
  input_double_vector("perturbk", E->convection.number_of_perturbations, E->convection.perturb_k);

  input_string("prevT", E->convection.old_T_file, "initialize");

  advection_diffusion_parameters(E);

  if (E->control.restart)
  {
    input_int("restart_timesteps", &(E->monitor.solution_cycles), "0");
    input_string("oldfile", E->convection.old_T_file, "initialize");
  }

  return;
}

/* =================================================================
   Any setup which relates only to the convection stuff goes in here
   ================================================================= */

void convection_derived_values(E) struct All_variables *E;

{

  return;
}

void convection_allocate_memory(E) struct All_variables *E;

{
  void advection_diffusion_allocate_memory();

  advection_diffusion_allocate_memory(E);

  return;
}

/* ============================================ */

void convection_initial_fields(E) struct All_variables *E;

{
  int i;
  int me = E->parallel.me;
  void convection_initial_temperature();
  void convection_initial_markers();
  const int nno = E->lmesh.nno;
  const int nel = E->lmesh.nel;
  if (E->control.composition)
  {
    E->Cdot = (double *)malloc((nno + 1) * sizeof(double));
    for (i = 1; i <= nno; i++)
      E->Cdot[i] = 0.0;
    if (!(strcmp(E->control.comp_adv_method, "field") == 0))
    {
      E->advection.markers = E->advection.markers_per_ele * E->mesh.nel;
      E->advection.markers = E->advection.markers * E->lmesh.volume / E->mesh.volume;
      E->advection.markers_uplimit = E->advection.markers * 2;
      fprintf(stderr, "CPU %d %d %.4e %.4e", E->parallel.me, E->advection.markers, E->lmesh.volume, E->mesh.volume);
      for (i = 1; i <= E->mesh.nsd; i++)
      {
        E->VO[i] = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
        E->XMC[i] = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
        E->XMCpred[i] = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
        E->Vpred[i] = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
      }
      E->C12 = (int *)malloc((E->advection.markers_uplimit + 1) * sizeof(int));
      E->CElement = (int *)malloc((E->advection.markers_uplimit + 1) * sizeof(int));
      E->traces_leave = (int *)malloc((E->advection.markers_uplimit + 1) * sizeof(int));
    }

    if (E->control.phasefile_C || E->control.phasefile_Complete)
    {
      if (me == 0)
        fprintf(stderr, "construct marker array\n");
      /* creat array for tracer to store  in marker size */
      /* note first two are int for flavor */
      /* the next two are for initial processer number and marker ID */
      /* the first n double are not used yet */

      E->C_phasefile_markers_int_num_store = 4;

      E->C_phasefile_markers_double_num_store = 0;
      if (E->control.phasevisc_d)
      {
        E->C_phasefile_d_start_double = E->C_phasefile_markers_double_num_store;
        E->C_phasefile_d_start_element = E->control.phasefile_C_num_element;
        E->C_phasefile_d_start_nno = E->control.phasefile_C_num_nno;
        E->C_phasefile_markers_double_num_store += 4;
        E->control.phasefile_C_num_element += 4;
        E->control.phasefile_C_num_nno += 4;
      }

      if (E->C_phasefile_markers_double_num_store == 0)
      {
        E->C_phasefile_markers_double_num_store = 1;
      }

      for (i = 0; i < E->C_phasefile_markers_int_num_store; i++)
      {
        E->C_phasefile_marker_int[i] = (int *)malloc((E->advection.markers_uplimit + 1) * sizeof(int));
      }

      for (i = 0; i < E->C_phasefile_markers_double_num_store; i++)
      {
        E->C_phasefile_marker_double[i] = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
      }
      /* creat array for tracer to store  in element size */
      if (me == 0)
        fprintf(stderr, "construct marker element array\n");
      /* the first n are for ratio of each type particle */

      for (i = 0; i < E->control.phasefile_C_num_element; i++)
      {
        E->C_phasefile_element[i] = (double *)malloc((nel + 1) * sizeof(double));
      }
      if (me == 0)
        fprintf(stderr, "construct marker nno array\n");

      /* creat array for tracer to store  in nno size */
      /* the first n are for C component at nodes */
      for (i = 0; i < E->control.phasefile_C_num_nno; i++)
      {
        E->C_phasefile_nno[i] = (double *)malloc((nno + 1) * sizeof(double));
      }
    }

    if (E->control.phasevisc_C)
    {
      E->Cphasedot = (double *)malloc((nno + 1) * sizeof(double));
      E->Cphasedotnum = (double *)malloc((nno + 1) * sizeof(double));
      E->Cphase_node = (double *)malloc((nno + 1) * sizeof(double));

      E->Tphase_node = (double *)malloc((nno + 1) * sizeof(double));
      E->Pphase_node = (double *)malloc((nno + 1) * sizeof(double));

      for (i = 1; i <= nno; i++)
      {
        E->Cphasedot[i] = 0.0;
        E->Cphasedotnum[i] = 0.0;
      }
      //      if (!(strcmp(E->control.comp_adv_method,"field")==0)) {
      E->Cphase_marker = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
      E->Cphase_marker_old = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
      E->Tphase_marker = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
      E->Pphase_marker = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
      if (E->control.phasevisc_d)
      {
        E->d_marker = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
        E->d_marker_old = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
        E->d_dotnum = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
        E->d_dot = (double *)malloc((E->advection.markers_uplimit + 1) * sizeof(double));
        E->d_node = (double *)malloc((nno + 1) * sizeof(double));
      }
      //      }
    }
  }
  report(E, "convection, initial temperature");
  convection_initial_temperature(E);

  return;
}

/* =========================================== */

void convection_boundary_conditions(E) struct All_variables *E;

{
  void velocity_boundary_conditions();
  void temperature_boundary_conditions();
  void temperatures_conform_bcs();
  void composition_boundary_conditions();

  velocity_boundary_conditions(E); /* universal */
  fprintf(stderr, "okvelo \n");
  temperature_boundary_conditions(E);
  fprintf(stderr, "oktemp\n");
  temperatures_conform_bcs(E);
  fprintf(stderr, "ok10\n");
  composition_boundary_conditions(E);

  return;
}

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature(E) struct All_variables *E;
{
  int i, j, k, p, node, ii, jj;
  double temp, base, radius, radius2;
  double drand48();
  FILE *fp;
  void remove_horiz_ave();
  void temperatures_conform_bcs();
  void thermal_buoyancy();
  void process_restart_mk();
  void process_restart_tc();
  void convection_initial_markers();
  void convection_initial_markers1();
  void convection_initial_markers2();
  void convection_initial_markers_phasechange();
  int in1, in2, in3, instance, nfz, ok, noz2, ll, mm;
  char output_file[255];
  double lithThick;
  double tbase, tbase1, t1, r1, weight, para1, plate_velocity, delta_temp, age;
  double x00, x01, x02, slope, con;
  double age_x, tempt, temp1, temp2, temp3, Ti, inter_x, inter_z;
  const int dims = E->mesh.nsd;
  const double e_5 = 1.0e-5;
  double tempdist, center_x, center_z, center_r;
  const int noy = E->lmesh.noy;
  const int noz = E->lmesh.noz;
  const int nox = E->lmesh.nox;
  const int nno = nox * noz;

  para1 = E->control.Ts * E->data.ref_temperature + 0.4 * E->data.ref_temperature;

  tbase = (para1 - E->control.Ts * E->data.ref_temperature) / E->data.ref_temperature;
  tbase1 = (para1 + 200 - E->control.Ts * E->data.ref_temperature) / E->data.ref_temperature;

  mm = E->convection.perturb_k[0];

  /*
        noz2 = (noz-1)/2+1;
        con = (noz-1);
*/
  con = E->convection.perturb_mag[0];

  if (E->control.restart == 0)
  {

    for (i = 1; i <= noy; i++)
      for (j = 1; j <= nox; j++)
        for (k = 1; k <= noz; k++)
        {
          node = k + (j - 1) * noz + (i - 1) * nox * noz;
          t1 = E->X[1][node];
          r1 = E->X[2][node];

          E->T[node] = 0.0;
          E->C[node] = 0.0;
          E->T[node] = 1 - r1;
          Ti = E->control.inter_temp;
          /*                if (k==10) */

          switch (E->control.initialTOption)
          {
          case -1:
            E->T[node] = 0.0;
            break;
          case 0:
            if (E->control.temperature_perturbation)
            {
              E->T[node] += con * cos(M_PI * 2.0 * t1);
            }
            break;
          case 1:
            if (E->control.temperature_perturbation)
            {
              E->T[node] += con * cos(M_PI * 2.0 * t1);
            }
            if (E->control.continent_lith)
            {
              if (t1 >= E->control.continent_loc_left && t1 <= E->control.continent_loc_right)
              {
                E->T[node] = Ti;
                if (r1 >= 1.0 - E->control.depth_continent_lith)
                {
                  E->T[node] = Ti * (1.0 - r1) / E->control.depth_continent_lith;
                } /*end of depth*/
              }   /*end of loc */
            }     /*end of whether continent temperature */

            if (E->control.ocean_lith)
            {
              if (t1 >= E->control.age_loc_left && t1 <= E->control.age_loc_right)
              {
                E->T[node] = Ti;
                if (r1 >= 1.0 - E->control.depth_ocean_lith)
                {
                  age_x = E->control.age_left + (t1 - E->control.age_loc_left) / (E->control.age_loc_right - E->control.age_loc_left) * (E->control.age_right - E->control.age_left);
                  tempt = (1.0 - r1) * 0.5 / sqrt(age_x) * E->data.layer_km * 1.0e3 / sqrt(365 * 24 * 3600 * 1e6 * E->data.therm_diff);
                  E->T[node] = Ti * erf(tempt);
                  if (E->control.platemodel)
                  {
                    temp1 = (1.0 - r1) / E->control.depth_ocean_lith;
                    temp2 = 365 * 24 * 3600 * 1e6 * E->data.therm_diff * age_x;
                    temp3 = temp2 / (E->control.depth_ocean_lith * E->data.layer_km * 1.0e3 * E->control.depth_ocean_lith * E->data.layer_km * 1.0e3);
                    E->T[node] = Ti * (temp1 + 2.0 / M_PI * exp(-M_PI * M_PI * temp3) * sin(M_PI * temp1) + 1.0 / M_PI * exp(-4.0 * M_PI * M_PI * temp3) * sin(2.0 * M_PI * temp1) + 2.0 / 3.0 / M_PI * exp(-9.0 * M_PI * M_PI * temp3) * sin(3.0 * M_PI * temp1) + 2.0 / 4.0 / M_PI * exp(-16.0 * M_PI * M_PI * temp3) * sin(4.0 * M_PI * temp1));
                    if (E->T[node] > Ti)
                      E->T[node] = Ti;
                  } /*end of plate model */
                }   /*end of depth */
              }     /*end of age_loc*/
              if (E->control.ocean_lith_margin)
              {
                tempdist = (t1 - E->control.age_loc_left) * tan(E->control.dip_margin * 3.14159265 / 180.0) + 1.0 - r1;
                /*tempdist *= cos(E->control.dip_margin*3.14159265/180.0); */
                if (tempdist >= 0 && tempdist <= E->control.depth_ocean_lith && r1 >= 1.0 - E->control.depth_lith_margin && t1 < E->control.age_loc_left)
                {
                  age_x = E->control.age_left;
                  tempt = tempdist * 0.5 / sqrt(age_x) * E->data.layer_km * 1.0e3 / sqrt(365 * 24 * 3600 * 1e6 * E->data.therm_diff);
                  E->T[node] = Ti * erf(tempt);
                  if (E->control.platemodel)
                  {
                    temp1 = tempdist / E->control.depth_ocean_lith;
                    temp2 = 365 * 24 * 3600 * 1e6 * E->data.therm_diff * age_x;
                    temp3 = temp2 / (E->control.depth_ocean_lith * E->data.layer_km * 1.0e3 * E->control.depth_ocean_lith * E->data.layer_km * 1.0e3);
                    E->T[node] = Ti * (temp1 + 2.0 / M_PI * exp(-M_PI * M_PI * temp3) * sin(M_PI * temp1) + 1.0 / M_PI * exp(-4.0 * M_PI * M_PI * temp3) * sin(2.0 * M_PI * temp1) + 2.0 / 3.0 / M_PI * exp(-9.0 * M_PI * M_PI * temp3) * sin(3.0 * M_PI * temp1) + 2.0 / 4.0 / M_PI * exp(-16.0 * M_PI * M_PI * temp3) * sin(4.0 * M_PI * temp1));
                    if (E->T[node] > Ti)
                      E->T[node] = Ti;
                  } /*end of plate model */
                }
                else if (tempdist > E->control.depth_ocean_lith && r1 >= 1.0 - E->control.depth_lith_margin && t1 < E->control.age_loc_left)
                {
                  E->T[node] = Ti;
                }
              } /*end of oceanic plate margin*/
            }   /*end of whether ocean temperature */
            break;
          case 2:
            if (E->control.continent_lith)
            {
              if (t1 >= E->control.continent_loc_left && t1 <= E->control.continent_loc_right)
              {
                E->T[node] = Ti;
                if (r1 >= 1.0 - E->control.depth_continent_lith)
                {
                  E->T[node] = Ti * (1.0 - r1) / E->control.depth_continent_lith;
                } /*end of depth*/
              }   /*end of loc */
            }     /*end of whether continent temperature */
            if (E->control.ocean_lith)
            {
              if (t1 >= E->control.age_loc_left && t1 <= E->control.age_loc_right)
              {
                E->T[node] = Ti;
                if (r1 >= 1.0 - E->control.depth_ocean_lith)
                {
                  age_x = E->control.age_left + (t1 - E->control.age_loc_left) / (E->control.age_loc_right - E->control.age_loc_left) * (E->control.age_right - E->control.age_left);
                  tempt = (1.0 - r1) * 0.5 / sqrt(age_x) * E->data.layer_km * 1.0e3 / sqrt(365 * 24 * 3600 * 1e6 * E->data.therm_diff);
                  E->T[node] = Ti * erf(tempt);

                  if (E->control.platemodel)
                  {
                    temp1 = (1.0 - r1) / E->control.depth_ocean_lith;
                    temp2 = 365 * 24 * 3600 * 1e6 * E->data.therm_diff * age_x;
                    temp3 = temp2 / (E->control.depth_ocean_lith * E->data.layer_km * 1.0e3 * E->control.depth_ocean_lith * E->data.layer_km * 1.0e3);
                    E->T[node] = Ti * (temp1 + 2.0 / M_PI * exp(-M_PI * M_PI * temp3) * sin(M_PI * temp1) + 1.0 / M_PI * exp(-4.0 * M_PI * M_PI * temp3) * sin(2.0 * M_PI * temp1) + 2.0 / 3.0 / M_PI * exp(-9.0 * M_PI * M_PI * temp3) * sin(3.0 * M_PI * temp1) + 2.0 / 4.0 / M_PI * exp(-16.0 * M_PI * M_PI * temp3) * sin(4.0 * M_PI * temp1));
                    if (E->T[node] > Ti)
                      E->T[node] = Ti;
                  } /*end of plate model */
                }   /*end of depth */
              }     /*end of age_loc*/
              center_x = E->control.age_loc_left + E->control.dip_center_x;
              center_z = 1.0 - E->control.dip_center_z;
              center_r = sqrt(E->control.dip_center_x * E->control.dip_center_x + E->control.dip_center_z * E->control.dip_center_z);
              inter_z = 1.0 - E->control.depth_ocean_lith;
              inter_x = center_x - sqrt((center_r - E->control.depth_ocean_lith) * (center_r - E->control.depth_ocean_lith) - (center_z - inter_z) * (center_z - inter_z));
              if (r1 >= 1.0 - E->control.depth_lith_margin && r1 <= 1.0 - E->viscosity.zcrust1 && (1.0 - inter_z) * t1 + (inter_x - 1.0) * r1 <= inter_x - inter_z)
              {
                tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
                /*if (t1 >= 0.9 && r1 >= 0.9)
                  fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf\n", node, t1, t1 - center_x, r1, center_z, r1 - center_z, tempdist);*/
                if (tempdist <= E->control.depth_ocean_lith && tempdist >= 0)
                {
                  age_x = E->control.age_left;
                  tempt = tempdist * 0.5 / sqrt(age_x) * E->data.layer_km * 1.0e3 / sqrt(365 * 24 * 3600 * 1e6 * E->data.therm_diff);
                  E->T[node] = Ti * erf(tempt);
                  if (E->control.platemodel)
                  {
                    temp1 = tempdist / E->control.depth_ocean_lith;
                    temp2 = 365 * 24 * 3600 * 1e6 * E->data.therm_diff * age_x;
                    temp3 = temp2 / (E->control.depth_ocean_lith * E->data.layer_km * 1.0e3 * E->control.depth_ocean_lith * E->data.layer_km * 1.0e3);
                    E->T[node] = Ti * (temp1 + 2.0 / M_PI * exp(-M_PI * M_PI * temp3) * sin(M_PI * temp1) + 1.0 / M_PI * exp(-4.0 * M_PI * M_PI * temp3) * sin(2.0 * M_PI * temp1) + 2.0 / 3.0 / M_PI * exp(-9.0 * M_PI * M_PI * temp3) * sin(3.0 * M_PI * temp1) + 2.0 / 4.0 / M_PI * exp(-16.0 * M_PI * M_PI * temp3) * sin(4.0 * M_PI * temp1));
                    if (E->T[node] > Ti)
                      E->T[node] = Ti;
                  } /*end of plate model */
                }
              } /*end of oceanic plate margin*/
            }   /*end of whether ocean temperature */

            break;
          case 10:
            if (E->control.continent_lith)
            {
              if (t1 >= E->control.continent_loc_left && t1 <= E->control.continent_loc_right)
              {
                E->T[node] = Ti;
                if (r1 >= 1.0 - E->control.depth_continent_lith)
                {
                  age_x = E->control.continent_age;

                  tempt = (1.0 - r1) * 0.5 / sqrt(age_x) * E->data.layer_km * 1.0e3 / sqrt(365 * 24 * 3600 * 1e6 * E->data.therm_diff);
                  E->T[node] = Ti * erf(tempt);
                  if (E->control.continent_platemodel)
                  {
                    temp1 = (1.0 - r1) / E->control.depth_continent_lith;
                    temp2 = 365 * 24 * 3600 * 1e6 * E->data.therm_diff * age_x;
                    temp3 = temp2 / (E->control.depth_continent_lith * E->data.layer_km * 1.0e3 * E->control.depth_continent_lith * E->data.layer_km * 1.0e3);
                    E->T[node] = Ti * (temp1 + 2.0 / M_PI * exp(-M_PI * M_PI * temp3) * sin(M_PI * temp1) + 1.0 / M_PI * exp(-4.0 * M_PI * M_PI * temp3) * sin(2.0 * M_PI * temp1) + 2.0 / 3.0 / M_PI * exp(-9.0 * M_PI * M_PI * temp3) * sin(3.0 * M_PI * temp1) + 2.0 / 4.0 / M_PI * exp(-16.0 * M_PI * M_PI * temp3) * sin(4.0 * M_PI * temp1));
                    if (E->T[node] > Ti)
                      E->T[node] = Ti;
                  } /*end of plate model */
                }   /*end of depth*/
              }     /*end of loc */
            }       /*end of whether continent temperature */
            if (E->control.ocean_lith)
            {
              if (t1 >= E->control.age_loc_left && t1 <= E->control.age_loc_right)
              {
                E->T[node] = Ti;
                if (r1 >= 1.0 - E->control.depth_ocean_lith)
                {
                  age_x = E->control.age_left + (t1 - E->control.age_loc_left) / (E->control.age_loc_right - E->control.age_loc_left) * (E->control.age_right - E->control.age_left);
                  tempt = (1.0 - r1) * 0.5 / sqrt(age_x) * E->data.layer_km * 1.0e3 / sqrt(365 * 24 * 3600 * 1e6 * E->data.therm_diff);
                  E->T[node] = Ti * erf(tempt);

                  if (E->control.platemodel)
                  {
                    temp1 = (1.0 - r1) / E->control.depth_ocean_lith;
                    temp2 = 365 * 24 * 3600 * 1e6 * E->data.therm_diff * age_x;
                    temp3 = temp2 / (E->control.depth_ocean_lith * E->data.layer_km * 1.0e3 * E->control.depth_ocean_lith * E->data.layer_km * 1.0e3);
                    E->T[node] = Ti * (temp1 + 2.0 / M_PI * exp(-M_PI * M_PI * temp3) * sin(M_PI * temp1) + 1.0 / M_PI * exp(-4.0 * M_PI * M_PI * temp3) * sin(2.0 * M_PI * temp1) + 2.0 / 3.0 / M_PI * exp(-9.0 * M_PI * M_PI * temp3) * sin(3.0 * M_PI * temp1) + 2.0 / 4.0 / M_PI * exp(-16.0 * M_PI * M_PI * temp3) * sin(4.0 * M_PI * temp1));
                    if (E->T[node] > Ti)
                      E->T[node] = Ti;
                  } /*end of plate model */
                }   /*end of depth */
              }     /*end of age_loc*/
              center_x = E->control.age_loc_left + E->control.dip_center_x;
              center_z = 1.0 - E->control.dip_center_z;
              center_r = sqrt(E->control.dip_center_x * E->control.dip_center_x + E->control.dip_center_z * E->control.dip_center_z);
              inter_z = 1.0 - E->control.depth_ocean_lith;
              inter_x = center_x - sqrt((center_r - E->control.depth_ocean_lith) * (center_r - E->control.depth_ocean_lith) - (center_z - inter_z) * (center_z - inter_z));
              if (r1 >= 1.0 - E->control.depth_lith_margin && r1 <= 1.0 - E->viscosity.zcrust1 && (1.0 - inter_z) * t1 + (inter_x - 1.0) * r1 <= inter_x - inter_z)
              {
                tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));

                /*if (t1 >= 0.9 && r1 >= 0.9)
                  fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf\n", node, t1, t1 - center_x, r1, center_z, r1 - center_z, tempdist);*/
                if (tempdist <= E->control.depth_ocean_lith && tempdist >= 0)
                {
                  age_x = E->control.age_left;
                  tempt = tempdist * 0.5 / sqrt(age_x) * E->data.layer_km * 1.0e3 / sqrt(365 * 24 * 3600 * 1e6 * E->data.therm_diff);
                  E->T[node] = Ti * erf(tempt);
                  if (E->control.platemodel)
                  {
                    temp1 = tempdist / E->control.depth_ocean_lith;
                    temp2 = 365 * 24 * 3600 * 1e6 * E->data.therm_diff * age_x;
                    temp3 = temp2 / (E->control.depth_ocean_lith * E->data.layer_km * 1.0e3 * E->control.depth_ocean_lith * E->data.layer_km * 1.0e3);
                    E->T[node] = Ti * (temp1 + 2.0 / M_PI * exp(-M_PI * M_PI * temp3) * sin(M_PI * temp1) + 1.0 / M_PI * exp(-4.0 * M_PI * M_PI * temp3) * sin(2.0 * M_PI * temp1) + 2.0 / 3.0 / M_PI * exp(-9.0 * M_PI * M_PI * temp3) * sin(3.0 * M_PI * temp1) + 2.0 / 4.0 / M_PI * exp(-16.0 * M_PI * M_PI * temp3) * sin(4.0 * M_PI * temp1));
                    if (E->T[node] > Ti)
                      E->T[node] = Ti;
                  } /*end of plate model */
                }
              } /*end of oceanic plate margin*/
            }   /*end of whether ocean temperature */
            break;

          } /* end of switch*/

          E->C[node] = 0.0;

          E->node[node] = E->node[node] | (INTX | INTZ | INTY);

        } /* close the loop for node */

    if (E->control.composition)
    {
      if (!(strcmp(E->control.comp_adv_method, "field") == 0))
      {
        convection_initial_markers(E);
      }
      else if ((strcmp(E->control.comp_adv_method, "field") == 0))
      {
        for (node = 1; node <= nno; node++)
        {
          t1 = E->X[1][node];
          r1 = E->X[2][node];
          switch (E->control.initialCOption)
          {
          case 1:
            if (r1 <= 1.0 - E->viscosity.zcrust1)
              E->C[node] = 1.0;
            else
              E->C[node] = 0.0;
            break;
          case 2:
            if (r1 >= 1.0 - E->viscosity.zcrust1)
              E->C[node] = 1.0;
            else
              E->C[node] = 0.0;
          default:
            break;
          }
        }
      }
    }
  }
  // end for restart=0

  else if (E->control.restart == 1)
  {
    process_restart_tc(E, E->mesh.levmax);

    if (E->control.composition && !(strcmp(E->control.comp_adv_method, "field") == 0))
      process_restart_mk(E);
  }

  temperatures_conform_bcs(E);

  thermal_buoyancy(E);

  return;
}

// Wei 2019 June 16
void convection_initial_markers(E) struct All_variables *E;
{
  int el, i, j, k, p, node, ii, jj;
  double dx, dX[4], dr, x, z;
  char input_s[100], output_file[255];
  FILE *fp;
  void get_C_from_markers();
  void get_C_from_markers_multi();
  void get_C_from_markers_double();
  double t1, r1, tempdist, loc_mid, center_x, center_z;
  node = 0;
  const int nel = E->lmesh.nel;
  const int nno = E->lmesh.nno;
  p = pow((double)E->advection.markers_per_ele, (double)(1.0 / E->mesh.dof));
  fprintf(stderr, "p = %d el = %d", p, nel);
  do
  {
    x = drand48() * (E->XG2[1] - E->XG1[1]);
    z = drand48() * (E->XG2[2] - E->XG1[2]);
    if ((x >= E->XP[1][1] && x <= E->XP[1][E->lmesh.nox]) && (z >= E->XP[2][1] && z <= E->XP[2][E->lmesh.noz]))
    {
      node++;
      E->XMC[1][node] = x;
      E->XMC[2][node] = z;
      el = get_element(E, E->XMC[1][node], E->XMC[2][node], dX);
      E->CElement[node] = el;
      E->C12[node] = 0;
      switch (E->control.initialCOption)
      {
      case 1:
        if (E->XMC[2][node] > (E->viscosity.zcrust1 + 0.02 * cos(M_PI * E->XMC[1][node] / E->mesh.layer[1])))
          E->C12[node] = 0;
        else
          E->C12[node] = 1;
        break;
      case 2:
        break;
      default:
        break;
      }

      if (E->control.phasefile_C || E->control.phasefile_Complete)
      {
        if (E->XMC[2][node] < 1.0 - E->viscosity.zcrust1)
          E->C12[node] = 0;
        else
          E->C12[node] = 1;
        E->C_phasefile_marker_int[0][node] = E->C12[node];
        t1 = E->XMC[1][node];
        r1 = E->XMC[2][node];
      }

      if (E->control.phasefile_C)
      {
        if (E->XMC[2][node] < 1.0 - E->viscosity.zcrust1 && E->XMC[2][node] > 1.0 - E->control.depth_harz)
        {
          E->C_phasefile_marker_int[0][node] = 2;
        }

        if (E->control.ocean_lith_margin)
        {
          switch (E->control.initialCOption)
          {
          case 1:
            tempdist = (t1 - E->control.age_loc_left) * tan(E->control.dip_margin * 3.14159265 / 180.0) + 1.0 - r1;
            tempdist *= cos(E->control.dip_margin * 3.14159265 / 180.0);
            break;
          case 2:
            center_x = E->control.age_loc_left + E->control.dip_center_x;
            center_z = 1.0 - E->control.dip_center_z;
            tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
            break;
          }

          if (tempdist >= 0 && tempdist <= E->viscosity.zcrust1 && r1 >= 1.0 - E->control.depth_lith_margin && t1 < E->control.age_loc_left)
          {
            E->C12[node] = 1;
            if (E->control.phasefile_C)
            {
              E->C_phasefile_marker_int[0][node] = 1;
            }
          }
          if (E->control.phasefile_C)
          {
            if (tempdist > E->viscosity.zcrust1 && tempdist < E->control.depth_harz && r1 >= 1.0 - E->control.depth_lith_margin && r1 <= 1.0 - E->viscosity.zcrust1 && t1 < E->control.age_loc_left)
            {
              E->C_phasefile_marker_int[0][node] = 2;
            }
          }
        }
      }

      if (E->control.phasefile_Complete)
      {
        if (r1 < E->viscosity.zlm)
        {
          E->C_phasefile_marker_int[0][node] = E->control.phasefile_C_flavor - 1;
        }
        if (t1 >= E->control.age_loc_left && t1 <= E->control.age_loc_right)
        {
          if (r1 < 1.0 - E->viscosity.zcrust1 && r1 > 1.0 - E->control.depth_harz)
          {
            E->C_phasefile_marker_int[0][node] = 2;
          }
          else if (r1 <= 1.0 - E->control.depth_harz && r1 > 1.0 - E->control.depth_ocean_lith)
          {
            E->C_phasefile_marker_int[0][node] = 3;
          }
        }
        switch (E->control.initialCOption)
        {
        case 1:
          tempdist = (t1 - E->control.age_loc_left) * tan(E->control.dip_margin * 3.14159265 / 180.0) + 1.0 - r1;
          //tempdist *= cos(E->control.dip_margin * 3.14159265 / 180.0);
          if (E->control.continent_lith)
          {
            if (t1 >= E->control.continent_loc_left && t1 <= E->control.continent_loc_right)
            {
              if (r1 >= 1.0 - E->control.depth_continent_crust)
              {
                E->C_phasefile_marker_int[0][node] = 4;
              }
              else if (r1 >= 1.0 - E->control.depth_continent_lith)
              {
                E->C_phasefile_marker_int[0][node] = 5;
              }
            }
          }
          break;
        case 2:
          center_x = E->control.age_loc_left + E->control.dip_center_x;
          center_z = 1.0 - E->control.dip_center_z;
          tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
          if (E->control.continent_lith)
          {
            if (t1 >= E->control.continent_loc_left && t1 <= E->control.continent_loc_right)
            {
              if (r1 >= 1.0 - E->control.depth_continent_crust)
              {
                E->C_phasefile_marker_int[0][node] = 4;
              }
              else if (r1 >= 1.0 - E->control.depth_continent_lith)
              {
                E->C_phasefile_marker_int[0][node] = 5;
              }
            }
          }
          break;
        case 10:
          center_x = E->control.age_loc_left + E->control.dip_center_x;
          center_z = 1.0 - E->control.dip_center_z;
          tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
          if (E->control.continent_lith)
          {
            if (t1 >= E->control.continent_loc_left && t1 <= E->control.continent_loc_right)
            {
              if (r1 >= 1.0 - E->viscosity.zcrust1)
              {
                E->C_phasefile_marker_int[0][node] = 4;
              }
              else if (r1 >= 1.0 - E->control.depth_harz)
              {
                E->C_phasefile_marker_int[0][node] = 5;
              }
              else if (r1 >= 1.0 - E->control.depth_continent_lith)
              {
                E->C_phasefile_marker_int[0][node] = 6;
              }
              else if (r1 >= E->viscosity.zlm)
              {
                E->C_phasefile_marker_int[0][node] = 0;
              }
            }
          }
          break;
        case 50:
          center_x = E->control.age_loc_left + E->control.dip_center_x;
          center_z = 1.0 - E->control.dip_center_z;
          tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
          if (r1 < E->viscosity.zlm)
          {
            E->C_phasefile_marker_int[0][node] = 0;
          }
          if (E->control.continent_lith)
          {
            if (t1 >= E->control.continent_loc_left && t1 <= E->control.continent_loc_right)
            {
              if (r1 >= 1.0 - E->viscosity.zcrust1)
              {
                E->C_phasefile_marker_int[0][node] = 4;
              }
              else if (r1 >= 1.0 - E->control.depth_harz)
              {
                E->C_phasefile_marker_int[0][node] = 5;
              }
              else if (r1 >= 1.0 - E->control.depth_continent_lith)
              {
                E->C_phasefile_marker_int[0][node] = 6;
              }
              else if (r1 >= E->viscosity.zlm)
              {
                E->C_phasefile_marker_int[0][node] = 0;
              }
            }
          }
          break;
        } // end of initialCOption

        if (r1 >= 1.0 - E->control.depth_lith_margin)
        {
          if (E->control.initialCOption == 1)
          {
            center_x = E->control.age_loc_left;
          }
          if (t1 < center_x)
          {
            if (tempdist >= 0 && tempdist <= E->viscosity.zcrust1)
            {
              E->C12[node] = 1;
              E->C_phasefile_marker_int[0][node] = 1;
            }
            else if (tempdist > E->viscosity.zcrust1 && tempdist <= E->control.depth_harz && r1 < 1.0 - E->viscosity.zcrust1)
            {
              E->C_phasefile_marker_int[0][node] = 2;
            }
            else if (tempdist > E->control.depth_harz && tempdist <= E->control.depth_ocean_lith && r1 < 1.0 - E->control.depth_harz)
            {
              E->C_phasefile_marker_int[0][node] = 3;
            }
          }
        }
      } // end of phase file Complete
      // here to store initial marker type
      E->C_phasefile_marker_int[1][node] = E->C_phasefile_marker_int[0][node];
      if (E->control.phasevisc_d)
      {
        E->C_phasefile_marker_double[E->C_phasefile_d_start_double][node] = E->control.phasevisc_d0;
        E->C_phasefile_marker_double[E->C_phasefile_d_start_double + 1][node] = E->control.phasevisc_d0;
        E->C_phasefile_marker_double[E->C_phasefile_d_start_double + 2][node] = 0.0;
        E->C_phasefile_marker_double[E->C_phasefile_d_start_double + 3][node] = 0.0;
      }
      E->C_phasefile_marker_int[2][node] = E->parallel.me;
      E->C_phasefile_marker_int[3][node] = node;
    } // end of coord confirm
  } while (node < E->advection.markers);

  E->advection.markers = node;
  //  fprintf(stderr, "before Cfrom markers\n");
  get_C_from_markers(E, E->C);
  //  fprintf(stderr, "finish Cfrom markers\n");

  if (E->control.phasefile_C || E->control.phasefile_Complete)
  {
    get_C_from_markers_multi(E, E->C_phasefile_marker_int[0], E->C_phasefile_nno, E->C_phasefile_element, 0, E->control.phasefile_C_flavor, E->CElement);
    if (E->control.phasevisc_d)
      get_C_from_markers_double(E, E->C_phasefile_marker_double[E->C_phasefile_d_start_double], E->C_phasefile_nno[E->C_phasefile_d_start_nno], E->C_phasefile_element[E->C_phasefile_d_start_element], E->CElement, 1);
    /* 1 store last step 0 store current step */
    for (i = 1; i <= E->advection.markers; i++)
    {
      /* give place for old markers */
      E->C_phasefile_marker_int[1][i] = E->C_phasefile_marker_int[0][i];
    }
  }

  return;
}

/* ====================================================================== */
// for both isochemical and thermochemical convection with field methods

void process_restart_tc(E, lev) struct All_variables *E;
int lev;
{
  int fileid[20];
  int i, j, k, ii, size2;
  char input_s[200], output_file[255], in_file[255];
  FILE *fp;
  double t1, t2, t3, t4, t5, t6;

  sprintf(output_file, "%s/temp_comp.%d.%d", E->convection.old_T_file, E->parallel.me, E->monitor.solution_cycles);

  fp = fopen(output_file, "r");

  fgets(input_s, 200, fp);
  sscanf(input_s, "%d %d %lf", &i, &E->advection.timesteps, &E->monitor.elapsed_time);
  if (E->control.composition)
  {
    for (i = 1; i <= E->lmesh.NNO[lev]; i++)
    {
      fgets(input_s, 200, fp);
      if (E->control.phasevisc_C)
      {
        if (E->control.phasevisc_d)
        {
          sscanf(input_s, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &E->T[i], &E->C[i], &E->V[1][i], &E->V[2][i], &E->Vi[i], &E->Cdot[i], &E->Cphasedotnum[i], &E->Cphasedot[i], &E->Cphase_node[i], &E->d_node[i]);
        }
      }
      else if (E->control.phasefile_C || E->control.phasefile_Complete)
      {
        if (E->control.phasevisc_d)
        {
          sscanf(input_s, "%lf %lf %lf %lf %lf %lf %lf %lf", &E->T[i], &E->C[i], &t1, &t2, &E->Vi[i], &t3, &t4, &t5, &E->C_phasefile_nno[E->C_phasefile_d_start_nno][i], &E->C_phasefile_nno[E->C_phasefile_d_start_nno + 1][i], &E->C_phasefile_nno[E->C_phasefile_d_start_nno + 2][i], &E->C_phasefile_nno[E->C_phasefile_d_start_nno + 3][i]);
        }
        else
          sscanf(input_s, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &E->T[i], &E->C[i], &t1, &t2, &E->Vi[i], &E->C_phasefile_nno[0][i], &E->C_phasefile_nno[1][i], &E->C_phasefile_nno[2][i], &E->C_phasefile_nno[3][i]);
      }
      else
      {
        sscanf(input_s, "%lf %lf %lf %lf %lf", &E->T[i], &E->C[i], &t1, &t2, &t3);
      }
      E->U[E->id[i].doff[1]] = t1;
      E->U[E->id[i].doff[2]] = t2;
    }
    for (i = 1; i <= E->lmesh.NEL[lev]; i++)
    {
      fgets(input_s, 200, fp);
      sscanf(input_s, "%lf", &t1);
      E->P[i] = t1;
    }
  }
  else
  {
    for (i = 1; i <= E->lmesh.NNO[lev]; i++)
    {
      fgets(input_s, 200, fp);
      sscanf(input_s, "%lf %lf %lf", &E->T[i], &t1, &t2);
      E->U[E->id[i].doff[1]] = t1;
      E->U[E->id[i].doff[2]] = t2;
      E->C[i] = 0;
    }
    for (i = 1; i <= E->lmesh.NEL[lev]; i++)
    {
      fgets(input_s, 200, fp);
      sscanf(input_s, "%lf", &t1);
      E->P[i] = t1;
    }
  }

  E->advection.timesteps = E->monitor.solution_cycles;

  fclose(fp);

  return;
}

/* ====================================================================== */

void process_restart_mk(E) struct All_variables *E;
{
  int fileid[20];
  int i, j, k, ii, size2;
  char input_s[200], output_file[255], in_file[255];
  FILE *fp;
  int temp_int_1, temp_int_2;
  double t1, temp_double_1;
  const int nel = E->lmesh.nel;
  void get_C_from_markers();
  void get_C_from_markers_multi();
  void get_C_from_markers_double();
  sprintf(output_file, "%s/traces.%d.%d", E->convection.old_T_file, E->parallel.me, E->monitor.solution_cycles);

  fp = fopen(output_file, "r");

  fgets(input_s, 200, fp);
  sscanf(input_s, "%d %d %d %lf", &E->advection.markers, &temp_int_1, &temp_int_2, &temp_double_1);
  for (i = 1; i <= E->advection.markers; i++)
  {
    fgets(input_s, 200, fp);
    if (E->control.phasevisc_C)
    {
      if (E->control.phasevisc_d)
      {
        sscanf(input_s, "%lf %lf %d %d %lf %lf %lf %lf %lf %lf\n", &E->XMC[1][i], &E->XMC[2][i], &E->CElement[i], &E->C12[i], &E->Cphase_marker[i], &E->Cphase_marker_old[i], &E->d_marker[i], &E->d_marker_old[i], &E->d_dotnum[i], &E->d_dot[i]);
      }
    }
    else if (E->control.phasefile_C || E->control.phasefile_Complete)
    {
      if (E->control.phasevisc_d)
      {
        sscanf(input_s, "%lf %lf %d %d %d %d %d %d %lf %lf %lf %lf\n", &E->XMC[1][i], &E->XMC[2][i], &E->CElement[i], &E->C12[i], &E->C_phasefile_marker_int[0][i], &E->C_phasefile_marker_int[1][i], &E->C_phasefile_marker_int[2][i], &E->C_phasefile_marker_int[3][i], &E->C_phasefile_marker_double[E->C_phasefile_d_start_double][i], &E->C_phasefile_marker_double[E->C_phasefile_d_start_double + 1][i], &E->C_phasefile_marker_double[E->C_phasefile_d_start_double + 2][i], &E->C_phasefile_marker_double[E->C_phasefile_d_start_double + 3][i]);
      }
      else
        sscanf(input_s, "%lf %lf %d %d %d %d %d %d", &E->XMC[1][i], &E->XMC[2][i], &E->CElement[i], &E->C12[i], &E->C_phasefile_marker_int[0][i], &E->C_phasefile_marker_int[1][i], &E->C_phasefile_marker_int[2][i], &E->C_phasefile_marker_int[3][i]);
    }
    else
      sscanf(input_s, "%lf %lf %d %d", &E->XMC[1][i], &E->XMC[2][i], &E->CElement[i], &E->C12[i]);
    if (E->XMC[2][i] < E->XP[2][1])
      E->XMC[2][i] = E->XP[2][1];
    else if (E->XMC[2][i] > E->XP[2][E->lmesh.noz])
      E->XMC[2][i] = E->XP[2][E->lmesh.noz];
    if (E->XMC[1][i] < E->XP[1][1])
      E->XMC[1][i] = E->XP[1][1];
    else if (E->XMC[1][i] > E->XP[1][E->lmesh.nox])
      E->XMC[1][i] = E->XP[1][E->lmesh.nox];
  }
  for (i = 1; i <= nel; i++)
  {
    fgets(input_s, 200, fp);
    if (E->control.phasevisc_C)
    {
      sscanf(input_s, "%lf %lf\n", E->CE[i], E->CphaseE[i]);
    }
    else
      sscanf(input_s, "%lf", &E->CE[i]);
  }
  fclose(fp);

  E->advection.timesteps = E->monitor.solution_cycles;

  get_C_from_markers(E, E->C);
  if (E->control.phasefile_C || E->control.phasefile_Complete)
  {
    get_C_from_markers_multi(E, E->C_phasefile_marker_int[0], E->C_phasefile_nno, E->C_phasefile_element, 0, E->control.phasefile_C_flavor, E->CElement);
    if (E->control.phasevisc_d)
    {

      get_C_from_markers_double(E, E->C_phasefile_marker_double[E->C_phasefile_d_start_double], E->C_phasefile_nno[E->C_phasefile_d_start_nno], E->C_phasefile_element[E->C_phasefile_d_start_element], E->CElement, 1);
    }
  }

  return;
}
