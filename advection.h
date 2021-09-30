struct ADVECTION
{
  int ADVECTION;

  double gamma;
  double timestep;
  double fine_tune_dt;
  double dt_reduced;
  double fixed_timestep;
  double max_dimensionless_time;

  int markers_g, markers, markers1, markers_uplimit, markers_per_ele;
  double marker_maxdist, marker_mindist;
  int markerIX, markerIZ;

  int min_timesteps;
  int max_timesteps;
  int max_total_timesteps;
  int timesteps;
  int total_timesteps;
  int temp_iterations;
  int max_substeps;
  int sub_iterations;
  int last_sub_iterations;

  double vel_substep_aggression;
  double temp_updatedness;
  double visc_updatedness;

  double lid_defining_velocity;
  double sub_layer_sample_level;

} advection;
