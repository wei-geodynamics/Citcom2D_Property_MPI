/*****************************************
 *   CC  III  TTTTT   CC   OO   MM MM    *
 *  C     I     T    C    O  O  M M M    *
 *  C     I     T    C    O  O  M   M    *
 *   CC  III    T     CC   OO   M   M    *
 *                                       *  
 * Developed at CIT for COnvection in    *
 * the Mantle by Louis Moresi 1992-today *
 *                                       *
 * You are free to use this code but it  * 
 * is distrubuted as BeWare i.e. it does *
 * not carry any guarantees or warranties *
 * of reliability.                       *
 *                                       *
 * Please respect all the time and work  *
 * that went into the development of the *
 * code.                                 *  
 *                                       *
 * LM                                    *
 *****************************************/
/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void get_viscosity_option(E) struct All_variables *E;
{
    void viscosity_for_system();

    /* general, essential default */

    E->viscosity.update_allowed = 0;
    E->viscosity.SDEPV = E->viscosity.TDEPV = E->viscosity.CHEMDEPV = 0;
    E->viscosity.EXPX = 0;

    input_string("Viscosity", E->viscosity.STRUCTURE, NULL); /* Which form of viscosity */

    input_boolean("VISC_EQUIVDD", &(E->viscosity.EQUIVDD), "off"); /* Whether to average it */
    input_int("equivdd_opt", &(E->viscosity.equivddopt), "1");
    input_int("equivdd_x", &(E->viscosity.proflocx), "1");
    input_int("equivdd_y", &(E->viscosity.proflocy), "1");

    input_boolean("VISC_SMOOTH", &(E->viscosity.SMOOTH), "off");
    input_int("visc_smooth_cycles", &(E->viscosity.smooth_cycles), "0");

    if (strcmp(E->viscosity.STRUCTURE, "system") == 0) /* Interpret */
    {
        fprintf(E->fp, "Viscosity derived from system state\n");
        E->viscosity.FROM_SYSTEM = 1;
        viscosity_for_system(E);
    }

    return;
}

/* ============================================ */

void viscosity_for_system(E) struct All_variables *E;
{
    void get_system_viscosity();
    void twiddle_thumbs();
    int l, i;

    /* default values .... */

    for (i = 0; i < 40; i++)
    {
        E->viscosity.N0[i] = 1.0;
        E->viscosity.T[i] = 0.0;
        E->viscosity.Z[i] = 0.0;
        E->viscosity.E[i] = 0.0;
        E->viscosity.T0[i] = 0.0;
        E->viscosity.sdepv_expt[i] = 0.0;
        E->viscosity.sdepv_trns[i] = 0.0;
    }

    /* read in information */
    input_int("rheol", &(E->viscosity.RHEOL), "essential");
    input_int("num_mat", &(E->viscosity.num_mat), "1");

    input_double_vector("viscT", E->viscosity.num_mat, (E->viscosity.T)); /* redundant */
    input_double_vector("viscT1", E->viscosity.num_mat, (E->viscosity.T));
    input_double_vector("viscZ", E->viscosity.num_mat, (E->viscosity.Z));
    input_double_vector("viscE", E->viscosity.num_mat, (E->viscosity.E));
    input_double_vector("viscT0", E->viscosity.num_mat, (E->viscosity.T0));
    input_double_vector("visc0", E->viscosity.num_mat, (E->viscosity.N0)); /* redundant */
    input_double_vector("viscsdepv_expt", E->viscosity.num_mat, (E->viscosity.sdepv_expt));
    input_double_vector("viscsdepv_trns", E->viscosity.num_mat, (E->viscosity.sdepv_trns));
    input_boolean("CHEMDEPV", &(E->viscosity.CHEMDEPV), "off");
    input_boolean("TDEPV", &(E->viscosity.TDEPV), "on");
    input_boolean("SDEPV", &(E->viscosity.SDEPV), "off");
    input_boolean("PDEPV", &(E->viscosity.PDEPV), "off");
    input_boolean("lower_diff", &(E->viscosity.lower_diff), "off");
    input_double("sdepv_misfit", &(E->viscosity.sdepv_misfit), "0.001");
    input_double_vector("sdepv_expt", E->viscosity.num_mat, (E->viscosity.sdepv_expt));
    input_double_vector("sdepv_trns", E->viscosity.num_mat, (E->viscosity.sdepv_trns));

    input_double("sdepv_iter_damp", &(E->viscosity.sdepv_iter_damp), "1.0");

    input_boolean("BDEPV", &(E->viscosity.BDEPV), "off");

    input_double_vector("abyerlee", E->viscosity.num_mat,
                        (E->viscosity.abyerlee));
    input_double_vector("bbyerlee", E->viscosity.num_mat,
                        (E->viscosity.bbyerlee));
    input_double_vector("lbyerlee", E->viscosity.num_mat,
                        (E->viscosity.lbyerlee));

    /* 1: transition 0: min/max transition */
    input_boolean("plasticity_trans", &(E->viscosity.plasticity_trans), "on");
    /* SRW */
    input_boolean("psrw", &(E->viscosity.psrw), "off");
    input_boolean("TDEPV_AVE", &(E->viscosity.TDEPV_AVE), "off");
    input_boolean("VFREEZE", &(E->viscosity.FREEZE), "off");
    input_boolean("VMAX", &(E->viscosity.MAX), "off");
    input_boolean("VMIN", &(E->viscosity.MIN), "off");
    input_boolean("VISC_UPDATE", &(E->viscosity.update_allowed), "on");

    input_double("freeze_thresh", &(E->viscosity.freeze_thresh), "0.0");
    input_double("freeze_value", &(E->viscosity.freeze_value), "1.0");
    input_double("visc_max", &(E->viscosity.max_value), "nodefault");
    input_double("visc_min", &(E->viscosity.min_value), "nodefault");
    input_double("weakzone_ratio", &(E->viscosity.weakzone_ratio), "1.0");
    input_boolean("VISC_GUESS", &(E->viscosity.guess), "off");
    input_string("visc_old_file", E->viscosity.old_file, " ");
    input_boolean("visc_platebond", &(E->viscosity.visc_platebond), "off");
    input_boolean("visc_platebond_const", &(E->viscosity.visc_platebond_const), "off");
    input_boolean("visc_platebond_selfadapt", &(E->viscosity.visc_platebond_selfadapt), "off");
    input_double("z_weakzone_platebond", &(E->viscosity.z_weakzone_platebond), "0.95644599");
    input_double("width_weakzone_platebond", &(E->viscosity.width_weakzone_platebond), "0.0348432");
    input_double("visc_reduce_platebond", &(E->viscosity.visc_reduce_platebond), "0.01");
    input_double("left_weakzone_platebond", &(E->viscosity.left_weakzone_platebond), "0.0348432");
    input_double("right_weakzone_platebond", &(E->viscosity.right_weakzone_platebond), "0.0");
    input_double("left_weakzone_platebond_ratio", &(E->viscosity.left_weakzone_platebond_ratio), "6.0");
    /*
    for (l=0;l<E->viscosity.num_mat;l++)  {
      E->viscosity.Z[l] = E->data.density*E->data.grav_acc*E->data.layer_meter*E->viscosity.Z[l]/(E->data.gas_const*E->data.ref_temperature);
      E->viscosity.E[l] = E->viscosity.E[l]/(E->data.gas_const*E->data.ref_temperature);
     fprintf(E->fp,"E & Z %lf %lf\n",E->viscosity.E[l],E->viscosity.Z[l]);
     }
*/

    if (!E->viscosity.update_allowed)
    {
        get_system_viscosity(E, 1, E->EVI[E->mesh.levmax], E->VI[E->mesh.levmax]);
    }

    return;
}

void get_system_viscosity(E, propogate, evisc, visc) struct All_variables *E;
int propogate;
double *visc, *evisc;
{
    void visc_from_mat();
    void visc_from_T();
    void visc_from_S();
    void visc_from_P();
    void apply_viscosity_smoother();
    void v_to_nodes();
    void visc_from_gint_to_nodes();
    int i, j;
    double *evisc_old0, *visc_old0, *evisc_old1, *visc_old1;

    const int vpts = vpoints[E->mesh.nsd];

    if (E->viscosity.TDEPV)
        visc_from_T(E, visc, evisc, propogate);
    else
        visc_from_mat(E, visc, evisc);

    if (E->viscosity.SDEPV)
        visc_from_S(E, visc, evisc, propogate);
    if (E->viscosity.PDEPV)
        visc_from_P(E, visc, evisc, propogate);
    if (E->viscosity.SMOOTH)
        apply_viscosity_smoother(E, visc, evisc);

    if (E->viscosity.MAX)
    {
        for (i = 1; i <= E->lmesh.nel; i++)
            for (j = 1; j <= vpts; j++)
            {
                if (evisc[(i - 1) * vpts + j] > E->viscosity.max_value)
                    evisc[(i - 1) * vpts + j] = E->viscosity.max_value;
            }
    }

    if (E->viscosity.MIN)
    {
        for (i = 1; i <= E->lmesh.nel; i++)
            for (j = 1; j <= vpts; j++)
                if (evisc[(i - 1) * vpts + j] < E->viscosity.min_value)
                    evisc[(i - 1) * vpts + j] = E->viscosity.min_value;
    }

    /*v_to_nodes(E, evisc, visc, E->mesh.levmax);*/
    visc_from_gint_to_nodes(E, evisc, visc, E->mesh.levmax);

    return;
}

void apply_viscosity_smoother(E, visc, evisc) struct All_variables *E;
double *visc, *evisc;

{
    void p_to_centres();
    void p_to_nodes();
    void element_markers();
    void get_C_from_markers_multi();

    double *ViscCentre;
    int i;

    ViscCentre = (double *)malloc((E->lmesh.nno + 10) * sizeof(double));

    for (i = 1; i <= E->viscosity.smooth_cycles; i++)
    {
        p_to_centres(E, visc, ViscCentre, E->mesh.levmax);
        p_to_nodes(E, ViscCentre, visc, E->mesh.levmax);
    }

    free((void *)ViscCentre);

    return;
}

void visc_from_mat(E, Eta, EEta) struct All_variables *E;
double *Eta, *EEta;
{

    int i, j, k, l, z, jj, kk;
    double temp1, temp[9], C1, C2, visc1, pres0, temp2, rii[9];
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    for (i = 1; i <= E->lmesh.nel; i++)
    {
        for (jj = 1; jj <= vpoints[E->mesh.nsd]; jj++)
        {
            EEta[(i - 1) * vpoints[E->mesh.nsd] + jj] = E->viscosity.N0[E->mat[i] - 1];
        }
    }

    return;
}

void visc_from_T(E, Eta, EEta, propogate) struct All_variables *E;
double *Eta, *EEta;
int propogate;

{
    void remove_horiz_ave();
    int i, j, k, l, z, jj, kk, imark, newnum, markvnum, m, left, right, mid, mark_num_depth;
    double c1, c2, c3, e_6, eta0, Tave, depth, temp[9], rii[9], tempa, TT[9], CC[9];
    double CCdot[9], Cratio, dd[9];
    double temp1, sumMat, aveMat, temp_visc_T, temp_visc_N0, temp_visc_E, temp_visc_all, temp_ratio;
    double temp0, temp2, temp3, temp5, temp7, Bvisc;
    static int visits = 0;
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    const int nno = E->lmesh.nno;
    const int noz = E->lmesh.noz;
    const int nox = E->lmesh.nox;
    const int elz = E->lmesh.elz;
    const int elx = E->lmesh.elx;

    double visc_weakzone, viscincreaseslab = 1.0, tempdist, velobegin, viscreducephase = 1.0, viscreducec = 1.0;
    double visc_r = 1.0;
    double tempdist_left, center_x, center_z, temp_loc_mid, temp_loc, temp_dist, temp_visc_log, temp_visc_log1, temp_visc_in;
    int nz, newelz;
    double e_pressure, loc_mid, t1, r1;
    double area, XMCold[4], dX[4], weigh1, weigh2, weigh3, weigh4;
    double x_ocean_front[1000], x_ocean_front_min[1000];
    void element_markers();
    void get_C_from_markers_multi();
    static const double one = 1.0;
    static const double zero = 0.0;
    static const double pt5 = 0.5;
    double global_fmin();
    void return_horiz_min();
    switch (E->viscosity.RHEOL)
    {
    case 1:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 1:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
        }

        if (E->control.phasevisc)
        {
            return_horiz_ave(E, E->T, E->Have.T);
        }

        for (i = 1; i <= nel; i++)
        {
            l = E->mat[i];
            tempa = E->viscosity.N0[l - 1];
            j = 0;

            for (kk = 1; kk <= ends; kk++)
            {
                TT[kk] = E->T[E->ien[i].node[kk]];
                CC[kk] = E->C[E->ien[i].node[kk]];
                if (E->control.phasevisc_C)
                {
                    CCdot[kk] = E->Cphasedot[E->ien[i].node[kk]] * E->control.phasevisc_C_time;
                }
                if (E->control.phasevisc_d)
                {
                    dd[kk] = E->d_node[E->ien[i].node[kk]];
                }
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                temp0 = 1.0e-32;
                temp3 = 1.0e-32;
                temp5 = 0.0;
                temp7 = 0.0;
                visc_weakzone = 1.0;
                viscincreaseslab = 1.0;
                viscreducephase = 1.0;
                viscreducec = 1.0;
                for (kk = 1; kk <= ends; kk++)
                {
                    temp0 += min(TT[kk], one) * E->N.vpt[GNVINDEX(kk, jj)];
                    if (E->control.phasevisc_C)
                    {
                        temp3 += CCdot[kk] * E->N.vpt[GNVINDEX(kk, jj)];
                    }
                    if (E->control.phasevisc_d)
                    {
                        temp5 += dd[kk] * E->N.vpt[GNVINDEX(kk, jj)];
                    }
                    if (E->control.Visc_C)
                    {
                        temp7 += CC[kk] * E->N.vpt[GNVINDEX(kk, jj)];
                    }
                }
                if (E->control.visc_leftcor)
                {
                    if (E->X[1][E->ien[i].node[jj]] <= E->control.x_weakzone_leftcor && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_leftcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_rightcor)
                {
                    if (E->X[1][E->ien[i].node[jj]] >= E->control.x_weakzone_rightcor && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_rightcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_mid)
                {
                    if (E->X[1][E->ien[i].node[jj]] >= E->control.x_weakzone_mid_left && E->X[1][E->ien[i].node[jj]] <= E->control.x_weakzone_mid_right && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_mid)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }

                if (E->control.visc_mid_dip)
                {
                    tempdist = (E->X[1][E->ien[i].node[jj]] - E->control.x_weakzone_mid_left) * tan(E->control.dip_weakzone_mid * 3.14159265 / 180.0) + 1.0 - E->X[2][E->ien[i].node[jj]];
                    if (tempdist >= 0.0 && tempdist <= (E->control.x_weakzone_mid_right - E->control.x_weakzone_mid_left) * tan(E->control.dip_weakzone_mid * 3.14159265 / 180.0) && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_mid)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }

                if (E->viscosity.visc_platebond)
                {
                    if (E->control.imposevelo)
                    {

                        velobegin = E->VB[1][noz * 7];
                        newnum = 8;
                        while (newnum <= nox - 1 && E->VB[1][noz * newnum] > velobegin - 910.0)
                        {
                            newnum++;
                        }
                        markvnum = newnum;
                        if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond && E->X[1][E->ien[i].node[jj]] - E->X[1][markvnum * noz] >= 0.0 - E->viscosity.width_weakzone_platebond && E->X[1][E->ien[i].node[jj]] - E->X[1][markvnum * noz] <= E->viscosity.width_weakzone_platebond)
                        {
                            visc_weakzone = E->control.visc_weakzone;
                        }

                    } /* end of control.imposevelo */
                    if (E->control.trechmigrate)
                    {
                        loc_mid = E->control.velo_surf_loc_mid;
                        loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
                        t1 = E->X[1][E->ien[i].node[jj]];
                        r1 = E->X[2][E->ien[i].node[jj]];
                        if (E->control.initialTOption == 1)
                        {
                            tempdist = (t1 - loc_mid) * tan(E->control.dip_margin * M_PI / 180.0) + 1.0 - r1;
                            tempdist_left = (t1 - loc_mid) * tan(E->control.dip_margin_left * M_PI / 180.0) + 1.0 - r1;
                            if (tempdist_left >= (0.0 - E->viscosity.left_weakzone_platebond) *
                                                     tan(E->control.dip_margin_left * M_PI / 180.0) &&
                                tempdist <= E->viscosity.right_weakzone_platebond *
                                                tan(E->control.dip_margin * M_PI / 180.0) &&
                                r1 >= E->viscosity.z_weakzone_platebond)
                            {
                                visc_weakzone = E->control.visc_weakzone;
                            }
                        }
                        if (E->control.initialTOption == 2)
                        {
                            center_x = loc_mid + E->control.dip_center_x;
                            center_z = 1.0 - E->control.dip_center_z;
                            if (r1 >= E->viscosity.z_weakzone_platebond && t1 < center_x)
                            {
                                tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
                                if (tempdist >= 0.0 - E->viscosity.left_weakzone_platebond && tempdist <= E->viscosity.right_weakzone_platebond)
                                {
                                    visc_weakzone = E->control.visc_weakzone;
                                }
                            }
                        } /*end of oceanic plate margin*/
                    }     /* end of trenchmigrate */

                } /* end of platebond*/

                if (E->control.slab_visc)
                {
                    if (E->X[2][E->ien[i].node[jj]] >= 1.0 - E->control.slab_visc_depth && E->X[2][E->ien[i].node[jj]] <= E->viscosity.zlith)
                    {
                        if (temp0 < E->control.inter_temp - E->control.temp_slabvisc)
                        {
                            viscincreaseslab = E->control.viscincreaseslab;
                        }
                    }
                }
                /*phase visc*/
                if (E->control.phasevisc)
                {
                    if (E->X[2][E->ien[i].node[jj]] < E->viscosity.zlm + E->control.phaseTop && E->X[2][E->ien[i].node[jj]] > E->viscosity.zlm - E->control.phaseBot)
                    {
                        newelz = ((i - 1) % elz) + 1;
                        nz = (E->ien[i].node[jj] % noz) + 1;

                        e_pressure = E->viscosity.zlm - E->X[2][E->ien[i].node[jj]] - E->control.clapeyron670 * (temp0 - E->control.transT670);
                        Bvisc = pt5 * (one + tanh(E->control.width670 * e_pressure));
                        /*printf("%lf %lf %lf %lf\n",E->X[2][E->ien[i].node[jj]],Bvisc,E->data.layer_meter/E->control.width670,e_pressure);*/
                        if (Bvisc < 0.5 && Bvisc > 0)
                        {
                            l = 4;
                        }
                        if (Bvisc < E->control.phaseffactorLa && Bvisc > E->control.phaseffactorSm)
                        {
                            if (E->control.phasevisc_slab)
                            {
                                if (temp0 <= E->Have.T[nz] + E->control.phasevisc_slab_T)
                                {
                                    viscreducephase = E->control.ViscReduce;
                                    l = 4;
                                }
                            }
                            else
                            {
                                viscreducephase = E->control.ViscReduce;
                                l = 4;
                                /*printf("%lf %lf\n",E->viscosity.N0[l-1], E->viscosity.N0[l-1]*viscreducephase*visc_weakzone*viscincreaseslab*exp(E->viscosity.E[l-1]*(E->viscosity.T[l-1]-temp0)));*/
                            }
                        }
                    }
                }
                /*end phase visc*/
                /*starting of phase viscC */
                if (E->control.phasevisc_d)
                {
                    if (E->X[2][E->ien[i].node[jj]] < E->viscosity.zlm + E->control.phaseTop && E->X[2][E->ien[i].node[jj]] > E->viscosity.zlm - E->control.phaseBot)
                    {
                        newelz = ((i - 1) % elz) + 1;
                        nz = (E->ien[i].node[jj] % noz) + 1;

                        e_pressure = E->viscosity.zlm - E->X[2][E->ien[i].node[jj]] - E->control.clapeyron670 * (temp0 - E->control.transT670);
                        Bvisc = pt5 * (one + tanh(E->control.width670 * e_pressure));
                        /*                        if (E->X[2][E->ien[i].node[jj]]-0.76655<0.002&&E->X[2][E->ien[i].node[jj]]-0.76655>-0.002)
                         printf("%lf %lf %lf %lf\n",E->X[2][E->ien[i].node[jj]],Bvisc,E->data.layer_meter/E->control.width670,e_pressure);
*/
                        if (Bvisc < 0.5 && Bvisc > 0)
                        {
                            l = 4;
                        }
                        if (Bvisc < E->control.phaseffactorLa && Bvisc > E->control.phaseffactorSm)
                        {
                            viscreducephase = powf(temp5 / E->control.phasevisc_d0, E->control.phasevisc_dp);
                            if (viscreducephase < E->control.ViscReduce)
                            {
                                viscreducephase = E->control.ViscReduce;
                            }
                            l = 5;
                        }
                    }
                }

                else if (E->control.phasevisc_C)
                {
                    if (E->X[2][E->ien[i].node[jj]] < E->viscosity.zlm + E->control.phaseTop && E->X[2][E->ien[i].node[jj]] > E->viscosity.zlm - E->control.phaseBot)
                    {
                        newelz = ((i - 1) % elz) + 1;
                        nz = (E->ien[i].node[jj] % noz) + 1;

                        e_pressure = E->viscosity.zlm - E->X[2][E->ien[i].node[jj]] - E->control.clapeyron670 * (temp0 - E->control.transT670);
                        Bvisc = pt5 * (one + tanh(E->control.width670 * e_pressure));
                        /*printf("%lf %lf %lf %lf\n",E->X[2][E->ien[i].node[jj]],Bvisc,E->data.layer_meter/E->control.width670,e_pressure);*/
                        if (Bvisc < 0.5 && Bvisc > 0)
                        {
                            l = 4;
                        }
                        if (Bvisc < E->control.phaseffactorLa && Bvisc > E->control.phaseffactorSm)
                        {
                            if (E->control.phasevisc_slab)
                            {
                                if (temp0 <= E->Have.T[nz] + E->control.phasevisc_slab_T)
                                {
                                    viscreducephase = E->control.ViscReduce;
                                    l = 4;
                                }
                            }
                            else
                            {
                                viscreducephase = E->control.ViscReduce / (0.001 + temp3);
                                l = 4;
                            }
                        }
                    }
                }
                /*  end of phase viscC */
                /*  start C_visc */
                if (E->control.Visc_C && E->X[2][E->ien[i].node[jj]] < 1.0 - E->control.depth_ocean_lith)
                {
                    if (temp7 >= 0.05)
                    {
                        viscreducec = E->control.ViscReduce_C;
                    }
                }
                /*  end of */

                temp2 = E->viscosity.E[l - 1] * (E->viscosity.T[l - 1] - temp0);
                tempa = E->viscosity.N0[l - 1];
                if (temp0 >= E->control.inter_temp && E->X[2][E->ien[i].node[jj]] > E->viscosity.zlith)
                {
                    /* fprintf(stderr,"%lf %lf %lf\n",E->X[1][E->ien[i].node[jj]]*2870,(1-E->X[2][E->ien[i].node[jj]])*2870,E->viscosity.N0[1]); */
                    tempa = E->viscosity.N0[1];
                }
                EEta[(i - 1) * vpts + jj] = tempa * viscreducephase * visc_weakzone * viscincreaseslab * viscreducec * exp(temp2);

                if (E->control.visc_const_cor)
                {
                    if (E->control.visc_leftcor)
                    {
                        if (E->X[1][E->ien[i].node[jj]] <= E->control.x_weakzone_leftcor && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_leftcor)
                        {
                            EEta[(i - 1) * vpts + jj] = E->control.visc_weakzone;
                        }
                    }
                    if (E->control.visc_rightcor)
                    {
                        if (E->X[1][E->ien[i].node[jj]] >= E->control.x_weakzone_rightcor && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_rightcor)
                        {
                            EEta[(i - 1) * vpts + jj] = E->control.visc_weakzone;
                        }
                    }
                    if (E->control.visc_mid)
                    {
                        if (E->X[1][E->ien[i].node[jj]] >= E->control.x_weakzone_mid_left && E->X[1][E->ien[i].node[jj]] <= E->control.x_weakzone_mid_right && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_mid)
                        {
                            EEta[(i - 1) * vpts + jj] = E->control.visc_weakzone;
                        }
                    }
                    if (E->control.visc_mid_dip)
                    {
                        tempdist = (E->X[1][E->ien[i].node[jj]] - E->control.x_weakzone_mid_left) * tan(E->control.dip_weakzone_mid * 3.14159265 / 180.0) + 1.0 - E->X[2][E->ien[i].node[jj]];
                        if (tempdist >= 0.0 && tempdist <= (E->control.x_weakzone_mid_right - E->control.x_weakzone_mid_left) * tan(E->control.dip_weakzone_mid * 3.14159265 / 180.0) && E->X[2][E->ien[i].node[jj]] >= E->control.z_weakzone_mid)
                        {
                            EEta[(i - 1) * vpts + jj] = E->control.visc_weakzone;
                        }
                    } //end visc_mid_dip

                    if (E->viscosity.visc_platebond)
                    {
                        if (E->control.imposevelo)
                        {
                            velobegin = E->VB[1][noz * 7];
                            newnum = 8;
                            while (newnum <= nox - 1 && E->VB[1][noz * newnum] > velobegin - 910.0)
                            {
                                newnum++;
                            }
                            markvnum = newnum;
                            if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond && E->X[1][E->ien[i].node[jj]] - E->X[1][markvnum * noz] >= 0.0 - E->viscosity.width_weakzone_platebond && E->X[1][E->ien[i].node[jj]] - E->X[1][markvnum * noz] <= E->viscosity.width_weakzone_platebond)
                            {
                                EEta[(i - 1) * vpts + jj] = E->control.visc_weakzone;
                            }
                        } /* end of imoposevelo */
                        if (E->control.trechmigrate)
                        {
                            loc_mid = E->control.velo_surf_loc_mid;
                            loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
                            t1 = E->X[1][E->ien[i].node[jj]];
                            r1 = E->X[2][E->ien[i].node[jj]];
                            if (E->control.initialTOption == 1)
                            {
                                tempdist = (t1 - loc_mid) * tan(E->control.dip_margin * 3.14159265 / 180.0) + 1.0 - r1;
                                if (tempdist >= (0.0 - E->viscosity.left_weakzone_platebond) * tan(E->control.dip_margin * 3.14159265 / 180.0) && tempdist <= E->viscosity.right_weakzone_platebond * tan(E->control.dip_margin * 3.14159265 / 180.0) && r1 >= E->viscosity.z_weakzone_platebond)
                                {
                                    EEta[(i - 1) * vpts + jj] = E->control.visc_weakzone;
                                }
                            }
                            if (E->control.initialTOption == 2)
                            {
                                center_x = loc_mid + E->control.dip_center_x;
                                center_z = 1.0 - E->control.dip_center_z;
                                if (r1 >= E->viscosity.z_weakzone_platebond && t1 < center_x)
                                {
                                    tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
                                    if (tempdist >= 0.0 - E->viscosity.left_weakzone_platebond && tempdist <= E->viscosity.right_weakzone_platebond)
                                    {
                                        EEta[(i - 1) * vpts + jj] = E->control.visc_weakzone;
                                    }
                                }
                            } /*end of oceanic plate margin*/
                        }     /* end of trenchmigrate */

                    } //end visc_platebond

                } //end control.visc_const_cor
            }     //end jj

        } //end i

        break;

    case 2:
        for (i = 1; i <= nel; i++)
        {
            l = E->mat[i];
            tempa = E->viscosity.N0[l - 1];

            j = 0;

            for (kk = 1; kk <= ends; kk++)
                TT[kk] = E->T[E->ien[i].node[kk]];

            for (jj = 1; jj <= vpts; jj++)
            {
                rii[jj] = 0;
                temp[jj] = 0;
                for (kk = 1; kk <= ends; kk++)
                {
                    rii[jj] += E->X[2][E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)];
                    temp[jj] += E->T[E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)];
                }
                temp2 = E->viscosity.E[l - 1] * (E->viscosity.T[l - 1] - 0.5);

                EEta[(i - 1) * vpoints[E->mesh.nsd] + jj] = E->viscosity.N0[E->mat[i] - 1];
            }
        }

        break;

    case 3:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 3:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
        }

        for (i = 1; i <= nel; i++)
        {
            l = E->mat[i];
            tempa = E->viscosity.N0[l - 1];

            j = 0;

            for (kk = 1; kk <= ends; kk++)
                TT[kk] = E->T[E->ien[i].node[kk]];

            for (jj = 1; jj <= vpts; jj++)
            {
                rii[jj] = 0;
                temp[jj] = 0;
                for (kk = 1; kk <= ends; kk++)
                {
                    rii[jj] += E->X[2][E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)];
                    temp[jj] += E->T[E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)];
                }
                temp1 = (E->viscosity.E[E->mat[i] - 1] +
                         E->viscosity.Z[E->mat[i] - 1] * (1.0 - rii[jj])) /
                            (temp[jj] + E->control.Ts) -
                        (E->viscosity.E[E->mat[i] - 1] +
                         E->viscosity.Z[E->mat[i] - 1]) /
                            (1.0 + E->control.Ts);
                EEta[(i - 1) * vpoints[E->mesh.nsd] + jj] = E->viscosity.N0[E->mat[i] - 1] * exp(temp1);
            }
        }
        break;
    case 10:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 10:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
        }

        if (E->control.phasefile_C || E->control.phasefile_Complete)
        {
            get_C_from_markers_multi(E, E->C_phasefile_marker_int[0], E->C_phasefile_nno, E->C_phasefile_element, 0, E->control.phasefile_C_flavor, E->CElement);
        }
        if (E->viscosity.visc_platebond_selfadapt)
        {
            fprintf(stderr, "try to get weak zone\n");
            for (i = 1; i <= noz; i++)
            {
                x_ocean_front[i] = E->X[1][nno];
                //x_ocean_front_lith[i] = E->X[1][E->mesh.nno];
            }
            for (i = 1; i <= E->advection.markers; i++)
            {
                if (E->XMC[2][i] < E->viscosity.z_weakzone_platebond - 0.002 || E->C_phasefile_marker_int[0][i] == 0 || E->C_phasefile_marker_int[0][i] > 2)
                {
                    continue;
                }
                left = 1;
                right = noz;
                mid = (left + right) / 2;
                if (E->XMC[2][i] > 1.0)
                {
                    E->XMC[2][i] = 1.0;
                }
                if (E->XMC[2][i] < 0)
                {
                    E->XMC[2][i] = 0;
                }
                while (left < right - 1)
                {
                    if (E->XMC[2][i] < E->X[2][mid])
                    {
                        right = mid;
                    }
                    else
                    {
                        left = mid;
                    }
                    mid = (left + right) / 2;
                }
                if (E->XMC[2][i] < E->X[2][mid])
                {
                    mark_num_depth = mid;
                }
                else
                {
                    mark_num_depth = mid + 1;
                }
                if (E->XMC[1][i] < x_ocean_front[mark_num_depth])
                {
                    x_ocean_front[mark_num_depth] = E->XMC[1][i];
                }
            }
        }

        x_ocean_front[0] = x_ocean_front[E->mesh.noz];
        for (i = E->mesh.noz; i >= 1; i--)
        {
            if (E->X[2][i] >= E->viscosity.z_weakzone_platebond)
            {
                if (i < E->mesh.noz)
                {
                    if (x_ocean_front[i - 1] - x_ocean_front[i] > 100.0 / 2870 || x_ocean_front[i - 1] - x_ocean_front[i] < -100.0 / 2870)
                    {
                        x_ocean_front[i - 1] = x_ocean_front[i];
                    }
                }
                fprintf(stderr, "ocean_front %d %lf %lf\n", i, E->X[2][i], x_ocean_front[i] * 2870 - 2870);
            }
        }
        for (i = 1; i <= nel; i++)
        {

            j = 0;
            //fprintf(stderr, "%d %lf %lf %lf\n", i, temp_visc_N0, temp_visc_E, temp_visc_T);

            for (kk = 1; kk <= ends; kk++)
            {
                TT[kk] = E->T[E->ien[i].node[kk]];
                CC[kk] = 1.0;
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                temp0 = 1.0e-32;
                temp3 = 1.0e-32;
                temp5 = 0.0;
                temp7 = 0.0;
                visc_weakzone = 1.0;
                viscincreaseslab = 1.0;
                viscreducephase = 1.0;
                viscreducec = 1.0;
                if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond)
                {
                    if (E->X[1][E->ien[i].node[jj]] <= x_ocean_front[E->ien[i].node[jj] % E->mesh.noz] &&
                        E->X[1][E->ien[i].node[jj]] >= x_ocean_front[E->ien[i].node[jj] % E->mesh.noz] - E->viscosity.left_weakzone_platebond)
                    {
                        visc_weakzone = E->viscosity.visc_reduce_platebond;
                    }
                }
                if (E->control.visc_leftcor)
                {
                    if (t1 <= E->control.x_weakzone_leftcor && r1 >= E->control.z_weakzone_leftcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_rightcor)
                {
                    if (t1 >= E->control.x_weakzone_rightcor && r1 >= E->control.z_weakzone_rightcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                for (kk = 1; kk <= ends; kk++)
                {
                    temp0 += min(TT[kk], one) * E->N.vpt[GNVINDEX(kk, jj)];
                }
                for (m = 0; m < E->control.phasefile_C_flavor; m++)
                {
                    if (m == 0)
                    {
                        temp_visc_all = 0;
                    }
                    temp2 = E->viscosity.E[m] * (E->viscosity.T[m] - temp0);
                    temp3 = E->viscosity.N0[m] * exp(temp2);
                    if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond && m == 1)
                    {
                        temp3 *= 100.0;
                    }
                    temp_visc_all += log(temp3) * E->C_phasefile_element[m][i];
                    //fprintf(stderr, "%d %lf %lf %lf\n", m, E->C_phasefile_element[m][i], E->viscosity.N0[m], E->viscosity.E[m]);
                }
                EEta[(i - 1) * vpts + jj] = visc_weakzone * exp(temp_visc_all);

                if (E->viscosity.visc_platebond_const)
                {
                    if (visc_weakzone < 0.9)
                    {
                        EEta[(i - 1) * vpts + jj] = E->viscosity.visc_reduce_platebond;
                    }
                }
            }
        }
        break;
    case 11:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 10:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
        }
        if (E->control.phasefile_C || E->control.phasefile_Complete)
        {
            get_C_from_markers_multi(E, E->C_phasefile_marker_int[0], E->C_phasefile_nno, E->C_phasefile_element, 0, E->control.phasefile_C_flavor, E->CElement);
        }
        if (E->viscosity.visc_platebond_selfadapt)
        {
            fprintf(stderr, "try to get weak zone\n");
            for (i = 1; i <= E->mesh.noz; i++)
            {
                x_ocean_front[i] = E->X[1][E->mesh.nno];
                //x_ocean_front_lith[i] = E->X[1][E->mesh.nno];
            }
            for (i = 1; i <= E->advection.markers; i++)
            {
                if (E->XMC[2][i] < E->viscosity.z_weakzone_platebond - 0.002 || E->C_phasefile_marker_int[0][i] == 0 || E->C_phasefile_marker_int[0][i] > 2)
                {
                    continue;
                }
                left = 1;
                right = E->mesh.noz;
                mid = (left + right) / 2;
                if (E->XMC[2][i] > 1.0)
                {
                    E->XMC[2][i] = 1.0;
                }
                if (E->XMC[2][i] < 0)
                {
                    E->XMC[2][i] = 0;
                }
                while (left < right - 1)
                {
                    if (E->XMC[2][i] < E->X[2][mid])
                    {
                        right = mid;
                    }
                    else
                    {
                        left = mid;
                    }
                    mid = (left + right) / 2;
                }
                if (E->XMC[2][i] < E->X[2][mid])
                {
                    mark_num_depth = mid;
                }
                else
                {
                    mark_num_depth = mid + 1;
                }
                if (E->XMC[1][i] < x_ocean_front[mark_num_depth])
                {
                    x_ocean_front[mark_num_depth] = E->XMC[1][i];
                }
            }
        }

        x_ocean_front[0] = x_ocean_front[E->mesh.noz];
        for (i = E->mesh.noz; i >= 1; i--)
        {
            if (E->X[2][i] >= E->viscosity.z_weakzone_platebond)
            {
                if (i < E->mesh.noz)
                {
                    if (x_ocean_front[i - 1] - x_ocean_front[i] > 100.0 / 2870 || x_ocean_front[i - 1] - x_ocean_front[i] < -100.0 / 2870)
                    {
                        x_ocean_front[i - 1] = x_ocean_front[i];
                    }
                }
                fprintf(stderr, "ocean_front %d %lf %lf\n", i, E->X[2][i], x_ocean_front[i] * 2870 - 2870);
            }
        }
        for (i = 1; i <= nel; i++)
        {

            j = 0;
            //fprintf(stderr, "%d %lf %lf %lf\n", i, temp_visc_N0, temp_visc_E, temp_visc_T);

            for (kk = 1; kk <= ends; kk++)
            {
                TT[kk] = E->T[E->ien[i].node[kk]];
                CC[kk] = 1.0;
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                temp0 = 1.0e-32;
                temp3 = 1.0e-32;
                temp5 = 0.0;
                temp7 = 0.0;
                visc_weakzone = 1.0;
                viscincreaseslab = 1.0;
                viscreducephase = 1.0;
                viscreducec = 1.0;
                if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond)
                {
                    if (E->X[1][E->ien[i].node[jj]] <= x_ocean_front[E->ien[i].node[jj] % E->mesh.noz] &&
                        E->X[1][E->ien[i].node[jj]] >= x_ocean_front[E->ien[i].node[jj] % E->mesh.noz] - E->viscosity.left_weakzone_platebond)
                    {
                        visc_weakzone = E->viscosity.visc_reduce_platebond;
                    }
                }
                if (E->control.visc_leftcor)
                {
                    if (t1 <= E->control.x_weakzone_leftcor && r1 >= E->control.z_weakzone_leftcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_rightcor)
                {
                    if (t1 >= E->control.x_weakzone_rightcor && r1 >= E->control.z_weakzone_rightcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                for (kk = 1; kk <= ends; kk++)
                {
                    temp0 += min(TT[kk], one) * E->N.vpt[GNVINDEX(kk, jj)];
                }
                for (m = 0; m < E->control.phasefile_C_flavor; m++)
                {
                    if (m == 0)
                    {
                        temp_visc_all = 0;
                    }
                    temp2 = E->viscosity.E[m] * (1 / (temp0 + E->viscosity.T[m]) - 1 / (E->viscosity.T[m] + 1));
                    temp3 = E->viscosity.N0[m] * exp(temp2);
                    if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond && m == 1)
                    {
                        temp3 *= 100.0;
                    }
                    temp_visc_all += log(temp3) * E->C_phasefile_element[m][i];
                    //fprintf(stderr, "%d %lf %lf %lf\n", m, E->C_phasefile_element[m][i], E->viscosity.N0[m], E->viscosity.E[m]);
                }
                EEta[(i - 1) * vpts + jj] = visc_weakzone * exp(temp_visc_all);

                if (E->viscosity.visc_platebond_const)
                {
                    if (visc_weakzone < 0.9)
                    {
                        EEta[(i - 1) * vpts + jj] = E->viscosity.visc_reduce_platebond;
                    }
                }
            }
        }
    case 20:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 20:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
        }
        if (E->control.phasefile_C || E->control.phasefile_Complete)
        {
            get_C_from_markers_multi(E, E->C_phasefile_marker_int[0], E->C_phasefile_nno, E->C_phasefile_element, 0, E->control.phasefile_C_flavor, E->CElement);
        }
        if (E->viscosity.visc_platebond_selfadapt)
        {
            fprintf(stderr, "try to get weak zone\n");
            for (i = 1; i <= E->mesh.noz; i++)
            {
                x_ocean_front[i] = E->X[1][E->mesh.nno];
                //x_ocean_front_lith[i] = E->X[1][E->mesh.nno];
            }
            for (i = 1; i <= E->advection.markers; i++)
            {
                if (E->XMC[2][i] < E->viscosity.z_weakzone_platebond - 0.002 || E->C_phasefile_marker_int[0][i] == 0 || E->C_phasefile_marker_int[0][i] > 2)
                {
                    continue;
                }
                left = 1;
                right = E->mesh.noz;
                mid = (left + right) / 2;
                if (E->XMC[2][i] > 1.0)
                {
                    E->XMC[2][i] = 1.0;
                }
                if (E->XMC[2][i] < 0)
                {
                    E->XMC[2][i] = 0;
                }
                while (left < right - 1)
                {
                    if (E->XMC[2][i] < E->X[2][mid])
                    {
                        right = mid;
                    }
                    else
                    {
                        left = mid;
                    }
                    mid = (left + right) / 2;
                }
                if (E->XMC[2][i] < E->X[2][mid])
                {
                    mark_num_depth = mid;
                }
                else
                {
                    mark_num_depth = mid + 1;
                }
                if (E->XMC[1][i] < x_ocean_front[mark_num_depth])
                {
                    x_ocean_front[mark_num_depth] = E->XMC[1][i];
                }
            }
        }

        x_ocean_front[0] = x_ocean_front[E->mesh.noz];
        for (i = E->mesh.noz; i >= 1; i--)
        {
            if (E->X[2][i] >= E->viscosity.z_weakzone_platebond)
            {
                if (i < E->mesh.noz)
                {
                    if (x_ocean_front[i - 1] - x_ocean_front[i] > 100.0 / 2870 || x_ocean_front[i - 1] - x_ocean_front[i] < -100.0 / 2870)
                    {
                        x_ocean_front[i - 1] = x_ocean_front[i];
                    }
                }
                fprintf(stderr, "ocean_front %d %lf %lf\n", i, E->X[2][i], x_ocean_front[i] * 2870 - 2870);
            }
        }
        for (i = 1; i <= nel; i++)
        {

            j = 0;
            //fprintf(stderr, "%d %lf %lf %lf\n", i, temp_visc_N0, temp_visc_E, temp_visc_T);

            for (kk = 1; kk <= ends; kk++)
            {
                TT[kk] = E->T[E->ien[i].node[kk]];
                CC[kk] = 1.0;
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                temp0 = 1.0e-32;
                temp3 = 1.0e-32;
                temp5 = 0.0;
                temp7 = 0.0;
                visc_weakzone = 1.0;
                viscincreaseslab = 1.0;
                viscreducephase = 1.0;
                viscreducec = 1.0;
                if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond)
                {
                    temp_loc = x_ocean_front[E->ien[i].node[jj] % E->mesh.noz];
                    temp_loc_mid = temp_loc - pt5 * E->viscosity.left_weakzone_platebond;
                    temp_dist = E->X[1][E->ien[i].node[jj]] - temp_loc_mid;
                    temp_visc_log = log(E->viscosity.visc_reduce_platebond);
                    temp_visc_in = cos(M_PI / 2.0 * tanh(temp_dist / E->viscosity.left_weakzone_platebond * E->viscosity.left_weakzone_platebond_ratio));
                    visc_weakzone = exp(temp_visc_in * temp_visc_log);
                }
                if (E->control.visc_leftcor)
                {
                    if (t1 <= E->control.x_weakzone_leftcor && r1 >= E->control.z_weakzone_leftcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_rightcor)
                {
                    if (t1 >= E->control.x_weakzone_rightcor && r1 >= E->control.z_weakzone_rightcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                for (kk = 1; kk <= ends; kk++)
                {
                    temp0 += min(TT[kk], one) * E->N.vpt[GNVINDEX(kk, jj)];
                }
                for (m = 0; m < E->control.phasefile_C_flavor; m++)
                {
                    if (m == 0)
                    {
                        temp_visc_all = 0;
                    }
                    temp2 = E->viscosity.E[m] * (E->viscosity.T[m] - temp0);
                    temp3 = E->viscosity.N0[m] * exp(temp2);
                    if (E->X[2][E->ien[i].node[jj]] >= 1.0 - E->viscosity.z_weakzone_platebond && m == 1)
                    {
                        temp3 *= 100.0;
                    }
                    temp_visc_all += log(temp3) * E->C_phasefile_element[m][i];
                    //fprintf(stderr, "%d %lf %lf %lf\n", m, E->C_phasefile_element[m][i], E->viscosity.N0[m], E->viscosity.E[m]);
                }
                EEta[(i - 1) * vpts + jj] = visc_weakzone * exp(temp_visc_all);

                if (E->viscosity.visc_platebond_const)
                {
                    if (visc_weakzone < 0.9)
                    {
                        EEta[(i - 1) * vpts + jj] = E->viscosity.visc_reduce_platebond;
                    }
                }
            }
        }
        break;

    case 30:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 30:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
            for (l = 0; l <= noz; l++)
            {
                x_ocean_front[i] = 0;
            }
        }

        //                fprintf(stderr, "viscosity inside%d \n",E->parallel.me);

        if (E->control.phasefile_C || E->control.phasefile_Complete)
        {
            /*
                fprintf(stderr, "before element marker%d \n",E->parallel.me);

            element_markers(E, 1);
                fprintf(stderr, "after element marker%d \n",E->parallel.me);

            get_C_from_markers_multi(E, E->C_phasefile_marker_int[0], E->C_phasefile_nno, E->C_phasefile_element, 0, E->control.phasefile_C_flavor, E->CElement);
*/
            if (E->viscosity.visc_platebond_selfadapt && E->parallel.me_loc[2] == E->parallel.nprocz - 1)
            {
                fprintf(stderr, "try to get weak zone\n");
                for (i = 1; i <= E->lmesh.noz; i++)
                {
                    x_ocean_front[i] = E->mesh.layer[1];
                    //x_ocean_front_lith[i] = E->X[1][E->mesh.nno];
                }
                for (i = 1; i <= E->advection.markers; i++)
                {
                    if (E->XMC[2][i] < E->viscosity.z_weakzone_platebond - 0.002 || E->C_phasefile_marker_int[0][i] == 0 || E->C_phasefile_marker_int[0][i] > 2)
                    {
                        continue;
                    }
                    left = 1;
                    right = E->lmesh.noz;
                    mid = (left + right) / 2;
                    if (E->XMC[2][i] > 1.0)
                    {
                        E->XMC[2][i] = 1.0;
                    }
                    if (E->XMC[2][i] < 0)
                    {
                        E->XMC[2][i] = 0;
                    }
                    while (left < right - 1)
                    {
                        if (E->XMC[2][i] < E->X[2][mid])
                        {
                            right = mid;
                        }
                        else
                        {
                            left = mid;
                        }
                        mid = (left + right) / 2;
                    }
                    if (E->XMC[2][i] < E->X[2][mid])
                    {
                        mark_num_depth = mid;
                    }
                    else
                    {
                        mark_num_depth = mid + 1;
                    }
                    if (E->XMC[1][i] < x_ocean_front[mark_num_depth])
                    {
                        x_ocean_front[mark_num_depth] = E->XMC[1][i];
                    }
                }
            }
            return_horiz_min(E, x_ocean_front, x_ocean_front_min, noz + 1);
            x_ocean_front_min[0] = x_ocean_front_min[noz];
            for (i = noz; i >= 1; i--)
            {
                if (E->X[2][i] >= E->viscosity.z_weakzone_platebond)
                {
                    if (i < noz)
                    {
                        if (x_ocean_front_min[i - 1] - x_ocean_front_min[i] > 100.0 / 2870 || x_ocean_front_min[i - 1] - x_ocean_front_min[i] < -100.0 / 2870)
                        {
                            x_ocean_front_min[i - 1] = x_ocean_front_min[i];
                        }
                    }
                    if (E->parallel.me_loc[1] == 1)
                        fprintf(stderr, "ocean_front %d %lf %lf\n", i, E->X[2][i], x_ocean_front_min[i] * 2870 - 2870);
                }
                else
                    break;
            }
        }

        for (i = 1; i <= nel; i++)
        {

            j = 0;
            //fprintf(stderr, "%d %lf %lf %lf\n", i, temp_visc_N0, temp_visc_E, temp_visc_T);

            for (kk = 1; kk <= ends; kk++)
            {
                TT[kk] = E->T[E->ien[i].node[kk]];
                CC[kk] = 1.0;
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                temp0 = 1.0e-32;
                temp3 = 1.0e-32;
                temp5 = 0.0;
                temp7 = 0.0;
                visc_weakzone = 1.0;
                viscincreaseslab = 1.0;
                viscreducephase = 1.0;
                viscreducec = 1.0;
                t1 = E->X[1][E->ien[i].node[jj]];
                r1 = E->X[2][E->ien[i].node[jj]];
                if (E->control.phasefile_C || E->control.phasefile_Complete)
                {
                    if (r1 >= E->viscosity.z_weakzone_platebond)
                    {
                        temp_loc = x_ocean_front_min[E->ien[i].node[jj] % E->lmesh.noz];
                        temp_ratio = (E->viscosity.weakzone_ratio - 1.0) * (r1 - E->viscosity.z_weakzone_platebond) / (1.0 - E->viscosity.z_weakzone_platebond) + 1.0;
                        temp_loc_mid = temp_loc - pt5 * E->viscosity.left_weakzone_platebond * temp_ratio;
                        temp_dist = t1 - temp_loc_mid;
                        temp_visc_log = log(E->viscosity.visc_reduce_platebond);
                        temp_visc_in = cos(M_PI / 2.0 * tanh(temp_dist / E->viscosity.left_weakzone_platebond * E->viscosity.left_weakzone_platebond_ratio / temp_ratio));
                        visc_weakzone = exp(temp_visc_in * temp_visc_log);
                    }
                }

                if (E->control.visc_leftcor)
                {
                    if (t1 <= E->control.x_weakzone_leftcor && r1 >= E->control.z_weakzone_leftcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_rightcor)
                {
                    if (t1 >= E->control.x_weakzone_rightcor && r1 >= E->control.z_weakzone_rightcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                for (kk = 1; kk <= ends; kk++)
                {
                    temp0 += min(TT[kk], one) * E->N.vpt[GNVINDEX(kk, jj)];
                }
                for (m = 0; m < E->control.phasefile_C_flavor; m++)
                {
                    if (m == 0)
                    {
                        temp_visc_all = 0;
                    }
                    temp2 = E->viscosity.E[m] * (E->viscosity.T[m] - temp0);
                    temp3 = E->viscosity.N0[m] * exp(temp2);
                    /*
                    if (r1 >= E->viscosity.zlm && r1 <= 1.0 - E->control.depth_ocean_lith && m >= 1)
                    {
                        temp3 *= E->viscosity.N0[0];
                    }
                    if (t1 <= temp_loc && r1 >= 1.0 - E->control.depth_continent_lith && m >= 1)
                    {
                        temp3 *= E->viscosity.N0[0];
                    }
*/
                    if (r1 >= 1.0 - E->viscosity.zcrust1 && m == 1)
                    {
                        temp3 *= E->viscosity.N0[0] / E->viscosity.N0[1];
                    }

                    temp_visc_all += log(temp3) * E->C_phasefile_element[m][i];
                    //fprintf(stderr, "%d %lf %lf %lf\n", m, E->C_phasefile_element[m][i], E->viscosity.N0[m], E->viscosity.E[m]);
                }
                EEta[(i - 1) * vpts + jj] = visc_weakzone * exp(temp_visc_all);
                //    fprintf(stderr, "visc CPU %d %d %.4e\n",E->parallel.me, (i - 1) * vpts + jj, EEta[(i - 1) * vpts + jj]);
                if (E->viscosity.visc_platebond_const && r1 > E->viscosity.z_weakzone_platebond)
                {
                    temp_loc = x_ocean_front_min[E->ien[i].node[jj] % E->lmesh.noz];
                    temp_ratio = (E->viscosity.weakzone_ratio - 1.0) * (r1 - E->viscosity.z_weakzone_platebond) / (1.0 - E->viscosity.z_weakzone_platebond) + 1.0;
                    temp_loc_mid = temp_loc - pt5 * E->viscosity.left_weakzone_platebond * temp_ratio;
                    temp_dist = t1 - temp_loc_mid;
                    if (r1 >= E->viscosity.z_weakzone_platebond && temp_dist >= -pt5 * E->viscosity.left_weakzone_platebond * temp_ratio && temp_dist <= pt5 * E->viscosity.left_weakzone_platebond * temp_ratio)
                    {
                        EEta[(i - 1) * vpts + jj] = E->viscosity.visc_reduce_platebond;
                    }
                }
            }
        }
        break;
    case -1:
    {
        for (i = 1; i <= nel; i++)
        {
            for (jj = 1; jj <= vpts; jj++)
            {
                EEta[(i - 1) * vpts + jj] = 1.0;
            }
        }
    }
    break;

    case -2:
    {
        for (i = 1; i <= nel; i++)
        {

            for (jj = 1; jj <= vpts; jj++)
            {
                t1 = E->X[1][E->ien[i].node[jj]];
                r1 = E->X[2][E->ien[i].node[jj]];
                if (r1 <= 0.2)
                    EEta[(i - 1) * vpts + jj] = 0.1;
                else
                    EEta[(i - 1) * vpts + jj] = 1.0;
            }
        }
    }
    break;
    case -3:
    {
        for (i = 1; i <= nel; i++)
        {

            for (jj = 1; jj <= vpts; jj++)
            {
                t1 = E->X[1][E->ien[i].node[jj]];
                r1 = E->X[2][E->ien[i].node[jj]];
                if (r1 <= 0.2)
                    EEta[(i - 1) * vpts + jj] = 0.01;
                else
                    EEta[(i - 1) * vpts + jj] = 1.0;
            }
        }
    }
    break;
    case 35:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 30:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
            for (l = 0; l <= noz; l++)
            {
                x_ocean_front[i] = 0;
            }
        }

        //                fprintf(stderr, "viscosity inside%d \n",E->parallel.me);

        if (E->control.phasefile_C || E->control.phasefile_Complete)
        {
            /*
                fprintf(stderr, "before element marker%d \n",E->parallel.me);

            element_markers(E, 1);
                fprintf(stderr, "after element marker%d \n",E->parallel.me);

            get_C_from_markers_multi(E, E->C_phasefile_marker_int[0], E->C_phasefile_nno, E->C_phasefile_element, 0, E->control.phasefile_C_flavor, E->CElement);
*/
            if (E->viscosity.visc_platebond_selfadapt && E->parallel.me_loc[2] == E->parallel.nprocz - 1)
            {
                fprintf(stderr, "try to get weak zone\n");
                for (i = 1; i <= E->lmesh.noz; i++)
                {
                    x_ocean_front[i] = E->mesh.layer[1];
                    //x_ocean_front_lith[i] = E->X[1][E->mesh.nno];
                }
                for (i = 1; i <= E->advection.markers; i++)
                {
                    if (E->XMC[2][i] < E->viscosity.z_weakzone_platebond - 0.002 || E->C_phasefile_marker_int[0][i] == 0 || E->C_phasefile_marker_int[0][i] > 2)
                    {
                        continue;
                    }
                    left = 1;
                    right = E->lmesh.noz;
                    mid = (left + right) / 2;
                    if (E->XMC[2][i] > 1.0)
                    {
                        E->XMC[2][i] = 1.0;
                    }
                    if (E->XMC[2][i] < 0)
                    {
                        E->XMC[2][i] = 0;
                    }
                    while (left < right - 1)
                    {
                        if (E->XMC[2][i] < E->X[2][mid])
                        {
                            right = mid;
                        }
                        else
                        {
                            left = mid;
                        }
                        mid = (left + right) / 2;
                    }
                    if (E->XMC[2][i] < E->X[2][mid])
                    {
                        mark_num_depth = mid;
                    }
                    else
                    {
                        mark_num_depth = mid + 1;
                    }
                    if (E->XMC[1][i] < x_ocean_front[mark_num_depth])
                    {
                        x_ocean_front[mark_num_depth] = E->XMC[1][i];
                    }
                }
            }
            return_horiz_min(E, x_ocean_front, x_ocean_front_min, noz + 1);
            x_ocean_front_min[0] = x_ocean_front_min[noz];
            for (i = noz; i >= 1; i--)
            {
                if (E->X[2][i] >= E->viscosity.z_weakzone_platebond)
                {
                    if (i < noz)
                    {
                        if (x_ocean_front_min[i - 1] - x_ocean_front_min[i] > 100.0 / 2870 || x_ocean_front_min[i - 1] - x_ocean_front_min[i] < -100.0 / 2870)
                        {
                            x_ocean_front_min[i - 1] = x_ocean_front_min[i];
                        }
                    }
                    if (E->parallel.me_loc[1] == 1)
                        fprintf(stderr, "ocean_front %d %lf %lf\n", i, E->X[2][i], x_ocean_front_min[i] * 2870 - 2870);
                }
                else
                    break;
            }
        }

        for (i = 1; i <= nel; i++)
        {

            j = 0;
            //fprintf(stderr, "%d %lf %lf %lf\n", i, temp_visc_N0, temp_visc_E, temp_visc_T);

            for (kk = 1; kk <= ends; kk++)
            {
                TT[kk] = E->T[E->ien[i].node[kk]];
                CC[kk] = 1.0;
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                temp0 = 1.0e-32;
                temp3 = 1.0e-32;
                temp5 = 0.0;
                temp7 = 0.0;
                visc_weakzone = 1.0;
                viscincreaseslab = 1.0;
                viscreducephase = 1.0;
                viscreducec = 1.0;
                t1 = E->X[1][E->ien[i].node[jj]];
                r1 = E->X[2][E->ien[i].node[jj]];
                if (E->control.phasefile_C || E->control.phasefile_Complete)
                {
                    if (r1 >= E->viscosity.z_weakzone_platebond)
                    {
                        temp_loc = x_ocean_front_min[E->ien[i].node[jj] % E->lmesh.noz];
                        temp_ratio = (E->viscosity.weakzone_ratio - 1.0) * (r1 - E->viscosity.z_weakzone_platebond) / (1.0 - E->viscosity.z_weakzone_platebond) + 1.0;
                        temp_loc_mid = temp_loc - pt5 * E->viscosity.left_weakzone_platebond * temp_ratio;
                        temp_dist = t1 - temp_loc_mid;
                        temp_visc_log = log(E->viscosity.visc_reduce_platebond);
                        temp_visc_in = cos(M_PI / 2.0 * tanh(temp_dist / E->viscosity.left_weakzone_platebond * E->viscosity.left_weakzone_platebond_ratio / temp_ratio));
                        visc_weakzone = exp(temp_visc_in * temp_visc_log);
                    }
                }

                if (E->control.visc_leftcor)
                {
                    if (t1 <= E->control.x_weakzone_leftcor && r1 >= E->control.z_weakzone_leftcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_rightcor)
                {
                    if (t1 >= E->control.x_weakzone_rightcor && r1 >= E->control.z_weakzone_rightcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                for (kk = 1; kk <= ends; kk++)
                {
                    temp0 += min(TT[kk], one) * E->N.vpt[GNVINDEX(kk, jj)];
                }
                for (m = 0; m < E->control.phasefile_C_flavor; m++)
                {
                    if (m == 0)
                    {
                        temp_visc_all = 0;
                    }
                    temp2 = E->viscosity.E[m] * (E->viscosity.T[m] - temp0);
                    temp3 = E->viscosity.N0[m] * exp(temp2);
                    /*
                    if (r1 >= E->viscosity.zlm && r1 <= 1.0 - E->control.depth_ocean_lith && m >= 1)
                    {
                        temp3 *= E->viscosity.N0[0];
                    }
                    if (t1 <= temp_loc && r1 >= 1.0 - E->control.depth_continent_lith && m >= 1)
                    {
                        temp3 *= E->viscosity.N0[0];
                    }
*/
                    if (r1 >= 1.0 - E->viscosity.zcrust1 && m == 1)
                    {
                        temp3 *= E->viscosity.N0[0] / E->viscosity.N0[1];
                    }

                    temp_visc_all += log(temp3) * E->C_phasefile_element[m][i];
                    //fprintf(stderr, "%d %lf %lf %lf\n", m, E->C_phasefile_element[m][i], E->viscosity.N0[m], E->viscosity.E[m]);
                }
                EEta[(i - 1) * vpts + jj] = visc_weakzone * exp(temp_visc_all);
                //    fprintf(stderr, "visc CPU %d %d %.4e\n",E->parallel.me, (i - 1) * vpts + jj, EEta[(i - 1) * vpts + jj]);
                if (E->viscosity.visc_platebond_const && r1 > E->viscosity.z_weakzone_platebond)
                {
                    temp_loc = x_ocean_front_min[E->ien[i].node[jj] % E->lmesh.noz];
                    temp_ratio = (E->viscosity.weakzone_ratio - 1.0) * (r1 - E->viscosity.z_weakzone_platebond) / (1.0 - E->viscosity.z_weakzone_platebond) + 1.0;
                    temp_loc_mid = temp_loc - pt5 * E->viscosity.left_weakzone_platebond * temp_ratio;
                    temp_dist = t1 - temp_loc_mid;
                    if (r1 >= E->viscosity.z_weakzone_platebond && temp_dist >= -pt5 * E->viscosity.left_weakzone_platebond * temp_ratio && temp_dist <= pt5 * E->viscosity.left_weakzone_platebond * temp_ratio)
                    {
                        EEta[(i - 1) * vpts + jj] = E->viscosity.visc_reduce_platebond;
                    }
                }
            }
        }
        break;

    case 100:
        if (propogate && visits == 0)
        {
            fprintf(E->fp, "\tRheological option 100:\n");

            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(E->fp, "\tlayer %d/%d: E=%lf T1=%lf \n",
                        l, E->viscosity.num_mat,
                        E->viscosity.E[l - 1], E->viscosity.T[l - 1]);
            }
            fflush(E->fp);
        }

        //                fprintf(stderr, "viscosity inside%d \n",E->parallel.me);

        for (i = 1; i <= nel; i++)
        {
            j = 0;
            //fprintf(stderr, "%d %lf %lf %lf\n", i, temp_visc_N0, temp_visc_E, temp_visc_T);
            for (kk = 1; kk <= ends; kk++)
            {
                TT[kk] = E->T[E->ien[i].node[kk]];
                CC[kk] = 1.0;
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                temp0 = 1.0e-32;
                temp3 = 1.0e-32;
                temp5 = 0.0;
                temp7 = 0.0;
                visc_weakzone = 1.0;
                viscincreaseslab = 1.0;
                viscreducephase = 1.0;
                viscreducec = 1.0;
                t1 = E->X[1][E->ien[i].node[jj]];
                r1 = E->X[2][E->ien[i].node[jj]];

                if (E->control.visc_leftcor)
                {
                    if (t1 <= E->control.x_weakzone_leftcor && r1 >= E->control.z_weakzone_leftcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->control.visc_rightcor)
                {
                    if (t1 >= E->control.x_weakzone_rightcor && r1 >= E->control.z_weakzone_rightcor)
                    {
                        visc_weakzone = E->control.visc_weakzone;
                    }
                }
                if (E->viscosity.visc_platebond)
                {
                    if (E->control.imposevelo)
                    {

                        velobegin = E->VB[1][noz * 7];
                        newnum = 8;
                        while (newnum <= nox - 1 && E->VB[1][noz * newnum] > velobegin - 910.0)
                        {
                            newnum++;
                        }
                        markvnum = newnum;
                        if (E->X[2][E->ien[i].node[jj]] >= E->viscosity.z_weakzone_platebond && E->X[1][E->ien[i].node[jj]] - E->X[1][markvnum * noz] >= 0.0 - E->viscosity.width_weakzone_platebond && E->X[1][E->ien[i].node[jj]] - E->X[1][markvnum * noz] <= E->viscosity.width_weakzone_platebond)
                        {
                            visc_weakzone = E->viscosity.visc_reduce_platebond;
                        }

                    } /* end of control.imposevelo */
                    if (E->control.trechmigrate)
                    {
                        loc_mid = E->control.velo_surf_loc_mid;
                        loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
                        t1 = E->X[1][E->ien[i].node[jj]];
                        r1 = E->X[2][E->ien[i].node[jj]];
                        switch (E->control.initialTOption)
                        {
                        case 1:
                            tempdist = (t1 - loc_mid) * tan(E->control.dip_margin * M_PI / 180.0) + 1.0 - r1;
                            tempdist_left = (t1 - loc_mid) * tan(E->control.dip_margin_left * M_PI / 180.0) + 1.0 - r1;
                            if (tempdist_left >= (0.0 - E->viscosity.left_weakzone_platebond) *
                                                     tan(E->control.dip_margin_left * M_PI / 180.0) &&
                                tempdist <= E->viscosity.right_weakzone_platebond *
                                                tan(E->control.dip_margin * M_PI / 180.0) &&
                                r1 >= E->viscosity.z_weakzone_platebond)
                            {
                                visc_weakzone = E->viscosity.visc_reduce_platebond;
                            }
                            break;
                        case 10:
                            center_x = loc_mid + E->control.dip_center_x;
                            center_z = 1.0 - E->control.dip_center_z;
                            if (r1 >= E->viscosity.z_weakzone_platebond && t1 < center_x)
                            {
                                tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
                                if (tempdist >= 0.0 - E->viscosity.left_weakzone_platebond && tempdist <= E->viscosity.right_weakzone_platebond)
                                {
                                    visc_weakzone = E->viscosity.visc_reduce_platebond;
                                }
                            }

                            break;

                        default:
                            break;
                        }
                    } /* end of trenchmigrate */
                }

                for (kk = 1; kk <= ends; kk++)
                {
                    temp0 += min(TT[kk], one) * E->N.vpt[GNVINDEX(kk, jj)];
                }
                for (m = 0; m < E->control.phasefile_C_flavor; m++)
                {
                    if (m == 0)
                    {
                        temp_visc_all = 0;
                    }
                    temp2 = E->viscosity.E[m] * (E->viscosity.T[m] - temp0);
                    temp3 = E->viscosity.N0[m] * exp(temp2);

                    temp_visc_all += log(temp3) * E->C_phasefile_element[m][i];
                    //fprintf(stderr, "%d %lf %lf %lf\n", m, E->C_phasefile_element[m][i], E->viscosity.N0[m], E->viscosity.E[m]);
                }

                EEta[(i - 1) * vpts + jj] = visc_weakzone * exp(temp_visc_all);
                //    fprintf(stderr, "visc CPU %d %d %.4e\n",E->parallel.me, (i - 1) * vpts + jj, EEta[(i - 1) * vpts + jj]);
                if (E->viscosity.visc_platebond_const && r1 > E->viscosity.z_weakzone_platebond)
                {
                    if (E->control.trechmigrate)
                    {
                        loc_mid = E->control.velo_surf_loc_mid;
                        loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
                        t1 = E->X[1][E->ien[i].node[jj]];
                        r1 = E->X[2][E->ien[i].node[jj]];
                        switch (E->control.initialTOption)
                        {
                        case 1:
                            tempdist = (t1 - loc_mid) * tan(E->control.dip_margin * M_PI / 180.0) + 1.0 - r1;
                            tempdist_left = (t1 - loc_mid) * tan(E->control.dip_margin_left * M_PI / 180.0) + 1.0 - r1;
                            if (tempdist_left >= (0.0 - E->viscosity.left_weakzone_platebond) *
                                                     tan(E->control.dip_margin_left * M_PI / 180.0) &&
                                tempdist <= E->viscosity.right_weakzone_platebond *
                                                tan(E->control.dip_margin * M_PI / 180.0) &&
                                r1 >= E->viscosity.z_weakzone_platebond)
                            {
                                EEta[(i - 1) * vpts + jj] = E->viscosity.visc_reduce_platebond;
                            }
                            break;
                        case 10:
                            center_x = loc_mid + E->control.dip_center_x;
                            center_z = 1.0 - E->control.dip_center_z;
                            if (r1 >= E->viscosity.z_weakzone_platebond && t1 < center_x)
                            {
                                tempdist = E->control.ocean_lith_margin_curve - sqrt((t1 - center_x) * (t1 - center_x) + (r1 - center_z) * (r1 - center_z));
                                if (tempdist >= 0.0 - E->viscosity.left_weakzone_platebond && tempdist <= E->viscosity.right_weakzone_platebond)
                                {
                                    EEta[(i - 1) * vpts + jj] = E->viscosity.visc_reduce_platebond;
                                }
                            }

                            break;

                        default:
                            break;
                        }
                    } /* end of trenchmigrate */
                }
            }
        }
        break;
    }

    fprintf(stderr, "visc finished\n");

    visits++;

    return;
}

void visc_from_S(E, Eta, EEta, propogate) struct All_variables *E;
double *Eta, *EEta;
int propogate;
{
    static int visits = 0;
    double one, two, scale, stress_magnitude, depth, exponent1, temp3;
    double *eedot;

    void strain_rate_2_inv();
    int e, l, z, jj, kk, m, i;
    double temp_visc_sdepv_expt, temp_visc_sdepv_trns, temp_visc_all;
    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->lmesh.nel;

    eedot = (double *)malloc((2 + nel) * sizeof(double));
    one = 1.0;
    two = 2.0;

    if (visits == 0)
    {
        for (e = 1; e <= nel; e++)
            eedot[e] = one;
    }
    else
        strain_rate_2_inv(E, eedot, 1);
    for (e = 1; e <= nel; e++)
    {
        eedot[e] = max(eedot[e], 1.0e-16);
    }
    for (e = 1; e <= nel; e++)
    {
        for (jj = 1; jj <= vpts; jj++)
        {
            for (m = 0; m < E->control.phasefile_C_flavor; m++)
            {
                if (m == 0)
                {
                    temp_visc_all = 0;
                }
                if (E->X[2][E->ien[e].node[jj]] < E->viscosity.zlm && E->viscosity.lower_diff)
                {
                    temp3 = EEta[(e - 1) * vpts + jj];
                }
                else
                {
                    exponent1 = one - one / E->viscosity.sdepv_expt[m];
                    scale = pow(two * eedot[e] / E->viscosity.sdepv_trns[m], exponent1);
                    temp3 = two * EEta[(e - 1) * vpts + jj] / (one + scale * pow(EEta[(e - 1) * vpts + jj], exponent1));
                }
                temp_visc_all += log(temp3) * E->C_phasefile_element[m][e];
            }
            EEta[(e - 1) * vpts + jj] = exp(temp_visc_all);
        }
    }
    visits++;
    free((void *)eedot);
    return;
}

void visc_from_P(E, Eta, EEta, propogate) struct All_variables *E;
double *Eta, *EEta;
int propogate;
{
    static int visits = 0;
    double one, two, scale, stress_magnitude, depth, exponent1;
    double *eedot;

    void strain_rate_2_inv();
    int e, l, z, jj, kk;

    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->mesh.nel;

    eedot = (double *)malloc((2 + nel) * sizeof(double));
    one = 1.0;
    two = 2.0;

    if (visits == 0)
    {
        for (e = 1; e <= nel; e++)
            eedot[e] = one;
    }
    else
        strain_rate_2_inv(E, eedot, 1);

    for (e = 1; e <= nel; e++)
    {

        for (jj = 1; jj <= vpts; jj++)
            EEta[(e - 1) * vpts + jj] = two * EEta[(e - 1) * vpts + jj] / (one + scale * pow(EEta[(e - 1) * vpts + jj], exponent1));
    }

    visits++;

    free((void *)eedot);
    return;
}

void strain_rate_2_inv_moreout(E, EEDOT, EEDOT11, EEDOT22, EEDOT12, SQRT) struct All_variables *E;
double *EEDOT, *EEDOT11, *EEDOT22, *EEDOT12;
int SQRT;
{
    void get_global_shape_fn();
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    double aaa, xk[3][5], edot[4][4], dudx[4][4];
    double VV[4][9];

    int e, i, p, q, n, k;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->lmesh.nno;
    const int vpts = vpoints[dims];

    const int nel = E->lmesh.nel;

    for (e = 1; e <= nel; e++)
    {

        for (i = 1; i <= ends; i++)
        {
            n = E->ien[e].node[i];
            VV[1][i] = E->V[1][n];
            VV[2][i] = E->V[2][n];
            if (dims == 3)
                VV[3][i] = E->V[3][n];
        }

        for (p = 1; p <= dims; p++)
            for (q = 1; q <= dims; q++)
                dudx[p][q] = 0.0;
        for (i = 1; i <= ends; i++)
            for (p = 1; p <= dims; p++)
                for (q = 1; q <= dims; q++)
                    dudx[p][q] += VV[p][i] * E->gNX[e].ppt[GNPXINDEX(q - 1, i, 1)];

        for (p = 1; p <= dims; p++)
            for (q = 1; q <= dims; q++)
                edot[p][q] = dudx[p][q] + dudx[q][p];

        if (dims == 2)
        {
            EEDOT[e] = edot[1][1] * edot[1][1] + edot[2][2] * edot[2][2] + edot[1][2] * edot[1][2] * 2.0;
            EEDOT11[e] = edot[1][1];
            EEDOT22[e] = edot[2][2];
            EEDOT12[e] = edot[1][2];
        }

        else if (dims == 3)
            EEDOT[e] = edot[1][1] * edot[1][1] + edot[1][2] * edot[1][2] * 2.0 + edot[2][2] * edot[2][2] + edot[2][3] * edot[2][3] * 2.0 + edot[3][3] * edot[3][3] + edot[1][3] * edot[1][3] * 2.0;
    }
    if (SQRT)
        for (e = 1; e <= nel; e++)
            EEDOT[e] = sqrt(0.5 * EEDOT[e]);
    else
        for (e = 1; e <= nel; e++)
            EEDOT[e] *= 0.5;

    return;
}

void strain_rate_2_inv(E, EEDOT, SQRT) struct All_variables *E;
double *EEDOT;
int SQRT;
{
    void get_global_shape_fn();
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    double aaa, xk[3][5], edot[4][4], dudx[4][4];
    double VV[4][9];

    int e, i, p, q, n, nel, k;
    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];

    nel = E->mesh.nel;

    for (e = 1; e <= nel; e++)
    {

        for (i = 1; i <= ends; i++)
        {
            n = E->ien[e].node[i];
            VV[1][i] = E->V[1][n];
            VV[2][i] = E->V[2][n];
            if (dims == 3)
                VV[3][i] = E->V[3][n];
        }
        for (p = 1; p <= dims; p++)
            for (q = 1; q <= dims; q++)
                dudx[p][q] = 0.0;

        for (i = 1; i <= ends; i++)
            for (p = 1; p <= dims; p++)
                for (q = 1; q <= dims; q++)
                    dudx[p][q] += VV[p][i] * E->gNX[e].ppt[GNPXINDEX(q - 1, i, 1)];

        for (p = 1; p <= dims; p++)
            for (q = 1; q <= dims; q++)
                edot[p][q] = dudx[p][q] + dudx[q][p];

        if (dims == 2)
            EEDOT[e] = edot[1][1] * edot[1][1] + edot[2][2] * edot[2][2] + edot[1][2] * edot[1][2] * 2.0;
        else if (dims == 3)
            EEDOT[e] = edot[1][1] * edot[1][1] + edot[1][2] * edot[1][2] * 2.0 + edot[2][2] * edot[2][2] + edot[2][3] * edot[2][3] * 2.0 + edot[3][3] * edot[3][3] + edot[1][3] * edot[1][3] * 2.0;
    }

    if (SQRT)
        for (e = 1; e <= nel; e++)
            EEDOT[e] = sqrt(0.5 * EEDOT[e]);
    else
        for (e = 1; e <= nel; e++)
            EEDOT[e] *= 0.5;

    return;
}

int layers(E, x2) struct All_variables *E;
double x2;
{
    int llayers = 0;

    if (x2 >= E->viscosity.zlith)
    {
        llayers = 1;
    }
    else if (x2 >= E->viscosity.z300)
    {
        llayers = 2;
    }
    else if (x2 >= E->viscosity.z410)
    {
        llayers = 3;
    }
    else if (x2 >= E->viscosity.zlm)
    {
        llayers = 4;
    }
    else if (x2 >= E->viscosity.z1000)
    {
        llayers = 5;
    }
    else
    {
        llayers = 6;
    }

    return (llayers);
}

static void visc_from_B(struct All_variables *E, double *Eta, double *EEta, int propogate)
{
    static int visited = 0;
    double scale, stress_magnitude, depth, exponent1, eta_old, eta_old2, eta_new;
    double *eedot;
    double zzz, zz[9];
    double tau, tau2, ettby, ettnew;
    int m, l, z, jj, kk, i;
    static double ndz_to_m;
#ifdef DEBUG
    FILE *out;
#endif
    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    const int ends = enodes[E->mesh.nsd];

    eedot = (double *)malloc((2 + nel) * sizeof(double));

#ifdef DEBUG
    out = fopen("tmp.visc", "w");
#endif
    if (!visited)
    {
        /* 
       scaling from nod-dim radius (0...1) to m 

       (only used for dimensional version)

    */
        ndz_to_m = E->monitor.length_scale;

        /*  */
        if (E->parallel.me == 0)
        { /* control output */
            for (l = 1; l <= E->viscosity.num_mat; l++)
            {
                fprintf(stderr, "Plasticity: %d/%d: a=%g b=%g p=%g\n",
                        l, E->viscosity.num_mat,
                        E->viscosity.abyerlee[l - 1],
                        E->viscosity.bbyerlee[l - 1],
                        E->viscosity.lbyerlee[l - 1]);
            }
            fprintf(stderr, "\tdim: %i trans: %i offset: %g\n",
                    E->viscosity.plasticity_dimensional,
                    E->viscosity.plasticity_trans,
                    E->viscosity.plasticity_viscosity_offset);
            fprintf(stderr, "\tpsrw: %i\n", E->viscosity.psrw);
        }
        /* 
       get strain rates for all elements 
    */
        for (i = 1; i <= nel; i++)
            eedot[i] = 1.0;
    }
    else
    {
        if (E->viscosity.psrw)
            strain_rate_2_inv(E, eedot, 0);
        else
            strain_rate_2_inv(E, eedot, 1);
    }
    if (E->viscosity.psrw)
    {
        /* strain-rate weakening */
        for (i = 1; i <= nel; i++)
        {
            l = E->mat[i] - 1; /* material of element */
            for (kk = 1; kk <= ends; kk++)
            { /* loop through integration points*/

                zz[kk] = (1.0 - E->X[3][E->ien[i].node[kk]]);
                if (E->viscosity.plasticity_dimensional)
                    zz[kk] *= ndz_to_m; /* scale to meters */
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                zzz = 0.0;
                for (kk = 1; kk <= ends; kk++)
                {
                    zzz += zz[kk] * E->N.vpt[GNVINDEX(kk, jj)];
                }
                if (E->viscosity.plasticity_dimensional)
                {
                    tau = (E->viscosity.abyerlee[l] * zzz +
                           E->viscosity.bbyerlee[l]) *
                          E->viscosity.lbyerlee[l];
                    tau /= E->monitor.tau_scale;
                }
                else
                {
                    tau = E->viscosity.abyerlee[l] * zzz + E->viscosity.bbyerlee[l];
                    tau = min(tau, E->viscosity.lbyerlee[l]);
                }
                if ((visited > 1) && (tau < 1e15))
                {
                    tau2 = tau * tau;
                    eta_old = EEta[(i - 1) * vpts + jj];
                    eta_old2 = eta_old * eta_old;
                    eta_new = (tau2 * eta_old) / (tau2 + 2.0 * eta_old2 * eedot[i]);
                    EEta[(i - 1) * vpts + jj] = ettnew;
                    //if(E->parallel.me==0)fprintf(stderr,"tau: %11g eII: %11g eta_old: %11g eta_new: %11g\n",tau, eedot[i],eta_old,eta_new);
                }
            }
        }
    }
    else
    {

        /* regular plasticity */

        for (i = 1; i <= nel; i++)
        {
            /* 
	 loop through all elements 
      */
            l = E->mat[i] - 1; /* material of element */
            for (kk = 1; kk <= ends; kk++)
            { /* loop through integration points*/
                /* depth in meters */

                zz[kk] = (1.0 - E->X[3][E->ien[i].node[kk]]);
                if (E->viscosity.plasticity_dimensional)
                    zz[kk] *= ndz_to_m; /* scale to meters */
            }
            for (jj = 1; jj <= vpts; jj++)
            {
                /* loop over nodes in element */
                zzz = 0.0;
                for (kk = 1; kk <= ends; kk++)
                {
                    /* 
	     depth  [m] 
	  */
                    zzz += zz[kk] * E->N.vpt[GNVINDEX(kk, jj)];
                }
                if (E->viscosity.plasticity_dimensional)
                {
                    /* byerlee type */

                    /* 
	     yield stress in [Pa] 
	  */
                    tau = (E->viscosity.abyerlee[l] * zzz + E->viscosity.bbyerlee[l]) * E->viscosity.lbyerlee[l];
                    /* 
	     scaled stress 
	  */
                    tau /= E->monitor.tau_scale;
                }
                else
                {

                    tau = E->viscosity.abyerlee[l] * zzz + E->viscosity.bbyerlee[l];

                    tau = min(tau, E->viscosity.lbyerlee[l]);
                }
                /* 
	   
	`byerlee viscosity' : tau = 2 eps eta, this is non-dim
	plus some offset as in Stein et al. 
	
	*/
                ettby = tau / (2.0 * (eedot[i] + 1e-7)) + E->viscosity.plasticity_viscosity_offset;
                /* 
	   
	decide on the plasticity transition
	
	
	*/
                if (E->viscosity.plasticity_trans)
                {
                    /* 
	     eta = 1./(1./eta(k)+1./ettby)  
	  */

                    ettnew = 1.0 / (1.0 / EEta[(i - 1) * vpts + jj] + 1.0 / ettby);
                    //fprintf(stderr,"a: %g %g %g\n",EEta[ (i-1)*vpts + jj ],ettby,ettnew);
                }
                else
                {
                    /* 
	     min(\eta_p, \eta_visc )
	  */
                    ettnew = min(EEta[(i - 1) * vpts + jj], ettby);
                    //fprintf(stderr,"m: %g %g %g\n",EEta[ (i-1)*vpts + jj ],ettby,ettnew);
                }
#ifdef DEBUG
                /* output format 
	   
	z[km] tau[MPa] eps[s^-1] eta_b[Pas] eta_T[Pas] eta_c[Pas]
	
	*/
                if (visited)
                    fprintf(out, "%10.2f %17.4e %17.4e %17.4e %17.4e %17.4e\n",
                            zzz / 1e3, tau * E->monitor.tau_scale / 1e6,
                            eedot[i] / E->monitor.time_scale,
                            ettby * E->data.ref_viscosity,
                            EEta[(i - 1) * vpts + jj] * E->data.ref_viscosity,
                            ettnew * E->data.ref_viscosity);
#endif
                //      if(visited)
                //	fprintf(stderr,"%11g %11g %11g %11g\n",ettnew,EEta[ (i-1)*vpts + jj ] ,ettby,eedot[i]);
                EEta[(i - 1) * vpts + jj] = ettnew;
            }
        } /* end regular plasticity */
    }
#ifdef DEBUG
    fclose(out);
#endif
    visited = 1;
    free((void *)eedot);
    return;
}
