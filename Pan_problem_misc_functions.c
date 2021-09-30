
#include "element_definitions.h"
#include "global_defs.h"

#include <stdio.h>

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <unistd.h>

#if (defined __sunos__)
#include <string.h>
#else
#include <strings.h>
#endif

#if defined(__sgi) || defined(__osf__)
#include <sys/types.h>
#endif

int get_process_identifier()
{
    int pid;

    pid = (int)getpid();
    return (pid);
}

void unique_copy_file(E, name, comment) struct All_variables *E;
char *name, *comment;
{
    char unique_name[500];
    char command[600];

    sprintf(unique_name, "%06d.%s-%s", E->control.PID, comment, name);
    sprintf(command, "cp -f %s %s\n", name, unique_name);
    system(command);
}

void thermal_buoyancy(E) struct All_variables *E;

{
    int i, m, n, mark, num;
    double coord_z, coord_z_km, loc_mid, t1, r1, tempdist;
    double *H, slope, temp1;
    void remove_horiz_ave();
    void phase_change();
    void phase_change_basalt();
    void get_phase_buoyancy();
    slope = (E->data.therm_exp_factor - 1.0);

    H = (double *)malloc((E->mesh.noz + 1) * sizeof(double));

    if (abs(E->control.Ra_670) > 0.0 || abs(E->control.Ra_410) > 0.0)
    {

        if (E->control.Ra_670_basalt > 0.0)
        {
            phase_change_basalt(E, E->Fas670, E->Fas670_b, E->Fas410, E->Fas410_b, E->Fas670_basalt, E->Fas670_basalt_b, E->Fas670_all_b);
        }
        else
        {
            phase_change(E, E->Fas670, E->Fas670_b, E->Fas410, E->Fas410_b);
        }
        for (i = 1; i <= E->mesh.nno; i++)
            E->buoyancy[i] = -E->control.Ra_670 * E->Fas670[i] - E->control.Ra_410 * E->Fas410[i];
    }
    else
        for (i = 1; i <= E->mesh.nno; i++)
            E->buoyancy[i] = 0.0;

    if (E->control.phasefile_buoyancy)
    {
        get_phase_buoyancy(E);
        fprintf(stderr, "buoyancy phasefile finished\n");
        if (E->control.phasefile_buoyancy_depth > 0)
        {
            /* first get transition depth */
            for (n = 1; n <= E->mesh.noz; n++)
            {
                coord_z = E->X[2][n];
                coord_z_km = (1.0 - coord_z) * E->data.layer_km;
                if (coord_z_km >= E->control.phasefile_buoyancy_depth)
                {
                    mark = n;
                    break;
                }
            }
            for (m = 1; m <= E->mesh.nox; m++)
            {
                for (n = mark; n <= E->mesh.noz; n++)
                {
                    num = n + E->mesh.noz * (m - 1);
                    E->buoyancy[num] = E->control.Ra_temp * E->T[num] * (slope * E->X[2][num] + 1.0) - E->control.Ra_comp * E->C[num];
                }
            }
        }
        /* considering continent buoyancy effects */
        if (E->control.phasefile_buoyancy_continent > 0.0)
        {
            if (E->control.trechmigrate)
            {
                loc_mid = E->control.velo_surf_loc_mid;
                loc_mid += E->control.velo_surf_loc_mid_rate * E->monitor.elapsed_time * E->control.timescale;
                for (n = 1; n <= E->mesh.noz; n++)
                {
                    coord_z = E->X[2][n];
                    if (coord_z >= E->viscosity.zlith)
                    {
                        mark = n;
                        break;
                    }
                }

                for (m = 1; m <= E->mesh.nox; m++)
                {
                    for (n = mark; n <= E->mesh.noz; n++)
                    {
                        num = n + E->mesh.noz * (m - 1);
                        t1 = E->X[1][num];
                        r1 = E->X[2][num];
                        tempdist = (t1 - loc_mid) * tan(E->control.dip_margin * 3.14159265 / 180.0) + 1.0 - r1;
                        if (tempdist <= (E->viscosity.right_weakzone_platebond) * tan(E->control.dip_margin * 3.14159265 / 180.0) && r1 >= E->viscosity.zlith)
                        {
                            E->buoyancy[num] += E->control.phasefile_buoyancy_continent * E->control.Ra_temp / (E->data.therm_exp * E->data.ref_temperature * E->data.density);
                        }
                    }
                }
            }
        }
        fprintf(stderr, "continent phasefile finished\n");

        /* end continent correction */
        /* correction for crust denstiy */
        if (E->control.phasefile_buoyancy_crust > 0.0)
        {
            for (n = 1; n <= E->mesh.noz; n++)
            {
                coord_z = E->X[2][n];
                coord_z_km = (1.0 - coord_z) * E->data.layer_km;
                if (coord_z_km >= E->control.phasefile_buoyancy_crust_depth)
                {
                    mark = n;
                    break;
                }
            }

            for (m = 1; m <= E->mesh.nox; m++)
            {
                for (n = mark; n <= E->mesh.noz; n++)
                {
                    num = n + E->mesh.noz * (m - 1);
                    coord_z = E->X[2][n];
                    coord_z_km = (1.0 - coord_z) * E->data.layer_km;
                    if (coord_z_km >= E->control.phasefile_buoyancy_crust_depth)
                    {
                        if (E->C[num] > 0.0)
                        {
                            E->buoyancy[num] += E->C[num] * E->control.phasefile_buoyancy_crust * E->control.Ra_temp / (E->data.therm_exp * E->data.ref_temperature * E->data.density);
                        }
                    }
                }
            }
        }
        fprintf(stderr, "crust phasefile finished\n");

        /* end crust correction */
    }
    else
    {

        for (i = 1; i <= E->mesh.nno; i++)
        {

            /*  for constant density planet with g decreasing linearly with depth */
            /*        E->buoyancy[i] = E->control.Ra_temp * E->T[i] * E->X[2][i]
                       + E->control.Ra_comp * E->C[i] * E->X[2][i];
*/
            /*  for constant g in the mantle */
            /*  first option is using density from phase diagram input files, the other is traditional way using Ra */
            E->buoyancy[i] += E->control.Ra_temp * E->T[i] * (slope * E->X[2][i] + 1.0) - E->control.Ra_comp * E->C[i];
        }
    }

    /*    for(i=1;i<=E->mesh.nno;i++)   {
        fprintf(stderr,"buoyancy %d %lf\n",i,E->buoyancy[i]);

    }
*/
    remove_horiz_ave(E, E->buoyancy, H, 0);

    fprintf(stderr, "buoyancy phasefile finished\n");

    free((void *)H);
    return;
}

double SIN_D(x) double x;
{
#if defined(__osf__)
    return sind(x);
#else
    return sin((x / 180.0) * M_PI);
#endif
}

double COT_D(x) double x;
{
#if defined(__osf__)
    return cotd(x);
#else
    return tan(((90.0 - x) / 180.0) * M_PI);
#endif
}

/* non-runaway malloc */

void *Malloc1(bytes, file, line) int bytes;
char *file;
int line;
{
    void *ptr;

    ptr = malloc((size_t)bytes);
    if (ptr == (void *)NULL)
    {
        fprintf(stderr, "Memory: cannot allocate another %d bytes \n(line %d of file %s)\n", bytes, line, file);
        exit(0);
    }

    return (ptr);
}

/* Copy one field to another on a different (but similarly structured) mesh.  The
   field TT (output) is on the standard mesh for this problem, T (input) is on the different
   mesh whose coordinates are described by X,Z,Y 


   */

void fcopy_interpolating(E, X, Z, Y, nx, nz, ny, T, TT) struct All_variables *E;
double *X, *Z, *Y;
int nx, nz, ny;
double *T, *TT;
{
    double cross2d();
    void p_to_nodes();
    void p_to_centres();
    void vcopy();
    double CPU_time0(), time;

    double *P;

    int i, j, found, not_found;
    int elX, elY, elZ, ex, ez;
    int old_ex, old_ez, old_ey;
    int node1, node2, node3, node4;
    int node5, node6, node7, node8;
    double inside1, inside2, inside3, inside4;
    double inside5, inside6, inside7, inside8, inside9, inside10, inside11, inside12;
    double distance1, distance2, distance3, distance4;
    double distance5, distance6, distance7, distance8;
    double d1, d2, d3, d4, d5, d6, d7, d8;

    const int dims = E->mesh.nsd;

    P = (double *)malloc((E->mesh.nno + 1) * sizeof(double));

    /* run over all the data points (take care to hit only one element), determine
       inside/outside of each element. The criterion for inside-ness is that the
       cross product of the vectors joining the point to the ends of each edge should
       have the same sign as the dot product of the centre to the ends of each edge.
       There are, undoubtedly, better ways to do this !!!

       Because a CITCOM node-ordering is assumed, we can also guess that the best place to start
       looking for a node is where you found the last one !

       */

    old_ex = old_ez = old_ey = 1;

    elX = nx - 1;
    elZ = nz - 1;

    if (E->control.print_convergence)
    {
        time = CPU_time0();
        fprintf(stderr, "Interpolating ...");
    }

    not_found = 0;

    if (2 == dims)
        for (i = 1; i <= E->mesh.nno; i++)
        {
            found = 0;
            for (ex = old_ex; ex <= elX && found == 0; ex++)
                for (ez = 1; ez <= elZ && found == 0; ez++)
                {
                    node1 = ez + (ex - 1) * nz;
                    node2 = node1 + offset[2].vector[1] + offset[2].vector[0] * nz;
                    node3 = node1 + offset[3].vector[1] + offset[3].vector[0] * nz;
                    node4 = node1 + offset[4].vector[1] + offset[4].vector[0] * nz;

                    if ((inside1 = cross2d(X[node1] - E->X[1][i], Z[node1] - E->X[2][i], X[node2] - E->X[1][i], Z[node2] - E->X[2][i], 3)) <= 0.0 &&
                        (inside4 = cross2d(X[node4] - E->X[1][i], Z[node4] - E->X[2][i], X[node1] - E->X[1][i], Z[node1] - E->X[2][i], 3)) <= 0.0 &&
                        (inside2 = cross2d(X[node2] - E->X[1][i], Z[node2] - E->X[2][i], X[node3] - E->X[1][i], Z[node3] - E->X[2][i], 3)) <= 0.0 &&
                        (inside3 = cross2d(X[node3] - E->X[1][i], Z[node3] - E->X[2][i], X[node4] - E->X[1][i], Z[node4] - E->X[2][i], 3)) <= 0.0)
                    {
                        found = node1;
                        old_ex = ex;
                    }
                }

            /* finish the loop if not found */
            for (ex = 1; ex <= old_ex && found == 0; ex++)
                for (ez = 1; ez <= elZ && found == 0; ez++)
                {
                    node1 = ez + (ex - 1) * nz;
                    node2 = node1 + offset[2].vector[1] + offset[2].vector[0] * nz;
                    node3 = node1 + offset[3].vector[1] + offset[3].vector[0] * nz;
                    node4 = node1 + offset[4].vector[1] + offset[4].vector[0] * nz;

                    if ((inside1 = cross2d(X[node1] - E->X[1][i], Z[node1] - E->X[2][i], X[node2] - E->X[1][i], Z[node2] - E->X[2][i], 3)) <= 0.0 &&
                        (inside4 = cross2d(X[node4] - E->X[1][i], Z[node4] - E->X[2][i], X[node1] - E->X[1][i], Z[node1] - E->X[2][i], 3)) <= 0.0 &&
                        (inside2 = cross2d(X[node2] - E->X[1][i], Z[node2] - E->X[2][i], X[node3] - E->X[1][i], Z[node3] - E->X[2][i], 3)) <= 0.0 &&
                        (inside3 = cross2d(X[node3] - E->X[1][i], Z[node3] - E->X[2][i], X[node4] - E->X[1][i], Z[node4] - E->X[2][i], 3)) <= 0.0)
                    {
                        found = node1;
                        old_ex = ex;
                    }
                }

            /* and having found the right node location, interpolate the appropriate value to it */

            if (!found)
                not_found++;
            else
            {
                distance1 = ((X[node1] - E->X[1][i]) * (X[node1] - E->X[1][i]) + (Z[node1] - E->X[2][i]) * (Z[node1] - E->X[2][i]));
                distance2 = ((X[node2] - E->X[1][i]) * (X[node2] - E->X[1][i]) + (Z[node2] - E->X[2][i]) * (Z[node2] - E->X[2][i]));
                distance3 = ((X[node3] - E->X[1][i]) * (X[node3] - E->X[1][i]) + (Z[node3] - E->X[2][i]) * (Z[node3] - E->X[2][i]));
                distance4 = ((X[node4] - E->X[1][i]) * (X[node4] - E->X[1][i]) + (Z[node4] - E->X[2][i]) * (Z[node4] - E->X[2][i]));

                d1 = distance2 * distance3 * distance4;
                d2 = distance1 * distance3 * distance4;
                d3 = distance2 * distance1 * distance4;
                d4 = distance2 * distance3 * distance1;

                TT[i] = (d1 * T[node1] + d2 * T[node2] + d3 * T[node3] + d4 * T[node4]) / (d1 + d2 + d3 + d4);
            }
        }

    else
    {
        elY = (3 == dims) ? ny - 1 : 1;
        fprintf(stderr, "3D interolator not yet implemented !!! \n");
        vcopy(TT, T, 1, E->mesh.nno);
    }

    if (E->control.print_convergence)
        fprintf(stderr, ". done (%lf secs)\n", CPU_time0() - time);

    if (not_found)
        fprintf(E->fp, "Warning: unable to interpolate old  data to %d nodes in the new mesh\n", not_found);
    else
    {
        p_to_centres(E, TT, P, E->mesh.levmax); /* if interpolated, apply slight smoothing */
        p_to_nodes(E, P, TT, E->mesh.levmax);
    }

    free((void *)P);

    return;
}

/* returns the out of plane component of the cross product of
   the two vectors assuming that one is looking AGAINST the 
   direction of the axis of D, anti-clockwise angles
   are positive (are you sure ?), and the axes are ordered 2,3 or 1,3 or 1,2 */

double cross2d(x11, x12, x21, x22, D) double x11, x12, x21, x22;
int D;
{
    if (1 == D)
        return (x11 * x22 - x12 * x21);
    if (2 == D)
        return (-x11 * x22 + x12 * x21);
    if (3 == D)
        return (x11 * x22 - x12 * x21);
}

void field_arbitrary_rectangle_file(E, parse_and_apply, RECT, name, field, BC, bcbitf, bcmask_on, bcmask_off) struct All_variables *E;
int parse_and_apply;
struct Rect *RECT;
char *name;
double *field;
int BC;
unsigned int *bcbitf;
unsigned int bcmask_on, bcmask_off;
{

    char read_string[500];
    double weight, radius2, weight2;
    double x1, y1, z1;
    int in1, in2, in3;
    int number, node;
    int combine_option;

    void field_arbitrary_rectangle();

    sprintf(read_string, "%s_rect", name);
    input_int(read_string, &(RECT->numb), "0");
    sprintf(read_string, "%s_rectx1", name);
    input_double_vector(read_string, RECT->numb, RECT->x1);
    sprintf(read_string, "%s_rectx2", name);
    input_double_vector(read_string, RECT->numb, RECT->x2);
    sprintf(read_string, "%s_rectz1", name);
    input_double_vector(read_string, RECT->numb, RECT->z1);
    sprintf(read_string, "%s_rectz2", name);
    input_double_vector(read_string, RECT->numb, RECT->z2);
    sprintf(read_string, "%s_recty1", name);
    input_double_vector(read_string, RECT->numb, RECT->y1);
    sprintf(read_string, "%s_recty2", name);
    input_double_vector(read_string, RECT->numb, RECT->y2);
    sprintf(read_string, "%s_recthw", name);
    input_double_vector(read_string, RECT->numb, RECT->halfw);
    sprintf(read_string, "%s_rectmag", name);
    input_double_vector(read_string, RECT->numb, RECT->mag);
    sprintf(read_string, "%s_rectovl", name);
    input_char_vector(read_string, RECT->numb, RECT->overlay);

    if (parse_and_apply)
        field_arbitrary_rectangle(E, RECT, field, BC, bcbitf, bcmask_on, bcmask_off);

    return;
}

void field_arbitrary_rectangle(E, RECT, field, BC, bcbitf, bcmask_on, bcmask_off) struct All_variables *E;
struct Rect *RECT;
double *field;
int BC;
unsigned int *bcbitf;
unsigned int bcmask_on, bcmask_off;
{

    double weight, radius2, weight2;
    double x1, y1, z1;
    int in1, in2, in3;
    int number, node;
    int combine_option;

    for (node = 1; node <= E->mesh.nno; node++)
    {
        x1 = E->X[1][node];
        z1 = E->X[2][node];
        y1 = (E->mesh.nsd != 3) ? 0.0 : E->X[3][node];

        for (number = 0; number < RECT->numb; number++)
        {

            switch (RECT->overlay[number])
            {
            case 'M':
                weight = 1.0;
                break;
            case 'R':
            case 'A':
                weight = 0.0;
                break;
            }

            in1 = (x1 >= RECT->x1[number] && x1 <= RECT->x2[number]);
            in2 = (z1 >= RECT->z1[number] && z1 <= RECT->z2[number]);
            in3 = (3 != E->mesh.nsd || y1 >= RECT->y1[number] && y1 <= RECT->y2[number]);

            if (in1 && in2 && in3)
            {
                weight = RECT->mag[number];
                radius2 = 0.0;
            }
            else
            {
                radius2 =
                    (in1 ? 0.0 : min(1.0e-10 + (x1 - RECT->x1[number]) * (x1 - RECT->x1[number]), 1.0e-10 + (x1 - RECT->x2[number]) * (x1 - RECT->x2[number]))) +
                    (in2 ? 0.0 : min(1.0e-10 + (z1 - RECT->z1[number]) * (z1 - RECT->z1[number]), 1.0e-10 + (z1 - RECT->z2[number]) * (z1 - RECT->z2[number]))) +
                    (in3 ? 0.0 : min(1.0e-10 + (y1 - RECT->y1[number]) * (y1 - RECT->y1[number]), 1.0e-10 + (y1 - RECT->y2[number]) * (y1 - RECT->y2[number])));

                weight += RECT->mag[number] * exp(-radius2 / (0.5 * RECT->halfw[number] * RECT->halfw[number]));
            }

            switch (RECT->overlay[number])
            {
            case 'R':
                if (radius2 > RECT->halfw[number] * RECT->halfw[number])
                    break;
                weight2 = 1.0 / (1.0 + radius2 / (1.0e-10 + RECT->halfw[number] * RECT->halfw[number]));
                field[node] = (1.0 - weight2) * field[node] + weight2 * weight;
                /*fprintf(stderr," %d (%lf/%lf) = %lf,%lf   /   %lf (%lf)\n",node,x1,z1,weight,weight2,field[node],radius2);*/
                break;
            case 'M':
                field[node] *= weight;
                break;
            case 'A':
                field[node] += weight;
                break;
            default:
                fprintf(E->fp, "RECTANGLE: %d can't work out how to combine new/old fields\n", number);
                break;
            }

            if (BC)
            {
                bcbitf[node] = (bcbitf[node] | bcmask_on);
                bcbitf[node] = (bcbitf[node] & (~bcmask_off));
            }
        }
    }

    return;
}

void field_arbitrary_circle_file(E, parse_and_apply, CIRC, name, field, BC, bcbitf, bcmask_on, bcmask_off) struct All_variables *E;
int parse_and_apply;
struct Circ *CIRC;
char *name;
double *field;
int BC;
unsigned int *bcbitf;
unsigned int bcmask_on, bcmask_off;
{

    char read_string[500];
    double weight, radius2, weight2;
    double x1, y1, z1;
    int in1;
    int number, node;
    int combine_option;
    void field_arbitrary_circle();

    sprintf(read_string, "%s_circ", name);
    input_int(read_string, &(CIRC->numb), "0");
    sprintf(read_string, "%s_circx", name);
    input_double_vector(read_string, CIRC->numb, CIRC->x);
    sprintf(read_string, "%s_circz", name);
    input_double_vector(read_string, CIRC->numb, CIRC->z);
    sprintf(read_string, "%s_circy", name);
    input_double_vector(read_string, CIRC->numb, CIRC->y);
    sprintf(read_string, "%s_circrad", name);
    input_double_vector(read_string, CIRC->numb, CIRC->rad);
    sprintf(read_string, "%s_circmag", name);
    input_double_vector(read_string, CIRC->numb, CIRC->mag);
    sprintf(read_string, "%s_circhw", name);
    input_double_vector(read_string, CIRC->numb, CIRC->halfw);
    sprintf(read_string, "%s_circovl", name);
    input_char_vector(read_string, CIRC->numb, CIRC->overlay);

    if (parse_and_apply)
        field_arbitrary_circle(E, CIRC, field, BC, bcbitf, bcmask_on, bcmask_off);

    return;
}
void field_arbitrary_circle(E, CIRC, field, BC, bcbitf, bcmask_on, bcmask_off) struct All_variables *E;
struct Circ *CIRC;
double *field;
int BC;
unsigned int *bcbitf;
unsigned int bcmask_on, bcmask_off;
{

    char read_string[500];
    double weight, radius2, weight2;
    double x1, y1, z1;
    int in1;
    int number, node;
    int combine_option;

    for (node = 1; node <= E->mesh.nno; node++)
    {
        x1 = E->X[1][node];
        z1 = E->X[2][node];
        y1 = (E->mesh.nsd != 3) ? 0.0 : E->X[3][node];

        for (number = 0; number < CIRC->numb; number++)
        {
            switch (CIRC->overlay[number])
            {
            case 'M':
                weight = 1.0;
                break;
            case 'R':
            case 'A':
                weight = 0.0;
                break;
            }

            radius2 =
                (x1 - CIRC->x[number]) * (x1 - CIRC->x[number]) +
                (z1 - CIRC->z[number]) * (z1 - CIRC->z[number]) +
                ((E->mesh.nsd != 3) ? 0.0 : (y1 - CIRC->y[number]) * (y1 - CIRC->y[number]));

            if (radius2 <= CIRC->rad[number] * CIRC->rad[number])
            {
                weight = CIRC->mag[number];
                radius2 = 0.0;
            }
            else
            {
                radius2 -= CIRC->rad[number] * CIRC->rad[number];
                weight += CIRC->mag[number] * exp(-2.0 * radius2 / (1.0e-10 + CIRC->halfw[number] * CIRC->halfw[number]));
            }

            switch (CIRC->overlay[number])
            {
            case 'R':
                if (radius2 > CIRC->halfw[number] * CIRC->halfw[number])
                    break;
                weight2 = 1.0 / (1.0 + radius2 / (1.0e-10 + CIRC->halfw[number] * CIRC->halfw[number]));
                field[node] = (1.0 - weight2) * field[node] + weight2 * weight;

                /*fprintf(stderr," %d (%lf/%lf) = %lf,%lf   /   %lf (%lf)\n",node,x1,z1,weight,weight2,field[node],radius2);*/
                break;
            case 'M':
                field[node] *= weight;
                break;
            case 'A':
                field[node] += weight;
                break;
            default:
                fprintf(E->fp, "CIRCLE: %d can't work out how to combine new/old fields\n", number);
                break;
            }

            if (BC)
            {
                bcbitf[node] = (bcbitf[node] | bcmask_on);
                bcbitf[node] = (bcbitf[node] & (~bcmask_off));
            }
        }
    }

    return;
}

void field_arbitrary_harmonic_file(E, parse_and_apply, HARM, name, field, BC, bcbitf, bcmask_on, bcmask_off) struct All_variables *E;
int parse_and_apply;
struct Harm *HARM;
char *name;
int BC;
double *field;
unsigned int *bcbitf;
unsigned int bcmask_on, bcmask_off;
{
    char read_string[500];
    void field_arbitrary_harmonic();
    int i;

    sprintf(read_string, "%s_harm", name);
    input_int(read_string, &(HARM->numb), "0");
    sprintf(read_string, "%s_harms", name);
    input_int(read_string, &(HARM->harms), "0,0,19");
    sprintf(read_string, "%s_harmoff", name);
    input_double_vector(read_string, HARM->numb, HARM->off);
    sprintf(read_string, "%s_harmx1", name);
    input_double_vector(read_string, HARM->numb, HARM->x1);
    sprintf(read_string, "%s_harmx2", name);
    input_double_vector(read_string, HARM->numb, HARM->x2);
    sprintf(read_string, "%s_harmz1", name);
    input_double_vector(read_string, HARM->numb, HARM->z1);
    sprintf(read_string, "%s_harmz2", name);
    input_double_vector(read_string, HARM->numb, HARM->z2);
    sprintf(read_string, "%s_harmy1", name);
    input_double_vector(read_string, HARM->numb, HARM->z1);
    sprintf(read_string, "%s_harmy2", name);
    input_double_vector(read_string, HARM->numb, HARM->z2);
    sprintf(read_string, "%s_harmovl", name);
    input_char_vector(read_string, HARM->numb, HARM->overlay);

    for (i = 0; i < HARM->harms; i++)
    {
        sprintf(read_string, "%s_harmkx%02d", name, i + 1);
        input_double_vector(read_string, HARM->numb, HARM->kx[i]);
        sprintf(read_string, "%s_harmkz%02d", name, i + 1);
        input_double_vector(read_string, HARM->numb, HARM->kz[i]);
        sprintf(read_string, "%s_harmky%02d", name, i + 1);
        input_double_vector(read_string, HARM->numb, HARM->ky[i]);
        sprintf(read_string, "%s_harmka%02d", name, i + 1);
        input_double_vector(read_string, HARM->numb, HARM->ka[i]);
        sprintf(read_string, "%s_harmphx%02d", name, i + 1);
        input_double_vector(read_string, HARM->numb, HARM->phx[i]);
        sprintf(read_string, "%s_harmphz%02d", name, i + 1);
        input_double_vector(read_string, HARM->numb, HARM->phz[i]);
        sprintf(read_string, "%s_harmphy%02d", name, i + 1);
        input_double_vector(read_string, HARM->numb, HARM->phy[i]);
    }

    if (parse_and_apply)
        field_arbitrary_harmonic(E, HARM, field, BC, bcbitf, bcmask_on, bcmask_off);

    return;
}

void field_arbitrary_harmonic(E, HARM, field, BC, bcbitf, bcmask_on, bcmask_off) struct All_variables *E;
struct Harm *HARM;
double *field;
int BC;
unsigned int *bcbitf;
unsigned int bcmask_on, bcmask_off;
{

    double weight, radius2, weight2;
    double x1, y1, z1;
    int in1, in2, in3;
    int number, node, l;
    int combine_option;

    for (node = 1; node <= E->mesh.nno; node++)
    {
        x1 = E->X[1][node];
        z1 = E->X[2][node];
        y1 = (E->mesh.nsd != 3) ? 0.0 : E->X[3][node];

        for (number = 0; number < HARM->numb; number++)
        {

            switch (HARM->overlay[number])
            {
            case 'M':
                weight = 1.0;
                break;
            case 'R':
            case 'A':
                weight = 0.0;
                break;
            }

            in1 = (x1 >= HARM->x1[number] && x1 <= HARM->x2[number]);
            in2 = (z1 >= HARM->z1[number] && z1 <= HARM->z2[number]);
            in3 = (3 != E->mesh.nsd || y1 >= HARM->y1[number] && y1 <= HARM->y2[number]);

            if (in1 && in2 && in3)
            {
                weight = HARM->off[number];
                for (l = 0; l < HARM->harms; l++)
                {
                    weight += HARM->ka[l][number] *
                              cos((HARM->kx[l][number] * x1 + HARM->phx[l][number]) * M_PI) *
                              cos((HARM->kz[l][number] * z1 + HARM->phz[l][number]) * M_PI) *
                              cos((HARM->ky[l][number] * y1 + HARM->phy[l][number]) * M_PI);
                }

                switch (HARM->overlay[number])
                {
                case 'R':
                    field[node] = weight;
                    break;
                case 'M':
                    field[node] *= weight;
                    break;
                case 'A':
                    field[node] += weight;
                    break;
                default:
                    fprintf(E->fp, "POLYNOMIAL: %d can't work out how to combine new/old fields\n", number);
                    break;
                }

                if (BC)
                {
                    bcbitf[node] = (bcbitf[node] | bcmask_on);
                    bcbitf[node] = (bcbitf[node] & (~bcmask_off));
                }
            }
        }
    }

    return;
}

/* ===================================  */
double sqrt_multis(jj, ii) int ii, jj;
{
    int i;
    double sqrt_multisa;

    sqrt_multisa = 1.0;
    if (jj > ii)
        for (i = jj; i > ii; i--)
            sqrt_multisa *= 1.0 / sqrt((double)i);

    return sqrt_multisa;
}

/* ===================================  */
double multis(ii) int ii;
{
    int i;
    double multisa;

    multisa = 1.0;
    if (ii)
        for (i = 2; i <= ii; i++)
            multisa *= (double)i;

    return multisa;
}

/* ===================================  */
int int_multis(ii) int ii;
{
    int i, multisa;

    multisa = 1;
    if (ii)
        for (i = 2; i <= ii; i++)
            multisa *= i;

    return multisa;
}

/* get buoyancy from phase diagram files */
void get_phase_buoyancy(E) struct All_variables *E;
{
    static FILE *fp0;
    int i, n1, n2;
    char output_file[255];
    static int been_here = 0;
    int m, n, l, num, number1, number2, number3, number4, j, count_mineral, count_mat, temp_num_mineral;
    int left, right, mid;
    double r1, t1, loc_mid, tempdist;
    const int nno = E->mesh.nno;
    int mark_P, mark_T, mark_P0_PREM = 160, mark_P_PREM;
    double coord_x, coord_z, temp, comp, coord_x_km, coord_z_km, T_C, P_GPa;
    double delta, delta1, delta2, delta_P, delta_T, delta_1, delta_2, delta_3, delta_4;
    double den_basa, Vp_basa, Vs_basa, den_pyro, Vp_pyro, Vs_pyro, den_harz, Vp_harz, Vs_harz;
    double den_mineral[100], Vp_mineral[100], Vs_mineral[100];
    fprintf(stderr, "start buoyancy phasefile\n");
    static int mark_P_store[1001], mark_P_PREM_store[1001];
    if (!been_here)
    {
        for (n = 1; n <= E->mesh.noz; n++)
        {
            mark_P_store[n] = E->control.phasefile_noP;
            mark_P_PREM_store[n] = mark_P0_PREM;
        }
        mark_P = E->control.phasefile_noP;
        mark_P_PREM = mark_P0_PREM;

        for (n = 1; n <= E->mesh.noz; n++)
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
            if (n == E->mesh.noz)
            {
                mark_P_PREM_store[n] = 2;
                mark_P_store[n] = 2;
            }

            mark_P = mark_P_store[n];
            fprintf(stderr, "%d %d %d\n", n, mark_P_PREM_store[n], mark_P_store[n]);
        }
        been_here++;
    }
    
    for (m = 1; m <= E->mesh.nox; m++)
    {
        mark_P = E->control.phasefile_noP;
        mark_T = E->control.phasefile_noT;
        mark_P_PREM = mark_P0_PREM;
        for (n = 1; n <= E->mesh.noz; n++)
        {
            mark_P = mark_P_store[n];
            mark_P_PREM = mark_P_PREM_store[n];
            num = n + E->mesh.noz * (m - 1);
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
            }
*/
            delta = E->control.phase_PREM_depth[mark_P_PREM] - E->control.phase_PREM_depth[mark_P_PREM - 1];
            delta1 = E->control.phase_PREM_depth[mark_P_PREM] - coord_z_km;
            delta2 = coord_z_km - E->control.phase_PREM_depth[mark_P_PREM - 1];
            P_GPa = E->control.phase_PREM_P[mark_P_PREM] * delta2 / delta + E->control.phase_PREM_P[mark_P_PREM - 1] * delta1 / delta;
            /*            for(l=mark_P;l>=2;l--) {
                if(P_GPa<=E->control.phase_P[l] && P_GPa>=E->control.phase_P[l-1]) {
                    mark_P = l;
                    break;
                }
            }
*/
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

            if (delta_2 < 0 || delta_4 < 0)
            {
                fprintf(stderr, "%d %d %d %d %d %d\n", m, n, number1, number2, number3, number4);
                fprintf(stderr, "%d %d %e num1 %e num2 %e num3 %e num4 %e\n", m, n, den_basa, E->control.phase_basa_density[number1], E->control.phase_basa_density[number2], E->control.phase_basa_density[number3], E->control.phase_basa_density[number4]);
                fprintf(stderr, "%d %d %e num1 %e num2 %e num3 %e num4 %e\n", m, n, den_pyro, E->control.phase_pyro_density[number1], E->control.phase_pyro_density[number2], E->control.phase_pyro_density[number3], E->control.phase_pyro_density[number4]);
            }

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
            E->buoyancy[num] -= E->density_phase[num] * E->control.Ra_temp / (E->data.therm_exp * E->data.ref_temperature * E->data.density);

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

        } // end of noz
    }     // end of nox
    return;
}
