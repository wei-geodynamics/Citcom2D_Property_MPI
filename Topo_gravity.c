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
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

#define c_re(a) a.real
#define c_im(a) a.imag
typedef struct compl
{
    double real;
    double imag;
}
COMPLEX;

/* ===================================================================
   Consistent boundary flux method for stress ... Zhong,Gurnis,Hulbert 

   Solve for the stress as the code defined it internally, rather than
   what was intended to be solved. This is more appropriate.

   Note also that the routine is dependent on the method 
   used to solve for the velocity in the first place.
   ===================================================================  */

/* call this only for top and bottom processors */
void get_CBF_topo(E, H, HB) /* call this only for top and bottom processors*/
    struct All_variables *E;
double *H, *HB;

{
    void get_elt_k();
    void get_elt_g();
    void get_elt_f();
    void matrix_transform_g();
    void get_global_1d_shape_fn();
    void exchange_node_f20();

    int el, elb, els, node, nodeb, nodes, i, j, k, l, m, n, count;
    int nodel, nodem, nodesl, nodesm, lnsf, nel2;

    struct Shape_function1 GM, GMb;
    struct Shape_function1_dA dGammax, dGammabx;

    double *temp, *eltTU, *eltTL, *SU, *SL, *RU, *RL;

    double eltk[24 * 24], eltf[24];
    double eltkb[24 * 24], eltfb[24];
    double res[24], resb[24], eu[24], eub[24];
    higher_precision eltg[24][1], eltgb[24][1];

    const int dims = E->mesh.nsd;
    const int Tsize = 5; /* maximum values, applicable to 3d, harmless for 2d */
    const int Ssize = 4;
    const int ends = enodes[dims];
    const int noz = E->lmesh.noz;
    const int noy = E->lmesh.noy;
    const int nno = E->lmesh.nno;
    const int onedv = onedvpoints[dims];
    const int snode1 = 1, snode2 = 4, snode3 = 5, snode4 = 8;
    const int elz = E->lmesh.elz;
    const int ely = E->lmesh.ely;
    const int lev = E->mesh.levmax;

    lnsf = E->lmesh.nsf;

    eltTU = (double *)malloc((1 + Tsize) * sizeof(double));
    eltTL = (double *)malloc((1 + Tsize) * sizeof(double));
    SU = (double *)malloc((1 + lnsf) * sizeof(double));
    SL = (double *)malloc((1 + lnsf) * sizeof(double));
    RU = (double *)malloc((1 + lnsf) * sizeof(double));
    RL = (double *)malloc((1 + lnsf) * sizeof(double));
    temp = (double *)malloc((1 + nno) * sizeof(double));

    for (i = 0; i <= nno; i++)
        temp[i] = 0;

    for (i = 0; i <= lnsf; i++)
        RU[i] = RL[i] = SU[i] = SL[i] = 0.0;

    /* calculate the element residuals */

    for (els = 1; els <= E->lmesh.snel; els++)
    {
        el = E->surf_element[els];
        elb = el + elz - 1;

        for (m = 0; m < ends; m++)
        { /* for bottom elements */
            nodeb = E->ien[elb].node[m + 1];
            eub[m * dims] = E->V[1][nodeb];
            eub[m * dims + 1] = E->V[2][nodeb];
            if (3 == dims)
                eub[m * dims + 2] = E->V[3][nodeb];
        }

        for (m = 0; m < ends; m++)
        {
            node = E->ien[el].node[m + 1];
            eu[m * dims] = E->V[1][node];
            eu[m * dims + 1] = E->V[2][node];
            if (3 == dims)
                eu[m * dims + 2] = E->V[3][node];
        }

        get_elt_k(E, el, eltk, lev, 0);
        get_elt_k(E, elb, eltkb, lev, 0);
        get_elt_f(E, el, eltf, 0, 0);
        get_elt_f(E, elb, eltfb, 0, 0);
        get_elt_g(E, el, eltg, lev);
        get_elt_g(E, elb, eltgb, lev);

        for (m = 0; m < dims * ends; m++)
        {
            res[m] = eltf[m] - E->elt_del[lev][el].g[m][0] * E->P[el];
            resb[m] = eltfb[m] - E->elt_del[lev][elb].g[m][0] * E->P[elb];
        }

        for (m = 0; m < dims * ends; m++)
            for (l = 0; l < dims * ends; l++)
            {
                res[m] -= eltk[ends * dims * m + l] * eu[l];
                resb[m] -= eltkb[ends * dims * m + l] * eub[l];
            }

        /* Put relevant (vertical & surface) parts of element residual into surface residual */

        for (m = 1; m <= ends; m++)
        { /* for bottom elements */
            switch (m)
            {
            case 2:
                RL[E->sien[els].node[1]] += resb[(m - 1) * dims + 1];
                break;
            case 3:
                RL[E->sien[els].node[2]] += resb[(m - 1) * dims + 1];
                break;
            case 7:
                RL[E->sien[els].node[3]] += resb[(m - 1) * dims + 1];
                break;
            case 6:
                RL[E->sien[els].node[4]] += resb[(m - 1) * dims + 1];
                break;
            }
        }

        for (m = 1; m <= ends; m++)
        {
            switch (m)
            {
            case 1:
                RU[E->sien[els].node[1]] += res[(m - 1) * dims + 1];
                break;
            case 4:
                RU[E->sien[els].node[2]] += res[(m - 1) * dims + 1];
                break;
            case 8:
                RU[E->sien[els].node[3]] += res[(m - 1) * dims + 1];
                break;
            case 5:
                RU[E->sien[els].node[4]] += res[(m - 1) * dims + 1];
                break;
            }
        }
    }

    /* calculate the LHS */

    for (els = 1; els <= E->lmesh.snel; els++)
    {
        el = E->surf_element[els];
        elb = el + elz - 1;

        get_global_1d_shape_fn(E, el, &GM, &dGammax, 1);
        get_global_1d_shape_fn(E, elb, &GMb, &dGammabx, 1);

        for (m = 1; m <= onedv; m++)
        {
            eltTU[m - 1] = 0.0;
            eltTL[m - 1] = 0.0;
            for (n = 1; n <= onedv; n++)
            {
                eltTU[m - 1] +=
                    dGammax.vpt[GMVGAMMA(1, n)] * l_1d[n].weight[dims - 1] * E->L.vpt[GMVINDEX(m, n)] * E->L.vpt[GMVINDEX(m, n)];
                eltTL[m - 1] +=
                    dGammabx.vpt[GMVGAMMA(1 + dims, n)] * l_1d[n].weight[dims - 1] * E->L.vpt[GMVINDEX(m, n)] * E->L.vpt[GMVINDEX(m, n)];
            }
        }

        for (m = 1; m <= onedv; m++) /* for bottom */
            SL[E->sien[els].node[m]] += eltTL[m - 1];

        for (m = 1; m <= onedv; m++)
            SU[E->sien[els].node[m]] += eltTU[m - 1];
    }

    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        temp[node] = RU[i];
    }
    exchange_node_f20(E, temp, E->mesh.levmax);
    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        RU[i] = temp[node];
    }
    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        temp[node] = SU[i];
    }
    exchange_node_f20(E, temp, E->mesh.levmax);
    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        SU[i] = temp[node];
    }

    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        temp[node] = RL[i];
    }
    exchange_node_f20(E, temp, E->mesh.levmax);
    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        RL[i] = temp[node];
    }
    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        temp[node] = SL[i];
    }
    exchange_node_f20(E, temp, E->mesh.levmax);
    for (i = 1; i <= E->lmesh.nsf; i++)
    {
        node = E->surf_node[i];
        SL[i] = temp[node];
    }

    if (E->parallel.me_loc[3] == 0)
    {
        for (i = 1; i <= E->lmesh.nsf; i++)
            H[i] = -RU[i] / SU[i];
    }

    if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
    {
        for (i = 1; i <= E->lmesh.nsf; i++)
            HB[i] = -RL[i] / SL[i];
    }

    free((void *)eltTU);
    free((void *)eltTL);
    free((void *)SU);
    free((void *)SL);
    free((void *)RU);
    free((void *)RL);
    free((void *)temp);
    return;
}

void get_STD_topo(E, tpg, tpgb, ii) struct All_variables *E;
double *tpg;
double *tpgb;
int ii;
{

    void get_surf_stress();
    void get_global_shape_fn();
    void exchange_node_f20();
    int i, j, k, e, nel2, snode, node;

    double *SZZ, *SXX, *SYY, *SXY, *SXZ, *SZY, VZ[9], VY[9], VX[9], Szz, Sxx, Syy, Sxy, Sxz, Szy;
    double Vzz[9], Vxx[9], Vyy[9], Vxy[9], Vxz[9], Vzy[9];
    double pre[9], el_volume, tww[9], Visc, a, b;

    const int dims = E->mesh.nsd, dofs = E->mesh.dof;
    const int vpts = vpoints[dims];
    const int ppts = ppoints[dims];
    const int ends = enodes[dims];
    const int nno = E->lmesh.nno;
    const int nel = E->lmesh.nel;
    const int lev = E->mesh.levmax;

    SXX = (double *)malloc((nno + 1) * sizeof(double));
    SYY = (double *)malloc((nno + 1) * sizeof(double));
    SXY = (double *)malloc((nno + 1) * sizeof(double));
    SXZ = (double *)malloc((nno + 1) * sizeof(double));
    SZY = (double *)malloc((nno + 1) * sizeof(double));
    SZZ = (double *)malloc((nno + 1) * sizeof(double));

    for (i = 1; i <= nno; i++)
    {
        SZZ[i] = 0.0;
        SXX[i] = 0.0;
        SYY[i] = 0.0;
        SXY[i] = 0.0;
        SXZ[i] = 0.0;
        SZY[i] = 0.0;
    }

    for (e = 1; e <= nel; e++)
    {
        Szz = 0.0;
        Sxx = 0.0;
        Syy = 0.0;
        Sxy = 0.0;
        Sxz = 0.0;
        Szy = 0.0;

        for (j = 1; j <= vpts; j++)
        {
            pre[j] = E->EVI[lev][(e - 1) * vpts + j] * E->gDA[e].vpt[j];
            Vzz[j] = 0.0;
            Vxx[j] = 0.0;
            Vyy[j] = 0.0;
            Vxy[j] = 0.0;
            Vxz[j] = 0.0;
            Vzy[j] = 0.0;
        }

        for (j = 1; j <= ends; j++)
        {
            VX[j] = E->V[1][E->ien[e].node[j]];
            VZ[j] = E->V[2][E->ien[e].node[j]];
            if (dims == 3)
                VY[j] = E->V[3][E->ien[e].node[j]];
        }

        for (i = 1; i <= vpts; i++)
        {
            for (j = 1; j <= ends; j++)
            {
                Vzz[i] += VZ[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)];
                Vxx[i] += VX[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)];
                Vxz[i] += (VX[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)] + VZ[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)]);
                if (dims == 3)
                {
                    Vyy[i] += VY[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)];
                    Vxy[i] += (VX[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)] + VY[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)]);
                    Vzy[i] += (VY[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)] + VZ[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)]);
                }
            }
            Sxx += 2.0 * pre[i] * Vxx[i];
            Syy += 2.0 * pre[i] * Vyy[i];
            Szz += 2.0 * pre[i] * Vzz[i];
            Sxy += pre[i] * Vxy[i];
            Sxz += pre[i] * Vxz[i];
            Szy += pre[i] * Vzy[i];
        }

        Sxx /= E->eco[e].area;
        Syy /= E->eco[e].area;
        Szz /= E->eco[e].area;
        Sxz /= E->eco[e].area;
        Sxy /= E->eco[e].area;
        Szy /= E->eco[e].area;

        Szz -= E->P[e]; /* add the pressure term */
        Sxx -= E->P[e]; /* add the pressure term */
        if (dims == 3)
            Syy -= E->P[e]; /* add the pressure term */

        for (j = 1; j <= ends; j++)
        {
            node = E->ien[e].node[j];
            SZZ[node] += E->TWW[E->mesh.levmax][e].node[j] * Szz;
            SXX[node] += E->TWW[E->mesh.levmax][e].node[j] * Sxx;
            SXZ[node] += E->TWW[E->mesh.levmax][e].node[j] * Sxz;
            if (dims == 3)
            {
                SYY[node] += E->TWW[E->mesh.levmax][e].node[j] * Syy;
                SXY[node] += E->TWW[E->mesh.levmax][e].node[j] * Sxy;
                SZY[node] += E->TWW[E->mesh.levmax][e].node[j] * Szy;
            }
        }
    }

    exchange_node_f20(E, SXX, lev);
    exchange_node_f20(E, SZZ, lev);
    exchange_node_f20(E, SXZ, lev);
    if (dims == 3)
    {
        exchange_node_f20(E, SXY, lev);
        exchange_node_f20(E, SYY, lev);
        exchange_node_f20(E, SZY, lev);
    }

    for (i = 1; i <= nno; i++)
    {
        SZZ[i] = SZZ[i] * E->Mass[i];
        SXX[i] = SXX[i] * E->Mass[i];
        SXZ[i] = SXZ[i] * E->Mass[i];
    }
    if (dims == 3)
        for (i = 1; i <= nno; i++)
        {
            SYY[i] = SYY[i] * E->Mass[i];
            SXY[i] = SXY[i] * E->Mass[i];
            SZY[i] = SZY[i] * E->Mass[i];
        }

    for (snode = 1; snode <= E->lmesh.nsf; snode++)
    {
        node = E->surf_node[snode];
        tpg[snode] = -2 * SZZ[node] + SZZ[node + 1];
        tpgb[snode] = 2 * SZZ[node + E->lmesh.noz - 1] - SZZ[node + E->lmesh.noz - 2];
    }

    free((void *)SXX);
    free((void *)SYY);
    free((void *)SXY);
    free((void *)SXZ);
    free((void *)SZY);
    free((void *)SZZ);

    return;
}
