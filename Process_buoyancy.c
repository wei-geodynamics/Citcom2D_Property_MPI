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

void process_temp_field(E,ii)
 struct All_variables *E;
    int ii;
{ 
    void heat_flux();
    void output_temp();

    if ( ((ii % E->control.record_every) == 0))    {
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

void heat_flux(E)
    struct All_variables *E;
{
    int e,i,j,node,lnode;
    double *mass,*flux,*SU,*RU;
    double VZ[9],u[9],T[9],dTdz[9],area,uT;

    double xk[3][5];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    void get_global_shape_fn();

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->mesh.nno;
    const int lev = E->mesh.levmax;

    flux = (double *) malloc((1+nno)*sizeof(double));

    for(i=1;i<=nno;i++)   {
      flux[i] = 0.0;
      }
    
    for(e=1;e<=E->mesh.nel;e++) {
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,2,E->mesh.levmax);

        for(j=1;j<=ends;j++)
          VZ[j] = E->V[2][E->ien[e].node[j]];

      for(i=1;i<=ppts;i++)   {
        u[i] = 0.0;
        T[i] = 0.0;
        dTdz[i] = 0.0;
        for(j=1;j<=ends;j++)  {
          u[i] += VZ[j]*E->N.ppt[GNPINDEX(j,i)];
          T[i] += E->T[E->ien[e].node[j]]*E->N.ppt[GNPINDEX(j,i)];
          dTdz[i] += -E->T[E->ien[e].node[j]]*GNx.ppt[GNPXINDEX(1,j,i)];
          }
        }

      uT = 0.0;
      area = 0.0;
      for(i=1;i<=ppts;i++)   {
        uT += u[i]*T[i]*dOmega.ppt[i] + dTdz[i]*dOmega.ppt[i];
        area += dOmega.ppt[i];
        }

      uT /= area;
      for(j=1;j<=ends;j++)  {
        flux[E->ien[e].node[j]] += uT*E->TWW[E->mesh.levmax][e].node[j];
        }
    }             /* end of e */

    for(i=1;i<=nno;i++)   {
      flux[i] = flux[i]*E->Mass[i];
      }

    for(i=1;i<=E->mesh.nsf;i++)   {
      E->slice.shflux[i] = 2*flux[E->surf_node[i]]
                           - flux[E->surf_node[i]-1];

      E->slice.bhflux[i] = 2*flux[E->surf_node[i]-E->mesh.noz+1]
                           - flux[E->surf_node[i]-E->mesh.noz+2];
      }
   free((void *)flux);

  return;  
  }
  
