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

void node_locations(E)
     struct All_variables *E;
{ 
  int lev,i,j,k,ijk[4],ii,d,node;
  double x00,*XX[4],*XG[4],dx[4],dxx[40],dx1,dx2,dc,rc,dr1,dr2,dr3,dr4,dr0,dr5,dr6;
  int n00,nox,noz,noy,fn,step,ncr;

  const int dims = E->mesh.nsd;
  int m;
  FILE *fp_read;
  char input_s[255],output_file[255];  


  input_int("z_grid_layers",&(E->segment.zlayers),"1");
  input_double_vector("zz",E->segment.zlayers,(E->segment.zzlayer));
  input_int_vector("nz",E->segment.zlayers,(E->segment.nzlayer));

  input_int("x_grid_layers",&(E->segment.xlayers),"1");
  input_double_vector("xx",E->segment.xlayers,(E->segment.xxlayer));
  input_int_vector("nx",E->segment.xlayers,(E->segment.nxlayer));
     for(d=1;d<=E->mesh.nsd;d++) {
       XX[d] = (double *)malloc((2+E->mesh.nnx[d])*sizeof(double));
       XG[d] = (double *)malloc ((2+1)*sizeof(double));
       }

  if(E->control.meshReadIn) {
  /*    input_string("fileMeshX",E->control.fileMeshX,"MeshX.txt"); */
      sprintf(output_file,"%s",E->control.fileMeshX);
      fp_read = fopen(output_file,"r");
      for (m=1;m<=E->segment.nxlayer[E->segment.xlayers-1];m++) {
          fgets(input_s,200,fp_read);
          sscanf(input_s,"%lf",&XX[1][m]);
          fprintf(stderr,"%d %lf\n",m,XX[1][m]);
      }
      fclose(fp_read);
/*      input_string("fileMeshZ",E->control.fileMeshZ,"MeshZ.txt"); */
      sprintf(output_file,"%s",E->control.fileMeshZ);
      fp_read = fopen(output_file,"r");
      for (m=1;m<=E->segment.nzlayer[E->segment.zlayers-1];m++) {
          fgets(input_s,200,fp_read);
          sscanf(input_s,"%lf",&XX[2][m]);
          fprintf(stderr,"%d %lf\n",m,XX[2][m]);
      }
      fclose(fp_read);



  }  else {


     dx[1] = E->mesh.layer[1]/(E->mesh.nnx[1]-1);
     XX[1][1] = 0.0;
     for(i=2;i<=E->mesh.nnx[1];i++)
  	      XX[1][i] = XX[1][i-1]+dx[1];

     dx[2] = E->mesh.layer[2]/(E->mesh.nnx[2]-1);
     XX[2][1] = 0.0;
     for(i=2;i<=E->mesh.nnx[2];i++)
  	      XX[2][i] = XX[2][i-1]+dx[2];

  for (j=1;j<E->segment.xlayers;j++)
    dxx[j] = (E->segment.xxlayer[j]-E->segment.xxlayer[j-1])
            /(E->segment.nxlayer[j]-E->segment.nxlayer[j-1]);
  j=1;
  for(i=2;i<E->mesh.nnx[1];i++)   {
    if (i<=E->segment.nxlayer[j])
       XX[1][i] = XX[1][i-1] + dxx[j];
    if (i==E->segment.nxlayer[j])
       j++;
    }


  for (j=1;j<E->segment.zlayers;j++)
    dxx[j] = (E->segment.zzlayer[j]-E->segment.zzlayer[j-1])
            /(E->segment.nzlayer[j]-E->segment.nzlayer[j-1]);
  j=1;
  for(i=2;i<E->mesh.nnx[2];i++)   {
    if (i<=E->segment.nzlayer[j])
       XX[2][i] = XX[2][i-1] + dxx[j];
    if (i==E->segment.nzlayer[j])
       j++;
    }

  }

   for(d=1;d<=E->mesh.nsd;d++)
     for(i=1;i<=E->mesh.nnx[d];i++)
       E->XP[d][i] = XX[d][i]; 
    
   for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {

     nox=E->mesh.NOX[lev]; 
     noy=E->mesh.NOY[lev];
     noz=E->mesh.NOZ[lev];

    if (E->control.NMULTIGRID||E->control.EMULTIGRID)
        step = (int) pow(2.0,(double)(E->mesh.levmax-lev));
    else
        step = 1;

     for(ijk[1]=1;ijk[1]<=nox;ijk[1]++)
       for(ijk[2]=1;ijk[2]<=noz;ijk[2]++)
         for(ijk[3]=1;ijk[3]<=noy;ijk[3]++)   {
           node=ijk[2]+(ijk[1]-1)*noz+(ijk[3]-1)*noz*nox;
           for(d=1;d<=E->mesh.nsd;d++)
             E->XX[lev][d][node] = XX[d][(ijk[d]-1)*step+1]; 
           }
      }

     for(d=1;d<=E->mesh.nsd;d++) {
       free((void *)XX[d]);
       free((void *)XG[d]);
       }

  if (E->control.verbose) 
    for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {
      fprintf(E->fp,"output_coordinates %d\n",lev);
      if (dims==2)  {
         for (i=1;i<=E->mesh.NNO[lev];i++)
             fprintf(E->fp,"%d %lf %lf\n",i,E->XX[lev][1][i],E->XX[lev][2][i]);
         }
      else if (dims==3)  {
         for (i=1;i<=E->mesh.NNO[lev];i++)
             fprintf(E->fp,"%d %lf %lf %lf\n",i,E->XX[lev][1][i],E->XX[lev][2][i],E->XX[lev][3][i]);
         }
      }

return;  }




void flogical_mesh_to_real(E,data,level)
     struct All_variables *E;
     double *data;
     int level;

{ int i,j,n1,n2;

  return;
}

void dp_to_nodes(E,P,PN,lev)
     struct All_variables *E;
     double *P;
     double *PN;
     int lev;

{ int e,element,node,j;

  for(node=1;node<=E->mesh.NNO[lev];node++)
    PN[node] =  0.0;
	  
  for(element=1;element<=E->mesh.NEL[lev];element++) {

      for(j=1;j<=enodes[E->mesh.nsd];j++)  {
     	  node = E->IEN[lev][element].node[j];
    	  PN[node] += P[element] * E->TW[lev][node] ; 
    	  }

      } 

     return; }

void p_to_nodes(E,P,PN,lev)
     struct All_variables *E;
     double *P,*PN;
     int lev;

{ int e,element,node,j;

  for(node=1;node<=E->mesh.NNO[lev];node++)
    PN[node] =  0.0;
	  
  for(element=1;element<=E->mesh.NEL[lev];element++) {

      for(j=1;j<=enodes[E->mesh.nsd];j++)  {
     	  node = E->IEN[lev][element].node[j];
    	  PN[node] += P[element] * E->TW[lev][node] ; 
    	  }

      } 

     return; }


void p_to_centres(E,PN,P,lev)
     struct All_variables *E;
     double *PN,*P;
     int lev;

{  int p,element,node,j;
   double weight;

   for(p=1;p<=E->mesh.NEL[lev];p++)
     P[p] = 0.0;

   weight=1.0/((double)enodes[E->mesh.nsd]) ;
   
   for(p=1;p<=E->mesh.NEL[lev];p++)
     for(j=1;j<=enodes[E->mesh.nsd];j++)
       P[p] +=  PN[E->IEN[lev][p].node[j]] * weight;

   return;  
   }


void v_to_intpts(E,VN,VE,lev)
  struct All_variables *E;
  double *VN,*VE;
  int lev;
  {

   int e,i,j,k;
   const int nsd=E->mesh.nsd;
   const int vpts=vpoints[nsd];
   const int ends=enodes[nsd];

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)                 {
        VE[(e-1)*vpts + i] = 0.0;
        for(j=1;j<=ends;j++)
          VE[(e-1)*vpts + i] += VN[E->IEN[lev][e].node[j]] *  E->N.vpt[GNVINDEX(j,i)];
        }

   return;
  }

void v_to_nodes(E,VE,VN,lev)
   struct All_variables *E;
   double *VE,*VN;
   int lev;
   {
    int e,i,j,k,n;
    const int nsd=E->mesh.nsd;
    const int vpts=vpoints[nsd];
    const int ends=enodes[nsd];

    for(i=1;i<=E->mesh.NNO[lev];i++)
	    VN[i] = 0.0;

    for(e=1;e<=E->mesh.NEL[lev];e++)
      for(j=1;j<=ends;j++) {
        n = E->IEN[lev][e].node[j];
        for(i=1;i<=vpts;i++)
          VN[n] += E->N.vpt[GNVINDEX(j,i)] * E->TW[lev][n] * VE[(e-1)*vpts + i];
        }

    return;
    }

void visc_to_intpts(E,VN,VE,lev)
   struct All_variables *E;
   double *VN,*VE;
   int lev;
   {

   int e,i,j,k;
   const int nsd=E->mesh.nsd;
   const int vpts=vpoints[nsd];
   const int ends=enodes[nsd];

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++) {
        VE[(e-1)*vpts + i] = 0.0;
	for(j=1;j<=ends;j++)
          VE[(e-1)*vpts + i] += log(VN[E->IEN[lev][e].node[j]]) *  E->N.vpt[GNVINDEX(j,i)];
        VE[(e-1)*vpts + i] = exp(VE[(e-1)*vpts + i]);
        }

  }


void visc_to_nodes(E,VE,VN,lev)
  struct All_variables *E;
  double *VE,*VN;
  int lev;
  {
  int e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

  for(i=1;i<=E->mesh.NNO[lev];i++)
    VN[i] = 0.0;

  for(e=1;e<=E->mesh.NEL[lev];e++)
    for(j=1;j<=ends;j++) {
      n = E->IEN[lev][e].node[j];
      temp_visc=0.0;
      for(i=1;i<=vpts;i++)
	temp_visc += E->TW[lev][n] * log(E->N.vpt[GNVINDEX(j,i)] * VE[(e-1)*vpts + i]);
      VN[n] += exp(temp_visc);
      }
   return;
}

void visc_from_gint_to_nodes(E,VE,VN,lev)
  struct All_variables *E;
  double *VE,*VN;
  int lev;
  {
  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(i=1;i<=E->mesh.NNO[lev];i++)
     VN[i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(j=1;j<=ends;j++)                {
       n = E->IEN[lev][e].node[j];
       temp_visc=0.0;
       for(i=1;i<=vpts;i++)
         temp_visc += E->TW[lev][n] * E->N.vpt[GNVINDEX(j,i)] * VE[(e-1)*vpts + i];
       VN[n] += temp_visc;
       }

   return;
}


 void visc_from_nodes_to_gint(E,VN,VE,lev)
  struct All_variables *E;
  double *VE,*VN;
  int lev;
  {

  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;
   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)
       VE[(e-1)*vpts+i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)      {
       temp_visc=0.0;
       for(j=1;j<=ends;j++)
         temp_visc += E->N.vpt[GNVINDEX(j,i)]*VN[E->IEN[lev][e].node[j]];

       VE[(e-1)*vpts+i] = temp_visc;
       }

   return;
   }

void visc_from_gint_to_ele(E,VE,VN,lev)
  struct All_variables *E;
  double *VE,*VN;
  int lev;
  {
  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(i=1;i<=E->mesh.NEL[lev];i++)
     VN[i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)   {
     temp_visc=0.0;
     for(i=1;i<=vpts;i++)
        temp_visc += VE[(e-1)*vpts + i];
     temp_visc = temp_visc/vpts;

     VN[e] = temp_visc;
    }

   return;
}


 void visc_from_ele_to_gint(E,VN,VE,lev)
  struct All_variables *E;
  double *VE,*VN;
  int lev;
  {

  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)
       VE[(e-1)*vpts+i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)      {

       VE[(e-1)*vpts+i] = VN[e];
       }

   return;
 }

