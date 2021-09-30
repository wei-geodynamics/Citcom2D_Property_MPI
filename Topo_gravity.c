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
  typedef struct compl {  double real;
			  double imag;    } COMPLEX;

extern double Zj0[1000],Zj1[1000];


/* ===================================================================
   Consistent boundary flux method for stress ... Zhong,Gurnis,Hulbert 

   Solve for the stress as the code defined it internally, rather than
   what was intended to be solved. This is more appropriate.

   Note also that the routine is dependent on the method 
   used to solve for the velocity in the first place.
   ===================================================================  */



void get_CBF_topo(E,H,HB)       /* call this only for top and bottom processors*/
    struct All_variables *E;
    double *H,*HB;
   
{
    void get_elt_k();
    void get_elt_g();
    void get_elt_f();
    void matrix_transform_g();
    void get_global_1d_shape_fn();

    int el,elb,els,node,nodeb,nodes,i,j,k,l,m,n,count;
    int nodel,nodem,nodesl,nodesm,lnsf,nel2;
    
    struct Shape_function1 GM,GMb;
    struct Shape_function1_dA dGammax,dGammabx;
 
    double *eltTU,*eltTL,*SU,*SL,*RU,*RL;

    double eltk[24*24],eltf[24];
    double eltkb[24*24],eltfb[24];
    double res[24],resb[24],eu[24],eub[24];
    higher_precision eltg[24][1],eltgb[24][1];

    const int dims=E->mesh.nsd;
    const int Tsize=5;   /* maximum values, applicable to 3d, harmless for 2d */ 
    const int Ssize=4;
    const int ends=enodes[dims];
    const int noz=E->mesh.noz;
    const int noy=E->mesh.noy;
    const int nno=E->mesh.nno;
    const int onedv=onedvpoints[dims];
    const int snode1=1,snode2=4,snode3=5,snode4=8;
    const int elz = E->mesh.elz;
    const int ely = E->mesh.ely;
    const int lev=E->mesh.levmax;
 
    lnsf=E->mesh.nsf;
 
    eltTU = (double *)malloc((1+Tsize)*sizeof(double)); 
    eltTL = (double *)malloc((1+Tsize)*sizeof(double));
    SU = (double *)malloc((1+lnsf)*sizeof(double));
    SL = (double *)malloc((1+lnsf)*sizeof(double));
    RU = (double *)malloc((1+lnsf)*sizeof(double));
    RL = (double *)malloc((1+lnsf)*sizeof(double));

    for(i=0;i<=lnsf;i++)
      RU[i] = RL[i] = SU[i] = SL[i] = 0.0;

    /* calculate the element residuals */

    for(els=1;els<=E->mesh.snel;els++) {
	    el = E->surf_element[els];
	    elb = el + elz-1;

	    for(m=0;m<ends;m++) {  /* for bottom elements */
          nodeb= E->ien[elb].node[m+1];
          eub[m*dims  ] = E->V[1][nodeb];
          eub[m*dims+1] = E->V[2][nodeb];
          if(3==dims) 
            eub[m*dims+2] = E->V[3][nodeb]; 
          }

	      for(m=0;m<ends;m++) {  
          node = E->ien[el].node[m+1];
          eu [m*dims  ] = E->V[1][node];
          eu [m*dims+1] = E->V[2][node];
          if(3==dims)
            eu [m*dims+2] = E->V[3][node];
          }

	    get_elt_k(E,el,eltk,lev,0);
	    get_elt_k(E,elb,eltkb,lev,0);
	    get_elt_f(E,el,eltf,0,0);
	    get_elt_f(E,elb,eltfb,0,0);
        get_elt_g(E,el,eltg,lev);
	    get_elt_g(E,elb,eltgb,lev);
	   
	    for(m=0;m<dims*ends;m++) {
          res[m]  = eltf[m]  - E->elt_del[lev][el].g[m][0]  * E->P[el];
          resb[m] = eltfb[m] - E->elt_del[lev][elb].g[m][0]* E->P[elb];
	      }
	   
	    for(m=0;m<dims*ends;m++)
		for(l=0;l<dims*ends;l++) {
		    res[m]  -= eltk[ends*dims*m+l]  * eu[l];
		    resb[m] -= eltkb[ends*dims*m+l] * eub[l];
		   }
	   
	    /* Put relevant (vertical & surface) parts of element residual into surface residual */
		
	    for(m=1;m<=ends;m++) {     /* for bottom elements */
		switch (m) {
		case 2:
		    RL[E->sien[els].node[1]] += resb[(m-1)*dims+1];  
		    break;
		case 3:
		    RL[E->sien[els].node[2]] += resb[(m-1)*dims+1];  
		    break;
		case 7:
		    RL[E->sien[els].node[3]] += resb[(m-1)*dims+1];  
		    break;
		case 6:
		    RL[E->sien[els].node[4]] += resb[(m-1)*dims+1];  
		    break;
		    }
		}


	    for(m=1;m<=ends;m++) {
		switch (m) {
		case 1:
		    RU[E->sien[els].node[1]] += res[(m-1)*dims+1];  
		    break;
		case 4:
		    RU[E->sien[els].node[2]] += res[(m-1)*dims+1];  
		    break;
		case 8:
		    RU[E->sien[els].node[3]] += res[(m-1)*dims+1];  
		    break;
		case 5:
		    RU[E->sien[els].node[4]] += res[(m-1)*dims+1];  
		    break;
		    }
		}
	}
    
    /* calculate the LHS */
 
    for(els=1;els<=E->mesh.snel;els++) {
	    el = E->surf_element[els];
	    elb = el + elz-1;

	    get_global_1d_shape_fn(E,el,&GM,&dGammax,1);
	    get_global_1d_shape_fn(E,elb,&GMb,&dGammabx,1);
   
	    for(m=1;m<=onedv;m++)        {
	      eltTU[m-1] = 0.0;
	      eltTL[m-1] = 0.0; 
	      for(n=1;n<=onedv;n++)          {
		     eltTU[m-1] += 
		         dGammax.vpt[GMVGAMMA(1,n)] * l_1d[n].weight[dims-1]
			 * E->L.vpt[GMVINDEX(m,n)] * E->L.vpt[GMVINDEX(m,n)];
		     eltTL[m-1] += 
			 dGammabx.vpt[GMVGAMMA(1+dims,n)]*l_1d[n].weight[dims-1]
			 * E->L.vpt[GMVINDEX(m,n)] * E->L.vpt[GMVINDEX(m,n)];
		     }
		  }

        for (m=1;m<=onedv;m++)     /* for bottom */
	      SL[E->sien[els].node[m]] += eltTL[m-1];

          for (m=1;m<=onedv;m++) 
	           SU[E->sien[els].node[m]] += eltTU[m-1];
	    }


      for(i=1;i<=E->mesh.nsf;i++)
        H[i] = -RU[i]/SU[i];
        
      for(i=1;i<=E->mesh.nsf;i++)
        HB[i] = -RL[i]/SL[i];

    free((void *)eltTU);
    free((void *)eltTL);
    free((void *)SU);
    free((void *)SL);
    free((void *)RU);
    free((void *)RL);
    return;

}

void get_STD_topo(E,tpg,tpgb,ii)
    struct All_variables *E;
    double *tpg;
    double *tpgb;
    int ii;
{
    int i,j,k,snode,node;
    double *Szz;
    void get_Szz();

    Szz=(double *) malloc((1+E->mesh.nno)*sizeof(double));
    
    get_Szz(E,E->P,Szz,ii);

    for(snode=1;snode<=E->mesh.nsf;snode++)   {
           node = E->surf_node[snode];
	   tpg[snode]  = -2*Szz[node] + Szz[node-1];
	   tpgb[snode] = 2*Szz[node-E->mesh.noz+1]-Szz[node-E->mesh.noz+2];
           //fprintf(stderr,"Topo %d  %.4e %.4e\n",snode, Szz[node],Szz[node-1]);
         
	   }

    free((void *)Szz);
    return;
}

void get_Szz(E,P,SZZ,file_number)
     struct All_variables *E;
     double *P;
     int file_number;
     double *SZZ;

{
    void get_surf_stress();
    void get_global_shape_fn();
    int i,j,k,e,node, nel2;
    
    double *SXX,*SYY,*SXY,*SXZ,*SZY,VZ[9],VY[9],VX[9],Szz,Sxx,Syy,Sxy,Sxz,Szy;
    double EXX[9],EZZ[9],EXZ[9],EYY[9];
    double el_volume,visc[9],Visc,a,b,xk[3][5];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;

/*
    SXX = (double *)malloc((E->mesh.nno+1)*sizeof(double));
    SYY = (double *)malloc((E->mesh.nno+1)*sizeof(double));
    SXY = (double *)malloc((E->mesh.nno+1)*sizeof(double));
    SXZ = (double *)malloc((E->mesh.nno+1)*sizeof(double));
    SZY = (double *)malloc((E->mesh.nno+1)*sizeof(double));
*/

    for(i=1;i<=E->mesh.nno;i++) {
      SZZ[i] = 0.0;
      }

    for(e=1;e<=E->mesh.nel;e++)  {
      Szz = 0.0;
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,0,E->mesh.levmax);
	
      for(i=1;i<=vpts;i++)   {
	  visc[i] =  E->EVI[E->mesh.levmax][(e-1)*vpts+i] * dOmega.vpt[i];
          EZZ[i] = 0.0;
          }
     
      for(j=1;j<=ends;j++) {
          VZ[j] = E->V[2][E->ien[e].node[j]];
          if(E->X[2][E->ien[e].node[j]]-1.0<=0.00001&&E->X[2][E->ien[e].node[j]]-1.0>=-0.00001)
            VZ[j] = 0.0;
	  }

      for(i=1;i<=vpts;i++)  {
        for(j=1;j<=ends;j++)  {
          EZZ[i] += VZ[j]* GNx.vpt[GNVXINDEX(1,j,i)];
	  }
        Szz += 2.0 * visc[i] * EZZ[i];
	}

      Szz /= E->eco[e].area;
      for(j=1;j<=ends;j++) {
        node = E->ien[e].node[j];
        //if(node%E->mesh.noz==0)
          //fprintf(stderr,"Szz1 %d %.4e %.4e %lf\n",node, E->X[1][node],E->X[2][node], SZZ);
      }

      Szz -= P[e];  /* add the pressure term */
      for(j=1;j<=ends;j++) {
        node = E->ien[e].node[j];
        //if(node%E->mesh.noz==0)
          //fprintf(stderr,"Szz2 %d %.4e %.4e %lf\n",node, E->X[1][node],E->X[2][node], SZZ);
      }

     
      for(j=1;j<=ends;j++) {
	    node = E->ien[e].node[j];
	    SZZ[node] += E->TWW[E->mesh.levmax][e].node[j] * Szz;  
        }

      }

   for(node=1;node<=E->mesh.nno;node++)  {
     SZZ[node] = SZZ[node]*E->Mass[node];
     //if(node%E->mesh.noz==0) 
       //fprintf(stderr,"Szz3 %d %.4e %.4e  %lf\n",node, E->X[1][node], E->X[2][node],SZZ[node]);

   }


/*
   get_surf_stress(E,SXX,SYY,SZZ,SXY,SXZ,SZY);

    free((void *)SXX);
    free((void *)SYY);
    free((void *)SXY);
    free((void *)SXZ);
    free((void *)SZY);
*/

    return; 
}   

void sphere_expansion(E,TG,sphc,sphs)
     struct All_variables *E;
     double *TG,*sphc,*sphs;
{
    int maxindice;
    int i,j,ll,mm,node;
    double sum,dx;
    int nx=E->mesh.nox;
    int nz=E->mesh.noz;
    double dx_all;
    dx_all = E->X[1][E->mesh.nno] -  E->X[1][1];
    maxindice = (E->control.llmax);
    for (i=0;i<maxindice;i++)    {
        sphc[i] = 0.0;
        sphs[i] = 0.0;
    }
    for (i=1;i<=nx;i++) {
        node = 1+(i-1)*nz;
        if(i==0) {
           dx = (E->X[1][node+nz] - E->X[1][1])/2.0;
           
        }
        else if(i==nx) {
           dx = (E->X[1][node] -  E->X[1][node-nz])/2.0;
        }
        else {
           dx = (E->X[1][node+nz] - E->X[1][node-nz])/2.0;
        }
        for (ll=0;ll<E->control.llmax;ll++) {
           sphc[ll] += TG[i]*2.0*dx/dx_all*sin(E->X[1][node]*2.0*M_PI/(dx_all/(double)(ll+1)));
           sphs[ll] += TG[i]*2.0*dx/dx_all*cos(E->X[1][node]*2.0*M_PI/(dx_all/(double)(ll+1)));
                    
           //if (ll==0) 
           //fprintf(stderr,"%.4e %.4e %.4e %d  %.4e\n",dx,sin(E->X[1][node]*2.0*M_PI/(dx_all/(double)(ll+1))),2.0*dx/dx_all*sin(E->X[1][node]*2.0*M_PI/(dx_all/(double)(ll+1)))
           //            *sin(E->X[1][node]*2.0*M_PI/(dx_all/(double)(ll+1))), i, sphc[ll]);

        }         
    }
return;
}
void expand_topo_sph_harm(struct All_variables *E,
                                 double *tpgt[2],
                                 double *tpgb[2])
{
    double scaling, stress_scaling, topo_scaling1,topo_scaling2;
    double den_contrast1, den_contrast2, grav1, grav2;
    int i, j;
    int maxindice;
    void sphere_expansion();
    maxindice = (E->control.llmax);
 
    stress_scaling = E->data.ref_viscosity*E->data.therm_diff/
        (E->data.layer_km*E->data.layer_km*1e6);
    den_contrast1 = E->data.density - E->data.density_above;
    den_contrast2 = E->data.density_below - E->data.density;
    grav1 = E->data.grav_acc;
    grav2 = E->data.grav_acc;
    topo_scaling1 = stress_scaling / (den_contrast1 * grav1);
    topo_scaling2 = stress_scaling / (den_contrast2 * grav2);
    scaling = 4.0 * M_PI * 1.0e3 * E->data.layer_km * E->data.grav_const
        / E->data.grav_acc * (E->X[1][E->mesh.nno] -  E->X[1][1]);
    sphere_expansion(E, E->slice.tpg, tpgt[0], tpgt[1]);
    sphere_expansion(E, E->slice.tpgb, tpgb[0], tpgb[1]);
    fprintf(stderr, "TopoScaling %.4e  %.4e %.4e\n",stress_scaling,topo_scaling1,topo_scaling2 );

  
    for (j=0; j<2; j++){
        for (i=0; i<maxindice; i++) {
            //fprintf(stderr,"Topo Scaling before %d  %.4e %.4e\n",i, tpgt[j][i],tpgb[j][i]);

            tpgt[j][i] *= topo_scaling1;
            tpgb[j][i] *= topo_scaling2; 
            //fprintf(stderr,"Topo Scaling %d  %.4e %.4e\n",i, tpgt[j][i],tpgb[j][i]);

        }
   }
return;
}


void geoid_from_buoyancy(struct All_variables *E,
double *harm_geoid[2], double *harm_geoidb[2])
{
    int m,k,ll,mm,node,i,j,p,noz,snode,nxnz;
    double *TT,radius,*geoid[2],dlayer,con1,grav,scaling2,scaling,radius_m;
    double cont, conb;
    double buoy2rho;
    int maxindice;
    void sphere_expansion();

    maxindice = (E->control.llmax); 
    radius_m = E->data.layer_km*1e3;
    /* scale for buoyancy */
    scaling2 = -E->data.therm_exp*E->data.ref_temperature*E->data.density
        / E->control.Ra_temp;
    /* scale for geoid */
    scaling = 4.0 * M_PI * 1.0e3 * E->data.layer_km * E->data.grav_const
        / E->data.grav_acc * (E->X[1][E->mesh.nno] -  E->X[1][1]);
    /* density of one layer */
    TT = (double *) malloc ((E->mesh.nox+1)*sizeof(double));
     
    geoid[0] = (double*)malloc((maxindice+1)*sizeof(double));
    geoid[1] = (double*)malloc((maxindice+1)*sizeof(double));
    for (p = 0; p < maxindice; p++) {
        harm_geoid[0][p] = 0;
        harm_geoid[1][p] = 0;
        harm_geoidb[0][p] = 0;
        harm_geoidb[1][p] = 0;
    }

    for(k=1;k<E->mesh.noz;k++)  {
        grav = 1.0;
        buoy2rho = scaling2 / grav;
        for(j=1;j<=E->mesh.nox;j++)  {
           node = k  +  (j-1)*E->mesh.noz;
           TT[j] = (E->buoyancy[node]+E->buoyancy[node+1])
                 * 0.5 * buoy2rho;

        }
        sphere_expansion(E,TT,geoid[0],geoid[1]);
        dlayer = (E->X[2][k+1]-E->X[2][k])*radius_m;
        radius = (1.0-(E->X[2][k+1]+E->X[2][k])*0.5) / (E->X[1][E->mesh.nno] -  E->X[1][1]);
        for (ll=0;ll<E->control.llmax;ll++) {
            con1 = scaling * dlayer / (2.0*(double)(ll+1));
            cont = exp(-(double)(ll+1)*radius);
            conb = 1.0/buoy2rho;
            //if(ll==0)
              // fprintf(stderr,"scaling %.4e %.4e  %.4e %.4e\n",scaling, dlayer,con1,cont);
         //   fprintf(stderr,"Buoyancy%d %.4e %.4e %.4e\n",ll,cont,geoid[0][ll],geoid[1][ll]);
            harm_geoid[0][ll] += con1*cont*geoid[0][ll];
            harm_geoid[1][ll] += con1*cont*geoid[1][ll];
            harm_geoidb[0][ll] += con1*conb*geoid[0][ll];
            harm_geoidb[1][ll] += con1*conb*geoid[1][ll];
            //if(ll==0)
            //fprintf(stderr,"Geoid buoyancy%d  %.4e %.4e\n",ll, geoid[0][ll],con1*cont*geoid[0][ll]);
        }

    }

    free ((void *)TT);
    free ((void *)geoid[0]);
    free ((void *)geoid[1]);


return;

}

void geoid_from_topography(struct All_variables *E,
                                  double *tpgt[2],
                                  double *tpgb[2],
                                  double *geoid_tpgt[2],
                                  double *geoid_tpgb[2])
{
    void sphere_expansion();
 
    double con1,con2,scaling,den_contrast1,den_contrast2;
    int i,j,k,ll,mm,s;
    int maxindice;
    maxindice = (E->control.llmax);
    den_contrast1 = E->data.density - E->data.density_above;
    den_contrast2 = E->data.density_below - E->data.density;
    for (i = 0; i < maxindice; i++) {
        geoid_tpgt[0][i] = 0;
        geoid_tpgt[1][i] = 0;
        geoid_tpgb[0][i] = 0;
        geoid_tpgb[1][i] = 0;
    }

    // note here *(E->X[1][E->mesh.nno] -  E->X[1][1]) because the x coordinate is (E->X[1][E->mesh.nno] -  E->X[1][1]) times d, and we expand that into (E->X[1][E->mesh.nno] -  E->X[1][1])d
    scaling = 4.0 * M_PI * 1.0e3 * E->data.layer_km * E->data.grav_const
            / E->data.grav_acc *(E->X[1][E->mesh.nno] -  E->X[1][1]) ;
    fprintf(stderr, "TopoGeoidScaling %.4e  %.4e %.4e\n",scaling,scaling*den_contrast1,scaling*den_contrast2 );

    for (j=0; j<2; j++) {
        for (ll=0; ll<maxindice; ll++)   {
            con1 = den_contrast1 * scaling / (2.0*(double)(ll+1));
            con2 = den_contrast2 * scaling / (2.0*(double)(ll+1)) * exp(-1.0/(E->X[1][E->mesh.nno] -  E->X[1][1])*(double)(ll+1));
            geoid_tpgt[j][ll] = tpgt[j][ll] * con1;
            geoid_tpgb[j][ll] = tpgb[j][ll] * con2;
        }
    }
 
return;
}


static void geoid_from_topography_self_g(struct All_variables *E,
                                         double *tpgt[2],
                                         double *tpgb[2],
                                         double *geoid_bncy[2],
                                         double *geoid_bncy_botm[2],
                                         double *geoid_tpgt[2],
                                         double *geoid_tpgb[2])
{
    void sphere_expansion();

    double den_contrast1,den_contrast2,grav1,grav2;
    double topo2stress1, topo2stress2;
    long double con4, ri;
    long double a1,b1,c1_0,c1_1,a2,b2,c2_0,c2_1,a11,a12,a21,a22,f1_0,f2_0,f1_1,f2_1,denom;
    int i,j,k,ll,mm,s;
    int maxindice;
    maxindice = (E->control.llmax);
    den_contrast1 = E->data.density - E->data.density_above;
    den_contrast2 = E->data.density_below - E->data.density;
    topo2stress1 = den_contrast1 * E->data.grav_acc;
    topo2stress2 = den_contrast2 * E->data.grav_acc;

    con4 = 4.0*M_PI*E->data.grav_const*E->data.layer_km*1000.0;
    for (i = 0; i < maxindice; i++) {
        geoid_tpgt[0][i] = 0;
        geoid_tpgt[1][i] = 0;
        geoid_tpgb[0][i] = 0;
        geoid_tpgb[1][i] = 0;
   }
    for (ll=1;ll<=E->control.llmax;ll++)   {
        a1 = 0.0;
   
    }
return;

}



void compute_geoid(E) 
     struct All_variables *E;
{
    void geoid_from_buoyancy();
    void expand_topo_sph_harm();
    void geoid_from_topography();
    int i, p;
    int maxindice;
    maxindice = (E->control.llmax);
    geoid_from_buoyancy(E, E->harm_geoid_from_bncy,
                        E->harm_geoid_from_bncy_botm);
    expand_topo_sph_harm(E, E->harm_tpgt, E->harm_tpgb);

    geoid_from_topography(E, E->harm_tpgt, E->harm_tpgb,
                             E->harm_geoid_from_tpgt,
                             E->harm_geoid_from_tpgb);
    for (i = 0; i < 2; i++) {
        for (p = 0; p < maxindice; p++) {
            E->harm_geoid[i][p]
                = E->harm_geoid_from_bncy[i][p]
                  + E->harm_geoid_from_tpgt[i][p]
                  + E->harm_geoid_from_tpgb[i][p];
        }
    }
}
