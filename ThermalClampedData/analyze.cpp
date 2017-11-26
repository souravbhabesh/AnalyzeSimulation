#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

double bendingEner[NMAX];
double bondHarmonicEner[NMAX];
double backbone;
double total_DHE,total_BHE;
double hgt_fluctuation[NMAX]; 
double h_avg_node[NMAX];
double h_width[MAXFRAMES][NXMAX];
double h_bb[MAXFRAMES][NXMAX];

double bond_length (int i,int j)
{
  return (sqrt(pow(position[3*i]-position[3*j],2)+pow(position[3*i+1]-position[3*j+1],2)+\
					pow(position[3*i+2]-position[3*j+2],2)));
}

double backbone_length(int frame)
{
  int node;
  backbone=0;
  for(int i=1;i<nx;i++)
  {
	node = (NY/2)*nx + i;
	backbone += bond_length(node,node-1);
	//printf("%d\t%d\t%lf\t%.8f\t%.8f\n",node,node-1,bond_length(node,node-1),position[3*node],position[3*(node-1)]);	  
  }
  //fprintf (fp,"%d\t%lf\n",frame,backbone);
  return (backbone);
}


/*	Cross product	*/
double* cross_product(double u[3],double v[3])
{
  static double u_cross_v[3];
  u_cross_v[0] = u[1]*v[2] - u[2]*v[1];
  u_cross_v[1] = u[2]*v[0] - u[0]*v[2];
  u_cross_v[2] = u[0]*v[1] - u[1]*v[0];

  double mod_u_cross_v = sqrt(u_cross_v[0]*u_cross_v[0] + u_cross_v[1]*u_cross_v[1] + u_cross_v[2]*u_cross_v[2]);

  u_cross_v[0] = u_cross_v[0]/mod_u_cross_v;
  u_cross_v[1] = u_cross_v[1]/mod_u_cross_v;
  u_cross_v[2] = u_cross_v[2]/mod_u_cross_v;

  //printf("CP:\t%lf,%lf,%lf\n",u_cross_v[0],u_cross_v[1],u_cross_v[2]);
  return u_cross_v ;
}


/*	Function evaluating Dihedral Harmonic Energy	*/
double bending_energy()
{
  double vec_cb[3],vec_ab[3],vec_dc[3];
  double *A,*B,V_A[3],V_B[3];
  double be,dot_AB;
  double total_DHE = 0;

  for(int i=0;i<Nd;i++)
  {
        for(int j=0;j<3;j++)
        {
                vec_cb[j] = position[3*dihedralGroup[4*i+2]+j] - position[3*dihedralGroup[4*i+1]+j];
                vec_ab[j] = position[3*dihedralGroup[4*i]+j] - position[3*dihedralGroup[4*i+1]+j];
                vec_dc[j] = position[3*dihedralGroup[4*i+3]+j] - position[3*dihedralGroup[4*i+2]+j];
        }
        //printf ("Dihedral %d:\t%d %d %d %d\n",i,dihedralGroup[4*i],dihedralGroup[4*i+1],dihedralGroup[4*i+2],dihedralGroup[4*i+3]);
        //printf("vec_cb:\t%lf,%lf,%lf\n",vec_cb[0],vec_cb[1],vec_cb[2]);
        //printf("vec_ab:\t%lf,%lf,%lf\n",vec_ab[0],vec_ab[1],vec_ab[2]);
        //printf("vec_dc:\t%lf,%lf,%lf\n",vec_dc[0],vec_dc[1],vec_dc[2]);

        A = cross_product(vec_cb,vec_ab);
	for(int k=0;k<3;k++)
	{
		V_A[k]=*(A+k);		
	}
        B = cross_product(vec_cb,vec_dc);
	for(int k=0;k<3;k++)
        {
                V_B[k]=*(B+k);
        }
	//printf("vec_cb X vec_ab : %.8f,%.8f,%.8f\n",V_A[0],V_A[1],V_A[2]);
	//printf("vec_cb X vec_dc : %.8f,%.8f,%.8f\n",V_B[0],V_B[1],V_B[2]);

        dot_AB = V_A[0]*V_B[0]+V_A[1]*V_B[1]+V_A[2]*V_B[2];
        //printf("dot_AB = %lf\n",dot_AB);
        be = 0.5 * KAPPA * (1+dot_AB);//Using HOOMD kappa

	total_DHE += be;
        //printf("BE = %lf\n",be);

  }
  return(total_DHE);
}

/*	Function evaluating Bond Harmonic Energy	*/
double bond_harmonic_energy()
{
  total_BHE = 0;
  double l;//current length of bond
  double se;
  for(int i=0;i<Nb;i++)
  {
	l=0;
  	for(int j=0;j<3;j++)
  	{
		l = l + (position[3*bondGroup[2*i]+j] - position[3*bondGroup[2*i+1]+j]) * (position[3*bondGroup[2*i]+j] - position[3*bondGroup[2*i+1]+j]);
	}
	l = sqrt(l);
	se = 0.5 * EPSILON * (l-a) * (l-a);
	//printf("Bond %d %d , se = %lf\n",bondGroup[2*i],bondGroup[2*i+1],se);
	total_BHE = total_BHE + se;
  } 
  return (total_BHE);
}

/*	Function evaluating Average Z height above z=0 plane	*/
/*	particle_id=1 clamped points on the left		*/
/*	particle_id=0 Normal Lattice Sites			*/
/*	particle_id=3 Nodes on the right constrained to X	*/
/*	particle_id=4 Backbone of the ribbon excluding two	*/
/*		 lattice sites at each boundary			*/

double avg_hgt()
{
   double hgt=0;
   int node_cnt=0;
   for(int i=0;i<N;i++)
   {
	//if(particleID[i]==0 || particleID[i]==4)
	//{
	hgt+=position[3*i+2];
	node_cnt++;
	//}
   }
   return (hgt/node_cnt);
}

/*	Average <h^2> = 1/N * Sum_i(z_i - <z>)^2	*/
double avg_hgt_sq()
{
   double hgtSq=0;
   double h_avg = avg_hgt();
   int node_cnt=0;
   for(int i=0;i<N;i++)
   {
	//if(particleID[i]==0 || particleID[i]==4)
	//{
        	hgtSq+=pow((position[3*i+2]-h_avg),2);
		node_cnt++;
	//}
   }
   return (hgtSq/node_cnt);
   
}

/*	Slider average position	*/
double avg_slider_pos()
{
   double slider_pos=0;
   int slider_node=0;
   for(int i=0;i<N;i++)
   {
	if(particleID[i]==3)
	{
		slider_pos+=position[3*i];
		slider_node++;
	}
   }
   return (slider_pos/slider_node);
}

/*	Initialize arrays	*/
int initialize1()
{
   for(int i=0;i<N;i++)
   {
	h_avg_node[i]=0;
   }
   return 0;
}

int initialize2()
{
   for(int i=0;i<N;i++)
   {
        hgt_fluctuation[i]=0;
   }
   return 0;
}

int initialize3()
{
   for(int i=0;i<FRAMES;i++)
   {
	for(int j=0;j<nx;j++)
	{
   	   h_width[i][j]=0;
	}
	for(int k=0;k<nx;k++)
        {
           h_bb[i][k]=0;
        }
   }
   return 0;
}

/*	Height sum at each node	over frames	*/
int sum_hgt_node()
{
   for(int i=0;i<N;i++)
   {
	h_avg_node[i]+=position[3*i+2];
   }
   return 0;
}

/*	Average Height at each node over all frames from all runs	*/	
int avg_hgt_node(int tot_frames)
{
   if (tot_frames == 0)
   {
        printf("Average Node Hgt computation is dividing by Zero\n");
	return 0;
   }
   for(int i=0;i<N;i++)
   {
	h_avg_node[i]=h_avg_node[i]/tot_frames;
   }
   return 0;
}

/*	Height fluctuation profile of the ribbon	*/
int hgt_profile()
{
   for(int i=0;i<N;i++)
   {
        if(particleID[i]==0 || particleID[i]==4)
                hgt_fluctuation[i]+=pow((position[3*i+2]-h_avg_node[i]),2);
   }
   return 0;
}

/*	Average Height Fluctuation	*/
int avg_hgt_profile(FILE *hgt,int tot_frames)
{
   if (tot_frames == 0)
   {
        printf("Average Hgt Fluctuation computation is dividing by Zero\n");
	return 0;
   }
   printf("Total Frames in last half of Simulation over all runs: %d\n",tot_frames);
   for(int i=0;i<N;i++)
   {
        if(particleID[i]==0 || particleID[i]==4)
                hgt_fluctuation[i]=hgt_fluctuation[i]/tot_frames;
        if(i%nx!=nx-1)
		fprintf(hgt,"%.8f ",sqrt(hgt_fluctuation[i]));
	else
		fprintf(hgt,"%.8f\n",sqrt(hgt_fluctuation[i]));
	
   }
   return 0;
}

/*	Backbone Height average		*/
int bb_hgt(int frame)
{
 int j=0; //count over the backbone nodes
 h_bb[frame][0] = 0; //first two nodes are fixed 
 h_bb[frame][1] = 0;
 j = 2;
 for(int i=0;i<N;i++)
   {
        if(particleID[i]==4)
        {
		h_bb[frame][j] = position[3*i+2];
		j++;
        }
   }
 h_bb[frame][j] = 0;
 h_bb[frame][j+1] = 0; 
/*
 if(frame == 0)
 {
	for(int i=0;i<nx;i++)
	{
		printf("%d\t%.8f\n",i,h_bb[frame][i]);
	}
 }
*/
 return 0;
}

/*	Width Height Average	*/
int width_hgt(int frame)
{
  int k=0,k_cnt;
  for(int i=0;i<nx;i++)
  {
	do {
	   //if (frame==0)
		//printf("%d\t%d\t%d\t%.8f\n",N,i,(i/2)+2*k*nx,position[3*((i/2)+2*k*nx)+2]);
	   h_width[frame][i] += position[3*(i+2*k*nx)+2];
	   k++;
	}while(((i/2)+2*k*nx) < N);
	k_cnt = k;
	k=0;
  }

  for(int i=0;i<nx;i++)
  {
	h_width[frame][i] = h_width[frame][i]/k_cnt;
	//if(frame == 0)
		//printf ("%d\t%.8f\t%d\t%d\n",i,h_width[frame][i],k_cnt,j_cnt);
  }
  return 0;
}

/*	Printing Width Height Average for all frames	*/
int print_width(FILE *wid, int total_frames)
{
  //fwrite(h_width, sizeof(double),MAXFRAMES*NXMAX,wid);//FRAMES*nx,wid);
  fwrite(&nx,sizeof(int),1,wid);
  fwrite(&total_frames,sizeof(int),1,wid);
  for(int i=0;i<total_frames;i++)
  {
        for(int j=0;j<nx;j++)
        {
                fwrite(&h_width[i][j],sizeof(double),1,wid);
        }
  }
   printf("#Frames in last half of simulation : %d\n\n",FRAMES/2);

  return 0;
}

/*	Printing Backbone height for all frames	*/
int print_bb(FILE *bb,int total_frames)
{
  //fwrite(h_bb, sizeof(double),MAXFRAMES*NXMAX,bb);//FRAMES*nx/2,bb);
  fwrite(&nx,sizeof(int),1,bb);
  fwrite(&total_frames,sizeof(int),1,bb);
  for(int i=0;i<total_frames;i++)
  {
        for(int j=0;j<nx;j++)
        {
                fwrite(&h_bb[i][j],sizeof(double),1,bb);
        }
  }
  return 0;
}

