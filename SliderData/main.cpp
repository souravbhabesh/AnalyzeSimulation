#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
//#include "gsd_tools.h"
//#include "gsd_fn.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"


int nx,NY,LEN,RUN;
double KAPPA,EPSILON;
int STEPS,FRAMES;



int main(int argc, char **argv)
{
   switch (argc){
     case 6:
       sscanf(argv[1],"%d",&nx);    
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%lf",&KAPPA);
       sscanf(argv[4],"%d",&RUN);
       sscanf(argv[5],"%d",&STEPS); 
       break;
     default:
       print_and_exit("Usage: %s nx NY KAPPA RUN STEPS\n",
           argv[0]);
  }
  
  FRAMES = STEPS/PERIOD;
  EPSILON = 720.0*KAPPA;

  FILE *fp,*hgt,*wid,*bb;
  char filepath[256],init_strip[256],trajectory_file[256],hgt_profile_file[256],hgt_width_file[256],hgt_bb_file[256];
  double dhe,bhe;
  double backbone_T0,slider_T0;;
  int frame_cnt=0;

  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",nx,NY);
  printf("Init_strip.gsd : %s\n",init_strip);

  // Avg, Height Squared ribbon profile path
/*  sprintf(hgt_profile_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/hgt_prof_real.dat",nx,NY,KAPPA);
  printf("Height Profile File: %s\n",hgt_profile_file);

  hgt = fopen(hgt_profile_file, "w");
  if (hgt == NULL)
   {
	print_and_exit("Could Not Open File to write height profile data");
   }
*/
  /* Initializing the arrays	*/
  //initialize();

  for(int run=1;run<=RUN;run++)
  {

	  // Output filepath 
	  sprintf(filepath,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/analyzeSlider.log",nx,NY,KAPPA,run);
	  printf("Filename of analyzed data: %s\n",filepath);
	  
	  // Trajectory.gsd filepath
	  sprintf(trajectory_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/traj.gsd",nx,NY,KAPPA,run);
	  printf("Trajectory File : %s\n",trajectory_file);

/*	  //Avg Width height of the ribbon
	  sprintf(hgt_width_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/width.bin",nx,NY,KAPPA,run);
	  printf("Height width File: %s\n",hgt_width_file);

	  //Height of the ribbon backbone
          sprintf(hgt_bb_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/backbone.bin",nx,NY,KAPPA,run);
          printf("Height width File: %s\n",hgt_bb_file);
*/

	  fp = fopen(filepath, "w");
	  if (fp == NULL)
	   {
		print_and_exit("Could Not Open File to write analyzed data");
	   }

/*	  wid = fopen(hgt_width_file, "wb");
  	  if (wid == NULL)
   	  {
        	print_and_exit("Could Not Open File to write height width data");
   	  }

	  bb = fopen(hgt_bb_file, "wb");
          if (wid == NULL)
          {
                print_and_exit("Could Not Open File to write backbone height data");
          }
*/
	  //printf("Reading GSD file: %s\n",argv[1]);
	  //load_gsd(argv[1],0);
	  
	  /*	T=0 evaluations		*/
	  load_gsd(init_strip,0);
	  backbone_T0 = backbone_length(0);
	  //backbone_length(0,fp);
	  
          int c=0;//count of frames > FRAMES/2
	  initialize1();
	  initialize3();
	  
	  slider_T0 = avg_slider_pos();
	  fprintf(fp,"Frames\tDihedral_Bending_Energy\tBond_Harmonic_Energy\tPotential_Energy\tBackbone\tAvg_hgt\tAvg_hgt_Sq\tAvg_Slider_Pos\tDelta_Slider\n");  
	  fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",0,bending_energy(),bond_harmonic_energy(),bending_energy()+bond_harmonic_energy(),backbone_length(0),avg_hgt(),avg_hgt_sq(),avg_slider_pos(),(slider_T0-avg_slider_pos())/avg_slider_pos());

	  for(int frames=1;frames<FRAMES;frames++)
	  {
		initialize1();
		//load_gsd(argv[2],frames);
		load_gsd(trajectory_file,frames);
		//backbone_length(frames,fp);
		//printf("%d\t%lf\t%lf\n",frames,position[3*(nx-1)],position[3*(LEN-nx)]);
		dhe = bending_energy();
		bhe = bond_harmonic_energy();
		
		fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",frames,dhe,bhe,dhe+bhe,backbone_length(frames),avg_hgt(),avg_hgt_sq(),avg_slider_pos(),(slider_T0-avg_slider_pos())/avg_slider_pos());
		
		if(frames>=FRAMES/2)
		{
			initialize1();
			frame_cnt++;
			sum_hgt_node();
			//width_hgt(c);
			//bb_hgt(c);
			c++;
			//printf("%d\t",frames - (FRAMES/2 + 1));
		}
		//if(frames == FRAMES/2 + 1)
		//width_hgt(0);			
	  }
	  //print_width(wid); 
	  //print_bb(bb);
	  fclose(fp);
	  //fclose(wid);
	  //fclose(bb);
  }
/*
  //Average Height of each node (averaged over last half of the frames)
  avg_hgt_node(frame_cnt);

  frame_cnt=0;
  for(int run=1;run<=RUN;run++)
  {
	// Trajectory.gsd filepath
	sprintf(trajectory_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/traj_thermal.gsd",nx,NY,KAPPA,run);
	initialize2(); // Initializing the hgt_fluctuation array 
	for(int frames=FRAMES/2;frames<FRAMES;frames++)
	{
		load_gsd(trajectory_file,frames);
		//Height fluctuation profile of the ribbon
		hgt_profile();
		frame_cnt++;
	}
  }

  //Average Height Fluctuation 
  avg_hgt_profile(hgt,frame_cnt);

  fclose(hgt);
*/
  return 0;
}
