#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"


int nx,NY,LEN,RUN;
double KAPPA,EPSILON;
int STEPS,FRAMES;
double cnode[MAXRUN][MAXFRAMES];


int main(int argc, char **argv)
{

   //File containing the valid runs for which Analysis will be performed
   char runfile[256];
   FILE *file;

   switch (argc){
     case 6:
       sscanf(argv[1],"%d",&nx);    
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%lf",&KAPPA);
       sscanf(argv[4],"%s",&runfile);
       sscanf(argv[5],"%d",&STEPS); 
       break;
     default:
       print_and_exit("Usage: %s nx NY KAPPA runfile STEPS\n",argv[0]);
  }
  
  if(NULL==(file=fopen(runfile,"r")))
        print_and_exit("I could not open file with simulation run numbers %s\n",runfile);

  FRAMES = STEPS/PERIOD;
  EPSILON = 720.0*KAPPA;

  FILE *fp,*hgt,*wid,*bb,*cn,*rf,*vf;
  char filepath[256],init_strip[256],trajectory_file[256],hgt_profile_file[256],hgt_width_file[256],hgt_bb_file[256],cnode_file[256],analyzedruns_file[256],numvalidruns_file[256];
  double dhe,bhe;
  double backbone_T0,slider_T0;;
  int frame_cnt=0,runnum,r;

  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",nx,NY);
  printf("Init_strip.gsd : %s\n",init_strip);

  // Avg, Height Squared ribbon profile path
  sprintf(hgt_profile_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/hgt_prof_real.dat",nx,NY,KAPPA);
  printf("Height Profile File: %s\n",hgt_profile_file);

  //Time series of central node of the ribbon
  sprintf(cnode_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/cnode.bin",nx,NY,KAPPA);
  printf("Central node time series File: %s\n",cnode_file);

  hgt = fopen(hgt_profile_file, "w");
  if (hgt == NULL)
   {
	print_and_exit("Could Not Open File to write height profile data");
   }

  cn = fopen(cnode_file, "wb");
  if (cn == NULL)
  {
  	print_and_exit("Could Not Open File to write central node time series");
  } 

  /* Initializing the arrays	*/
  //initialize();

  int c; //counter for frames inside each run 

  //File containing filepaths of runs
  sprintf(analyzedruns_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/validruns.log",nx,NY,KAPPA);

  rf = fopen(analyzedruns_file, "w");
  if (rf == NULL)
   {
        print_and_exit("Could Not Open File to write validruns.log");
   }

  //File containing number of valid runs
  sprintf(numvalidruns_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/numvalidruns.log",nx,NY,KAPPA);
  vf = fopen(numvalidruns_file, "w");
  if (vf == NULL)
  {
        print_and_exit("Could Not Open File to write numvalidruns.log");
  }

  r=0;
  while (fscanf(file, "%d", &runnum) == 1)// 1 is returned if fscanf reads a number
  {

	  // Output filepath 
	  sprintf(filepath,"../Sim_dump_stretched/L%d/W%d/k%.1f/r%d/analyze.log",nx,NY,KAPPA,runnum);
	  printf("Filename of analyzed data: %s\n",filepath);
	  
	  // Trajectory.gsd filepath
	  sprintf(trajectory_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/r%d/traj_stretched.gsd",nx,NY,KAPPA,runnum);
	  printf("Trajectory File : %s\n",trajectory_file);

	  //Avg Width height of the ribbon
	  sprintf(hgt_width_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/r%d/hgt_widthavg.bin",nx,NY,KAPPA,runnum);
	  printf("Height width File: %s\n",hgt_width_file);

	  //Height of the ribbon backbone
          sprintf(hgt_bb_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/r%d/backbone.bin",nx,NY,KAPPA,runnum);
          printf("Backbone Height File: %s\n",hgt_bb_file);

	  fp = fopen(filepath, "w");
	  if (fp == NULL)
	   {
		print_and_exit("Could Not Open File to write analyzed data");
	   }

	  wid = fopen(hgt_width_file, "wb");
  	  if (wid == NULL)
   	  {
        	print_and_exit("Could Not Open File to write height width data");
   	  }

	  bb = fopen(hgt_bb_file, "wb");
          if (wid == NULL)
          {
                print_and_exit("Could Not Open File to write backbone height data");
          }
	  
	  /*	T=0 evaluations		*/
	  load_gsd(init_strip,0);
	  backbone_T0 = backbone_length(0);
	
          c=0;//count of frames > FRAMES/2
	  initialize1();
	  initialize3();
	  
	  slider_T0 = avg_slider_pos();
	  fprintf(fp,"Frames\tDihedral_Bending_Energy\tBond_Harmonic_Energy\tPotential_Energy\tDelta_Backbone\tAvg_hgt\tAvg_hgt_Sq\tAvg_Slider_Pos\tDelta_Slider\n");  
	  fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",0,bending_energy(),bond_harmonic_energy(),bending_energy()+bond_harmonic_energy(),backbone_length(0),avg_hgt(),avg_hgt_sq(),avg_slider_pos(),(slider_T0-avg_slider_pos())/avg_slider_pos());

	  for(int frames=1;frames<FRAMES;frames++)
	  {
		load_gsd(trajectory_file,frames);
		dhe = bending_energy();
		bhe = bond_harmonic_energy();
		
		fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",frames,dhe,bhe,dhe+bhe,backbone_length(frames)-backbone_T0,avg_hgt(),avg_hgt_sq(),avg_slider_pos(),(slider_T0-avg_slider_pos())/avg_slider_pos());
		// Complete time series of central node
		cnode[r][frames-1]=position[3*((N+1)/2)+2];

		// Height of the ribbon averaged over the width
                width_hgt(frames-1);

                // Height of ribbon backbone
                bb_hgt(frames-1);


		// Average Height of each node in the last half of the simulation
		if(frames>=FRAMES/2)
		{
			frame_cnt++;//counting over all runs no reset required
			sum_hgt_node();
		}
	  }

	  r++; //counting the number of valid runs

	  print_width(wid,FRAMES); 
	  print_bb(bb,FRAMES);
	  fclose(fp);
	  fclose(wid);
	  fclose(bb);
  }

  //Number of valid runs written to file
  fprintf(vf,"%d\n",r);

  fclose(rf);
  fclose(vf);
  fclose(file);
  
  /*    writing central node height time series         */
  fwrite(&r,sizeof(int),1,cn);//Number of runs
  for(int i=0;i<r;i++)
  {
        for(int j=0;j<FRAMES;j++)
        {
                fwrite(&cnode[i][j],sizeof(double),1,cn);
        }
  }
  fclose(cn);

  // Output files are generated only for the valid runs

  //Average Height of each node (averaged over last half of the frames)
  avg_hgt_node(frame_cnt);

  if(NULL==(file=fopen(runfile,"r")))
        print_and_exit("I could not open file with simulation run numbers %s\n",runfile);

  initialize2(); // Initializing the hgt_fluctuation array 
  frame_cnt=0;

  while (fscanf(file, "%d", &runnum) == 1)// 1 is returned if fscanf reads a number
  {
	// Trajectory.gsd filepath
	sprintf(trajectory_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/r%d/traj_stretched.gsd",nx,NY,KAPPA,runnum);
	for(int frames=FRAMES/2;frames<FRAMES;frames++)
	{
		load_gsd(trajectory_file,frames);
		//Height fluctuation profile of the ribbon
		hgt_profile();
		frame_cnt++;
	}
  }

  fclose(file);

  //Average Height Fluctuation 
  avg_hgt_profile(hgt,frame_cnt);

  fclose(hgt);

  return 0;
}
