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
       sscanf(argv[4],"%s",runfile);
       sscanf(argv[5],"%d",&STEPS); 
       break;
     default:
       print_and_exit("Usage: %s nx NY KAPPA runfile STEPS\n",
           argv[0]);
  }

  if(NULL==(file=fopen(runfile,"r")))
        print_and_exit("I could not open file with simulation run numbers %s\n",runfile);
  
  FRAMES = STEPS/PERIOD;
  EPSILON = 720.0*KAPPA;

  FILE *fp,*rf;
  char filepath[256],init_strip[256],trajectory_file[256],analyzedruns_file[256];
  double dhe,bhe;
  double backbone_T0,slider_T0;;
  int runnum;

  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",nx,NY);
  printf("Init_strip.gsd : %s\n",init_strip);

  //File containing filepaths of analyzed log files
  sprintf(analyzedruns_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/valid_slider_runs.log",nx,NY,KAPPA);

  rf = fopen(analyzedruns_file, "w");
  if (rf == NULL)
   {
	print_and_exit("Could Not Open File to write validruns.log");
   }

//  for(int run=1;run<=RUN;run++)
  while (fscanf(file, "%d", &runnum) == 1)// 1 is returned if fscanf reads a number
  {
	  // Output filepath 
	  sprintf(filepath,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/analyzeSlider.log",nx,NY,KAPPA,runnum);
	  printf("Filename of analyzed data: %s\n",filepath);
	  
	  // Trajectory.gsd filepath
	  sprintf(trajectory_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/traj.gsd",nx,NY,KAPPA,runnum);
	  printf("Trajectory File : %s\n",trajectory_file);

	  fp = fopen(filepath, "w");
	  if (fp == NULL)
	   {
		print_and_exit("Could Not Open File to write analyzed data");
	   }

	  fprintf(rf,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d\n",nx,NY,KAPPA,runnum);
	  
	  /*	T=0 evaluations		*/
	  load_gsd(init_strip,0);
	  backbone_T0 = backbone_length(0);
	  //backbone_length(0,fp);
	  slider_T0 = avg_slider_pos();
	  
	  initialize1();
	  initialize3();
	  
	  fprintf(fp,"#Timestep\tDihedral_Bending_Energy\tBond_Harmonic_Energy\tPotential_Energy\tBackbone\tAvg_Slider_Pos\tDelta_Slider/Avg_Slider_Pos\n");  
	  fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",0,bending_energy(),bond_harmonic_energy(),bending_energy()+bond_harmonic_energy(),backbone_length(0),avg_slider_pos(),(slider_T0-avg_slider_pos())/avg_slider_pos());

	  for(int frames=1;frames<FRAMES;frames++)
	  {
		load_gsd(trajectory_file,frames);
		dhe = bending_energy();
		bhe = bond_harmonic_energy();
		
		fprintf(fp,"%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",frames*PERIOD,dhe,bhe,dhe+bhe,backbone_length(frames),avg_slider_pos(),(slider_T0-avg_slider_pos())/avg_slider_pos());
  	  }

  	fclose(fp);
  }

  fclose(rf);
  fclose(file);

  return 0;
}
