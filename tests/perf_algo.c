/*
Copyright 2016-2018 Frédéric Magniette
This file is part of Libgem.

Libgem is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Libgem is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Libgem.  If not, see <http://www.gnu.org/licenses/>
*/

#include <libgem.h>
#include <sys/time.h>

//check performances

int main(int argc,char **argv) {

  struct gem_ws *result;
  struct dataset *ds;
  double convcrit;
  double noise_level;
  int size;
  FILE *f;
  int i,j;
  int nb_lines;
  int nb_iter;
  int seed;
  double scalecrit;
  double pc;
  struct point *pmin;
  struct point *pmax;
  float uerr;
  float serr;
  //struct graphics *gws;
  //struct point *gmin;
  //struct point *gmax;

  //error histograms
  #define NB_BINS 8
  int histo[NB_BINS];
  
  if (argc!=3) {
    printf("usage : %s seed convcrit\n",argv[0]);
    return -1;
  }

  
  //parameters
  seed=atoi(argv[1]);
  convcrit=atof(argv[2]);
  srand(seed);
  noise_level=0.1;
  size=100;
  scalecrit=150;
  nb_iter=5000;
  pmin=new_2D_point(-50,-50);
  pmax=new_2D_point(50,50);
  //gmin=new_2D_point(-200,-200);
  //gmax=new_2D_point(200,200);
  //gws=new_graphics(2,1000,1000,gmin,gmax,param2verb("summary"));
  f=fopen("/tmp/perf_algo.txt","w");
  uerr=0.0;

  //initialize error histogram
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++)
    histo[nb_lines-1]=0;

    
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++) {
    printf("nb_lines=%d\n",nb_lines);
    for(i=0;i<nb_iter;i++) {
      //create dataset
      ds=new_empty_dataset(2);
      for(j=0;j<nb_lines;j++)
        add_random_line_dataset(ds,20,size,pmin,pmax);
      noise_dataset(ds,noise_level);
      
      //exec the multifit
      result=multifit(ds,convcrit,scalecrit,0,0,NULL);
      
      //check for error
      if (result->nb_objects!=nb_lines) {
        histo[nb_lines-1]++;
        //plot_gem(result,gws);
        //getchar();
      }
      
      //free data
      free_gem(result);
      free_dataset(ds);
    }

    pc=(float)histo[nb_lines-1]/(float)nb_iter*100.0;
    if (nb_lines==2)
      uerr=pc;
    serr=uerr*((float)nb_lines)*((float)nb_lines-1)/2.0;
    fprintf(f,"%d %d %f %f %f\n",nb_lines,histo[nb_lines-1],pc,serr,pc-serr);
    printf("%d %d %f %f %f\n",nb_lines,histo[nb_lines-1],pc,serr,pc-serr);

  }
 
  fclose(f);
  return 1;
}

