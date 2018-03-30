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
  

  //initialize error histogram
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++)
    histo[nb_lines-1]=0;

    
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++) {
    for(i=0;i<nb_iter;i++) {
      
      //create dataset
      ds=new_empty_dataset(2);
      for(j=0;j<nb_lines;j++)
        add_random_line_dataset(ds,20,size,pmin,pmax);
      noise_dataset(ds,noise_level);
      
      //exec the multifit
      result=multifit(ds,convcrit,scalecrit,0,0,NULL);
      remove_degenerated_objects_gem(result,3);
      
      //check for error
      if (result->nb_objects!=nb_lines) {
        histo[nb_lines-1]++;
      }
      
      //free data
      free_gem(result);
      free_dataset(ds);
    }
  }
 

  //save result 
  f=fopen("/tmp/perf_algo.txt","w");
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++) {
    pc=(float)histo[nb_lines-1]/(float)nb_iter*100.0;
    fprintf(f,"%d %d %f\n",nb_lines,histo[nb_lines-1],pc);
    printf("%d %d %f\n",nb_lines,histo[nb_lines-1],pc);
  }
  fclose(f);


  return 1;
}

