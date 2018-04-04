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
  FILE *fc;
  FILE *fcp;
  FILE *fcpd;
  int i,j,k;
  int nb_lines;
  int nb_iter;
  double scalecrit;
  struct point *pmin;
  struct point *pmax;
  int stop;
  float abl;
  float angres;

  //graphical debug environment
  #ifdef DEBUG
  struct graphics *gws;
  struct point *gmin;
  struct point *gmax;
#endif

  
  //error histograms
  #define NB_BINS 8
  #define MAX_SIZE NB_BINS+5
  int histo_wc[NB_BINS][MAX_SIZE];
  int histo_wcp[NB_BINS][MAX_SIZE];
  int histo_wcpd[NB_BINS][MAX_SIZE];
  struct line *targets[NB_BINS];

  
  //parameters
  convcrit=0.00001;
  noise_level=0.1;
  size=80;
  scalecrit=150;
  nb_iter=2000;
  pmin=new_2D_point(-50,-50);
  pmax=new_2D_point(50,50);
  angres=7.0;

  //initialize random generator
  srand(0);
  
  //initialize graphical debug environment
#ifdef DEBUG
  gmin=new_2D_point(-200,-200);
  gmax=new_2D_point(200,200);
  gws=new_graphics(2,1000,1000,gmin,gmax,param2verb("summary"));
#endif
  
  //initialize error histograms
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++)
    for(i=0;i<=MAX_SIZE;i++) {
      histo_wc[nb_lines-1][i]=0;
      histo_wcp[nb_lines-1][i]=0;
      histo_wcpd[nb_lines-1][i]=0;
    }

  //open the files
  fc=fopen("/tmp/perf_algo_wc.txt","w");
  fcp=fopen("/tmp/perf_algo_wcp.txt","w");
  fcpd=fopen("/tmp/perf_algo_wcpd.txt","w");
  
  //iterate on mixtures
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++) {
    for(i=0;i<nb_iter;i++) {
      printf("nb_lines %d iter %d\n",nb_lines,i);
      
      //create dataset : nb_lines without colinearity at more that angres
      ds=new_empty_dataset(2);
      for(j=0;j<nb_lines;j++) {
        do {
          targets[j]=random_line(2,pmin,pmax);
          stop=1;
          for(k=0;k<j;k++) {
            abl=angle_between_line(targets[k],targets[j]);
            if ((abl<angres) || (abl>180-angres)) {
              //printf("abl=%f\n",abl);
              stop=0;
            }
          }
          if (stop==0)
            free_line(targets[j]);
        } while(stop==0);
        add_line_dataset(ds,targets[j],20,size);
      }
      noise_dataset(ds,noise_level);
      
      //exec the multifit
      result=multifit(ds,convcrit,scalecrit,0,0,NULL);

      //fill histogram
      if (result->nb_objects>MAX_SIZE)
	histo_wc[nb_lines-1][MAX_SIZE]++;
      else
	histo_wc[nb_lines-1][result->nb_objects]++;

      //apply degenerated objects post-processing
      remove_degenerated_objects_gem(result,3);

      //fill histogram
      if (result->nb_objects>MAX_SIZE)
	histo_wcp[nb_lines-1][MAX_SIZE]++;
      else
	histo_wcp[nb_lines-1][result->nb_objects]++;

      //apply duplicated objects post-processing
      remove_dup_objects_gem(result,0.1);

      //fill histogram
      if (result->nb_objects>MAX_SIZE)
	histo_wcpd[nb_lines-1][MAX_SIZE]++;
      else
	histo_wcpd[nb_lines-1][result->nb_objects]++;


      //print if debug
#ifdef DEBUG
      if (result->nb_objects!=nb_lines) {
        plot_gem(result,gws);
        getchar();
      }
#endif
      
      //free data
      free_gem(result);
      free_dataset(ds);
      for(j=0;j<nb_lines;j++) {
        free_line(targets[j]);
      }
    }

    //dump result in files
    fprintf(fc,"%d ",nb_lines);
    fprintf(fcp,"%d ",nb_lines);
    fprintf(fcpd,"%d ",nb_lines);
    for(i=0;i<=MAX_SIZE;i++) {
      fprintf(fc,"%f ",(float)histo_wc[nb_lines-1][i]/(float)nb_iter*100.0);
      fprintf(fcp,"%f ",(float)histo_wcp[nb_lines-1][i]/(float)nb_iter*100.0);
      fprintf(fcpd,"%f ",(float)histo_wcpd[nb_lines-1][i]/(float)nb_iter*100.0);
    }
    fprintf(fc,"%f\n",100.0-(float)histo_wc[nb_lines-1][nb_lines]/(float)nb_iter*100.0);
    fprintf(fcp,"%f\n",100.0-(float)histo_wcp[nb_lines-1][nb_lines]/(float)nb_iter*100.0);
    fprintf(fcpd,"%f\n",100.0-(float)histo_wcpd[nb_lines-1][nb_lines]/(float)nb_iter*100.0);
  }
 
  fclose(fc);
  fclose(fcp);
  fclose(fcpd);
  return 1;
}

