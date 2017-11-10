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

#include "libgem.h"

//*********************************** STANDARD GAUSSIAN

struct gauss *new_gauss() {
  struct gauss *result=malloc(sizeof(struct gauss));
  result->nb_data=0.0;
  result->mean=0.0;
  result->old_mean=0.0;
  result->variance=0.0;
  result->stdev=0.0;
  result->s=0.0;
  result->old_s=0.0;
  return result;
}

void add_data_gauss(struct gauss *dist,double x) {
  printf("adding %f\n",x);
  dist->nb_data+=1.0;
  printf("nb data=%f\n",dist->nb_data);
  //printf("old mean=%f\n",dist->old_mean);
  //printf("old s=%f\n",dist->old_s);
  if (dist->nb_data==1.0) {
    dist->mean=x;
    dist->s=0;
  } else {
    dist->mean=dist->old_mean+((x-dist->old_mean)/dist->nb_data);
    dist->s=dist->old_s+((x-dist->old_mean)*(x-dist->mean));
  }

  dist->variance=dist->s/dist->nb_data;
  dist->stdev=sqrt(dist->variance);
  
  dist->old_s=dist->s;
  dist->old_mean=dist->mean;
  dist->old_s=dist->s;
  printf("mean=%f\n",dist->mean);
  printf("s=%f\n",dist->s);
}

double likelyhood_gauss(struct gauss *dist,double x) {
  return 1/(dist->stdev*sqrt(2*PI))*exp(-1*pow(x-dist->mean,2))/(2*pow(dist->stdev,2));
}

void dump_gauss(struct gauss *dist,double xmin,double xmax,int nb_steps,char *filename) {
  double x;
  FILE *f=fopen(filename,"w");
  for(x=xmin;x<=xmax;x+=(xmax-xmin)/nb_steps) {
    fprintf(f,"%f %f\n",x,likelyhood_gauss(dist,x));
  }
  fclose(f);
}

//*************************** WEIGHTED GAUSSIAN


struct wgauss *new_wgauss() {
  struct wgauss *result=malloc(sizeof(struct wgauss));
  result->sum_weights=0.0;
  result->mean=0.0;
  result->old_mean=0.0;
  result->s=0.0;
  result->old_s=0.0;
  result->stdev=0.0;
  result->variance=0;
  return result;
}

struct wgauss *cp_wgauss(struct wgauss *dist) {
  struct wgauss *result=malloc(sizeof(struct wgauss));
  result->sum_weights=dist->sum_weights;
  result->mean=dist->mean;
  result->old_mean=dist->old_mean;
  result->s=dist->s;
  result->old_s=dist->old_s;
  result->stdev=dist->stdev;
  result->variance=dist->variance;
  return result;
}

void free_wgauss(struct wgauss *dist) {
  free(dist);
}

void reinit_wgauss(struct wgauss *dist) {
  dist->sum_weights=0.0;
  dist->mean=0.0;
  dist->old_mean=0.0;
  dist->s=0.0;
  dist->old_s=0.0;
  dist->stdev=0.0;
  dist->variance=0;
}

void print_wgauss(struct wgauss *dist) {
  printf("mean=%f stdev=%f\n",dist->mean,dist->stdev);
}


void add_data_wgauss(struct wgauss *dist,double x,double weight) {
  //printf("adding %f with weight %f\n",x,weight);
  //printf("old mean=%f\n",dist->old_mean);
  //printf("old s=%f\n",dist->old_s);

  //calculate recusive mean and variance precursor (s)
  if (dist->sum_weights==0.0) {
    dist->sum_weights=weight;
    dist->mean=x;
    dist->s=0;
  } else {
    dist->sum_weights+=weight;
    dist->mean=dist->old_mean+((x-dist->old_mean)*weight/dist->sum_weights);
    dist->s=dist->old_s+(weight*(x-dist->old_mean)*(x-dist->mean));
  }

  //calculate variance and standard deviation
  dist->variance=dist->s/dist->sum_weights;
  dist->stdev=sqrt(dist->variance);

  dist->old_mean=dist->mean;
  dist->old_s=dist->s;
  //printf("sum weight=%f\n",dist->sum_weights);
  //printf("mean=%f\n",dist->mean);
  //printf("s=%f\n",dist->s);
}

double likelyhood_wgauss(struct wgauss *dist,double x) {
  return 1/(dist->stdev*sqrt(2*PI))*exp(-1*pow(x-dist->mean,2))/(2*pow(dist->stdev,2));
}

double normalized_likelyhood_wgauss(struct wgauss *dist,double x) {
  double result;

  if (dist->stdev==0.0) {
    if (x==dist->mean)
      return 1/sqrt(2*PI);
    else
      return 0.0;
  }
  result=1/sqrt(2*PI)*exp(-1*pow((x-dist->mean)/dist->stdev,2)/2);

  if (isnan(result)) {
    printf("ERROR : likelyhood is nan : \n");
    printf("mean=%f\n",dist->mean);
    printf("stdev=%f\n",dist->stdev);
    printf("x=%f\n",x);
    exit(1);
  }
  return result;
}


void dump_wgauss(struct wgauss *dist,double xmin,double xmax,int nb_steps,char *filename) {
  double x;
  FILE *f=fopen(filename,"w");
  for(x=xmin;x<=xmax;x+=(xmax-xmin)/nb_steps) {
    fprintf(f,"%f %f\n",x,likelyhood_wgauss(dist,x));
  }
  fclose(f);
}

//*************************** CENTERED WEIGHTED GAUSSIAN


struct cwgauss *new_cwgauss() {
  struct cwgauss *result=malloc(sizeof(struct cwgauss));
  result->sum_weights=0.0;
  result->mean=0.0;
  result->old_stdev=0.0;
  result->stdev=0.0;
  return result;
}

struct cwgauss *cp_cwgauss(struct cwgauss *dist) {
  struct cwgauss *result=malloc(sizeof(struct cwgauss));
  result->sum_weights=dist->sum_weights;
  result->mean=dist->mean;
  result->old_stdev=dist->old_stdev;
  result->stdev=dist->stdev;
  return result;
}

void free_cwgauss(struct cwgauss *dist) {
  free(dist);
}

void reinit_cwgauss(struct cwgauss *dist) {
  dist->sum_weights=0.0;
  dist->mean=0.0;
  dist->old_stdev=0.0;
  dist->stdev=0.0;
}

void print_cwgauss(struct cwgauss *dist) {
  printf("mean=%f stdev=%f\n",dist->mean,dist->stdev);
}


void add_data_cwgauss(struct cwgauss *dist,double x,double weight) {

  //calculate recusive mean and variance precursor (s)
  if (dist->sum_weights==0.0) {
    dist->sum_weights=weight;
    dist->mean=0;
    dist->stdev=x*weight;
   } else {
    dist->sum_weights+=weight;
    dist->mean=0;
    dist->stdev=dist->old_stdev+((x-dist->old_stdev)*weight/dist->sum_weights);
    if (dist->stdev<0.001)
      dist->stdev=0.001;
  }
  dist->old_stdev=dist->stdev;
}

double likelyhood_cwgauss(struct cwgauss *dist,double x) {
  return 1/(dist->stdev*sqrt(2*PI))*exp(-1*pow(x-dist->mean,2))/(2*pow(dist->stdev,2));
}

double normalized_likelyhood_cwgauss(struct cwgauss *dist,double x) {
  double result;
  
  if (dist->stdev==0.0) {
    if (x==dist->mean)
      return 1/sqrt(2*PI);
    else
      return 0.0;
  }
  result=1/sqrt(2*PI)*exp(-1*pow((x-dist->mean)/dist->stdev,2)/2);
  
  if (isnan(result)) {
    printf("ERROR : likelyhood is nan : \n");
    printf("mean=%f\n",dist->mean);
    printf("stdev=%f\n",dist->stdev);
    printf("x=%f\n",x);
    return 0;
  }
  return result;
}

void dump_cwgauss(struct cwgauss *dist,double xmin,double xmax,int nb_steps,char *filename) {
  double x;
  FILE *f=fopen(filename,"w");
  for(x=xmin;x<=xmax;x+=(xmax-xmin)/nb_steps) {
    fprintf(f,"%f %f\n",x,normalized_likelyhood_cwgauss(dist,x));
  }
  fclose(f);
}





