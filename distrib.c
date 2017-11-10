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

struct distrib *new_distrib(int type) {
  struct distrib *result=malloc(sizeof(struct distrib));
  memset(result,0,sizeof(struct distrib));
  result->type=type;
  switch(type) {
  case WGAUSS:
    //printf("new wgauss\n");
    result->data=(void *)new_wgauss();
    return result;
    break;
  case CWGAUSS:
    //printf("new cwgauss\n");
    result->data=(void *)new_cwgauss();
    return result;
    break;  
  default:
    printf("error unknown distribution type %d\n",type);
    return NULL;
  }
}

void add_data_distrib(struct distrib *d,double value,double weight) {
  switch(d->type) {
  case WGAUSS:
    add_data_wgauss((struct wgauss *)d->data,value,weight);
    break;
  case CWGAUSS:
    add_data_cwgauss((struct cwgauss *)d->data,value,weight);
    break;  
  default:
    break;
  }
}

void reinit_distrib(struct distrib *d) {
  switch(d->type) {
  case WGAUSS:
    reinit_wgauss((struct wgauss *)d->data);
    break;
  case CWGAUSS:
    reinit_cwgauss((struct cwgauss *)d->data);
    break;  
  default:
    break;
  }
}

void print_distrib(struct distrib *d) {
  switch(d->type) {
  case WGAUSS:
    printf("weighted gaussian distribution : \n");
    print_wgauss((struct wgauss *)d->data);
    break;
  case CWGAUSS:
    printf("centred weighted gaussian distribution : \n");
    print_cwgauss((struct cwgauss *)d->data);
    break;
  default:
    printf("unknown distribution\n");
    break;
  }
}

struct distrib *cp_distrib(struct distrib *d) {
  struct distrib *result=malloc(sizeof(struct distrib));
  memset(result,0,sizeof(struct distrib));
  result->type=d->type;
  switch(d->type) {
  case WGAUSS:
    result->data=(void *)cp_wgauss((struct wgauss *)d->data);
    break;
  case CWGAUSS:
    result->data=(void *)cp_cwgauss((struct cwgauss *)d->data);
    break;  
  default:
    break;
  }
  return result;
}


void free_distrib(struct distrib *d) {
  switch(d->type) {
  case WGAUSS:
    free_wgauss((struct wgauss *)d->data);
    break;
  case CWGAUSS:
    free_cwgauss((struct cwgauss *)d->data);
    break;
  default:
    break;
  }
  free(d);
}

double normalized_likelyhood_distrib(struct distrib *d,double x) {
  switch(d->type) {
  case WGAUSS:
    return normalized_likelyhood_wgauss((struct wgauss *)d->data,x);
    break;
  case CWGAUSS:
    return normalized_likelyhood_cwgauss((struct cwgauss *)d->data,x);
    break; 
  default:
    return 0.0;
    break;
  }
}

double mean_distrib(struct distrib *d) {
switch(d->type) {
 case WGAUSS:
   return ((struct wgauss *)d->data)->mean;
   break;
 case CWGAUSS:
   return ((struct cwgauss *)d->data)->mean;
   break;
 default:
   return 0.0;
   break;
 }
}


double stdev_distrib(struct distrib *d) {
switch(d->type) {
 case WGAUSS:
   return ((struct wgauss *)d->data)->stdev;
   break;
 case CWGAUSS:
   return ((struct cwgauss *)d->data)->stdev;
   break;
 default:
   return 0.0;
   break;
 }
}

void set_stdev_distrib(struct distrib *d,double stdev)  {
switch(d->type) {
 case WGAUSS:
   ((struct wgauss *)d->data)->stdev=stdev;
   break;
 case CWGAUSS:
   ((struct cwgauss *)d->data)->stdev=stdev;
   break;
 default:
   break;
 }
}

void dump_distrib(struct distrib *d,double center,double xmin,double xmax,int nb_steps,char *filename) {
  double x;
  double value;
  FILE *f=fopen(filename,"w");
  for(x=xmin;x<=xmax;x+=(xmax-xmin)/nb_steps) {
    value=normalized_likelyhood_distrib(d,fabs(x-center));
    fprintf(f,"%f %f\n",x,value);
  }
  fclose(f);
}

