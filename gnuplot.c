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

struct gplot *new_gplot(int w,int l) {
  struct gplot *result=malloc(sizeof(struct gplot));
  result->fn_index=0;
  result->pipe=popen("gnuplot","w");
  result->first_plot=1;
  
  
  if (result->pipe==NULL) {
    printf("Error : cant open gnuplot...exiting\n");
    exit(1);
  }
  //fprintf(result->pipe,"set term x11 noraise\n");
  fprintf(result->pipe,"set term wxt noraise size %d,%d\n",w,l);
  fprintf(result->pipe,"set size .721,1.\n");
  system("if test ! -d /tmp/libgem_gp/; then mkdir /tmp/libgem_gp/; fi");
  return result;
}

void plot_histo_file(struct gplot *gp,double xmin,double xmax,double box_large,char *filename,char replot,char *title) {

  int real_replot=replot;
  if (gp->first_plot && replot) {
    real_replot=0;
  }
  gp->first_plot=0;

  fprintf(gp->pipe,"set style fill solid 1.00 border lt -1\n");
  fprintf(gp->pipe,"set boxwidth %f relative\n",box_large);
  if (real_replot) {
    fprintf(gp->pipe,"replot \"%s\" with boxes title \"%s\"\n",filename,title);
  } else {
    fprintf(gp->pipe,"set xrange[%f:%f]\n",xmin,xmax);
    fprintf(gp->pipe,"plot \"%s\" with boxes title \"%s\"\n",filename,title);
  }
  
  fflush(gp->pipe);
}

void plot_file(struct gplot *gp,int dim,char *filename,char withlines,char replot,char *title) {
  
  int real_replot=replot;
  if (gp->first_plot && replot) {
    real_replot=0;
  }
  gp->first_plot=0;

  if (dim==1) {
    fprintf(gp->pipe,"set style fill solid 1.00 border lt -1\n");
    if (real_replot)
      fprintf(gp->pipe,"replot \"%s\" with boxes title \"%s\"\n",filename,title);
    else
      fprintf(gp->pipe,"plot \"%s\" with boxes title \"%s\"\n",filename,title);
  }
  
  if (dim==2) {
    if (withlines) {
      if (real_replot)
        fprintf(gp->pipe,"replot \"%s\" with lines title \"%s\"\n",filename,title);
      else
        fprintf(gp->pipe,"plot \"%s\" with lines title \"%s\"\n",filename,title);
    } else {
      if (real_replot) 
        fprintf(gp->pipe,"replot \"%s\" title \"%s\"\n",filename,title);
      else 
        fprintf(gp->pipe,"plot \"%s\" title \"%s\"\n",filename,title);
    }
  }

  if (dim==3) {
    if (withlines) {
      if (real_replot)
        fprintf(gp->pipe,"replot \"%s\" with lines title \"%s\"\n",filename,title);
      else
        fprintf(gp->pipe,"splot \"%s\" with lines title \"%s\"\n",filename,title);
    } else {
      if (real_replot) 
        fprintf(gp->pipe,"replot \"%s\" title \"%s\"\n",filename,title);
      else 
        fprintf(gp->pipe,"splot \"%s\" title \"%s\"\n",filename,title);
    }
  }

  fflush(gp->pipe);
}

void get_plotfile(struct gplot *gp,char *filename) {
  sprintf(filename,"/tmp/libgem_gp/gplot_%d.gdat",gp->fn_index);
  gp->fn_index++;
}

void replot(struct gplot *gp) {
  fprintf(gp->pipe,"replot\n");
  fflush(gp->pipe);
}

void free_gplot(struct gplot *gp) {
  if (!gp)
    return;
  if (gp->pipe!=NULL)
    fclose(gp->pipe);
  free(gp);
  system("rm -rf /tmp/libgem_gp");
}
