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

//*************************** DATASET *******************************

struct dataset * new_dataset_zero(int nb_points,int dim) {
  struct dataset *ds=malloc(sizeof(struct dataset));
  int i,j;
  ds->nb_points=nb_points;
  ds->dim=dim;
  ds->points=malloc(nb_points*sizeof(struct point *));
  for (i=0;i<nb_points;i++) {
    ds->points[i]=malloc(sizeof(struct point));
    ds->points[i]->dim=dim;
    ds->points[i]->coords=malloc(dim*sizeof(double));
    for(j=0;j<dim;j++) {
      ds->points[i]->coords[j]=0;
    }
  }
  ds->point_min=new_point(dim);
  ds->point_max=new_point(dim);
  return ds;
}

struct dataset * new_empty_dataset(int dim) {
  struct dataset *ds=malloc(sizeof(struct dataset));
  ds->nb_points=0;
  ds->points=NULL;
  ds->dim=dim;
  ds->point_min=new_point(dim);
  ds->point_max=new_point(dim);
  return ds;
}

struct dataset * new_empty_boxed_dataset(int dim,struct point *point_min,struct point *point_max) {
  struct dataset *ds=malloc(sizeof(struct dataset));
  ds->nb_points=0;
  ds->points=NULL;
  ds->dim=dim;
  ds->point_min=cp_point(point_min);
  ds->point_max=cp_point(point_max);
  return ds;
}

struct dataset * new_dataset_fromfile(int dim,char *filename) {
  FILE *f;
  struct dataset *ds;
  int i,j;
  int nb_points=0;
  int res=1;
  double tmpf;

  f=fopen(filename,"r");
  if (f==NULL)
    return NULL;
  
  while(res==1) {
    for(j=0;j<dim;j++) {
      res=fscanf(f,"%lf",&tmpf);
      if (res!=1)
        break;
      //printf("%lf ",tmpf);
    }
    //printf("\n");
    if (res!=1)
      break;
    nb_points++;
    //printf("incrementing nb_points result=%d\n",nb_points);
  }
  //printf("nb_points=%d\n",nb_points);
  rewind(f);

  ds=new_dataset_zero(nb_points,dim);
  for(i=0;i<nb_points;i++) {    
    for(j=0;j<dim;j++) {
      fscanf(f,"%lf",&(ds->points[i]->coords[j]));
    }
  }
  free_point(ds->point_min);
  free_point(ds->point_max);
  ds->point_min=new_point(dim);
  ds->point_max=new_point(dim);
  calc_min_point_dataset(ds);
  calc_max_point_dataset(ds);
  return ds;
}

void print_dataset(struct dataset *ds) {
  int i;
  printf("dataset(dim=%d)[%d]=\n",ds->dim,ds->nb_points);
  for(i=0;i<ds->nb_points;i++) {
    printf("point %d ",i);
    print_point(ds->points[i]);
  }
    //print max and min coords
  printf("min and max coords : \n");
  print_point(ds->point_min);
  print_point(ds->point_max);
}

void dump_dataset(struct dataset *ds,char *filename) {
  int d,p;
  FILE *f=fopen(filename,"w");
  for(p=0;p<ds->nb_points;p++) {
    for(d=0;d<ds->dim;d++) {
      fprintf(f,"%f ",ds->points[p]->coords[d]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
}

void calc_min_point_dataset(struct dataset *ds) {
  int d,p;
  if (ds->nb_points==0) {
    for(d=0;d<ds->dim;d++) {
      ds->point_min->coords[d]=0;
    }
    return;
  }
   
  for(d=0;d<ds->dim;d++) {
    ds->point_min->coords[d]=ds->points[0]->coords[d];
  }
  for(p=0;p<ds->nb_points;p++) {
    for(d=0;d<ds->dim;d++) {
      if (ds->points[p]->coords[d]<ds->point_min->coords[d])
        ds->point_min->coords[d]=ds->points[p]->coords[d];
    }
  }
}

void calc_max_point_dataset(struct dataset *ds) {
  int d,p;
  if (ds->nb_points==0) {
    for(d=0;d<ds->dim;d++) {
      ds->point_max->coords[d]=0;
    }
    return;
  }
   
  for(d=0;d<ds->dim;d++) {
    ds->point_max->coords[d]=ds->points[0]->coords[d];
  }
  for(p=0;p<ds->nb_points;p++) {
    for(d=0;d<ds->dim;d++) {
      if (ds->points[p]->coords[d]>ds->point_max->coords[d])
        ds->point_max->coords[d]=ds->points[p]->coords[d];
    }
  }
}

void free_dataset(struct dataset *ds) {
  int i;
  for(i=0;i<ds->nb_points;i++) {
    free_point(ds->points[i]);
  }
  free(ds->points);
  free_point(ds->point_min);
  free_point(ds->point_max);
  free(ds);
}

#define HISTO_SIZE 20

void plot_weighted_dataset(struct dataset *ds,struct weights **w,int nb_lines,struct graphics *gws) {
  char name[1024];
  int l,p;
  int rref,gref,bref;
  int r,g,b;
  int histo[HISTO_SIZE];
  double histovals[HISTO_SIZE];
  int i;
  FILE *f;
  double xmin,xmax;
  double step;

  if (ds->dim==1) {
    if (gws->gp) {
      get_plotfile(gws->gp,name);
      xmin=ds->point_min->coords[0];
      xmax=ds->point_max->coords[0];
      step=(xmax-xmin)/HISTO_SIZE;
      for(i=0;i<20;i++) {
        histovals[i]=xmin+i*step;
        histo[i]=0;
      }
      for(p=0;p<ds->nb_points;p++) {
        for(i=0;i<HISTO_SIZE-1;i++) {
          if((ds->points[p]->coords[0]>=histovals[i]) && (ds->points[p]->coords[0]<histovals[i+1]))
            histo[i]++;
        }
      }
      f=fopen(name,"w");
      for(i=0;i<HISTO_SIZE;i++) 
        fprintf(f,"%f %d\n",histovals[i]+(step/2),histo[i]);
      fclose(f);
      plot_histo_file(gws->gp,xmin,xmax,step,name,0,"data");
    }
  }

  if (ds->dim==2) {
    if (gws->sp) {
      for(p=0;p<ds->nb_points;p++) {
        r=255;
        g=255;
        b=255;
        for(l=0;l<nb_lines;l++) {
          attribute_color_splot(l,&rref,&gref,&bref);
          //r+=(int)(w[l]->coefs[p]*rref);
          //g+=(int)(w[l]->coefs[p]*gref);
          //b+=(int)(w[l]->coefs[p]*bref);
        }
        plot_real_point_splot(ds->points[p],gws->sp,r,g,b);
      }
    }
  }

  if (ds->dim==3) {
    if (gws->gp) {
      get_plotfile(gws->gp,name);
      dump_dataset(ds,name);
      plot_file(gws->gp,ds->dim,name,0,0,"data");
    }
  }
}

void plot_dataset(struct dataset *ds,struct graphics *gws) {
  char name[1024];
  int p;
  int histo[HISTO_SIZE];
  double histovals[HISTO_SIZE];
  int i;
  FILE *f;
  double xmin,xmax;
  double step;

  if (ds->dim==1) {
    if (gws->gp) {
      //TODO convert to real value histogram in a new class
      get_plotfile(gws->gp,name);
      xmin=ds->point_min->coords[0];
      xmax=ds->point_max->coords[0];
      step=(xmax-xmin)/HISTO_SIZE;
      for(i=0;i<20;i++) {
        histovals[i]=xmin+i*step;
        histo[i]=0;
      }
      for(p=0;p<ds->nb_points;p++) {
        for(i=0;i<HISTO_SIZE-1;i++) {
          if((ds->points[p]->coords[0]>=histovals[i]) && (ds->points[p]->coords[0]<histovals[i+1]))
            histo[i]++;
        }
      }
      f=fopen(name,"w");
      for(i=0;i<HISTO_SIZE;i++) 
        fprintf(f,"%f %d\n",histovals[i]+(step/2),histo[i]);
      fclose(f);
      plot_histo_file(gws->gp,xmin,xmax,step,name,0,"data");
    }
  }

  if (ds->dim==2) {
    if (gws->sp) {
      //clear_splot(gws->sp);
      for(p=0;p<ds->nb_points;p++) {
        plot_real_point_splot(ds->points[p],gws->sp,255,255,0);
      }
      apply_splot(gws->sp);
    }
  }

  if (ds->dim==3) {
    if (gws->gp) {
      get_plotfile(gws->gp,name);
      dump_dataset(ds,name);
      plot_file(gws->gp,ds->dim,name,0,0,"data");
    }
  }
}

void add_line_dataset(struct dataset *ds,struct line *l,int nb_steps,double length) {
  int i;
  double tmp;
  struct vector *norm_dir=cp_vector(l->dir_vect);
  struct point **points;
  struct point *point_min=cp_point(l->ref);
  struct point *point_max=cp_point(l->ref);

  //calc the limits
  mult_lambda_vect(norm_dir,(length/2)/norme_vector(norm_dir));
  add_vector_point(point_max,norm_dir);
  mult_lambda_vect(norm_dir,-1);
  add_vector_point(point_min,norm_dir);
  for(i=0;i<2;i++) {
    if (point_min->coords[i]>point_max->coords[i]) {
      tmp=point_max->coords[i];
      point_max->coords[i]=point_min->coords[i];
      point_min->coords[i]=tmp;
    }
      
  }

  points=boxed_line(l,point_min,point_max,nb_steps);

  for(i=0;i<nb_steps;i++) {
    if (points[i]!=NULL) {
      add_point_dataset(ds,points[i]);
    }
  }
  free(points);
  free_point(point_min);
  free_point(point_max);
  free_vector(norm_dir);
}

void add_angle_line_dataset(struct dataset *ds,int angle,int nb_steps,double x,double y,double length) {
  struct line *angline=angle_line(angle,x,y);
  add_line_dataset(ds,angline,nb_steps,length);
  free_line(angline);
}

void add_random_line_dataset(struct dataset *ds,int nb_steps,double length) {
  struct line *rline=random_line(ds->dim,ds->point_min,ds->point_max);
  add_line_dataset(ds,rline,nb_steps,length);
  free_line(rline);
}

void add_point_dataset(struct dataset *ds,struct point *p) {
  int newsize=ds->nb_points+1;
  ds->points=realloc(ds->points,newsize*sizeof(struct point *));
  ds->points[ds->nb_points]=p;
  ds->nb_points=newsize;
  calc_max_point_dataset(ds);
  calc_min_point_dataset(ds);
}

void add_random_noise_dataset(struct dataset *ds,int nb_points) {
  int i,d;
  struct point *p;
  double unif;
  for(i=0;i<nb_points;i++) {
    p=new_point(ds->dim);
    for(d=0;d<ds->dim;d++) {
      unif=rand()/((double)RAND_MAX);
      p->coords[d]=ds->point_min->coords[d]+unif*(ds->point_max->coords[d]-ds->point_min->coords[d]);
    }
    add_point_dataset(ds,p);
  }
}

void noise_dataset(struct dataset *ds,double dispersion) {
  int d,p;
  struct normal_gene *gene=new_normal_gene();
  for(d=0;d<ds->dim;d++) {
    for(p=0;p<ds->nb_points;p++) {
      ds->points[p]->coords[d]+=gener_normal(gene)*dispersion;
    }
  }
  free_normal_gene(gene);
  calc_max_point_dataset(ds);
  calc_min_point_dataset(ds);
}

struct dataset *ortho_projection_dataset(struct dataset *ds,struct matrix *base,struct vector *dir,struct point *ref,double *long_dist) {
  struct dataset *result=new_empty_dataset(ds->dim-1);
  struct point *p_bp;
  int i;
  p_bp=change_base_point(ref,base);
  printf("projection of ref\n");
  print_point(p_bp);
  free_point(p_bp);
  //printf("new base : \n");
  //print_matrix(base);
  for(i=0;i<ds->nb_points;i++) {
    //project point
    p_bp=ortho_project_point(ds->points[i],base,dir,ref);
    //save the longitudinal distance
    if (long_dist!=NULL)
      long_dist[i]=p_bp->coords[0];
    //add projected point to dataset
    add_point_dataset(result,cut_first_coord_point(p_bp));
    free_point(p_bp);
    //free_point(p);
  }
  return result;
}
