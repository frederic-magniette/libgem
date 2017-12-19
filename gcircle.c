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


double belonging_proba_gcircle(struct gcircle *gc,struct point *p) {
  return normalized_likelyhood_distrib(gc->pgauss,dist_circle_point(gc->circle,p));
}

double dist_gcircle(struct gcircle *gc,struct point *p) {
  return dist_circle_point(gc->circle,p);
}


//weighted cicular fit
struct circle *iRCP(struct dataset *ds,struct weights *w) {

  struct graphics *gws=NULL;
  struct circle *result=new_circle(new_point(ds->dim),0);
  struct vector *move;
  double *dists;
  struct vector *unimove;
  int p;
  double old_sd=10000;
  double sum_weights=0;
  int iter=0;
  struct distrib *dd;
  double lambda;
  int nb_iter=0;
  struct point *new_center;

  //uncomment this if you want graphical hint about circle fit
  /*
    struct point *pmin, *pmax;
    pmin=cp_point(ds->point_min);
    pmax=cp_point(ds->point_max);
    for(p=0;p<ds->dim;p++) {
    pmin->coords[p]-=2;
    pmax->coords[p]+=2;
    }
    gws=new_graphics(ds->dim,1000,1000,pmin,pmax,param2verb("detail"));
  */
 
  if (gws) {
    plot_weighted_dataset(ds,&w,1,gws);
    printf("calculating new gcircle for weights : ");
    print_weights(w);
  }
  
  for(p=0;p<ds->nb_points;p++) {
    sum_weights+=w->coefs[p];
  }
  if (sum_weights==0) {
    printf("cant evaluate a circle with all weights to zero\n");
    return result;
  }
  
  //allocate necessary memory
  dists=malloc(sizeof(double)*ds->nb_points);
  dd=new_distrib(WGAUSS);
  move=new_vector(ds->dim);
  
  //initialize the circle center to weighted barycenter 
  free_point(result->center);
  result->center=weighted_barycenter_point(ds,w);
  if (gws)
    plot_point(result->center,2,gws);
  
  //initialize the weighted distribution of distance from points to center
  reinit_distrib(dd);
  for(p=0;p<ds->nb_points;p++) {
    dists[p]=dist_points(result->center,ds->points[p]);
    add_data_distrib(dd,dists[p],1.0);
  }
  
  while(1) {
    
    nb_iter++;
    if (nb_iter>50)
      break; //TODO make a parameter
    
    //calculate the new center
    zero_vector(move,0.0);
    new_center=cp_point(result->center);
    for(p=0;p<ds->nb_points;p++) {
      unimove=vector_by2points(ds->points[p],result->center);
      normalize_vector(unimove);
      lambda=(mean_distrib(dd)-dists[p])*w->coefs[p]/sum_weights;
      if (!isnan(lambda) && !isnan_vector(unimove)) {
        mult_lambda_vect(unimove,lambda);
        add_vector(move,unimove);
      } 
      if (gws) {
	printf("coef for point %d = %f (weight=%f)\n",p,(mean_distrib(dd)-dists[p]),w->coefs[p]);
	//plot_vector(unimove,ds->points[p],iter+1,gws);
	printf("move for point %d of norme %f : ",p,norme_vector(unimove));
	print_vector(unimove);
      }
      free_vector(unimove);
    }
    
    if (gws) {
      printf("move : ");
      print_vector(move);
    }
    if (!isnan_vector(move))
      add_vector_point(new_center,move);
    else
      printf("move is nan\n");
    
    if (gws) {
      printf("new center : ");
      print_point(new_center);
      plot_vector(move,result->center,iter+1,gws);
    }
    
    free_point(result->center);
    result->center=new_center;
    
    
    //calculate the distance of points to circle center
    reinit_distrib(dd);
    for(p=0;p<ds->nb_points;p++) {
      dists[p]=dist_points(result->center,ds->points[p]);
      if (isnan(dists[p])) {
	printf("dist nan for %d between point : \n",p);
        print_point(ds->points[p]);
        printf("and circle : \n");
        print_circle(result);
        return NULL;
      }
      add_data_distrib(dd,dists[p],w->coefs[p]);
      if (gws) {
        printf("distance for point %d = %f\n",p,dists[p]);
      }
    }
    if (gws) {
      printf("mean dist=%f\n",mean_distrib(dd));
      printf("standard deviation = %f\n",stdev_distrib(dd));
    }
    
    //calculate the radius
    result->radius=mean_distrib(dd);
    
    if (gws)
      plot_circle(result,iter,gws); 
    
    //check for ending
    if (fabs(stdev_distrib(dd)-old_sd)<0.01)
      break;
    if (iter>1)
      old_sd=stdev_distrib(dd);
    
    //printf("continuing\n");
    iter++;
    if (gws) {
      apply_splot(gws->sp);
      getchar();
    }
  } //end of while(1)
  

  if (gws) {
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("!!!!!!!!!!!!!ENDED!!!!!!!!!!!!!!!\n");
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    getchar();
  }

  //free allocated memory
  free(dists);
  free_distrib(dd);
  free_vector(move);
  return result;
}

struct gcircle *new_gcircle(struct dataset *ds,struct weights *w,struct graphics *gws) {
  int p;
  struct gcircle *result=malloc(sizeof(struct gcircle));
  result->pgauss=new_distrib(CWGAUSS);
  result->circle=iRCP(ds,w);
  if (result->circle==NULL) {
    free_distrib(result->pgauss);
    free(result);
    return NULL;
  }

  //calculate the centred weighted gaussian distrib for new circle
  for(p=0;p<ds->nb_points;p++) {
    add_data_distrib(result->pgauss,dist_circle_point(result->circle,ds->points[p]),w->coefs[p]);
  }

  return result;
}

struct gcircle *cp_gcircle(struct gcircle *gc) {
  struct gcircle *result=malloc(sizeof(struct gcircle));
  result->circle=cp_circle(gc->circle);
  result->pgauss=cp_distrib(gc->pgauss);
  return result;
}

void print_gcircle(struct gcircle *gc) {
  print_circle(gc->circle);
  print_distrib(gc->pgauss);
}

void dump_gcircle(struct gcircle *gc,FILE *f) {
  dump_circle(gc->circle,f);
}

void free_gcircle(struct gcircle *gc) {
  free_circle(gc->circle);
  free_distrib(gc->pgauss);
  free(gc);
}

void plot_gcircle(struct gcircle *gc,int id,struct graphics *gws) {
  plot_circle(gc->circle,id,gws);
}

void plot_field_gcircle(struct gcircle *gc,int id,struct graphics *gws) {
  int rref,gref,bref;
  double distance;
  int x,y;
  int r,g,b;
  double proba;
  double llh;
  struct point *pt;
  if (gws->sp) {
    //plot the field
    for(x=0;x<gws->w;x++) {
      for(y=0;y<gws->h;y++) {

        //calc the proba
	pt=sdl2real_splot(x,y,gws->sp);
	distance=dist_gcircle(gc,pt);
	free_point(pt);
        llh=normalized_likelyhood_distrib(gc->pgauss,distance);
        if (isnan(llh)) {
          printf("error : likelihood is nan\n");
          printf("dont print point %d,%d\n",x,y);
          continue;
        }
	proba=log(1+llh*sqrt(2*PI));
        
	//calc the color : proportional and thresholded
        get_pixel_splot(gws->sp,x,y,&r,&g,&b);
        attribute_color_splot(id,&rref,&gref,&bref);
        r+=(int)(proba*rref);
        g+=(int)(proba*gref);
        b+=(int)(proba*bref);

        set_pixel_splot(gws->sp,x,y,r,g,b);
      }
    }
  }
}

int same_gcircle(struct gcircle *gc1,struct gcircle *gc2,double distlim) {
  if ((dist_points(gc1->circle->center,gc2->circle->center)>distlim) || (fabs(gc1->circle->radius-gc2->circle->radius)>distlim))
    return 0;
  else
    return 1;
}

