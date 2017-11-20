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

//************************* EM ALGO ****************************

struct gem_ws *init_gem(struct dataset *ds,double convcrit) { 
  struct point *pmin=cp_point(ds->point_min);
  struct point *pmax=cp_point(ds->point_max);
  struct gem_ws *result=malloc(sizeof(struct gem_ws));
  result->ds=ds;
  result->nb_objects=0;
  result->nb_iter=0;
  result->global_dist=0;
  result->scale=dist_points(pmin,pmax);
  result->convergence_criterion=convcrit;
  result->fit_weights=NULL;
  result->distance=NULL;
  result->fit_objects=NULL;
  free_point(pmin);
  free_point(pmax);
  return result;
}

int add_object_gem(struct gem_ws *ws,enum otype type,struct graphics *gws) {

  //this function creates a new object to get a better fit
  
  int p,l,o;
  double max_dist=0.0;
  int wf=0;
  double *dist_to_ref=malloc(sizeof(double)*ws->ds->nb_points);
  int newo_idx=ws->nb_objects;
  struct neighbouring *wfn;
  double pc_out=0.000001; //weight for excluded objects : to avoid 0 
  int obj_type;
  int dsindex;
  ws->nb_objects++;

  //reallocate weights
  ws->fit_weights=realloc(ws->fit_weights,sizeof(struct weights *)*ws->nb_objects);
  ws->fit_weights[newo_idx]=new_weights(ws->ds->nb_points,0.0);


  //reallocate distances
  ws->distance=realloc(ws->distance,sizeof(double *)*ws->nb_objects);
  ws->distance[newo_idx]=malloc(sizeof(double)*ws->ds->nb_points);
  for(p=0;p<ws->ds->nb_points;p++) {
    ws->distance[newo_idx][p]=0.0;
  }

  //if there is only one object : use uniform weights
  if (ws->nb_objects==1) {
    for(p=0;p<ws->ds->nb_points;p++) {
      ws->fit_weights[0]->coefs[p]=1.0;
      if (gws) {
        plot_point(ws->ds->points[p],newo_idx,gws);
      }
    }
  } else {

    //we calculate the weights for the new object : ws->fit_weights[newo_idx]
    //procedure is to find the worst fitted point (wfp) in the dataset
    //i.e. the point which have the longuest distance to its reference object
    
    //then find the neighbour of wfp that have a distance to wfp less than the fifth of the data scale

    //for the points that cope this we change weights : 1-pc_out in new object and pc_out for other objects
    //for the others points : pc_out for this object and unchanged for other objects
    
    //get the worst fitted point (wfp)
    for(p=0;p<ws->ds->nb_points;p++) {
      dist_to_ref[p]=ws->distance[0][p];
      for(l=1;l<newo_idx;l++) {
        if (ws->distance[l][p]<dist_to_ref[p]) {
          dist_to_ref[p]=ws->distance[l][p];
        }
      }
      if (dist_to_ref[p]>max_dist) {
        max_dist=dist_to_ref[p];
	wf=p;
      }
    }
    //printf("worst fitted point is %d : ",wf);
    //print_point(ws->ds->points[wf]);

    //get the neighbouring of wfp
    wfn=new_neighbouring(ws->ds,ws->ds->points[wf]);

    //wfp new weights
    ws->fit_weights[newo_idx]->coefs[wf]=1-pc_out;
    //affect the weight for all other objects
    for(l=0;l<newo_idx;l++) {
      ws->fit_weights[l]->coefs[wf]=pc_out/(float)newo_idx;
      printf("removing point %d from object %d\n",wf,l);
    }
    if (gws) {
      plot_point(ws->ds->points[wf],newo_idx+1,gws);
    }

    printf("choosing with respect to scale/2=%f\n",ws->scale/2);
    //printf("choosing with respect to distance to reference\n");
    //calculate the new weights
    for(p=0;p<wfn->size;p++) {
      dsindex=wfn->neighbours[p]->ds_index;
      //choose the point not too far from wfp
      //TODO : should be a parameter
      //if (wfn->neighbours[p]->dist<ws->scale/2) {
      if (wfn->neighbours[p]->dist<dist_to_ref[wf]) {

	//affect the weight for the new object
        ws->fit_weights[newo_idx]->coefs[dsindex]=1-pc_out;
        if (gws) {
          plot_point(wfn->neighbours[p]->p,newo_idx,gws);
        }
	
	//affect the weight for all other objects
        for(l=0;l<newo_idx;l++) {
          ws->fit_weights[l]->coefs[dsindex]=pc_out/(float)newo_idx;
          //printf("removing point %d from object %d\n",dsindex,l);
        }
      } else {
	//for non choosen points
        ws->fit_weights[newo_idx]->coefs[dsindex]=pc_out;
        for(l=0;l<newo_idx;l++) {
          ws->fit_weights[l]->coefs[dsindex]-=pc_out/(float)newo_idx;
        }
      }
    }
    free_neighbouring(wfn);
  }

  if (gws && gws->verbosity==VERB_DEBUG) {
    for(p=0;p<ws->ds->nb_points;p++) {
      printf("point %d ",p);
      for(o=0;o<ws->nb_objects;o++) {
        printf("%f ",ws->fit_weights[o]->coefs[p]);
      }
      printf("\n");
    }
    printf("points choosen\n");
    apply_graphics(gws);
    getchar();
    }
  

  //create the new object
  ws->fit_objects=realloc(ws->fit_objects,sizeof(struct object *)*ws->nb_objects);
  ws->fit_objects[newo_idx]=new_object(type,ws->ds,ws->fit_weights[newo_idx],gws);
  //print_object(ws->fit_objects[newo_idx]);

  //modify the stdev of the new object : TODO use the parameter
  if (stdev_object(ws->fit_objects[newo_idx])<(ws->scale/20))
    set_stdev_object(ws->fit_objects[newo_idx],ws->scale/20);

  //recreate the other objects with their new weights
  for(o=0;o<newo_idx;o++) {
    obj_type=ws->fit_objects[o]->type;
    free_object(ws->fit_objects[o]);
    ws->fit_objects[o]=new_object(obj_type,ws->ds,ws->fit_weights[o],gws);
  }
  
  //reset the nb of iterations
  ws->nb_iter=0;
  
  free(dist_to_ref);
  return newo_idx;
}

int add_random_object_gem(struct gem_ws *ws,enum otype type) {
  int p,l;
  int newo_idx=ws->nb_objects;
  double pc_out=0.001; //weight for excluded objects : to avoid 0 
  ws->nb_objects++;

  //reallocate weights
  ws->fit_weights=realloc(ws->fit_weights,sizeof(struct weights *)*ws->nb_objects);
  ws->fit_weights[newo_idx]=new_weights(ws->ds->nb_points,0.0);


  //reallocate distances
  ws->distance=realloc(ws->distance,sizeof(double *)*ws->nb_objects);
  ws->distance[newo_idx]=malloc(sizeof(double)*ws->ds->nb_points);
  for(p=0;p<ws->ds->nb_points;p++) {
    ws->distance[newo_idx][p]=0.0;
  }

  //if there is only one object : use uniform weights
  if (ws->nb_objects==1) {
    for(p=0;p<ws->ds->nb_points;p++) {
      ws->fit_weights[0]->coefs[p]=1.0;
    }
  } else {   
    for(p=0;p<ws->ds->nb_points;p++) {
      if (rand()/((double)RAND_MAX)>0.5) {
	//affect the weight for the new object
        ws->fit_weights[newo_idx]->coefs[p]=1-pc_out;
	//affect the weight for all other objects
	for(l=0;l<newo_idx;l++) {
          ws->fit_weights[l]->coefs[p]=pc_out;
        }
      } else {
	//for non choosen points
        ws->fit_weights[newo_idx]->coefs[p]=pc_out;
      }
    }
  }
  

  //create the new object
  ws->fit_objects=realloc(ws->fit_objects,sizeof(struct object *)*ws->nb_objects);
  ws->fit_objects[newo_idx]=new_object(type,ws->ds,ws->fit_weights[newo_idx],NULL);
  //print_object(ws->fit_objects[newo_idx]);

  //modify the stdev of the new object : TODO use the parameter
  if (stdev_object(ws->fit_objects[newo_idx])<(ws->scale/20))
    set_stdev_object(ws->fit_objects[newo_idx],ws->scale/20);

  //recreate the other objects with their new weights
  int o;
  for(o=0;o<newo_idx;o++) {
    free_object(ws->fit_objects[o]);
    ws->fit_objects[o]=new_object(type,ws->ds,ws->fit_weights[o],NULL);
  }
  
  //reset the nb of iterations
  ws->nb_iter=0;
  
  return newo_idx;
}

void remove_object_gem(struct gem_ws *ws,int obj_idx) {
  int l;
  for(l=obj_idx;l<ws->nb_objects-1;l++)
    ws->fit_objects[l]=ws->fit_objects[l+1];
  for(l=obj_idx;l<ws->nb_objects-1;l++)
    ws->fit_weights[l]=ws->fit_weights[l+1];
  ws->nb_objects--;
  //printf("%d objects left\n",ws->nb_objects);
}


void print_gem(struct gem_ws *ws) {
  int l;
  
  printf("print em workspace : \n");

  //print the fit objects
  for(l=0;l<ws->nb_objects;l++) {
    printf("object %d : ",l);
    print_object(ws->fit_objects[l]);
  }
  printf("iterations number %d\n",ws->nb_iter);
  printf("global distance=%f\n",ws->global_dist);
  
}

void free_gem(struct gem_ws *ws) {
  int l;

  if (!ws)
    return;

  //free fit weights
  if (ws->fit_weights) {
    for(l=0;l<ws->nb_objects;l++) {
      free_weights(ws->fit_weights[l]);
    }
    free(ws->fit_weights);
  }

  //free fit objects
  if (ws->fit_objects) {
    for(l=0;l<ws->nb_objects;l++) {
      free_object(ws->fit_objects[l]);
    }
    free(ws->fit_objects);
    }

  //free distance
  if (ws->distance) {
    for(l=0;l<ws->nb_objects;l++) {
      free(ws->distance[l]);
    }
    free(ws->distance);
  }

  //free(ws->non_fitted);
  free(ws);
}

struct gem_ws *cp_gem(struct gem_ws *src) {
  int l,p;
  struct gem_ws *result=malloc(sizeof(struct gem_ws));
  result->ds=src->ds;
  result->nb_objects=src->nb_objects;
  result->nb_iter=src->nb_iter;
  result->scale=src->scale;
  result->convergence_criterion=src->convergence_criterion;
  result->global_dist=src->global_dist;

  result->fit_weights=malloc(sizeof(struct weights *)*src->nb_objects);
  for(l=0;l<src->nb_objects;l++) {
    result->fit_weights[l]=cp_weights(src->fit_weights[l]);
  }

  //allocate distances
  result->distance=malloc(sizeof(double *)*src->nb_objects);
  for(l=0;l<src->nb_objects;l++) {
    result->distance[l]=malloc(sizeof(double)*src->ds->nb_points);
    for(p=0;p<src->ds->nb_points;p++) {
      result->distance[l][p]=src->distance[l][p];
    }
  }
  
  //allocate objects
  result->fit_objects=malloc(sizeof(struct object *)*src->nb_objects);
  for(l=0;l<src->nb_objects;l++) {
    result->fit_objects[l]=cp_object(src->fit_objects[l]);
  }
    
  return result;
}

void plot_gem(struct gem_ws *ws,struct graphics *gws) {
  int l;

  //print_gem(ws);
  
  if (!gws)
    return;
  
  if (gws->sp) {
    clear_splot(gws->sp);
  }
  if (gws->gp) {
    plot_weighted_dataset(ws->ds,ws->fit_weights,ws->nb_objects,gws); 
  }

  //plot fit proba field
  for(l=0;l<ws->nb_objects;l++) {
    plot_field_object(ws->fit_objects[l],l,gws);
  }
  
  //plot fit objects
  for(l=0;l<ws->nb_objects;l++) {
    plot_object(ws->fit_objects[l],l,gws);
  }

  if (gws->sp) {
    //plot points
    plot_weighted_dataset(ws->ds,ws->fit_weights,ws->nb_objects,gws);
    apply_splot(gws->sp);
  }
  
}

void update_distances_gem(struct gem_ws *ws,struct graphics *gws) {
  int p,o;
  //calculate the distance between any point and any line
  for(p=0;p<ws->ds->nb_points;p++) {
    for(o=0;o<ws->nb_objects;o++) {
      ws->distance[o][p]=dist_object(ws->fit_objects[o],ws->ds->points[p]);
    }
  }
  if (gws) {
    for(p=0;p<ws->ds->nb_points;p++) {
      printf("distance for point %d : ",p);
      for(o=0;o<ws->nb_objects;o++) {
        printf("%f ", ws->distance[o][p]);
      }
      printf("\n");
    }
  }
}

int algo_gem(struct gem_ws *ws,struct graphics *gws) {
  int p; //point iterator
  int o; //object iterator
  enum otype type;

  int hvar=1;
  struct distrib *global_distrib=new_distrib(CWGAUSS);

  //alias
  double nb_points=ws->ds->nb_points;
  double nb_objects=ws->nb_objects;

  //probability of points
  double **bprobas; //probability of belonging for all points to every objects
  double *total_proba; // total of probability from every gaussian for every points
  double *pk;
  
  //convergence criterion
  int converged=0;
  double diffcrit;
  double oldcrit=1.0;
  double newcrit;
  //double loglik; //log-likelyhood

  int max_iter=100;

  //non fitted points
  //int nf;
  
  //FILE *floglik=fopen("/tmp/loglik.gdat","w");

  if (nb_objects==0) {
    printf("error : no objects for em\n");
    return 0;
  }

  //allocate probas and probas sums
  bprobas=new_dtab(nb_points,nb_objects);
  total_proba=malloc(sizeof(double)*nb_points);
  pk=malloc(sizeof(double)*nb_objects);
    
  do {
    //printf("iteration %d\n",ws->nb_iter);
    ws->nb_iter++;

    update_distances_gem(ws,gws);

    //global distribution for homogenous variance option
    if (hvar) {
      //printf("homogenous variance option active\n");
      reinit_distrib(global_distrib);
      for(p=0;p<nb_points;p++) {
	for(o=0;o<nb_objects;o++) {
	  add_data_distrib(global_distrib,ws->distance[o][p],ws->fit_weights[o]->coefs[p]);
	}
      }
    }
      
    //calculate the probas and their sum for every point
    for(p=0;p<nb_points;p++) {
      total_proba[p]=0;
      for(o=0;o<nb_objects;o++) {
	if (hvar) {
	  bprobas[p][o]=normalized_likelyhood_distrib(global_distrib,ws->distance[o][p]);
	  total_proba[p]+=bprobas[p][o];
	} else {
	  bprobas[p][o]=belonging_proba_object(ws->fit_objects[o],ws->ds->points[p]);
	  total_proba[p]+=bprobas[p][o];
	}
      }
    }

    //calculate the pk
    for(o=0;o<nb_objects;o++) {
      pk[o]=0;
      for(p=0;p<nb_points;p++) {
        pk[o]+=ws->fit_weights[o]->coefs[p];
      }
      pk[o]/=nb_points;
    }

    //print the probas if necessary
    if (gws) {
      for(p=0;p<nb_points;p++) {
	printf("probas for point %d : ",p);
	for(o=0;o<nb_objects;o++) {
	  printf("%f ",bprobas[p][o]);
	}
	printf(" (total = %f)\n",total_proba[p]);
      }
    }

    //reevaluate the weights with pk
    //printf("reevaluate the weights\n");
    
      int j;
      double spk;
      for(p=0;p<nb_points;p++) {
      spk=0;
      for(j=0;j<nb_objects;j++) {
	spk+=pk[j]*bprobas[p][j];
      }
      for(o=0;o<nb_objects;o++) {
	if (spk==0) {
	  printf("equirepartition\n");
          ws->fit_weights[o]->coefs[p]=1.0/nb_objects;
        } else {
	  ws->fit_weights[o]->coefs[p]=pk[o]*bprobas[p][o]/spk;
	}
      }
      }

    //reevaluate the weights without pk
    //printf("reevaluate the weights\n");
    /*for(o=0;o<nb_objects;o++) {
      for(p=0;p<nb_points;p++) {
        if (total_proba[p]==0)
          ws->fit_weights[o]->coefs[p]=1.0/nb_objects;
        else
	  ws->fit_weights[o]->coefs[p]=bprobas[p][o]/total_proba[p];
      }
      }*/

    //print the weights if necessary
    if (gws) {
      for(p=0;p<nb_points;p++) {
	printf("weight for point %d : ",p);
	for(o=0;o<nb_objects;o++) {
	  printf("%f ",ws->fit_weights[o]->coefs[p]);
	}
	printf("\n");
      }
    }

    //recompute the objects
    for(o=0;o<nb_objects;o++) {
      type=ws->fit_objects[o]->type;
      free_object(ws->fit_objects[o]);
      ws->fit_objects[o]=new_object(type,ws->ds,ws->fit_weights[o],gws);
    }

    //calculate the global distance measure
    ws->global_dist=0;
    int best;
    double best_weight=0;
    for(p=0;p<nb_points;p++) {
      best_weight=0;
      best=-1;
      for(o=0;o<nb_objects;o++) {
        if (ws->fit_weights[o]->coefs[p]>best_weight) {
          best_weight=ws->fit_weights[o]->coefs[p];
          best=o;
        }
      }
      //printf("for point %d best is %d with weight %f\n",p,best,best_weight);
      ws->global_dist+=dist_object(ws->fit_objects[best],ws->ds->points[p]);
    }
    //printf("global distance=%f\n",ws->global_dist);

    //choice of the convergence criterion
    newcrit=ws->global_dist;

    
    //calculate the convergence criterion
    diffcrit=fabs(oldcrit-newcrit);
    //printf("oldcrit=%f newcrit=%f diffcrit=%f\n",oldcrit,newcrit,diffcrit);
    if (isnan(diffcrit)) {
      printf("criteria is nan->converged\n");
      converged=1;
    }
    if (!converged) { 
      //printf("convergence criterion=%f (vs %f)\n",diffcrit,ws->convergence_criterion); 
      if (diffcrit<=ws->convergence_criterion) {
        converged=1;
      } else {
        converged=0;
      }
    }

    
    oldcrit=newcrit;
    
    if (converged==0) {
      if (gws) {
        printf("not converged\n");
        print_gem(ws);
        plot_gem(ws,gws);
	printf("Press a key to continue\n");
	getchar();
      }
    } else {
      if (gws) {
	printf("CONVERGED!!!!!!!!!!!!!\n");
        //print_gem(ws);
	plot_gem(ws,gws);
        printf("Press a key to continue\n");
        getchar();
      }
    }

    if (ws->nb_iter>max_iter) {
      if (gws)
	printf("em has diverged (nb_iter>%d)\n",max_iter);
      break;
    }
    
  } while(converged==0);

  //free local arrays
  for(p=0;p<nb_points;p++)
    free(bprobas[p]);
  free(bprobas);
  free(total_proba);
  free(pk);
  free_distrib(global_distrib);
  return converged;
}

void dump_gem(struct gem_ws *ws,char *filename) {
  FILE *f;
  int i;
  f=fopen(filename,"w");
  fprintf(f,"%d %d\n",ws->nb_objects,ws->ds->dim);
  for(i=0;i<ws->nb_objects;i++) {
    dump_object(ws->fit_objects[i],f);
  }
  fclose(f);
}


int master_object(struct point *p,struct object **objects,int nb_objects) {
  int i;
  double distance=100000000;
  int master=-1;
  double dist;
  for(i=0;i<nb_objects;i++) {
    dist=dist_object(objects[i],p);
    if (dist<distance) {
      distance=dist;
      master=i;
    }
  }
  return master;
}

int * calc_obj_score_gem(struct gem_ws *ws) {
  int o;
  int master;
  int *np=malloc(sizeof(int)*ws->nb_objects);
  int p;

  //init result to zero
  for(o=0;o<ws->nb_objects;o++) {
    np[o]=0;
  }

  //for each point search the master
  for(p=0;p<ws->ds->nb_points;p++) {
    master=master_object(ws->ds->points[p],ws->fit_objects,ws->nb_objects);
    if (master!=-1)
      np[master]++;
  }
  return np;
}

int is_scaled_gem(struct gem_ws *ws,double scalecrit,int output) {
  double stdd;
  double limit;
  int l;
  for(l=0;l<ws->nb_objects;l++) {
    stdd=stdev_object(ws->fit_objects[l]);
    limit=ws->scale/scalecrit;
    if (output>1)
      printf("stdev_distrib=%f compared to scale/%f=%f\n",stdd,scalecrit,limit);
    if (stdd>limit) {
      if (output)
        printf("not scaled : object %d stdev %f exceed scale/%f=%f\n",l,stdd,scalecrit,limit);
      return 0;
    }
  }
  if (output)
    printf("scaled\n");
  return 1;
}


void remove_degenerated_objects_gem(struct gem_ws *ws,int seuil) {
  int o;
  int mod=1;
  int nbo;
  int *np;
  int iter=1;
  while(mod==1) {
    //printf("iteration %d\n",iter);
    mod=0;
    //update_distances_gem(ws,NULL);
    //print_degenerated_objects_gem(ws);
    np=calc_obj_score_gem(ws);
    nbo=ws->nb_objects;
    for(o=0;o<nbo;o++) {
      if (np[o]<=seuil) {
        //printf("removing object %d with score %d\n",o,np[o]);
        mod=1;
        remove_object_gem(ws,o);
        break;
      }
    }
    free(np);
    iter++;
  }
}

void print_degenerated_objects_gem(struct gem_ws *ws) {
  int o;
  int *np=calc_obj_score_gem(ws);
  for(o=0;o<ws->nb_objects;o++) {
    printf("score for object %d is %d\n",o,np[o]);
  }
  free(np);
}

void remove_dup_objects_gem(struct gem_ws *ws,double distlim) {
  int i,j;
  int rem=1;
  while(rem==1) {
    rem=0;
    for(i=0;i<ws->nb_objects;i++) {
      for(j=i+1;j<ws->nb_objects;j++) {
	if (same_object(ws->fit_objects[i],ws->fit_objects[j],distlim)) {
	  //printf("object %d is a duplicate of %d\n",j,i);
	  remove_object_gem(ws,j);
	  i=ws->nb_objects;
	  j=ws->nb_objects;
	  rem=1;
	}
      }
    }
  }
}

struct line *get_line_from_gem(struct gem_ws *ws,int objid) {

  struct gline *gl;

  if (ws->nb_objects<=objid) {
    printf("error no object with id %d\n",objid);
    return NULL;
  }
  if (ws->fit_objects[objid]->type!=LINE) {
    printf("error : object %d is not a line\n",objid);
    return NULL;
  }
  gl=(struct gline *)ws->fit_objects[objid]->container;
  return gl->line;
}

struct circle *get_circle_from_gem(struct gem_ws *ws,int objid) {

  struct gcircle *gc;

  if (ws->nb_objects<=objid) {
    printf("error no object with id %d\n",objid);
    return NULL;
  }
  if (ws->fit_objects[objid]->type!=CIRCLE) {
    printf("error : object %d is not a circle\n",objid);
    return NULL;
  }
  gc=(struct gcircle *)ws->fit_objects[objid]->container;
  return gc->circle;
}


