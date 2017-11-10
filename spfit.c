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

struct vpoint {
  int id;
  double ldist;
};


void sort_vpoints(struct vpoint **vps,int nbp) {
  //sort the data structure on longitudinal distance
  struct vpoint *tmp;
  int ch=1;
  int i;
  while(ch==1) {
    ch=0;
    for(i=0;i<nbp-1;i++) {
      if (vps[i]->ldist>vps[i+1]->ldist) {
        tmp=vps[i];
        vps[i]=vps[i+1];
        vps[i+1]=tmp;
        ch=1;
      }      
    }
  }
}

int main(int argc,char **argv) {
  struct dataset *ds;
  struct dataset *ds_proj;
  struct gem_ws *emws_3d;
  struct gem_ws *emws_2d;
  struct graphics *gws_3d; 
  struct graphics *gws_2d; 
  struct line *support;
  struct circle *fit;
  int i;
  struct vpoint **vps;
  double *ldist;
  double cbase[2]={1,0};
  struct vector *vbase=new_valued_vector(2,cbase);
  struct vector *pvp;
  struct vector *pvp_old;
  struct matrix *base;
  int nbp;
  double ang_offset;
  double ang_speed;
  struct spiral *result;
  double angle;
  double asl;
  int res;
  struct distrib *aspeed;
  char name[1024];
  char title[1024];

  if (argc!=2) {
    printf("usage %s data_filename\n",argv[0]);
    exit(1);
  }

  //0. import data
  printf("step 0 : getting spiral data\n");
  ds=new_dataset_fromfile(3,argv[1]);
  nbp=ds->nb_points;

  //plot data
  gws_3d=new_graphics(3,700,700,ds->point_min,ds->point_max);
  plot_dataset(ds,gws_3d);
  getchar();
  

  //1. obtain support by gem
  printf("step 1 : obtain spiral support by gem\n");
  emws_3d=init_gem(ds,0.001);
  add_object_gem(emws_3d,LINE,NULL);
  res=algo_gem(emws_3d,NULL);
  if (res==0)
    printf("gem algo has diverged\n");
  support=get_line_from_gem(emws_3d,0);

  //plot support 
  plot_gem(emws_3d,gws_3d);
  //plot_point(support->ref,0,gws_3d);
  print_line(support);
  for(i=0;i<ds->nb_points;i++) {
    //printf("for point %d dist=%f\n",i,dist_line_point(support,ds->points[i]));
  }
  getchar(); 

  //2. project on orthogonal plan
  printf("step 2 : project points on orthogonal plan\n");
  normalize_vector(support->dir_vect);
  ldist=malloc(sizeof(double)*nbp);
  base=gram_schmidt(support->dir_vect);
  ds_proj=ortho_projection_dataset(ds,base,support->dir_vect,support->ref,ldist);

  //plot projection
  printf("proj base=\n");
  print_matrix(base);
  gws_2d=new_graphics(2,500,500,ds_proj->point_min,ds_proj->point_max);
  plot_dataset(ds_proj,gws_2d);
  get_plotfile(gws_3d->gp,name);
  dump_base_matrix(base,support->ref,name);
  sprintf(title,"ortho base");
  plot_file(gws_3d->gp,3,name,0,1,title);
  getchar();

  //3. circular fit of circle
  printf("step 3 : fit circle on projection with gem\n");
  emws_2d=init_gem(ds_proj,0.001);
  add_object_gem(emws_2d,CIRCLE,NULL);
  algo_gem(emws_2d,NULL);
  fit=get_circle_from_gem(emws_2d,0);

  //plot the circle fit
  print_circle(fit);
  plot_gem(emws_2d,gws_2d);
  plot_point(fit->center,1,gws_2d);
  apply_splot(gws_2d->sp);
  getchar();

  //4.fill and sort data structure 
  printf("step 4. fill and sort the data structure\n");
  vps=malloc(sizeof(struct vpoint *)*nbp);
  for(i=0;i<ds_proj->nb_points;i++) {
    vps[i]=malloc(sizeof(struct vpoint));
    vps[i]->id=i;
    vps[i]->ldist=ldist[i];
  }
  free(ldist);
  sort_vpoints(vps,nbp);
  printf("\n\n");

  //5. calc angular speed
  printf("step 5. calculating angular speed\n");
  aspeed=new_distrib(WGAUSS);
  pvp_old=vector_by2points(fit->center,ds_proj->points[vps[0]->id]);
  for(i=1;i<nbp;i++) {
    pvp=vector_by2points(fit->center,ds_proj->points[vps[i]->id]);
    angle=angle_between_2_vector(pvp_old,pvp);
    free_vector(pvp_old);
    pvp_old=pvp;
    asl=angle/(vps[i]->ldist-vps[i-1]->ldist);
    printf("for point %d, aspeed=%f\n",i,asl);
    add_data_distrib(aspeed,asl,1); //replace 1 by weight in multifit context
  }
  free_vector(pvp_old);
  ang_speed=mean_distrib(aspeed);

  //plot angular speed
  printf("angular speed is %f\n\n\n",ang_speed);
  getchar();

  //6. calc angular offset
  printf("step 6. calculating angular offset\n");
  pvp=vector_by2points(fit->center,ds_proj->points[vps[0]->id]);
  ang_offset=fmod(angle_between_2_vector(vbase,pvp)-ang_speed*vps[0]->ldist,360);
  free_vector(pvp);

  //plot angular offset
  printf("angular offset is %f\n\n\n",ang_offset);
  getchar();

  //7. create the spiral structure and plot it
  printf("step 7. plot the resulting spiral\n");
  result=new_based_spiral(support,fit->radius,ang_offset,ang_speed,base);
  print_spiral(result);
  plot_spiral(result,0,gws_3d);
  getchar();

  //free memory
  free_spiral(result);
  for(i=0;i<nbp;i++)
    free(vps[i]);
  free(vps);
  free_matrix(base);
  free_graphics(gws_2d);
  free_graphics(gws_3d);
  free_gem(emws_2d);
  free_gem(emws_3d);
  free_dataset(ds_proj);
  free_dataset(ds);
  return 0;
}
