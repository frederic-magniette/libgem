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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <SDL/SDL.h>

#define PI 3.14159265359

#define max(a,b) ((a) > (b) ? (a) : (b)) 
#define min(a,b) ((a) < (b) ? (a) : (b)) 

struct gem_ws;
struct graphics;
struct weights;
struct matrix;
struct vector;
struct line;
struct dataset;

//points

struct point {
  int dim;
  double *coords;
};

struct point *new_point(int dim);
struct point *new_valued_point(int dim,double *values);
struct point *new_random_point(int dim,struct point *pmin,struct point *pmax);
struct point *cp_point(struct point *p);
void print_point(struct point *p);
void free_point(struct point *p);
double dist_points(struct point *p1,struct point *p2);
void plot_point(struct point *p,int id,struct graphics *gws);
void add_vector_point(struct point *p,struct vector *v);
int is_valid_point(struct point *p,struct point *pmin,struct point *pmax);
struct point *ortho_project_point(struct point *p,struct matrix *base,struct vector *dir,struct point *ref);
struct point *cut_first_coord_point(struct point *p);
struct point *change_base_point(struct point *p,struct matrix *base);
struct point *barycenter_point(struct dataset *ds);
struct point *weighted_barycenter_point(struct dataset *ds,struct weights *weights);
void dump_point(struct point *p,FILE *f);

//datasets

struct dataset {
  int dim;
  struct point **points;
  struct point *point_min;
  struct point *point_max;
  int nb_points;
};

struct dataset * new_dataset_fromfile(int dim,char *filename);
void print_dataset(struct dataset *ds);
void dump_dataset(struct dataset *ds,char *filename);
void calc_min_point_dataset(struct dataset *ds);
void calc_max_point_dataset(struct dataset *ds);
void free_dataset(struct dataset *ds);
struct dataset * new_zero_dataset(int nb_points,int dim);
struct dataset * new_empty_dataset(int dim);
void plot_weighted_dataset(struct dataset *ds,struct weights **w,int nb_lines,struct graphics *gws);
void plot_dataset(struct dataset *ds,struct graphics *gws);
void noise_dataset(struct dataset *ds,double dispersion);
void add_point_dataset(struct dataset *ds,struct point *p);
void add_line_dataset(struct dataset *ds,struct line *l,int nb_steps,double length);
void add_angle_line_dataset(struct dataset *ds,int angle,int nb_steps,double x,double y,double length);
void add_random_line_dataset(struct dataset *ds,int nb_steps,double length);
void add_random_noise_dataset(struct dataset *ds,int nb_points);
struct dataset *ortho_projection_dataset(struct dataset *ds,struct matrix *base,struct vector *dir,struct point *ref,double *long_dist);

//neighbouring

struct neighbour {
  struct point *p;
  double dist;
  int ds_index;
};

struct neighbouring {
  struct point *center;
  struct neighbour **neighbours;
  int size;
};

struct neighbouring *new_neighbouring(struct dataset *ds,struct point *center);
void print_neighbouring(struct neighbouring *n);
void free_neighbouring(struct neighbouring *n);

//angles

double deg2rad(double angle);
double rad2deg(double angle);

//weights

struct weights {
  int nb_points;
  double *coefs;
};

struct weights *new_weights(int nb_points,double init_value);
void free_weights(struct weights *p);
struct weights *cp_weights(struct weights *src);
void print_weights(struct weights *p);

//vectors

struct vector {
  int dim;
  double *coords;
};

struct vector *new_vector(int dim);
struct vector *new_orthonormal_vector(int dim,int one_index);
struct vector *new_valued_vector(int dim,double *values);
struct vector *new_random_vector(int dim,double max_value);
void free_vector(struct vector *v);
void print_vector(struct vector *v);
void mult_lambda_vect(struct vector *v,double lambda);
void add_vector(struct vector *dst,struct vector *src);
void sub_vector(struct vector *dst,struct vector *src);
struct vector *sum_vector(struct vector *v1,struct vector *v2);
struct vector *project_vector(struct vector *src,struct vector *support);
double coef_project_vector(struct vector *src,struct vector *support);
int eq_vector(struct vector *v1,struct vector *v2);
struct vector *cp_vector(struct vector *v);
struct vector *vector_by2points(struct point *origin,struct point *dest);
struct vector *vector_abs_by2points(struct point *origin,struct point *dest);
double scalar_vector(struct vector *v1,struct vector *v2);
double norme_vector(struct vector *v);
void normalize_vector(struct vector *v);
struct vector *weighted_mean_vector(struct dataset *ds,struct point *ref,struct weights *w);
void plot_vector(struct vector *v,struct point *origin,int id,struct graphics *gws);
void zero_vector(struct vector *v,double value);
int is_zero_vector(struct vector *v,double limit);
int isnan_vector(struct vector *v);
struct vector *change_base_vector(struct vector *v,struct matrix *base);
double angle_between_2_vector(struct vector *v1,struct vector *v2);
struct vector *cut_first_coord_vector(struct vector *v);

//matrices

struct matrix {
  int nb_cols;
  int nb_lines;
  double **coefs;
};

double **new_dtab(int nb_cols,int nb_lines);
struct matrix *new_matrix(int nb_cols,int nb_lines);
void free_matrix(struct matrix *m);
void print_matrix(struct matrix *m);
void set_col_matrix(struct matrix *m,struct vector *v,int column);
struct vector *mult_vector_matrix(struct matrix *m,struct vector *v);
struct point *mult_point_matrix(struct matrix *m,struct point *p);
struct vector *col_vector_matrix(struct matrix *m,int col);
void dump_base_matrix(struct matrix *base,struct point *orig,char *filename);

//algebra

struct matrix *gram_schmidt(struct vector *v);
struct point **orthonormal_projection(struct point **points,int nb_points,struct vector *v,struct point *ref);

//lines

struct line {
  int dim;
  struct point *ref;
  struct vector *dir_vect;
};

struct line *new_line(struct point *ref,struct vector *dir_vect);
struct line *cp_line(struct line *l);
void free_line(struct line *l);
void print_line(struct line *l);
double dist_line_point(struct line *l,struct point *p);
void sort_line_array(struct line **la,int la_size);
void dump_vect_line(struct line *l,char *filename);
struct point **boxed_line(struct line *l,struct point *point_min,struct point *point_max,int nb_steps);
void dump_boxed_line(struct line *l,char *filename,struct point *point_min,struct point *point_max,int nb_steps);
void dump_line(struct line *l,FILE *f);
struct point *get_coord_line(struct line *l,double val,int axis);
void plot_line(struct line *l,int id,struct graphics *gws);
struct line *angle_line(int angle,double x,double y);
struct line *random_line(int dim,struct point *pmin,struct point *pmax);

//circles

struct circle {
  int dim;
  struct point *center;
  double radius;
};

struct circle *new_circle(struct point *center,double radius);
struct circle *cp_circle(struct circle *src);
void free_circle(struct circle *c);
void print_circle(struct circle *c);
double dist_circle_point(struct circle *c,struct point *p);
void plot_circle(struct circle *c,int id,struct graphics *gws);
void dump_boxed_circle(struct circle *c,char *filename,struct point *point_min,struct point *point_max,int nb_steps);
void dump_circle(struct circle *c,FILE *f);

//spirals

struct spiral {
  int dim;
  struct line *support;
  struct matrix *proj_base;
  double radius;
  double angular_offset;
  double angular_speed;
};


struct spiral *new_spiral(struct line *support,double radius,double angular_offset,double angular_speed);
struct spiral *new_based_spiral(struct line *support,double radius,double angular_offset,double angular_speed,struct matrix *base);
void free_spiral(struct spiral *s);
void print_spiral(struct spiral *s);
void dump_boxed_spiral(struct spiral *s,char *filename,struct point *point_min,struct point *point_max,int nb_steps);
void plot_spiral(struct spiral *s,int id,struct graphics *gws);

//statistic distributions

struct gauss {
  double nb_data;
  double old_mean;
  double mean;
  double s;
  double old_s;
  double variance;
  double stdev;
};

struct gauss *new_gauss();
void add_data_gauss(struct gauss *dist,double x);
void dump_gauss(struct gauss *dist,double xmin,double xmax,int nb_steps,char *filename);

struct wgauss {
  double sum_weights;
  double old_mean;
  double mean;
  double old_s;
  double s;
  double variance;
  double stdev;
};

struct wgauss *new_wgauss();
struct wgauss *cp_wgauss(struct wgauss *dist);
void free_wgauss(struct wgauss *dist);
void add_data_wgauss(struct wgauss *dist,double x,double weight);
void print_wgauss(struct wgauss *dist);
void dump_wgauss(struct wgauss *dist,double xmin,double xmax,int nb_steps,char *filename);
void reinit_wgauss(struct wgauss *dist);
double likelyhood_wgauss(struct wgauss *dist,double x);
double normalized_likelyhood_wgauss(struct wgauss *dist,double x);
double normalized_centred_likelyhood_wgauss(double stdev,double x);

struct cwgauss {
  double sum_weights;
  double old_stdev;
  double mean;
  double stdev;
};

struct cwgauss *new_cwgauss();
struct cwgauss *cp_cwgauss(struct cwgauss *dist);
void free_cwgauss(struct cwgauss *dist);
void reinit_cwgauss(struct cwgauss *dist);
void print_cwgauss(struct cwgauss *dist);
void add_data_cwgauss(struct cwgauss *dist,double x,double weight);
double normalized_likelyhood_cwgauss(struct cwgauss *dist,double x);
void dump_cwgauss(struct cwgauss *dist,double xmin,double xmax,int nb_steps,char *filename);

#define WGAUSS 1
#define CWGAUSS 2

struct distrib {
  int type;
  void *data;  
};

struct distrib *new_distrib(int type);
void add_data_distrib(struct distrib *d,double value,double weight);
struct distrib *cp_distrib(struct distrib *d);
void free_distrib(struct distrib *d);
double normalized_likelyhood_distrib(struct distrib *d,double x);
void print_distrib(struct distrib *d);
double mean_distrib(struct distrib *d);
double stdev_distrib(struct distrib *d);
void set_stdev_distrib(struct distrib *d,double stdev);
void dump_distrib(struct distrib *d,double center,double xmin,double xmax,int nb_steps,char *filename);
void reinit_distrib(struct distrib *d);

//generator of gaussian data

struct normal_gene {
  char dispo;
  double memory;
};

struct normal_gene *new_normal_gene();
double gener_normal(struct normal_gene *gene);
void free_normal_gene(struct normal_gene *gene);

//gnuplot

struct gplot {
  int first_plot;
  int fn_index;
  FILE *pipe;
};

struct gplot *new_gplot(int w,int l);
void plot_file(struct gplot *gp,int dim,char *filename,char withlines,char replot,char *title);
void plot_histo_file(struct gplot *gp,double xmin,double xmax,double box_large,char *filename,char replot,char *title);
void replot(struct gplot *gp);
void get_plotfile(struct gplot *gp,char *filename);
void free_gplot(struct gplot *gp);

//SDL

struct splot {
  struct point *point_min;
  struct point *point_max;
  int w,h;
  double margin;
  SDL_Surface *screen;
};

struct splot *new_splot(int w,int h,struct point *point_min,struct point*point_max);
void attribute_color_splot(int id,int *r,int *v,int *b);
void plot_real_point_splot(struct point *src,struct splot *ws,int r,int g,int b);
void plot_file_splot(struct splot *ws,char *filename,int r,int g,int b);
void apply_splot(struct splot *ws);
void clear_splot(struct splot *ws);
void free_splot(struct splot *ws);
struct point *sdl2real_splot(int i,int j,struct splot *ws);
struct point *real2sdl_splot(double x,double y,struct splot *ws);
void set_pixel_splot(struct splot *sp,int x,int y,int r,int g,int b);
void get_pixel_splot(struct splot *sp,int x,int y,int *r,int *g,int *b);

struct graphics {
  int dim;
  struct point *point_min;
  struct point *point_max;
  int w;
  int h;
  struct gplot *gp;
  struct splot *sp;
};

struct graphics *new_graphics(int dim,int w,int h,struct point *point_min,struct point*point_max);
void free_graphics(struct graphics *gws);
void apply_graphics(struct graphics *gws);

struct gline {
  struct line *line; //the fit lines
  struct distrib *pgauss; //a weightsated gaussian
};

//probabilized lines

struct gline *new_gline(struct dataset *ds,struct weights *w);
//struct gline *new_precomputed_gline(int dim,struct point *ref,struct vector *dir_vect,struct dataset *ds,struct weights *w);
struct gline *cp_gline(struct gline *gl);
void print_gline(struct gline *gl);
void dump_gline(struct gline *gl,FILE *f);
void dump_boxed_gline(struct gline *gl,char *filename,struct point *point_min,struct point*point_max,int nb_steps);
void free_gline(struct gline *gl);
double belonging_proba_gline(struct gline *gl,struct point *p);
double dist_gline(struct gline *gl,struct point *p);
void plot_gline(struct gline *gl,int id,struct graphics *gws);
void plot_field_gline(struct gline *gl,int id,struct graphics *gws);
int same_gline(struct gline *gl1,struct gline *gl2,double distlim);

//probabilized circles

struct gcircle {
  struct circle *circle;
  struct distrib *pgauss;
};

struct gcircle *new_gcircle(struct dataset *ds,struct weights *w);
double belonging_proba_gcircle(struct gcircle *gc,struct point *p);
void free_gcircle(struct gcircle *gc);
void print_gcircle(struct gcircle *gc);
void dump_gcircle(struct gcircle *gc,FILE *f);
struct gcircle *cp_gcircle(struct gcircle *gc);
double dist_gcircle(struct gcircle *gc,struct point *p);
void plot_gcircle(struct gcircle *gc,int id,struct graphics *gws);
void plot_field_gcircle(struct gcircle *gc,int id,struct graphics *gws);
int same_gcircle(struct gcircle *gc1,struct gcircle *gc2,double distlim);

//objects

enum otype {LINE,CIRCLE};

struct object {
  enum otype type;
  void *container;
};

double belonging_proba_object(struct object *o,struct point *p);
double dist_object(struct object *o,struct point *p);
struct object *new_object(enum otype type,struct dataset *ds,struct weights *w);
struct object *cp_object(struct object *src);
void free_object(struct object *o);
void plot_object(struct object *o,int id,struct graphics *gws);
void plot_field_object(struct object *o,int id,struct graphics *gws);
void print_object(struct object *o);
void dump_object(struct object *o,FILE *f);
double stdev_object(struct object *o);
void set_stdev_object(struct object *o,double stdev);

//em algorithm

struct gem_ws {
  int nb_objects;
  double scale;
  struct weights **fit_weights; //weights of points
  struct object **fit_objects; //the fit gaussian objects
  struct dataset *ds;
  int nb_iter;
  double convergence_criterion;
  double **distance;
  double global_dist;
};

void free_gem(struct gem_ws *ws);
struct gem_ws *init_gem(struct dataset *ds,double convcrit);
int algo_gem(struct gem_ws *ws,struct graphics *gws);
void em(struct dataset *ds,int nb_lines,double convcrit);
int add_object_gem(struct gem_ws *ws,enum otype type,struct graphics *gws);
int add_random_object_gem(struct gem_ws *ws,enum otype type);
void remove_object_gem(struct gem_ws *ws,int obj_idx);
int same_object(struct object *o1,struct object *o2,double distlim);
void print_gem(struct gem_ws *ws);
void plot_gem(struct gem_ws *ws,struct graphics *gws);
struct gem_ws *cp_gem(struct gem_ws *src);
void dump_gem(struct gem_ws *ws,char *filename);
int is_scaled_gem(struct gem_ws *ws,double scalecrit,int output);
void remove_degenerated_objects_gem(struct gem_ws *ws,int seuil);
void print_degenerated_objects_gem(struct gem_ws *ws);
void remove_dup_objects_gem(struct gem_ws *ws,double distlim);
struct line *get_line_from_gem(struct gem_ws *ws,int objid);
struct circle *get_circle_from_gem(struct gem_ws *ws,int objid);


//mfit algorithm

struct gem_ws *multifit(struct dataset *ds,double convcrit,double scalecrit,enum otype method,int output,struct graphics *gws);

//ufit algorithm

struct ufit_tree {
  struct ufit_elem *root;
  struct ufit_elem *current;
  struct ufit_elem *best;
  struct ufit_elem *first_son;
  struct ufit_elem *chain;
  int tfs;
  int level;
};

struct ufit_elem {
  struct gem_ws *gem;
  char done;
  double score;
  int method;
  struct ufit_elem *parent;
  struct ufit_elem *next;
  struct ufit_elem *sons;
};

#define NB_METHODS 2

struct ufit_tree *ufit(struct dataset *ds,double convcrit,double scalecrit,int output,struct graphics *gws);
void free_ufit_tree(struct ufit_tree *tree);
void print_ufit_tree(struct ufit_tree *tree);
