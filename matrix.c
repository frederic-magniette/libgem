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

//*********************** MATRIX ***************


double **new_dtab(int nb_cols,int nb_lines) {
  int i;
  double **result=malloc(sizeof(double *)*nb_cols);
  for(i=0;i<nb_cols;i++)
    result[i]=malloc(sizeof(double)*nb_lines);
  return result;
}

void free_dtab(double **dtab,int nb_cols,int nb_lines) {
  int i;
  for(i=0;i<nb_cols;i++) {
    free(dtab[i]);
  }
  free(dtab);
}

struct matrix *new_matrix(int nb_cols,int nb_lines) {
  struct matrix *result=malloc(sizeof(struct matrix));
  result->nb_cols=nb_cols;
  result->nb_lines=nb_lines;
  result->coefs=new_dtab(nb_cols,nb_lines);
  return result;
}

void free_matrix(struct matrix *m) {
  free_dtab(m->coefs,m->nb_cols,m->nb_lines);
  free(m);
}

void print_matrix(struct matrix *m) {
  int i,j;
  printf("matrix[%d][%d]=\n",m->nb_cols,m->nb_lines);
  for(i=0;i<m->nb_cols;i++) {
    for(j=0;j<m->nb_lines;j++) {
      printf("%f ",m->coefs[i][j]);
    }
    printf("\n");
  }
}

void set_col_matrix(struct matrix *m,struct vector *v,int column) {
  int i;
  for(i=0;i<m->nb_lines;i++) {
    m->coefs[column][i]=v->coords[i];
  }
}

struct vector *mult_vector_matrix(struct matrix *m,struct vector *v) {
  struct vector *result=new_vector(v->dim);
  int i,j;
  for(i=0;i<m->nb_lines;i++) {
    for(j=0;j<m->nb_cols;j++) {
      result->coords[i]+=v->coords[j]*m->coefs[j][i];
    }
  }
  return result;
}

struct point *mult_point_matrix(struct matrix *m,struct point *p) {
  struct point *result=new_point(p->dim);
  int i,j;
  for(i=0;i<m->nb_lines;i++) {
    for(j=0;j<m->nb_cols;j++) {
      result->coords[i]+=p->coords[j]*m->coefs[j][i];
    }
  }
  return result;
}

struct vector *col_vector_matrix(struct matrix *m,int col) {
  int i;
  struct vector *result=new_vector(m->nb_lines);
  for(i=0;i<m->nb_lines;i++) {
    result->coords[i]=m->coefs[col][i];
  }
  return result;
}

void dump_base_matrix(struct matrix *base,struct point *orig,char *filename) {
  int i;
  struct vector *v;
  struct point *p;
  FILE *f=fopen(filename,"w");
  dump_point(orig,f);
  for(i=0;i<base->nb_cols;i++) {
    p=cp_point(orig);
    v=col_vector_matrix(base,i);
    add_vector_point(p,v);
    dump_point(p,f);
    free_vector(v);
    free_point(p);
  }
  fclose(f);
}


