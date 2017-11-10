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

//*********************** WEIGHTS *****************************

struct weights *new_weights(int nb_points,double init_value) {
  int i;
  struct weights *result=malloc(sizeof(struct weights));
  result->nb_points=nb_points;
  result->coefs=malloc(nb_points*sizeof(double));
  for(i=0;i<nb_points;i++)
    result->coefs[i]=init_value;
  return result;
}

void free_weights(struct weights *p) {
  if (p) {
    free(p->coefs);
    free(p);
  }
}

struct weights *cp_weights(struct weights *src) {
  int i;
  struct weights *result=malloc(sizeof(struct weights));
  result->nb_points=src->nb_points;
  result->coefs=malloc(src->nb_points*sizeof(double));
  for(i=0;i<src->nb_points;i++)
    result->coefs[i]=src->coefs[i];
  return result;
}

void print_weights(struct weights *p) {
  int i;
  printf("weights[%d]=\n",p->nb_points);
  for(i=0;i<p->nb_points;i++) {
    printf("%2lf\t",p->coefs[i]);
    if ((i+1)%5==0)
      printf("\n");
  }
  printf("\n");
}
