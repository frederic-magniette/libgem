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

//***************** BOX-MULLER GENERATOR *************************

struct normal_gene *new_normal_gene() {
  struct normal_gene *result=malloc(sizeof(struct normal_gene));
  result->dispo=0;
  //srand(time(NULL));
  return result;
}

void free_normal_gene(struct normal_gene *gene) {
  free(gene);
}

double gener_normal(struct normal_gene *gene) {
  double x,y;
  double s=0;
  double result;
  if (gene->dispo) {
    gene->dispo=0;
    return gene->memory;
  } else {
    while ((s==0.0) || (s>=1.0)) {
      x=rand()/((double)RAND_MAX);
      y=rand()/((double)RAND_MAX);
      s=pow(x,2)+pow(y,2);
      //printf("x=%f y=%f s=%f\n",x,y,s);
    }
    result=x*sqrt(-2*log(s)/s);
    gene->memory=y*sqrt(-2*log(s)/s);
    gene->dispo=1;
    return result;
  }
}
