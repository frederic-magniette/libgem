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

void print_neighbouring(struct neighbouring *n) {
  int i;
  printf("neighbouring of size %d center on ",n->size);
  print_point(n->center);
  for(i=0;i<n->size;i++) {
    printf("index %d with distance %f ",n->neighbours[i]->ds_index,n->neighbours[i]->dist);
    print_point(n->neighbours[i]->p);
  }
  printf("\n");
}

struct neighbouring *new_neighbouring(struct dataset *ds,struct point *center) {
  int i;
  int j;
  int k;
  double dist;
  struct neighbouring *result=malloc(sizeof(struct neighbouring));
  result->center=center;
  result->size=ds->nb_points-1;
  result->neighbours=malloc(result->size*sizeof(struct neighbour *));
  for(i=0;i<result->size;i++) {
    result->neighbours[i]=NULL;
  }
  for(i=0;i<ds->nb_points;i++) {
    if (ds->points[i]==center) {
      //printf("center\n");
      continue;
    }
    dist=dist_points(center,ds->points[i]);
    for(j=0;j<result->size;j++) {
      if (result->neighbours[j]==NULL) {
	//printf("filling slot %d\n",j);
	result->neighbours[j]=malloc(sizeof(struct neighbour));
	result->neighbours[j]->p=ds->points[i];
	result->neighbours[j]->dist=dist;
	result->neighbours[j]->ds_index=i;
	break;
      } else {
	if (result->neighbours[j]->dist<dist)
	  continue;
	//this is the place we want to place the new neighbour thus moving everybody on the right
	//printf("wanting the slot %d\n",j);
	for(k=result->size-1;k>j;k--) {
	  //printf("moving slot %d to %d\n",k-1,k);
	  result->neighbours[k]=result->neighbours[k-1];
	}
	result->neighbours[j]=malloc(sizeof(struct neighbour));
	result->neighbours[j]->p=ds->points[i];
	result->neighbours[j]->dist=dist;
	result->neighbours[j]->ds_index=i;
	break;
      }
    }
  }
  return result;
}

void free_neighbouring(struct neighbouring *n) {
  int i;
  for(i=0;i<n->size;i++) {
    free(n->neighbours[i]);
  }
  free(n->neighbours);
  free(n);
}
