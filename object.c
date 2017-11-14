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

double belonging_proba_object(struct object *o,struct point *p) {
  switch(o->type) {
  case LINE:
    return belonging_proba_gline((struct gline *)o->container,p);
  case CIRCLE:
    return belonging_proba_gcircle((struct gcircle *)o->container,p);
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}

double dist_object(struct object *o,struct point *p) {
  switch(o->type) {
  case LINE:
    return dist_gline((struct gline *)o->container,p);
  case CIRCLE:
    return dist_gcircle((struct gcircle *)o->container,p);
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}

struct object *new_object(enum otype type,struct dataset *ds,struct weights *w,struct graphics *gws) {
  struct object *result=malloc(sizeof(struct object));
  result->type=type;
  switch(type) {
  case LINE:
    result->container=(void *)new_gline(ds,w,gws);
    break;
  case CIRCLE:
    result->container=(void *)new_gcircle(ds,w,gws);
    break;
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
  return result;
}

struct object *cp_object(struct object *src) {
  struct object *result=malloc(sizeof(struct object));
  result->type=src->type;
  switch(src->type) {
  case LINE:
    result->container=(void *)cp_gline((struct gline *)src->container);
    break;
  case CIRCLE:
    result->container=(void *)cp_gcircle((struct gcircle *)src->container);
    break;
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
  return result;
}

void free_object(struct object *o) {
  if (!o)
    return;
  switch(o->type) {
  case LINE:
    free_gline((struct gline *)o->container);
    break;
  case CIRCLE:
    free_gcircle((struct gcircle *)o->container);
    break;
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
  free(o);
}

void plot_object(struct object *o,int id,struct graphics *gws) {
   switch(o->type) {
  case LINE:
    return plot_gline((struct gline *)o->container,id,gws);
  case CIRCLE:
    return plot_gcircle((struct gcircle *)o->container,id,gws);
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}

void plot_field_object(struct object *o,int id,struct graphics *gws) {
   switch(o->type) {
  case LINE:
    return plot_field_gline((struct gline *)o->container,id,gws);
  case CIRCLE:
    return plot_field_gcircle((struct gcircle *)o->container,id,gws);
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}

void print_object(struct object *o) {
  switch(o->type) {
  case LINE:
    return print_gline((struct gline *)o->container);
  case CIRCLE:
    return print_gcircle((struct gcircle *)o->container);
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}

void dump_object(struct object *o,FILE *f) {
  switch(o->type) {
  case LINE:
    fprintf(f,"l ");
    dump_gline((struct gline *)o->container,f);
    break;
  case CIRCLE:
    fprintf(f,"c ");
    dump_gcircle((struct gcircle *)o->container,f);
    break;
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}

double stdev_object(struct object *o) {
  switch(o->type) {
  case LINE:
    return stdev_distrib(((struct gline *)o->container)->pgauss);
  case CIRCLE:
    return stdev_distrib(((struct gcircle *)o->container)->pgauss);
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}

void set_stdev_object(struct object *o,double stdev)  {
  switch(o->type) {
  case LINE:
    set_stdev_distrib(((struct gline *)o->container)->pgauss,stdev);
    break;
  case CIRCLE:
    set_stdev_distrib(((struct gcircle *)o->container)->pgauss,stdev);
    break;
  default:
    printf("unknown object type %d for set_stdev...exiting\n",o->type);
    exit(1);
  }
}

int same_object(struct object *o1,struct object *o2,double distlim) {
  if (o1->type!=o2->type)
    return 0;
  switch(o1->type) {
  case LINE:
    return same_gline(((struct gline *)o1->container),((struct gline *)o2->container),distlim);
  case CIRCLE:
    return same_gcircle(((struct gcircle *)o1->container),((struct gcircle *)o2->container),distlim);
  default:
    printf("unknown object type...exiting\n");
    exit(1);
  }
}
