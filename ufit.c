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

//**************************** MAIN ****************************

int main(int argc,char **argv) {


  struct ufit_tree *result;
  if (argc!=6) {
    printf("usage : %s datafile dim convcrit scalecrit output\n",argv[0]);
    return -1;
  }
  int dim=atoi(argv[2]);
  double convcrit=atof(argv[3]);
  double scalecrit=atof(argv[4]);
  struct dataset *ds=new_dataset_fromfile(dim,argv[1]);
  struct graphics *gws=NULL;
  struct gem_ws *result_gem;

  int o=param2verb(argv[5]);

  if (o!=0)
    gws=new_graphics(dim,700,700,ds->point_min,ds->point_max,o); 

  if (ds==NULL) {
    printf("error reading file\n");
    return -1;
  }
  result=ufit(ds,convcrit,scalecrit,0,gws);
  printf("removing degenerated objects\n");
  print_ufit_tree(result);

  result_gem=result->best->gem;
  dump_gem(result_gem,"/tmp/ufit_res.txt");
  remove_degenerated_objects_gem(result_gem,4);
  if (o!=0) {
    plot_gem(result_gem,gws);
    getchar();
  }
  free_ufit_tree(result);
  free_dataset(ds);
  free_graphics(gws);
  return 1;
}

