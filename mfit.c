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
  int o=-1;
  struct gem_ws *result;

  if (argc!=7) {
    printf("%d args given (6 wanted)\n",argc-1);
    printf("usage : %s datafile dim type[line|circle] convcrit scalecrit output[detail|summary|none]\n",argv[0]);
    return -1;
  }
  int dim=atoi(argv[2]);
  struct dataset *ds=new_dataset_fromfile(dim,argv[1]);
  struct graphics *gws=NULL;
  int type=-1;
  double convcrit=atof(argv[4]);
  double scalecrit=atof(argv[5]);

  if (!strcmp(argv[3],"circle"))
    type=CIRCLE;
  if (!strcmp(argv[3],"line"))
    type=LINE;
  if (type==-1) {
    printf("unknown type %s\n",argv[3]);
    exit(1);
  }
  
  o=param2verb(argv[6]);

  if (o!=0)
    gws=new_graphics(dim,700,700,ds->point_min,ds->point_max,o);
  

  if (ds==NULL) {
    printf("error reading file\n");
    return -1;
  }
  result=multifit(ds,convcrit,scalecrit,type,o,gws);
  printf("removing degenerated objects\n");
  remove_degenerated_objects_gem(result,4);
  print_gem(result);
  if (o!=0) {
    plot_gem(result,gws);
    getchar();
  }
  dump_gem(result,"/tmp/mfit_res.txt");
  free_dataset(ds);
  free_graphics(gws);
  free_gem(result);
  return 1;
}

