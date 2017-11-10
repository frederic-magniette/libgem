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
  int i;
  struct gem_ws *emws;
  struct graphics *gws=NULL;
  if (argc!=7) {
    printf("usage : %s datafile dim nb_lines seed convcrit output\n",argv[0]);
    return -1;
  }
  int dim=atoi(argv[2]);
  int nb_lines=atoi(argv[3]);
  struct dataset *ds=new_dataset_fromfile(dim,argv[1]);
  int seed=atoi(argv[4]);
  double convcrit=atof(argv[5]);

  int o=-1;
  if (!strcmp(argv[6],"none"))
    o=0;
  if (!strcmp(argv[6],"summary"))
    o=1;
  if (!strcmp(argv[6],"detail"))
    o=2;
  if (o==-1) {
    printf("unknown output : %s\n",argv[3]);
    exit(1);
  }

  if (ds==NULL) {
    printf("error reading file\n");
    return -1;
  }
  if (o!=0)
    gws=new_graphics(dim,700,700,ds->point_min,ds->point_max);
  srand(seed);
  emws=init_gem(ds,convcrit);
  for(i=0;i<nb_lines;i++)
    add_random_object_gem(emws,LINE);
  if (o==2)
    algo_gem(emws,gws);
  else 
    algo_gem(emws,NULL);

  if (o!=0) {
    printf("finalizing\n");
    plot_gem(emws,gws);
    getchar();
  }

  dump_gem(emws,"/tmp/lfit_res.txt");
  free_gem(emws);
  free_dataset(ds);
  free_graphics(gws);
  return 1;
}

