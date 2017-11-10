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


int main(int argc,char **argv) {
  double dispersion;
  
  if (argc<4) {
    printf("usage : %s datafile dim dispersion resultfilename\n",argv[0]);
    return -1;
  }
  dispersion=atof(argv[3]);
  struct dataset *ds=new_dataset_fromfile(atoi(argv[2]),argv[1]);
  if (ds==NULL) {
    printf("error reading file\n");
    return -1;
  }
  noise_dataset(ds,dispersion);
  dump_dataset(ds,argv[4]);
  free_dataset(ds);
  return 1;
}
