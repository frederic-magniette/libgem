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

int main() {
  struct dataset*ds;
  struct neighbouring *n;

  ds=new_dataset_fromfile(2,"data/test_neigh.txt");
  if (ds==NULL) {
    printf("unable to load data\n");
    exit(1);
  }
  print_dataset(ds);
  n=new_neighbouring(ds,ds->points[0]);
  print_neighbouring(n);
  free_neighbouring(n);
  return 1;
}
