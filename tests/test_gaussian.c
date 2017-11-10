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
  struct distrib *pdist=new_distrib(CWGAUSS);
  add_data_distrib(pdist,1,4);
  add_data_distrib(pdist,2,2);
  add_data_distrib(pdist,0,1);

  printf("cpg : \n");
  printf("mean=%f\n",mean_distrib(pdist));
  printf("std dev=%f\n",stdev_distrib(pdist));
  dump_distrib(pdist,0,-4,4,200,"/tmp/test_pdist.dat");
  print_distrib(pdist);
  return 0;
}
