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

//build a plot of angular resolution depending on scalecrit

int main(int argc,char **argv) {

  struct gem_ws *result;
  int angle_iter;
  struct dataset *ds;
  double noise_level=0.1;
  FILE *f;
  int scalecrit;
  int size=100;

  f=fopen("/tmp/angres.txt","w");
  
  for(scalecrit=10;scalecrit<=180;scalecrit++) {
    for(angle_iter=0;angle_iter<180;angle_iter++) {

      //generate data
      ds=new_empty_dataset(2);
      add_angle_line_dataset(ds,0,20,0,0,size);
      add_angle_line_dataset(ds,angle_iter,20,0,0,size);
      noise_dataset(ds,noise_level);
      
      //exec the multifit
      result=multifit(ds,0.0001,scalecrit,LINE,0,NULL);
      //remove_dup_objects_gem(result);
      remove_degenerated_objects_gem(result,3);

      free_dataset(ds);
      
      if (result->nb_objects>1) {
        printf("%d %d\n",scalecrit,angle_iter);
        fprintf(f,"%d %d\n",scalecrit,angle_iter);
	free_gem(result);
        break;
      }
      
      //free data
      free_gem(result);
    }
  }
  return 1;
}

