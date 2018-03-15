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
  char cmd[4096];
  struct dataset *ds;
  double noise_level=0.1;
  FILE *f;
  int scalecrit;

  f=fopen("/tmp/angul_resol.txt","w");
  
  for(scalecrit=10;scalecrit<=100;scalecrit++) {
    for(angle_iter=0;angle_iter<180;angle_iter++) {
      //printf("\nscalecrit=%d angle=%d\n",scalecrit,angle_iter);
      
      //generate data
      system("../data/angle_line.exe 2 10 l1.txt 0 -100 -100 100 100");
      sprintf(cmd,"../data/noise.exe l1.txt 2 %f n1.txt",noise_level);
      system(cmd);
      sprintf(cmd,"../data/angle_line.exe 2 10 l2.txt %d -100 -100 100 100",angle_iter);
      system(cmd);
      sprintf(cmd,"../data/noise.exe l2.txt 2 %f n2.txt",noise_level);
      system(cmd);
      system("cat n1.txt n2.txt > test.txt");
      
      //read the dataset
      ds=new_dataset_fromfile(2,"test.txt");
      if (ds==NULL) {
        printf("error reading file\n");
        return -1;
      }
      
      //exec the multifit
      result=multifit(ds,0.0001,scalecrit,LINE,0,NULL);
      //remove_dup_objects_gem(result);
      remove_degenerated_objects_gem(result,3);

      free_dataset(ds);
      system("rm -f l1.txt l2.txt n1.txt n2.txt test.txt");
      
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

