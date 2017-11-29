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

//**************************** MULTIFIT ****************************

struct gem_ws *multifit(struct dataset *ds,double convcrit,double scalecrit,enum otype method,int output,struct graphics *gws) {

  struct gem_ws *gemws;
  gemws=init_gem(ds,convcrit);
  int res;
  int nb_objects=0;
  
  while(1) {

    //adding a new object
    if (gws)
      printf("adding an object with method %d\n",method);
    nb_objects++;
    if (nb_objects>10) {
      //printf("enough objects : stopping\n");
      break;
    }
    switch(method) {
    case 0:
      add_object_gem(gemws,LINE,gws);
      break;
    case 1:
      add_object_gem(gemws,CIRCLE,gws);
      break;
    default:
      printf("unknown method %d\n",method);
      return NULL;
    }

    if (output!=0) {
      plot_gem(gemws,gws);
      getchar();
    }
    
    //apply the gem algorithm
    if (output==2)
      res=algo_gem(gemws,gws);
    else
      res=algo_gem(gemws,NULL);

    //check if gem has diverged
    if (res==0) {
      if (output!=0) {
	printf("gem has diverged\n");
        plot_gem(gemws,gws);
        getchar();
      }
    }

    //continue if gem is not scaled
    if (!is_scaled_gem(gemws,scalecrit,output))
      continue;
    else
      break;
  }

  if (output!=0) {
    printf("finalizing\n");
    plot_gem(gemws,gws);
    getchar();
  }
  return gemws;
}


//******************************* UFIT *******************************

void print_method_tree(struct ufit_elem *current) {
  if (current->parent)
    print_method_tree(current->parent);
  switch(current->method) {
  case 0:
    printf("L");
    break;
  case 1:
    printf("C");
    break;
  default:
    break;
  }
}

int add_object_ufit(struct gem_ws *gemws,int method) {
  int res;
  switch(method) {
  case 0:
    //printf("ufit add a line\n");
    res=add_object_gem(gemws,LINE,NULL);
    break;
  case 1:
    //printf("ufit add a circle\n");
    res=add_object_gem(gemws,CIRCLE,NULL);
    break;
  default:
    printf("unknown method %d\n",method);
    res=-1;
  }
  return res;
}

void print_elem_ufit(struct ufit_tree *tree,struct ufit_elem *el,int level) {
  int i;
  if (el->done==0)
    return;
  for(i=0;i<level;i++) {
    printf("-|");
  }
  print_method_tree(el);
  printf("(%.2f)",el->score);
  
  if (el==tree->best) {
    printf("  !!BEST!!\n");
  }  else {
    printf("\n");
  }

}

void print_eu_score(struct ufit_tree *tree,struct ufit_elem *el,int level) {
  struct ufit_elem *sons;
  print_elem_ufit(tree,el,level);
  sons=el->sons;
  while(sons!=NULL && sons->parent==el) {
    print_eu_score(tree,sons,level+1);
    sons=sons->next;
  }
}

void print_ufit_tree(struct ufit_tree *tree) {
  printf("choices(score)\n");
  printf("root");
  print_eu_score(tree,tree->root,0);
}

void new_elem_ufit(struct ufit_tree *tree,int method) {
  struct ufit_elem *ures;

  //create a new node in the tree
  ures=malloc(sizeof(struct ufit_elem));
  ures->gem=cp_gem(tree->current->gem);
  ures->score=10000000;
  ures->done=0;
  ures->parent=tree->current;
  if (tree->current->sons==NULL) 
    tree->current->sons=ures;
  ures->next=NULL;
  ures->sons=NULL;
  ures->method=method;
  
  //add an object to this gem
  add_object_ufit(ures->gem,method);

  //chain the new elem
  if (tree->tfs==1) {
    //printf("first son\n");
    tree->first_son=ures;
    tree->chain=ures;
    tree->tfs=0;
  } else {
    tree->chain->next=ures;
    tree->chain=ures;
  }
}

void next_current(struct ufit_tree *tree,int found_one,int directcut) {
  //printf("\nsearching next current elem\n");
  if (tree->current->next!=NULL) {
    if (found_one && directcut) {
      tree->current=NULL;
    } else {
      //printf("going to next gem in the current level %d\n",tree->level);
      tree->current=tree->current->next;
    }
  } else {
    if (found_one) {
      tree->current=NULL;
    } else {
      tree->level++;
      //printf("going to next level %d\n",tree->level);
      tree->tfs=1;
      tree->current=tree->first_son;
      tree->first_son=NULL;
    }
  }
}

struct ufit_tree *new_ufit_tree(struct dataset *ds,double convcrit) {
  struct ufit_tree *result=malloc(sizeof(struct ufit_tree));
  memset(result,0,sizeof(struct ufit_tree));
  result->root=malloc(sizeof(struct ufit_elem));
  result->root->gem=init_gem(ds,convcrit);
  result->root->score=1000000;
  result->root->parent=NULL;
  result->root->method=-1;
  result->root->next=NULL;
  result->root->sons=NULL;
  result->best=result->root;
  result->current=result->root;
  result->first_son=NULL;
  result->chain=NULL;
  result->tfs=1;
  result->level=1;
  return result;
}

void free_elem_ufit(struct ufit_elem *root) {
  struct ufit_elem *iter;
  struct ufit_elem *todelete;
  if (root->sons)
    free_elem_ufit(root->sons);
  iter=root;
  while(iter!=NULL) {
    todelete=iter;
    iter=iter->next;
    free_gem(todelete->gem);
    free(todelete);
  }
}

void free_ufit_tree(struct ufit_tree *tree) {
  free_elem_ufit(tree->root);
  free(tree);
}

struct ufit_tree *ufit(struct dataset *ds,double convcrit,double scalecrit,int output,struct graphics *gws) {
  struct ufit_tree *tree=new_ufit_tree(ds,convcrit);
  int method;
  int gem_res;
  int found_one=0;
  
  do { 

    //apply gem algo on current
    if (tree->current->gem->nb_objects>0) {
      if (output==2)
        gem_res=algo_gem(tree->current->gem,gws);
      else
        gem_res=algo_gem(tree->current->gem,NULL);
      
      //if algo has converged and score is better
      tree->current->score=tree->current->gem->global_dist;
      if ((gem_res!=0) && (tree->current->score<tree->best->score)) {
        printf("best result so far %f ",tree->current->gem->global_dist);
        tree->best=tree->current;
        if (is_scaled_gem(tree->current->gem,scalecrit,output)) {
          found_one=1;
          printf("scaled\n");
        } else {
          printf("not scaled\n");
        }
      } 
    }
    tree->current->done=1;

    print_elem_ufit(tree,tree->current,tree->level);

    //plot gem result
    if (output!=0) {
      if (tree->current!=tree->root) {
        plot_gem(tree->current->gem,gws);
        getchar();
      }
    }

  
    //if necessary create the sons of the current state
    for(method=0;method<NB_METHODS;method++) {
      new_elem_ufit(tree,method);
    }
 
    //go to next current
    //last param is directcut : 1 if you want to stop after the first best scaled, 0 otherwise
    next_current(tree,found_one,1);
    
    //getchar();

  } while (tree->current);
  
  printf("\n\n\nFinalizing\n");
  if (tree->best!=tree->root) {
    if (output!=0) {
      plot_gem(tree->best->gem,gws);
      print_gem(tree->best->gem);
      getchar();
    }
  } else {
    printf("no best elem selectioned by ufit\n");

  }
  return tree;
}
