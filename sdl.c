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
#include <signal.h>

void init_splot() {
}

struct splot *new_splot(int w,int h,struct point *point_min,struct point*point_max) {
  struct splot *result=malloc(sizeof(struct splot));
  SDL_Rect rect;
  Uint32 color;
  int i;
  double diffx,diffy;

  //init SDL 1.2
  SDL_Init(SDL_INIT_VIDEO|SDL_INIT_NOPARACHUTE);
  signal(SIGINT,SIG_DFL);

  result->point_min=cp_point(point_min);
  result->point_max=cp_point(point_max);
  
  //orthonormalize the frame
  diffx=point_max->coords[0]-point_min->coords[0];
  diffy=point_max->coords[1]-point_min->coords[1];
  if (diffx>diffy) {
    result->point_max->coords[1]=result->point_min->coords[1]+diffx;
  } else {
    result->point_max->coords[0]=result->point_min->coords[0]+diffy;
  }

  //add a margin 
  result->margin=dist_points(result->point_min,result->point_max)/20;
  for(i=0;i<point_min->dim;i++) {
    result->point_min->coords[i]-=result->margin;
    result->point_max->coords[i]+=result->margin;
  }

  //set window size
  result->screen=SDL_SetVideoMode(w,h,32,SDL_DOUBLEBUF|SDL_HWSURFACE);
  result->w=w;
  result->h=h;

  //fill the background with black
  color=SDL_MapRGB(result->screen->format,255,255,255);
  rect.x=0;
  rect.y=0;
  rect.w=w;
  rect.h=h;
  SDL_FillRect(result->screen,&rect,color);

  //return the workspace
  return result;
}

void free_splot(struct splot *ws) {
  if (!ws)
    return;
  free_point(ws->point_min);
  free_point(ws->point_max);
  free(ws);
}

struct point *sdl2real_splot(int i,int j,struct splot *ws) {
  struct point *pt=new_point(2);
  double xmin=ws->point_min->coords[0];
  double ymin=ws->point_min->coords[1];
  double xmax=ws->point_max->coords[0];
  double ymax=ws->point_max->coords[1];
  pt->coords[0]=xmin+((double)(i-ws->margin)/(double)(ws->w-2*ws->margin))*(xmax-xmin);
  pt->coords[1]=ymin+((double)(ws->h-j-ws->margin)/(double)(ws->h-2*ws->margin))*(ymax-ymin);
  return pt;
}

struct point *real2sdl_splot(double x,double y,struct splot *ws) {
  struct point *pt=new_point(2);
  double xmin=ws->point_min->coords[0];
  double ymin=ws->point_min->coords[1];
  double xmax=ws->point_max->coords[0];
  double ymax=ws->point_max->coords[1];
  pt->coords[0]=(int)((x-xmin)/(xmax-xmin)*(ws->w-2*ws->margin))+ws->margin;
  pt->coords[1]=ws->h-(int)((y-ymin)/(ymax-ymin)*(ws->h-2*ws->margin)+ws->margin);
  //printf("sdl coords=%d,%d\n",(int)pt->coords[0],(int)pt->coords[1]);
  return pt;
}

void plot_real_point_splot(struct point *src,struct splot *ws,int r,int g,int b) {
  struct point *pt=real2sdl_splot(src->coords[0],src->coords[1],ws);
  Uint32 color=SDL_MapRGB(ws->screen->format,r,g,b);
  SDL_Rect rect;
  rect.x=(int)pt->coords[0]-1;
  rect.y=(int)pt->coords[1]-1;
  //printf("plotting point at %d,%d\n",rect.x,rect.y);
  free_point(pt);
  rect.w=3;
  rect.h=3;
  SDL_FillRect(ws->screen,&rect,color);
}

void plot_file_splot(struct splot *ws,char *filename,int r,int g,int b) {
  FILE *f;
  struct point *pt=new_point(2);
  f=fopen(filename,"r");
  while(!feof(f)) {
    if (fscanf(f,"%lf %lf",&(pt->coords[0]),&(pt->coords[1]))==2)
      plot_real_point_splot(pt,ws,r,g,b);
  }
  fclose(f);
  free_point(pt);
}

void apply_splot(struct splot *ws) {
  SDL_Flip(ws->screen);
}

void clear_splot(struct splot *ws) {
  Uint32 color;
  SDL_Rect rect;
  color=SDL_MapRGB(ws->screen->format,0,0,0);
  rect.x=0;
  rect.y=0;
  rect.w=ws->w;
  rect.h=ws->h;
  SDL_FillRect(ws->screen,&rect,color);
}

void set_pixel_splot(struct splot *ws,int x,int y,int r,int g,int b) {
  Uint32 color;
  SDL_Surface *surface=ws->screen;
  int bpp=surface->format->BytesPerPixel;
  Uint8 *p=(Uint8 *)surface->pixels+y*surface->pitch+x*bpp;
  color=SDL_MapRGB(surface->format,r,g,b);
  *(Uint32 *)p=color;
}

void get_pixel_splot(struct splot *ws,int x,int y,int *r,int *g,int *b) {
  Uint8 rp,gp,bp;
  SDL_Surface *surface=ws->screen;
  int bpp=surface->format->BytesPerPixel;
  Uint8 *p=(Uint8 *)surface->pixels+y*surface->pitch+x*bpp;
  SDL_GetRGB(*(Uint32 *)p,surface->format,&rp,&gp,&bp);
  *r=rp;
  *g=gp;
  *b=bp;
}

void attribute_color_splot(int id,int *r,int *v,int *b) {
  if (id%3==0) {
    *r=255;
    *v=0;
    *b=0;
  }
  if (id%3==1) {
    *r=0;
    *v=255;
    *b=0;
  }
  if (id%3==2) {
    *r=0;
    *v=0;
    *b=255;
  }
}



