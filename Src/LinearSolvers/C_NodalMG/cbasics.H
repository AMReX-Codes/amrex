#ifndef _CBASICS_H_
#define _CBASICS_H_
/*
  Louis Howell, 4/19/90
  Version:  10/12/92
  Header file for using cbasics.c from C++
*/

extern "C" void gopen(int);
extern "C" void newdev(int);
extern "C" void gclose();
extern "C" void newpage();
extern "C" void linewidth(int);
extern "C" void black();
extern "C" void color(int,int,int);
extern "C" void getcolor(int*,int*,int*);
extern "C" void colorscale(Real);
extern "C" void dim(Real,Real,Real,Real);
extern "C" void fitbox(Real,Real,Real,Real,Real);
extern "C" void clipxy(Real,Real,Real,Real);
extern "C" void mkbox(Real,Real,Real,Real);
extern "C" void pm(Real,Real);
extern "C" void pc(Real,Real);
extern "C" void gflush(int);
extern "C" void getpt(Real*,Real*,int*);

#endif
