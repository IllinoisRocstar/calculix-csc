/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

void writeHomogenizedVals(double *stn, double *een, ITG *nk, ITG *ne, double *co,
                          ITG *kon, ITG *ipkon, double tStep, double *elemvols) {

  // Open File
  FILE *fp;
  fp=fopen("homogenization.csv","a+");

  // Storing all stresses and strains element-wise
  double *E11=NULL;
  NNEW(E11,double,*ne);
  double *E22=NULL;
  NNEW(E22,double,*ne);
  double *E33=NULL;
  NNEW(E33,double,*ne);
  double *E23=NULL;
  NNEW(E23,double,*ne);
  double *E13=NULL;
  NNEW(E13,double,*ne);
  double *E12=NULL;
  NNEW(E12,double,*ne);

  double *S11=NULL;
  NNEW(S11,double,*ne);
  double *S22=NULL;
  NNEW(S22,double,*ne);
  double *S33=NULL;
  NNEW(S33,double,*ne);
  double *S23=NULL;
  NNEW(S23,double,*ne);
  double *S13=NULL;
  NNEW(S13,double,*ne);
  double *S12=NULL;
  NNEW(S12,double,*ne);

  int ne_ = (int) *ne;
  float sumVols = 0.0;

  for (int k=0; k<ne_; k++) {
    float elemE11 = 0;
    float elemE22 = 0;
    float elemE33 = 0;
    float elemE23 = 0;
    float elemE13 = 0;
    float elemE12 = 0;
    float elemS11 = 0;
    float elemS22 = 0;
    float elemS33 = 0;
    float elemS23 = 0;
    float elemS13 = 0;
    float elemS12 = 0;

    for (int j=0; j<4; j++) {
      elemE11 = elemE11 + een[0+6*(kon[ipkon[k]+j]-1)];
      elemE22 = elemE22 + een[1+6*(kon[ipkon[k]+j]-1)];
      elemE33 = elemE33 + een[2+6*(kon[ipkon[k]+j]-1)];
      elemE23 = elemE23 + een[3+6*(kon[ipkon[k]+j]-1)];
      elemE13 = elemE13 + een[4+6*(kon[ipkon[k]+j]-1)];
      elemE12 = elemE12 + een[5+6*(kon[ipkon[k]+j]-1)];
      elemS11 = elemS11 + stn[0+6*(kon[ipkon[k]+j]-1)];
      elemS22 = elemS22 + stn[1+6*(kon[ipkon[k]+j]-1)];
      elemS33 = elemS33 + stn[1+6*(kon[ipkon[k]+j]-1)];
      elemS23 = elemS23 + stn[1+6*(kon[ipkon[k]+j]-1)];
      elemS13 = elemS13 + stn[1+6*(kon[ipkon[k]+j]-1)];
      elemS12 = elemS12 + stn[1+6*(kon[ipkon[k]+j]-1)];
    }

    E11[k] = (elemE11/4)*elemvols[k];
    E22[k] = (elemE22/4)*elemvols[k];
    E33[k] = (elemE33/4)*elemvols[k];
    E23[k] = (elemE23/4)*elemvols[k];
    E13[k] = (elemE13/4)*elemvols[k];
    E12[k] = (elemE12/4)*elemvols[k];
    S11[k] = (elemS11/4)*elemvols[k];
    S22[k] = (elemS22/4)*elemvols[k];
    S33[k] = (elemS33/4)*elemvols[k];
    S23[k] = (elemS23/4)*elemvols[k];
    S13[k] = (elemS13/4)*elemvols[k];
    S12[k] = (elemS12/4)*elemvols[k];

    sumVols = sumVols + elemvols[k];
  }

  // Sum all stresses and strains followed by division with total volume of domain
  float totalE11 = 0.0;
  float totalE22 = 0.0;
  float totalE33 = 0.0;
  float totalE23 = 0.0;
  float totalE13 = 0.0;
  float totalE12 = 0.0;
  float totalS11 = 0.0;
  float totalS22 = 0.0;
  float totalS33 = 0.0;
  float totalS23 = 0.0;
  float totalS13 = 0.0;
  float totalS12 = 0.0;
  for (int i=0; i<ne_; i++) {
    totalE11 = totalE11 + E11[i];
    totalE22 = totalE22 + E22[i];
    totalE33 = totalE33 + E33[i];
    totalE23 = totalE23 + E23[i];
    totalE13 = totalE13 + E13[i];
    totalE12 = totalE12 + E12[i];
    totalS11 = totalS11 + S11[i];
    totalS22 = totalS22 + S22[i];
    totalS33 = totalS33 + S33[i];
    totalS23 = totalS23 + S23[i];
    totalS13 = totalS13 + S13[i];
    totalS12 = totalS12 + S12[i];
  }

  totalE11 = totalE11/(sumVols);
  totalE22 = totalE22/(sumVols);
  totalE33 = totalE33/(sumVols);
  totalE23 = totalE23/(sumVols);
  totalE13 = totalE13/(sumVols);
  totalE12 = totalE12/(sumVols);
  totalS11 = totalS11/(sumVols);
  totalS22 = totalS22/(sumVols);
  totalS33 = totalS33/(sumVols);
  totalS23 = totalS23/(sumVols);
  totalS13 = totalS13/(sumVols);
  totalS12 = totalS12/(sumVols);

  fprintf(fp, "%0.10f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f\n",
    tStep,totalE11,totalE22,totalE33,totalE23,totalE13,totalE12,totalS11,totalS22,totalS33,totalS23,totalS13,totalS12);
  fclose(fp);
  printf(" Writing homogenized stress and strain values....\n\n");
}