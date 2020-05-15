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

void calcElmVols(ITG *nk, ITG *ne, double *co, ITG *kon, ITG *ipkon, char *lakon,
                 double *elemvols) {

  // Calculate 4-node tetrahedral element volumes
  for (int i=0; i<*ne; i++) {
    double A[3], B[3], C[3], D[3];

    // checking element type
    if (lakon[8*i+3]!='4')
    {
        printf("Element %i is not supported yet.\n",i);
        continue;
    }

    if (ipkon[i] > -1) {
      A[0] = co[0+(kon[ipkon[i]+0]-1)*3];
      A[1] = co[1+(kon[ipkon[i]+0]-1)*3];
      A[2] = co[2+(kon[ipkon[i]+0]-1)*3];

      B[0] = co[0+(kon[ipkon[i]+1]-1)*3];
      B[1] = co[1+(kon[ipkon[i]+1]-1)*3];
      B[2] = co[2+(kon[ipkon[i]+1]-1)*3];

      C[0] = co[0+(kon[ipkon[i]+2]-1)*3];
      C[1] = co[1+(kon[ipkon[i]+2]-1)*3];
      C[2] = co[2+(kon[ipkon[i]+2]-1)*3];

      D[0] = co[0+(kon[ipkon[i]+3]-1)*3];
      D[1] = co[1+(kon[ipkon[i]+3]-1)*3];
      D[2] = co[2+(kon[ipkon[i]+3]-1)*3];
    }

    // Calculating volume now
    double a[3], b[3], c[3];
    for (int j=0; j<3; j++) {
      a[j] = (A[j] - B[j]);
      b[j] = (B[j] - C[j]);
      c[j] = (C[j] - D[j]);
    }

    elemvols[i] = fabs(((a[0]*((b[1]*c[2])-(c[1]*b[2])))
                  - (a[1]*((b[0]*c[2])-(c[0]*b[2])))
                  + (a[2]*((b[0]*c[1])-(c[0]*b[1]))))/6.0);
  }
}
