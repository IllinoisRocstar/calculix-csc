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
    if (lakon[8*i+3]=='4') {
      // Find volume for tet
      double A[3], B[3], C[3], D[3];

      if (ipkon[i] > -1) {
        A[0] = co[0+(kon[ipkon[i]+0]-1)*3];  // 0
        A[1] = co[1+(kon[ipkon[i]+0]-1)*3];  // 0
        A[2] = co[2+(kon[ipkon[i]+0]-1)*3];  // 0

        B[0] = co[0+(kon[ipkon[i]+1]-1)*3];  // 1
        B[1] = co[1+(kon[ipkon[i]+1]-1)*3];  // 1
        B[2] = co[2+(kon[ipkon[i]+1]-1)*3];  // 1

        C[0] = co[0+(kon[ipkon[i]+2]-1)*3];  // 2
        C[1] = co[1+(kon[ipkon[i]+2]-1)*3];  // 2
        C[2] = co[2+(kon[ipkon[i]+2]-1)*3];  // 2

        D[0] = co[0+(kon[ipkon[i]+3]-1)*3];  // 3
        D[1] = co[1+(kon[ipkon[i]+3]-1)*3];  // 3
        D[2] = co[2+(kon[ipkon[i]+3]-1)*3];  // 3
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

    } else if (lakon[8*i+3] == '8') {
      // Find volume for hex
      double A[3], B[3], C[3], D[3], E[3], F[3], G[3], H[3];

      // Get all the vertices
      if (ipkon[i] > -1) {
        A[0] = co[0+(kon[ipkon[i]+0]-1)*3];  // 0
        A[1] = co[1+(kon[ipkon[i]+0]-1)*3];  // 0
        A[2] = co[2+(kon[ipkon[i]+0]-1)*3];  // 0

        B[0] = co[0+(kon[ipkon[i]+1]-1)*3];  // 1
        B[1] = co[1+(kon[ipkon[i]+1]-1)*3];  // 1
        B[2] = co[2+(kon[ipkon[i]+1]-1)*3];  // 1

        C[0] = co[0+(kon[ipkon[i]+2]-1)*3];  // 2
        C[1] = co[1+(kon[ipkon[i]+2]-1)*3];  // 2
        C[2] = co[2+(kon[ipkon[i]+2]-1)*3];  // 2

        D[0] = co[0+(kon[ipkon[i]+3]-1)*3];  // 3
        D[1] = co[1+(kon[ipkon[i]+3]-1)*3];  // 3
        D[2] = co[2+(kon[ipkon[i]+3]-1)*3];  // 3

        E[0] = co[0+(kon[ipkon[i]+4]-1)*3];  // 4
        E[1] = co[1+(kon[ipkon[i]+4]-1)*3];  // 4
        E[2] = co[2+(kon[ipkon[i]+4]-1)*3];  // 4

        F[0] = co[0+(kon[ipkon[i]+5]-1)*3];  // 5
        F[1] = co[1+(kon[ipkon[i]+5]-1)*3];  // 5
        F[2] = co[2+(kon[ipkon[i]+5]-1)*3];  // 5

        G[0] = co[0+(kon[ipkon[i]+6]-1)*3];  // 6
        G[1] = co[1+(kon[ipkon[i]+6]-1)*3];  // 6
        G[2] = co[2+(kon[ipkon[i]+6]-1)*3];  // 6

        H[0] = co[0+(kon[ipkon[i]+7]-1)*3];  // 7
        H[1] = co[1+(kon[ipkon[i]+7]-1)*3];  // 7
        H[2] = co[2+(kon[ipkon[i]+7]-1)*3];  // 7
      }

      // Initialize
      double totalVol = 0.0;

      // Find volume for tet 1
      double a1[3], b1[3], c1[3];
      for (int j=0; j<3; j++) {
        a1[j] = (A[j] - B[j]);
        b1[j] = (B[j] - D[j]);
        c1[j] = (D[j] - E[j]);
      }

      totalVol = totalVol + fabs(((a1[0]*((b1[1]*c1[2])-(c1[1]*b1[2])))
                          - (a1[1]*((b1[0]*c1[2])-(c1[0]*b1[2])))
                          + (a1[2]*((b1[0]*c1[1])-(c1[0]*b1[1]))))/6.0);

      // Find volume for tet 2
      double a2[3], b2[3], c2[3];
      for (int j=0; j<3; j++) {
        a2[j] = (B[j] - C[j]);
        b2[j] = (C[j] - D[j]);
        c2[j] = (D[j] - G[j]);
      }

      totalVol = totalVol + fabs(((a2[0]*((b2[1]*c2[2])-(c2[1]*b2[2])))
                          - (a2[1]*((b2[0]*c2[2])-(c2[0]*b2[2])))
                          + (a2[2]*((b2[0]*c2[1])-(c2[0]*b2[1]))))/6.0);

      // Find volume for tet 3
      double a3[3], b3[3], c3[3];
      for (int j=0; j<3; j++) {
        a3[j] = (B[j] - E[j]);
        b3[j] = (E[j] - F[j]);
        c3[j] = (F[j] - G[j]);
      }

      totalVol = totalVol + fabs(((a3[0]*((b3[1]*c3[2])-(c3[1]*b3[2])))
                          - (a3[1]*((b3[0]*c3[2])-(c3[0]*b3[2])))
                          + (a3[2]*((b3[0]*c3[1])-(c3[0]*b3[1]))))/6.0);

      // Find volume for tet 4
      double a4[3], b4[3], c4[3];
      for (int j=0; j<3; j++) {
        a4[j] = (D[j] - E[j]);
        b4[j] = (E[j] - G[j]);
        c4[j] = (G[j] - H[j]);
      }

      totalVol = totalVol + fabs(((a4[0]*((b4[1]*c4[2])-(c4[1]*b4[2])))
                          - (a4[1]*((b4[0]*c4[2])-(c4[0]*b4[2])))
                          + (a4[2]*((b4[0]*c4[1])-(c4[0]*b4[1]))))/6.0);

      // Find volume for tet 5
      double a5[3], b5[3], c5[3];
      for (int j=0; j<3; j++) {
        a5[j] = (B[j] - D[j]);
        b5[j] = (D[j] - E[j]);
        c5[j] = (E[j] - G[j]);
      }

      totalVol = totalVol + fabs(((a5[0]*((b5[1]*c5[2])-(c5[1]*b5[2])))
                          - (a5[1]*((b5[0]*c5[2])-(c5[0]*b5[2])))
                          + (a5[2]*((b5[0]*c5[1])-(c5[0]*b5[1]))))/6.0);

      elemvols[i] = 1*totalVol;
    } else {
      printf("Element %i is not supported yet.\n",i);
      continue;
    }
  }
}
