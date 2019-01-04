// Backwoard procedure for multiple change point detection
// Seung Jun Shin, Yichao Wu, and Ning Hao (2014+)
// "(Tentative) Backwoard procedure for multiple change point detection"
// The code here is copyright (c) 2014 Seung Jun Shin
// and may be used under the terms of the GNU Public Licsense
// version 2.0 or a later version which can be found at
// http://www.gnu.org/licenses/gpl.html
//  This code comes with no warranty.  Use at your own risk.

#include <R.h>
#include <Rmath.h>

void backward(
 double *m,
 int *z,
 double *rss,
 int *n,
 int *l,
 double *value,
 int *index,
 double *sigma,
 double *cutoff,
 int *ngroup,
 int *kmin
)
{
   int i, j, u, z1, z2;
   double m1, m2, v;
// double cutoffvalue = (*sigma) * (*sigma) * (*cutoff) * (*cutoff);
    for (i = 0; i <= (*n - *ngroup - 1); i++) {
      u = 0;
      v = rss[0];
      // find minimum
      for (j = 1; j <= (*l-2); j++)  {
                      if (v > rss[j]) {
            v = rss[j];
            u = j;
          }
      }
      value[i] = v;
      index[i] = u;

      z1 = z[u];
      z2 = z[u + 1];

      m1 = m[u];
      m2 = m[u + 1];

      //      if (v > cutoffvalue && m1 > *kmin && m2 > *kmin) return;

      // Update "m"ean vector
      m[u] = (m1 * z1 + m2 * z2) / (z1 + z2);
      if (u < (*l-2)) {
        for (j = (u+1); j <= (*l-2); j++) {
          m[j] = m[j+1];
          }
      }

      // update si"z"e vector
      z[u] = z1 + z2;

      if (u < (*l-2)) {
        for (j = (u+1); j <= (*l-2); j++) {
          z[j] = z[j+1];
          }
       }
      // update rss vector
      if (u > 0) { //left
         rss[u-1] = ((double)z[u-1] * (double)z[u])/(z[u-1] + z[u]) * (m[u-1] - m[u]) * (m[u-1] - m[u]);
      }

      if (u < (*l-2)) { // right
         rss[u] = ((double)z[u+1] * (double)z[u])/(z[u+1] + z[u]) * (m[u+1] - m[u]) * (m[u+1] - m[u]);
      }
      if (u < (*l-3))
      {
         for (j = (u+1); j <= (*l-3); j++) {
           rss[j] = rss[j+1];
         }
      }

     // update "l"ength
     *l = *l - 1;
     }
}
