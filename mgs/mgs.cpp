# include <cstdlib>
# include <cmath>

# include "mgs.hpp"

void mgs ( int m, int n, float **c, float **r, float **q )
{
  int i;
  int j;
  int k;
  float *x;
  float xn;

  x = new float[m];
  for ( k = 0; k < n; k++ )
  {
  for ( j = 0; j < m; j++ )
  {
  x[j] = c[j][k];
  }
  xn = 0;
  for ( j = 0; j < m; j++ )
  {
  xn = xn + x[j]*x[j];
  }
  r[k][k] = sqrt(xn);
  if ( 0.0 < r[k][k] )
  {
  for ( j = 0; j < m; j++ )
  {
  q[j][k] = c[j][k]/r[k][k];
  }
  }
  else
  {
  for ( j = 0; j < m; j++ )
  {
  q[j][k] = 0.0;
  }
  }
  for ( j = k + 1; j < n; j++ )
  {
  r[k][j] = 0;
  for ( i = 0; i < m; i++ )
  {
  r[k][j] = r[k][j] + q[i][k]*c[i][j];
  }
  for ( i = 0; i < m; i++ )
  {
  c[i][j] = c[i][j] - q[i][k]*r[k][j];
  }
  }
  }
  delete [] x;
  return;
}
