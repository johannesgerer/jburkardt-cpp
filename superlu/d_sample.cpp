# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "slu_ddefs.h"

int main ( );
double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] );
void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  string title );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    D_SAMPLE tests the SUPERLU solver with a 5x5 double precision real matrix.
//
//  Discussion:
//
//    The general (GE) representation of the matrix is:
//
//      [ 19  0 21 21  0
//        12 21  0  0  0
//         0 12 16  0  0 
//         0  0  0  5 21
//        12 12  0  0 18 ]
//
//    The (0-based) compressed column (CC) representation of this matrix is:
//
//      I  CC   A
//     --  --  --
//      0   0  19
//      1      12
//      4      12
//
//      1   3  21
//      2      12
//      4      12
//
//      0   6  21
//      2      16
//
//      0   8  21
//      3       5
//
//      3  10  21
//      4      18
//
//      *  12   *
//
//    The right hand side B and solution X are
//
//      #   B     X
//     --  --  ----------
//      0   1  -0.03125
//      1   1   0.0654762
//      2   1   0.0133929
//      3   1   0.0625
//      4   1   0.0327381 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Demmel, John Gilbert, Xiaoye Li,
//    SuperLU Users's Guide.
//
{
# define N 5
# define NCC 12

  SuperMatrix A;
  double acc[NCC] = { 
    19.0, 12.0, 12.0, 
    21.0, 12.0, 12.0,
    21.0, 16.0, 
    21.0,  5.0, 
    21.0, 18.0 };
  double *b;
  double *b2;
  SuperMatrix B;
  int ccc[N+1] = { 0, 3, 6, 8, 10, 12 };
  int i;
  int icc[NCC] = { 
    0, 1, 4,
    1, 2, 4,
    0, 2, 
    0, 3, 
    3, 4 };
  int info;
  int j;
  SuperMatrix L;
  int m = N;
  int n = N;
  int nrhs = 1;
  int ncc = NCC;
  superlu_options_t options;
  int *perm_c;
  int permc_spec;
  int *perm_r;
  SuperLUStat_t stat;
  SuperMatrix U;

  timestamp ( );
  cout << "\n";
  cout << "D_SAMPLE:\n";
  cout << "  C++ version\n";
  cout << "  SUPERLU solves a double precision real linear system.\n";
//
//  Print the matrix.
//
  cc_print ( m, n, ncc, icc, ccc, acc, "  CC matrix:" );
//
//  Convert the compressed column (CC) matrix into a SuperMatrix A. 
//
  dCreate_CompCol_Matrix ( &A, m, n, ncc, acc, icc, ccc, SLU_NC, SLU_D, SLU_GE );    
//
//  Create the right-hand side matrix.
//
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    b[i] = 1.0;
  }
  cout << "\n";
  cout << "  Right hand side:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << b[i] << "\n";
  }
//
//  Create Super Right Hand Side.
//
  dCreate_Dense_Matrix ( &B, m, nrhs, b, m, SLU_DN, SLU_D, SLU_GE );
//
//  Set space for the permutations.
//
  perm_r = new int[m];
  perm_c = new int[n];
//
//  Set the input options. 
//
  set_default_options ( &options );
  options.ColPerm = NATURAL;
//
//  Initialize the statistics variables. 
//
  StatInit ( &stat );
//
//  Solve the linear system. 
//
  dgssv ( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info );
    
  dPrint_CompCol_Matrix ( ( char * ) "A", &A );
  dPrint_CompCol_Matrix ( ( char * ) "U", &U );
  dPrint_SuperNode_Matrix ( ( char * ) "L", &L );
  print_int_vec ( ( char * ) "\nperm_r", m, perm_r );
//
//  By some miracle involving addresses, 
//  the solution has been put into the B vector.
//
  cout << "\n";
  cout << "  Computed solution:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << b[i] << "\n";
  }
//
//  Demonstrate that RHS is really the solution now.
//  Multiply it by the matrix.
//
  b2 = cc_mv ( m, n, ncc, icc, ccc, acc, b );
  cout << "\n";
  cout << "  Product A*X:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << b2[i] << "\n";
  }
//
//  Free memory.
//
  free ( b );
  free ( b2 );
  free ( perm_c );
  free ( perm_r );

  Destroy_SuperMatrix_Store ( &A );
  Destroy_SuperMatrix_Store ( &B );
  Destroy_SuperNode_Matrix ( &L );
  Destroy_CompCol_Matrix ( &U );
  StatFree ( &stat );
//
//  Terminate.
//
  cout << "\n";
  cout << "D_SAMPLE:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;

# undef N
# undef NCC
}
//****************************************************************************80

double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CC_MV multiplies a CC matrix by a vector
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, int NCC, the number of CC values.
//
//    Input, int ICC[NCC], the CC rows.
//
//    Input, int CCC[N+1], the compressed CC columns
//
//    Input, double ACC[NCC], the CC values.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Output, double CC_MV[M], the product A*X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = ccc[j]; k < ccc[j+1]; k++ )
    {
      i = icc[k];
      b[i] = b[i] + acc[k] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    CC_PRINT prints a sparse matrix in CC format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in the matrix.
//
//    Input, int N, the number of columns in the matrix.
//
//    Input, int NCC, the number of CC elements.
//
//    Input, int ICC[NCC], the CC rows.
//
//    Input, int CCC[N+1], the compressed CC columns.
//
//    Input, double ACC[NCC], the CC values.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int j;
  int jnext;
  int k;

  cout << "\n";
  cout << title << "\n";
  cout << "     #     I     J       A\n";
  cout << "  ----  ----  ----  --------------\n";
  cout << "\n";

  j = 0;
  jnext = ccc[1];

  for ( k = 0; k < ncc; k++ )
  {
    i = icc[k];
    while ( jnext <= k )
    {
      j = j + 1;
      jnext = ccc[j+1];
    }
 
    cout << setw(4) << k << "  "
         << setw(4) << i << "  "
         << setw(4) << j << "  "
         << setw(16) << acc[k] << "\n";
  }

  return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
