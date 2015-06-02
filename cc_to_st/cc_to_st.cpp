# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "cc_to_st.hpp"

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
//    23 July 2014
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
  int k;

  cout << "\n";
  cout << title << "\n";
  cout << "     #     I     J         A\n";
  cout << "  ----  ----  ----  ----------------\n";
  cout << "\n";

  if ( ccc[0] == 0 )
  {
    j = 0;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j+1] <= k )
      {
        j = j + 1;
      }
      cout << setw(4) << k << "  "
           << setw(4) << i << "  "
           << setw(4) << j << "  "
           << setw(16) << acc[k] << "\n";
    }
  }
  else
  {
    j = 1;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j] <= k + 1 )
      {
        j = j + 1;
      }
      cout << setw(4) << k + 1 << "  "
           << setw(4) << i << "  "
           << setw(4) << j << "  "
          << setw(16) << acc[k] << "\n";
    }
  }

  return;
}
//****************************************************************************80

void cc_to_st ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  int &nst, int ist[], int jst[], double ast[] )

//****************************************************************************80
//
//  Purpose:
//
//    CC_TO_ST converts sparse matrix information from CC to ST format.
//
//  Discussion:
//
//    Only JST actually needs to be computed.  The other three output 
//    quantities are simply copies.  
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, int NCC, the number of CC elements.
//
//    Input, int ICC[NCC], the CC rows.
//
//    Input, int CCC[N+1], the CC compressed columns.
//
//    Input, double ACC[NCC], the CC values.
//
//    Output, int &NST, the number of ST elements.
//
//    Output, int IST[NST], JST[NST], the ST rows and columns.
//
//    Output, double AST[NST], the ST values.
//
{
  int j;
  int jhi;
  int jlo;
  int k;
  int khi;
  int klo;

  nst = 0;

  if ( ccc[0] == 0 )
  {
    jlo = 0;
    jhi = n - 1;
  
    for ( j = jlo; j <= jhi; j++ )
    {
      klo = ccc[j];
      khi = ccc[j+1] - 1;

      for ( k = klo; k <= khi; k++ )
      {
        ist[nst] = icc[k];
        jst[nst] = j;
        ast[nst] = acc[k];
        nst = nst + 1;
      }
    }
  }
  else
  {
    jlo = 1;
    jhi = n;
  
    for ( j = jlo; j <= jhi; j++ )
    {
      klo = ccc[j-1];
      khi = ccc[j] - 1;

      for ( k = klo; k <= khi; k++ )
      {
        ist[nst] = icc[k-1];
        jst[nst] = j;
        ast[nst] = acc[k-1];
        nst = nst + 1;
      }
    }
  }

  return;
}
//****************************************************************************80

void st_print ( int m, int n, int nst, int ist[], int jst[], double ast[], 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    ST_PRINT prints a sparse matrix in ST format.
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
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, int NST, the number of ST elements.
//
//    Input, int IST[NST], JST[NST], the ST rows and columns.
//
//    Input, double AST[NST], the ST values.
//
//    Input, string TITLE, a title.
//
{
  int k;

  cout << "\n";
  cout << title << "\n";
  cout << "     #     I     J       A\n";
  cout << "  ----  ----  ----  --------------\n";
  cout << "\n";
  for ( k = 0; k < nst; k++ )
  {
    cout << setw(4) << k << "  "
         << setw(4) << ist[k] << "  "
         << setw(4) << jst[k] << "  "
         << setw(16) << ast[k] << "\n";
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
