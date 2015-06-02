# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "monomial.hpp"

int main ( );
void mono_between_enum_test ( );
void mono_between_next_grevlex_test ( );
void mono_between_next_grlex_test ( );
void mono_between_random_test ( );
void mono_next_grevlex_test ( );
void mono_next_grlex_test ( );
void mono_print_test ( );
void mono_rank_grlex_test ( );
void mono_total_enum_test ( );
void mono_total_next_grevlex_test ( );
void mono_total_next_grlex_test ( );
void mono_total_random_test ( );
void mono_unrank_grlex_test ( );
void mono_upto_enum_test ( );
void mono_upto_next_grlex_test ( );
void mono_upto_next_grevlex_test ( );
void mono_upto_random_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MONOMIAL_PRB.
//
//  Discussion:
//
//    MONOMIAL_PRB tests the MONOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "MONOMIAL_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the MONOMIAL library.\n";

  mono_between_enum_test ( );
  mono_between_next_grevlex_test ( );
  mono_between_next_grlex_test ( );
  mono_between_random_test ( );
  mono_next_grevlex_test ( );
  mono_next_grlex_test ( );
  mono_print_test ( );
  mono_rank_grlex_test ( );
  mono_total_enum_test ( );
  mono_total_next_grevlex_test ( );
  mono_total_next_grlex_test ( );
  mono_total_random_test ( );
  mono_unrank_grlex_test ( );
  mono_upto_enum_test ( );
  mono_upto_next_grevlex_test ( );
  mono_upto_next_grlex_test ( ); 
  mono_upto_random_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "MONOMIAL_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void mono_between_enum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_ENUM_TEST tests MONO_BETWEEN_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n1;
  int n2;
  int v;

  cout << "\n";
  cout << "MONO_BETWEEN_ENUM_TEST\n";
  cout << "  MONO_BETWEEN_ENUM can enumerate the number of monomials\n";
  cout << "  in M variables, of total degree between N1 and N2.\n";

  m = 3;
  cout << "\n";
  cout << "  Using spatial dimension m = " << m << "\n";
  cout << "\n";
  cout << "   N2:";
  for ( n2 = 0; n2 <= 8; n2++ )
  {
    cout << "  " << setw(4) << n2;
  }
  cout << "\n";
  cout << "  N1 +------------------------------------------------------\n";
  for ( n1 = 0; n1 <= 8; n1++ )
  {
    cout << "  " << setw(2) << n1 << " |";
    for ( n2 = 0; n2 <= 8; n2++ )
    {
      v = mono_between_enum ( m, n1, n2 );
      cout << "  " << setw(4) << v;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void mono_between_next_grevlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_NEXT_GREVLEX_TEST tests MONO_BETWEEN_NEXT_GREVLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n1;
  int n2;
  int x[3];

  cout << "\n";
  cout << "MONO_BETWEEN_NEXT_GREVLEX_TEST\n";
  cout << "  MONO_BETWEEN_NEXT_GREVLEX can list the monomials\n";
  cout << "  in M variables, of total degree N between N1 and N2,\n";
  cout << "  in grevlex order, one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,N1).\n";
  cout << "  The process ends with (N2,0,...,0,0)\n";

  n1 = 2;
  n2 = 3;

  cout << "\n";
  cout << "  Let M =  " << m << "\n";
  cout << "      N1 = " << n1 << "\n";
  cout << "      N2 = " << n2 << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = n1; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n2 )
    {
      break;
    }
 
    mono_between_next_grevlex ( m, n1, n2, x );
    i = i + 1;
  }

  return;
}
//****************************************************************************80

void mono_between_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_NEXT_GRLEX_TEST tests MONO_BETWEEN_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n1;
  int n2;
  int x[3];

  cout << "\n";
  cout << "MONO_BETWEEN_NEXT_GRLEX_TEST\n";
  cout << "  MONO_BETWEEN_NEXT_GRLEX can list the monomials\n";
  cout << "  in M variables, of total degree N between N1 and N2,\n";
  cout << "  in grlex order, one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,N1).\n";
  cout << "  The process ends with (N2,0,...,0,0)\n";

  n1 = 2;
  n2 = 3;

  cout << "\n";
  cout << "  Let M =  " << m << "\n";
  cout << "      N1 = " << n1 << "\n";
  cout << "      N2 = " << n2 << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = n1; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n2 )
    {
      break;
    }
 
    mono_between_next_grlex ( m, n1, n2, x );
    i = i + 1;
  }

  return;
}
//****************************************************************************80

void mono_between_random_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_RANDOM_TEST tests MONO_BETWEEN_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int j;
  int n1;
  int n2;
  int rank;
  int seed;
  int test;
  int test_num;
  int *x;

  cout << "\n";
  cout << "MONO_BETWEEN_RANDOM_TEST\n";
  cout << "  MONO_BETWEEN_RANDOM selects at random a monomial\n";
  cout << "  in M dimensions of total degree between N1 and N2.\n";

  n1 = 2;
  n2 = 3;

  cout << "\n";
  cout << "  Let M =  " << m << "\n";
  cout << "      N1 = " << n1 << "\n";
  cout << "      N2 = " << n2 << "\n";
  cout << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    x = mono_between_random ( m, n1, n2, seed, rank );
    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
    delete [] x;
  }

  return;
}
//****************************************************************************80

void mono_next_grevlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_NEXT_GREVLEX_TEST tests MONO_NEXT_GREVLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int k;
  int m = 4;
  int *x;

  cout << "\n";
  cout << "MONO_NEXT_GREVLEX_TEST\n";
  cout << "  MONO_NEXT_GREVLEX computes the next monomial\n";
  cout << "  in M variables, in grevlex order.\n";
  cout << "\n";
  cout << "  Let M =  " << m << "\n";

  k = 0;
  x = new int[m];
  for ( i = 0; i < m; i++ )
  {
    x[i] = 0;
  }

  for ( ;; )
  {
    d = i4vec_sum ( m, x );
    cout << "  " << setw(2) << k
         << "  " << setw(2) << d
         << "  |";
    for ( i = 0; i < m; i++ )
    {
      cout << setw(2) << x[i];
    }
    cout << "\n";
    if ( x[0] == 3 )
    {
      break;
    }
    k = k + 1;
    mono_next_grevlex ( m, x );
  }

  delete [] x;

  return;
}
//****************************************************************************80

void mono_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_NEXT_GRLEX_TEST tests MONO_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int k;
  int m = 4;
  int *x;

  cout << "\n";
  cout << "MONO_NEXT_GRLEX_TEST\n";
  cout << "  MONO_NEXT_GRLEX computes the next monomial\n";
  cout << "  in M variables, in grlex order.\n";
  cout << "\n";
  cout << "  Let M =  " << m << "\n";

  k = 0;
  x = new int[m];
  for ( i = 0; i < m; i++ )
  {
    x[i] = 0;
  }

  for ( ;; )
  {
    d = i4vec_sum ( m, x );
    cout << "  " << setw(2) << k
         << "  " << setw(2) << d
         << "  |";
    for ( i = 0; i < m; i++ )
    {
      cout << setw(2) << x[i];
    }
    cout << "\n";
    if ( x[0] == 3 )
    {
      break;
    }
    k = k + 1;
    mono_next_grlex ( m, x );
  }

  delete [] x;

  return;
}
//****************************************************************************80

void mono_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_PRINT_TEST tests MONO_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  int f1[1] = { 5 };
  int f2[1] = { -5 };
  int f3[4] = { 2, 1, 0, 3 };
  int f4[3] = { 17, -3, 199 };
  int m;

  cout << "\n";
  cout << "MONO_PRINT_TEST\n";
  cout << "  MONO_PRINT can print out a monomial.\n";
  cout << "\n";

  m = 1;
  mono_print ( m, f1, "  Monomial [5]:" );

  m = 1;
  mono_print ( m, f2, "  Monomial [5]:" );

  m = 4;
  mono_print ( m, f3, "  Monomial [2,1,0,3]:" );

  m = 3;
  mono_print ( m, f4, "  Monomial [17,-3,199]:" );

  return;
}
//****************************************************************************80

void mono_rank_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_RANK_GRLEX_TEST tests MONO_RANK_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int rank;
  int test;
  int test_num = 8;
  int x[3];
  int x_test[3*8] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 0, 1, 
    0, 2, 0, 
    1, 0, 2, 
    0, 3, 1, 
    3, 2, 1, 
    5, 2, 1 };

  cout << "\n";
  cout << "MONO_RANK_GRLEX_TEST\n";
  cout << "  MONO_RANK_GRLEX returns the rank of a monomial in the sequence\n";
  cout << "  of all monomials in M dimensions, in grlex order.\n";

  cout << "\n";
  cout << "  Print a monomial sequence with ranks assigned.\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = 0; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(3) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  cout << "\n";
  cout << "  Now, given a monomial, retrieve its rank in the sequence:\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    for ( j = 0; j < m; j++ )
    {
      x[j] = x_test[j+test*m];
    }
    rank = mono_rank_grlex ( m, x );

    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void mono_total_enum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_ENUM_TEST tests MONO_TOTAL_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int v;

  cout << "\n";
  cout << "MONO_TOTAL_ENUM_TEST\n";
  cout << "  MONO_TOTAL_ENUM can enumerate the number of monomials\n";
  cout << "  in M variables, of total degree N.\n";

  cout << "\n";
  cout << "    N:";
  for ( n = 0; n <= 8; n++ )
  {
    cout << "  " << setw(4) << n;
  }
  cout << "\n";
  cout << "   m +------------------------------------------------------\n";
  for ( m = 1; m <= 8; m++ )
  {
    cout << "  " << setw(2) << m << " |";
    for ( n = 0; n <= 8; n++ )
    {
      v = mono_total_enum ( m, n );
      cout << "  " << setw(4) << v;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void mono_total_next_grevlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_NEXT_GREVLEX_TEST tests MONO_TOTAL_NEXT_GREVLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int x[3];

  cout << "\n";
  cout << "MONO_TOTAL_NEXT_GREVLEX_TEST\n";
  cout << "  MONO_TOTAL_NEXT_GREVLEX can list the monomials\n";
  cout << "  in M variables, of total degree N,\n";
  cout << "  in grevlex order, one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,N).\n";
  cout << "  The process ends with (N,0,...,0,0)\n";

  n = 3;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = n; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_total_next_grevlex ( m, n, x );
    i = i + 1;
  }

  return;
}
//****************************************************************************80

void mono_total_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_NEXT_GRLEX_TEST tests MONO_TOTAL_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int x[3];

  cout << "\n";
  cout << "MONO_TOTAL_NEXT_GRLEX_TEST\n";
  cout << "  MONO_TOTAL_NEXT_GRLEX can list the monomials\n";
  cout << "  in M variables, of total degree N,\n";
  cout << "  in grlex order, one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,N).\n";
  cout << "  The process ends with (N,0,...,0,0)\n";

  n = 3;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = n; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_total_next_grlex ( m, n, x );
    i = i + 1;
  }

  return;
}
//****************************************************************************80

void mono_total_random_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_RANDOM_TEST tests MONO_TOTAL_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int j;
  int n;
  int rank;
  int seed;
  int test;
  int test_num;
  int *x;

  cout << "\n";
  cout << "MONO_TOTAL_RANDOM_TEST\n";
  cout << "  MONO_TOTAL_RANDOM selects at random a monomial\n";
  cout << "  in M dimensions of total degree N.\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    x = mono_total_random ( m, n, seed, rank );
    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
    delete [] x;
  }

  return;
}
//****************************************************************************80

void mono_unrank_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UNRANK_GRLEX_TEST tests MONO_UNRANK_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int rank;
  int rank_max;
  int seed;
  int test;
  int test_num;
  int *x;

  cout << "\n";
  cout << "MONO_UNRANK_GRLEX_TEST\n";
  cout << "  MONO_UNRANK_GRLEX is given a rank, and returns the corresponding\n";
  cout << "  monomial in the sequence of all monomials in M dimensions\n";
  cout << "  in grlex order.\n";

  cout << "\n";
  cout << "  For reference, print a monomial sequence with ranks.\n";

  n = 4;
  rank_max = mono_upto_enum ( m, n );

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x = new int[n];
  for ( i = 0; i < m; i++ )
  {
    x[i] = 0;
  }

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(3) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  cout << "\n";
  cout << "  Now choose random ranks between 1 and " << rank_max << "\n";
  cout << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    rank = i4_uniform_ab ( 1, rank_max, seed );
    x = mono_unrank_grlex ( m, rank );
    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
    free ( x );
  }

  return;
}
//****************************************************************************80

void mono_upto_enum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_ENUM_TEST tests MONO_UPTO_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int v;

  cout << "\n";
  cout << "MONO_UPTO_ENUM_TEST\n";
  cout << "  MONO_UPTO_ENUM can enumerate the number of monomials\n";
  cout << "  in M variables, of total degree 0 up to N.\n";

  cout << "\n";
  cout << "    N:\n";
  for ( n = 0; n <= 8; n++ )
  {
    cout << "  " << setw(4) << n;
  }
  cout << "\n";
  cout << "   m +------------------------------------------------------\n";
  for ( m = 1; m <= 8; m++ )
  {
    cout << "  " << setw(2) << m << "  |";
    for ( n = 0; n <= 8; n++ )
    {
      v = mono_upto_enum ( m, n );
      cout << " " << setw(5) << v;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void mono_upto_next_grevlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_NEXT_GREVLEX_TEST tests MONO_UPTO_NEXT_GREVLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int x[3];

  cout << "\n";
  cout << "MONO_UPTO_NEXT_GREVLEX_TEST\n";
  cout << "  MONO_UPTO_NEXT_GREVLEX can list the monomials\n";
  cout << "  in M variables, of total degree up to N,\n";
  cout << "  in grevlex order, one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,0).\n";
  cout << "  The process ends with (N,0,...,0,0)\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = 0; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";


    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grevlex ( m, n, x );
    i = i + 1;
  }

  return;
}
//****************************************************************************80

void mono_upto_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_NEXT_GRLEX_TEST tests MONO_UPTO_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int i;
  int j;
  int n;
  int x[3];

  cout << "\n";
  cout << "MONO_UPTO_NEXT_GRLEX_TEST\n";
  cout << "  MONO_UPTO_NEXT_GRLEX can list the monomials\n";
  cout << "  in M variables, of total degree up to N,\n";
  cout << "  in grlex order, one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,0).\n";
  cout << "  The process ends with (N,0,...,0,0)\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = 0; 

  i = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << i << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  return;
}
//****************************************************************************80

void mono_upto_random_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_RANDOM_TEST tests MONO_UPTO_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int j;
  int n;
  int rank;
  int seed;
  int test;
  int test_num;
  int *x;

  cout << "\n";
  cout << "MONO_UPTO_RANDOM_TEST\n";
  cout << "  MONO_UPTO_RANDOM selects at random a monomial\n";
  cout << "  in M dimensions of total degree no greater than N.\n";

  n = 4;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";
  cout << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    x = mono_upto_random ( m, n, seed, rank );
    cout << "  " << setw(3) << rank << "    ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(2) << x[j];
    }
    cout << "\n";
    delete [] x;
  }

  return;
}
//****************************************************************************80

void mono_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_VALUE_TEST tests MONO_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int *f;
  int j;
  int m = 3;
  int n;
  int nx = 2;
  int rank;
  int seed;
  int test;
  int test_num;
  double *v;
  double x[3*2] = {
     1.0, 2.0, 3.0, 
    -2.0, 4.0, 1.0 };

  cout << "\n";
  cout << "MONO_VALUE_TEST\n";
  cout << "  MONO_VALUE evaluates a monomial.\n";

  n = 6;

  cout << "\n";
  cout << "  Let M = " << m << "\n";
  cout << "      N = " << n << "\n";

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    f = mono_upto_random ( m, n, seed, rank );
    cout << "\n";
    mono_print ( m, f, "  M(X) = " );
    v = mono_value ( m, nx, f, x );
    for ( j = 0; j < nx; j++ )
    {
      cout << "  M(" << x[0+j*m] 
           << "," << x[1+j*m]
           << "," << x[2+j*m]
           << ") = " << v[j] << "\n";
    }
    delete [] f;
    delete [] v;
  }

  return;
}

