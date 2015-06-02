# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>

using namespace std;

# include "hammersley.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HAMMERSLEY_PRB.
//
//  Discussion:
//
//    HAMMERSLEY_PRB tests the HAMMERSLEY library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HAMMERSLEY_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the HAMMERSLEY library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HAMMERSLEY_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests I4_TO_HAMMERSLEY_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  int nmax = 1000;
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of \n";
  cout << "  a Hammersley sequence on a single call.\n";
  cout << "  All arguments are specified explicitly.\n";
  cout << "\n";
  cout << "  In this example, we compute the first 10 elements\n";
  cout << "  of a \"classical\" Hammersley sequence, and then\n";
  cout << "  the \"last\" 10 elements.\n";

  step = 1;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  base[0] = -nmax;
  for ( i = 1; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i );
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  We can jump ahead in the sequence by changing STEP:\n";

  step = nmax - N + 1;

  cout << "\n";
  cout << "  STEP = " << step << "\n";

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                        << "  "
         << setw(6) << step+j-1 << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests I4_TO_HAMMERSLEY_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 12
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of \n";
  cout << "  a Hammersley sequence on a single call.\n";
  cout << "  All arguments are specified explicitly.\n";
  cout << "\n";
  cout << "  We are free to choose the values of BASE.\n";
  cout << "  Any negative value indicates a sequence of\n";
  cout << "  J/(-BASE) in that coordinate.\n";
  cout << "\n";
  cout << "  In this example, that is the only kind of base we use.\n";

  step = 0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = - ( int ) pow ( ( double ) 10, i+1 );
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests I4_TO_HAMMERSLEY_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 12
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of \n";
  cout << "  a Hammersley sequence on a single call.\n";
  cout << "  All arguments are specified explicitly.\n";
  cout << "\n";
  cout << "  The SEED vector allows us to define the zeroth\n";
  cout << "  element of the coordinate subsequence.\n";
  cout << "  That is, if we ask for the STEP=0 entry of the\n";
  cout << "  subsequence, we will get the SEED(I)th entry\n";
  cout << "  of the full sequence.\n";
  cout << "\n";
  cout << "  In this example, we use a fixed base for simplicity.\n";

  step = 0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 10 * ( i + 1 );
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = -100;
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests I4_TO_HAMMERSLEY_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 12
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of \n";
  cout << "  a Hammersley sequence on a single call.\n";
  cout << "  All arguments are specified explicitly.\n";
  cout << "\n";
  cout << "  The LEAP vector allows us to define the distance\n";
  cout << "  (in the original sequence) between successive\n";
  cout << "  subsequence elements.\n";
  cout << "\n";
  cout << "  A LEAP of 1 means that, once we start sampling\n";
  cout << "  the sequence, we are taking every element.\n";
  cout << "  A LEAP of 2 takes every other sequence element,\n";
  cout << "  and so on.\n";
  cout << "\n";
  cout << "  In this example, we use a fixed base for simplicity.\n";

  step = 0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = ( int ) pow ( ( double ) 2, i );
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = -100;
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests I4_TO_HAMMERSLEY_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int k;
  int leap[DIM_NUM];
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of \n";
  cout << "  a Hammersley sequence on a single call.\n";
  cout << "  All arguments are specified explicitly.\n";
  cout << "\n";
  cout << "  Any entry of the Hammersley sequence can be computed\n";
  cout << "  immediately, without having to compute the previous\n";
  cout << "  entries.  This is also true of the entries of the\n";
  cout << "  leaped Hammersley subsequences we generate.\n";
  cout << "\n";
  cout << "  The value of a component of the Hammersley sequence\n";
  cout << "  is computed directly from its index.  But there\n";
  cout << "  should not be much difficulty handling indices\n";
  cout << "  that go as high as a million or a billion.\n";
  cout << "\n";
  cout << "  In this example, we look at high index entries,\n";
  cout << "  attained by large values of STEP, or SEED or LEAP.\n";
  cout << "\n";
  cout << "  In this example, we use the default bases.\n";
  cout << "\n";
  cout << "\n";
  cout << "  BIG VALUES OF STEP:\n";
  cout << "\n";

  for ( k = 1; k <= 4; k++ )
  {
    step = ( int ) pow ( ( double ) 100, k );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      seed[i] = 0;
    }
    for ( i = 0; i < DIM_NUM; i++ )
    {
      leap[i] = 1;
    }
    base[0] = -(step+N-1);
    for ( i = 1; i < DIM_NUM; i++ )
    {
      base[i] = prime ( i );
    }

    cout << "\n";
    cout << "  DIM_NUM = " << DIM_NUM << "\n";
    cout << "  N =    " << N    << "\n";
    cout << "  STEP = " << step << "\n";
    i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
    i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
    i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

    i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

    cout << "\n";
    cout << "    STEP   Hammersley\n";
    cout << "\n";
    for ( j = 0; j < N; j++ )
    {
      cout                      << "  "
           << setw(6) << step+j << "  ";
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << setw(12) << r[i+j*DIM_NUM] << "  ";
      }
      cout << "\n";
    }
  }

  cout << "\n";
  cout << "\n";
  cout << "  BIG VALUES OF SEED:\n";
  cout << "\n";

  step = 0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = pow ( ( double ) 100, i+1 );
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  base[0] = -(100+N-1);
  for ( i = 1; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i );
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "\n";
  cout << "  BIG VALUES OF LEAP:\n";
  cout << "\n";

  step = 0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = ( int ) pow ( ( double ) 100, i+1 );
  }
  base[0] = - ( N * 100 );
  for ( i = 1; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i );
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests HAMMERSLEY_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  int nmax = 1000;
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  HAMMERSLEY_SEQUENCE computes N elements of\n";
  cout << "  a Hammersley sequence on a single call.\n";
  cout << "  All arguments are specified externally, by calling\n";
  cout << "  various setup routines.\n";
  cout << "\n";
  cout << "  In this example, we compute the first 10 elements\n";
  cout << "  of a \"classical\" Hammersley sequence, and then\n";
  cout << "  the \"last\" 10 elements.\n";

  step = 1;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  base[0] = -nmax;
  for ( i = 1; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i );
  }

  hammersley_dim_num_set ( DIM_NUM );
  hammersley_step_set ( step );
  hammersley_seed_set ( seed );
  hammersley_leap_set ( leap );
  hammersley_base_set ( base );

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  hammersley_sequence ( N, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  We can jump ahead in the sequence by changing STEP:\n";

  step = nmax - N + 1;
  hammersley_step_set ( step );

  cout << "\n";
  cout << "  STEP = " << step << "\n";

  hammersley_sequence ( N, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests HALHAM_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 4

  int base[DIM_NUM];
  char file_name[80] = "hammersley_04_00010.txt";
  int i;
  int j;
  int leap[DIM_NUM];
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  HALHAM_WRITE writes a Halton or Hammersley dataset to a file.\n";

  step = 0;
  seed[0] = 1;
  for ( i = 1; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  base[0] = -N;
  for ( i = 1; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i );
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  halham_write ( DIM_NUM, N, step, seed, leap, base, r, file_name );

  cout << "\n";
  cout << "  The data was written to \"" << file_name << "\".\n";

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests I4_TO_HAMMERSLEY_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 20
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  int nmax = 10;
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of \n";
  cout << "  a Hammersley sequence on a single call.\n";
  cout << "  All arguments are specified explicitly.\n";
  cout << "\n";
  cout << "  In this example, we demonstrate that any coordinate of\n";
  cout << "  the generalized Hammersley sequence that is generated\n";
  cout << "  as a fractional sequence J/|BASE(I)| will\n";
  cout << "  \"wrap around\".\n";

  step = 1;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  base[0] = -nmax;
  for ( i = 1; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i );
  }

  cout << "\n";
  cout << "  DIM_NUM = " << DIM_NUM << "\n";
  cout << "  N =    " << N    << "\n";
  cout << "  STEP = " << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_hammersley_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Hammersley\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                      << "  "
         << setw(6) << step+j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
