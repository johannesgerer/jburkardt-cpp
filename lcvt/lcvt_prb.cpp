# include <cstdlib>
# include <iostream>
# include <iomanip>

# include "lcvt.hpp"

using namespace std;

int main ( );
void test01 ( int sample_function_cvt );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    LCVT_PRB calls a set of problems for LCVT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int sample_function_cvt;

  timestamp ( );

  cout << "\n";
  cout << "LCVT_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LCVT library.\n";
 
  for ( i = -1; i <= 2; i++ )
  {
    sample_function_cvt = i;
    test01 ( sample_function_cvt );
  }

  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LCVT_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int sample_function_cvt )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CVT, R8MAT_LATINIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define N 25

  double generator[M*N];
  int i;
  int latin_steps = 3;
  int sample_function_init = 0;
  int sample_num_cvt = 100000;
  int sample_num_steps = 50;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  R8MAT_LATINIZE makes it a Latin Hypersquare.\n";
  cout << "\n";
  cout << "  In this test, we vary the sampling used during the\n";
  cout << "  CVT Latin iteration.\n";
//
//  GET_SEED can be used to produce a different seed on each run.
//  But using a fixed seed is useful for debugging.
//
  seed = get_seed ( );

  seed = 123456789;

  cout << "\n";
  cout << "  Spatial dimension M =        " << M << "\n";
  cout << "  Number of generators =       " << N << "\n";
  cout << "  Initial random number seed = " << seed << "\n";;
  cout << "\n";

  if ( sample_function_init == -1 )
  {
    cout << "  Initialize using RANDOM_NUMBER (C++ STDLIB intrinsic).\n";
  }
  else if ( sample_function_init == 0 )
  {
    cout << "  Initialize using UNIFORM.\n";
  }
  else if ( sample_function_init == 1 )
  {
    cout << "  Initialize using HALTON.\n";
  }
  else if ( sample_function_init == 2 )
  {
    cout << "  Initialize using GRID.\n";
  }
  else if ( sample_function_init == 3 )
  {
    cout << "  USER will initialize data.\n";
  }

  if ( sample_function_cvt == -1 )
  {
    cout << "  Sample using RANDOM_NUMBER (C++ STDLIB intrinsic).\n";
  }
  else if ( sample_function_cvt == 0 )
  {
    cout << "  Sample using UNIFORM.\n";
  }
  else if ( sample_function_cvt == 1 )
  {
    cout << "  Sample using HALTON.\n";
  }
  else if ( sample_function_cvt == 2 )
  {
    cout << "  Sample using GRID.\n";
  }

  cout << "  Number of sample points = " << sample_num_cvt << "\n";
  cout << "  Number of sample steps =  " << sample_num_steps << "\n";

  for ( i = 1; i <= latin_steps; i++ )
  {
    cvt ( M, N, sample_function_init, sample_function_cvt,
      sample_num_cvt, sample_num_steps, &seed, generator );

    r8mat_transpose_print ( M, N, generator, "  After CVT steps:" );

    r8mat_latinize ( M, N, generator );

    r8mat_transpose_print ( M, N, generator, "  After Latin step:" );

  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CVT. R8MAT_LATINIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define N 25

  double generator[M*N];
  int i;
  int j;
  int latin_steps = 3;
  int ngrid;
  int rank;
  int sample_function_cvt = 0;
  int sample_function_init = 3;
  int sample_num_cvt = 100000;
  int sample_num_steps = 50;
  int seed;
  int tuple[M];
//
  cout << "\n";
  cout << "TEST02\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  R8MAT_LATINIZE makes it a Latin Hypersquare.\n";
  cout << "\n";
  cout << "  In this test, we initialize the generators to\n";
  cout << "  grid points; this is an unstable CVT solution.\n";
//
//  GET_SEED can be used to produce a different seed on each run.
//  But using a fixed seed is useful for debugging.
//
  seed = get_seed ( );
 
  seed = 123456789;

  cout << "\n";
  cout << "  Spatial dimension M =        " << M << "\n";
  cout << "  Number of generators =       " << N << "\n";
  cout << "  Initial random number seed = " << seed << "\n";;
  cout << "\n";

  if ( sample_function_init == -1 )
  {
    cout << "  Initialize using RANDOM_NUMBER (C++ STDLIB intrinsic).\n";
  }
  else if ( sample_function_init == 0 )
  {
    cout << "  Initialize using UNIFORM.\n";
  }
  else if ( sample_function_init == 1 )
  {
    cout << "  Initialize using HALTON.\n";
  }
  else if ( sample_function_init == 2 )
  {
    cout << "  Initialize using GRID.\n";
  }
  else if ( sample_function_init == 3 )
  {
    cout << "  USER will initialize data.\n";
  }

  if ( sample_function_cvt == -1 )
  {
    cout << "  Sample using RANDOM_NUMBER (C++ STDLIB intrinsic).\n";
  }
  else if ( sample_function_cvt == 0 )
  {
    cout << "  Sample using UNIFORM.\n";
  }
  else if ( sample_function_cvt == 1 )
  {
    cout << "  Sample using HALTON.\n";
  }
  else if ( sample_function_cvt == 2 )
  {
    cout << "  Sample using GRID.\n";
  }

  cout << "  Number of sample points = " << sample_num_cvt << "\n";
  cout << "  Number of sample steps =  " << sample_num_steps << "\n";

  ngrid = 5;

  for ( rank = 0; rank <= N-1; rank++ )
  {
    tuple_next_fast ( ngrid, M, rank, tuple );
    for ( i = 0; i < M; i++ )
    {
      generator[i+rank*M] = ( double ) ( 2 * tuple[i] - 1 ) 
                          / ( double ) ( 2 * ngrid        );
    }
  }

  r8mat_transpose_print ( M, N, generator, "  Initial generators (rows):" );

  for ( i = 1; i <=latin_steps; i++ )
  {
    cvt ( M, N, sample_function_init, sample_function_cvt, 
      sample_num_cvt, sample_num_steps, &seed, generator );

    r8mat_transpose_print ( M, N, generator, "  After CVT steps:" );

    r8mat_latinize ( M, N, generator );

    r8mat_transpose_print ( M, N, generator, "  After Latin step:" );

  }

  return;
# undef M
# undef N
}
