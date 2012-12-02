# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_tri_int.H"

//****************************************************************************80

int get_problem_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_PROBLEM_NUM returns the number of test integration problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int PROBLEM_NUM, the number of test integration problems.
//
{
  int problem_num;

  problem_num = 22;

  return problem_num;
}
//****************************************************************************80

unsigned long get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, unsigned long GET_SEED, a random seed value.
//
{
# define UNSIGNED_LONG_MAX 4294967295UL

  time_t clock;
  int i;
  int hours;
  int minutes;
  int seconds;
  struct tm *lt;
  static unsigned long seed = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed == 0 )
  {
    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Extract HOURS.
//
    hours = lt->tm_hour;
//
//  In case of 24 hour clocks, shift so that HOURS is between 1 and 12.
//
    if ( 12 < hours )
    {
      hours = hours - 12;
    }
//
//  Move HOURS to 0, 1, ..., 11
//
    hours = hours - 1;

    minutes = lt->tm_min;

    seconds = lt->tm_sec;

    seed = seconds + 60 * ( minutes + 60 * hours );
//
//  We want values in [1,43200], not [0,43199].
//
    seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,UNSIGNED_LONG_MAX].
//
    seed = ( unsigned long ) 
      ( ( ( double ) seed )
      * ( ( double ) UNSIGNED_LONG_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
  }
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;

# undef UNSIGNED_LONG_MAX
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" I4.
//
{
  return 2147483647;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

double *p00_fun ( int problem, int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_FUN evaluates the integrand for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the number of the desired test problem.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2,N], the evaluation points.
//
//    Output, double P00_FUN[N], the function values.
//
{
  double *value;

  if ( problem == 1 )
  {
    value = p01_fun ( n, p );
  }
  else if ( problem == 2 )
  {
    value = p02_fun ( n, p );
  }
  else if ( problem == 3 )
  {
    value = p03_fun ( n, p );
  }
  else if ( problem == 4 )
  {
    value = p04_fun ( n, p );
  }
  else if ( problem == 5 )
  {
    value = p05_fun ( n, p );
  }
  else if ( problem == 6 )
  {
    value = p06_fun ( n, p );
  }
  else if ( problem == 7 )
  {
    value = p07_fun ( n, p );
  }
  else if ( problem == 8 )
  {
    value = p08_fun ( n, p );
  }
  else if ( problem == 9 )
  {
    value = p09_fun ( n, p );
  }
  else if ( problem == 10 )
  {
    value = p10_fun ( n, p );
  }
  else if ( problem == 11 )
  {
    value = p11_fun ( n, p );
  }
  else if ( problem == 12 )
  {
    value = p12_fun ( n, p );
  }
  else if ( problem == 13 )
  {
    value = p13_fun ( n, p );
  }
  else if ( problem == 14 )
  {
    value = p14_fun ( n, p );
  }
  else if ( problem == 15 )
  {
    value = p15_fun ( n, p );
  }
  else if ( problem == 16 )
  {
    value = p16_fun ( n, p );
  }
  else if ( problem == 17 )
  {
    value = p17_fun ( n, p );
  }
  else if ( problem == 18 )
  {
    value = p18_fun ( n, p );
  }
  else if ( problem == 19 )
  {
    value = p19_fun ( n, p );
  }
  else if ( problem == 20 )
  {
    value = p20_fun ( n, p );
  }
  else if ( problem == 21 )
  {
    value = p21_fun ( n, p );
  }
  else if ( problem == 22 )
  {
    value = p22_fun ( n, p );
  }
  else
  {
    cout << "\n";
    cout << "P00_FUN - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

double p00_monte_carlo ( int problem, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    P00_MONTE_CARLO applies the Monte Carlo rule to integrate a function.
//
//  Discussion:
//
//    The function f(x,y) is to be integrated over a triangle T.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Input, int N, the number of sample points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P00_MONTE_CARLO, the approximate integral.
//
{
  double area;
  double *f;
  double *p;
  double result;
  double *t;

  t = p00_vertices ( problem );

  p = triangle_sample ( t, n, seed );

  f = p00_fun ( problem, n, p );

  area = triangle_area ( t );

  result = area * r8vec_sum ( n, f ) / ( double ) ( n );

  delete [] f;
  delete [] p;
  delete [] t;

  return result;
}
//****************************************************************************80

int p00_singularity ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_SINGULARITY warns of common singularities for any problem.
//
//  Discussion:
//
//    This routine can be used to check whether the integrand function
//    for a given problem has singularities at the vertices or along
//    the edges of the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Output, int P00_SINGULARITY.
//    0, there are no vertex or edge singularities.
//    1, there are singularities at one or more vertices, but not on edges.
//    2, there are singularities on one or more edges, possibly 
//       including vertices.
//    3, there are singularities somewhere inside or on the triangle.
//
{
  int singularity;

  if ( problem == 1 )
  {
    singularity = 0;
  }
  else if ( problem == 2 )
  {
    singularity = 0;
  }
  else if ( problem == 3 )
  {
    singularity = 0;
  }
  else if ( problem == 4 )
  {
    singularity = 0;
  }
  else if ( problem == 5 )
  {
    singularity = 0;
  }
  else if ( problem == 6 )
  {
    singularity = 0;
  }
  else if ( problem == 7 )
  {
    singularity = 0;
  }
  else if ( problem == 8 )
  {
    singularity = 0;
  }
  else if ( problem == 9 )
  {
    singularity = 0;
  }
  else if ( problem == 10 )
  {
    singularity = 2;
  }
  else if ( problem == 11 )
  {
    singularity = 1;
  }
  else if ( problem == 12 )
  {
    singularity = 2;
  }
  else if ( problem == 13 )
  {
    singularity = 2;
  }
  else if ( problem == 14 )
  {
    singularity = 2;
  }
  else if ( problem == 15 )
  {
    singularity = 2;
  }
  else if ( problem == 16 )
  {
    singularity = 2;
  }
  else if ( problem == 17 )
  {
    singularity = 3;
  }
  else if ( problem == 18 )
  {
    singularity = 1;
  }
  else if ( problem == 19 )
  {
    singularity = 0;
  }
  else if ( problem == 20 )
  {
    singularity = 0;
  }
  else if ( problem == 21 )
  {
    singularity = 1;
  }
  else if ( problem == 22 )
  {
    singularity = 1;
  }
  else
  {
    cout << "\n";
    cout << "P00_SINGULARITY - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return singularity;
}
//****************************************************************************80

char *p00_title ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns the title for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the number of the desired test problem.
//
//    Output, char *P00_TITLE, the title.
//
{
  char *value;

  if ( problem == 1 )
  {
    value = p01_title ( );
  }
  else if ( problem == 2 )
  {
    value = p02_title ( );
  }
  else if ( problem == 3 )
  {
    value = p03_title ( );
  }
  else if ( problem == 4 )
  {
    value = p04_title ( );
  }
  else if ( problem == 5 )
  {
    value = p05_title ( );
  }
  else if ( problem == 6 )
  {
    value = p06_title ( );
  }
  else if ( problem == 7 )
  {
    value = p07_title ( );
  }
  else if ( problem == 8 )
  {
    value = p08_title ( );
  }
  else if ( problem == 9 )
  {
    value = p09_title ( );
  }
  else if ( problem == 10 )
  {
    value = p10_title ( );
  }
  else if ( problem == 11 )
  {
    value = p11_title ( );
  }
  else if ( problem == 12 )
  {
    value = p12_title ( );
  }
  else if ( problem == 13 )
  {
    value = p13_title ( );
  }
  else if ( problem == 14 )
  {
    value = p14_title ( );
  }
  else if ( problem == 15 )
  {
    value = p15_title ( );
  }
  else if ( problem == 16 )
  {
    value = p16_title ( );
  }
  else if ( problem == 17 )
  {
    value = p17_title ( );
  }
  else if ( problem == 18 )
  {
    value = p18_title ( );
  }
  else if ( problem == 19 )
  {
    value = p19_title ( );
  }
  else if ( problem == 20 )
  {
    value = p20_title ( );
  }
  else if ( problem == 21 )
  {
    value = p21_title ( );
  }
  else if ( problem == 22 )
  {
    value = p22_title ( );
  }
  else
  {
    cout << "\n";
    cout << "P00_TITLE - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

void p00_vertex_sub ( int problem, int level, int *n, double *result )

//****************************************************************************80
//
//  Purpose:
//
//   P00_VERTEX_SUB approximates an integral in a triangle by subdivision.
//
//  Discussion:
//
//    The function f(x,y) is to be integrated over a triangle T.
//
//    The first approximation averages the values at the vertices.
//
//    If a second approximation is requested, the routine subdivides each
//    existing triangle into 4, evaluates the function at the new vertices,
//    and returns an improved estimate.
//
//    The routine may be called repeatedly in this way, to get an improved
//    estimate of the integral.
//
//    Note that this routine will fail in the case that there
//    are singularities at the vertices or along the sides of the triangle.
//
//    Moreover, since the number of new vertices grows as a power of 4,
//    the use of an automatic array to store all the new vertices at one
//    time may fail when a memory limit is reached.
//
//    Finally, note that the rule has a very low order of convergence.
//    We're just exhibiting this routine as an EXAMPLE of how one
//    would use the test integrands.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Input, int LEVEL, the level of subdivision.  The first call
//    should be with LEVEL = 0.  For successive refinement, the routine
//    may be called repeatedly.  Each time, the user should increase the
//    value of LEVEL by 1, and also input the value of RESULT that was
//    output on the previous call.
//
//    Input/output, int *N, the number of function evaluations used.
//    If LEVEL = 0, the input value is ignored.  Otherwise, the input
//    value is assumed to be the output value from the previous call.
//
//    Input/output, double *RESULT, the approximate integral.
//    If LEVEL = 0, then the input value is ignored.  Otherwise, the
//    input value is assumed to be the result from the previous call,
//    at the previous level.  The output value is based on the input
//    value, adjusted by information determined at the new level.
//
{
  double area;
  double *f;
  int i;
  int j;
  int n_new;
  int order_max_1d;
  double *p;
  double result_new;
  double *t;
  double x;
  double xsi[3];
  double y;

  t = p00_vertices ( problem );
  area = triangle_area ( t );
//
//  Compute first level.
//
  if ( level == 0 )
  {
    *n = 3;
    p = new double[2 * (*n)];
    r8vec_copy ( 2 * (*n), t, p );

    f = p00_fun ( problem, *n, p );

    *result = r8vec_sum ( *n, f ) * area / 3.0;
  }
//
//  Compute next level.
//
  else
  {
    if ( level == 1 )
    {
      n_new = 3;
    }
    else
    {
      n_new = ( i4_power ( 2, level - 1 ) + 1 ) 
              * i4_power ( 2, level - 1 ) * 3;
    }

    p = new double[2*n_new];

    order_max_1d = i4_power ( 2, level );

    n_new = 0;

    for ( i = 0; i <= order_max_1d - 1; i = i + 2 )
    {
      for ( j = 0; j <= order_max_1d - 1 - i; j = j + 2 )
      {
        xsi[0] = ( double ) ( order_max_1d - i - 1 - j ) 
               / ( double ) ( order_max_1d             );
        xsi[1] = ( double ) (                i + 1     ) 
               / ( double ) ( order_max_1d             );
        xsi[2] = ( double ) (                        j ) 
               / ( double ) ( order_max_1d             );

        r8mat_mxv ( 2, 3, t, xsi, p+n_new*2 );

        n_new = n_new + 1;

        xsi[0] = ( double ) ( order_max_1d - i - j - 1 ) 
               / ( double ) ( order_max_1d             );
        xsi[1] = ( double ) (                i         )
               / ( double ) ( order_max_1d             );
        xsi[2] = ( double ) (                    j + 1 )
               / ( double ) ( order_max_1d	       );

        r8mat_mxv ( 2, 3, t, xsi, p+n_new*2 );

        n_new = n_new + 1;

        xsi[0] = ( double ) ( order_max_1d - i - 1 - j - 1 )
               / ( double ) ( order_max_1d                 );
        xsi[1] = ( double ) (                i + 1         )
               / ( double ) ( order_max_1d		   );
        xsi[2] = ( double ) (                        j + 1 )
               / ( double ) ( order_max_1d                 );

        r8mat_mxv ( 2, 3, t, xsi, p+n_new*2 );

        n_new = n_new + 1;
      }
    }

    f = p00_fun ( problem, n_new, p );

    result_new = r8vec_sum ( n_new, f ) * area / ( double ) ( n_new );

    *result = ( 3.0 * result_new + *result ) / 4.0;

    *n = *n + n_new;
  }

  delete [] f;
  delete [] p;
  delete [] t;

  return;
}
//****************************************************************************80

double *p00_vertices ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_VERTICES returns the vertices for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer PROBLEM, the problem number.
//
//    Output, real ( kind = 8 ) P00_VERTICES[2*3], the vertices.
//
{
  double *t;

  if ( problem == 1 )
  {
    t = p01_vertices ( );
  }
  else if ( problem == 2 )
  {
    t = p02_vertices ( );
  }
  else if ( problem == 3 )
  {
    t = p03_vertices ( );
  }
  else if ( problem == 4 )
  {
    t = p04_vertices ( );
  }
  else if ( problem == 5 )
  {
    t = p05_vertices ( );
  }
  else if ( problem == 6 )
  {
    t = p06_vertices ( );
  }
  else if ( problem == 7 )
  {
    t = p07_vertices ( );
  }
  else if ( problem == 8 )
  {
    t = p08_vertices ( );
  }
  else if ( problem == 9 )
  {
    t = p09_vertices ( );
  }
  else if ( problem == 10 )
  {
    t = p10_vertices ( );
  }
  else if ( problem == 11 )
  {
    t = p11_vertices ( );
  }
  else if ( problem == 12 )
  {
    t = p12_vertices ( );
  }
  else if ( problem == 13 )
  {
    t = p13_vertices ( );
  }
  else if ( problem == 14 )
  {
    t = p14_vertices ( );
  }
  else if ( problem == 15 )
  {
    t = p15_vertices ( );
  }
  else if ( problem == 16 )
  {
    t = p16_vertices ( );
  }
  else if ( problem == 17 )
  {
    t = p17_vertices ( );
  }
  else if ( problem == 18 )
  {
    t = p18_vertices ( );
  }
  else if ( problem == 19 )
  {
    t = p19_vertices ( );
  }
  else if ( problem == 20 )
  {
    t = p20_vertices ( );
  }
  else if ( problem == 21 )
  {
    t = p21_vertices ( );
  }
  else if ( problem == 22 )
  {
    t = p22_vertices ( );
  }
  else
  {
    cout << "\n";
    cout << "P00_VERTICES - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return t;
}
//****************************************************************************80

double p00_wandzura05_sub ( int problem, int level, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    P00_WANDZURA05_SUB uses subdivision and a Wandzura rule.
//
//  Discussion:
//
//    The Wandzura rule is a seven point rule of polynomial exactness 5.
//
//    The function f(x,y) is to be integrated over a triangle T.
//
//    The triangle is subdivided by subdividing each side into LEVEL sections,
//    which produces LEVEL*LEVEL subtriangles.  The Wandzura rule is then
//    applied to each subtriangle, and the result is summed.
//
//    The abscissas of this Wandzura rule do not lie on the vertices
//    or sides of the reference triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, pages 1829-1840, 2003.
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Input, integer LEVEL, the level of subdivision.  This indicates the
//    number of equally spaced subedges into which each edge of the triangle
//    is to be divided.  This will result in a total of LEVEL*LEVEL subtriangles
//    being used.
//
//    Output, int *N, the number of function evaluations used.
//
//    Output, double P00_WANDZURA05_SUB, the approximate integral.
//
{
# define ORDER 7

  double area;
  double *f;
  int i1;
  int i2;
  int i3;
  int j1;
  int j2;
  int j3;
  bool more;
  double result;
  double sub_tri_phys[2*3];
  double *tri_phys;
  double w[ORDER] = {
    0.22500000000000, 
    0.13239415278851, 
    0.13239415278851, 
    0.13239415278851, 
    0.12593918054483, 
    0.12593918054483, 
    0.12593918054483  
  };
  double xsi[3*3];
  double xy_phys[2*ORDER];
  double xy_ref[3*ORDER] = {
      0.33333333333333, 0.33333333333333, 0.33333333333333, 
      0.05971587178977, 0.47014206410512, 0.47014206410512, 
      0.47014206410512, 0.05971587178977, 0.47014206410512, 
      0.47014206410512, 0.47014206410512, 0.05971587178977, 
      0.79742698535309, 0.10128650732346, 0.10128650732346, 
      0.10128650732346, 0.79742698535309, 0.10128650732346, 
      0.10128650732346, 0.10128650732346, 0.79742698535309  
  };

  result = 0.0;

  tri_phys = p00_vertices ( problem );

  area = triangle_area ( tri_phys );

  more = false;

  for ( ; ; )
  {
//
//  Get the integer indices of the next reference subtriangle.
//
    subtriangle_next ( level, &more, &i1, &j1, &i2, &j2, &i3, &j3 );
//
//  Get the barycentric coordinates of the vertices of the reference subtriangle.
//  
    xsi[0+0*3] = ( double ) (         i1      ) / ( double ) ( level );
    xsi[1+0*3] = ( double ) (              j1 ) / ( double ) ( level );
    xsi[2+0*3] = ( double ) ( level - i1 - j1 ) / ( double ) ( level );

    xsi[0+1*3] = ( double ) (         i2      ) / ( double ) ( level );
    xsi[1+1*3] = ( double ) (              j2 ) / ( double ) ( level );
    xsi[2+1*3] = ( double ) ( level - i2 - j2 ) / ( double ) ( level );

    xsi[0+2*3] = ( double ) (         i3      ) / ( double ) ( level );
    xsi[1+2*3] = ( double ) (              j3 ) / ( double ) ( level );
    xsi[2+2*3] = ( double ) ( level - i3 - j3 ) / ( double ) ( level );
//
//  Map the reference subtriangle to the physical subtriangle.
//
    r8mat_mxm ( 2, 3, 3, tri_phys, xsi, sub_tri_phys );
//
//  Now map the integration abscissas to the physical subtriangle.
//
    r8mat_mxm ( 2, 3, ORDER, sub_tri_phys, xy_ref, xy_phys );
//
//  Evaluate the function.
//
    f = p00_fun ( problem, ORDER, xy_phys );
//
//  Update the quadrature estimate.
//
    result = result + r8vec_dot ( ORDER, w, f );

    delete [] f;

    if ( !more )
    {
      break;
    }
  }
//
//  Scale by area and number of subtriangles.
//
  *n = level * level * ORDER;

  result = result * area / ( double ) ( level * level );

  delete [] tri_phys;

  return result;
# undef ORDER
}
//****************************************************************************80

double *p01_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_FUN evaluates the integrand for problem 1.
//
//  Integrand:
//
//    f(x,y) = 2
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 2.0;
  }

  return f;
}
//****************************************************************************80

char *p01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns the name of problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P01_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[11];

  strcpy ( value, "f(x,y) = 2" );

  return value;
}
//****************************************************************************80

double *p01_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_VERTICES returns the vertices for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p02_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_FUN evaluates the integrand for problem 2.
//
//  Integrand:
//
//    f(x,y) = 6 * x
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P02_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 6.0 * p[0+i*2];
  }

  return f;
}
//****************************************************************************80

char *p02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_TITLE returns the name of problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P02_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[15];

  strcpy ( value, "f(x,y) = 6 * x" );

  return value;
}
//****************************************************************************80

double *p02_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_VERTICES returns the vertices for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p03_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_FUN evaluates the integrand for problem 3.
//
//  Integrand:
//
//    f(x,y) = 6 * y
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P03_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 6.0 * p[1+i*2];
  }

  return f;
}
//****************************************************************************80

char *p03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_TITLE returns the name of problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P03_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[15];

  strcpy ( value, "f(x,y) = 6 * y" );

  return value;
}
//****************************************************************************80

double *p03_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_VERTICES returns the vertices for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p04_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_FUN evaluates the integrand for problem 4.
//
//  Integrand:
//
//    f(x,y) = 12 * x^2
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P04_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 12.0 * pow ( p[0+i*2], 2 );
  }

  return f;
}
//****************************************************************************80

char *p04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_TITLE returns the name of problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P04_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[18];

  strcpy ( value, "f(x,y) = 12 * x^2" );

  return value;
}
//****************************************************************************80

double *p04_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_VERTICES returns the vertices for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p05_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_FUN evaluates the integrand for problem 5.
//
//  Integrand:
//
//    f(x,y) = 24 * x*y
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P05_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 24.0 * p[0+i*2] * p[1+i*2];
  }

  return f;
}
//****************************************************************************80

char *p05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_TITLE returns the name of problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P01_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[18];

  strcpy ( value, "f(x,y) = 24 * x*y" );

  return value;
}
//****************************************************************************80

double *p05_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_VERTICES returns the vertices for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p06_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_FUN evaluates the integrand for problem 6.
//
//  Integrand:
//
//    f(x,y) = 12 * y^2
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P06_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 12.0 * pow ( p[1+i*2], 2 );
  }

  return f;
}
//****************************************************************************80

char *p06_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_TITLE returns the name of problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P06_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[18];

  strcpy ( value, "f(x,y) = 12 * y^2" );

  return value;
}
//****************************************************************************80

double *p06_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_VERTICES returns the vertices for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p07_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_FUN evaluates the integrand for problem 7.
//
//  Integrand:
//
//    f(x,y) = 20 * x^3
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P07_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 20.0 * pow ( p[0+i*2], 3 );
  }

  return f;
}
//****************************************************************************80

char *p07_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_TITLE returns the name of problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P07_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[18];

  strcpy ( value, "f(x,y) = 20 * x^3" );

  return value;
}
//****************************************************************************80

double *p07_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_VERTICES returns the vertices for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p08_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_FUN evaluates the integrand for problem 8.
//
//  Integrand:
//
//    f(x,y) = 30 * x^4
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P08_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 30.0 * pow ( p[0+i*2], 4 );
  }

  return f;
}
//****************************************************************************80

char *p08_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_TITLE returns the name of problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P08_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[18];

  strcpy ( value, "f(x,y) = 30 * x^4" );

  return value;
}
//****************************************************************************80

double *p08_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_VERTICES returns the vertices for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p09_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P09_FUN evaluates the integrand for problem 9.
//
//  Integrand:
//
//    f(x,y) = 42 * x^5
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P09_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 42.0 * pow ( p[0+i*2], 5 );
  }

  return f;
}
//****************************************************************************80

char *p09_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_TITLE returns the name of problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P01_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[18];

  strcpy ( value, "f(x,y) = 42 * x^5" );

  return value;
}
//****************************************************************************80

double *p09_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_VERTICES returns the vertices for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p10_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P10_FUN evaluates the integrand for problem 10.
//
//  Discussion:
//
//    The integral has been transformed from the integral of G(X,Y)
//    over the unit reference triangle.
//
//  Integrand:
//
//    PA = 1
//    PB = 5
//    PC = 0
//    PD = 0
//    PG = 0.25
//    PH = -0.25
//    D = PB - PA
//    U(X) = ( X - PA ) / D
//    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
//    V(X,Y) = V1(X) * ( Y - PC * X - PD )
//
//    G(X,Y) = X^(-0.2)
//
//    f(x,y) = g ( u(x), v(x,y) ) * v1(x) / d
//
//  Vertices:
//
//    (1,0), (5,0), (5,1)
//
//  Integral:
//
//    0.6944444444444444D+00
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P10_FUN[N], the function values.
//
{
  double d;
  double *f;
  int i;
  double pa = 1.0;
  double pb = 5.0;
  double pc = 0.0;
  double pd = 0.0;
  double pg = 0.25;
  double ph = -0.25;
  double power = -0.2;
  double u;
  double v;
  double v1;

  f = new double[n];

  d = pb - pa;

  for ( i = 0; i < n; i++ )
  {
    u = ( p[0+i*2] - pa ) / d;
    v1 = ( 1.0 - u ) / ( ( pg - pc ) * p[0+i*2] + ph - pd );
    v = v1 * ( p[1+i*2] - pc * p[0+i*2] - pd );
    f[i] = ( 36.0 / 25.0 ) * pow ( u, power ) * v1 / d;
  }

  return f;
}
//****************************************************************************80

char *p10_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_TITLE returns the name of problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P10_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[41];

  strcpy ( value, "f(x,y) = x^(-0.2) on ((1,0),(5,0),(5,1))" );

  return value;
}
//****************************************************************************80

double *p10_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_VERTICES returns the vertices for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 1.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 5.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 5.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p11_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P11_FUN evaluates the integrand for problem 11.
//
//  Discussion:
//
//    The integral has been transformed from the integral of G(X,Y)
//    over the unit reference triangle.
//
//  Integrand:
//
//    PA = 0
//    PB = 1
//    PC = 0
//    PD = 0
//    PG = -1.0
//    PH = 1.0
//    D = PB - PA
//    U(X) = ( X - PA ) / D
//    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
//    V(X,Y) = V1(X) * ( Y - PC * X - PD )
//
//    G(X,Y) = (X+Y)^(-0.2)
//
//    f(x,y) = g ( u(x), v(x,y) ) * v1(x) / d
//
//  Vertices:
//
//    (0,1), (1,0), (0,1)
//
//  Integral:
//
//    0.5555555555555556D+00
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P11_FUN[N], the function values.
//
{
  double d;
  double *f;
  int i;
  double pa = 0.0;
  double pb = 1.0;
  double pc = 0.0;
  double pd = 0.0;
  double pg = -1.00;
  double ph =  1.00;
  double power = -0.2;
  double u;
  double v;
  double v1;

  f = new double[n];

  d = pb - pa;

  for ( i = 0; i < n; i++ )
  {
    u = ( p[0+i*2] - pa ) / d;
    v1 = ( 1.0 - u ) / ( ( pg - pc ) * p[0+i*2] + ph - pd );
    v = v1 * ( p[1+i*2] - pc * p[0+i*2] - pd );
    f[i] = ( 9.0 / 5.0 ) * pow ( u + v, power ) * v1 / d;
  }

  return f;
}
//****************************************************************************80

char *p11_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_TITLE returns the name of problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P11_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[22];

  strcpy ( value, "f(x,y) = (x+y)^(-0.2)" );

  return value;
}
//****************************************************************************80

double *p11_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_VERTICES returns the vertices for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p12_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P12_FUN evaluates the integrand for problem 12.
//
//  Discussion:
//
//    The integral has been transformed from the integral of G(X,Y)
//    over the unit reference triangle.
//
//  Integrand:
//
//    PA = -1
//    PB = 3
//    PC = 0.25
//    PD = -2.75
//    PG = -1.0
//    PH = 1.0
//    D = PB - PA
//    U(X) = ( X - PA ) / D
//    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
//    V(X,Y) = V1(X) * ( Y - PC * X - PD )
//
//    G(X,Y) = (1-X-Y)^(-0.2)
//
//    f(x,y) = g ( u(x), v(x,y) ) * v1(x) / d
//
//  Vertices:
//
//    (-1,-3), (3,-2), (-1,2)
//
//  Integral:
//
//    0.6944444444444444D+00
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P12_FUN[N], the function values.
//
{
  double d;
  double *f;
  int i;
  double pa = -1.0;
  double pb = 3.0;
  double pc = 0.25;
  double pd = -2.75;
  double pg = -1.00;
  double ph =  1.00;
  double power = -0.2;
  double u;
  double v;
  double v1;

  f = new double[n];

  d = pb - pa;

  for ( i = 0; i < n; i++ )
  {
    u = ( p[0+i*2] - pa ) / d;
    v1 = ( 1.0 - u ) / ( ( pg - pc ) * p[0+i*2] + ph - pd );
    v = v1 * ( p[1+i*2] - pc * p[0+i*2] - pd );
    f[i] = ( 36.0 / 25.0 ) * pow ( 1.0 - u - v, power ) * v1 / d;
  }

  return f;
}
//****************************************************************************80

char *p12_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P12_TITLE returns the name of problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P12_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[51];

  strcpy ( value, "f(x,y) = (1-x-y)^(-0.2) on ((-1,-3),(3,-2),(-1,2))" );

  return value;
}
//****************************************************************************80

double *p12_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P12_VERTICES returns the vertices for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = -1.0;
  t[1+0*2] = -3.0;
  t[0+1*2] =  3.0;
  t[1+1*2] = -2.0;
  t[0+2*2] = -1.0;
  t[1+2*2] =  2.0;

  return t;
}
//****************************************************************************80

double *p13_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P13_FUN evaluates the integrand for problem 13.
//
//  Discussion:
//
//    The integral has been transformed from the integral of G(X,Y)
//    over the unit reference triangle.
//
//  Integrand:
//
//    PA = 0
//    PB = -7
//    PC = 0
//    PD = 0
//    PG = -3/7
//    PH = -3
//    D = PB - PA
//    U(X) = ( X - PA ) / D
//    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
//    V(X,Y) = V1(X) * ( Y - PC * X - PD )
//
//    G(X,Y) = (X*Y)^(-0.2)
//
//    f(x,y) = g ( u(x), v(x,y) ) * v1(x) / d
//
//  Vertices:
//
//    (0,0), (-7,0), (0,-3)
//
//  Integral:
//
//    0.94810264549557699446
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P13_FUN[N], the function values.
//
{
  double c = 0.94810264549557699446;
  double d;
  double *f;
  int i;
  double pa = 0.0;
  double pb = -7.0;
  double pc = 0.0;
  double pd = 0.0;
  double pg = ( -3.00 / 7.00 );
  double ph =  -3.00;
  double power = -0.2;
  double u;
  double v;
  double v1;

  f = new double[n];

  d = pb - pa;

  for ( i = 0; i < n; i++ )
  {
    u = ( p[0+i*2] - pa ) / d;
    v1 = ( 1.0 - u ) / ( ( pg - pc ) * p[0+i*2] + ph - pd );
    v = v1 * ( p[1+i*2] - pc * p[0+i*2] - pd );
    f[i] = pow ( u * v, power ) * v1 / d / c;
  }

  return f;
}
//****************************************************************************80

char *p13_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P13_TITLE returns the name of problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P13_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[47];

  strcpy ( value, "f(x,y) = (x*y)^(-0.2) on ((0,0),(-7,0),(0,-3))" );

  return value;
}
//****************************************************************************80

double *p13_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P13_VERTICES returns the vertices for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] =  0.0;
  t[1+0*2] =  0.0;
  t[0+1*2] = -7.0;
  t[1+1*2] =  0.0;
  t[0+2*2] =  0.0;
  t[1+2*2] = -3.0;

  return t;
}
//****************************************************************************80

double *p14_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P14_FUN evaluates the integrand for problem 14.
//
//  Integrand:
//
//    f(x,y) = 1 / sqrt ( X ) + 1 / sqrt ( Y ) + 1 / sqrt ( X + Y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P14_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 3.0 / 10.0 * (      
           1.0 / sqrt ( p[0+i*2] ) 
         + 1.0 / sqrt ( p[1+i*2] ) 
         + 1.0 / sqrt ( p[0+i*2] + p[1+i*2] ) );
  }

  return f;
}
//****************************************************************************80

char *p14_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P14_TITLE returns the name of problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P14_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[45];

  strcpy ( value, "f(x,y) = 1/sqrt(x) + 1/sqrt(y) + 1/sqrt(x+y)" );

  return value;
}
//****************************************************************************80

double *p14_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P14_VERTICES returns the vertices for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p15_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P15_FUN evaluates the integrand for problem 15.
//
//  Integrand:
//
//    f(x,y) = 1 / sqrt ( 1 - X - Y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P15_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 3.0 / 4.0 / sqrt ( 1.0 - p[0+i*2] - p[1+i*2] );
  }

  return f;
}
//****************************************************************************80

char *p15_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P15_TITLE returns the name of problem 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P15_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[23];

  strcpy ( value, "f(x,y) = 1/sqrt(1-x-y)" );

  return value;
}
//****************************************************************************80

double *p15_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P15_VERTICES returns the vertices for problem 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p16_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P16_FUN evaluates the integrand for problem 16.
//
//  Integrand:
//
//    f(x,y) = log ( x * y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P16_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = ( - 2.0 / 3.0 ) * log ( p[0+i*2] * p[1+i*2] );
  }

  return f;
}
//****************************************************************************80

char *p16_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P16_TITLE returns the name of problem 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P16_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[18];

  strcpy ( value, "f(x,y) = log(x*y)" );

  return value;
}
//****************************************************************************80

double *p16_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P16_VERTICES returns the vertices for problem 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p17_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P17_FUN evaluates the integrand for problem 17.
//
//  Integrand:
//
//    f(x,y) = 1/sqrt(|x-1/4|) + 1/sqrt(|y-1/2|)
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P17_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = ( 1.0 / 3.11357229949 ) * 
    ( 1.0 / sqrt ( r8_abs ( p[0+i*2] - 0.25 ) ) 
    + 1.0 / sqrt ( r8_abs ( p[1+i*2] - 0.50 ) ) );
  }

  return f;
}
//****************************************************************************80

char *p17_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P17_TITLE returns the name of problem 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P17_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[43];

  strcpy ( value, "f(x,y) = 1/sqrt(|x-1/4|) + 1/sqrt(|y-1/2|)" );

  return value;
}
//****************************************************************************80

double *p17_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P17_VERTICES returns the vertices for problem 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p18_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P18_FUN evaluates the integrand for problem 18.
//
//  Integrand:
//
//    f(x,y) = log ( x + y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P18_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = -4.0 * log ( p[0+i*2] + p[1+i*2] );
  }

  return f;
}
//****************************************************************************80

char *p18_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P18_TITLE returns the name of problem 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P18_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[23];

  strcpy ( value, "f(x,y) = log ( x + y )" );

  return value;
}
//****************************************************************************80

double *p18_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P18_VERTICES returns the vertices for problem 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p19_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P19_FUN evaluates the integrand for problem 19.
//
//  Integrand:
//
//    f(x,y) = sin ( x ) * cos ( 5 * y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P19_FUN[N], the function values.
//
{
  double c = 0.043052326655855175018;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = sin ( p[0+i*2] ) * cos ( 5.0 * p[1+i*2] ) / c;
  }

  return f;
}
//****************************************************************************80

char *p19_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P19_TITLE returns the name of problem 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P19_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[31];

  strcpy ( value, "f(x,y) = sin ( x ) cos ( 5 y )" );

  return value;
}
//****************************************************************************80

double *p19_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P19_VERTICES returns the vertices for problem 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p20_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P20_FUN evaluates the integrand for problem 20.
//
//  Integrand:
//
//    f(x,y) = sin ( 11 x ) * cos ( y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P20_FUN[N], the function values.
//
{
  double c = 0.085468091995313041919;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = sin ( 11.0 * p[0+i*2] ) * cos ( p[1+i*2] ) / c;
  }

  return f;
}
//****************************************************************************80

char *p20_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P20_TITLE returns the name of problem 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P20_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[32];

  strcpy ( value, "f(x,y) = sin ( 11 x ) cos ( y )" );

  return value;
}
//****************************************************************************80

double *p20_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P20_VERTICES returns the vertices for problem 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p21_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P21_FUN evaluates the integrand for problem 21.
//
//  Discussion:
//
//    To do this integral by hand, convert to polar coordinates:
//
//    Integral ( 0 <= t <= Pi/2 )
//      Integral ( 0 <= r <= 1/(cos(t)+sin(t)) ) 1/r * r dr dt
//
//  Integrand:
//
//    f(x,y) = 1 / sqrt ( x * x + y * y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P21_FUN[N], the function values.
//
{
  double c = 1.2464504802804610268;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0 / sqrt ( pow ( p[0+i*2], 2 ) + pow ( p[1+i*2], 2 ) ) / c;
  }

  return f;
}
//****************************************************************************80

char *p21_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P21_TITLE returns the name of problem 21.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P21_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[40];

  strcpy ( value, "f(x,y) = 1 / r = 1 / sqrt ( x^2 + y^2 )" );

  return value;
}
//****************************************************************************80

double *p21_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P21_VERTICES returns the vertices for problem 21.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double *p22_fun ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    P22_FUN evaluates the integrand for problem 22.
//
//  Discussion:
//
//    To do this integral by hand, convert to polar coordinates:
//
//    Integral ( 0 <= t <= Pi/2 )
//      Integral ( 0 <= r <= 1/(cos(t)+sin(t))) Log(r)/r * r dr dt
//
//  Integrand:
//
//    f(x,y) = log ( sqrt ( x * x + y * y ) ) / sqrt ( x * x + y * y )
//
//  Vertices:
//
//    (0,0), (1,0), (0,1)
//
//  Integral:
//
//    1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double P22_FUN[N], the function values.
//
{
  double c = -1.5280234546641884580;
  double *f;
  int i;
  double r;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    r = sqrt ( pow ( p[0+i*2], 2 ) + pow ( p[1+i*2], 2 ) );
    f[i] = log ( r ) / r / c;
  }

  return f;
}
//****************************************************************************80

char *p22_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P22_TITLE returns the name of problem 22.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *P22_TITLE, the title of the problem.
//
{
  char *value;

  value = new char[23];

  strcpy ( value, "f(x,y) = log ( r ) / r" );

  return value;
}
//****************************************************************************80

double *p22_vertices ( )

//****************************************************************************80
//
//  Purpose:
//
//    P22_VERTICES returns the vertices for problem 22.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double T[2*3], the vertices.
//
{
  double *t;

  t = new double[2*3];

  t[0+0*2] = 0.0;
  t[1+0*2] = 0.0;
  t[0+1*2] = 1.0;
  t[1+1*2] = 0.0;
  t[0+2*2] = 0.0;
  t[1+2*2] = 1.0;

  return t;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_mxm ( int n1, int n2, int n3, double a[], double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MXM multiplies two matrices.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double C[N1*N3], the product matrix C = A * B.
//
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return;
}
//****************************************************************************80

void r8mat_mxv ( int m, int n, double a[], double x[], double ax[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MXV multiplies a matrix times a vector.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double AX[M], the product A*X.
//
{
  int i;
  int j;

  for ( i = 0; i < m; i++ )
  {
    ax[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      ax[i] = ax[i] + a[i+j*m] * x[j];
    }
  }

  return;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Input, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }

  return value;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double sum;

  sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

double *r8vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
  int i;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge ( );
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

unsigned long random_initialize ( unsigned long seed )

//****************************************************************************80
//
//  Purpose:
//
//    RANDOM_INITIALIZE initializes the RANDOM random number generator.
//
//  Discussion:
//
//    If you don't initialize RANDOM, the random number generator, 
//    it will behave as though it were seeded with value 1.  
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRANDOM (the appropriate routine 
//    to call when setting the seed for RANDOM).  The seed is also
//    returned to the user as the value of the function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned long SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RANDOM_INITIALIZE, is the value of the seed
//    passed to SRANDOM, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
# define DEBUG 0

  if ( seed != 0 )
  {
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with user SEED = " << seed << "\n";
    }
  }
  else
  {
    seed = get_seed ( );
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with arbitrary SEED = " << seed << "\n";
    }
  }
//
//  Now set the seed.
//
  srandom ( seed );

  return seed;
# undef DEBUG
}
//****************************************************************************80

void subtriangle_next ( int n, bool *more, int *i1, int *j1, int *i2, int *j2, 
  int *i3, int *j3 )

//****************************************************************************80
//
//  Purpose:
//
//    SUBTRIANGLE_NEXT computes the next subtriangle of a triangle.
//
//  Discussion:
//
//    The three sides of a triangle have been subdivided into N segments,
//    inducing a natural subdivision of the triangle into N*N subtriangles.
//    It is desired to consider each subtriangle, one at a time, in some
//    definite order.  This routine can produce information defining each 
//    of the subtriangles, one after another.
//
//    The subtriangles are described in terms of the integer coordinates 
//    (I,J) of their vertices.  These coordinates both range from 0 to N,
//    with the additional restriction that I + J <= N.
//
//    The vertices of each triangle are listed in counterclockwise order.
//
//  Example:
//
//    N = 4
//
//    4  *
//       |\
//       16\
//    3  *--*
//       |14|\
//       13\15\
//    2  *--*--*
//       |\9|11|\
//       |8\10\12\
//    1  *--*--*--*
//       |\2|\4|\6|\
//       |1\|3\|5\|7\
//   0   *--*--*--*--*
//
//       0  1  2  3  4
//
//    Rank  I1 J1  I2 J2  I3 J3
//    ----  -----  -----  ----- 
//       1   0  0   1  0   0  1
//       2   1  1   0  1   1  0
//       3   1  0   2  0   1  1
//       4   2  1   1  1   2  0
//       5   2  0   3  0   2  1
//       6   3  1   1  1   3  0
//       7   3  0   4  0   3  1
//       8   0  1   1  1   0  2
//       9   1  2   0  2   1  1
//      10   1  1   2  1   1  2
//      11   2  2   1  2   2  1
//      12   2  1   3  1   2  2
//      13   0  2   1  2   0  3
//      14   1  3   0  3   1  2
//      15   1  2   2  2   1  3
//      16   0  3   1  3   0  4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, indicates the number of subdivisions of each side
//    of the original triangle.
//
//    Input/output, bool *MORE.
//    On first call, set MORE to FALSE.  Thereafter, the output value of MORE
//    will be TRUE if there are more subtriangles that can be generated by
//    further calls.  However, if MORE is returned as FALSE, the accompanying
//    subtriangle information refers to the last subtriangle that can be
//    generated.
//
//    Input/output, int *I1, *J1, *I2, *J2, *I3, *J3, the indices of the 
//    vertices of the subtriangle.
//
{
  if ( n <= 0 )
  {
    *more = false;
    return;
  }

  if ( ! ( *more ) )
  {
    *i1 = 0;
    *j1 = 0;
    *i2 = 1;
    *j2 = 0;
    *i3 = 0;
    *j3 = 1;

    if ( n == 1 )
    {
      *more = false;
    }
    else
    {
      *more = true;
    }
  }
//
//  We last generated a triangle like:
//
//    2---1
//     \  |
//      \ |
//       \|
//        3
//
  else if ( *i2 < *i3 )
  {
    *i1 = *i3;
    *j1 = *j3;
    *i2 = *i1 + 1;
    *j2 = *j1;
    *i3 = *i1;
    *j3 = *j1 + 1;
  }
//
//  We last generated a triangle like
//
//    3
//    |\
//    | \
//    |  \
//    1---2
//
  else if ( *i1 + 1 + *j1 + 1 <= n )
  {
    *i1 = *i1 + 1;
    *j1 = *j1 + 1;
    *i2 = *i1 - 1;
    *j2 = *j1;
    *i3 = *i1;
    *j3 = *j1 - 1;
  }
//
//  We must be at the end of a row.
//
  else
  {
    *i1 = 0;
    *j1 = *j1 + 1;
    *i2 = *i1 + 1;
    *j2 = *j1;
    *i3 = *i1;
    *j3 = *j1 + 1;

    if ( n <= *j1 + 1 )
    {
      *more = false;
    }
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
//    24 September 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

double triangle_area ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA computes the area of a triangle.
//
//  Discussion:
//
//    If the triangle's vertices are given in counter clockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed 
//    area of a triangle, it is possible to easily compute the area 
//    of a nonconvex polygon as the sum of the (possibly negative) 
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA, the area of the triangle.
//
{
  double area;

  area = 0.5 * ( 
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) + 
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) + 
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );
 
  return area;
}
//****************************************************************************80

double *triangle_sample ( double t[2*3], int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE returns random points in a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Input, integer N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double TRIANGLE_SAMPLE[2*N], a random point in the triangle.
//
{
# define DIM_NUM 2

  double alpha;
  double beta;
  int j;
  double r;
  double *p;
  double p12[DIM_NUM];
  double p13[DIM_NUM];

  p = new double[2*n];

  for ( j = 0; j < n; j++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    alpha = sqrt ( r );
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
    p12[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+1*2];
    p12[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+1*2];

    p13[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+2*2];;
    p13[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+2*2];;
//
//  Now choose, uniformly at random, a point on the line L.
//
    beta = r8_uniform_01 ( seed );

    p[0+j*2] = ( 1.0 - beta ) * p12[0] + beta * p13[0];
    p[1+j*2] = ( 1.0 - beta ) * p12[1] + beta * p13[1];
  }

  return p;
# undef DIM_NUM
}
