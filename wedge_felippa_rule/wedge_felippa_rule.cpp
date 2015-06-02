# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>
# include <cmath>

using namespace std;

# include "wedge_felippa_rule.hpp"

//****************************************************************************80

void comp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed 
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

  return;
}
//****************************************************************************80

void line_o01 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_O01 returns a 1 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[1], the weights.
//
//    Output, double X[1], the abscissas.
//
{
  int i;
  int order = 1;
  double w_save[1] = { 1.0 };
  double x_save[1] = { 0.0 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_o02 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_O02 returns a 2 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[2], the weights.
//
//    Output, double X[2], the abscissas.
//
{
  int i;
  int order = 2;
  double w_save[2] = {
    0.5,
    0.5 };
  double x_save[2] = {
    -0.57735026918962576451,
     0.57735026918962576451 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_o03 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_O03 returns a 3 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[3], the weights.
//
//    Output, double X[3], the abscissas.
//
{
  int i;
  int order = 3;
  double w_save[3] = {
    0.27777777777777777777, 
    0.44444444444444444444, 
    0.27777777777777777777  };
  double x_save[3] = {
    -0.77459666924148337704,
     0.00000000000000000000,
     0.77459666924148337704 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_o04 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_O04 returns a 4 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[4], the weights.
//
//    Output, double X[4], the abscissas.
//
{
  int i;
  int order = 4;
  double w_save[4] = {
    0.173927422568727, 
    0.326072577431273, 
    0.326072577431273, 
    0.173927422568727 };
  double x_save[4] = {
    -0.86113631159405257522,
    -0.33998104358485626480,
     0.33998104358485626480,
     0.86113631159405257522 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_o05 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_O05 returns a 5 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[5], the weights.
//
//    Output, double X[5], the abscissas.
//
{
  int i;
  int order = 5;
  double w_save[5] = {
    0.118463442528095, 
    0.239314335249683, 
    0.284444444444444, 
    0.239314335249683, 
    0.118463442528095 };
  double x_save[5] = {
    -0.90617984593866399280,
    -0.53846931010568309104,
     0.00000000000000000000,
     0.53846931010568309104,
     0.90617984593866399280 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

double *monomial_value ( int m, int n, int e[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= i <= m ) x(i)^e(i)
//
//    The combination 0.0^0 is encountered is treated as 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
//    Input, int E[M], the exponents.
//
//    Input, double X[M*N], the point coordinates.
//
//    Output, double MONOMIAL_VALUE[N], the monomial values.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];
  for ( j = 0; j < n; j++)
  {
    v[j] = 1.0;
  }
//v = r8vec_ones_new ( n );

  for ( i = 0; i < m; i++ )
  {
    if ( 0 != e[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], e[i] );
      }
    }
  }

  return v;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
//  Discussion:
//
//    An R8VEC is a vector of R8's.
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

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2007
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
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
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

double *r8vec_ones_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ONES_NEW creates a vector of 1's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ONES_NEW[N], a vector of 1's.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 1.0;
  }
  return a;
}
//****************************************************************************80

void subcomp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    SUBCOMP_NEXT computes the next subcomposition of N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to a value of N.
//
//    A subcomposition of the integer N into K parts is a composition
//    of M into K parts, where 0 <= M <= N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose subcompositions are desired.
//
//    Input, int K, the number of parts in the subcomposition.
//
//    Input/output, int A[K], the parts of the subcomposition.
//
//    Input/output, bool *MORE, set by the user to start the computation,
//    and by the routine to terminate it.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;
  static bool more2 = false;
  static int n2 = 0;
//
//  The first computation.
//
  if ( !( *more ) )
  {
    n2 = 0;

    for ( i = 0; i < k; i++ )
    {
      a[i] = 0;
    }
    more2 = false;
    *h = 0;
    *t = 0;

    *more = true;
  }
//
//  Do the next element at the current value of N.
//
  else if ( more2 )
  {
    comp_next ( n2, k, a, &more2, h, t );
  }
  else
  {
    more2 = false;
    n2 = n2 + 1;

    comp_next ( n2, k, a, &more2, h, t );
  }
//
//  Termination occurs if MORE2 = FALSE and N2 = N.
//
  if ( !more2 && n2 == n )
  {
    *more = false;
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
//****************************************************************************80

void triangle_o01 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_O01 returns a 1 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 1.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Output, double W[1], the weights.
//
//    Output, double XY[2*1], the abscissas.
//
{
  int order = 1;

  double w_save[1] = {
    1.0 };
  double xy_save[2*1] = {
    0.33333333333333333333,  0.33333333333333333333 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_o03 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_O03 returns a 3 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 2.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Output, double W[3], the weights.
//
//    Output, double XY[2*3], the abscissas.
//
{
  int order = 3;

  double w_save[3] = {
    0.33333333333333333333,
    0.33333333333333333333,
    0.33333333333333333333 };
  double xy_save[2*3] = {
    0.66666666666666666667,  0.16666666666666666667,
    0.16666666666666666667,  0.66666666666666666667,
    0.16666666666666666667,  0.16666666666666666667 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_o03b ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_O03B returns a 3 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 2.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Output, double W[3], the weights.
//
//    Output, double XY[2*3], the abscissas.
//
{
  int order = 3;

  double w_save[3] = {
    0.33333333333333333333,
    0.33333333333333333333,
    0.33333333333333333333 };
  double xy_save[2*3] = {
    0.0,  0.5,
    0.5,  0.0,
    0.5,  0.5 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_o06 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_O06 returns a 6 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 4.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Output, double W[6], the weights.
//
//    Output, double XY[2*6], the abscissas.
//
{
  int order = 6;

  double w_save[6] = {
    0.22338158967801146570,
    0.22338158967801146570,
    0.22338158967801146570,
    0.10995174365532186764,
    0.10995174365532186764,
    0.10995174365532186764 };
  double xy_save[2*6] = {
    0.10810301816807022736,  0.44594849091596488632,
    0.44594849091596488632,  0.10810301816807022736,
    0.44594849091596488632,  0.44594849091596488632,
    0.81684757298045851308,  0.091576213509770743460,
    0.091576213509770743460,  0.81684757298045851308,
    0.091576213509770743460,  0.091576213509770743460 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_o06b ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_O06B returns a 6 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 3.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Output, double W[6], the weights.
//
//    Output, double XY[2*6], the abscissas.
//
{
  int order = 6;

  double w_save[6] = {
    0.30000000000000000000,
    0.30000000000000000000,
    0.30000000000000000000,
    0.033333333333333333333,
    0.033333333333333333333,
    0.033333333333333333333 };
  double xy_save[2*6] = {
    0.66666666666666666667,  0.16666666666666666667,
    0.16666666666666666667,  0.66666666666666666667,
    0.16666666666666666667,  0.16666666666666666667,
    0.0,  0.5,
    0.5,  0.0,
    0.5,  0.5 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_o07 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_O07 returns a 7 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 5.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Output, double W[7], the weights.
//
//    Output, double XY[2*7], the abscissas.
//
{
  int order = 7;

  double w_save[7] = {
    0.12593918054482715260,
    0.12593918054482715260,
    0.12593918054482715260,
    0.13239415278850618074,
    0.13239415278850618074,
    0.13239415278850618074,
    0.22500000000000000000 };
  double xy_save[2*7] = {
    0.79742698535308732240,  0.10128650732345633880,
    0.10128650732345633880,  0.79742698535308732240,
    0.10128650732345633880,  0.10128650732345633880,
    0.059715871789769820459,  0.47014206410511508977,
    0.47014206410511508977,  0.059715871789769820459,
    0.47014206410511508977,  0.47014206410511508977,
    0.33333333333333333333,  0.33333333333333333333 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_o12 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_O12 returns a 12 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 6.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Output, double W[12], the weights.
//
//    Output, double XY[2*12], the abscissas.
//
{
  int order = 12;

  double w_save[12] = {
     0.050844906370206816921,
     0.050844906370206816921,
     0.050844906370206816921,
     0.11678627572637936603,
     0.11678627572637936603,
     0.11678627572637936603,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194 };
  double xy_save[2*12] = {
    0.87382197101699554332,  0.063089014491502228340,
    0.063089014491502228340,  0.87382197101699554332,
    0.063089014491502228340,  0.063089014491502228340,
    0.50142650965817915742,  0.24928674517091042129,
    0.24928674517091042129,  0.50142650965817915742,
    0.24928674517091042129,  0.24928674517091042129,
    0.053145049844816947353,  0.31035245103378440542,
    0.31035245103378440542,  0.053145049844816947353,
    0.053145049844816947353,  0.63650249912139864723,
    0.31035245103378440542,  0.63650249912139864723,
    0.63650249912139864723,  0.053145049844816947353,
    0.63650249912139864723,  0.31035245103378440542 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

double wedge_integral ( int expon[3] )

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE_INTEGRAL: monomial integral in a unit wedge.
//
//  Discussion:
//
//    This routine returns the integral of
//
//      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
//
//    over the unit wedge.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int EXPON[3], the exponents.
//
//    Output, double WEDGE_INTEGRAL, the integral of the monomial.
//
{
  int i;
  int k;
  double value;
//
//  The first computation ends with VALUE = 1.0;
//
  value = 1.0;

  k = expon[0];

  for ( i = 1; i <= expon[1]; i++ )
  {
    k = k + 1;
    value = value * ( double ) ( i ) / ( double ) ( k );
  }

  k = k + 1;
  value = value / ( double ) ( k );

  k = k + 1;
  value = value / ( double ) ( k );
//
//  Now account for integration in Z.
//
  if ( expon[2] == - 1 )
  {
    cerr << "\n";
    cerr << "WEDGE_INTEGRAL - Fatal error!\n";
    cerr << "  EXPON[2] = -1 is not a legal input.\n";
    exit ( 1 );
  }
  else if ( ( expon[2] % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = value * 2.0 / ( double ) ( expon[2] + 1 );
  }

  return value;
}
//****************************************************************************80

void wedge_rule ( int line_order, int triangle_order, double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE_RULE returns a quadrature rule for the unit wedge.
//
//  Discussion:
//
//    It is usually sensible to take LINE_ORDER and TRIG_ORDER so that
//    the line and triangle rules are roughly the same precision.  For that
//    criterion, we recommend the following combinations:
//
//      TRIANGLE_ORDER  LINE_ORDER  Precision
//      --------------  ----------  ---------
//          1               1       1 x 1
//          3               2       2 x 3
//         -3               2       2 x 3
//          6               3       4 x 5
//         -6               2       3 x 3
//          7               3       5 x 5
//         12               4       6 x 7
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, int LINE_ORDER, the index of the line rule.
//    The index of the rule is equal to the order of the rule.
//    1 <= LINE_ORDER <= 5.
//
//    Input, int TRIANGLE_ORDER, the indes of the triangle rule.
//    The index of the rule is 1, 3, -3, 6, -6, 7 or 12.
//
//    Output, double W[LINE_ORDER*abs(TRIANGLE_ORDER)], the weights.
//
//    Output, double XYZ[3*LINE_ORDER*abs(TRIANGLE_ORDER)], the abscissas.
//
{
  int i;
  int j;
  int k;
  double *line_w;
  double *line_x;
  double *triangle_w;
  double *triangle_xy;

  line_w = new double[line_order];
  line_x = new double[line_order];

  if ( line_order == 1 )
  {
    line_o01 ( line_w, line_x );
  }
  else if ( line_order == 2 )
  {
    line_o02 ( line_w, line_x );
  }
  else if ( line_order == 3 )
  {
    line_o03 ( line_w, line_x );
  }
  else if ( line_order == 4 )
  {
    line_o04 ( line_w, line_x );
  }
  else if ( line_order == 5 )
  {
    line_o05 ( line_w, line_x );
  }
  else
  {
    cerr << "\n";
    cerr << "WEDGE_RULE - Fatal error!\n";
    cerr << "  Illegal value of LINE_ORDER.\n";
    exit ( 1 );
  }

  triangle_w = new double[abs(triangle_order)];
  triangle_xy = new double[2 * abs(triangle_order)];

  if ( triangle_order == 1 )
  {
    triangle_o01 ( triangle_w, triangle_xy );
  }
  else if ( triangle_order == 3 )
  {
    triangle_o03 ( triangle_w, triangle_xy );
  }
  else if ( triangle_order == - 3 )
  {
    triangle_o03b ( triangle_w, triangle_xy );
  }
  else if ( triangle_order == 6 )
  {
    triangle_o06 ( triangle_w, triangle_xy );
  }
  else if ( triangle_order == - 6 )
  {
    triangle_o06b ( triangle_w, triangle_xy );
  }
  else if ( triangle_order == 7 )
  {
    triangle_o07 ( triangle_w, triangle_xy );
  }
  else if ( triangle_order == 12 )
  {
    triangle_o12 ( triangle_w, triangle_xy );
  }
  else
  {
    cerr << "\n";
    cerr << "WEDGE_RULE - Fatal error!\n";
    cerr << "  Illegal value of TRIANGLE_ORDER.\n";
    exit ( 1 );
  }

  k = 0;
  for ( i = 0; i < line_order; i++ )
  {
    for ( j = 0; j < abs ( triangle_order ); j++ )
    {
      w[k] = line_w[i] * triangle_w[j];
      xyz[0+k*3] = triangle_xy[0+j*2];
      xyz[1+k*3] = triangle_xy[1+j*2];
      xyz[2+k*3] = line_x[i];
      k = k + 1;
    }
  }

  delete [] line_w;
  delete [] line_x;
  delete [] triangle_w;
  delete [] triangle_xy;

  return;
}
//****************************************************************************80

double wedge_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE_VOLUME: volume of a unit wedge.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double WEDGE_VOLUME, the volume.
//
{
  double value;

  value = 1.0;

  return value;
}
