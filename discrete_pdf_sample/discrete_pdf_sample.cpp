# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
double *discrete_cdf_to_xy ( double cdf[20*20], int n, double u[], int *seed );
double *get_discrete_pdf ( );
double *r8mat_copy_new ( int m, int n, double a1[] );
void r8mat_scale ( int m, int n, double s, double a[] );
double r8mat_sum ( int m, int n, double a[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double *r8vec_uniform_01_new ( int n, int *seed );
double *set_discrete_cdf ( double pdf[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DISCRETE_PDF_SAMPLE.
//
//  Discussion:
//
//    This program is an example of how discrete sample or density data
//    can be used to define a PDF (probability density function). 
//
//    In this function and the functions it calls, we assume that we have
//    data for an array of 20 by 20 square subcells of the unit square.
//    We wish to derive a PDF that will allow us to sample an arbitrary
//    number of points from this region.
//
//    In particular, we intend to use the discrete data to generate a PDF
//    which we will then use to generate sample points.
//
//    Roughly speaking, we have kept track of how many fish we caught in
//    each part of a lake, and now we want to simulate catching N fish
//    under the same conditions.
//
//    The statistics for each simulation should be governed by the discrete
//    PDF, but with random variation.  In other words, the actual number
//    of points taken from each subregion is random, and the actual location of
//    each point in a subregion is random, but over many simulations, the
//    statistics of the sample points should reproduce the statistics of
//    the original discrete sample that defined the PDF.
//
//  Usage:
//
//    discrete_pdf_sample n
//
//    where
//
//    * n is the number of sample points desired;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *cdf;
  string filename;
  int n;
  double *pdf;
  int seed;
  double *u;
  double *xy;

  timestamp ( );
  cout << "\n";
  cout << "DISCRETE_PDF_SAMPLE:\n";
  cout << "  Generate sample data using a discrete PDF.\n";
//
//  Get the value of N.
//
  if ( 1 < argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter N, the number of samples to generate: ";
    cin >> n;
  }
//
//  Construct a PDF from the data.
//
  pdf = get_discrete_pdf ( );
//
//  "Integrate" the data over rows and columns of the region to get the CDF.
//
  cdf = set_discrete_cdf ( pdf );
//
//  Choose N CDF values at random.
//
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );
//
//  Find the cell corresponding to each CDF value,
//  and choose a random point in that cell.
//
  xy = discrete_cdf_to_xy ( cdf, n, u, &seed );
//
//  Write data to a file for examination, plotting, or analysis.
//
  filename = "discrete_pdf_sample.txt";
  r8mat_write ( filename, 2, n, xy );

  cout << "\n";
  cout << "  Wrote sample data to file \"" << filename << "\".\n";
//
//  Free memory.
//
  delete [] cdf;
  delete [] pdf;
  delete [] u;
  delete [] xy;
//
//  Terminate.
//
  cout << "\n";
  cout << "DISCRETE_PDF_SAMPLE:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double *discrete_cdf_to_xy ( double cdf[20*20], int n, double u[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_CDF_TO_XY finds XY points corresponding to discrete CDF values.
//
//  Discussion:
//
//    This program is given a discrete CDF function and a set of N random
//    values U.  Each value of U corresponds to a particular (I,J) subregion
//    whose CDF value just exceeds the value of U.  Inside that subregion,
//    we pick a point at random - this is equivalent to assuming the PDF
//    is constant over the subregion.
//
//    This function is part of an example program, for which various
//    assumptions have been made.  In particular, the region is the unit
//    square, and the subregions are formed by a 20 by 20 grid of equal
//    subsquares.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF[20*20], the CDF values associated with each 
//    subcell.  A particular ordering has been given to the subcells so that the
//    CDF is a monotonoe function when the subcells are listed in that order.
//
//    Input, int N, the number of sample points.
//
//    Input, double U[N], N random values.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double DISCRETE_CDF_TO_XY[2*N], the sample points.
//
{
  double high;
  int i;
  int j;
  int k;
  double low;
  double *r;
  double *xy;

  xy = new double[2*n];

  low = 0.0;
  for ( j = 0; j < 20; j++ )
  {
    for ( i = 0; i < 20; i++ )
    {
      high = cdf[i+j*20];
      for ( k = 0; k < n; k++ )
      {
        if ( low <= u[k] && u[k] <= high )
        {
          r = r8vec_uniform_01_new ( 2, seed );
          xy[0+k*2] = ( ( double ) ( i ) + r[0] ) / 20.0;
          xy[1+k*2] = ( ( double ) ( j ) + r[1] ) / 20.0;
          delete [] r;
        }
      }
      low = high;
    }
  }
  return xy;
}
//****************************************************************************80

double *get_discrete_pdf ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_DISCRETE_PDF returns the value of the discrete PDF function in each cell.
//
//  Discussion:
//
//    Cell (I,J) extends from 
//
//      (I-1) * H < Y < I * H
//      (J-1) * H < X < J * H
//
//    We have data for each cell, representing the integral of some PDF
//    over that cell.  The function pdf(x,y) must be nonnegative.  However,
//    we don't impose any other conditions on it.
//
//    The array PDF(:,:) contains the integral of pdf(x,y) over each cell,
//    or, almost as good, simply a sample or average value.
//
//    We load the array PDF, and then we normalize it so that the sum of
//    all the entries is 1.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double GET_DISCRETE_PDF[20*20].  PDF(I,J) is the discrete PDF 
//    for the cell (I,J), normalized so that the sum over all cells is 1.
//
{
  double *pdf;
  double pdf_save[20*20] = {
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
    0.0000, 0.0000, 0.0001, 0.0001, 0.0002, 
    0.0002, 0.0002, 0.0003, 0.0003, 0.0003, 
    0.0003, 0.0003, 0.0002, 0.0002, 0.0002, 
    0.0002, 0.0001, 0.0001, 0.0001, 0.0000, 
    0.0000, 0.0001, 0.0002, 0.0003, 0.0004, 
    0.0004, 0.0005, 0.0006, 0.0006, 0.0006, 
    0.0006, 0.0006, 0.0005, 0.0005, 0.0004, 
    0.0003, 0.0003, 0.0002, 0.0001, 0.0000, 
    0.0000, 0.0002, 0.0003, 0.0005, 0.0006, 
    0.0008, 0.0009, 0.0009, 0.0010, 0.0010, 
    0.0010, 0.0009, 0.0008, 0.0008, 0.0007, 
    0.0006, 0.0004, 0.0003, 0.0002, 0.0000, 
    0.0000, 0.0003, 0.0005, 0.0008, 0.0010, 
    0.0012, 0.0014, 0.0015, 0.0015, 0.0015, 
    0.0015, 0.0014, 0.0013, 0.0011, 0.0010, 
    0.0008, 0.0006, 0.0005, 0.0003, 0.0000, 
    0.0000, 0.0004, 0.0009, 0.0013, 0.0016, 
    0.0019, 0.0021, 0.0023, 0.0023, 0.0023, 
    0.0021, 0.0020, 0.0018, 0.0016, 0.0013, 
    0.0011, 0.0009, 0.0007, 0.0004, 0.0000, 
    0.0000, 0.0007, 0.0014, 0.0020, 0.0025, 
    0.0030, 0.0033, 0.0034, 0.0034, 0.0033, 
    0.0031, 0.0028, 0.0025, 0.0022, 0.0018, 
    0.0015, 0.0012, 0.0009, 0.0006, 0.0000, 
    0.0000, 0.0011, 0.0021, 0.0031, 0.0039, 
    0.0045, 0.0049, 0.0051, 0.0050, 0.0047, 
    0.0043, 0.0039, 0.0034, 0.0029, 0.0024, 
    0.0019, 0.0015, 0.0011, 0.0007, 0.0000, 
    0.0000, 0.0017, 0.0033, 0.0048, 0.0060, 
    0.0069, 0.0074, 0.0074, 0.0072, 0.0066, 
    0.0059, 0.0052, 0.0045, 0.0037, 0.0031, 
    0.0025, 0.0019, 0.0014, 0.0009, 0.0000, 
    0.0000, 0.0025, 0.0050, 0.0073, 0.0091, 
    0.0104, 0.0109, 0.0107, 0.0101, 0.0091, 
    0.0080, 0.0068, 0.0057, 0.0047, 0.0038, 
    0.0030, 0.0023, 0.0017, 0.0011, 0.0000, 
    0.0000, 0.0038, 0.0075, 0.0110, 0.0136, 
    0.0153, 0.0157, 0.0151, 0.0138, 0.0121, 
    0.0104, 0.0087, 0.0071, 0.0058, 0.0046, 
    0.0036, 0.0027, 0.0019, 0.0012, 0.0000, 
    0.0000, 0.0055, 0.0110, 0.0160, 0.0198, 
    0.0218, 0.0219, 0.0205, 0.0182, 0.0155, 
    0.0129, 0.0106, 0.0085, 0.0068, 0.0053, 
    0.0041, 0.0031, 0.0022, 0.0014, 0.0000, 
    0.0000, 0.0077, 0.0154, 0.0224, 0.0276, 
    0.0299, 0.0293, 0.0266, 0.0229, 0.0190, 
    0.0154, 0.0123, 0.0098, 0.0077, 0.0059, 
    0.0045, 0.0034, 0.0024, 0.0015, 0.0000, 
    0.0000, 0.0100, 0.0202, 0.0295, 0.0362, 
    0.0385, 0.0368, 0.0324, 0.0271, 0.0219, 
    0.0174, 0.0137, 0.0107, 0.0082, 0.0063, 
    0.0048, 0.0035, 0.0025, 0.0016, 0.0000, 
    0.0000, 0.0120, 0.0244, 0.0356, 0.0432, 
    0.0455, 0.0426, 0.0366, 0.0298, 0.0236, 
    0.0184, 0.0143, 0.0110, 0.0084, 0.0064, 
    0.0048, 0.0035, 0.0025, 0.0016, 0.0000, 
    0.0000, 0.0134, 0.0266, 0.0382, 0.0461, 
    0.0480, 0.0445, 0.0376, 0.0301, 0.0235, 
    0.0181, 0.0139, 0.0106, 0.0081, 0.0061, 
    0.0046, 0.0033, 0.0023, 0.0015, 0.0000, 
    0.0000, 0.0151, 0.0261, 0.0362, 0.0436, 
    0.0447, 0.0412, 0.0347, 0.0276, 0.0214, 
    0.0164, 0.0125, 0.0095, 0.0072, 0.0054, 
    0.0041, 0.0029, 0.0021, 0.0013, 0.0000, 
    0.0000, 0.0174, 0.0220, 0.0295, 0.0349, 
    0.0361, 0.0333, 0.0281, 0.0225, 0.0175, 
    0.0134, 0.0102, 0.0078, 0.0059, 0.0044, 
    0.0033, 0.0024, 0.0017, 0.0010, 0.0000, 
    0.0000, 0.0097, 0.0152, 0.0200, 0.0235, 
    0.0244, 0.0227, 0.0193, 0.0156, 0.0122, 
    0.0094, 0.0072, 0.0055, 0.0041, 0.0031, 
    0.0023, 0.0017, 0.0012, 0.0007, 0.0000, 
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000 };
  double scale;
  double total;

  pdf = r8mat_copy_new ( 20, 20, pdf_save );
//
//  Normalize to get an integral of 1.
//
  total = r8mat_sum ( 20, 20, pdf );
  scale = 1.0 / total;

  r8mat_scale ( 20, 20, scale, pdf );

  return pdf;
}
//****************************************************************************80

double *r8mat_copy_new ( int m, int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
//
{
  double *a2;
  int i;
  int j;

  a2 = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
//****************************************************************************80

void r8mat_scale ( int m, int n, double s, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SCALE multiplies an R8MAT by a scalar.
//
//  Discussion:
//
//    An R8MAT is an array of R8 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double S, the scale factor.
//
//    Input/output, double A[M*N], the matrix to be scaled.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * s;
    }
  }
  return;
}
//****************************************************************************80

double r8mat_sum ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SUM returns the sum of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], the array.
//
//    Output, double R8MAT_SUM, the sum of the entries.
//
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m];
    }
  }
  return value;
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

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
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
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *set_discrete_cdf ( double pdf[] )

//****************************************************************************80
//
//  Purpose:
//
//    SET_DISCRETE_PDF sets a CDF from a discrete PDF.
//
//  Discussion:
//
//    Here, we proceed from cell (1,1) to (2,1) to (20,1), (1,2), (2,2)...(20,20).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double PDF[20*20], the discrete PDF for the cell (I,J),
//    normalized so that the sum over all cells is 1.
//
//    Output, double CDF[20*20], the discrete CDF for the cell (I,J).
//    CDF(20,20) should be 1.
//
{
  double *cdf;
  int i;
  int j;
  double total;

  cdf = new double[20*20];

  total = 0.0;
  for ( j = 0; j < 20; j++ )
  {
    for ( i = 0; i < 20; i++ )
    {
      total = total + pdf[i+j*20];
      cdf[i+j*20] = total;
    }
  }
  return cdf;
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
