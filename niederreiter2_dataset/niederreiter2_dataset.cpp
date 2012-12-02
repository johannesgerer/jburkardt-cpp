# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>

# define MAXDEG 50
# define DIM_MAX 20
# define NBITS 31

using namespace std;

int main ( int argc, char *argv[] );
static double RECIP = 1.0 / ( double ) ( 1 << NBITS );
void calcc2 ( int dimen, int cj[DIM_MAX][NBITS] );
void calcv2 ( int maxv, int px_deg, int px[MAXDEG+1], int add[2][2], 
  int mul[2][2], int sub[2][2], int *b_deg, int b[MAXDEG+1], 
  int v[] );
int i4_power ( int i, int j );
void niederreiter2 ( int dim, int *seed, double quasi[] );
double *niederreiter2_generate ( int dim_num, int n, int *seed );
void plymul2 ( int add[2][2], int mul[2][2], int pa_deg, 
  int pa[MAXDEG+1], int pb_deg, int pb[MAXDEG+1], 
  int *pc_deg, int pc[MAXDEG+1] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void setfld2 ( int add[2][2], int mul[2][2], int sub[2][2] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for NIEDERREITER2_DATASET.
//
//  Discussion:
//
//    NIEDERREITER2_DATASET generates a Niederreiter2 dataset and writes it out.
//
//    These program assumes that your computer's word length
//    is 31 bits, excluding sign.  If this is not the case,
//    modify the parameter NBITS throughout accordingly.
//
//  Usage:
//
//    niederreiter2_dataset m n skip
//
//    where
//
//    * M, the spatial dimension,
//    * N, the number of points to generate,
//    * SKIP, number of initial values to skip.
//
//    creates an M by N dataset and writes it to the
//    file "niederreiter2_M_N.txt".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 December 2009
//
//  Author:
//
//    John Burkardt
//
{
  int dim_max = 20;
  int i;
  int m;
  ostringstream m_ostring;
  int n;
  ostringstream n_ostring;
  string output_filename;
  double *r;
  int skip;

  timestamp ( );

  cout << "\n";
  cout << "NIEDERREITER2_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Generate a Niederreiter2 dataset,\n";
  cout << "  of spatial dimension M and BASE = 2.\n";
//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    m = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of M\n";
    cin >> m;
  }

  if ( dim_max < m )
  {
    cout << "\n";
    cout << "NIEDERREITER2_DATASET - Fatal error!\n";
    cout << "  The dimension may not exceed " << dim_max << "\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";
//
//  Get N.
//
//  The sequence length is the number of quasi-random points used to 
//  estimate the integral, excluding warm-up.
//
//  Some users may wish to rewrite the driver to test a [heuristic] 
//  "convergence" criterion, stopping the generation of points
//  when it is passed or when a specified number of points have been 
//  generated, whichever occurs first.
//
  if ( 2 < argc )
  {
    cout << "\n";
    cout << "  Choose the sequence length:\n";
    cout << "\n";
    cout << "  If you do not have a preference, we\n";
    cout << "  suggest using a large power of two, such as:\n";
    cout << "\n";
    cout << "  2^10 = " << i4_power ( 2, 10 ) << "\n";
    cout << "  2^15 = " << i4_power ( 2, 15 ) << "\n";
    cout << "  2^20 = " << i4_power ( 2, 20 ) << "\n";
    cout << "\n";
    cout << "  Enter the sequence length:\n";
    n = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the number of points N\n";
    cin >> n;
  }

  cout << "  Number of points N = " << n << "\n";
//
//  Get SKIP.
//
  if ( 3 < argc )
  {
    cout << "\n";
    cout << "  Choose the number of values to skip:\n";
    cout << "\n";
    cout << "  There is reason to believe that BASE^E,\n";
    cout << "  where E is defined for example in\n";
    cout << "  Bratley, Fox, and Niederreiter [1991],\n";
    cout << "  would be a good choice.  Our formula has\n";
    cout << "  has the form:\n";
    cout << "    SKIP = 2^POWER,\n";
    cout << "  where POWER is chosen so that SKIP comes nowhere\n";
    cout << "  near the largest possible machine-representable\n";
    cout << "  integer ( kind = 4 ).  It does not hurt to choose\n";
    cout << "  POWER larger than E, because skipping is\n";
    cout << "  done implicitly in O(1) time.\n";
    cout << "\n";
    cout << "  The maximum value of E for a fixed dimension\n";
    cout << "  S grows like log S.  We allow some 'fat' for\n";
    cout << "  the implicit constant.\n";
    cout << "\n";
    cout << "  Numerically, 2^POWER = " << i4_power ( 2, 12 ) << "\n";
    cout << "\n";
    cout << "  Enter SKIP (not necessarily the value above)\n";
    skip = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of SKIP\n";
    cin >> skip;
  }

  cout << "  SKIP = " << skip << "\n";
//
//  Compute the data.
//
  r = niederreiter2_generate ( m, n, &skip );
//
//  Write it to a file.
//
  m_ostring << m;
  n_ostring << n;

  output_filename = "niederreiter2_" + m_ostring.str ( ) + "_" 
    + n_ostring.str ( ) + ".txt";

  r8mat_write ( output_filename, m, n, r );

  cout << "\n";
  cout << "  The data was written to the file \"" 
    << output_filename << "\".\n";
//
//  Terminate.
//
  delete [] r;

  cout << "\n";
  cout << "NIEDERREITER2_DATASET:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void calcc2 ( int dim_num, int cj[DIM_MAX][NBITS] )

//****************************************************************************80
//
//  Purpose:
//
//    CALCC2 computes values of the constants C(I,J,R).
//
//  Discussion:
//
//    This program calculates the values of the constants C(I,J,R).
//
//    As far as possible, Niederreiter's notation is used.
//
//    For each value of I, we first calculate all the corresponding
//    values of C.  These are held in the array CI.  All these
//    values are either 0 or 1.  
//
//    Next we pack the values into the
//    array CJ, in such a way that CJ(I,R) holds the values of C
//    for the indicated values of I and R and for every value of
//    J from 1 to NBITS.  The most significant bit of CJ(I,R)
//    (not counting the sign bit) is C(I,1,R) and the least
//    significant bit is C(I,NBITS,R).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    R Lidl, Harald Niederreiter, 
//    Finite Fields,
//    Cambridge University Press, 1984, page 553.
//
//    Harald Niederreiter,
//    Low-discrepancy and low-dispersion sequences,
//    Journal of Number Theory,
//    Volume 30, 1988, pages 51-70.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the sequence to be generated.
//
//    Output, int CJ[DIM_MAX][NBITS], the packed values of 
//    Niederreiter's C(I,J,R)
//
//  Local Parameters:
//
//    Local, int MAXE; we need DIM_MAX irreducible polynomials over Z2.
//    MAXE is the highest degree among these.
//
//    Local, int MAXV, the maximum possible index used in V.
//
{
# define MAXE 6

  int add[2][2];
  int b[MAXDEG+1];
  int b_deg;
  int ci[NBITS][NBITS];
  int count;
  int e;
  int i;
  static int irred[DIM_MAX][MAXE+1] =
  {
    { 0,1,0,0,0,0,0 },
    { 1,1,0,0,0,0,0 },
    { 1,1,1,0,0,0,0 },
    { 1,1,0,1,0,0,0 },
    { 1,0,1,1,0,0,0 },
    { 1,1,0,0,1,0,0 },
    { 1,0,0,1,1,0,0 },
    { 1,1,1,1,1,0,0 },
    { 1,0,1,0,0,1,0 },
    { 1,0,0,1,0,1,0 },
    { 1,1,1,1,0,1,0 },
    { 1,1,1,0,1,1,0 },
    { 1,1,0,1,1,1,0 },
    { 1,0,1,1,1,1,0 },
    { 1,1,0,0,0,0,1 },
    { 1,0,0,1,0,0,1 },
    { 1,1,1,0,1,0,1 },
    { 1,1,0,1,1,0,1 },
    { 1,0,0,0,0,1,1 },
    { 1,1,1,0,0,1,1 }
  };
  int irred_deg[DIM_MAX] = 
    { 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6 };
  int j;
  int maxv = NBITS + MAXE;
  int mul[2][2];
  int nextq[DIM_MAX];
  int p;
  int px[MAXDEG+1];
  int px_deg;
  int q;
  int r;
  int sub[2][2];
  int term;
  int u;
  int v[NBITS+MAXE+1];
//
//  Prepare to work in Z2.
//
  setfld2 ( add, mul, sub );

  for ( i = 0; i < dim_num; i++ )
  {
//
//  For each dimension, we need to calculate powers of an
//  appropriate irreducible polynomial:  see Niederreiter
//  page 65, just below equation (19).
//
//  Copy the appropriate irreducible polynomial into PX,
//  and its degree into E.  Set polynomial B = PX ** 0 = 1.
//  M is the degree of B.  Subsequently B will hold higher
//  powers of PX.
//
    e = irred_deg[i];

    px_deg = irred_deg[i];

    for ( j = 0; j <= px_deg; j++ )
    {
      px[j] = irred[i][j];
    }

    b_deg = 0;
    b[0] = 1;
//
//  Niederreiter (page 56, after equation (7), defines two
//  variables Q and U.  We do not need Q explicitly, but we do need U.
//
    u = 0;

    for ( j = 0; j < NBITS; j++ )
    {
//
//  If U = 0, we need to set B to the next power of PX
//  and recalculate V.  This is done by subroutine CALCV.
//
      if ( u == 0 )
      {
        calcv2 ( maxv, px_deg, px, add, mul, sub, &b_deg, b, v );
      }
//
//  Now C is obtained from V.  Niederreiter obtains A from V (page 65, 
//  near the bottom), and then gets C from A (page 56, equation (7)).  
//  However this can be done in one step.  Here CI(J,R) corresponds to
//  Niederreiter's C(I,J,R).
//
      for ( r = 0; r < NBITS; r++ )
      {
        ci[j][r] = v[r+u];
      }
//
//  Increment U.  
//
//  If U = E, then U = 0 and in Niederreiter's
//  paper Q = Q + 1.  Here, however, Q is not used explicitly.
//
      u = u + 1;
      if ( u == e )
      {
        u = 0;
      }

    }
//
//  The array CI now holds the values of C(I,J,R) for this value
//  of I.  We pack them into array CJ so that CJ(I,R) holds all
//  the values of C(I,J,R) for J from 1 to NBITS.
//
    for ( r = 0; r < NBITS; r++ )
    {
      term = 0;
      for ( j = 0; j < NBITS; j ++ )
      {
        term = 2 * term + ci[j][r];
      }
      cj[i][r] = term;
    }

  }

  return;
# undef MAXE
}
//****************************************************************************80

void calcv2 ( int maxv, int px_deg, int px[MAXDEG+1], int add[2][2], 
  int mul[2][2], int sub[2][2], int *b_deg, int b[MAXDEG+1], 
  int v[] )

//****************************************************************************80
//
//  Purpose:
//
//    CALCV2 calculates the value of the constants V(J,R).
//
//  Discussion:
//
//    This program calculates the values of the constants V(J,R) as
//    described in the reference (BFN) section 3.3.  It is called from CALCC2.  
//
//    Polynomials stored as arrays have the coefficient of degree N 
//    in POLY(N).  
//
//    A polynomial which is identically 0 is given degree -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, pages 494-495, 1994.
//
//  Parameters:
//
//    Input, int MAXV, the dimension of the array V.
//
//    Input, int PX_DEG, the degree of PX.
//
//    Input, int PX[MAXDEG+1], the appropriate irreducible polynomial 
//    for the dimension currently being considered.  
//
//    Input, int ADD[2][2], MUL[2][2], SUB[2][2], the addition, multiplication, 
//    and subtraction tables, mod 2.
//
//    Input/output, int *B_DEG, the degree of the polynomial B.
//
//    Input/output, int B[MAXDEG+1].  On input, B is the polynomial 
//    defined in section 2.3 of BFN.  The degree of B implicitly defines 
//    the parameter J of section 3.3, by degree(B) = E*(J-1).  On output,
//    B has been multiplied by PX, so its degree is now E * J.
//
//    Output, int V[MAXV+1], the computed V array.
//
//  Local Parameters:
//
//    Local, int ARBIT, indicates where the user can place
//    an arbitrary element of the field of order 2.  This means 
//    0 <= ARBIT < 2.  
//
//    Local, int BIGM, is the M used in section 3.3.
//    It differs from the [little] m used in section 2.3,
//    denoted here by M.
//
//    Local, int NONZER, shows where the user must put an arbitrary 
//    non-zero element of the field.  For the code, this means 
//    0 < NONZER < 2.
//
{
  static int arbit = 1;
  int bigm;
  int e;
  int h[MAXDEG+1];
  int h_deg;
  int i;
  int j;
  int kj;
  int m;
  static int nonzer = 1;
  static int p = 2;
  int pb_deg;
  static int q = 2;
  int r;
  int term;
//
  e = px_deg;
//
//  The polynomial H is PX**(J-1), which is the value of B on arrival.
//
//  In section 3.3, the values of Hi are defined with a minus sign:
//  don't forget this if you use them later!
//
  h_deg = *b_deg;

  for ( i = 0; i <= h_deg; i++ )
  {
    h[i] = b[i];
  }

  bigm = h_deg;
//
//  Multiply B by PX so B becomes PX**J.
//  In section 2.3, the values of Bi are defined with a minus sign:
//  don't forget this if you use them later!
//
  pb_deg = *b_deg;

  plymul2 ( add, mul, px_deg, px, pb_deg, b, &pb_deg, b );

  *b_deg = pb_deg;
  m = *b_deg;
//
//  We don't use J explicitly anywhere, but here it is just in case.
//
  j = m / e;
//
//  Now choose a value of Kj as defined in section 3.3.
//  We must have 0 <= Kj < E*J = M.
//  The limit condition on Kj does not seem very relevant
//  in this program.
//
  kj = bigm;
//
//  Choose values of V in accordance with the conditions in section 3.3.
//
  for ( r = 0; r < kj; r++ )
  {
    v[r] = 0;
  }
  v[kj] = 1;

  if ( kj < bigm )
  {
    term = sub [ 0 ] [ h[kj] ];

    for ( r = kj+1; r <= bigm-1; r++ )
    {
      v[r] = arbit;
//
//  Check the condition of section 3.3,
//  remembering that the H's have the opposite sign.
//
      term = sub [ term ] [ mul [ h[r] ] [ v[r] ] ];

    }
//
//  Now V(BIGM) is anything but TERM.
//
    v[bigm] = add [ nonzer] [ term ];

    for ( r = bigm+1; r <= m-1; r++ )
    {
      v[r] = arbit;
    }
  }
  else
  {
    for ( r = kj+1; r <= m-1; r++ )
    {
      v[r] = arbit;
    }

  }
//
//  Calculate the remaining V's using the recursion of section 2.3,
//  remembering that the B's have the opposite sign.
//
  for ( r = 0; r <= maxv - m; r++ )
  {
    term = 0;
    for ( i = 0; i <= m-1; i++ )
    {
      term = sub [ term] [ mul [ b[i] ] [ v[r+i] ] ];
    }
    v[r+m] = term;
  }

  return;
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

void niederreiter2 ( int dim_num, int *seed, double quasi[] )

//****************************************************************************80
//
//  Purpose:
//
//    NIEDERREITER2 returns an element of the Niederreiter sequence base 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Harald Niederreiter,
//    Low-discrepancy and low-dispersion sequences,
//    Journal of Number Theory,
//    Volume 30, 1988, pages 51-70.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the sequence to be generated.
//
//    Input/output, int *SEED, the index of the element entry to
//    compute.  On output, SEED is typically reset by this routine
//    to SEED+1.
//
//    Output, double QUASI[DIM_NUM], the next quasirandom vector.
//
//  Local Parameters:
//
//    Local, int CJ(DIM_MAX,0:NBITS-1), the packed values of 
//    Niederreiter's C(I,J,R).
//
//    Local, int DIM_SAVE, the spatial dimension of the sequence
//    as specified on an initialization call.
//
//    Local, int COUNT, the index of the current item in the sequence,
//    expressed as an array of bits.  COUNT(R) is the same as Niederreiter's
//    AR(N) (page 54) except that N is implicit.
//
//    Local, int NEXTQ[DIM_MAX], the numerators of the next item in the
//    series.  These are like Niederreiter's XI(N) (page 54) except that
//    N is implicit, and the NEXTQ are integers.  To obtain
//    the values of XI(N), multiply by RECIP.
//
{
  static int cj[DIM_MAX][NBITS];
  static int dim_save = 0;
  int gray;
  int i;
  static int nextq[DIM_MAX];
  int r;
  int skip;
  static int seed_save = 0;
//
//  Initialization.
//
  if ( dim_save < 1 || dim_num != dim_save || *seed <= 0 )
  {
    if ( dim_num <= 0 || DIM_MAX < dim_num )
    {
      cout << "\n";
      cout << "NIEDERREITER2 - Fatal error!\n";
      cout << "  Bad spatial dimension.\n";
      exit ( 1 );
    }

    dim_save = dim_num;

    if ( *seed < 0 )
    {
      *seed = 0;
    }

    seed_save = *seed;
//
//  Calculate the C array.
//
    calcc2 ( dim_save, cj );
  }
//
//  Set up NEXTQ appropriately, depending on the Gray code of SEED.
//
//  You can do this every time, starting NEXTQ back at 0,
//  or you can do it once, and then carry the value of NEXTQ
//  around from the previous computation.
//
  if ( *seed != seed_save + 1 )
  {
    gray = ( *seed ) ^ ( *seed / 2 );

    for ( i = 0; i < dim_save; i++ )
    {
      nextq[i] = 0;
    }

    r = 0;

    while ( gray != 0 )
    {
      if ( ( gray % 2 ) != 0 )
      {
        for ( i = 0; i < dim_save; i++ )
        {
          nextq[i] = ( nextq[i] ) ^ ( cj[i][r] );
        }
      }
      gray = gray / 2;
      r = r + 1;
    }
  }
//
//  Multiply the numerators in NEXTQ by RECIP to get the next
//  quasi-random vector.
//
  for ( i = 0; i < dim_save; i++ )
  {
    quasi[i] = ( ( double ) nextq[i] ) * RECIP;
  }
//
//  Find the position of the right-hand zero in SEED.  This
//  is the bit that changes in the Gray-code representation as
//  we go from SEED to SEED+1.
//
  r = 0;
  i = *seed;

  while ( ( i % 2 ) != 0 )
  {
    r = r + 1;
    i = i / 2;
  }
//
//  Check that we have not passed 2**NBITS calls.
//
  if ( NBITS <= r )
  {
    cout << "\n";
    cout << "NIEDERREITER2 - Fatal error!\n";
    cout << "  Too many calls!\n";
    exit ( 1 );
  }
//
//  Compute the new numerators in vector NEXTQ.
//
  for ( i = 0; i < dim_save; i++ )
  {
    nextq[i] = ( nextq[i] ) ^ ( cj[i][r] );
  }

  seed_save = *seed;
  *seed = *seed + 1;

  return;
}
//****************************************************************************80

double *niederreiter2_generate ( int dim_num, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    NIEDERREITER2_GENERATE generates a set of Niederreiter values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points desired.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double R[DIM_NUM*N], the points.
//
{
  int j;
  double *r;

  r = new double[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
    niederreiter2 ( dim_num, seed, r+j*dim_num );
  }

  return r;
}
//****************************************************************************80

void plymul2 ( int add[2][2], int mul[2][2], int pa_deg, 
  int pa[MAXDEG+1], int pb_deg, int pb[MAXDEG+1], 
  int *pc_deg, int pc[MAXDEG+1] )

//****************************************************************************80
//
//  Purpose:
//
//    PLYMUL2 multiplies two polynomials in the field of order 2
//
//  Discussion:
//
//    Polynomials stored as arrays have the coefficient of degree N in 
//    POLY(N), and the degree of the polynomial in POLY(-1).  
//
//    A polynomial which is identically 0 is given degree -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int ADD[2][2], MUL[2][2], 
//    the addition and multiplication tables, mod 2.
//
//    Input, int PA_DEG, the degree of PA.
//
//    Input, int PA[MAXDEG+1], the first polynomial factor.
//
//    Input, int PB_DEG, the degree of PB.
//
//    Input, int PB[MAXDEG+1], the second polynomial factor.
//
//    Output, int *PC_DEG, the degree of the product.
//
//    Output, int PC[MAXDEG+1], the product polynomial.
//
{
  int i;
  int j;
  int jhi;
  int jlo;
  int pt[MAXDEG+1];
  int term;

  if ( pa_deg == -1 || pb_deg == -1 )
  {
    *pc_deg = -1;
  }
  else
  {
    *pc_deg = pa_deg + pb_deg;
  }

  if ( MAXDEG < *pc_deg )
  {
    cout << "\n";
    cout << "PLYMUL2 - Fatal error!\n";
    cout << "  Degree of the product exceeds MAXDEG.\n";
    exit ( 1 );
  }

  for ( i = 0; i <= *pc_deg; i++ )
  {

    jlo = i - pa_deg;
    if ( jlo < 0 ) 
    {
      jlo = 0;
    }

    jhi = pb_deg;
    if ( i < jhi ) 
    {
      jhi = i;
    }

    term = 0;

    for ( j = jlo; j <= jhi; j++ ) 
    {
      term = add [ term ] [ mul [ pa[i-j] ] [ pb[j] ] ];
    }
    pt[i] = term;
  }

  for ( i = 0; i <= *pc_deg; i++ )
  {
    pc[i] = pt[i];
  }

  for ( i = *pc_deg + 1; i <= MAXDEG; i++ )
  {
    pc[i] = 0;
  }

  return;
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
//    Input, double TABLE[M*N], the table data.
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
    return;
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

void setfld2 ( int add[2][2], int mul[2][2], int sub[2][2] )

//****************************************************************************80
//
//  Purpose:
//
//    SETFLD2 sets up arithmetic tables for the finite field of order 2.
//
//  Discussion:
//
//    SETFLD2 sets up addition, multiplication, and subtraction tables 
//    for the finite field of order QIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int ADD[2][2], MUL[2][2], SUB[2][2], the addition, multiplication, 
//    and subtraction tables, mod 2.
//
{
  int i;
  int j;
  int p = 2;
  int q = 2;
//
  for ( i = 0; i < q; i++ )
  {
    for ( j = 0; j < q; j++ )
    {
      add[i][j] = ( i + j ) % p;
      mul[i][j] = ( i * j ) % p;
    }
  }
//
//  Use the addition table to set the subtraction table.
//
  for ( i = 0; i < q; i++ )
  {
    for ( j = 0; j < q; j++ )
    {
      sub[ add[i][j] ] [i] = j;
    }
  }

  return;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
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
