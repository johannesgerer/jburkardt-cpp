# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>

using namespace std;
//
//  GLOBAL DATA.
//
//    The following GLOBAL data, used by many functions,
//    gives the order Q of a field, its characteristic P, and its
//    addition, multiplication, and subtraction tables.
//
//    Global, int DEG_MAX, the maximum degree of the polynomials
//    to be considered.
//
//    Global, int P, the characteristic of the field.
//
//    Global, int Q, the order of the field.
//
//    Global, int Q_MAX, the order of the largest field to
//    be handled.
//
//    Global, int ADD[Q_MAX][Q_MAX], the field addition table. 
//
//    Global, int MUL[Q_MAX][Q_MAX], the field multiplication table. 
//
//    Global, int SUB[Q_MAX][Q_MAX], the field subtraction table.
//
const int DEG_MAX = 50;
int P;
int Q;
const int Q_MAX = 50;

int add[Q_MAX][Q_MAX];
int mul[Q_MAX][Q_MAX];
int sub[Q_MAX][Q_MAX];

int main ( );
int find ( int n, int tab[], int i, int tab_max );
int i4_characteristic ( int q );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void irred ( ofstream &output, int q_init );
int *itop ( int in );
int *plymul ( int pa[], int pb[] );
int ptoi ( int poly[] );
void setfld ( int q_init );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GFPLYS.
//
//  Discussion:
//
//    GFPLYS writes out data about irreducible polynomials.
//
//    The program calculates irreducible polynomials for various
//    finite fields, and writes them out to the file "gfplys.txt".
//
//    Finite field arithmetic is carried out with the help of
//    precalculated addition and multiplication tables found on
//    the file "gfarit.txt".  This file should have been computed
//    and written by the program GFARIT.
//
//    The format of the irreducible polynomials on the output file is
//
//      Q
//      d1   a(1)  a(2) ... a(d1)
//      d2   b(1)  b(2) ... b(d2)
//      ...
//
//    where 
//
//      Q is the order of the field, 
//      d1 is the degree of the first irreducible polynomial, 
//      a(1), a(2), ..., a(d1) are its coefficients.
//
//    Polynomials stored as arrays have the coefficient of degree N in 
//    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
//    DEG is just to remind us of this last fact.  A polynomial which is
//    identically 0 is given degree -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2007
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
{
  char *output_filename = "gfplys.txt";
  ofstream output;
  int q_init;

  timestamp ( );

  cout << "\n";
  cout << "GFPLYS:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  A program to compute a set of irreducible\n";
  cout << "  polynomials over fields of certain orders Q.\n";
  cout << "\n";

  output.open ( output_filename );

  for ( q_init = 2; q_init <= Q_MAX; q_init++ )
  {
    irred ( output, q_init );
  }

  output.close ( );

  cout << "\n";
  cout << "GFPLYS:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int find ( int n, int tab[], int i, int tab_max )

//****************************************************************************80
//
//  Purpose:
//
//    FIND seeks the value N in the range TAB(I) to TAB(TAB_MAX).
//
//  Discussion:
//
//    The vector TAB does not have to be sorted or have any other
//    special properties.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2007
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
//  Parameters:
//
//    Input, int N, the value being sought.
//
//    Input, int TAB[], the table to be searched.
//
//    Input, int I, TAB_MAX, the first and last entries of
//    TAB to be examined.
//
//    Output, int FIND, is the index ( between I and TAB_MAX) of the 
//    entry in TAB that is equal to N, or else -1 if no such value
//    was found.
//
{
  int j;
  int value;
  
  value = -1;

  if ( tab[tab_max-1]  <  n )
  {
    return value;
  }

  for ( j = i; j <= tab_max; j++ )
  {
    if ( tab[j-1] == n )
    {
      value = j;
      return value;
    }
  }

  return value;
}
//****************************************************************************80

int i4_characteristic ( int q )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHARACTERISTIC gives the characteristic for an integer.
//
//  Discussion:
//
//    For any positive integer Q, the characteristic is:
//
//    Q, if Q is a prime;
//    P, if Q = P^N for some prime P and some integer N;
//    0, otherwise, that is, if Q is negative, 0, 1, or the product
//       of more than one distinct prime.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2007
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
//  Parameters:
//
//    Input, int Q, the value to be tested.
//
//    Output, int I4_CHARACTERISTIC, the characteristic of Q.
//
{
  int i;
  int i_max;
  int q_copy;
  int value;
  
  if ( q <= 1 )
  {
    value = 0;
    return value;
  }
//
//  If Q is not prime, then there is at least one prime factor
//  of Q no greater than SQRT(Q)+1.
//
//  A faster code would only consider prime values of I,
//  but that entails storing a table of primes and limiting the 
//  size of Q.  Simplicity and flexibility for now//
//
  i_max = ( int ) ( sqrt ( ( double ) ( q ) ) ) + 1;
  q_copy = q;

  for ( i = 2; i <= i_max; i++ )
  {
    if ( ( q_copy % i ) == 0 )
    {
      while ( ( q_copy % i ) == 0 )
	  {
        q_copy = q_copy / i;
      }

      if ( q_copy == 1 )
	  {
        value = i;
	  }
      else
	  {
        value = 0;
      }
      return value;
    }
  }
//
//  If no factor was found, then Q is prime.
//
  value = q;

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

void irred ( ofstream &output, int q_init )

//****************************************************************************80
//
//  Purpose:
//
//    IRRED computes and writes out a set of irreducible polynomials.
//
//  Discussion:
//
//    We find the irreducible polynomials using a sieve.  
//
//    Polynomials stored as arrays have the coefficient of degree n in 
//    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
//    DEG is just to remind us of this last fact.  A polynomial which is
//    identically 0 is given degree -1.
//
//    Note that the value of NPOL controls the number of polynomials
//    computed, and hence the maximum spatial dimension for the
//    subsequence Niederreiter sequences, JVB, 07 June 2010.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
//  Parameters:
//
//    Input, ofstream &OUTPUT, a reference to the output stream.
//
//    Input, int Q_INIT, the order of the field.
//
//  Local Parameters:
//
//    Local, int SIEVE_MAX, the size of the sieve.  
//
//    Array MONPOL holds monic polynomials.
//
//    Array SIEVE says whether the polynomial is still OK.
//
//    Local, int NPOLS, the number of irreducible polynomials to
//    be calculated for a given field.
//
{
# define SIEVE_MAX 400

  int i;
  int j;
  int k;
  int l;
  int monpol[SIEVE_MAX];
  int n;
  int npols = 50;
  int *pi;
  int *pj;
  int *pk;
  bool sieve[SIEVE_MAX];

  if ( q_init <= 1 || Q_MAX < q_init )
  {
    cerr << "\n";
    cerr << "IRRED - Fatal error!\n";
    cerr << "  Bad value of Q = " << q_init << "\n";
    exit ( 1 );
  }

  P = i4_characteristic ( q_init );
//
//  If no field of order Q_INIT exists, there is nothing to do.
//
  if ( P <= 0 )
  {
    return;
  }

  cout << "  IRRED setting up case for Q = " << q_init << "\n";
//
//  Set up the field arithmetic tables.
//  (Note that SETFLD sets Q = q_init!)
//
  setfld ( q_init );
//
//  Set up the sieve containing only monic polynomials.
//
  i = 0;
  j = 1;
  k = Q;

  for ( n = 1; n <= SIEVE_MAX; n++ )
  {
    i = i + 1;

    if ( i == j )
    {
      i = k;
      j = 2 * k;
      k = Q * k;
    }

    monpol[n-1] = i;
    sieve[n-1] = true;
  }
//
//  Write out the irreducible polynomials as they are found.
//
  n = 0;
  output << setw(3) << Q << "\n";

  for ( i = 1; i <= SIEVE_MAX; i++ )
  {
    if ( sieve[i-1] )
    {
      pi = itop ( monpol[i-1] );
      k = pi[0];
      output << setw(3) << k;
      for ( l = 0; l <= k; l++ )
      {
        output << setw(3) << pi[l+1];
      }
      output << "\n";
      n = n + 1;

      if ( n == npols )
      {
        delete [] pi;
        return;
      }

      for ( j = i; j <= SIEVE_MAX; j++ )
      {
        pj = itop ( monpol[j-1] );
        pk = plymul ( pi, pj );

        k = find ( ptoi ( pk ), monpol, j, SIEVE_MAX );

        if ( k != -1 )
        {
          sieve[k-1] = false;
        }
        delete [] pj;
        delete [] pk;
      }
      delete [] pi;
    }
  }

  cerr << "\n";
  cerr << "IRRED - Warning!\n";
  cerr << "  The sieve size SIEVE_MAX is too small.\n";
  cerr << "  Number of irreducible polynomials found: " << n << "\n";
  cerr << "  Number needed: " << npols << "\n";

  return;
# undef SIEVE_MAX
}
//****************************************************************************80

int *itop ( int in )

//****************************************************************************80
//
//  Purpose:
//
//    ITOP converts an integer to a polynomial in the field of order Q.
//
//  Discussion:
//
//    A nonnegative integer IN can be decomposed into a polynomial in
//    powers of Q, with coefficients between 0 and P-1, by setting:
//
//      J = 0
//      do while ( 0 < IN )
//        POLY(J) = mod ( IN, Q )
//        J = J + 1
//        IN = IN / Q
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2007
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
//  Parameters:
//
//    Input, int IN, the (nonnegative) integer containing the 
//    polynomial information.
//
//    Output, int ITOP[DEG_MAX+2], the polynomial information.
//    ITOP[0] contains the degree of the polynomial.  ITOP[I+1] contains
//    the coefficient of degree I.  Each coefficient is an element of
//    the field of order Q; in other words, each coefficient is
//    between 0 and Q-1.
//
{
  int i;
  int j;
  int *poly;
  
  poly = new int[DEG_MAX+2];

  for ( j = 0; j < DEG_MAX + 2; j++ )
  {
    poly[j] = 0;
  }

  i = in;
  j = -1;

  while ( 0 < i )
  {
    j = j + 1;

    if ( DEG_MAX < j )
	{
      cerr << "\n";
      cerr << "ITOP - Fatal error!\n";
      cerr << "  The polynomial degree exceeds DEG_MAX.\n";
      exit ( 2 );
    }
    poly[j+1] = ( i % Q );
    i = i / Q;
  }

  poly[0] = j;

  return poly;
}
//****************************************************************************80

int *plymul ( int pa[], int pb[] )

//****************************************************************************80
//
//  Purpose:
//
//    PLYMUL multiplies one polynomial by another.
//
//  Discussion:
//
//    Polynomial coefficients are elements of the field of order Q.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2007
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
//  Parameters:
//
//    Input, int PA[DEG_MAX+2], the first polynomial.
//
//    Input, int PB[DEG_MAX+2], the second polynomial.
//
//    Output, int PLYMUL[DEG_MAX+2], the product polynomial.
//
{
  int dega;
  int degb;
  int degc;
  int i;
  int j;
  int p;
  int *pc;
  int term;

  pc = new int[DEG_MAX+2];
  
  dega = pa[0];
  degb = pb[0];

  if ( dega == -1 || degb == -1 )
  {
    degc = -1;
  }
  else
  {
    degc = dega + degb;
  }

  if ( DEG_MAX < degc )
  {
    cerr << "\n";
    cerr << "PLYMUL - Fatal error!\n";
    cerr << "  The degree of the product exceeds DEG_MAX.\n";
    exit ( 1 );
  }

  for ( i = 0; i <= degc; i++ )
  {
    term = 0;
    for ( j = i4_max ( 0, i-dega ); j <= i4_min ( degb, i ); j++ )
	{
      term = add [ term ] [ mul [ pa[i-j+1] ] [ pb[j+1] ] ];
    }
    pc[i+1] = term;
  }

  pc[0] = degc;

  for ( i = degc + 1; i <= DEG_MAX; i++ )
  {
    pc[i+1] = 0;
  }

  return pc;
}
//****************************************************************************80

int ptoi ( int poly[] )

//****************************************************************************80
//
//  Purpose:
//
//    PTOI converts a polynomial in the field of order Q to an integer.
//
//  Discussion:
//
//    A polynomial with coefficients A(*) in the field of order Q
//    can also be stored in an integer I, with
//
//      I = AN*Q**N + ... + A0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2007
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
//  Parameters:
//
//    Input, int POLY[DEG_MAX+2], the polynomial information.
//    POLY[0] contains the degree of the polynomial.  POLY[I] contains
//    the coefficient of degree I-1.  Each coefficient is an element of
//    the field of order Q; in other words, each coefficient is
//    between 0 and Q-1.
//
//    Output, int PTOI, the (nonnegative) integer containing the 
//    polynomial information.
//
{
  int degree;
  int i;
  int j;

  degree = poly[0];
  
  i = 0;
  for ( j = degree; 0 <= j; j-- )
  {
    i = i * Q + poly[j+1];
  }

  return i;
}
//****************************************************************************80

void setfld ( int q_init )

//****************************************************************************80
//
//  Purpose:
//
//    SETFLD sets up arithmetic tables for the finite field of order QIN.
//
//  Discussion:
//
//    This subroutine sets up addition, multiplication, and
//    subtraction tables for the finite field of order QIN.
//
//    If necessary, it reads precalculated tables from the file
//    "gfarit.txt", which are supposed to have been created by GFARIT.
//
//    Polynomials stored as arrays have the coefficient of degree n in 
//    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
//    DEG is just to remind us of this last fact.  A polynomial which is
//    identically 0 is given degree -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2007
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
//    Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, 1994, pages 494-495.
//
//  Parameters:
//
//    Input, Q_INIT, the order of the field.
//
{
  int i;
  ifstream input;
  char *input_filename = "gfarit.txt";
  int j;
  int n;

  if ( q_init <= 1 || Q_MAX < q_init )
  {
    cerr << "\n";
    cerr << "SETFLD - Fatal error!\n";
    cerr << "  Bad value of Q = " << q_init << "\n";
    exit ( 1 );
  }

  Q = q_init;
  P = i4_characteristic ( Q );

  if ( P == 0 )
  {
    cerr << "\n";
    cerr << "SETFLD - Fatal error!\n";
    cerr << "  There is no field of order Q = " << Q << "\n";
    exit ( 1 );
  }
//
//  Handle a field of prime order:  calculate ADD and MUL.
//
  if ( P == Q )
  {
    for ( i = 0; i < Q; i++ )
    {
      for ( j = 0; j < Q; j++ )
      {
        add[i][j] = ( i + j ) % P;
        mul[i][j] = ( i * j ) % P;
      }
    }
  }
//
//  Handle a field of prime-power order:  tables for
//  ADD and MUL are in the file "gfarit.txt".
//
  else
  {
    input.open ( input_filename );

    if ( !input )
    {
      cerr << "\n";
      cerr << "SETFLD - Fatal error!\n";
      cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
      exit ( 1 );
    }

    for ( ; ; )
    {
      input >> n;

      if ( input.eof ( ) )
      {
        cerr << "\n";
        cerr << "SETFLD - Fatal error!\n";
        cerr << "  Could not find tables for Q = " << Q << "\n";
        exit ( 1 );
      }

      for ( i = 0; i < n; i++ )
      {
        for ( j = 0; j < n; j++)
        {
          input >> add[i][j];
        }
      }

      for ( i = 0; i < n; i++ )
      {
        for ( j = 0; j < n; j++)
        {
          input >> mul[i][j];
        }
      }

      if ( n == Q )
      {
        break;
      }
    }
    input.close ( );
  }
//
//  Use the addition table to set the subtraction table.
//
  for ( i = 0; i < Q; i++ )
  {
    for ( j = 0; j < Q; j++ )
    {
      sub[ add[i][j] ] [ i ] = j;
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
