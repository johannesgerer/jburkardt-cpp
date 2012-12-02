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

int main ( void );
void gftab ( ofstream &output, int q_init );
int i4_characteristic ( int q );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *itop ( int in, int p );
int *plyadd ( int pa[], int pb[] );
void plydiv ( int pa[], int pb[], int pq[], int pr[] );
int *plymul ( int pa[], int pb[] );
int ptoi ( int poly[], int q );
void setfld ( int q );
void timestamp ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GFARIT.
//
//  Discussion:
//
//    GFARIT writes the arithmetic tables called "gfarit.txt".
//
//    The program calculates addition and multiplication tables
//    for arithmetic in finite fields, and writes them out to
//    the file "gfarit.txt".  Tables are only calculated for fields
//    of prime-power order Q, the other cases being trivial.
//
//    For each value of Q, the file contains first Q, then the
//    addition table, and lastly the multiplication table.
//
//    After "gfarit.txt" has been set up, run GFPLYS to set up 
//    the file "gfplys.txt".  That operation requires reading 
//    "gfarit.txt".  
//
//    The files "gfarit.txt" and "gfplys.txt" should be saved 
//    for future use.  
//
//    Thus, a user needs to run GFARIT and GFPLYS just once,
//    before running the set of programs associated with GENIN.  
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
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
  char *output_filename = "gfarit.txt";
  ofstream output;
  int q_init;

  timestamp ( );

  cout << "\n";
  cout << "GFARIT:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  A program which computes a set of arithmetic\n";
  cout << "  tables, and writes them to a file.\n";
  cout << "\n";
  cout << "  Tables will be created for fields of prime or prime power order\n";
  cout << "  Q between 2 and " << Q_MAX << ".\n";
  cout << "\n";

  output.open ( output_filename );

  if ( !output )
  {
    cout << "\n";
    cout << "GFARIT - Fatal error!\n";
    cout << "  Could not open the output file: \"" << output_filename << "\"\n";
    return 1;
  }
  
  for ( q_init = 2; q_init <= Q_MAX; q_init++ )
  {
    gftab ( output, q_init );
  }

  output.close ( );

  cout << "\n";
  cout << "GFARIT:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void gftab ( ofstream &output, int q_init )

//****************************************************************************80
//
//  Purpose:
//
//    GFTAB computes and writes data for a particular field size Q_INIT.
//
//  Discussion:
//
//    A polynomial with coefficients A(*) in the field of order Q
//    can also be stored in an integer I, with
//
//      I = AN*Q**N + ... + A0.
//
//    Polynomials stored as arrays have the
//    coefficient of degree n in POLY(N), and the degree of the
//    polynomial in POLY(-1).  The parameter DEG is just to remind
//    us of this last fact.  A polynomial which is identically 0
//    is given degree -1.
//
//    IRRPLY holds irreducible polynomials for constructing
//    prime-power fields.  IRRPLY(-2,I) says which field this
//    row is used for, and then the rest of the row is a
//    polynomial (with the degree in IRRPLY(-1,I) as usual).
//    The chosen irreducible poly is copied into MODPLY for use.
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
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
//    Input, int Q_INIT, the order of the field for which the
//    addition and multiplication tables are needed.
//
{
  int gfadd[Q_MAX][Q_MAX];
  int gfmul[Q_MAX][Q_MAX];
  int i;
  static int irrply[8][8] = {
      {  4, 2, 1, 1, 1, 0, 0, 0 },
	{  8, 3, 1, 1, 0, 1, 0, 0 },
	{  9, 2, 1, 0, 1, 0, 0, 0 },
	{ 16, 4, 1, 1, 0, 0, 1, 0 },
	{ 25, 2, 2, 0, 1, 0, 0, 0 },
	{ 27, 3, 1, 2, 0, 1, 0, 0 },
	{ 32, 5, 1, 0, 1, 0, 0, 1 },
	{ 49, 2, 1, 0, 1, 0, 0, 0 } };
  int j;
  int modply[DEG_MAX+2];
  int *pi;
  int *pj;
  int *pk;
  int *pl;

  if ( q_init <= 1 || Q_MAX < q_init )
  {
    cout << "\n";
    cout << "GFTAB - Fatal error!\n";
    cout << "  Bad value of Q_INIT.\n";
    exit ( 1 );
  }

  P = i4_characteristic ( q_init );
//
//  If QIN is not a prime power, we are not interested.
//
  if ( P == 0 || P == q_init )
  {
    return;
  }

  cout << "  GFTAB computing table for Q = " << q_init
       << "  with characteristic P = " << P << ".\n";
//
//  Otherwise, we set up the elements of the common /FIELD/
//  ready to do arithmetic mod P, the characteristic of Q_INIT.
//
  setfld ( q_init );
//
//  Next find a suitable irreducible polynomial and copy it to array MODPLY.
//
  i = 1;

  while ( irrply[i-1][-2+2] != q_init )
  {
    i = i + 1;
  }

  for ( j = -1; j <= irrply[i-1][-1+2]; j++ )
  {
    modply[j+1] = irrply[i-1][j+2];
  }

  for ( j = irrply[i-1][-1+2]+1; j <= DEG_MAX; j++ )
  {
    modply[j+1] = 0;
  }
//
//  Deal with the trivial cases.
//
  for ( i = 0; i < q_init; i++ )
  {
    gfadd[i][0] = i;
    gfadd[0][i] = i;
    gfmul[i][0] = 0;
    gfmul[0][i] = 0;
  }

  for ( i = 1; i < q_init; i++ )
  {
    gfmul[i][1] = i;
    gfmul[1][i] = i;
  }
//
//  Now deal with the rest.  Each integer from 1 to Q-1
//  is treated as a polynomial with coefficients handled mod P.
//  Multiplication of polynomials is mod MODPLY.
//
  pl = new int[DEG_MAX+2];
 
  for ( i = 1; i < q_init; i++ )
  {
    pi = itop ( i, P );

    for ( j = 1; j <= i; j++ )
    {
      pj = itop ( j, P );
      pk = plyadd ( pi, pj );
      gfadd[i][j] = ptoi ( pk, P );
      gfadd[j][i] = gfadd[i][j];
      delete [] pk;

      if ( 1 < i && 1 < j )
      {
        pk = plymul ( pi, pj );
        plydiv ( pk, modply, pj, pl );
        gfmul[i][j] = ptoi ( pl, P );
        gfmul[j][i] = gfmul[i][j];
        delete [] pk;
      }
      delete [] pj;
    }
    delete [] pi;
  }
  delete [] pl;
//
//  Write out the tables.
//
  output << " " << q_init << "\n";

  for ( i = 0; i < q_init; i++ )
  {
    for ( j = 0; j < q_init; j++ )
    {
      output << " " << gfadd[i][j];
    }
    output << "\n";
  }

  for ( i = 0; i < q_init; i++ )
  {
    for ( j = 0; j < q_init; j++ )
    {
      output << " " << gfmul[i][j];
    }
    output << "\n";
  }

  return;
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
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
//  size of Q.  Simplicity and flexibility for now!
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

int *itop ( int in, int p )

//****************************************************************************80
//
//  Purpose:
//
//    ITOP converts an integer to a polynomial in the field of order P.
//
//  Discussion:
//
//    A nonnegative integer IN can be decomposed into a polynomial in
//    powers of P, with coefficients between 0 and P-1, by setting:
//
//      J = 0
//      do while ( 0 < IN )
//        POLY(J) = mod ( IN, P )
//        J = J + 1
//        IN = IN / P
//      end do
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
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
//    Input, int P, the order of the field.
//
//    Output, int ITOP[DEG_MAX+2], the polynomial information.
//    ITOP[0] contains the degree of the polynomial.  ITOP[I+1] contains
//    the coefficient of degree I.  Each coefficient is an element of
//    the field of order P; in other words, each coefficient is
//    between 0 and P-1.
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
      cout << "\n";
      cout << "ITOP - Fatal error!\n";
      cout << "  The polynomial degree exceeds DEG_MAX.\n";
      exit ( 2 );
    }
    poly[j+1] = ( i % p );
    i = i / p;
  }

  poly[0] = j;

  return poly;
}
//****************************************************************************80

int *plyadd ( int pa[], int pb[] )

//****************************************************************************80
//
//  Purpose:
//
//    PLYADD adds two polynomials.
//
//  Discussion:
//
//    POLY[0] contains the degree of the polynomial.  POLY[I+1] contains
//    the coefficient of degree I.  Each coefficient is an element of
//    the field of order Q; in other words, each coefficient is
//    between 0 and Q-1.
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
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
//    Output, int PLYADD[DEG_MAX+2], the sum polynomial.
//
{
  int degc;
  int i;
  int maxab;
  int *pc;

  pc = new int[DEG_MAX+2];

  maxab = i4_max ( pa[0], pb[0] );

  degc = -1;

  for ( i = 0; i <= maxab; i++ )
  {
    pc[i+1] = add [ pa[i+1] ] [ pb[i+1] ];

    if ( pc[i+1] != 0 )
    {
      degc = i;
    }
  }

  pc[0] = degc;

  for ( i = maxab+1; i <= DEG_MAX; i++ )
  {
    pc[i+1] = 0;
  }

  return pc;
}
//****************************************************************************80

void plydiv ( int pa[], int pb[], int pq[], int pr[] )

//****************************************************************************80
//
//  Purpose:
//
//    PLYDIV divides one polynomial by another.
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
//    06 September 2007
//
//  Author:
//
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
//    Output, int PQ[DEG_MAX+2], the quotient polynomial.
//
//    Output, int PR[DEG_MAX+2], the remainder polynomial.
//
{
  int binv;
  int d;
  int degb;
  int degq;
  int degr;
  int i;
  int j;
  int m;

  if ( pb[0] == -1 )
  {
    cout << "\n";
    cout << "PLYDIV -  Fatal error!\n";
    cout << "  Division by zero polynomial.\n";
    exit ( 1 );
  }

  for ( i = -1; i <= DEG_MAX; i++ )
  {
    pq[i+1] = 0;
    pr[i+1] = pa[i+1];
  }

  degr = pa[0];
  degb = pb[0];
  degq = degr - degb;

  if ( degq < 0 )
  {
    degq = -1;
  }
//
//  Find the inverse of the leading coefficient of PB.
//
  j = pb[degb+1];
  
  for ( i = 1; i <= P - 1; i++ )
  {
    if ( mul[i][j] == 1 )
    {
      binv = i;
    }
  }

  for ( d = degq; 0 <= d; d-- )
  {
    m = mul [ pr[degr+1] ] [ binv ];
	for ( i = degb; 0 <= i; i-- )
	{
      pr[degr+i-degb+1] = sub [ pr[degr+i-degb+1] ] [ mul[m][pb[i+1]] ];
    }
    degr = degr - 1;
    pq[d+1] = m;
  }

  pq[0] = degq;

  while ( pr[degr+1] == 0 && 0 <= degr ) 
  {
    degr = degr - 1;
  }

  pr[0] = degr;

  return;
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
//    06 September 2007
//
//  Author:
//
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
    cout << "\n";
    cout << "PLYMUL - Fatal error!\n";
    cout << "  The degree of the product exceeds DEG_MAX.\n";
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

int ptoi ( int poly[], int q )

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
//    06 September 2007
//
//  Author:
//
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
//    Input, int Q, the order of the field.
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
    i = i * q + poly[j+1];
  }

  return i;
}
//****************************************************************************80

void setfld ( int q_init )

//****************************************************************************80
//
//  Purpose: 
//
//    SETFLD sets up the arithmetic tables for a finite field.
//
//  Discussion:
//
//    This subroutine sets up addition, multiplication, and
//    subtraction tables for the finite field of order QIN.
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
//    06 September 2007
//
//  Author:
//
//    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
//    Input, int Q_INIT, the order of the field.
//
{
  int i;
  ifstream input;
  char *input_filename = "gfarit.txt";
  int j;
  int n;

  if ( q_init <= 1 || Q_MAX < q_init )
  {
    cout << "\n";
    cout << "SETFLD - Fatal error!\n";
    cout << "  Bad value of Q = " << q_init << "\n";
    exit ( 1 );
  }

  Q = q_init;
  P = i4_characteristic ( Q );

  if ( P == 0 )
  {
    cout << "\n";
    cout << "SETFLD - Fatal error!\n";
    cout << "  There is no field of order Q = " << Q << "\n";
    exit ( 1 );
  }
//
//  Set up to handle a field of prime or prime-power order.
//  Calculate the addition and multiplication tables.
//
  for ( i = 0; i < P; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      add[i][j] = ( i + j ) % P;
    }
  }

  for ( i = 0; i < P; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      mul[i][j] = ( i * j ) % P;
    }
  }
//
//  Use the addition table to set the subtraction table.
//
  for ( i = 0; i < P; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      sub [ add[i][j] ] [i] = j;
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
