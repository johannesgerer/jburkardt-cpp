# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "quadrule.hpp"

//****************************************************************************80

void bashforth_set ( int n, double x[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASHFORTH_SET sets abscissas and weights for Adams-Bashforth quadrature.
//
//  Discussion:
//
//    Adams-Bashforth quadrature formulas are normally used in solving
//    ordinary differential equations, and are not really suitable for
//    general quadrature computations.  However, an Adams-Bashforth formula
//    is equivalent to approximating the integral of F(Y(X)) between X(M)
//    and X(M+1), using an explicit formula that relies only on known values
//    of F(Y(X)) at X(M-N+1) through X(M).  For this reason, the formulas
//    have been included here.
//
//    Suppose the unknown function is denoted by Y(X), with derivative
//    F(Y(X)), and that approximate values of the function are known at a
//    series of X values, which we write as X(1), X(2), ..., X(M).  We write
//    the value Y(X(1)) as Y(1) and so on.
//
//    Then the solution of the ODE Y'=F(X,Y) at the next point X(M+1) is
//    computed by:
//
//      Y(M+0] = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
//             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+1-I)) 
//               approximately.
//
//    In the documentation that follows, we replace F(Y(X)) by F(X).
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) F(X) dX.
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) weight[I) * F ( 1 - I ),
//
//    The Adams-Bashforth formulas require equally spaced data.
//
//    Here is how the formula is applied in the case with non-unit spacing:
//
//      Integral ( A <= X <= A+H ) F(X) dX =
//      H * Sum ( 1 <= I <= N ) weight[I) * F ( A - (I-1)*H ),
//      approximately.
//
//    The reference lists the second coefficient of the order 8 Adams-Bashforth
//    formula as
//      weight[1] =  -1162169.0 / 120960.0
//    but this should be
//      weight[1] =  -1152169.0 / 120960.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N should be between 1 and 10, or 12, 14, 16, 18 or 20.
//
//    Output, double X[N], the abscissas.
//
//    Output, double WEIGHT[N], the weights.
//    weight[0] is the weight at X = 0, weight[1] the weight at X = -1,
//    and so on.  The weights are rational, and should sum to 1.  Some
//    weights may be negative.
//
{
  double d;
  int i;

  if ( n == 1 )
  {
    weight[0] = 1.0;
  }
  else if ( n == 2 )
  {
    d = 2.0;

    weight[0] =   3.0 / d;
    weight[1] = - 1.0 / d;
  }
  else if ( n == 3 )
  {
    d = 12.0;

    weight[0] =   23.0 / d;
    weight[1] = - 16.0 / d;
    weight[2] =    5.0 / d;
  }
  else if ( n == 4 )
  {
    d = 24.0;

    weight[0] =   55.0 / d;
    weight[1] = - 59.0 / d;
    weight[2] =   37.0 / d;
    weight[3] =  - 9.0 / d;
  }
  else if ( n == 5 )
  {
    d = 720.0;

    weight[0] =   1901.0 / d;
    weight[1] = - 2774.0 / d;
    weight[2] =   2616.0 / d;
    weight[3] = - 1274.0 / d;
    weight[4] =    251.0 / d;
  }
  else if ( n == 6 )
  {
    d = 1440.0;

    weight[0] =   4277.0 / d;
    weight[1] = - 7923.0 / d;
    weight[2] =   9982.0 / d;
    weight[3] = - 7298.0 / d;
    weight[4] =   2877.0 / d;
    weight[5] =  - 475.0 / d;
  }
  else if ( n == 7 )
  {
    d = 60480.0;

    weight[0] =    198721.0 / d;
    weight[1] =  - 447288.0 / d;
    weight[2] =    705549.0 / d;
    weight[3] =  - 688256.0 / d;
    weight[4] =    407139.0 / d;
    weight[5] =  - 134472.0 / d;
    weight[6] =     19087.0 / d;
  }
  else if ( n == 8 )
  {
    d = 120960.0;

    weight[0] =     434241.0 / d;
    weight[1] =  - 1152169.0 / d;
    weight[2] =    2183877.0 / d;
    weight[3] =  - 2664477.0 / d;
    weight[4] =    2102243.0 / d;
    weight[5] =  - 1041723.0 / d;
    weight[6] =     295767.0 / d;
    weight[7] =    - 36799.0 / d;
  }
  else if ( n == 9 )
  {
    d = 3628800.0;

    weight[0] =   14097247.0 / d;
    weight[1] =  -43125206.0 / d;
    weight[2] =   95476786.0 / d;
    weight[3] = -139855262.0 / d;
    weight[4] =  137968480.0 / d;
    weight[5] =  -91172642.0 / d;
    weight[6] =   38833486.0 / d;
    weight[7] =   -9664106.0 / d;
    weight[8] =    1070017.0 / d;
  }
  else if ( n == 10 )
  {
    d = 7257600.0;

    weight[ 0] =   30277247.0 / d;
    weight[ 1] = -104995189.0 / d;
    weight[ 2] =  265932680.0 / d;
    weight[ 3] = -454661776.0 / d;
    weight[ 4] =  538363838.0 / d;
    weight[ 5] = -444772162.0 / d;
    weight[ 6] =  252618224.0 / d;
    weight[ 7] =  -94307320.0 / d;
    weight[ 8] =   20884811.0 / d;
    weight[ 9] =   -2082753.0 / d;
  }
  else if ( n == 12 )
  {
    d = 958003200.0;

    weight[ 0] =    4527766399.0 / d;
    weight[ 1] =  -19433810163.0 / d;
    weight[ 2] =   61633227185.0 / d;
    weight[ 3] = -135579356757.0 / d;
    weight[ 4] =  214139355366.0 / d;
    weight[ 5] = -247741639374.0 / d;
    weight[ 6] =  211103573298.0 / d;
    weight[ 7] = -131365867290.0 / d;
    weight[ 8] =   58189107627.0 / d;
    weight[ 9] =  -17410248271.0 / d;
    weight[10] =    3158642445.0 / d;
    weight[11] =    -262747265.0 / d;
  }
  else if ( n == 14 )
  {
    d = 5230697472000.0;

    weight[ 0] =    27511554976875.0 / d;
    weight[ 1] =  -140970750679621.0 / d;
    weight[ 2] =   537247052515662.0 / d;
    weight[ 3] = -1445313351681906.0 / d;
    weight[ 4] =  2854429571790805.0 / d;
    weight[ 5] = -4246767353305755.0 / d;
    weight[ 6] =  4825671323488452.0 / d;
    weight[ 7] = -4204551925534524.0 / d;
    weight[ 8] =  2793869602879077.0 / d;
    weight[ 9] = -1393306307155755.0 / d;
    weight[10] =   505586141196430.0 / d;
    weight[11] =  -126174972681906.0 / d;
    weight[12] =    19382853593787.0 / d;
    weight[13] =    -1382741929621.0 / d;
  }
  else if ( n == 16 )
  {
    d = 62768369664000.0;

    weight[ 0] =     362555126427073.0 / d;
    weight[ 1] =   -2161567671248849.0 / d;
    weight[ 2] =    9622096909515337.0 / d;
    weight[ 3] =  -30607373860520569.0 / d;
    weight[ 4] =   72558117072259733.0 / d;
    weight[ 5] = -131963191940828581.0 / d;
    weight[ 6] =  187463140112902893.0 / d;
    weight[ 7] = -210020588912321949.0 / d;
    weight[ 8] =  186087544263596643.0 / d;
    weight[ 9] = -129930094104237331.0 / d;
    weight[10] =   70724351582843483.0 / d;
    weight[11] =  -29417910911251819.0 / d;
    weight[12] =    9038571752734087.0 / d;
    weight[13] =   -1934443196892599.0 / d;
    weight[14] =     257650275915823.0 / d;
    weight[15] =     -16088129229375.0 / d;
  }
  else if ( n == 18 )
  {
    d = 64023737057280000.0;

    weight[ 0] =     401972381695456831.0 / d;
    weight[ 1] =   -2735437642844079789.0 / d;
    weight[ 2] =   13930159965811142228.0 / d;
    weight[ 3] =  -51150187791975812900.0 / d;
    weight[ 4] =  141500575026572531760.0 / d;
    weight[ 5] = -304188128232928718008.0 / d;
    weight[ 6] =  518600355541383671092.0 / d;
    weight[ 7] = -710171024091234303204.0 / d;
    weight[ 8] =  786600875277595877750.0 / d;
    weight[ 9] = -706174326992944287370.0 / d;
    weight[10] =  512538584122114046748.0 / d;
    weight[11] = -298477260353977522892.0 / d;
    weight[12] =  137563142659866897224.0 / d;
    weight[13] =  -49070094880794267600.0 / d;
    weight[14] =   13071639236569712860.0 / d;
    weight[15] =   -2448689255584545196.0 / d;
    weight[16] =     287848942064256339.0 / d;
    weight[17] =     -15980174332775873.0 / d;
  }
  else if ( n == 20 )
  {
    d = 102181884343418880000.0;

    weight[ 0] =      691668239157222107697.0 / d;
    weight[ 1] =    -5292843584961252933125.0 / d;
    weight[ 2] =    30349492858024727686755.0 / d;
    weight[ 3] =  -126346544855927856134295.0 / d;
    weight[ 4] =   399537307669842150996468.0 / d;
    weight[ 5] =  -991168450545135070835076.0 / d;
    weight[ 6] =  1971629028083798845750380.0 / d;
    weight[ 7] = -3191065388846318679544380.0 / d;
    weight[ 8] =  4241614331208149947151790.0 / d;
    weight[ 9] = -4654326468801478894406214.0 / d;
    weight[10] =  4222756879776354065593786.0 / d;
    weight[11] = -3161821089800186539248210.0 / d;
    weight[12] =  1943018818982002395655620.0 / d;
    weight[13] =  -970350191086531368649620.0 / d;
    weight[14] =   387739787034699092364924.0 / d;
    weight[15] =  -121059601023985433003532.0 / d;
    weight[16] =    28462032496476316665705.0 / d;
    weight[17] =    -4740335757093710713245.0 / d;
    weight[18] =      498669220956647866875.0 / d;
    weight[19] =      -24919383499187492303.0 / d;
  }
  else
  {
    cerr << "\n";
    cerr << "BASHFORTH_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1 through 10, 12, 14, 16, 18 or 20.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( -i );
  }

  return;
}
//****************************************************************************80

void bdf_set ( int n, double alpha[], double *beta, double *gamma )

//****************************************************************************80
//
//  Purpose:
//
//    BDF_SET sets weights for backward differentiation ODE weights.
//
//  Discussion:
//
//    GAMMA * Y(N+1) = Sum ( 1 <= I <= N ) ALPHA(I) * Y(N+1-I)
//                     + dX * BETA * Y'(X(N+1),Y(N+1))
//
//    This is equivalent to the backward differentiation corrector formulas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order, between 1 and 6.
//
//    Output, double ALPHA[N], *BETA, *GAMMA, the weights.
//
{
  if ( n == 1 )
  {
    *beta =      1.0;
    *gamma =     1.0;
    alpha[0] = 1.0;
  }
  else if ( n == 2 )
  {
    *beta =        2.0;
    *gamma =       3.0;
    alpha[0] =   4.0;
    alpha[1] = - 1.0;
  }
  else if ( n == 3 )
  {
    *beta =        6.0;
    *gamma =      11.0;
    alpha[0] =  18.0;
    alpha[1] = - 9.0;
    alpha[2] =   2.0;
  }
  else if ( n == 4 )
  {
    *beta =        12.0;
    *gamma =       25.0;
    alpha[0] =   48.0;
    alpha[1] = - 36.0;
    alpha[2] =   16.0;
    alpha[3] =  - 3.0;
  }
  else if ( n == 5 )
  {
    *beta =         60.0;
    *gamma =       137.0;
    alpha[0] =   300.0;
    alpha[1] = - 300.0;
    alpha[2] =   200.0;
    alpha[3] =  - 75.0;
    alpha[4] =    12.0;
  }
  else if ( n == 6 )
  {
    *beta =         60.0;
    *gamma =       147.0;
    alpha[0] =   360.0;
    alpha[1] = - 450.0;
    alpha[2] =   400.0;
    alpha[3] = - 225.0;
    alpha[4] =    72.0;
    alpha[5] =  - 10.0;
  }
  else
  {
    cerr << "\n";
    cerr << "BDF_SET - Fatal error!\n";
    cerr << "  Illegal N requested = " << n << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void bdfc_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    BDFC_SET sets weights for backward differentiation corrector quadrature.
//
//  Discussion:
//
//    A backward differentiation corrector formula is defined for a set
//    of evenly spaced abscissas X(I) with X(0] = 1 and X(1] = 0.  Assuming
//    that the values of the function to be integrated are known at the
//    abscissas, the formula is written in terms of the function value at
//    X(1), and the backward differences at X(1) that approximate the
//    derivatives there.
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * BD^(I-1) F ( 1 ).
//
//    Here, "BD^(I-1) F ( 1 )" denotes the (I-1)st backward difference
//    of F at X = 1, using a spacing of 1.  In particular,
//
//    BD^0 F(0] = F(1)
//    BD^1 F(0] = F(1) - F(0)
//    BD^2 F(0] = F(1) - 2 * F(0) + F(-1 )
//
//    The relationship between a backward difference corrector and the
//    corresponding Adams-Moulton formula may be illustrated for the
//    BDF corrector of order 4:
//
//      BD^0 F(1) - 1/2 * BD^1 F(1) - 1/12 * BD^2 F(1) - 1/24 * BDF^3 F(1)
//      =            F(1)
//        -  1/2 * ( F(1) -         F(0) )
//        - 1/12 * ( F(1) - 2     * F(0) +        F(-1) )
//        - 1/24 * ( F(1) - 3     * F(0) + 3    * F(-1) -        F(-2) )
//      =   9/24 *   F(1) + 19/24 * F(0) - 5/24 * F(-1) + 1/24 * F(-2)
//
//    which is the Adams-Moulton formula of order 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Simeon Fatunla,
//    Numerical Methods for Initial Value Problems in Ordinary Differential
//    Equations,
//    Academic Press, 1988.
//
//  Parameters:
//
//    Input, int N, the order.
//    N can be any value from 1 to 19.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
# define N_MAX 19

  int i;
  static double w_save[N_MAX] = {
                  1.0,
                 -1.0 /               2.0,
                 -1.0 /              12.0,
                 -1.0 /              24.0,
                -19.0 /             720.0,
                 -3.0 /             160.0,
               -863.0 /           60480.0,
               -275.0 /           24792.0,
             -33953.0 /         3628800.0,
              -8183.0 /         1036800.0,
           -3250433.0 /       479001600.0,
              -4671.0 /          788480.0,
       -13695779093.0 /   2615348736000.0,
        -2224234463.0 /    475517952000.0,
      -132282840127.0 /  31384184832000.0,
        -2639651053.0 /    689762304000.0,
    111956703448001.0 /   3201186852864.0,
           50188465.0 /     15613165568.0,
   2334028946344463.0 / 786014494949376.0 };

  if ( N_MAX < n )
  {
    cerr << "\n";
    cerr << "BDFC_SET - Fatal error!\n";
    cerr << "  Input order N = " << n << "\n";
    cerr << "  exceeds maximum order N_MAX = " << N_MAX << "\n";
    exit ( 1 );
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = w_save[i];
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( 1 - i );
  }

  return;
#  undef N_MAX
}
//****************************************************************************80

void bdfp_set ( int order, double x[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    BDFP_SET sets weights for backward differentiation predictor quadrature.
//
//  Discussion:
//
//    A backward differentiation predictor formula is defined for a set
//    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
//    that the values of the function to be integrated are known at the
//    abscissas, the formula is written in terms of the function value at
//    X(2), and the backward differences at X(2) that approximate the
//    derivatives there.  A backward differentiation predictor formula
//    is equivalent to an Adams-Bashforth formula of the same order.
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * BD^(I-1) F ( 0 ),
//
//    Here, "BD^(I-1) F ( 0 )" denotes the (I-1)st backward difference
//    of F at X = 0, using a spacing of 1.  In particular,
//
//    BD^0 F(0) = F(0)
//    BD^1 F(0) = F(0) - F(-1)
//    BD^2 F(0) = F(0) - 2 * F(-1) + F(-2 )
//
//    The relationship between a backward difference predictor and the
//    corresponding Adams-Bashforth formula may be illustrated for the
//    BDF predictor of order 3:
//
//      BD^0 F(0) + 0.5 * BD^1 F(0) + 5/12 * BD^2 F(0)
//      =            F(0)
//        + 1/2  * ( F(0) -         F(1) )
//        + 5/12 * ( F(0) - 2     * F(-1) +      F(-2) )
//      =  23/12 *   F(0) - 16/12 * F(-1) + 5/12 F(-2)
//
//    which is the Adams-Bashforth formula of order 3.
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Simeon Fatunla,
//    Numerical Methods for Initial Value Problems in Ordinary Differential
//    Equations,
//    Academic Press, 1988.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER can be any value from 1 to 19.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weight.
//
{
# define ORDER_MAX 19

  int i;
  double w[ORDER_MAX] = {
                        1.0,
                        1.0 /                2.0,
                        5.0 /               12.0,
                        3.0 /                8.0,
                      251.0 /              720.0,
                       95.0 /              288.0,
                    19087.0 /            60480.0,
                     5257.0 /            17280.0,
                  1070017.0 /          3628800.0,
                    25713.0 /            89600.0,
                 26842253.0 /         95800320.0,
                  4777223.0 /         17418240.0,
             703604254357.0 /    2615348736000.0,
             106364763817.0 /     402361344000.0,
            1166309819657.0 /    4483454976000.0,
                 25221445.0 /         98402304.0,
         8092989203533249.0 /    3201186852864.0,
           85455477715379.0 /      34237292544.0,
     12600467236042756559.0 / 5109094217170944.0 };

  if ( ORDER_MAX < order )
  {
    cerr << "\n";
    cerr << "BDFP_SET - Fatal error!\n";
    cerr << "  Input order ORDER = " << order << "\n";
    cerr << "  exceeds maximum order ORDER_MAX = " << ORDER_MAX << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    weight[i] = w[i];
  }
  for ( i = 0; i < order; i++ )
  {
    x[i] = ( double ) ( - i );
  }

  return;
# undef ORDER_MAX
}
//****************************************************************************80

double bdf_sum ( double func ( double x ), int order, double x[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    BDF_SUM carries out explicit backward difference quadrature on [0,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      RESULT = Sum ( 1 <= I <= ORDER ) WEIGHT(I) * BDF**(I-1) FUNC ( 0 )
//
//    The integral from 0 to 1 is approximated using data at X = 0,
//    -1, -2, ..., -ORDER+1.  This is a form of extrapolation, and
//    the approximation can become poor as ORDER increases.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double X ), the name of the function which 
//    evaluates the integrand.
//
//    Input, int ORDER, the order.
//
//    Input, double X[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights.
//
//    Output, double BDF_SUM, the approximate value of the integral.
//
{
  double *diftab;
  int i;
  int j;
  double result;

  diftab = new double[order];

  for ( i = 0; i < order; i++ )
  {
    diftab[i] = func ( x[i] );
  }

  for ( i = 2; i <= order; i++ )
  {
    for ( j = i; j <= order; j++ )
    {
      diftab[order+i-j-1] = ( diftab[order+i-j-2] - diftab[order+i-j-1] );
    }
  }

  result = r8vec_dot_product ( order, weight, diftab );

  delete [] diftab;

  return result;
}
//****************************************************************************80

char ch_cap ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
//****************************************************************************80

void cheb_set ( int n, double xtab[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEB_SET sets abscissas and weights for Chebyshev quadrature.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( xtab[I) )
//
//    The Chebyshev rule is distinguished by using equal weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Hermann Engels,
//    Numerical Quadrature and Cubature,
//    Academic Press, 1980.
//
//    Zdenek Kopal,
//    Numerical Analysis,
//    John Wiley, 1955.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    N may only have the values 1, 2, 3, 4, 5, 6, 7 or 9.
//    There are NO other Chebyshev rules with real abscissas.
//
//    Output, double XTAB[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;

  if ( n == 1 )
  {
    xtab[0] = 0.0;
  }
  else if ( n == 2 )
  {
    xtab[0] = - 1.0 / sqrt ( 3.0 );
    xtab[1] =   1.0 / sqrt ( 3.0 );
  }
  else if ( n == 3 )
  {
    xtab[0] = - 1.0 / sqrt ( 2.0 );
    xtab[1] =   0.0;
    xtab[2] =   1.0 / sqrt ( 2.0 );
  }
  else if ( n == 4 )
  {
    xtab[0] =   - sqrt ( ( 1.0 + 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab[1] =   - sqrt ( ( 1.0 - 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab[2] =     sqrt ( ( 1.0 - 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab[3] =     sqrt ( ( 1.0 + 2.0/ sqrt ( 5.0 ) ) / 3.0 );
  }
  else if ( n == 5 )
  {
    xtab[0] = - sqrt ( ( 5.0 + sqrt ( 11.0 ) ) / 12.0 );
    xtab[1] = - sqrt ( ( 5.0 - sqrt ( 11.0 ) ) / 12.0 );
    xtab[2] =   0.0;
    xtab[3] =   sqrt ( ( 5.0 - sqrt ( 11.0 ) ) / 12.0 );
    xtab[4] =   sqrt ( ( 5.0 + sqrt ( 11.0 ) ) / 12.0 );
  }
  else if ( n == 6 )
  {
    xtab[0] = - 0.866246818107820591383598;
    xtab[1] = - 0.422518653761111529118546;
    xtab[2] = - 0.266635401516704720331534;
    xtab[3] =   0.266635401516704720331534;
    xtab[4] =   0.422518653761111529118546;
    xtab[5] =   0.866246818107820591383598;
  }
  else if ( n == 7 )
  {
    xtab[0] = - 0.883861700758049035704224;
    xtab[1] = - 0.529656775285156811385048;
    xtab[2] = - 0.323911810519907637519673;
    xtab[3] =   0.0;
    xtab[4] =   0.323911810519907637519673;
    xtab[5] =   0.529656775285156811385048;
    xtab[6] =   0.883861700758049035704224;
  }
  else if ( n == 9 )
  {
    xtab[0] = - 0.911589307728434473664949;
    xtab[1] = - 0.601018655380238071428128;
    xtab[2] = - 0.528761783057879993260181;
    xtab[3] = - 0.167906184214803943068031;
    xtab[4] =   0.0;
    xtab[5] =   0.167906184214803943068031;
    xtab[6] =   0.528761783057879993260181;
    xtab[7] =   0.601018655380238071428128;
    xtab[8] =   0.911589307728434473664949;
  }
  else
  {
    cerr << "\n";
    cerr << "CHEB_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1 through 7, and 9.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( double ) ( n );
  }

  return;
}
//****************************************************************************80

void chebyshev1_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) / sqrt ( 1 - x^2 ) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be greater than 0.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  double pi = 3.141592653589793;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CHEBYSHEV1_COMPUTE - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = pi / ( double ) ( n );
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( pi * ( double ) ( 2 * n - 1 - 2 * i ) 
                    / ( double ) ( 2 * n ) );
  }

  return;
}
//****************************************************************************80

double chebyshev1_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x^2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Output, double CHEBYSHEV1_INTEGRAL, the value of the exact integral.
//
{
  double bot;
  double exact;
  int i;
  double pi = 3.141592653589793;
  double top;
//
//  Get the exact value of the integral.
//
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }
	
    exact = pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;	
  }

  return exact;
}
//****************************************************************************80

void chebyshev2_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X)  sqrt ( 1 - x^2 )  dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be greater than 0.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CHEBYSHEV2_COMPUTE - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    angle = pi * ( double ) ( n - i ) / ( double ) ( n + 1 );
    w[i] = pi / ( double ) ( n + 1 ) * pow ( sin ( angle ), 2 );
    x[i] = cos ( angle );
  }

  return;
}
//****************************************************************************80

double chebyshev2_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x^2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Output, double CHEBYSHEV2_INTEGRAL, the value of the exact integral.
//
{
  double bot;
  double exact;
  int i;
  double pi = 3.141592653589793;
  double top;
//
//  Get the exact value of the integral.
//
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }

	bot = bot * ( double ) ( expon + 2 );

    exact = pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;
  }
  return exact;
}
//****************************************************************************80

void chebyshev3_compute ( int n, double xtab[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV3_COMPUTE computes a Gauss-Chebyshev type 3 quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) f(x) / sqrt ( 1 - x * x ) dx
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( XTAB(I) )
//
//    The N = 1 rule is exceptional.  It consists of a single
//    point at 0, with weight PI.
//
//    For rules with N = 2 or greater, the following remarks apply:
//
//    If N points are used, then Gauss-Chebyshev quadrature
//    will compute the integral exactly, whenever F(X) is a polynomial
//    of degree 2*N-3 or less.
//
//    The abscissas include -1 and 1.
//
//    The first and last weights are 0.5 * PI / ( N - 1),
//    and all other weights are PI / ( N - 1 ).
//
//    If the order is doubled, the abscissas of the new rule include
//    all the points of the old rule.  This fact can be used to
//    efficiently implement error estimation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, double XTAB[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CHEBYSHEV3_COMPUTE - Fatal error!\n";
    cerr << "  N must be at least 1.\n";
    cerr << "  The input value was N = " << n << "\n";
    exit ( 1 );
  }
//
//  Take care of the special case N = 1.
//
  if ( n == 1 )
  {
    xtab[0] = 0.0;
    w[0] = pi;
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( n - 1 - i ) * pi / ( double ) ( n - 1 );
    xtab[i] = cos ( angle );
  }

  w[0] = pi / ( double ) ( 2 * ( n - 1 ) );
  for ( i = 1; i < n - 1; i++ )
  {
    w[i] = pi / ( double ) ( n - 1 );
  }
  w[n-1] = pi / ( double ) ( 2 * ( n - 1 ) );

  return;
}
//****************************************************************************80

void clenshaw_curtis_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
//
//  Discussion:
//
//    This method uses a direct approach.  The paper by Waldvogel
//    exhibits a more efficient approach using Fourier transforms.
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//    The abscissas for the rule of order ORDER can be regarded 
//    as the cosines of equally spaced angles between 180 and 0 degrees:
//
//      X(I) = cos ( ( I - 1 ) * PI / ( ORDER - 1 ) )
//
//    except for the basic case ORDER = 1, when
//
//      X(1) = 0.
//
//    A Clenshaw-Curtis rule that uses ORDER points will integrate
//    exactly all polynomials of degrees 0 through ORDER-1.  If ORDER
//    is odd, then by symmetry the polynomial of degree ORDER will
//    also be integrated exactly.
//
//    If the value of ORDER is increased in a sensible way, then
//    the new set of abscissas will include the old ones.  One such
//    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double b;
  int i;
  int j;
  static double pi = 3.141592653589793;
  double *theta;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CLENSHAW_CURTIS_COMPUTE - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  theta = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    theta[i-1] = ( double ) ( i - 1 ) * pi 
               / ( double ) ( n - 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( theta[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0;

    for ( j = 1; j <= ( ( n - 1 ) / 2 ); j++ )
    {
      if ( 2 * j == ( n - 1 ) )
      {
        b = 1.0;
      }
      else
      {
        b = 2.0;
      }

      w[i] = w[i] - b * cos ( 2.0 * ( double ) ( j ) * theta[i] ) 
           / ( double ) ( 4 * j * j - 1 );
    }
  }

  w[0] = w[0] / ( double ) ( n - 1 );
  for ( i = 2; i <= n-1; i++ )
  {
    w[i-1] = 2.0 * w[i-1] / ( double ) ( n - 1 );
  }
  w[n-1] = w[n-1] / ( double ) ( n - 1 );

  delete [] theta;

  return;
}
//****************************************************************************80

void clenshaw_curtis_set ( int order, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_SET sets a Clenshaw-Curtis quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//    The abscissas for the rule of order ORDER can be regarded 
//    as the cosines of equally spaced angles between 180 and 0 degrees:
//
//      X(I) = cos ( ( I - 1 ) * PI / ( ORDER - 1 ) )
//
//    except for the basic case ORDER = 1, when
//
//      X(1) = 0.
//
//    A Clenshaw-Curtis rule that uses ORDER points will integrate
//    exactly all polynomials of degrees 0 through ORDER-1.  If ORDER
//    is odd, then by symmetry the polynomial of degree ORDER will
//    also be integrated exactly.
//
//    If the value of ORDER is increased in a sensible way, then
//    the new set of abscissas will include the old ones.  One such
//    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...
//    Thus, in the table below, the abscissas for order 9 include
//    those for order 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 17, 33, 65 or 129.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double W[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    x[0] =  0.00000000000000000000;
    w[0] =  2.00000000000000000000;
  }
  else if ( order == 2 )
  {
    x[0] = -1.00000000000000000000;
    x[1] =  1.00000000000000000000;

    w[0] =  1.00000000000000000000;
    w[1] =  1.00000000000000000000;
  }
  else if ( order == 3 )
  {
    x[0] = -1.00000000000000000000;
    x[1] =  0.00000000000000000000;
    x[2] =  1.00000000000000000000;

    w[0] =  0.33333333333333333333;
    w[1] =  1.33333333333333333333;
    w[2] =  0.33333333333333333333;
  }
  else if ( order == 4 )
  {
    x[0] = -1.00000000000000000000;
    x[1] = -0.50000000000000000000;
    x[2] =  0.50000000000000000000;
    x[3] =  1.00000000000000000000;

    w[0] =  0.11111111111111111111;
    w[1] =  0.88888888888888888889;
    w[2] =  0.88888888888888888889;
    w[3] =  0.11111111111111111111;
  }
  else if ( order == 5 )
  {
    x[0] = -1.00000000000000000000;
    x[1] = -0.70710678118654752440;
    x[2] =  0.00000000000000000000;
    x[3] =  0.70710678118654752440;
    x[4] =  1.00000000000000000000;

    w[0] =  0.06666666666666666667;
    w[1] =  0.53333333333333333333;
    w[2] =  0.80000000000000000000;
    w[3] =  0.53333333333333333333;
    w[4] =  0.06666666666666666667;
  }
  else if ( order == 6 )
  {
    x[0] = -1.00000000000000000000;
    x[1] = -0.80901699437494742410;
    x[2] = -0.30901699437494742410;
    x[3] =  0.30901699437494742410;
    x[4] =  0.80901699437493732410;
    x[5] =  1.00000000000000000000;

    w[0] =  0.04000000000000000000;
    w[1] =  0.36074304120001121619;
    w[2] =  0.59925695879998878381;
    w[3] =  0.59925695879998878381;
    w[4] =  0.36074304120001121619;
    w[5] =  0.04000000000000000000;
  }
  else if ( order == 7 )
  {
    x[0] = -1.00000000000000000000;
    x[1] = -0.86602540378443864676;
    x[2] = -0.50000000000000000000;
    x[3] =  0.00000000000000000000;
    x[4] =  0.50000000000000000000;
    x[5] =  0.86602540378443864676;
    x[6] =  1.00000000000000000000;
   
    w[0] =  0.02857142857142857143;
    w[1] =  0.25396825396825396825;
    w[2] =  0.45714285714285714286;
    w[3] =  0.52063492063492063492;
    w[4] =  0.45714285714285714286;
    w[5] =  0.25396825396825396825;
    w[6] =  0.02857142857142857143;
  }
  else if ( order == 8 )
  {
    x[0] = -1.00000000000000000000;
    x[1] = -0.90096886790241912624;
    x[2] = -0.62348980185873353053;
    x[3] = -0.22252093395631440429;
    x[4] =  0.22252093395631440429;
    x[5] =  0.62348980185873353053;
    x[6] =  0.90096886790241910624;
    x[7] =  1.00000000000000000000;
   
    w[0] =  0.02040816326530612245;
    w[1] =  0.19014100721820835178;
    w[2] =  0.35224242371815911533;
    w[3] =  0.43720840579832641044;
    w[4] =  0.43720840579832641044;
    w[5] =  0.35224242371815911533;
    w[6] =  0.19014100721820835178;
    w[7] =  0.02040816326530612245;
  }
  else if ( order == 9 )
  {
    x[0] = -1.00000000000000000000;
    x[1] = -0.92387953251128675613;
    x[2] = -0.70710678118654752440;
    x[3] = -0.38268343236508977173;
    x[4] =  0.00000000000000000000;
    x[5] =  0.38268343236508977173;
    x[6] =  0.70710678118654752440;
    x[7] =  0.92387953251128675613;
    x[8] =  1.00000000000000000000;

    w[0] =  0.01587301587301587302;
    w[1] =  0.14621864921601815501;
    w[2] =  0.27936507936507936508;
    w[3] =  0.36171785872048978150;
    w[4] =  0.39365079365079365079;
    w[5] =  0.36171785872048978150;
    w[6] =  0.27936507936507936508;
    w[7] =  0.14621864921601815501;
    w[8] =  0.01587301587301587302;
  }
  else if ( order == 10 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.93969262078590838405;
    x[2]  = -0.76604444311897903520;
    x[3]  = -0.50000000000000000000;
    x[4]  = -0.17364817766693034885;
    x[5]  =  0.17364817766693034885;
    x[6]  =  0.50000000000000000000;
    x[7]  =  0.76604444311897903520;
    x[8]  =  0.93969262078590838405;
    x[9] =  1.00000000000000000000;

    w[0]  =  0.01234567901234567901;
    w[1]  =  0.11656745657203712296;
    w[2]  =  0.22528432333810440813;
    w[3]  =  0.30194003527336860670;
    w[4]  =  0.34386250580414418320;
    w[5]  =  0.34386250580414418320;
    w[6]  =  0.30194003527336860670;
    w[7]  =  0.22528432333810440813;
    w[8]  =  0.11656745657203712296;
    w[9]  =  0.01234567901234567901;
  }
  else if ( order == 11 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.95105651629515357212;
    x[2]  = -0.80901699437494742410;
    x[3]  = -0.58778525229247312917;
    x[4]  = -0.30901699437494742410;
    x[5]  =  0.00000000000000000000;
    x[6]  =  0.30901699437494742410;
    x[7]  =  0.58778525229247312917;
    x[8]  =  0.80901699437494742410;
    x[9]  =  0.95105651629515357212;
    x[10] =  1.00000000000000000000;

    w[0]  =  0.01010101010101010101;
    w[1]  =  0.09457905488370156116;
    w[2]  =  0.18563521442424776529;
    w[3]  =  0.25358833328368660623;
    w[4]  =  0.29921327042423708320;
    w[5]  =  0.31376623376623376623;
    w[6]  =  0.29921327042423708320;
    w[7]  =  0.25358833328368660623;
    w[8]  =  0.18563521442424776529;
    w[9]  =  0.09457905488370156116;
    w[10] =  0.01010101010101010101;
  }
  else if ( order == 12 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.95949297361449738989;
    x[2]  = -0.84125353283118116886;
    x[3]  = -0.65486073394528506406;
    x[4]  = -0.41541501300188642553;
    x[5]  = -0.14231483827328514044;
    x[6]  =  0.14231483827328514044;
    x[7]  =  0.41541501300188642553;
    x[8]  =  0.65486073394528506406;
    x[9]  =  0.84125353283118116886;
    x[10] =  0.95949297361449738989;
    x[11] =  1.00000000000000000000;

    w[0]  =  0.00826446280991735537;
    w[1]  =  0.07856015374620000543;
    w[2]  =  0.15504045508256136552;
    w[3]  =  0.21556254600086858099;
    w[4]  =  0.25991734106691617602;
    w[5]  =  0.28265504129353651666;
    w[6]  =  0.28265504129353651666;
    w[7]  =  0.25991734106691617602;
    w[8]  =  0.21556254600086858099;
    w[9]  =  0.15504045508256136552;
    w[10] =  0.07856015374620000543;
    w[11] =  0.00826446280991735537;
  }
  else if ( order == 13 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.96592582628906828675;
    x[2]  = -0.86602540378443864676;
    x[3]  = -0.70710678118654752440;
    x[4]  = -0.50000000000000000000;
    x[5]  = -0.25881904510252076235;
    x[6]  =  0.00000000000000000000;
    x[7]  =  0.25881904510252076235;
    x[8]  =  0.50000000000000000000;
    x[9]  =  0.70710678118654752440;
    x[10] =  0.86602540378443864676;
    x[11] =  0.96592582628906828675;
    x[12] =  1.00000000000000000000;

    w[0]  =  0.00699300699300699301;
    w[1]  =  0.06605742495207439452;
    w[2]  =  0.13154253154253154253;
    w[3]  =  0.18476338476338476338;
    w[4]  =  0.22697302697302697303;
    w[5]  =  0.25267569378104433860;
    w[6]  =  0.26198986198986198986;
    w[7]  =  0.25267569378104433860;
    w[8]  =  0.22697302697302697303;
    w[9]  =  0.18476338476338476338;
    w[10] =  0.13154253154253154253;
    w[11] =  0.06605742495207439452;
    w[12] =  0.00699300699300699301;
  }
  else if ( order == 14 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.97094181742605202716;
    x[2]  = -0.88545602565320989590;
    x[3]  = -0.74851074817110109863;
    x[4]  = -0.56806474673115580251;
    x[5]  = -0.35460488704253562597;
    x[6]  = -0.12053668025532305335;
    x[7]  =  0.12053668025532305335;
    x[8]  =  0.35460488704253562597;
    x[9]  =  0.56806474673115580251;
    x[10] =  0.74851074817110109863;
    x[11] =  0.88545602565320989590;
    x[12] =  0.97094181742605202716;
    x[13] =  1.00000000000000000000;

    w[0]  =  0.00591715976331360947;
    w[1]  =  0.05646531376341444627;
    w[2]  =  0.11276867248985655881;
    w[3]  =  0.16003802611671868523;
    w[4]  =  0.19899241036578321848;
    w[5]  =  0.22590304977856444935;
    w[6]  =  0.23991536772234903239;
    w[7]  =  0.23991536772234903239;
    w[8]  =  0.22590304977856444935;
    w[9]  =  0.19899241036578321848;
    w[10] =  0.16003802611671868523;
    w[11] =  0.11276867248985655881;
    w[12] =  0.05646531376341444627;
    w[13] =  0.00591715976331360947;
  }
  else if ( order == 15 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.97492791218182360702;
    x[2]  = -0.90096886790241912624;
    x[3]  = -0.78183148246802980871;
    x[4]  = -0.62348980185873353053;
    x[5]  = -0.43388373911755812048;
    x[6]  = -0.22252093395631440429;
    x[7]  =  0.00000000000000000000;
    x[8]  =  0.22252093395631440429;
    x[9]  =  0.43388373911755812048;
    x[10] =  0.62348980185873353053;
    x[11] =  0.78183148246802980871;
    x[12] =  0.90096886790241912624;
    x[13] =  0.97492791218182360702;
    x[14] =  1.00000000000000000000;

    w[0]  =  0.00512820512820512821;
    w[1]  =  0.04869938729508823855;
    w[2]  =  0.09782039167605215913;
    w[3]  =  0.13966507849560431803;
    w[4]  =  0.17560578900106674677;
    w[5]  =  0.20205146748238357364;
    w[6]  =  0.21888151163057340180;
    w[7]  =  0.22429633858205286777;
    w[8]  =  0.21888151163057340180;
    w[9]  =  0.20205146748238357364;
    w[10] =  0.17560578900106674677;
    w[11] =  0.13966507849560431803;
    w[12] =  0.09782039167605215913;
    w[13] =  0.04869938729508823855;
    w[14] =  0.00512820512820512821;
  }
  else if ( order == 16 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.97814760073380563793;
    x[2]  = -0.91354545764260089550;
    x[3]  = -0.80901699437494742410;
    x[4]  = -0.66913060635885821383;
    x[5]  = -0.50000000000000000000;
    x[6]  = -0.30901699437494742410;
    x[7]  = -0.10452846326765347140;
    x[8]  =  0.10452846326765347140;
    x[9]  =  0.30901699437494742410;
    x[10] =  0.50000000000000000000;
    x[11] =  0.66913060635885821383;
    x[12] =  0.80901699437494742410;
    x[13] =  0.91354545764260089550;
    x[14] =  0.97814760073380563793;
    x[15] =  1.00000000000000000000;

    w[0]  =  0.00444444444444444444;
    w[1]  =  0.04251476624752508988;
    w[2]  =  0.08553884025933288291;
    w[3]  =  0.12294010082849361533;
    w[4]  =  0.15573317603967369176;
    w[5]  =  0.18132978132978132978;
    w[6]  =  0.19921478132638853955;
    w[7]  =  0.20828410952436040635;
    w[8]  =  0.20828410952436040635;
    w[9]  =  0.19921478132638853955;
    w[10] =  0.18132978132978132978;
    w[11] =  0.15573317603967369176;
    w[12] =  0.12294010082849361533;
    w[13] =  0.08553884025933288291;
    w[14] =  0.04251476624752508988;
    w[15] =  0.00444444444444444444;
  }
  else if ( order == 17 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.98078528040323044913;
    x[2]  = -0.92387953251128675613;
    x[3]  = -0.83146961230254523708;
    x[4]  = -0.70710678118654752440;
    x[5]  = -0.55557023301960222474;
    x[6]  = -0.38268343236508977173;
    x[7]  = -0.19509032201612826785;
    x[8]  =  0.00000000000000000000;
    x[9]  =  0.19509032201612826785;
    x[10] =  0.38268343236508977173;
    x[11] =  0.55557023301960222474;
    x[12] =  0.70710678118654752440;
    x[13] =  0.83146961230254523708;
    x[14] =  0.92387953251128675613;
    x[15] =  0.98078528040323044913;
    x[16] =  1.00000000000000000000;

    w[0]  =  0.00392156862745098039;
    w[1]  =  0.03736870283720561032;
    w[2]  =  0.07548233154315183441;
    w[3]  =  0.10890555258189093044;
    w[4]  =  0.13895646836823307412;
    w[5]  =  0.16317266428170330256;
    w[6]  =  0.18147378423649335700;
    w[7]  =  0.19251386461292564687;
    w[8]  =  0.19641012582189052777;
    w[9]  =  0.19251386461292564687;
    w[10] =  0.18147378423649335700;
    w[11] =  0.16317266428170330256;
    w[12] =  0.13895646836823307412;
    w[13] =  0.10890555258189093044;
    w[14] =  0.07548233154315183441;
    w[15] =  0.03736870283720561032;
    w[16] =  0.00392156862745098039;
  }
  else if ( order == 33 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.99518472667219688624;
    x[2]  = -0.98078528040323044913;
    x[3]  = -0.95694033573220886494;
    x[4]  = -0.92387953251128675613;
    x[5]  = -0.88192126434835502971;
    x[6]  = -0.83146961230254523708;
    x[7]  = -0.77301045336273696081;
    x[8]  = -0.70710678118654752440;
    x[9]  = -0.63439328416364549822;
    x[10] = -0.55557023301960222474;
    x[11] = -0.47139673682599764856;
    x[12] = -0.38268343236508977173;
    x[13] = -0.29028467725446236764;
    x[14] = -0.19509032201612826785;
    x[15] = -0.098017140329560601994;
    x[16] =  0.000000000000000000000;
    x[17] =  0.098017140329560601994;
    x[18] =  0.19509032201612826785;
    x[19] =  0.29028467725446236764;
    x[20] =  0.38268343236508977173;
    x[21] =  0.47139673682599764856;
    x[22] =  0.55557023301960222474;
    x[23] =  0.63439328416364549822;
    x[24] =  0.70710678118654752440;
    x[25] =  0.77301045336273696081;
    x[26] =  0.83146961230254523708;
    x[27] =  0.88192126434835502971;
    x[28] =  0.92387953251128675613;
    x[29] =  0.95694033573220886494;
    x[30] =  0.98078528040323044913;
    x[31] =  0.99518472667219688624;
    x[32] =  1.00000000000000000000;

    w[0]  =  0.00097751710654936461;
    w[1]  =  0.00939319796295501470;
    w[2]  =  0.01923424513268114918;
    w[3]  =  0.02845791667723369009;
    w[4]  =  0.03759434191404720602;
    w[5]  =  0.04626276283775174949;
    w[6]  =  0.05455501630398031044;
    w[7]  =  0.06227210954529400455;
    w[8]  =  0.06942757563043545090;
    w[9]  =  0.07588380044138847048;
    w[10] =  0.08163481765493851023;
    w[11] =  0.08657753844182743544;
    w[12] =  0.09070611286772099874;
    w[13] =  0.09394324443876873573;
    w[14] =  0.09629232594548817919;
    w[15] =  0.09769818820805558182;
    w[16] =  0.09817857778176829677;
    w[17] =  0.09769818820805558182;
    w[18] =  0.09629232594548817919;
    w[19] =  0.09394324443876873573;
    w[20] =  0.09070611286772099874;
    w[21] =  0.08657753844182743544;
    w[22] =  0.08163481765493851023;
    w[23] =  0.07588380044138847048;
    w[24] =  0.06942757563043545090;
    w[25] =  0.06227210954529400455;
    w[26] =  0.05455501630398031044;
    w[27] =  0.04626276283775174949;
    w[28] =  0.03759434191404720602;
    w[29] =  0.02845791667723369009;
    w[30] =  0.01923424513268114918;
    w[31] =  0.00939319796295501470;
    w[32] =  0.00097751710654936461;
  }
  else if ( order == 65 )
  {
    x[0]  = -1.00000000000000000000;
    x[1]  = -0.99879545620517239271;
    x[2]  = -0.99518472667219688624;
    x[3]  = -0.98917650996478097345;
    x[4]  = -0.98078528040323044913;
    x[5]  = -0.97003125319454399260;
    x[6]  = -0.95694033573220886494;
    x[7]  = -0.94154406518302077841;
    x[8]  = -0.92387953251128675613;
    x[9]  = -0.90398929312344333159;
    x[10] = -0.88192126434835502971;
    x[11] = -0.85772861000027206990;
    x[12] = -0.83146961230254523708;
    x[13] = -0.80320753148064490981;
    x[14] = -0.77301045336273696081;
    x[15] = -0.74095112535495909118;
    x[16] = -0.70710678118654752440;
    x[17] = -0.67155895484701840063;
    x[18] = -0.63439328416364549822;
    x[19] = -0.59569930449243334347;
    x[20] = -0.55557023301960222474;
    x[21] = -0.51410274419322172659;
    x[22] = -0.47139673682599764856;
    x[23] = -0.42755509343028209432;
    x[24] = -0.38268343236508977173;
    x[25] = -0.33688985339222005069;
    x[26] = -0.29028467725446236764;
    x[27] = -0.24298017990326388995;
    x[28] = -0.19509032201612826785;
    x[29] = -0.14673047445536175166;
    x[30] = -0.098017140329560601994;
    x[31] = -0.049067674327418014255;
    x[32] =  0.000000000000000000000;
    x[33] =  0.049067674327418014255;
    x[34] =  0.098017140329560601994;
    x[35] =  0.14673047445536175166;
    x[36] =  0.19509032201612826785;
    x[37] =  0.24298017990326388995;
    x[38] =  0.29028467725446236764;
    x[39] =  0.33688985339222005069;
    x[40] =  0.38268343236508977173;
    x[41] =  0.42755509343028209432;
    x[42] =  0.47139673682599764856;
    x[43] =  0.51410274419322172659;
    x[44] =  0.55557023301960222474;
    x[45] =  0.59569930449243334347;
    x[46] =  0.63439328416364549822;
    x[47] =  0.67155895484701840063;
    x[48] =  0.70710678118654752440;
    x[49] =  0.74095112535495909118;
    x[50] =  0.77301045336273696081;
    x[51] =  0.80320753148064490981;
    x[52] =  0.83146961230254523708;
    x[53] =  0.85772861000027206990;
    x[54] =  0.88192126434835502971;
    x[55] =  0.90398929312344333159;
    x[56] =  0.92387953251128675613;
    x[57] =  0.94154406518302077841;
    x[58] =  0.95694033573220886494;
    x[59] =  0.97003125319454399260;
    x[60] =  0.98078528040323044913;
    x[61] =  0.98917650996478097345;
    x[62] =  0.99518472667219688624;
    x[63] =  0.99879545620517239271;
    x[64] =  1.00000000000000000000;

    w[0]  =  0.00024420024420024420;
    w[1]  =  0.00235149067531170332;
    w[2]  =  0.00483146544879091264;
    w[3]  =  0.00719269316173611402;
    w[4]  =  0.00958233879528379039;
    w[5]  =  0.01192339471421277160;
    w[6]  =  0.01425206043235199679;
    w[7]  =  0.01653498765728958965;
    w[8]  =  0.01878652974179578354;
    w[9]  =  0.02098627442973743378;
    w[10] =  0.02314069493435819848;
    w[11] =  0.02523506498175476590;
    w[12] =  0.02727225714146838686;
    w[13] =  0.02924065319746833770;
    w[14] =  0.03114129710406762447;
    w[15] =  0.03296454656997632997;
    w[16] =  0.03471049818092511427;
    w[17] =  0.03637092028663918309;
    w[18] =  0.03794545992128481711;
    w[19] =  0.03942698871295609976;
    w[20] =  0.04081501340035783384;
    w[21] =  0.04210333111141810203;
    w[22] =  0.04329151496169082935;
    w[23] =  0.04437417923925731580;
    w[24] =  0.04535110955166067221;
    w[25] =  0.04621766751092557684;
    w[26] =  0.04697395904661414870;
    w[27] =  0.04761604458525019296;
    w[28] =  0.04814443257251220341;
    w[29] =  0.04855584485714105274;
    w[30] =  0.04885125664306609371;
    w[31] =  0.04902801843102555294;
    w[32] =  0.04908762351494245585;
    w[33] =  0.04902801843102555294;
    w[34] =  0.04885125664306609371;
    w[35] =  0.04855584485714105274;
    w[36] =  0.04814443257251220341;
    w[37] =  0.04761604458525019296;
    w[38] =  0.04697395904661414870;
    w[39] =  0.04621766751092557684;
    w[40] =  0.04535110955166067221;
    w[41] =  0.04437417923925731580;
    w[42] =  0.04329151496169082935;
    w[43] =  0.04210333111141810203;
    w[44] =  0.04081501340035783384;
    w[45] =  0.03942698871295609976;
    w[46] =  0.03794545992128481711;
    w[47] =  0.03637092028663918309;
    w[48] =  0.03471049818092511427;
    w[49] =  0.03296454656997632997;
    w[50] =  0.03114129710406762447;
    w[51] =  0.02924065319746833770;
    w[52] =  0.02727225714146838686;
    w[53] =  0.02523506498175476590;
    w[54] =  0.02314069493435819848;
    w[55] =  0.02098627442973743378;
    w[56] =  0.01878652974179578354;
    w[57] =  0.01653498765728958965;
    w[58] =  0.01425206043235199679;
    w[59] =  0.01192339471421277160;
    w[60] =  0.00958233879528379039;
    w[61] =  0.00719269316173611402;
    w[62] =  0.00483146544879091264;
    w[63] =  0.00235149067531170332;
    w[64] =  0.00024420024420024420;
  }
  else if ( order == 129 )
  {
    x[0]   = -1.00000000000000000000;
    x[1]   = -0.99969881869620422012;
    x[2]   = -0.99879545620517239271;
    x[3]   = -0.99729045667869021614;
    x[4]   = -0.99518472667219688624;
    x[5]   = -0.99247953459870999816;
    x[6]   = -0.98917650996478097345;
    x[7]   = -0.98527764238894124477;
    x[8]   = -0.98078528040323044913;
    x[9]   = -0.97570213003852854446;
    x[10]  = -0.97003125319454399260;
    x[11]  = -0.96377606579543986669;
    x[12]  = -0.95694033573220886494;
    x[13]  = -0.94952818059303666720;
    x[14]  = -0.94154406518302077841;
    x[15]  = -0.93299279883473888771;
    x[16]  = -0.92387953251128675613;
    x[17]  = -0.91420975570353065464;
    x[18]  = -0.90398929312344333159;
    x[19]  = -0.89322430119551532034;
    x[20]  = -0.88192126434835502971;
    x[21]  = -0.87008699110871141865;
    x[22]  = -0.85772861000027206990;
    x[23]  = -0.84485356524970707326;
    x[24]  = -0.83146961230254523708;
    x[25]  = -0.81758481315158369650;
    x[26]  = -0.80320753148064490981;
    x[27]  = -0.78834642762660626201;
    x[28]  = -0.77301045336273696081;
    x[29]  = -0.75720884650648454758;
    x[30]  = -0.74095112535495909118;
    x[31]  = -0.72424708295146692094;
    x[32]  = -0.70710678118654752440;
    x[33]  = -0.68954054473706692462;
    x[34]  = -0.67155895484701840063;
    x[35]  = -0.65317284295377676408;
    x[36]  = -0.63439328416364549822;
    x[37]  = -0.61523159058062684548;
    x[38]  = -0.59569930449243334347;
    x[39]  = -0.57580819141784530075;
    x[40]  = -0.55557023301960222474;
    x[41]  = -0.53499761988709721066;
    x[42]  = -0.51410274419322172659;
    x[43]  = -0.49289819222978403687;
    x[44]  = -0.47139673682599764856;
    x[45]  = -0.44961132965460660005;
    x[46]  = -0.42755509343028209432;
    x[47]  = -0.40524131400498987091;
    x[48]  = -0.38268343236508977173;
    x[49]  = -0.35989503653498814878;
    x[50]  = -0.33688985339222005069;
    x[51]  = -0.31368174039889147666;
    x[52]  = -0.29028467725446236764;
    x[53]  = -0.26671275747489838633;
    x[54]  = -0.24298017990326388995;
    x[55]  = -0.21910124015686979723;
    x[56]  = -0.19509032201612826785;
    x[57]  = -0.17096188876030122636;
    x[58]  = -0.14673047445536175166;
    x[59]  = -0.12241067519921619850;
    x[60]  = -0.098017140329560601994;
    x[61]  = -0.073564563599667423529;
    x[62]  = -0.049067674327418014255;
    x[63]  = -0.024541228522912288032;
    x[64]  =  0.00000000000000000000;
    x[65]  =  0.024541228522912288032;
    x[66]  =  0.049067674327418014255;
    x[67]  =  0.073564563599667423529;
    x[68]  =  0.098017140329560601994;
    x[69]  =  0.12241067519921619850;
    x[70]  =  0.14673047445536175166;
    x[71]  =  0.17096188876030122636;
    x[72]  =  0.19509032201612826785;
    x[73]  =  0.21910124015686979723;
    x[74]  =  0.24298017990326388995;
    x[75]  =  0.26671275747489838633;
    x[76]  =  0.29028467725446236764;
    x[77]  =  0.31368174039889147666;
    x[78]  =  0.33688985339222005069;
    x[79]  =  0.35989503653498814878;
    x[80]  =  0.38268343236508977173;
    x[81]  =  0.40524131400498987091;
    x[82]  =  0.42755509343028209432;
    x[83]  =  0.44961132965460660005;
    x[84]  =  0.47139673682599764856;
    x[85]  =  0.49289819222978403687;
    x[86]  =  0.51410274419322172659;
    x[87]  =  0.53499761988709721066;
    x[88]  =  0.55557023301960222474;
    x[89]  =  0.57580819141784530075;
    x[90]  =  0.59569930449243334347;
    x[91]  =  0.61523159058062684548;
    x[92]  =  0.63439328416364549822;
    x[93]  =  0.65317284295377676408;
    x[94]  =  0.67155895484701840063;
    x[95]  =  0.68954054473706692462;
    x[96]  =  0.70710678118654752440;
    x[97]  =  0.72424708295146692094;
    x[98]  =  0.74095112535495909118;
    x[99]  =  0.75720884650648454758;
    x[100] =  0.77301045336273696081;
    x[101] =  0.78834642762660626201;
    x[102] =  0.80320753148064490981;
    x[103] =  0.81758481315158369650;
    x[104] =  0.83146961230254523708;
    x[105] =  0.84485356524970707326;
    x[106] =  0.85772861000027206990;
    x[107] =  0.87008699110871141865;
    x[108] =  0.88192126434835502971;
    x[109] =  0.89322430119551532034;
    x[110] =  0.90398929312344333159;
    x[111] =  0.91420975570353065464;
    x[112] =  0.92387953251128675613;
    x[113] =  0.93299279883473888771;
    x[114] =  0.94154406518302077841;
    x[115] =  0.94952818059303666720;
    x[116] =  0.95694033573220886494;
    x[117] =  0.96377606579543986669;
    x[118] =  0.97003125319454399260;
    x[119] =  0.97570213003852854446;
    x[120] =  0.98078528040323044913;
    x[121] =  0.98527764238894124477;
    x[122] =  0.98917650996478097345;
    x[123] =  0.99247953459870999816;
    x[124] =  0.99518472667219688624;
    x[125] =  0.99729045667869021614;
    x[126] =  0.99879545620517239271;
    x[127] =  0.99969881869620422012;
    x[128] =  1.00000000000000000000;

    w[0]   =  0.00006103888176768602;
    w[1]   =  0.00058807215382869754;
    w[2]   =  0.00120930061875273991;
    w[3]   =  0.00180308126695362360;
    w[4]   =  0.00240715327877140915;
    w[5]   =  0.00300345869904497128;
    w[6]   =  0.00360197835812614147;
    w[7]   =  0.00419553798718534675;
    w[8]   =  0.00478862143341336763;
    w[9]   =  0.00537724746840184621;
    w[10]  =  0.00596388034730799521;
    w[11]  =  0.00654590843862298928;
    w[12]  =  0.00712483332325489785;
    w[13]  =  0.00769875778896082811;
    w[14]  =  0.00826865154203087108;
    w[15]  =  0.00883303867470133581;
    w[16]  =  0.00939256583934814871;
    w[17]  =  0.00994602784923457905;
    w[18]  =  0.01049386202576892125;
    w[19]  =  0.01103504877427254184;
    w[20]  =  0.01156988348290849967;
    w[21]  =  0.01209748052807164113;
    w[22]  =  0.01261803597977743271;
    w[23]  =  0.01313076516693974630;
    w[24]  =  0.01363579321293772047;
    w[25]  =  0.01413241437853094133;
    w[26]  =  0.01462070254634350205;
    w[27]  =  0.01510001572479266783;
    w[28]  =  0.01557039073899425960;
    w[29]  =  0.01603123858745057916;
    w[30]  =  0.01648256956220377909;
    w[31]  =  0.01692383985846499368;
    w[32]  =  0.01735504125411394958;
    w[33]  =  0.01777566938875279997;
    w[34]  =  0.01818570377926339481;
    w[35]  =  0.01858467519566908661;
    w[36]  =  0.01897255587067948426;
    w[37]  =  0.01934890842392451844;
    w[38]  =  0.01971370183700155725;
    w[39]  =  0.02006652805198357604;
    w[40]  =  0.02040735612003867863;
    w[41]  =  0.02073580533490147816;
    w[42]  =  0.02105184759002011131;
    w[43]  =  0.02135512797425970725;
    w[44]  =  0.02164562356712882440;
    w[45]  =  0.02192300400598756892;
    w[46]  =  0.02218725355897195088;
    w[47]  =  0.02243806539722630184;
    w[48]  =  0.02267543270456671718;
    w[49]  =  0.02289907134390605882;
    w[50]  =  0.02310898491627407168;
    w[51]  =  0.02330491126131143273;
    w[52]  =  0.02348686571193163505;
    w[53]  =  0.02365460746057766523;
    w[54]  =  0.02380816473024258975;
    w[55]  =  0.02394731750476901502;
    w[56]  =  0.02407210792327850000;
    w[57]  =  0.02418233623893147567;
    w[58]  =  0.02427805942075745923;
    w[59]  =  0.02435909748927643184;
    w[60]  =  0.02442552306156708690;
    w[61]  =  0.02447717542743444284;
    w[62]  =  0.02451414358881568292;
    w[63]  =  0.02453628559651495473;
    w[64]  =  0.02454370750551418263;
    w[65]  =  0.02453628559651495473;
    w[66]  =  0.02451414358881568292;
    w[67]  =  0.02447717542743444284;
    w[68]  =  0.02442552306156708690;
    w[69]  =  0.02435909748927643184;
    w[70]  =  0.02427805942075745923;
    w[71]  =  0.02418233623893147567;
    w[72]  =  0.02407210792327850000;
    w[73]  =  0.02394731750476901502;
    w[74]  =  0.02380816473024258975;
    w[75]  =  0.02365460746057766523;
    w[76]  =  0.02348686571193163505;
    w[77]  =  0.02330491126131143273;
    w[78]  =  0.02310898491627407168;
    w[79]  =  0.02289907134390605882;
    w[80]  =  0.02267543270456671718;
    w[81]  =  0.02243806539722630184;
    w[82]  =  0.02218725355897195088;
    w[83]  =  0.02192300400598756892;
    w[84]  =  0.02164562356712882440;
    w[85]  =  0.02135512797425970725;
    w[86]  =  0.02105184759002011131;
    w[87]  =  0.02073580533490147816;
    w[88]  =  0.02040735612003867863;
    w[89]  =  0.02006652805198357604;
    w[90]  =  0.01971370183700155725;
    w[91]  =  0.01934890842392451844;
    w[92]  =  0.01897255587067948426;
    w[93]  =  0.01858467519566908661;
    w[94]  =  0.01818570377926339481;
    w[95]  =  0.01777566938875279997;
    w[96]  =  0.01735504125411394958;
    w[97]  =  0.01692383985846499368;
    w[98]  =  0.01648256956220377909;
    w[99]  =  0.01603123858745057916;
    w[100] =  0.01557039073899425960;
    w[101] =  0.01510001572479266783;
    w[102] =  0.01462070254634350205;
    w[103] =  0.01413241437853094133;
    w[104] =  0.01363579321293772047;
    w[105] =  0.01313076516693974630;
    w[106] =  0.01261803597977743271;
    w[107] =  0.01209748052807164113;
    w[108] =  0.01156988348290849967;
    w[109] =  0.01103504877427254184;
    w[110] =  0.01049386202576892125;
    w[111] =  0.00994602784923457905;
    w[112] =  0.00939256583934814871;
    w[113] =  0.00883303867470133581;
    w[114] =  0.00826865154203087108;
    w[115] =  0.00769875778896082811;
    w[116] =  0.00712483332325489785;
    w[117] =  0.00654590843862298928;
    w[118] =  0.00596388034730799521;
    w[119] =  0.00537724746840184621;
    w[120] =  0.00478862143341336763;
    w[121] =  0.00419553798718534675;
    w[122] =  0.00360197835812614147;
    w[123] =  0.00300345869904497128;
    w[124] =  0.00240715327877140915;
    w[125] =  0.00180308126695362360;
    w[126] =  0.00120930061875273991;
    w[127] =  0.00058807215382869754;
    w[128] =  0.00006103888176768602;
  }
  else
  {
    cerr << "\n";
    cerr << "CLENSHAW_CURTIS_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are 1 to 17, 33, 65 or 129.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void fejer1_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER1_COMPUTE computes a Fejer type 1 quadrature rule.
//
//  Discussion:
//
//    This method uses a direct approach.  The paper by Waldvogel
//    exhibits a more efficient approach using Fourier transforms.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  int j;
  static double pi = 3.141592653589793;
  double *theta;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "FEJER1_COMPUTE - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  theta = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    theta[i-1] = ( double ) ( 2 * ( n - i ) + 1 ) * pi
               / ( double ) ( 2 * n     );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( theta[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0;
    for ( j = 1; j <= ( n / 2 ); j++ )
    {
      w[i] = w[i] - 2.0
        * cos ( 2.0 * ( double ) ( j ) * theta[i] ) 
        / ( double ) ( 4 * j * j - 1 );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( n );
  }

  delete [] theta;

  return;
}
//****************************************************************************80

void fejer1_set ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER1_SET sets abscissas and weights for Fejer type 1 quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int ORDER, the order.  ORDER should be
//    between 1 and 10.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0]   =  0.000000000000000;
    weight[0] =  2.000000000000000;
  }
  else if ( order == 2 )
  {
    xtab[0] =   -0.7071067811865475;
    xtab[1] =    0.7071067811865475;

    weight[0] =  1.000000000000000;
    weight[1] =  1.000000000000000;
  }
  else if ( order == 3 )
  {
    xtab[0] =  -0.8660254037844387;
    xtab[1] =   0.0000000000000000;
    xtab[2] =   0.8660254037844387;

    weight[0] =  0.4444444444444444;
    weight[1] =  1.111111111111111;
    weight[2] =  0.4444444444444444;
  }                  
  else if ( order == 4 )
  {
    xtab[0] =  -0.9238795325112867;
    xtab[1] =  -0.3826834323650897;
    xtab[2] =   0.3826834323650898;
    xtab[3] =   0.9238795325112867;

    weight[0] = 0.2642977396044841;
    weight[1] = 0.7357022603955158;
    weight[2] = 0.7357022603955158;
    weight[3] = 0.2642977396044841;
  }
  else if ( order == 5 )
  {
    xtab[0] =  -0.9510565162951535;
    xtab[1] =  -0.5877852522924730;
    xtab[2] =   0.0000000000000000;
    xtab[3] =   0.5877852522924731;
    xtab[4] =   0.9510565162951535;  

    weight[0] = 0.1677812284666835;
    weight[1] = 0.5255521048666498;
    weight[2] = 0.6133333333333333;
    weight[3] = 0.5255521048666498;
    weight[4] = 0.1677812284666835;
  }
  else if ( order == 6 )
  {
    xtab[0] =  -0.9659258262890682;
    xtab[1] =  -0.7071067811865475;
    xtab[2] =  -0.2588190451025206;
    xtab[3] =   0.2588190451025207;
    xtab[4] =   0.7071067811865476;
    xtab[5] =   0.9659258262890683;
      
    weight[0] = 0.1186610213812358;
    weight[1] = 0.3777777777777778;
    weight[2] = 0.5035612008409863;
    weight[3] = 0.5035612008409863;
    weight[4] = 0.3777777777777778;
    weight[5] = 0.1186610213812358;
  }
  else if ( order == 7 )
  {
    xtab[0] =  -0.9749279121818237;
    xtab[1] =  -0.7818314824680295;
    xtab[2] =  -0.4338837391175581;
    xtab[3] =   0.0000000000000000;
    xtab[4] =   0.4338837391175582;
    xtab[5] =   0.7818314824680298;
    xtab[6] =   0.9749279121818236;    
                    
    weight[0] = 0.08671618072672234;
    weight[1] = 0.2878313947886919;
    weight[2] = 0.3982415401308441;
    weight[3] = 0.4544217687074830;
    weight[4] = 0.3982415401308441;
    weight[5] = 0.2878313947886919;
    weight[6] = 0.08671618072672234;
  }
  else if ( order == 8 )
  {
    xtab[0] =  -0.9807852804032304;
    xtab[1] =  -0.8314696123025453;
    xtab[2] =  -0.5555702330196020;
    xtab[3] =  -0.1950903220161282;
    xtab[4] =   0.1950903220161283;
    xtab[5] =   0.5555702330196023;
    xtab[6] =   0.8314696123025452;
    xtab[7] =   0.9807852804032304;
                        
    weight[0] = 0.06698294569858981;
    weight[1] = 0.2229879330145788;
    weight[2] = 0.3241525190645244;
    weight[3] = 0.3858766022223071;
    weight[4] = 0.3858766022223071;
    weight[5] = 0.3241525190645244;
    weight[6] = 0.2229879330145788;
    weight[7] = 0.06698294569858981;
 }
 else if ( order == 9 )
 {
    xtab[0] =  -0.9848077530122080;
    xtab[1] =  -0.8660254037844385;
    xtab[2] =  -0.6427876096865394;
    xtab[3] =  -0.3420201433256685;
    xtab[4] =   0.0000000000000000;
    xtab[5] =   0.3420201433256688;
    xtab[6] =   0.6427876096865394;
    xtab[7] =   0.8660254037844387;
    xtab[8] =   0.9848077530122080;
                       
    weight[0] = 0.05273664990990676;
    weight[1] = 0.1791887125220458;
    weight[2] = 0.2640372225410044;
    weight[3] = 0.3308451751681364;
    weight[4] = 0.3463844797178130;
    weight[5] = 0.3308451751681364;
    weight[6] = 0.2640372225410044;
    weight[7] = 0.1791887125220458;
    weight[8] = 0.05273664990990676;
  }
  else if ( order == 10 )
  {
    xtab[0] =  -0.9876883405951377;
    xtab[1] =  -0.8910065241883678;
    xtab[2] =  -0.7071067811865475;
    xtab[3] =  -0.4539904997395467;
    xtab[4] =  -0.1564344650402306;
    xtab[5] =   0.1564344650402309;
    xtab[6] =   0.4539904997395468;
    xtab[7] =   0.7071067811865476;
    xtab[8] =   0.8910065241883679;
    xtab[9] =   0.9876883405951378;
  
    weight[0] = 0.04293911957413078;
    weight[1] = 0.1458749193773909; 
    weight[2] = 0.2203174603174603;
    weight[3] = 0.2808792186638755;
    weight[4] = 0.3099892820671425;
    weight[5] = 0.3099892820671425;
    weight[6] = 0.2808792186638755;
    weight[7] = 0.2203174603174603;
    weight[8] = 0.1458749193773909;
    weight[9] = 0.04293911957413078;
  }
  else
  {
    cerr << "\n";
    cerr << "FEJER1_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are 1 through 10.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void fejer2_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_COMPUTE computes a Fejer type 2 quadrature rule.
//
//  Discussion:
//
//    This method uses a direct approach.  The paper by Waldvogel
//    exhibits a more efficient approach using Fourier transforms.
//
//    Thanks to Mike Neish, University of Toronto, for pointing out a
//    memory leak in a previous version of this code, 13 October 2009.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  int j;
  double p;
  static double pi = 3.141592653589793;
  double *theta;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "FEJER2_COMPUTE - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  theta = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    theta[i-1] = ( double ) ( n + 1 - i ) * pi
               / ( double ) ( n + 1     );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( theta[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0;

    for ( j = 1; j <= ( ( n - 1 ) / 2 ); j++ )
    {
      w[i] = w[i] - 2.0 * cos ( 2.0 * ( double ) ( j ) * theta[i] ) 
        / ( double ) ( 4 * j * j - 1 );
    }

    if ( 2 < n )
    {
      p = 2.0 * ( double ) ( ( n + 1 ) / 2 ) - 1.0;
      w[i] = w[i] - cos ( ( p + 1.0 ) * theta[i] ) / p;
    }

  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( n + 1 );
  }

  delete [] theta;

  return;
}
//****************************************************************************80

void fejer2_set ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_SET sets abscissas and weights for Fejer type 2 quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int ORDER, the order.  
//    ORDER should be between 1 and 10.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0]   =  0.000000000000000;
    weight[0] =  2.000000000000000;
  }
  else if ( order == 2 )
  {
    xtab[0] =   -0.5000000000000000;
    xtab[1] =    0.5000000000000000;

    weight[0] =  1.0000000000000000;
    weight[1] =  1.0000000000000000;
  }
  else if ( order == 3 )
  {
    xtab[0] =  -0.7071067811865476;
    xtab[1] =   0.0000000000000000;
    xtab[2] =   0.7071067811865476;

    weight[0] =  0.6666666666666666;
    weight[1] =  0.6666666666666666;
    weight[2] =  0.6666666666666666;
  }                  
  else if ( order == 4 )
  {
    xtab[0] =  -0.8090169943749475;
    xtab[1] =  -0.3090169943749475;
    xtab[2] =   0.3090169943749475;
    xtab[3] =   0.8090169943749475;

    weight[0] = 0.4254644007500070;
    weight[1] = 0.5745355992499930;
    weight[2] = 0.5745355992499930;
    weight[3] = 0.4254644007500070;
  }
  else if ( order == 5 )
  {
    xtab[0] =  -0.8660254037844387;
    xtab[1] =  -0.5000000000000000;
    xtab[2] =   0.0000000000000000;
    xtab[3] =   0.5000000000000000;
    xtab[4] =   0.8660254037844387;  

    weight[0] = 0.3111111111111111;
    weight[1] = 0.4000000000000000;
    weight[2] = 0.5777777777777777;
    weight[3] = 0.4000000000000000;
    weight[4] = 0.3111111111111111;
  }
  else if ( order == 6 )
  {
    xtab[0] =  -0.9009688679024191;
    xtab[1] =  -0.6234898018587336;
    xtab[2] =  -0.2225209339563144;
    xtab[3] =   0.2225209339563144;
    xtab[4] =   0.6234898018587336;
    xtab[5] =   0.9009688679024191;
      
    weight[0] = 0.2269152467244296;
    weight[1] = 0.3267938603769863;
    weight[2] = 0.4462908928985841;
    weight[3] = 0.4462908928985841;
    weight[4] = 0.3267938603769863;
    weight[5] = 0.2269152467244296;
  }
  else if ( order == 7 )
  {
    xtab[0] =  -0.9238795325112867;
    xtab[1] =  -0.7071067811865476;
    xtab[2] =  -0.3826834323650898;
    xtab[3] =   0.0000000000000000;
    xtab[4] =   0.3826834323650898;
    xtab[5] =   0.7071067811865476;
    xtab[6] =   0.9238795325112867;    
                    
    weight[0] = 0.1779646809620499;
    weight[1] = 0.2476190476190476;
    weight[2] = 0.3934638904665215;
    weight[3] = 0.3619047619047619;
    weight[4] = 0.3934638904665215;
    weight[5] = 0.2476190476190476;
    weight[6] = 0.1779646809620499;
  }
  else if ( order == 8 )
  {
    xtab[0] =  -0.9396926207859084;
    xtab[1] =  -0.7660444431189780;
    xtab[2] =  -0.5000000000000000;
    xtab[3] =  -0.1736481776669304;
    xtab[4] =   0.1736481776669304;
    xtab[5] =   0.5000000000000000;
    xtab[6] =   0.7660444431189780;
    xtab[7] =   0.9396926207859084;
                        
    weight[0] = 0.1397697435050225;
    weight[1] = 0.2063696457302284;
    weight[2] = 0.3142857142857143;
    weight[3] = 0.3395748964790348;
    weight[4] = 0.3395748964790348;
    weight[5] = 0.3142857142857143;
    weight[6] = 0.2063696457302284;
    weight[7] = 0.1397697435050225;
  }
  else if ( order == 9 )
  {
    xtab[0] =  -0.9510565162951535;
    xtab[1] =  -0.8090169943749475;
    xtab[2] =  -0.5877852522924731;
    xtab[3] =  -0.3090169943749475;
    xtab[4] =   0.0000000000000000;
    xtab[5] =   0.3090169943749475;
    xtab[6] =   0.5877852522924731;
    xtab[7] =   0.8090169943749475;
    xtab[8] =   0.9510565162951535;
                       
    weight[0] = 0.1147810750857217;
    weight[1] = 0.1654331942222276;
    weight[2] = 0.2737903534857068;
    weight[3] = 0.2790112502222170;
    weight[4] = 0.3339682539682539;
    weight[5] = 0.2790112502222170;
    weight[6] = 0.2737903534857068;
    weight[7] = 0.1654331942222276;
    weight[8] = 0.1147810750857217;
  }
  else if ( order == 10 )
  {
    xtab[0] =  -0.9594929736144974;
    xtab[1] =  -0.8412535328311812;
    xtab[2] =  -0.6548607339452851;
    xtab[3] =  -0.4154150130018864;
    xtab[4] =  -0.1423148382732851;
    xtab[5] =   0.1423148382732851;
    xtab[6] =   0.4154150130018864;
    xtab[7] =   0.6548607339452851;
    xtab[8] =   0.8412535328311812;
    xtab[9] =   0.9594929736144974;
  
    weight[0] = 0.09441954173982806;
    weight[1] = 0.1411354380109716; 
    weight[2] = 0.2263866903636005;
    weight[3] = 0.2530509772156453;
    weight[4] = 0.2850073526699544;
    weight[5] = 0.2850073526699544;
    weight[6] = 0.2530509772156453;
    weight[7] = 0.2263866903636005;
    weight[8] = 0.1411354380109716;
    weight[9] = 0.09441954173982806;
  }
  else
  {
    cerr << "\n";
    cerr << "FEJER2_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are 1 through 10.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void gegenbauer_compute ( int order, double alpha, double xtab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_COMPUTE computes a Gauss-Gegenbauer quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//    Thanks to Janiki Raman for pointing out a problem in an earlier
//    version of the code that occurred when ALPHA was -0.5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the quadrature rule.
//
//    Input, double ALPHA, the exponent of (1-X^2) in the weight.  
//    -1.0 < ALPHA is required.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double an;
  double *c;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double temp;
  double x;
//
//  Check ORDER.
//
  if ( order < 1 )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_COMPUTE - Fatal error!\n";
    cerr << "  1 <= ORDER is required.\n";
    exit ( 1 );
  }
  
  c = new double[order];
//
//  Check ALPHA.
//
  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_COMPUTE - Fatal error!\n";
    cerr << "  -1.0 < ALPHA is required.\n";
    exit ( 1 );
  }
//
//  Set the recursion coefficients.
//
  c[0] = 0.0;
  if ( 2 <= order )
  {
    c[1] = 1.0 / ( 2.0 * alpha + 3.0 );
  }

  for ( i = 3; i <= order; i++ )
  {
    c[i-1] = ( double ) ( i - 1 ) 
          * ( alpha + alpha + ( double ) ( i - 1 ) ) / 
          ( ( alpha + alpha + ( double ) ( 2 * i - 1 ) ) 
          * ( alpha + alpha + ( double ) ( 2 * i - 3 ) ) );
  }

  delta = gamma ( alpha         + 1.0 ) 
        * gamma (         alpha + 1.0 ) 
        / gamma ( alpha + alpha + 2.0 );

  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + alpha + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );

      r1 = ( 1.0 + alpha ) 
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) ) 
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 2.44 * an + 1.282 * an * an;

      x = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) / 
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) * 
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * alpha * 
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x = x - r1 * r2 * r3 * ( 1.0 - x );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) );

      x = x - r1 * r2 * r3 * ( xtab[0] - x );
    }
    else if ( i < order - 1 )
    {
      x = 3.0 * xtab[i-2] - 3.0 * xtab[i-3] + xtab[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * alpha ) / ( 0.766 + 0.119 * alpha );

      r2 = 1.0 / ( 1.0 + 0.639 
        * ( ( double ) ( order ) - 4.0 ) 
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * 
        ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * alpha ) / ( 1.67 + 0.28 * alpha );

      r2 = 1.0 / 
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }

    gegenbauer_root ( &x, order, alpha, &dp2, &p1, c );

    xtab[i-1] = x;
    weight[i-1] = cc / ( dp2 * p1 );
  }
//
//  Reverse the order of the values.
//
  for ( i = 1; i <= order/2; i++ )
  {
    temp          = xtab[i-1];
    xtab[i-1]     = xtab[order-i];
    xtab[order-i] = temp;
  }

  for ( i = 1; i <=order/2; i++ )
  {
    temp            = weight[i-1];
    weight[i-1]     = weight[order-i];
    weight[order-i] = temp;
  }

  delete [] c;

  return;
}
//****************************************************************************80

double gegenbauer_integral ( int expon, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_INTEGRAL evaluates the integral of a monomial with Gegenbauer weight.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Input, double ALPHA, the exponent of (1-X^2) in the weight factor.
//
//    Output, double GEGENBAUER_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  c = ( double ) ( expon );

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = gamma ( 1.0 + c ) * 2.0 
    * gamma ( 1.0 + alpha  ) * value1 
    / gamma ( 2.0 + alpha  + c );

  return value;
}
//****************************************************************************80

void gegenbauer_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_RECUR finds the value and derivative of a Gegenbauer polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of J(ORDER)(X).
//
//    Output, double *DP2, the value of J'(ORDER)(X).
//
//    Output, double *P1, the value of J(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponents of (1-X^2).
//
//    Input, double C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x;
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = x *  ( *p1 ) - c[i-1] * p0;
    *dp2 = x * dp1 + ( *p1 ) - c[i-1] * dp0;
  }
  return;
}
//****************************************************************************80

void gegenbauer_root ( double *x, int order, double alpha,  double *dp2, 
  double *p1, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_ROOT improves an approximate root of a Gegenbauer polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponents of (1-X^2).
//
//    Output, double *DP2, the value of J'(ORDER)(X).
//
//    Output, double *P1, the value of J(ORDER-1)(X).
//
//    Input, double C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gegenbauer_recur ( &p2, dp2, p1, *x, order, alpha, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
//****************************************************************************80

void gen_hermite_dr_compute ( int order, double alpha, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_DR_COMPUTE computes a generalized Gauss-Hermite rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double W[ORDER], the weights.
//
{
  double alpha_laguerre;
  double arg;
  int i;
  int order_laguerre;
  double *w_laguerre;
  double *x_laguerre;

  if ( order == 1 )
  {
    arg = ( alpha + 1.0 ) / 2.0;
    x[0] = 0.0;
    w[0] = gamma ( arg );
    return;
  }

  if ( ( order % 2 ) == 0 ) 
  {
    order_laguerre = order / 2;
    alpha_laguerre = ( alpha - 1.0 ) / 2.0;
  }
  else
  {
    order_laguerre = ( order - 1 ) / 2;
    alpha_laguerre = ( alpha + 1.0 ) / 2.0;
  }
  
  w_laguerre = new double[order_laguerre];
  x_laguerre = new double[order_laguerre];

  gen_laguerre_ss_compute ( order_laguerre, alpha_laguerre, 
    x_laguerre, w_laguerre );

  if ( ( order % 2 ) == 0 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i];
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+i] = 0.5 * w_laguerre[i];
    }
  }
  else if ( ( order % 2 ) == 1 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    x[order_laguerre] = 0.0;
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+1+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i] / x_laguerre[order_laguerre-1-i];
    }

    arg = ( alpha + 1.0 ) / 2.0;
    w[order_laguerre] = gamma ( arg );
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre] = w[order_laguerre] - w_laguerre[i] / x_laguerre[i];
    }

    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+1+i] = 0.5 * w_laguerre[i] / x_laguerre[i];
    }
  }
  delete [] w_laguerre;
  delete [] x_laguerre;

  return;
}
//****************************************************************************80

void gen_hermite_ek_compute ( int n, double alpha, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_EK_COMPUTE computes a generalized Gauss-Hermite quadrature rule.
//
//  Discussion:
//
//    The code uses an algorithm by Elhay and Kautsky.
//
//    The integral:
//
//      integral ( -oo < x < +oo ) |x|^alpha exp(-x^2) f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the number of abscissas.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *bj;
  int i;
  double i_r8;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = gamma ( ( alpha + 1.0 ) / 2.0 );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];

  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    if ( ( i % 2 ) == 0 )
    {
      bj[i] = ( i_r8 + alpha ) / 2.0;
    }
    else
    {
      bj[i] = i_r8 / 2.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    bj[i] = sqrt ( bj[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double gen_hermite_integral ( int expon, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_INTEGRAL evaluates a monomial generalized Hermite integral.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x^2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent of the monomial.
//    0 <= EXPON.
//
//    Input, double ALPHA, the exponent of |X| in the weight function.
//    -1.0 < ALPHA.
//
//    Output, double GEN_HERMITE_INTEGRAL, the value of the integral.
//
{
  double a;
  double arg;
  double value;

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    a = alpha + ( double ) ( expon );
    if ( a <= - 1.0 )
    {
      value = - r8_huge ( );
    }
    else
    {
      arg = ( a + 1.0 ) / 2.0;
      value = gamma ( arg );
    }
  }
  return value;
}
//****************************************************************************80

void gen_laguerre_ek_compute ( int n, double alpha, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_EK_COMPUTE: generalized Gauss-Laguerre quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( 0 <= x < +oo ) exp ( - x ) * x^alpha * f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    The integral:
//
//      integral ( 0 <= x < +oo ) x^alpha * f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * exp ( x(i) ) * f ( x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, double ALPHA, the exponent of the X factor.
//    ALPHA must be nonnegative.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *bj;
  int i;
  double i_r8;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = gamma ( alpha + 1.0 );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];

  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    bj[i] = sqrt ( i_r8 * ( i_r8 + alpha ) );
  }

  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    x[i] = 2.0 * i_r8 - 1.0 + alpha;
  }

  w[0] = sqrt ( zemu );

  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double gen_laguerre_integral ( int expon, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_INTEGRAL evaluates a monomial generalized Laguerre integral.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= x < +oo ) x^n * x^alpha exp(-x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent of the monomial.
//    0 <= EXPON.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, double GEN_LAGUERRE_INTEGRAL, the value of the integral.
//
{
  double arg;
  double value;

  arg = alpha + ( double ) ( expon + 1.0 );
  value = gamma ( arg );

  return value;
}
//****************************************************************************80

void gen_laguerre_ss_compute ( int order, double alpha, double xtab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_SS_COMPUTE computes a generalized Gauss-Laguerre quadrature rule.
//
//  Discussion:
//
//    In the simplest case, ALPHA is 0, and we are approximating the
//    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
//    it is easy to modify the rule to approximate the integral from
//    A to +oo as well.
//
//    If ALPHA is nonzero, then there is no simple way to extend the
//    rule to approximate the integral from A to +oo.  The simplest
//    procedures would be to approximate the integral from 0 to A.
//
//    If the integral is:
//
//        Integral ( A <= X < +oo ) EXP ( - X ) * F(X) dX
//      or
//        Integral ( 0 <= X < +oo ) EXP ( - X ) * X**ALPHA * F(X) dX
//
//    then the quadrature rule is:
//
//      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
//    or
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//
//    If the integral is:
//
//        Integral ( A <= X < +oo ) F(X) dX
//      or
//        Integral ( 0 <= X < +oo ) X**ALPHA * F(X) dX
//
//    then the quadrature rule is:
//
//      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) 
//        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
//    or
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the quadrature rule to be computed.
//    ORDER must be at least 1.
//
//    Input, double ALPHA, the exponent of the X factor.
//    Set ALPHA = 0.0 for the simplest rule.
//    ALPHA must be nonnegative.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double *b;
  double *c;
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double ratio;
  double xval;

  b = new double[order];
  c = new double[order];
//
//  Set the recursion coefficients.
//
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( alpha + ( double ) ( 2 * i + 1 ) );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i ) * ( alpha + ( double ) ( i ) );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = gamma ( alpha + 1.0 ) * prod;

  for ( i = 0; i < order; i++ )
  {
//
//  Compute an estimate for the root.
//
    if ( i == 0 )
    {
      xval = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) / 
        ( 1.0 + 2.4 * ( double ) ( order ) + 1.8 * alpha );
    }
    else if ( i == 1 )
    {
      xval = xval + ( 15.0 + 6.25 * alpha ) / 
        ( 1.0 + 0.9 * alpha + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) ) 
        / ( 1.9 * ( double ) ( i - 1 ) );

      r2 = 1.26 * ( double ) ( i - 1 ) * alpha / 
        ( 1.0 + 3.5 * ( double ) ( i - 1 ) );

      ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );

      xval = xval + ratio * ( xval - xtab[i-2] );
    }
//
//  Use iteration to find the root.
//
    gen_laguerre_ss_root ( &xval, order, alpha, &dp2, &p1, b, c );
//
//  Set the abscissa and weight.
//
    xtab[i] = xval;
    weight[i] = ( cc / dp2 ) / p1;
  }

  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void gen_laguerre_ss_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_SS_RECUR evaluates a generalized Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of L(ORDER)(X).
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponent of the X factor in the
//    integrand.
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - alpha - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
//****************************************************************************80

void gen_laguerre_ss_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_SS_ROOT improves a root of a generalized Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponent of the X factor.
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gen_laguerre_ss_recur ( &p2, dp2, p1, *x, order, alpha, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void hermite_ek_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_EK_COMPUTE computes a Gauss-Hermite quadrature rule.
//
//  Discussion:
//
//    The code uses an algorithm by Elhay and Kautsky.
//
//    The abscissas are the zeros of the N-th order Hermite polynomial.
//
//    The integral:
//
//      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2012
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the number of abscissas.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double arg;
  double *bj;
  int i;
  double zemu;
//
//  Define the zero-th moment.
//
  arg = 0.5;
  zemu = gamma ( arg );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];

  for ( i = 0; i < n; i++ )
  {
    bj[i] = sqrt ( ( double ) ( i + 1 ) / 2.0 );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );
//
//  If N is odd, force the center X to be 0.
//
  if ( ( n % 2 ) == 1 )
  {
    x[(n-1)/2] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

void hermite_gk16_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK16_SET sets a Hermite Genz-Keister 16 rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    A nested family of rules for the Hermite integration problem
//    was produced by Genz and Keister.  The structure of the nested
//    family was denoted by 1+2+6+10+16, that is, it comprised rules 
//    of successive orders O = 1, 3, 9, 19, and 35.
//
//    The precisions of these rules are P = 1, 5, 15, 29, and 51.
//
//    Consider, however, the special cases where a rule of precision 
//    at least 7, 17, 31 or 33 is desired.  Ordinarily, this would
//    suggest using the nested rule of order 9, 19, 51 or 51 respectively.
//    In these cases, however, the order of the rule that is used exceeds
//    the precision requested.  Hence, it is possible simply to select
//    a subset of the points in the higher precision rule and get a 
//    rule of lower order and the desired precision.  This accounts for
//    the four extra rules in this family.
//
//    The entire list of rules is therefore:
//
//    L   P   N
//
//    0   1   1  <-- Full rule
//    1   5   3  <-- Full rule
//    2   7   7  <-- Partial rule
//    3  15   9  <-- Full rule
//    4  17  17  <-- Partial rule
//    5  29  19  <-- Full rule
//    6  31  31  <-- Partial rule
//    7  33  33  <-- Partial rule
//    8  51  35  <-- Full rule 4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz, Bradley Keister,
//    Fully symmetric interpolatory rules for multiple integrals
//    over infinite regions with Gaussian weight,
//    Journal of Computational and Applied Mathematics,
//    Volume 71, 1996, pages 299-309
//
//    Florian Heiss, Viktor Winschel,
//    Likelihood approximation by numerical integration on sparse grids,
//    Journal of Econometrics,
//    Volume 144, 2008, pages 62-80.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 7, 9, 17, 19, 31, 33 or 35.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[ 0] =   0.0000000000000000E+00;

    w[ 0] =   1.7724538509055159E+00;
  }
  else if ( n == 3 )
  {
    x[ 0] =  -1.2247448713915889E+00;
    x[ 1] =   0.0000000000000000E+00;
    x[ 2] =   1.2247448713915889E+00;

    w[ 0] =   2.9540897515091930E-01;
    w[ 1] =   1.1816359006036772E+00;
    w[ 2] =   2.9540897515091930E-01;
  }
  else if ( n == 7 )
  {
    x[ 0] =  -2.9592107790638380E+00;
    x[ 1] =  -1.2247448713915889E+00;
    x[ 2] =  -5.2403354748695763E-01;
    x[ 3] =   0.0000000000000000E+00;
    x[ 4] =   5.2403354748695763E-01;
    x[ 5] =   1.2247448713915889E+00;
    x[ 6] =   2.9592107790638380E+00;

    w[ 0] =   1.2330680655153448E-03;
    w[ 1] =   2.4557928535031393E-01;
    w[ 2] =   2.3286251787386100E-01;
    w[ 3] =   8.1310410832613500E-01;
    w[ 4] =   2.3286251787386100E-01;
    w[ 5] =   2.4557928535031393E-01;
    w[ 6] =   1.2330680655153448E-03;
  }
  else if ( n == 9 )
  {
    x[ 0] =  -2.9592107790638380E+00;
    x[ 1] =  -2.0232301911005157E+00;
    x[ 2] =  -1.2247448713915889E+00;
    x[ 3] =  -5.2403354748695763E-01;
    x[ 4] =   0.0000000000000000E+00;
    x[ 5] =   5.2403354748695763E-01;
    x[ 6] =   1.2247448713915889E+00;
    x[ 7] =   2.0232301911005157E+00;
    x[ 8] =   2.9592107790638380E+00;

    w[ 0] =   1.6708826306882348E-04;
    w[ 1] =   1.4173117873979098E-02;
    w[ 2] =   1.6811892894767771E-01;
    w[ 3] =   4.7869428549114124E-01;
    w[ 4] =   4.5014700975378197E-01;
    w[ 5] =   4.7869428549114124E-01;
    w[ 6] =   1.6811892894767771E-01;
    w[ 7] =   1.4173117873979098E-02;
    w[ 8] =   1.6708826306882348E-04;
  }
  else if ( n == 17 )
  {
    x[ 0] =  -4.4995993983103881E+00;
    x[ 1] =  -3.6677742159463378E+00;
    x[ 2] =  -2.9592107790638380E+00;
    x[ 3] =  -2.0232301911005157E+00;
    x[ 4] =  -1.8357079751751868E+00;
    x[ 5] =  -1.2247448713915889E+00;
    x[ 6] =  -8.7004089535290285E-01;
    x[ 7] =  -5.2403354748695763E-01;
    x[ 8] =   0.0000000000000000E+00;
    x[ 9] =   5.2403354748695763E-01;
    x[10] =   8.7004089535290285E-01;
    x[11] =   1.2247448713915889E+00;
    x[12] =   1.8357079751751868E+00;
    x[13] =   2.0232301911005157E+00;
    x[14] =   2.9592107790638380E+00;
    x[15] =   3.6677742159463378E+00;
    x[16] =   4.4995993983103881E+00;

    w[ 0] =   3.7463469943051758E-08;
    w[ 1] =  -1.4542843387069391E-06;
    w[ 2] =   1.8723818949278350E-04;
    w[ 3] =   1.2466519132805918E-02;
    w[ 4] =   3.4840719346803800E-03;
    w[ 5] =   1.5718298376652240E-01;
    w[ 6] =   2.5155825701712934E-02;
    w[ 7] =   4.5119803602358544E-01;
    w[ 8] =   4.7310733504965385E-01;
    w[ 9] =   4.5119803602358544E-01;
    w[10] =   2.5155825701712934E-02;
    w[11] =   1.5718298376652240E-01;
    w[12] =   3.4840719346803800E-03;
    w[13] =   1.2466519132805918E-02;
    w[14] =   1.8723818949278350E-04;
    w[15] =  -1.4542843387069391E-06;
    w[16] =   3.7463469943051758E-08;
  }
  else if ( n == 19 )
  {
    x[ 0] =  -4.4995993983103881E+00;
    x[ 1] =  -3.6677742159463378E+00;
    x[ 2] =  -2.9592107790638380E+00;
    x[ 3] =  -2.2665132620567876E+00;
    x[ 4] =  -2.0232301911005157E+00;
    x[ 5] =  -1.8357079751751868E+00;
    x[ 6] =  -1.2247448713915889E+00;
    x[ 7] =  -8.7004089535290285E-01;
    x[ 8] =  -5.2403354748695763E-01;
    x[ 9] =   0.0000000000000000E+00;
    x[10] =   5.2403354748695763E-01;
    x[11] =   8.7004089535290285E-01;
    x[12] =   1.2247448713915889E+00;
    x[13] =   1.8357079751751868E+00;
    x[14] =   2.0232301911005157E+00;
    x[15] =   2.2665132620567876E+00;
    x[16] =   2.9592107790638380E+00;
    x[17] =   3.6677742159463378E+00;
    x[18] =   4.4995993983103881E+00;

    w[ 0] =   1.5295717705322357E-09;
    w[ 1] =   1.0802767206624762E-06;
    w[ 2] =   1.0656589772852267E-04;
    w[ 3] =   5.1133174390883855E-03;
    w[ 4] =  -1.1232438489069229E-02;
    w[ 5] =   3.2055243099445879E-02;
    w[ 6] =   1.1360729895748269E-01;
    w[ 7] =   1.0838861955003017E-01;
    w[ 8] =   3.6924643368920851E-01;
    w[ 9] =   5.3788160700510168E-01;
    w[10] =   3.6924643368920851E-01;
    w[11] =   1.0838861955003017E-01;
    w[12] =   1.1360729895748269E-01;
    w[13] =   3.2055243099445879E-02;
    w[14] =  -1.1232438489069229E-02;
    w[15] =   5.1133174390883855E-03;
    w[16] =   1.0656589772852267E-04;
    w[17] =   1.0802767206624762E-06;
    w[18] =   1.5295717705322357E-09;
  }
  else if ( n == 31 )
  {
    x[ 0] =  -6.3759392709822356E+00;
    x[ 1] =  -5.6432578578857449E+00;
    x[ 2] =  -5.0360899444730940E+00;
    x[ 3] =  -4.4995993983103881E+00;
    x[ 4] =  -3.6677742159463378E+00;
    x[ 5] =  -2.9592107790638380E+00;
    x[ 6] =  -2.5705583765842968E+00;
    x[ 7] =  -2.2665132620567876E+00;
    x[ 8] =  -2.0232301911005157E+00;
    x[ 9] =  -1.8357079751751868E+00;
    x[10] =  -1.5794121348467671E+00;
    x[11] =  -1.2247448713915889E+00;
    x[12] =  -8.7004089535290285E-01;
    x[13] =  -5.2403354748695763E-01;
    x[14] =  -1.7606414208200893E-01;
    x[15] =   0.0000000000000000E+00;
    x[16] =   1.7606414208200893E-01;
    x[17] =   5.2403354748695763E-01;
    x[18] =   8.7004089535290285E-01;
    x[19] =   1.2247448713915889E+00;
    x[20] =   1.5794121348467671E+00;
    x[21] =   1.8357079751751868E+00;
    x[22] =   2.0232301911005157E+00;
    x[23] =   2.2665132620567876E+00;
    x[24] =   2.5705583765842968E+00;
    x[25] =   2.9592107790638380E+00;
    x[26] =   3.6677742159463378E+00;
    x[27] =   4.4995993983103881E+00;
    x[28] =   5.0360899444730940E+00;
    x[29] =   5.6432578578857449E+00;
    x[30] =   6.3759392709822356E+00;

    w[ 0] =   2.2365645607044459E-15;
    w[ 1] =  -2.6304696458548942E-13;
    w[ 2] =   9.0675288231679823E-12;
    w[ 3] =   1.4055252024722478E-09;
    w[ 4] =   1.0889219692128120E-06;
    w[ 5] =   1.0541662394746661E-04;
    w[ 6] =   2.6665159778939428E-05;
    w[ 7] =   4.8385208205502612E-03;
    w[ 8] =  -9.8566270434610019E-03;
    w[ 9] =   2.9409427580350787E-02;
    w[10] =   3.1210210352682834E-03;
    w[11] =   1.0939325071860877E-01;
    w[12] =   1.1594930984853116E-01;
    w[13] =   3.5393889029580544E-01;
    w[14] =   4.9855761893293160E-02;
    w[15] =   4.5888839636756751E-01;
    w[16] =   4.9855761893293160E-02;
    w[17] =   3.5393889029580544E-01;
    w[18] =   1.1594930984853116E-01;
    w[19] =   1.0939325071860877E-01;
    w[20] =   3.1210210352682834E-03;
    w[21] =   2.9409427580350787E-02;
    w[22] =  -9.8566270434610019E-03;
    w[23] =   4.8385208205502612E-03;
    w[24] =   2.6665159778939428E-05;
    w[25] =   1.0541662394746661E-04;
    w[26] =   1.0889219692128120E-06;
    w[27] =   1.4055252024722478E-09;
    w[28] =   9.0675288231679823E-12;
    w[29] =  -2.6304696458548942E-13;
    w[30] =   2.2365645607044459E-15;
  }
  else if ( n == 33 )
  {
    x[ 0] =  -6.3759392709822356E+00;
    x[ 1] =  -5.6432578578857449E+00;
    x[ 2] =  -5.0360899444730940E+00;
    x[ 3] =  -4.4995993983103881E+00;
    x[ 4] =  -4.0292201405043713E+00;
    x[ 5] =  -3.6677742159463378E+00;
    x[ 6] =  -2.9592107790638380E+00;
    x[ 7] =  -2.5705583765842968E+00;
    x[ 8] =  -2.2665132620567876E+00;
    x[ 9] =  -2.0232301911005157E+00;
    x[10] =  -1.8357079751751868E+00;
    x[11] =  -1.5794121348467671E+00;
    x[12] =  -1.2247448713915889E+00;
    x[13] =  -8.7004089535290285E-01;
    x[14] =  -5.2403354748695763E-01;
    x[15] =  -1.7606414208200893E-01;
    x[16] =   0.0000000000000000E+00;
    x[17] =   1.7606414208200893E-01;
    x[18] =   5.2403354748695763E-01;
    x[19] =   8.7004089535290285E-01;
    x[20] =   1.2247448713915889E+00;
    x[21] =   1.5794121348467671E+00;
    x[22] =   1.8357079751751868E+00;
    x[23] =   2.0232301911005157E+00;
    x[24] =   2.2665132620567876E+00;
    x[25] =   2.5705583765842968E+00;
    x[26] =   2.9592107790638380E+00;
    x[27] =   3.6677742159463378E+00;
    x[28] =   4.0292201405043713E+00;
    x[29] =   4.4995993983103881E+00;
    x[30] =   5.0360899444730940E+00;
    x[31] =   5.6432578578857449E+00;
    x[32] =   6.3759392709822356E+00;

    w[ 0] =  -1.7602932805372496E-15;
    w[ 1] =   4.7219278666417693E-13;
    w[ 2] =  -3.4281570530349562E-11;
    w[ 3] =   2.7547825138935901E-09;
    w[ 4] =  -2.3903343382803510E-08;
    w[ 5] =   1.2245220967158438E-06;
    w[ 6] =   9.8710009197409173E-05;
    w[ 7] =   1.4753204901862772E-04;
    w[ 8] =   3.7580026604304793E-03;
    w[ 9] =  -4.9118576123877555E-03;
    w[10] =   2.0435058359107205E-02;
    w[11] =   1.3032872699027960E-02;
    w[12] =   9.6913444944583621E-02;
    w[13] =   1.3726521191567551E-01;
    w[14] =   3.1208656194697448E-01;
    w[15] =   1.8411696047725790E-01;
    w[16] =   2.4656644932829619E-01;
    w[17] =   1.8411696047725790E-01;
    w[18] =   3.1208656194697448E-01;
    w[19] =   1.3726521191567551E-01;
    w[20] =   9.6913444944583621E-02;
    w[21] =   1.3032872699027960E-02;
    w[22] =   2.0435058359107205E-02;
    w[23] =  -4.9118576123877555E-03;
    w[24] =   3.7580026604304793E-03;
    w[25] =   1.4753204901862772E-04;
    w[26] =   9.8710009197409173E-05;
    w[27] =   1.2245220967158438E-06;
    w[28] =  -2.3903343382803510E-08;
    w[29] =   2.7547825138935901E-09;
    w[30] =  -3.4281570530349562E-11;
    w[31] =   4.7219278666417693E-13;
    w[32] =  -1.7602932805372496E-15;
  }
  else if ( n == 35 )
  {
    x[ 0] =  -6.3759392709822356E+00;
    x[ 1] =  -5.6432578578857449E+00;
    x[ 2] =  -5.0360899444730940E+00;
    x[ 3] =  -4.4995993983103881E+00;
    x[ 4] =  -4.0292201405043713E+00;
    x[ 5] =  -3.6677742159463378E+00;
    x[ 6] =  -3.3491639537131945E+00;
    x[ 7] =  -2.9592107790638380E+00;
    x[ 8] =  -2.5705583765842968E+00;
    x[ 9] =  -2.2665132620567876E+00;
    x[10] =  -2.0232301911005157E+00;
    x[11] =  -1.8357079751751868E+00;
    x[12] =  -1.5794121348467671E+00;
    x[13] =  -1.2247448713915889E+00;
    x[14] =  -8.7004089535290285E-01;
    x[15] =  -5.2403354748695763E-01;
    x[16] =  -1.7606414208200893E-01;
    x[17] =   0.0000000000000000E+00;
    x[18] =   1.7606414208200893E-01;
    x[19] =   5.2403354748695763E-01;
    x[20] =   8.7004089535290285E-01;
    x[21] =   1.2247448713915889E+00;
    x[22] =   1.5794121348467671E+00;
    x[23] =   1.8357079751751868E+00;
    x[24] =   2.0232301911005157E+00;
    x[25] =   2.2665132620567876E+00;
    x[26] =   2.5705583765842968E+00;
    x[27] =   2.9592107790638380E+00;
    x[28] =   3.3491639537131945E+00;
    x[29] =   3.6677742159463378E+00;
    x[30] =   4.0292201405043713E+00;
    x[31] =   4.4995993983103881E+00;
    x[32] =   5.0360899444730940E+00;
    x[33] =   5.6432578578857449E+00;
    x[34] =   6.3759392709822356E+00;

    w[ 0] =   1.8684014894510604E-18;
    w[ 1] =   9.6599466278563243E-15;
    w[ 2] =   5.4896836948499462E-12;
    w[ 3] =   8.1553721816916897E-10;
    w[ 4] =   3.7920222392319532E-08;
    w[ 5] =   4.3737818040926989E-07;
    w[ 6] =   4.8462799737020461E-06;
    w[ 7] =   6.3328620805617891E-05;
    w[ 8] =   4.8785399304443770E-04;
    w[ 9] =   1.4515580425155904E-03;
    w[10] =   4.0967527720344047E-03;
    w[11] =   5.5928828911469180E-03;
    w[12] =   2.7780508908535097E-02;
    w[13] =   8.0245518147390893E-02;
    w[14] =   1.6371221555735804E-01;
    w[15] =   2.6244871488784277E-01;
    w[16] =   3.3988595585585218E-01;
    w[17] =   9.1262675363737921E-04;
    w[18] =   3.3988595585585218E-01;
    w[19] =   2.6244871488784277E-01;
    w[20] =   1.6371221555735804E-01;
    w[21] =   8.0245518147390893E-02;
    w[22] =   2.7780508908535097E-02;
    w[23] =   5.5928828911469180E-03;
    w[24] =   4.0967527720344047E-03;
    w[25] =   1.4515580425155904E-03;
    w[26] =   4.8785399304443770E-04;
    w[27] =   6.3328620805617891E-05;
    w[28] =   4.8462799737020461E-06;
    w[29] =   4.3737818040926989E-07;
    w[30] =   3.7920222392319532E-08;
    w[31] =   8.1553721816916897E-10;
    w[32] =   5.4896836948499462E-12;
    w[33] =   9.6599466278563243E-15;
    w[34] =   1.8684014894510604E-18;
  }
  else
  {
    cerr << "\n";
    cerr << "HERMITE_GK16_SET - Fatal error!\n";
    cerr << "  Illegal input value of N.\n";
    cerr << "  N must be 1, 3, 7, 9, 17, 19, 31, 33 or 35.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void hermite_gk18_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK18_SET sets a Hermite Genz-Keister 18 rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    A nested family of rules for the Hermite integration problem
//    was produced by Genz and Keister.  The structure of the nested
//    family was denoted by 1+2+6+10+18, that is, it comprised rules 
//    of successive orders O = 1, 3, 9, 19, and 37.
//
//    The precisions of these rules are P = 1, 5, 15, 29, and 55.
//
//    Some of the data in this function was kindly supplied directly by
//    Alan Genz on 24 April 2011.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz, Bradley Keister,
//    Fully symmetric interpolatory rules for multiple integrals
//    over infinite regions with Gaussian weight,
//    Journal of Computational and Applied Mathematics,
//    Volume 71, 1996, pages 299-309
//
//    Florian Heiss, Viktor Winschel,
//    Likelihood approximation by numerical integration on sparse grids,
//    Journal of Econometrics,
//    Volume 144, 2008, pages 62-80.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 9, 19, or 37.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[ 0] =   0.0000000000000000E+00;

    w[ 0] =   1.7724538509055159E+00;
  }
  else if ( n == 3 )
  {
    x[ 0] =  -1.2247448713915889E+00;
    x[ 1] =   0.0000000000000000E+00;
    x[ 2] =   1.2247448713915889E+00;

    w[ 0] =   2.9540897515091930E-01;
    w[ 1] =   1.1816359006036772E+00;
    w[ 2] =   2.9540897515091930E-01;
  }
  else if ( n == 9 )
  {
    x[ 0] =  -2.9592107790638380E+00;
    x[ 1] =  -2.0232301911005157E+00;
    x[ 2] =  -1.2247448713915889E+00;
    x[ 3] =  -5.2403354748695763E-01;
    x[ 4] =   0.0000000000000000E+00;
    x[ 5] =   5.2403354748695763E-01;
    x[ 6] =   1.2247448713915889E+00;
    x[ 7] =   2.0232301911005157E+00;
    x[ 8] =   2.9592107790638380E+00;

    w[ 0] =   1.6708826306882348E-04;
    w[ 1] =   1.4173117873979098E-02;
    w[ 2] =   1.6811892894767771E-01;
    w[ 3] =   4.7869428549114124E-01;
    w[ 4] =   4.5014700975378197E-01;
    w[ 5] =   4.7869428549114124E-01;
    w[ 6] =   1.6811892894767771E-01;
    w[ 7] =   1.4173117873979098E-02;
    w[ 8] =   1.6708826306882348E-04;
  }
  else if ( n == 19 )
  {
    x[ 0] =  -4.4995993983103881E+00;
    x[ 1] =  -3.6677742159463378E+00;
    x[ 2] =  -2.9592107790638380E+00;
    x[ 3] =  -2.2665132620567876E+00;
    x[ 4] =  -2.0232301911005157E+00;
    x[ 5] =  -1.8357079751751868E+00;
    x[ 6] =  -1.2247448713915889E+00;
    x[ 7] =  -8.7004089535290285E-01;
    x[ 8] =  -5.2403354748695763E-01;
    x[ 9] =   0.0000000000000000E+00;
    x[10] =   5.2403354748695763E-01;
    x[11] =   8.7004089535290285E-01;
    x[12] =   1.2247448713915889E+00;
    x[13] =   1.8357079751751868E+00;
    x[14] =   2.0232301911005157E+00;
    x[15] =   2.2665132620567876E+00;
    x[16] =   2.9592107790638380E+00;
    x[17] =   3.6677742159463378E+00;
    x[18] =   4.4995993983103881E+00;

    w[ 0] =   1.5295717705322357E-09;
    w[ 1] =   1.0802767206624762E-06;
    w[ 2] =   1.0656589772852267E-04;
    w[ 3] =   5.1133174390883855E-03;
    w[ 4] =  -1.1232438489069229E-02;
    w[ 5] =   3.2055243099445879E-02;
    w[ 6] =   1.1360729895748269E-01;
    w[ 7] =   1.0838861955003017E-01;
    w[ 8] =   3.6924643368920851E-01;
    w[ 9] =   5.3788160700510168E-01;
    w[10] =   3.6924643368920851E-01;
    w[11] =   1.0838861955003017E-01;
    w[12] =   1.1360729895748269E-01;
    w[13] =   3.2055243099445879E-02;
    w[14] =  -1.1232438489069229E-02;
    w[15] =   5.1133174390883855E-03;
    w[16] =   1.0656589772852267E-04;
    w[17] =   1.0802767206624762E-06;
    w[18] =   1.5295717705322357E-09;
  }
  else if ( n == 37 )
  {
    x[ 0] =  -6.853200069757519;
    x[ 1] =  -6.124527854622158;
    x[ 2] =  -5.521865209868350;
    x[ 3] =  -4.986551454150765;
    x[ 4] =  -4.499599398310388;
    x[ 5] =  -4.057956316089741;
    x[ 6] =  -3.667774215946338;
    x[ 7] =  -3.315584617593290;
    x[ 8] =  -2.959210779063838;
    x[ 9] =  -2.597288631188366;
    x[10] =  -2.266513262056788;
    x[11] =  -2.023230191100516;
    x[12] =  -1.835707975175187;
    x[13] =  -1.561553427651873;
    x[14] =  -1.224744871391589;
    x[15] =  -0.870040895352903;
    x[16] =  -0.524033547486958;
    x[17] =  -0.214618180588171;
    x[18] =   0.000000000000000;
    x[19] =   0.214618180588171;
    x[20] =   0.524033547486958;
    x[21] =   0.870040895352903;
    x[22] =   1.224744871391589;
    x[23] =   1.561553427651873;
    x[24] =   1.835707975175187;
    x[25] =   2.023230191100516;
    x[26] =   2.266513262056788;
    x[27] =   2.597288631188366;
    x[28] =   2.959210779063838;
    x[29] =   3.315584617593290;
    x[30] =   3.667774215946338;
    x[31] =   4.057956316089741;
    x[32] =   4.499599398310388;
    x[33] =   4.986551454150765;
    x[34] =   5.521865209868350;
    x[35] =   6.124527854622158;
    x[36] =   6.853200069757519;

    w[ 0] = 0.19030350940130498E-20;
    w[ 1] = 0.187781893143728947E-16;
    w[ 2] = 0.182242751549129356E-13;
    w[ 3] = 0.45661763676186859E-11;
    w[ 4] = 0.422525843963111041E-09;
    w[ 5] = 0.16595448809389819E-07;
    w[ 6] = 0.295907520230744049E-06;
    w[ 7] = 0.330975870979203419E-05;
    w[ 8] = 0.32265185983739747E-04;
    w[ 9] = 0.234940366465975222E-03;
    w[10] = 0.985827582996483824E-03;
    w[11] = 0.176802225818295443E-02;
    w[12] = 0.43334988122723492E-02;
    w[13] = 0.15513109874859354E-01;
    w[14] = 0.442116442189845444E-01;
    w[15] = 0.937208280655245902E-01;
    w[16] = 0.143099302896833389E+00;
    w[17] = 0.147655710402686249E+00;
    w[18] = 0.968824552928425499E-01;
    w[19] = 0.147655710402686249E+00;
    w[20] = 0.143099302896833389E+00;
    w[21] = 0.937208280655245902E-01;
    w[22] = 0.442116442189845444E-01;
    w[23] = 0.15513109874859354E-01;
    w[24] = 0.43334988122723492E-02;
    w[25] = 0.176802225818295443E-02;
    w[26] = 0.985827582996483824E-03;
    w[27] = 0.234940366465975222E-03;
    w[28] = 0.32265185983739747E-04;
    w[29] = 0.330975870979203419E-05;
    w[30] = 0.295907520230744049E-06;
    w[31] = 0.16595448809389819E-07;
    w[32] = 0.422525843963111041E-09;
    w[33] = 0.45661763676186859E-11;
    w[34] = 0.182242751549129356E-13;
    w[35] = 0.187781893143728947E-16;
    w[36] = 0.19030350940130498E-20;
  }
  else
  {
    cerr << "\n";
    cerr << "HERMITE_GK18_SET - Fatal error!\n";
    cerr << "  Illegal input value of N.\n";
    cerr << "  N must be 1, 3, 9, 19, or 37.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void hermite_gk22_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK22_SET sets a Hermite Genz-Keister 22 rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    A nested family of rules for the Hermite integration problem
//    was produced by Genz and Keister.  The structure of the nested
//    family was denoted by 1+2+6+10+16, that is, it comprised rules 
//    of successive orders O = 1, 3, 9, 19, and 41.
//
//    The precisions of these rules are P = 1, 5, 15, 29, and 63.
//
//    Some of the data in this function was kindly supplied directly by
//    Alan Genz on 24 April 2011.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz, Bradley Keister,
//    Fully symmetric interpolatory rules for multiple integrals
//    over infinite regions with Gaussian weight,
//    Journal of Computational and Applied Mathematics,
//    Volume 71, 1996, pages 299-309
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 9, 19, or 41.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[ 0] =   0.0000000000000000E+00;

    w[ 0] =   1.7724538509055159E+00;
  }
  else if ( n == 3 )
  {
    x[ 0] =  -1.2247448713915889E+00;
    x[ 1] =   0.0000000000000000E+00;
    x[ 2] =   1.2247448713915889E+00;

    w[ 0] =   2.9540897515091930E-01;
    w[ 1] =   1.1816359006036772E+00;
    w[ 2] =   2.9540897515091930E-01;
  }
  else if ( n == 9 )
  {
    x[ 0] =  -2.9592107790638380E+00;
    x[ 1] =  -2.0232301911005157E+00;
    x[ 2] =  -1.2247448713915889E+00;
    x[ 3] =  -5.2403354748695763E-01;
    x[ 4] =   0.0000000000000000E+00;
    x[ 5] =   5.2403354748695763E-01;
    x[ 6] =   1.2247448713915889E+00;
    x[ 7] =   2.0232301911005157E+00;
    x[ 8] =   2.9592107790638380E+00;

    w[ 0] =   1.6708826306882348E-04;
    w[ 1] =   1.4173117873979098E-02;
    w[ 2] =   1.6811892894767771E-01;
    w[ 3] =   4.7869428549114124E-01;
    w[ 4] =   4.5014700975378197E-01;
    w[ 5] =   4.7869428549114124E-01;
    w[ 6] =   1.6811892894767771E-01;
    w[ 7] =   1.4173117873979098E-02;
    w[ 8] =   1.6708826306882348E-04;
  }
  else if ( n == 19 )
  {
    x[ 0] =  -4.4995993983103881E+00;
    x[ 1] =  -3.6677742159463378E+00;
    x[ 2] =  -2.9592107790638380E+00;
    x[ 3] =  -2.2665132620567876E+00;
    x[ 4] =  -2.0232301911005157E+00;
    x[ 5] =  -1.8357079751751868E+00;
    x[ 6] =  -1.2247448713915889E+00;
    x[ 7] =  -8.7004089535290285E-01;
    x[ 8] =  -5.2403354748695763E-01;
    x[ 9] =   0.0000000000000000E+00;
    x[10] =   5.2403354748695763E-01;
    x[11] =   8.7004089535290285E-01;
    x[12] =   1.2247448713915889E+00;
    x[13] =   1.8357079751751868E+00;
    x[14] =   2.0232301911005157E+00;
    x[15] =   2.2665132620567876E+00;
    x[16] =   2.9592107790638380E+00;
    x[17] =   3.6677742159463378E+00;
    x[18] =   4.4995993983103881E+00;

    w[ 0] =   1.5295717705322357E-09;
    w[ 1] =   1.0802767206624762E-06;
    w[ 2] =   1.0656589772852267E-04;
    w[ 3] =   5.1133174390883855E-03;
    w[ 4] =  -1.1232438489069229E-02;
    w[ 5] =   3.2055243099445879E-02;
    w[ 6] =   1.1360729895748269E-01;
    w[ 7] =   1.0838861955003017E-01;
    w[ 8] =   3.6924643368920851E-01;
    w[ 9] =   5.3788160700510168E-01;
    w[10] =   3.6924643368920851E-01;
    w[11] =   1.0838861955003017E-01;
    w[12] =   1.1360729895748269E-01;
    w[13] =   3.2055243099445879E-02;
    w[14] =  -1.1232438489069229E-02;
    w[15] =   5.1133174390883855E-03;
    w[16] =   1.0656589772852267E-04;
    w[17] =   1.0802767206624762E-06;
    w[18] =   1.5295717705322357E-09;
  }
  else if ( n == 41 )
  {
    x[ 0] =  -7.251792998192644;
    x[ 1] =  -6.547083258397540;
    x[ 2] =  -5.961461043404500;
    x[ 3] =  -5.437443360177798;
    x[ 4] =  -4.953574342912980;
    x[ 5] =  -4.4995993983103881;
    x[ 6] =  -4.070919267883068;
    x[ 7] =  -3.6677742159463378;
    x[ 8] =  -3.296114596212218;
    x[ 9] =  -2.9592107790638380;
    x[10] =  -2.630415236459871;
    x[11] =  -2.2665132620567876;
    x[12] =  -2.043834754429505;
    x[13] =  -2.0232301911005157;
    x[14] =  -1.8357079751751868;
    x[15] =  -1.585873011819188;
    x[16] =  -1.2247448713915889;
    x[17] =  -0.87004089535290285;
    x[18] =  -0.52403354748695763;
    x[19] =  -0.195324784415805;
    x[20] =   0.0000000000000000;
    x[21] =   0.195324784415805;
    x[22] =   0.52403354748695763;
    x[23] =   0.87004089535290285;
    x[24] =   1.2247448713915889;
    x[25] =   1.585873011819188;
    x[26] =   1.8357079751751868;
    x[27] =   2.0232301911005157;
    x[28] =   2.043834754429505;
    x[29] =   2.2665132620567876;
    x[30] =   2.630415236459871;
    x[31] =   2.9592107790638380;
    x[32] =   3.296114596212218;
    x[33] =   3.6677742159463378;
    x[34] =   4.070919267883068;
    x[35] =   4.4995993983103881;
    x[36] =   4.953574342912980;
    x[37] =   5.437443360177798;
    x[38] =   5.961461043404500;
    x[39] =   6.547083258397540;
    x[40] =   7.251792998192644;

    w[ 0] =   0.664195893812757801E-23;
    w[ 1] =   0.860427172512207236E-19;
    w[ 2] =   0.1140700785308509E-15;
    w[ 3] =   0.408820161202505983E-13;
    w[ 4] =   0.581803393170320419E-11;
    w[ 5] =   0.400784141604834759E-09;
    w[ 6] =   0.149158210417831408E-07;
    w[ 7] =   0.315372265852264871E-06;
    w[ 8] =   0.381182791749177506E-05;
    w[ 9] =   0.288976780274478689E-04;
    w[10] =   0.189010909805097887E-03;
    w[11] =   0.140697424065246825E-02;
    w[12] = - 0.144528422206988237E-01;
    w[13] =   0.178852543033699732E-01;
    w[14] =   0.705471110122962612E-03;
    w[15] =   0.165445526705860772E-01;
    w[16] =   0.45109010335859128E-01;
    w[17] =   0.928338228510111845E-01;
    w[18] =   0.145966293895926429E+00;
    w[19] =   0.165639740400529554E+00;
    w[20] =   0.562793426043218877E-01;
    w[21] =   0.165639740400529554E+00;
    w[22] =   0.145966293895926429E+00;
    w[23] =   0.928338228510111845E-01;
    w[24] =   0.45109010335859128E-01;
    w[25] =   0.165445526705860772E-01;
    w[26] =   0.705471110122962612E-03;
    w[27] =   0.178852543033699732E-01;
    w[28] = - 0.144528422206988237E-01;
    w[29] =   0.140697424065246825E-02;
    w[30] =   0.189010909805097887E-03;
    w[31] =   0.288976780274478689E-04;
    w[32] =   0.381182791749177506E-05;
    w[33] =   0.315372265852264871E-06;
    w[34] =   0.149158210417831408E-07;
    w[35] =   0.400784141604834759E-09;
    w[36] =   0.581803393170320419E-11;
    w[37] =   0.408820161202505983E-13;
    w[38] =   0.1140700785308509E-15;
    w[39] =   0.860427172512207236E-19;
    w[40] =   0.664195893812757801E-23;
  }
  else
  {
    cerr << "\n";
    cerr << "HERMITE_GK22_SET - Fatal error!\n";
    cerr << "  Illegal input value of N.\n";
    cerr << "  N must be 1, 3, 9, 19, or 41.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void hermite_gk24_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK24_SET sets a Hermite Genz-Keister 24 rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    A nested family of rules for the Hermite integration problem
//    was produced by Genz and Keister.  The structure of the nested
//    family was denoted by 1+2+6+10+16, that is, it comprised rules 
//    of successive orders O = 1, 3, 9, 19, and 43.
//
//    The precisions of these rules are P = 1, 5, 15, 29, and 67.
//
//    Some of the data in this function was kindly supplied directly by
//    Alan Genz on 24 April 2011.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz, Bradley Keister,
//    Fully symmetric interpolatory rules for multiple integrals
//    over infinite regions with Gaussian weight,
//    Journal of Computational and Applied Mathematics,
//    Volume 71, 1996, pages 299-309
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 9 19, or 43.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[ 0] =   0.0000000000000000E+00;

    w[ 0] =   1.7724538509055159E+00;
  }
  else if ( n == 3 )
  {
    x[ 0] =  -1.2247448713915889E+00;
    x[ 1] =   0.0000000000000000E+00;
    x[ 2] =   1.2247448713915889E+00;

    w[ 0] =   2.9540897515091930E-01;
    w[ 1] =   1.1816359006036772E+00;
    w[ 2] =   2.9540897515091930E-01;
  }
  else if ( n == 9 )
  {
    x[ 0] =  -2.9592107790638380E+00;
    x[ 1] =  -2.0232301911005157E+00;
    x[ 2] =  -1.2247448713915889E+00;
    x[ 3] =  -5.2403354748695763E-01;
    x[ 4] =   0.0000000000000000E+00;
    x[ 5] =   5.2403354748695763E-01;
    x[ 6] =   1.2247448713915889E+00;
    x[ 7] =   2.0232301911005157E+00;
    x[ 8] =   2.9592107790638380E+00;

    w[ 0] =   1.6708826306882348E-04;
    w[ 1] =   1.4173117873979098E-02;
    w[ 2] =   1.6811892894767771E-01;
    w[ 3] =   4.7869428549114124E-01;
    w[ 4] =   4.5014700975378197E-01;
    w[ 5] =   4.7869428549114124E-01;
    w[ 6] =   1.6811892894767771E-01;
    w[ 7] =   1.4173117873979098E-02;
    w[ 8] =   1.6708826306882348E-04;
  }
  else if ( n == 19 )
  {
    x[ 0] =  -4.4995993983103881E+00;
    x[ 1] =  -3.6677742159463378E+00;
    x[ 2] =  -2.9592107790638380E+00;
    x[ 3] =  -2.2665132620567876E+00;
    x[ 4] =  -2.0232301911005157E+00;
    x[ 5] =  -1.8357079751751868E+00;
    x[ 6] =  -1.2247448713915889E+00;
    x[ 7] =  -8.7004089535290285E-01;
    x[ 8] =  -5.2403354748695763E-01;
    x[ 9] =   0.0000000000000000E+00;
    x[10] =   5.2403354748695763E-01;
    x[11] =   8.7004089535290285E-01;
    x[12] =   1.2247448713915889E+00;
    x[13] =   1.8357079751751868E+00;
    x[14] =   2.0232301911005157E+00;
    x[15] =   2.2665132620567876E+00;
    x[16] =   2.9592107790638380E+00;
    x[17] =   3.6677742159463378E+00;
    x[18] =   4.4995993983103881E+00;

    w[ 0] =   1.5295717705322357E-09;
    w[ 1] =   1.0802767206624762E-06;
    w[ 2] =   1.0656589772852267E-04;
    w[ 3] =   5.1133174390883855E-03;
    w[ 4] =  -1.1232438489069229E-02;
    w[ 5] =   3.2055243099445879E-02;
    w[ 6] =   1.1360729895748269E-01;
    w[ 7] =   1.0838861955003017E-01;
    w[ 8] =   3.6924643368920851E-01;
    w[ 9] =   5.3788160700510168E-01;
    w[10] =   3.6924643368920851E-01;
    w[11] =   1.0838861955003017E-01;
    w[12] =   1.1360729895748269E-01;
    w[13] =   3.2055243099445879E-02;
    w[14] =  -1.1232438489069229E-02;
    w[15] =   5.1133174390883855E-03;
    w[16] =   1.0656589772852267E-04;
    w[17] =   1.0802767206624762E-06;
    w[18] =   1.5295717705322357E-09;
  }
  else if ( n == 43 )
  {
    x[ 0] = -10.167574994881873;
    x[ 1] =  -7.231746029072501;
    x[ 2] =  -6.535398426382995;
    x[ 3] =  -5.954781975039809;
    x[ 4] =  -5.434053000365068;
    x[ 5] =  -4.952329763008589;
    x[ 6] =  -4.4995993983103881;
    x[ 7] =  -4.071335874253583;
    x[ 8] =  -3.6677742159463378;
    x[ 9] =  -3.295265921534226;
    x[10] =  -2.9592107790638380;
    x[11] =  -2.633356763661946;
    x[12] =  -2.2665132620567876;
    x[13] =  -2.089340389294661;
    x[14] =  -2.0232301911005157;
    x[15] =  -1.8357079751751868;
    x[16] =  -1.583643465293944;
    x[17] =  -1.2247448713915889;
    x[18] =  -0.87004089535290285;
    x[19] =  -0.52403354748695763;
    x[20] =  -0.196029453662011;
    x[21] =   0.0000000000000000;
    x[22] =   0.196029453662011;
    x[23] =   0.52403354748695763;
    x[24] =   0.87004089535290285;
    x[25] =   1.2247448713915889;
    x[26] =   1.583643465293944;
    x[27] =   1.8357079751751868;
    x[28] =   2.0232301911005157;
    x[29] =   2.089340389294661;
    x[30] =   2.2665132620567876;
    x[31] =   2.633356763661946;
    x[32] =   2.9592107790638380;
    x[33] =   3.295265921534226;
    x[34] =   3.6677742159463378;
    x[35] =   4.071335874253583;
    x[36] =   4.4995993983103881;
    x[37] =   4.952329763008589;
    x[38] =   5.434053000365068;
    x[39] =   5.954781975039809;
    x[40] =   6.535398426382995;
    x[41] =   7.231746029072501;
    x[42] =  10.167574994881873;

    w[ 0] =   0.546191947478318097E-37;
    w[ 1] =   0.87544909871323873E-23;
    w[ 2] =   0.992619971560149097E-19;
    w[ 3] =   0.122619614947864357E-15;
    w[ 4] =   0.421921851448196032E-13;
    w[ 5] =   0.586915885251734856E-11;
    w[ 6] =   0.400030575425776948E-09;
    w[ 7] =   0.148653643571796457E-07;
    w[ 8] =   0.316018363221289247E-06;
    w[ 9] =   0.383880761947398577E-05;
    w[10] =   0.286802318064777813E-04;
    w[11] =   0.184789465688357423E-03;
    w[12] =   0.150909333211638847E-02;
    w[13] = - 0.38799558623877157E-02;
    w[14] =   0.67354758901013295E-02;
    w[15] =   0.139966252291568061E-02;
    w[16] =   0.163616873493832402E-01;
    w[17] =   0.450612329041864976E-01;
    w[18] =   0.928711584442575456E-01;
    w[19] =   0.145863292632147353E+00;
    w[20] =   0.164880913687436689E+00;
    w[21] =   0.579595986101181095E-01;
    w[22] =   0.164880913687436689E+00;
    w[23] =   0.145863292632147353E+00;
    w[24] =   0.928711584442575456E-01;
    w[25] =   0.450612329041864976E-01;
    w[26] =   0.163616873493832402E-01;
    w[27] =   0.139966252291568061E-02;
    w[28] =   0.67354758901013295E-02;
    w[29] = - 0.38799558623877157E-02;
    w[30] =   0.150909333211638847E-02;
    w[31] =   0.184789465688357423E-03;
    w[32] =   0.286802318064777813E-04;
    w[33] =   0.383880761947398577E-05;
    w[34] =   0.316018363221289247E-06;
    w[35] =   0.148653643571796457E-07;
    w[36] =   0.400030575425776948E-09;
    w[37] =   0.586915885251734856E-11;
    w[38] =   0.421921851448196032E-13;
    w[39] =   0.122619614947864357E-15;
    w[40] =   0.992619971560149097E-19;
    w[41] =   0.87544909871323873E-23;
    w[42] =   0.546191947478318097E-37;
  }
  else
  {
    cerr << "\n";
    cerr << "HERMITE_GK24_SET - Fatal error!\n";
    cerr << "  Illegal input value of N.\n";
    cerr << "  N must be 1, 3, 9, 19, or 43.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double hermite_integral ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_INTEGRAL evaluates a monomial Hermite integral.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -oo < x < +oo ) x^n exp(-x^2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the integral.  
//    0 <= N.
//
//    Output, double VALUE, the value of the integral.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( n < 0 )
  {
    value = - r8_huge ( );
  }
  else if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial2 ( n - 1 ) * sqrt ( pi ) / pow ( 2.0, n / 2 );
  }

  return value;
}
//****************************************************************************80

void hermite_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_SET sets abscissas and weights for Hermite quadrature.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) ).
//
//    Mathematica can numerically estimate the abscissas 
//    of order N to P digits by the command:
//
//      NSolve [ HermiteH [ n, x ] == 0, x, p ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20, or 31/32/33, 63/64/65, 127/128/129.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[0] = 0.0;

    w[0] = 1.77245385090551602729816748334;
  }
  else if ( n == 2 )
  {
    x[0] = - 0.707106781186547524400844362105;
    x[1] =   0.707106781186547524400844362105;

    w[0] = 0.886226925452758013649083741671;
    w[1] = 0.886226925452758013649083741671;
  }
  else if ( n == 3 )
  {
    x[0] = - 0.122474487139158904909864203735E+01;
    x[1] =   0.0;
    x[2] =   0.122474487139158904909864203735E+01;

    w[0] = 0.295408975150919337883027913890;
    w[1] = 0.118163590060367735153211165556E+01;
    w[2] = 0.295408975150919337883027913890;
  }
  else if ( n == 4 )
  {
    x[0] = - 0.165068012388578455588334111112E+01;
    x[1] = - 0.524647623275290317884060253835;
    x[2] =   0.524647623275290317884060253835;
    x[3] =   0.165068012388578455588334111112E+01;

    w[0] = 0.813128354472451771430345571899E-01;
    w[1] = 0.804914090005512836506049184481;
    w[2] = 0.804914090005512836506049184481;
    w[3] = 0.813128354472451771430345571899E-01;
  }
  else if ( n == 5 )
  {
    x[0] = - 0.202018287045608563292872408814E+01;
    x[1] = - 0.958572464613818507112770593893;
    x[2] =   0.0;
    x[3] =   0.958572464613818507112770593893;
    x[4] =   0.202018287045608563292872408814E+01;

    w[0] = 0.199532420590459132077434585942E-01;
    w[1] = 0.393619323152241159828495620852;
    w[2] = 0.945308720482941881225689324449;
    w[3] = 0.393619323152241159828495620852;
    w[4] = 0.199532420590459132077434585942E-01;
  }
  else if ( n == 6 )
  {
    x[0] = - 0.235060497367449222283392198706E+01;
    x[1] = - 0.133584907401369694971489528297E+01;
    x[2] = - 0.436077411927616508679215948251;
    x[3] =   0.436077411927616508679215948251;
    x[4] =   0.133584907401369694971489528297E+01;
    x[5] =   0.235060497367449222283392198706E+01;

    w[0] = 0.453000990550884564085747256463E-02;
    w[1] = 0.157067320322856643916311563508;
    w[2] = 0.724629595224392524091914705598;
    w[3] = 0.724629595224392524091914705598;
    w[4] = 0.157067320322856643916311563508;
    w[5] = 0.453000990550884564085747256463E-02;
  }
  else if ( n == 7 )
  {
    x[0] = - 0.265196135683523349244708200652E+01;
    x[1] = - 0.167355162876747144503180139830E+01;
    x[2] = - 0.816287882858964663038710959027;
    x[3] =   0.0;
    x[4] =   0.816287882858964663038710959027;
    x[5] =   0.167355162876747144503180139830E+01;
    x[6] =   0.265196135683523349244708200652E+01;

    w[0] = 0.971781245099519154149424255939E-03;
    w[1] = 0.545155828191270305921785688417E-01;
    w[2] = 0.425607252610127800520317466666;
    w[3] = 0.810264617556807326764876563813;
    w[4] = 0.425607252610127800520317466666;
    w[5] = 0.545155828191270305921785688417E-01;
    w[6] = 0.971781245099519154149424255939E-03;
  }
  else if ( n == 8 )
  {
    x[0] = - 0.293063742025724401922350270524E+01;
    x[1] = - 0.198165675669584292585463063977E+01;
    x[2] = - 0.115719371244678019472076577906E+01;
    x[3] = - 0.381186990207322116854718885584;
    x[4] =   0.381186990207322116854718885584;
    x[5] =   0.115719371244678019472076577906E+01;
    x[6] =   0.198165675669584292585463063977E+01;
    x[7] =   0.293063742025724401922350270524E+01;

    w[0] = 0.199604072211367619206090452544E-03;
    w[1] = 0.170779830074134754562030564364E-01;
    w[2] = 0.207802325814891879543258620286;
    w[3] = 0.661147012558241291030415974496;
    w[4] = 0.661147012558241291030415974496;
    w[5] = 0.207802325814891879543258620286;
    w[6] = 0.170779830074134754562030564364E-01;
    w[7] = 0.199604072211367619206090452544E-03;
  }
  else if ( n == 9 )
  {
    x[0] = - 0.319099320178152760723004779538E+01;
    x[1] = - 0.226658058453184311180209693284E+01;
    x[2] = - 0.146855328921666793166701573925E+01;
    x[3] = - 0.723551018752837573322639864579;
    x[4] =   0.0;
    x[5] =   0.723551018752837573322639864579;
    x[6] =   0.146855328921666793166701573925E+01;
    x[7] =   0.226658058453184311180209693284E+01;
    x[8] =   0.319099320178152760723004779538E+01;

    w[0] = 0.396069772632643819045862946425E-04;
    w[1] = 0.494362427553694721722456597763E-02;
    w[2] = 0.884745273943765732879751147476E-01;
    w[3] = 0.432651559002555750199812112956;
    w[4] = 0.720235215606050957124334723389;
    w[5] = 0.432651559002555750199812112956;
    w[6] = 0.884745273943765732879751147476E-01;
    w[7] = 0.494362427553694721722456597763E-02;
    w[8] = 0.396069772632643819045862946425E-04;
  }
  else if ( n == 10 )
  {
    x[0] =  - 0.343615911883773760332672549432E+01;
    x[1] =  - 0.253273167423278979640896079775E+01;
    x[2] =  - 0.175668364929988177345140122011E+01;
    x[3] =  - 0.103661082978951365417749191676E+01;
    x[4] =  - 0.342901327223704608789165025557;
    x[5] =    0.342901327223704608789165025557;
    x[6] =    0.103661082978951365417749191676E+01;
    x[7] =    0.175668364929988177345140122011E+01;
    x[8] =    0.253273167423278979640896079775E+01;
    x[9] =   0.343615911883773760332672549432E+01;

    w[0] =  0.764043285523262062915936785960E-05;
    w[1] =  0.134364574678123269220156558585E-02;
    w[2] =  0.338743944554810631361647312776E-01;
    w[3] =  0.240138611082314686416523295006;
    w[4] =  0.610862633735325798783564990433;
    w[5] =  0.610862633735325798783564990433;
    w[6] =  0.240138611082314686416523295006;
    w[7] =  0.338743944554810631361647312776E-01;
    w[8] =  0.134364574678123269220156558585E-02;
    w[9] = 0.764043285523262062915936785960E-05;
  }
  else if ( n == 11 )
  {
    x[0] =  - 0.366847084655958251845837146485E+01;
    x[1] =  - 0.278329009978165177083671870152E+01;
    x[2] =  - 0.202594801582575533516591283121E+01;
    x[3] =  - 0.132655708449493285594973473558E+01;
    x[4] =  - 0.656809566882099765024611575383;
    x[5] =    0.0;
    x[6] =    0.656809566882099765024611575383;
    x[7] =    0.132655708449493285594973473558E+01;
    x[8] =    0.202594801582575533516591283121E+01;
    x[9] =   0.278329009978165177083671870152E+01;
    x[10] =   0.366847084655958251845837146485E+01;

    w[0] =  0.143956039371425822033088366032E-05;
    w[1] =  0.346819466323345510643413772940E-03;
    w[2] =  0.119113954449115324503874202916E-01;
    w[3] =  0.117227875167708503381788649308;
    w[4] =  0.429359752356125028446073598601;
    w[5] =  0.654759286914591779203940657627;
    w[6] =  0.429359752356125028446073598601;
    w[7] =  0.117227875167708503381788649308;
    w[8] =  0.119113954449115324503874202916E-01;
    w[9] = 0.346819466323345510643413772940E-03;
    w[10] = 0.143956039371425822033088366032E-05;
  }
  else if ( n == 12 )
  {
    x[0] =  - 0.388972489786978191927164274724E+01;
    x[1] =  - 0.302063702512088977171067937518E+01;
    x[2] =  - 0.227950708050105990018772856942E+01;
    x[3] =  - 0.159768263515260479670966277090E+01;
    x[4] =  - 0.947788391240163743704578131060;
    x[5] =  - 0.314240376254359111276611634095;
    x[6] =    0.314240376254359111276611634095;
    x[7] =    0.947788391240163743704578131060;
    x[8] =    0.159768263515260479670966277090E+01;
    x[9] =   0.227950708050105990018772856942E+01;
    x[10] =   0.302063702512088977171067937518E+01;
    x[11] =   0.388972489786978191927164274724E+01;

    w[0] =  0.265855168435630160602311400877E-06;
    w[1] =  0.857368704358785865456906323153E-04;
    w[2] =  0.390539058462906185999438432620E-02;
    w[3] =  0.516079856158839299918734423606E-01;
    w[4] =  0.260492310264161129233396139765;
    w[5] =  0.570135236262479578347113482275;
    w[6] =  0.570135236262479578347113482275;
    w[7] =  0.260492310264161129233396139765;
    w[8] =  0.516079856158839299918734423606E-01;
    w[9] = 0.390539058462906185999438432620E-02;
    w[10] = 0.857368704358785865456906323153E-04;
    w[11] = 0.265855168435630160602311400877E-06;
  }
  else if ( n == 13 )
  {
    x[0] =  - 0.410133759617863964117891508007E+01;
    x[1] =  - 0.324660897837240998812205115236E+01;
    x[2] =  - 0.251973568567823788343040913628E+01;
    x[3] =  - 0.185310765160151214200350644316E+01;
    x[4] =  - 0.122005503659074842622205526637E+01;
    x[5] =  - 0.605763879171060113080537108602;
    x[6] =    0.0;
    x[7] =    0.605763879171060113080537108602;
    x[8] =    0.122005503659074842622205526637E+01;
    x[9] =   0.185310765160151214200350644316E+01;
    x[10] =   0.251973568567823788343040913628E+01;
    x[11] =   0.324660897837240998812205115236E+01;
    x[12] =   0.410133759617863964117891508007E+01;

    w[0] =  0.482573185007313108834997332342E-07;
    w[1] =  0.204303604027070731248669432937E-04;
    w[2] =  0.120745999271938594730924899224E-02;
    w[3] =  0.208627752961699392166033805050E-01;
    w[4] =  0.140323320687023437762792268873;
    w[5] =  0.421616296898543221746893558568;
    w[6] =  0.604393187921161642342099068579;
    w[7] =  0.421616296898543221746893558568;
    w[8] =  0.140323320687023437762792268873;
    w[9] = 0.208627752961699392166033805050E-01;
    w[10] = 0.120745999271938594730924899224E-02;
    w[11] = 0.204303604027070731248669432937E-04;
    w[12] = 0.482573185007313108834997332342E-07;
  }
  else if ( n == 14 )
  {
    x[0] =  - 0.430444857047363181262129810037E+01;
    x[1] =  - 0.346265693360227055020891736115E+01;
    x[2] =  - 0.274847072498540256862499852415E+01;
    x[3] =  - 0.209518325850771681573497272630E+01;
    x[4] =  - 0.147668273114114087058350654421E+01;
    x[5] =  - 0.878713787329399416114679311861;
    x[6] =  - 0.291745510672562078446113075799;
    x[7] =    0.291745510672562078446113075799;
    x[8] =    0.878713787329399416114679311861;
    x[9] =   0.147668273114114087058350654421E+01;
    x[10] =   0.209518325850771681573497272630E+01;
    x[11] =   0.274847072498540256862499852415E+01;
    x[12] =   0.346265693360227055020891736115E+01;
    x[13] =   0.430444857047363181262129810037E+01;

    w[0] =  0.862859116812515794532041783429E-08;
    w[1] =  0.471648435501891674887688950105E-05;
    w[2] =  0.355092613551923610483661076691E-03;
    w[3] =  0.785005472645794431048644334608E-02;
    w[4] =  0.685055342234652055387163312367E-01;
    w[5] =  0.273105609064246603352569187026;
    w[6] =  0.536405909712090149794921296776;
    w[7] =  0.536405909712090149794921296776;
    w[8] =  0.273105609064246603352569187026;
    w[9] = 0.685055342234652055387163312367E-01;
    w[10] = 0.785005472645794431048644334608E-02;
    w[11] = 0.355092613551923610483661076691E-03;
    w[12] = 0.471648435501891674887688950105E-05;
    w[13] = 0.862859116812515794532041783429E-08;
  }
  else if ( n == 15 )
  {
    x[0] =  - 0.449999070730939155366438053053E+01;
    x[1] =  - 0.366995037340445253472922383312E+01;
    x[2] =  - 0.296716692790560324848896036355E+01;
    x[3] =  - 0.232573248617385774545404479449E+01;
    x[4] =  - 0.171999257518648893241583152515E+01;
    x[5] =  - 0.113611558521092066631913490556E+01;
    x[6] =  - 0.565069583255575748526020337198;
    x[7] =    0.0;
    x[8] =    0.565069583255575748526020337198;
    x[9] =   0.113611558521092066631913490556E+01;
    x[10] =   0.171999257518648893241583152515E+01;
    x[11] =   0.232573248617385774545404479449E+01;
    x[12] =   0.296716692790560324848896036355E+01;
    x[13] =   0.366995037340445253472922383312E+01;
    x[14] =   0.449999070730939155366438053053E+01;

    w[0] =  0.152247580425351702016062666965E-08;
    w[1] =  0.105911554771106663577520791055E-05;
    w[2] =  0.100004441232499868127296736177E-03;
    w[3] =  0.277806884291277589607887049229E-02;
    w[4] =  0.307800338725460822286814158758E-01;
    w[5] =  0.158488915795935746883839384960;
    w[6] =  0.412028687498898627025891079568;
    w[7] =  0.564100308726417532852625797340;
    w[8] =  0.412028687498898627025891079568;
    w[9] = 0.158488915795935746883839384960;
    w[10] = 0.307800338725460822286814158758E-01;
    w[11] = 0.277806884291277589607887049229E-02;
    w[12] = 0.100004441232499868127296736177E-03;
    w[13] = 0.105911554771106663577520791055E-05;
    w[14] = 0.152247580425351702016062666965E-08;
  }
  else if ( n == 16 )
  {
    x[0] =  - 0.468873893930581836468849864875E+01;
    x[1] =  - 0.386944790486012269871942409801E+01;
    x[2] =  - 0.317699916197995602681399455926E+01;
    x[3] =  - 0.254620215784748136215932870545E+01;
    x[4] =  - 0.195178799091625397743465541496E+01;
    x[5] =  - 0.138025853919888079637208966969E+01;
    x[6] =  - 0.822951449144655892582454496734;
    x[7] =  - 0.273481046138152452158280401965;
    x[8] =    0.273481046138152452158280401965;
    x[9] =   0.822951449144655892582454496734;
    x[10] =   0.138025853919888079637208966969E+01;
    x[11] =   0.195178799091625397743465541496E+01;
    x[12] =   0.254620215784748136215932870545E+01;
    x[13] =   0.317699916197995602681399455926E+01;
    x[14] =   0.386944790486012269871942409801E+01;
    x[15] =   0.468873893930581836468849864875E+01;

    w[0] =  0.265480747401118224470926366050E-09;
    w[1] =  0.232098084486521065338749423185E-06;
    w[2] =  0.271186009253788151201891432244E-04;
    w[3] =  0.932284008624180529914277305537E-03;
    w[4] =  0.128803115355099736834642999312E-01;
    w[5] =  0.838100413989858294154207349001E-01;
    w[6] =  0.280647458528533675369463335380;
    w[7] =  0.507929479016613741913517341791;
    w[8] =  0.507929479016613741913517341791;
    w[9] = 0.280647458528533675369463335380;
    w[10] = 0.838100413989858294154207349001E-01;
    w[11] = 0.128803115355099736834642999312E-01;
    w[12] = 0.932284008624180529914277305537E-03;
    w[13] = 0.271186009253788151201891432244E-04;
    w[14] = 0.232098084486521065338749423185E-06;
    w[15] = 0.265480747401118224470926366050E-09;
  }
  else if ( n == 17 )
  {
    x[0] =  - 0.487134519367440308834927655662E+01;
    x[1] =  - 0.406194667587547430689245559698E+01;
    x[2] =  - 0.337893209114149408338327069289E+01;
    x[3] =  - 0.275776291570388873092640349574E+01;
    x[4] =  - 0.217350282666662081927537907149E+01;
    x[5] =  - 0.161292431422123133311288254454E+01;
    x[6] =  - 0.106764872574345055363045773799E+01;
    x[7] =  - 0.531633001342654731349086553718;
    x[8] =    0.0;
    x[9] =   0.531633001342654731349086553718;
    x[10] =   0.106764872574345055363045773799E+01;
    x[11] =   0.161292431422123133311288254454E+01;
    x[12] =   0.217350282666662081927537907149E+01;
    x[13] =   0.275776291570388873092640349574E+01;
    x[14] =   0.337893209114149408338327069289E+01;
    x[15] =   0.406194667587547430689245559698E+01;
    x[16] =   0.487134519367440308834927655662E+01;

    w[0] =  0.458057893079863330580889281222E-10;
    w[1] =  0.497707898163079405227863353715E-07;
    w[2] =  0.711228914002130958353327376218E-05;
    w[3] =  0.298643286697753041151336643059E-03;
    w[4] =  0.506734995762753791170069495879E-02;
    w[5] =  0.409200341495762798094994877854E-01;
    w[6] =  0.172648297670097079217645196219;
    w[7] =  0.401826469470411956577635085257;
    w[8] =  0.530917937624863560331883103379;
    w[9] = 0.401826469470411956577635085257;
    w[10] = 0.172648297670097079217645196219;
    w[11] = 0.409200341495762798094994877854E-01;
    w[12] = 0.506734995762753791170069495879E-02;
    w[13] = 0.298643286697753041151336643059E-03;
    w[14] = 0.711228914002130958353327376218E-05;
    w[15] = 0.497707898163079405227863353715E-07;
    w[16] = 0.458057893079863330580889281222E-10;
  }
  else if ( n == 18 )
  {
    x[0] =  - 0.504836400887446676837203757885E+01;
    x[1] =  - 0.424811787356812646302342016090E+01;
    x[2] =  - 0.357376906848626607950067599377E+01;
    x[3] =  - 0.296137750553160684477863254906E+01;
    x[4] =  - 0.238629908916668600026459301424E+01;
    x[5] =  - 0.183553160426162889225383944409E+01;
    x[6] =  - 0.130092085838961736566626555439E+01;
    x[7] =  - 0.776682919267411661316659462284;
    x[8] =  - 0.258267750519096759258116098711;
    x[9] =   0.258267750519096759258116098711;
    x[10] =   0.776682919267411661316659462284;
    x[11] =   0.130092085838961736566626555439E+01;
    x[12] =   0.183553160426162889225383944409E+01;
    x[13] =   0.238629908916668600026459301424E+01;
    x[14] =   0.296137750553160684477863254906E+01;
    x[15] =   0.357376906848626607950067599377E+01;
    x[16] =   0.424811787356812646302342016090E+01;
    x[17] =   0.504836400887446676837203757885E+01;

    w[0] =  0.782819977211589102925147471012E-11;
    w[1] =  0.104672057957920824443559608435E-07;
    w[2] =  0.181065448109343040959702385911E-05;
    w[3] =  0.918112686792940352914675407371E-04;
    w[4] =  0.188852263026841789438175325426E-02;
    w[5] =  0.186400423875446519219315221973E-01;
    w[6] =  0.973017476413154293308537234155E-01;
    w[7] =  0.284807285669979578595606820713;
    w[8] =  0.483495694725455552876410522141;
    w[9] = 0.483495694725455552876410522141;
    w[10] = 0.284807285669979578595606820713;
    w[11] = 0.973017476413154293308537234155E-01;
    w[12] = 0.186400423875446519219315221973E-01;
    w[13] = 0.188852263026841789438175325426E-02;
    w[14] = 0.918112686792940352914675407371E-04;
    w[15] = 0.181065448109343040959702385911E-05;
    w[16] = 0.104672057957920824443559608435E-07;
    w[17] = 0.782819977211589102925147471012E-11;
  }
  else if ( n == 19 )
  {
    x[0] =  - 0.522027169053748216460967142500E+01;
    x[1] =  - 0.442853280660377943723498532226E+01;
    x[2] =  - 0.376218735196402009751489394104E+01;
    x[3] =  - 0.315784881834760228184318034120E+01;
    x[4] =  - 0.259113378979454256492128084112E+01;
    x[5] =  - 0.204923170985061937575050838669E+01;
    x[6] =  - 0.152417061939353303183354859367E+01;
    x[7] =  - 0.101036838713431135136859873726E+01;
    x[8] =  - 0.503520163423888209373811765050;
    x[9] =   0.0;
    x[10] =   0.503520163423888209373811765050;
    x[11] =   0.101036838713431135136859873726E+01;
    x[12] =   0.152417061939353303183354859367E+01;
    x[13] =   0.204923170985061937575050838669E+01;
    x[14] =   0.259113378979454256492128084112E+01;
    x[15] =   0.315784881834760228184318034120E+01;
    x[16] =   0.376218735196402009751489394104E+01;
    x[17] =   0.442853280660377943723498532226E+01;
    x[18] =   0.522027169053748216460967142500E+01;

    w[0] =  0.132629709449851575185289154385E-11;
    w[1] =  0.216305100986355475019693077221E-08;
    w[2] =  0.448824314722312295179447915594E-06;
    w[3] =  0.272091977631616257711941025214E-04;
    w[4] =  0.670877521407181106194696282100E-03;
    w[5] =  0.798886677772299020922211491861E-02;
    w[6] =  0.508103869090520673569908110358E-01;
    w[7] =  0.183632701306997074156148485766;
    w[8] =  0.391608988613030244504042313621;
    w[9] = 0.502974888276186530840731361096;
    w[10] = 0.391608988613030244504042313621;
    w[11] = 0.183632701306997074156148485766;
    w[12] = 0.508103869090520673569908110358E-01;
    w[13] = 0.798886677772299020922211491861E-02;
    w[14] = 0.670877521407181106194696282100E-03;
    w[15] = 0.272091977631616257711941025214E-04;
    w[16] = 0.448824314722312295179447915594E-06;
    w[17] = 0.216305100986355475019693077221E-08;
    w[18] = 0.132629709449851575185289154385E-11;
  }
  else if ( n == 20 )
  {
    x[0] =  - 0.538748089001123286201690041068E+01;
    x[1] =  - 0.460368244955074427307767524898E+01;
    x[2] =  - 0.394476404011562521037562880052E+01;
    x[3] =  - 0.334785456738321632691492452300E+01;
    x[4] =  - 0.278880605842813048052503375640E+01;
    x[5] =  - 0.225497400208927552308233334473E+01;
    x[6] =  - 0.173853771211658620678086566214E+01;
    x[7] =  - 0.123407621539532300788581834696E+01;
    x[8] =  - 0.737473728545394358705605144252;
    x[9] = - 0.245340708300901249903836530634;
    x[10] =   0.245340708300901249903836530634;
    x[11] =   0.737473728545394358705605144252;
    x[12] =   0.123407621539532300788581834696E+01;
    x[13] =   0.173853771211658620678086566214E+01;
    x[14] =   0.225497400208927552308233334473E+01;
    x[15] =   0.278880605842813048052503375640E+01;
    x[16] =   0.334785456738321632691492452300E+01;
    x[17] =   0.394476404011562521037562880052E+01;
    x[18] =   0.460368244955074427307767524898E+01;
    x[19] =   0.538748089001123286201690041068E+01;

    w[0] =  0.222939364553415129252250061603E-12;
    w[1] =  0.439934099227318055362885145547E-09;
    w[2] =  0.108606937076928169399952456345E-06;
    w[3] =  0.780255647853206369414599199965E-05;
    w[4] =  0.228338636016353967257145917963E-03;
    w[5] =  0.324377334223786183218324713235E-02;
    w[6] =  0.248105208874636108821649525589E-01;
    w[7] =  0.109017206020023320013755033535;
    w[8] =  0.286675505362834129719659706228;
    w[9] = 0.462243669600610089650328639861;
    w[10] = 0.462243669600610089650328639861;
    w[11] = 0.286675505362834129719659706228;
    w[12] = 0.109017206020023320013755033535;
    w[13] = 0.248105208874636108821649525589E-01;
    w[14] = 0.324377334223786183218324713235E-02;
    w[15] = 0.228338636016353967257145917963E-03;
    w[16] = 0.780255647853206369414599199965E-05;
    w[17] = 0.108606937076928169399952456345E-06;
    w[18] = 0.439934099227318055362885145547E-09;
    w[19] = 0.222939364553415129252250061603E-12;
  }
  else if ( n == 30 )
  {
    x[ 0] =   -6.86334529352989158106110835756;
    x[ 1] =   -6.13827922012393462039499237854;
    x[ 2] =   -5.53314715156749572511833355558;
    x[ 3] =   -4.98891896858994394448649710633;
    x[ 4] =   -4.48305535709251834188703761971;
    x[ 5] =   -4.00390860386122881522787601332;
    x[ 6] =   -3.54444387315534988692540090217;
    x[ 7] =   -3.09997052958644174868873332237;
    x[ 8] =   -2.66713212453561720057110646422;
    x[9] =   -2.24339146776150407247297999483;
    x[10] =   -1.82674114360368803883588048351;
    x[11] =   -1.41552780019818851194072510555;
    x[12] =   -1.00833827104672346180498960870;
    x[13] =   -0.603921058625552307778155678757;
    x[14] =   -0.201128576548871485545763013244;
    x[15] =    0.201128576548871485545763013244;
    x[16] =    0.603921058625552307778155678757;
    x[17] =    1.00833827104672346180498960870;
    x[18] =    1.41552780019818851194072510555;
    x[19] =    1.82674114360368803883588048351;
    x[20] =    2.24339146776150407247297999483;
    x[21] =    2.66713212453561720057110646422;
    x[22] =    3.09997052958644174868873332237;
    x[23] =    3.54444387315534988692540090217;
    x[24] =    4.00390860386122881522787601332;
    x[25] =    4.48305535709251834188703761971;
    x[26] =    4.98891896858994394448649710633;
    x[27] =    5.53314715156749572511833355558;
    x[28] =    6.13827922012393462039499237854;
    x[29] =    6.86334529352989158106110835756;

    w[ 0] =   0.290825470013122622941102747365E-20;
    w[ 1] =   0.281033360275090370876277491534E-16;
    w[ 2] =   0.287860708054870606219239791142E-13;
    w[ 3] =   0.810618629746304420399344796173E-11;
    w[ 4] =   0.917858042437852820850075742492E-09;
    w[ 5] =   0.510852245077594627738963204403E-07;
    w[ 6] =   0.157909488732471028834638794022E-05;
    w[ 7] =   0.293872522892298764150118423412E-04;
    w[ 8] =   0.348310124318685523420995323183E-03;
    w[9] =   0.273792247306765846298942568953E-02;
    w[10] =   0.147038297048266835152773557787E-01;
    w[11] =   0.551441768702342511680754948183E-01;
    w[12] =   0.146735847540890099751693643152;
    w[13] =   0.280130930839212667413493211293;
    w[14] =   0.386394889541813862555601849165;
    w[15] =   0.386394889541813862555601849165;
    w[16] =   0.280130930839212667413493211293;
    w[17] =   0.146735847540890099751693643152;
    w[18] =   0.551441768702342511680754948183E-01;
    w[19] =   0.147038297048266835152773557787E-01;
    w[20] =   0.273792247306765846298942568953E-02;
    w[21] =   0.348310124318685523420995323183E-03;
    w[22] =   0.293872522892298764150118423412E-04;
    w[23] =   0.157909488732471028834638794022E-05;
    w[24] =   0.510852245077594627738963204403E-07;
    w[25] =   0.917858042437852820850075742492E-09;
    w[26] =   0.810618629746304420399344796173E-11;
    w[27] =   0.287860708054870606219239791142E-13;
    w[28] =   0.281033360275090370876277491534E-16;
    w[29] =   0.290825470013122622941102747365E-20;
  }
  else if ( n == 31 )
  {
    x[0] = -6.99568012371854027532485214732E+00;
    x[1] = -6.27507870494286014270365678125E+00;
    x[2] = -5.67396144461858832963325587893E+00;
    x[3] = -5.13359557711238070458629689140E+00;
    x[4] = -4.63155950631285994206679976543E+00;
    x[5] = -4.15627175581814517248313523153E+00;
    x[6] = -3.70074340323146942244971645897E+00;
    x[7] = -3.26032073231354081046454015096E+00;
    x[8] = -2.83168045339020545570156401514E+00;
    x[9] = -2.41231770548042010517401845821E+00;
    x[10] = -2.00025854893563896579755625986E+00;
    x[11] = -1.59388586047213982613884194556E+00;
    x[12] = -1.19182699835004642608213586492E+00;
    x[13] = -0.792876976915308939685930329988E+00;
    x[14] = -0.395942736471423110946700416634E+00;
    x[15] = 0.0E+00;
    x[16] = 0.395942736471423110946700416634E+00;
    x[17] = 0.792876976915308939685930329988E+00;
    x[18] = 1.19182699835004642608213586492E+00;
    x[19] = 1.59388586047213982613884194556E+00;
    x[20] = 2.00025854893563896579755625986E+00;
    x[21] = 2.41231770548042010517401845821E+00;
    x[22] = 2.83168045339020545570156401514E+00;
    x[23] = 3.26032073231354081046454015096E+00;
    x[24] = 3.70074340323146942244971645897E+00;
    x[25] = 4.15627175581814517248313523153E+00;
    x[26] = 4.63155950631285994206679976543E+00;
    x[27] = 5.13359557711238070458629689140E+00;
    x[28] = 5.67396144461858832963325587893E+00;
    x[29] = 6.27507870494286014270365678125E+00;
    x[30] = 6.99568012371854027532485214732E+00;

    w[0] = 4.61896839446420502132944426974E-22;
    w[1] = 5.11060900792715640739422641166E-18;
    w[2] = 5.89955649875387299038431589378E-15;
    w[3] = 1.86037352145214652437380892603E-12;
    w[4] = 2.35249200320864163398597795323E-10;
    w[5] = 1.46119883449105307352780323055E-08;
    w[6] = 5.04371255893979974253745671633E-07;
    w[7] = 0.0000104986027576756063228123279208E+00;
    w[8] = 0.000139520903950470433823653754396E+00;
    w[9] = 0.00123368330730688826551750402319E+00;
    w[10] = 0.00748279991403519848345678003016E+00;
    w[11] = 0.0318472307313003327772087235339E+00;
    w[12] = 0.0967179481608704535580338478886E+00;
    w[13] = 0.212132788668764779877735637343E+00;
    w[14] = 0.338772657894107724675701919878E+00;
    w[15] = 0.395778556098609545141783810611E+00;
    w[16] = 0.338772657894107724675701919878E+00;
    w[17] = 0.212132788668764779877735637343E+00;
    w[18] = 0.0967179481608704535580338478886E+00;
    w[19] = 0.0318472307313003327772087235339E+00;
    w[20] = 0.00748279991403519848345678003016E+00;
    w[21] = 0.00123368330730688826551750402319E+00;
    w[22] = 0.000139520903950470433823653754396E+00;
    w[23] = 0.0000104986027576756063228123279208E+00;
    w[24] = 5.04371255893979974253745671633E-07;
    w[25] = 1.46119883449105307352780323055E-08;
    w[26] = 2.35249200320864163398597795323E-10;
    w[27] = 1.86037352145214652437380892603E-12;
    w[28] = 5.89955649875387299038431589378E-15;
    w[29] = 5.11060900792715640739422641166E-18;
    w[30] = 4.61896839446420502132944426974E-22;
  }
  else if ( n == 32 )
  {
    x[0] = -7.12581390983072757279520760342E+00;
    x[1] = -6.40949814926966041217376374153E+00;
    x[2] = -5.81222594951591383276596615366E+00;
    x[3] = -5.27555098651588012781906048140E+00;
    x[4] = -4.77716450350259639303579405689E+00;
    x[5] = -4.30554795335119844526348653193E+00;
    x[6] = -3.85375548547144464388787292109E+00;
    x[7] = -3.41716749281857073587392729564E+00;
    x[8] = -2.99249082500237420628549407606E+00;
    x[9] = -2.57724953773231745403092930114E+00;
    x[10] = -2.16949918360611217330570559502E+00;
    x[11] = -1.76765410946320160462767325853E+00;
    x[12] = -1.37037641095287183816170564864E+00;
    x[13] = -0.976500463589682838484704871982E+00;
    x[14] = -0.584978765435932448466957544011E+00;
    x[15] = -0.194840741569399326708741289532E+00;
    x[16] = 0.194840741569399326708741289532E+00;
    x[17] = 0.584978765435932448466957544011E+00;
    x[18] = 0.976500463589682838484704871982E+00;
    x[19] = 1.37037641095287183816170564864E+00;
    x[20] = 1.76765410946320160462767325853E+00;
    x[21] = 2.16949918360611217330570559502E+00;
    x[22] = 2.57724953773231745403092930114E+00;
    x[23] = 2.99249082500237420628549407606E+00;
    x[24] = 3.41716749281857073587392729564E+00;
    x[25] = 3.85375548547144464388787292109E+00;
    x[26] = 4.30554795335119844526348653193E+00;
    x[27] = 4.77716450350259639303579405689E+00;
    x[28] = 5.27555098651588012781906048140E+00;
    x[29] = 5.81222594951591383276596615366E+00;
    x[30] = 6.40949814926966041217376374153E+00;
    x[31] = 7.12581390983072757279520760342E+00;

    w[0] = 7.31067642738416239327427845506E-23;
    w[1] = 9.23173653651829223349442007207E-19;
    w[2] = 1.19734401709284866582868189951E-15;
    w[3] = 4.21501021132644757296944521183E-13;
    w[4] = 5.93329146339663861451156821558E-11;
    w[5] = 4.0988321647708966182350410138E-09;
    w[6] = 1.57416779254559402926869257929E-07;
    w[7] = 3.65058512956237605737032418746E-06;
    w[8] = 0.0000541658406181998255800193939267E+00;
    w[9] = 0.000536268365527972045970238101501E+00;
    w[10] = 0.00365489032665442807912565712241E+00;
    w[11] = 0.017553428831573430303437844611E+00;
    w[12] = 0.0604581309559126141865857607833E+00;
    w[13] = 0.151269734076642482575147114614E+00;
    w[14] = 0.277458142302529898137698918542E+00;
    w[15] = 0.375238352592802392866818388907E+00;
    w[16] = 0.375238352592802392866818388907E+00;
    w[17] = 0.277458142302529898137698918542E+00;
    w[18] = 0.151269734076642482575147114614E+00;
    w[19] = 0.0604581309559126141865857607833E+00;
    w[20] = 0.017553428831573430303437844611E+00;
    w[21] = 0.00365489032665442807912565712241E+00;
    w[22] = 0.000536268365527972045970238101501E+00;
    w[23] = 0.0000541658406181998255800193939267E+00;
    w[24] = 3.65058512956237605737032418746E-06;
    w[25] = 1.57416779254559402926869257929E-07;
    w[26] = 4.0988321647708966182350410138E-09;
    w[27] = 5.93329146339663861451156821558E-11;
    w[28] = 4.21501021132644757296944521183E-13;
    w[29] = 1.19734401709284866582868189951E-15;
    w[30] = 9.23173653651829223349442007207E-19;
    w[31] = 7.31067642738416239327427845506E-23;
  }
  else if ( n == 33 )
  {
    x[0] = -7.25385182201520064607977244465E+00;
    x[1] = -6.54165544573807726095826608811E+00;
    x[2] = -5.94807118208714447981366477584E+00;
    x[3] = -5.41492900261419253992709076454E+00;
    x[4] = -4.92002852059500829241139910265E+00;
    x[5] = -4.45191114883282719009473876206E+00;
    x[6] = -4.00367160995693141451378357174E+00;
    x[7] = -3.57072198023271828561890330658E+00;
    x[8] = -3.14979668170382538461281438786E+00;
    x[9] = -2.73844582435135490694887052899E+00;
    x[10] = -2.33475115152951517708536069773E+00;
    x[11] = -1.93715458182220661643452220908E+00;
    x[12] = -1.54434826124312180914304288754E+00;
    x[13] = -1.15520020412678961356412482063E+00;
    x[14] = -0.768701379758868598107224561306E+00;
    x[15] = -0.383926014508409083771145488401E+00;
    x[16] = 0.0E+00;
    x[17] = 0.383926014508409083771145488401E+00;
    x[18] = 0.768701379758868598107224561306E+00;
    x[19] = 1.15520020412678961356412482063E+00;
    x[20] = 1.54434826124312180914304288754E+00;
    x[21] = 1.93715458182220661643452220908E+00;
    x[22] = 2.33475115152951517708536069773E+00;
    x[23] = 2.73844582435135490694887052899E+00;
    x[24] = 3.14979668170382538461281438786E+00;
    x[25] = 3.57072198023271828561890330658E+00;
    x[26] = 4.00367160995693141451378357174E+00;
    x[27] = 4.45191114883282719009473876206E+00;
    x[28] = 4.92002852059500829241139910265E+00;
    x[29] = 5.41492900261419253992709076454E+00;
    x[30] = 5.94807118208714447981366477584E+00;
    x[31] = 6.54165544573807726095826608811E+00;
    x[32] = 7.25385182201520064607977244465E+00;

    w[0] = 1.15331621854588454082208890757E-23;
    w[1] = 1.65709474153369051048226040291E-19;
    w[2] = 2.40778567955799442824587707068E-16;
    w[3] = 9.43481415901497503451931586527E-14;
    w[4] = 1.47398093709248867676655543441E-11;
    w[5] = 1.12892224710833129085848357165E-09;
    w[6] = 4.8077456763231909801575826594E-08;
    w[7] = 1.23769336720121013593677278301E-06;
    w[8] = 0.000020423684051423773240532989618E+00;
    w[9] = 0.000225442770596327415479556999963E+00;
    w[10] = 0.00171845463776092445684897958755E+00;
    w[11] = 0.00926568997068523330372581072023E+00;
    w[12] = 0.0359879823185769744486866448437E+00;
    w[13] = 0.102069079995541500792505520921E+00;
    w[14] = 0.213493931150291836488258605194E+00;
    w[15] = 0.331552000750741282288352789764E+00;
    w[16] = 0.383785266519863801349608543622E+00;
    w[17] = 0.331552000750741282288352789764E+00;
    w[18] = 0.213493931150291836488258605194E+00;
    w[19] = 0.102069079995541500792505520921E+00;
    w[20] = 0.0359879823185769744486866448437E+00;
    w[21] = 0.00926568997068523330372581072023E+00;
    w[22] = 0.00171845463776092445684897958755E+00;
    w[23] = 0.000225442770596327415479556999963E+00;
    w[24] = 0.000020423684051423773240532989618E+00;
    w[25] = 1.23769336720121013593677278301E-06;
    w[26] = 4.8077456763231909801575826594E-08;
    w[27] = 1.12892224710833129085848357165E-09;
    w[28] = 1.47398093709248867676655543441E-11;
    w[29] = 9.43481415901497503451931586527E-14;
    w[30] = 2.40778567955799442824587707068E-16;
    w[31] = 1.65709474153369051048226040291E-19;
    w[32] = 1.15331621854588454082208890757E-23;
  }
  else if ( n == 63 )
  {
    x[0] = -10.4354998778541680534681154273E+00;
    x[1] = -9.80287599129749636352239352865E+00;
    x[2] = -9.27920195430503913194047455065E+00;
    x[3] = -8.81185814372845464425266282756E+00;
    x[4] = -8.38076834518632193430106510438E+00;
    x[5] = -7.97559508014203731815418062985E+00;
    x[6] = -7.59013951986410667624797831945E+00;
    x[7] = -7.22031670788896784611613242225E+00;
    x[8] = -6.86325443317953685273532858761E+00;
    x[9] = -6.51683481068211606052733958540E+00;
    x[10] = -6.17943799227059698624184617873E+00;
    x[11] = -5.84978840008106734625265829615E+00;
    x[12] = -5.52685725264030314250475751228E+00;
    x[13] = -5.20979798304083548615751364163E+00;
    x[14] = -4.89790186449757423507450992149E+00;
    x[15] = -4.59056657444351902292712945691E+00;
    x[16] = -4.28727333528244040317276161995E+00;
    x[17] = -3.98756991041971574852270520681E+00;
    x[18] = -3.69105770009634651173228105598E+00;
    x[19] = -3.39738177133039118527559418063E+00;
    x[20] = -3.10622302792825663291386167460E+00;
    x[21] = -2.81729196728379777507471356574E+00;
    x[22] = -2.53032363047120109268552217185E+00;
    x[23] = -2.24507346048120662989959181793E+00;
    x[24] = -1.96131385830814852939220084113E+00;
    x[25] = -1.67883127917201375208028006226E+00;
    x[26] = -1.39742374860496251075707520637E+00;
    x[27] = -1.11689870509964626905109702778E+00;
    x[28] = -0.837071095589476159777377954613E+00;
    x[29] = -0.557761664279082216687636652538E+00;
    x[30] = -0.278795385671152239866876286272E+00;
    x[31] = 0.0E+00;
    x[32] = 0.278795385671152239866876286272E+00;
    x[33] = 0.557761664279082216687636652538E+00;
    x[34] = 0.837071095589476159777377954613E+00;
    x[35] = 1.11689870509964626905109702778E+00;
    x[36] = 1.39742374860496251075707520637E+00;
    x[37] = 1.67883127917201375208028006226E+00;
    x[38] = 1.96131385830814852939220084113E+00;
    x[39] = 2.24507346048120662989959181793E+00;
    x[40] = 2.53032363047120109268552217185E+00;
    x[41] = 2.81729196728379777507471356574E+00;
    x[42] = 3.10622302792825663291386167460E+00;
    x[43] = 3.39738177133039118527559418063E+00;
    x[44] = 3.69105770009634651173228105598E+00;
    x[45] = 3.98756991041971574852270520681E+00;
    x[46] = 4.28727333528244040317276161995E+00;
    x[47] = 4.59056657444351902292712945691E+00;
    x[48] = 4.89790186449757423507450992149E+00;
    x[49] = 5.20979798304083548615751364163E+00;
    x[50] = 5.52685725264030314250475751228E+00;
    x[51] = 5.84978840008106734625265829615E+00;
    x[52] = 6.17943799227059698624184617873E+00;
    x[53] = 6.51683481068211606052733958540E+00;
    x[54] = 6.86325443317953685273532858761E+00;
    x[55] = 7.22031670788896784611613242225E+00;
    x[56] = 7.59013951986410667624797831945E+00;
    x[57] = 7.97559508014203731815418062985E+00;
    x[58] = 8.38076834518632193430106510438E+00;
    x[59] = 8.81185814372845464425266282756E+00;
    x[60] = 9.27920195430503913194047455065E+00;
    x[61] = 9.80287599129749636352239352865E+00;
    x[62] = 10.4354998778541680534681154273E+00;

    w[0] = 3.70992064349030055823376157823E-48;
    w[1] = 1.04007786152246672212559599908E-42;
    w[2] = 1.97968047083199197900260998813E-38;
    w[3] = 8.46874781919035663281042885251E-35;
    w[4] = 1.30713059308206243904769877879E-31;
    w[5] = 9.34378371756582396450246862195E-29;
    w[6] = 3.60274266352851638202340658522E-26;
    w[7] = 8.29638631162099766157527065317E-24;
    w[8] = 1.22666299091434557721622529775E-21;
    w[9] = 1.22884356288353036990240371039E-19;
    w[10] = 8.69255369584585252225619256428E-18;
    w[11] = 4.48570586893158184069444097978E-16;
    w[12] = 1.73358179557891044383064226749E-14;
    w[13] = 5.1265062385197846998384009333E-13;
    w[14] = 1.18089218445696923817995132237E-11;
    w[15] = 2.15086982978749617679069862879E-10;
    w[16] = 3.13719295353830786449435629291E-09;
    w[17] = 3.70416259848969809883356560995E-08;
    w[18] = 3.57347329499908777461505032558E-07;
    w[19] = 2.83931144984692884712301165567E-06;
    w[20] = 0.0000187091130037887216027832755405E+00;
    w[21] = 0.000102848808006856425543062213642E+00;
    w[22] = 0.000474117026103206754395975199216E+00;
    w[23] = 0.0018409222622442103760124297917E+00;
    w[24] = 0.00604360445513757113209247151533E+00;
    w[25] = 0.0168292991996521044559098701555E+00;
    w[26] = 0.0398582640278170328649908688578E+00;
    w[27] = 0.0804670879942008323850873860195E+00;
    w[28] = 0.138719508176584635072239096351E+00;
    w[29] = 0.204486953468973988225911656103E+00;
    w[30] = 0.25799889943138332612723393346E+00;
    w[31] = 0.278766948849251654365527505911E+00;
    w[32] = 0.25799889943138332612723393346E+00;
    w[33] = 0.204486953468973988225911656103E+00;
    w[34] = 0.138719508176584635072239096351E+00;
    w[35] = 0.0804670879942008323850873860195E+00;
    w[36] = 0.0398582640278170328649908688578E+00;
    w[37] = 0.0168292991996521044559098701555E+00;
    w[38] = 0.00604360445513757113209247151533E+00;
    w[39] = 0.0018409222622442103760124297917E+00;
    w[40] = 0.000474117026103206754395975199216E+00;
    w[41] = 0.000102848808006856425543062213642E+00;
    w[42] = 0.0000187091130037887216027832755405E+00;
    w[43] = 2.83931144984692884712301165567E-06;
    w[44] = 3.57347329499908777461505032558E-07;
    w[45] = 3.70416259848969809883356560995E-08;
    w[46] = 3.13719295353830786449435629291E-09;
    w[47] = 2.15086982978749617679069862879E-10;
    w[48] = 1.18089218445696923817995132237E-11;
    w[49] = 5.1265062385197846998384009333E-13;
    w[50] = 1.73358179557891044383064226749E-14;
    w[51] = 4.48570586893158184069444097978E-16;
    w[52] = 8.69255369584585252225619256428E-18;
    w[53] = 1.22884356288353036990240371039E-19;
    w[54] = 1.22666299091434557721622529775E-21;
    w[55] = 8.29638631162099766157527065317E-24;
    w[56] = 3.60274266352851638202340658522E-26;
    w[57] = 9.34378371756582396450246862195E-29;
    w[58] = 1.30713059308206243904769877879E-31;
    w[59] = 8.46874781919035663281042885251E-35;
    w[60] = 1.97968047083199197900260998813E-38;
    w[61] = 1.04007786152246672212559599908E-42;
    w[62] = 3.70992064349030055823376157823E-48;
  }
  else if ( n == 64 )
  {
    x[0] = -10.5261231679605458833268262838E+00;
    x[1] = -9.89528758682953902120446147716E+00;
    x[2] = -9.37315954964672116254565243972E+00;
    x[3] = -8.90724909996476975729597288564E+00;
    x[4] = -8.47752908337986309056416634482E+00;
    x[5] = -8.07368728501022522585879114076E+00;
    x[6] = -7.68954016404049682844780422987E+00;
    x[7] = -7.32101303278094920118956936372E+00;
    x[8] = -6.96524112055110752924264219349E+00;
    x[9] = -6.62011226263602737903666010894E+00;
    x[10] = -6.28401122877482823541809319507E+00;
    x[11] = -5.95566632679948604534456718098E+00;
    x[12] = -5.63405216434997214724992048331E+00;
    x[13] = -5.31832522463327085732364951520E+00;
    x[14] = -5.00777960219876819644370262718E+00;
    x[15] = -4.70181564740749981609753801581E+00;
    x[16] = -4.39991716822813764776793253544E+00;
    x[17] = -4.10163447456665671497098123846E+00;
    x[18] = -3.80657151394536046116597200046E+00;
    x[19] = -3.51437593574090621153995058647E+00;
    x[20] = -3.22473129199203572584817111019E+00;
    x[21] = -2.93735082300462180968533902619E+00;
    x[22] = -2.65197243543063501100545778600E+00;
    x[23] = -2.36835458863240140411151126534E+00;
    x[24] = -2.08627287988176202083256330236E+00;
    x[25] = -1.80551717146554491890377357419E+00;
    x[26] = -1.52588914020986366294897013315E+00;
    x[27] = -1.24720015694311794069356453069E+00;
    x[28] = -0.969269423071178016743541489019E+00;
    x[29] = -0.691922305810044577268219287596E+00;
    x[30] = -0.414988824121078684576929129200E+00;
    x[31] = -0.138302244987009724115049767967E+00;
    x[32] = 0.138302244987009724115049767967E+00;
    x[33] = 0.414988824121078684576929129200E+00;
    x[34] = 0.691922305810044577268219287596E+00;
    x[35] = 0.969269423071178016743541489019E+00;
    x[36] = 1.24720015694311794069356453069E+00;
    x[37] = 1.52588914020986366294897013315E+00;
    x[38] = 1.80551717146554491890377357419E+00;
    x[39] = 2.08627287988176202083256330236E+00;
    x[40] = 2.36835458863240140411151126534E+00;
    x[41] = 2.65197243543063501100545778600E+00;
    x[42] = 2.93735082300462180968533902619E+00;
    x[43] = 3.22473129199203572584817111019E+00;
    x[44] = 3.51437593574090621153995058647E+00;
    x[45] = 3.80657151394536046116597200046E+00;
    x[46] = 4.10163447456665671497098123846E+00;
    x[47] = 4.39991716822813764776793253544E+00;
    x[48] = 4.70181564740749981609753801581E+00;
    x[49] = 5.00777960219876819644370262718E+00;
    x[50] = 5.31832522463327085732364951520E+00;
    x[51] = 5.63405216434997214724992048331E+00;
    x[52] = 5.95566632679948604534456718098E+00;
    x[53] = 6.28401122877482823541809319507E+00;
    x[54] = 6.62011226263602737903666010894E+00;
    x[55] = 6.96524112055110752924264219349E+00;
    x[56] = 7.32101303278094920118956936372E+00;
    x[57] = 7.68954016404049682844780422987E+00;
    x[58] = 8.07368728501022522585879114076E+00;
    x[59] = 8.47752908337986309056416634482E+00;
    x[60] = 8.90724909996476975729597288564E+00;
    x[61] = 9.37315954964672116254565243972E+00;
    x[62] = 9.89528758682953902120446147716E+00;
    x[63] = 10.5261231679605458833268262838E+00;

    w[0] = 5.53570653585694282057546330099E-49;
    w[1] = 1.67974799010815921866628833063E-43;
    w[2] = 3.42113801125574050432722182815E-39;
    w[3] = 1.55739062462976380230933538026E-35;
    w[4] = 2.54966089911299925660476658044E-32;
    w[5] = 1.92910359546496685030196877907E-29;
    w[6] = 7.86179778892591036909999149628E-27;
    w[7] = 1.91170688330064282995845696553E-24;
    w[8] = 2.98286278427985115447870070202E-22;
    w[9] = 3.15225456650378141612134668341E-20;
    w[10] = 2.35188471067581911695767591556E-18;
    w[11] = 1.28009339132243804163956329526E-16;
    w[12] = 5.21862372659084752295780851305E-15;
    w[13] = 1.62834073070972036208430708124E-13;
    w[14] = 3.95917776694772392723644586425E-12;
    w[15] = 7.61521725014545135331529567532E-11;
    w[16] = 1.17361674232154934354250646708E-09;
    w[17] = 1.4651253164761093549266220038E-08;
    w[18] = 1.49553293672724706110246169293E-07;
    w[19] = 1.25834025103118457615784218002E-06;
    w[20] = 8.7884992308503591814440474067E-06;
    w[21] = 0.0000512592913578627466082191141274E+00;
    w[22] = 0.000250983698513062486082362017982E+00;
    w[23] = 0.00103632909950757766345674174628E+00;
    w[24] = 0.00362258697853445876066812537162E+00;
    w[25] = 0.0107560405098791370494651727867E+00;
    w[26] = 0.0272031289536889184538348212615E+00;
    w[27] = 0.0587399819640994345496889462518E+00;
    w[28] = 0.108498349306186840633025845506E+00;
    w[29] = 0.171685842349083702000727970124E+00;
    w[30] = 0.232994786062678046650566029333E+00;
    w[31] = 0.271377424941303977945606508418E+00;
    w[32] = 0.271377424941303977945606508418E+00;
    w[33] = 0.232994786062678046650566029333E+00;
    w[34] = 0.171685842349083702000727970124E+00;
    w[35] = 0.108498349306186840633025845506E+00;
    w[36] = 0.0587399819640994345496889462518E+00;
    w[37] = 0.0272031289536889184538348212615E+00;
    w[38] = 0.0107560405098791370494651727867E+00;
    w[39] = 0.00362258697853445876066812537162E+00;
    w[40] = 0.00103632909950757766345674174628E+00;
    w[41] = 0.000250983698513062486082362017982E+00;
    w[42] = 0.0000512592913578627466082191141274E+00;
    w[43] = 8.7884992308503591814440474067E-06;
    w[44] = 1.25834025103118457615784218002E-06;
    w[45] = 1.49553293672724706110246169293E-07;
    w[46] = 1.4651253164761093549266220038E-08;
    w[47] = 1.17361674232154934354250646708E-09;
    w[48] = 7.61521725014545135331529567532E-11;
    w[49] = 3.95917776694772392723644586425E-12;
    w[50] = 1.62834073070972036208430708124E-13;
    w[51] = 5.21862372659084752295780851305E-15;
    w[52] = 1.28009339132243804163956329526E-16;
    w[53] = 2.35188471067581911695767591556E-18;
    w[54] = 3.15225456650378141612134668341E-20;
    w[55] = 2.98286278427985115447870070202E-22;
    w[56] = 1.91170688330064282995845696553E-24;
    w[57] = 7.86179778892591036909999149628E-27;
    w[58] = 1.92910359546496685030196877907E-29;
    w[59] = 2.54966089911299925660476658044E-32;
    w[60] = 1.55739062462976380230933538026E-35;
    w[61] = 3.42113801125574050432722182815E-39;
    w[62] = 1.67974799010815921866628833063E-43;
    w[63] = 5.53570653585694282057546330099E-49;
  }
  else if ( n == 65 )
  {
    x[0] = -10.6160229818782811890575602362E+00;
    x[1] = -9.98694169167668475289516351866E+00;
    x[2] = -9.46632932015538456209230219583E+00;
    x[3] = -9.00182332295913301957373193078E+00;
    x[4] = -8.57344474441790920512961590649E+00;
    x[5] = -8.17090617805258532129873946832E+00;
    x[6] = -7.78803908298957078257446963177E+00;
    x[7] = -7.42077883436632423648540822493E+00;
    x[8] = -7.06626794030689283695626648500E+00;
    x[9] = -6.72239982016573443713269171234E+00;
    x[10] = -6.38756373978709108766264932629E+00;
    x[11] = -6.06049177883150521431844958114E+00;
    x[12] = -5.74016182369022552144785845004E+00;
    x[13] = -5.42573329769734933377879520007E+00;
    x[14] = -5.11650300472141247403405938903E+00;
    x[15] = -4.81187385202746476469375701958E+00;
    x[16] = -4.51133211136821338520053781440E+00;
    x[17] = -4.21443050997195460766227433230E+00;
    x[18] = -3.92077540444472380734866143333E+00;
    x[19] = -3.63001687763289533380743256995E+00;
    x[20] = -3.34184096844683005348576537026E+00;
    x[21] = -3.05596348432867100471992904351E+00;
    x[22] = -2.77212500515709167382401243242E+00;
    x[23] = -2.49008679530393550665251461372E+00;
    x[24] = -2.20962741516918436363972079745E+00;
    x[25] = -1.93053987597722550177796442890E+00;
    x[26] = -1.65262921904032541019598694267E+00;
    x[27] = -1.37571042772366833613088581439E+00;
    x[28] = -1.09960660005694495033699221150E+00;
    x[29] = -0.824147324402412861055989047706E+00;
    x[30] = -0.549167211221599184571872835161E+00;
    x[31] = -0.274504541753944755855051087074E+00;
    x[32] = 0.0E+00;
    x[33] = 0.274504541753944755855051087074E+00;
    x[34] = 0.549167211221599184571872835161E+00;
    x[35] = 0.824147324402412861055989047706E+00;
    x[36] = 1.09960660005694495033699221150E+00;
    x[37] = 1.37571042772366833613088581439E+00;
    x[38] = 1.65262921904032541019598694267E+00;
    x[39] = 1.93053987597722550177796442890E+00;
    x[40] = 2.20962741516918436363972079745E+00;
    x[41] = 2.49008679530393550665251461372E+00;
    x[42] = 2.77212500515709167382401243242E+00;
    x[43] = 3.05596348432867100471992904351E+00;
    x[44] = 3.34184096844683005348576537026E+00;
    x[45] = 3.63001687763289533380743256995E+00;
    x[46] = 3.92077540444472380734866143333E+00;
    x[47] = 4.21443050997195460766227433230E+00;
    x[48] = 4.51133211136821338520053781440E+00;
    x[49] = 4.81187385202746476469375701958E+00;
    x[50] = 5.11650300472141247403405938903E+00;
    x[51] = 5.42573329769734933377879520007E+00;
    x[52] = 5.74016182369022552144785845004E+00;
    x[53] = 6.06049177883150521431844958114E+00;
    x[54] = 6.38756373978709108766264932629E+00;
    x[55] = 6.72239982016573443713269171234E+00;
    x[56] = 7.06626794030689283695626648500E+00;
    x[57] = 7.42077883436632423648540822493E+00;
    x[58] = 7.78803908298957078257446963177E+00;
    x[59] = 8.17090617805258532129873946832E+00;
    x[60] = 8.57344474441790920512961590649E+00;
    x[61] = 9.00182332295913301957373193078E+00;
    x[62] = 9.46632932015538456209230219583E+00;
    x[63] = 9.98694169167668475289516351866E+00;
    x[64] = 10.6160229818782811890575602362E+00;

    w[0] = 8.25161081325244640518686536873E-50;
    w[1] = 2.70767584528327632245086261566E-44;
    w[2] = 5.89628446597893219238447711362E-40;
    w[3] = 2.8541849032786262808377028501E-36;
    w[4] = 4.95258625502059879210418105309E-33;
    w[5] = 3.96328698707468682361835959189E-30;
    w[6] = 1.70591158107580273148997822331E-27;
    w[7] = 4.37697419487184691809226004173E-25;
    w[8] = 7.20161078913500757836854034749E-23;
    w[9] = 8.0222187354240312838311535001E-21;
    w[10] = 6.30789104558609987896303941119E-19;
    w[11] = 3.61819961904286485492939434525E-17;
    w[12] = 1.55466357223809604941702812296E-15;
    w[13] = 5.11391748171652449009988302839E-14;
    w[14] = 1.31125161063902569430172028735E-12;
    w[15] = 2.66086534779295548413319751434E-11;
    w[16] = 4.32865615344850974821379264835E-10;
    w[17] = 5.70758293277877491250362877931E-09;
    w[18] = 6.15779622145053848599380659292E-08;
    w[19] = 5.48045603501799498244047819842E-07;
    w[20] = 4.05224939102373376093012342174E-06;
    w[21] = 0.0000250453428904958321201946621231E+00;
    w[22] = 0.000130082916298451204382435919638E+00;
    w[23] = 0.000570398967523771524725931177791E+00;
    w[24] = 0.00211998163203684165580510045255E+00;
    w[25] = 0.00670140453800573713948633700424E+00;
    w[26] = 0.018069433112703589006399924887E+00;
    w[27] = 0.0416611087624784398909512089873E+00;
    w[28] = 0.0823001633697352251543326980867E+00;
    w[29] = 0.139526139482843953007755621004E+00;
    w[30] = 0.203250574154441897747728738409E+00;
    w[31] = 0.254628811852790103887643365928E+00;
    w[32] = 0.274478226559263167375288621205E+00;
    w[33] = 0.254628811852790103887643365928E+00;
    w[34] = 0.203250574154441897747728738409E+00;
    w[35] = 0.139526139482843953007755621004E+00;
    w[36] = 0.0823001633697352251543326980867E+00;
    w[37] = 0.0416611087624784398909512089873E+00;
    w[38] = 0.018069433112703589006399924887E+00;
    w[39] = 0.00670140453800573713948633700424E+00;
    w[40] = 0.00211998163203684165580510045255E+00;
    w[41] = 0.000570398967523771524725931177791E+00;
    w[42] = 0.000130082916298451204382435919638E+00;
    w[43] = 0.0000250453428904958321201946621231E+00;
    w[44] = 4.05224939102373376093012342174E-06;
    w[45] = 5.48045603501799498244047819842E-07;
    w[46] = 6.15779622145053848599380659292E-08;
    w[47] = 5.70758293277877491250362877931E-09;
    w[48] = 4.32865615344850974821379264835E-10;
    w[49] = 2.66086534779295548413319751434E-11;
    w[50] = 1.31125161063902569430172028735E-12;
    w[51] = 5.11391748171652449009988302839E-14;
    w[52] = 1.55466357223809604941702812296E-15;
    w[53] = 3.61819961904286485492939434525E-17;
    w[54] = 6.30789104558609987896303941119E-19;
    w[55] = 8.0222187354240312838311535001E-21;
    w[56] = 7.20161078913500757836854034749E-23;
    w[57] = 4.37697419487184691809226004173E-25;
    w[58] = 1.70591158107580273148997822331E-27;
    w[59] = 3.96328698707468682361835959189E-30;
    w[60] = 4.95258625502059879210418105309E-33;
    w[61] = 2.8541849032786262808377028501E-36;
    w[62] = 5.89628446597893219238447711362E-40;
    w[63] = 2.70767584528327632245086261566E-44;
    w[64] = 8.25161081325244640518686536873E-50;
  }
  else if ( n == 127 )
  {
    x[0] = -15.2283381481673509782469544335E+00;
    x[1] = -14.6695951588339726327463541129E+00;
    x[2] = -14.2090859952848707551682442509E+00;
    x[3] = -13.7997222902116766346452467467E+00;
    x[4] = -13.4235185900709500624382583219E+00;
    x[5] = -13.0712086604746019015839954396E+00;
    x[6] = -12.7372356524156863381380039241E+00;
    x[7] = -12.4179393788697158054458796241E+00;
    x[8] = -12.1107490209477476001321235081E+00;
    x[9] = -11.8137721982677271951345841362E+00;
    x[10] = -11.5255651125726965991678885886E+00;
    x[11] = -11.2449945837855434451943841943E+00;
    x[12] = -10.9711505698402474234230402639E+00;
    x[13] = -10.7032882010274813476709407447E+00;
    x[14] = -10.4407879577727728677425917980E+00;
    x[15] = -10.1831274734503438886241264504E+00;
    x[16] = -9.92986104951142507368470042737E+00;
    x[17] = -9.68060444124747280381507127327E+00;
    x[18] = -9.43502333898816501350195985063E+00;
    x[19] = -9.19282449884603057157741950525E+00;
    x[20] = -8.95374881085654043238078901700E+00;
    x[21] = -8.71756580870763073638339995485E+00;
    x[22] = -8.48406926898324733260971803400E+00;
    x[23] = -8.25307364544571565796941242439E+00;
    x[24] = -8.02441115147033755785947397968E+00;
    x[25] = -7.79792935138701054208291204556E+00;
    x[26] = -7.57348915560834540228349607633E+00;
    x[27] = -7.35096313922690527019612580437E+00;
    x[28] = -7.13023412203507106680640257134E+00;
    x[29] = -6.91119396154657131974656331094E+00;
    x[30] = -6.69374252087582941900744173817E+00;
    x[31] = -6.47778678116453654481449038215E+00;
    x[32] = -6.26324007427373543456097238571E+00;
    x[33] = -6.05002141614198456944654744824E+00;
    x[34] = -5.83805492487741873866016908078E+00;
    x[35] = -5.62726931054648166594234557949E+00;
    x[36] = -5.41759742592432407228484258729E+00;
    x[37] = -5.20897586931539835875702583722E+00;
    x[38] = -5.00134463203863600385208091074E+00;
    x[39] = -4.79464678437649250097485099309E+00;
    x[40] = -4.58882819476983729516064850312E+00;
    x[41] = -4.38383727784647362942537444075E+00;
    x[42] = -4.17962476753520313494211898924E+00;
    x[43] = -3.97614351206733559160358141959E+00;
    x[44] = -3.77334828812505267210046784001E+00;
    x[45] = -3.57119563177821804471997564852E+00;
    x[46] = -3.36964368417173978966436292400E+00;
    x[47] = -3.16865205019536301918577982615E+00;
    x[48] = -2.96818166859559102677616495215E+00;
    x[49] = -2.76819469218240588012265459589E+00;
    x[50] = -2.56865437694735017231440130224E+00;
    x[51] = -2.36952497904904010800124746457E+00;
    x[52] = -2.17077165874115068794984980837E+00;
    x[53] = -1.97236039041950200793247432276E+00;
    x[54] = -1.77425787805167915846764421037E+00;
    x[55] = -1.57643147532678013155195976219E+00;
    x[56] = -1.37884910992617780914415570537E+00;
    x[57] = -1.18147921137006858486785835984E+00;
    x[58] = -0.984290641940272777265689842138E+00;
    x[59] = -0.787252630218250341515968318790E+00;
    x[60] = -0.590334706809421021422304393461E+00;
    x[61] = -0.393506641851301365680378262002E+00;
    x[62] = -0.196738383924232519642722397371E+00;
    x[63] = 0.0E+00;
    x[64] = 0.196738383924232519642722397371E+00;
    x[65] = 0.393506641851301365680378262002E+00;
    x[66] = 0.590334706809421021422304393461E+00;
    x[67] = 0.787252630218250341515968318790E+00;
    x[68] = 0.984290641940272777265689842138E+00;
    x[69] = 1.18147921137006858486785835984E+00;
    x[70] = 1.37884910992617780914415570537E+00;
    x[71] = 1.57643147532678013155195976219E+00;
    x[72] = 1.77425787805167915846764421037E+00;
    x[73] = 1.97236039041950200793247432276E+00;
    x[74] = 2.17077165874115068794984980837E+00;
    x[75] = 2.36952497904904010800124746457E+00;
    x[76] = 2.56865437694735017231440130224E+00;
    x[77] = 2.76819469218240588012265459589E+00;
    x[78] = 2.96818166859559102677616495215E+00;
    x[79] = 3.16865205019536301918577982615E+00;
    x[80] = 3.36964368417173978966436292400E+00;
    x[81] = 3.57119563177821804471997564852E+00;
    x[82] = 3.77334828812505267210046784001E+00;
    x[83] = 3.97614351206733559160358141959E+00;
    x[84] = 4.17962476753520313494211898924E+00;
    x[85] = 4.38383727784647362942537444075E+00;
    x[86] = 4.58882819476983729516064850312E+00;
    x[87] = 4.79464678437649250097485099309E+00;
    x[88] = 5.00134463203863600385208091074E+00;
    x[89] = 5.20897586931539835875702583722E+00;
    x[90] = 5.41759742592432407228484258729E+00;
    x[91] = 5.62726931054648166594234557949E+00;
    x[92] = 5.83805492487741873866016908078E+00;
    x[93] = 6.05002141614198456944654744824E+00;
    x[94] = 6.26324007427373543456097238571E+00;
    x[95] = 6.47778678116453654481449038215E+00;
    x[96] = 6.69374252087582941900744173817E+00;
    x[97] = 6.91119396154657131974656331094E+00;
    x[98] = 7.13023412203507106680640257134E+00;
    x[99] = 7.35096313922690527019612580437E+00;
    x[100] = 7.57348915560834540228349607633E+00;
    x[101] = 7.79792935138701054208291204556E+00;
    x[102] = 8.02441115147033755785947397968E+00;
    x[103] = 8.25307364544571565796941242439E+00;
    x[104] = 8.48406926898324733260971803400E+00;
    x[105] = 8.71756580870763073638339995485E+00;
    x[106] = 8.95374881085654043238078901700E+00;
    x[107] = 9.19282449884603057157741950525E+00;
    x[108] = 9.43502333898816501350195985063E+00;
    x[109] = 9.68060444124747280381507127327E+00;
    x[110] = 9.92986104951142507368470042737E+00;
    x[111] = 10.1831274734503438886241264504E+00;
    x[112] = 10.4407879577727728677425917980E+00;
    x[113] = 10.7032882010274813476709407447E+00;
    x[114] = 10.9711505698402474234230402639E+00;
    x[115] = 11.2449945837855434451943841943E+00;
    x[116] = 11.5255651125726965991678885886E+00;
    x[117] = 11.8137721982677271951345841362E+00;
    x[118] = 12.1107490209477476001321235081E+00;
    x[119] = 12.4179393788697158054458796241E+00;
    x[120] = 12.7372356524156863381380039241E+00;
    x[121] = 13.0712086604746019015839954396E+00;
    x[122] = 13.4235185900709500624382583219E+00;
    x[123] = 13.7997222902116766346452467467E+00;
    x[124] = 14.2090859952848707551682442509E+00;
    x[125] = 14.6695951588339726327463541129E+00;
    x[126] = 15.2283381481673509782469544335E+00;

    w[0] = 1.25044975770895101066558695394E-101;
    w[1] = 1.72727980594728851329952877284E-94;
    w[2] = 8.93216815722645216635320162557E-89;
    w[3] = 7.7306185241134158744827181222E-84;
    w[4] = 2.01439576527109443920782513994E-79;
    w[5] = 2.15037147336771602203551878273E-75;
    w[6] = 1.13419242086298913875376620343E-71;
    w[7] = 3.34891390118992716444169809114E-68;
    w[8] = 6.04865489642049179016214753843E-65;
    w[9] = 7.13750929465743002965122123123E-62;
    w[10] = 5.78845633750656959788340019085E-59;
    w[11] = 3.3581166223962736386929935773E-56;
    w[12] = 1.4394641949298720336141068619E-53;
    w[13] = 4.68218083833618292793410025836E-51;
    w[14] = 1.18170544407210392716367665268E-48;
    w[15] = 2.35816591560823143778744566357E-46;
    w[16] = 3.78144279409152203964384313149E-44;
    w[17] = 4.9411031115925407477456893331E-42;
    w[18] = 5.32553037755907921458489847863E-40;
    w[19] = 4.78543906802804099967221020647E-38;
    w[20] = 3.61918834460649868835433546523E-36;
    w[21] = 2.3232083386415854084664074623E-34;
    w[22] = 1.27533314110484056196532640642E-32;
    w[23] = 6.02777538509463291699314327193E-31;
    w[24] = 2.4679773241854004762148469348E-29;
    w[25] = 8.8019567691972403392314252914E-28;
    w[26] = 2.74824892121260880467531987939E-26;
    w[27] = 7.54682189033203465872349657723E-25;
    w[28] = 1.83031346363374264415878982576E-23;
    w[29] = 3.93559908609832906838466602268E-22;
    w[30] = 7.52931616388155067444192947319E-21;
    w[31] = 1.28579977867628696999762170542E-19;
    w[32] = 1.96593268885070384943390296306E-18;
    w[33] = 2.69865119072980851232572568063E-17;
    w[34] = 3.33444143033026256341061235315E-16;
    w[35] = 3.71733031252663248624409938613E-15;
    w[36] = 3.74739544729563577089986076081E-14;
    w[37] = 3.42300944935037851188976963928E-13;
    w[38] = 2.83853037250817094975750489262E-12;
    w[39] = 2.14069202905212884993201956606E-11;
    w[40] = 1.47063312734774830028408333227E-10;
    w[41] = 9.21739409677215086782446989876E-10;
    w[42] = 5.27816639371369729333040255118E-09;
    w[43] = 2.76504970450371674155194812923E-08;
    w[44] = 1.32678558425807549298485884004E-07;
    w[45] = 5.83809442762947462901022315301E-07;
    w[46] = 2.35815617248490159838145978859E-06;
    w[47] = 8.75244680345528247507614056972E-06;
    w[48] = 0.0000298767905360019901790649251988E+00;
    w[49] = 0.0000938744357203646866361259710004E+00;
    w[50] = 0.000271707626280157286781639661883E+00;
    w[51] = 0.000724939297427239633212185817821E+00;
    w[52] = 0.0017841208326818955520088211458E+00;
    w[53] = 0.00405248551861722466559241860023E+00;
    w[54] = 0.00850002630418086349941683729112E+00;
    w[55] = 0.0164711422416609467530350356258E+00;
    w[56] = 0.0294992962483054353948393364098E+00;
    w[57] = 0.0488473871144520262535428484316E+00;
    w[58] = 0.074807989768816537216026182806E+00;
    w[59] = 0.10598520508123912472195529192E+00;
    w[60] = 0.138939453090947794093360848265E+00;
    w[61] = 0.168562360742603870987330592834E+00;
    w[62] = 0.189278495801793364889704841035E+00;
    w[63] = 0.196733406888845140995323677102E+00;
    w[64] = 0.189278495801793364889704841035E+00;
    w[65] = 0.168562360742603870987330592834E+00;
    w[66] = 0.138939453090947794093360848265E+00;
    w[67] = 0.10598520508123912472195529192E+00;
    w[68] = 0.074807989768816537216026182806E+00;
    w[69] = 0.0488473871144520262535428484316E+00;
    w[70] = 0.0294992962483054353948393364098E+00;
    w[71] = 0.0164711422416609467530350356258E+00;
    w[72] = 0.00850002630418086349941683729112E+00;
    w[73] = 0.00405248551861722466559241860023E+00;
    w[74] = 0.0017841208326818955520088211458E+00;
    w[75] = 0.000724939297427239633212185817821E+00;
    w[76] = 0.000271707626280157286781639661883E+00;
    w[77] = 0.0000938744357203646866361259710004E+00;
    w[78] = 0.0000298767905360019901790649251988E+00;
    w[79] = 8.75244680345528247507614056972E-06;
    w[80] = 2.35815617248490159838145978859E-06;
    w[81] = 5.83809442762947462901022315301E-07;
    w[82] = 1.32678558425807549298485884004E-07;
    w[83] = 2.76504970450371674155194812923E-08;
    w[84] = 5.27816639371369729333040255118E-09;
    w[85] = 9.21739409677215086782446989876E-10;
    w[86] = 1.47063312734774830028408333227E-10;
    w[87] = 2.14069202905212884993201956606E-11;
    w[88] = 2.83853037250817094975750489262E-12;
    w[89] = 3.42300944935037851188976963928E-13;
    w[90] = 3.74739544729563577089986076081E-14;
    w[91] = 3.71733031252663248624409938613E-15;
    w[92] = 3.33444143033026256341061235315E-16;
    w[93] = 2.69865119072980851232572568063E-17;
    w[94] = 1.96593268885070384943390296306E-18;
    w[95] = 1.28579977867628696999762170542E-19;
    w[96] = 7.52931616388155067444192947319E-21;
    w[97] = 3.93559908609832906838466602268E-22;
    w[98] = 1.83031346363374264415878982576E-23;
    w[99] = 7.54682189033203465872349657723E-25;
    w[100] = 2.74824892121260880467531987939E-26;
    w[101] = 8.8019567691972403392314252914E-28;
    w[102] = 2.4679773241854004762148469348E-29;
    w[103] = 6.02777538509463291699314327193E-31;
    w[104] = 1.27533314110484056196532640642E-32;
    w[105] = 2.3232083386415854084664074623E-34;
    w[106] = 3.61918834460649868835433546523E-36;
    w[107] = 4.78543906802804099967221020647E-38;
    w[108] = 5.32553037755907921458489847863E-40;
    w[109] = 4.9411031115925407477456893331E-42;
    w[110] = 3.78144279409152203964384313149E-44;
    w[111] = 2.35816591560823143778744566357E-46;
    w[112] = 1.18170544407210392716367665268E-48;
    w[113] = 4.68218083833618292793410025836E-51;
    w[114] = 1.4394641949298720336141068619E-53;
    w[115] = 3.3581166223962736386929935773E-56;
    w[116] = 5.78845633750656959788340019085E-59;
    w[117] = 7.13750929465743002965122123123E-62;
    w[118] = 6.04865489642049179016214753843E-65;
    w[119] = 3.34891390118992716444169809114E-68;
    w[120] = 1.13419242086298913875376620343E-71;
    w[121] = 2.15037147336771602203551878273E-75;
    w[122] = 2.01439576527109443920782513994E-79;
    w[123] = 7.7306185241134158744827181222E-84;
    w[124] = 8.93216815722645216635320162557E-89;
    w[125] = 1.72727980594728851329952877284E-94;
    w[126] = 1.25044975770895101066558695394E-101;
  }
  else if ( n == 128 )
  {
    x[0] = -15.2918197668827409717467886552E+00;
    x[1] = -14.7338424735892990556131447177E+00;
    x[2] = -14.2739813047878355625094356564E+00;
    x[3] = -13.8652069844762415897768433203E+00;
    x[4] = -13.4895564126231418263791177750E+00;
    x[5] = -13.1377747880276511010650586719E+00;
    x[6] = -12.8043120820671312950137141654E+00;
    x[7] = -12.4855125853494481606990566084E+00;
    x[8] = -12.1788086198312463132740693095E+00;
    x[9] = -11.8823101188783115808359168093E+00;
    x[10] = -11.5945750547414517467820845908E+00;
    x[11] = -11.3144716442899779172120028451E+00;
    x[12] = -11.0410909760196333842428986719E+00;
    x[13] = -10.7736891151614406713116609896E+00;
    x[14] = -10.5116473299148686173941369279E+00;
    x[15] = -10.2544439284709307170245436605E+00;
    x[16] = -10.0016337989301228460111363000E+00;
    x[17] = -9.75283321343916867454942614151E+00;
    x[18] = -9.50770832327905657695490182674E+00;
    x[19] = -9.26596630029617592185364037517E+00;
    x[20] = -9.02734841339478834482665573280E+00;
    x[21] = -8.79162454488868640635040291427E+00;
    x[22] = -8.55858879506450828896030380072E+00;
    x[23] = -8.32805592079014664500802003672E+00;
    x[24] = -8.09985842150789607545794348110E+00;
    x[25] = -7.87384413353543446678710891884E+00;
    x[26] = -7.64987422768100656113184995327E+00;
    x[27] = -7.42782152995230111565739552073E+00;
    x[28] = -7.20756910338733385441779947109E+00;
    x[29] = -6.98900904264477401185223744438E+00;
    x[30] = -6.77204144325592885820688621877E+00;
    x[31] = -6.55657351526448288962578894289E+00;
    x[32] = -6.34251881700177947172938955573E+00;
    x[33] = -6.12979658942216202462059597877E+00;
    x[34] = -5.91833117508581167511681743446E+00;
    x[35] = -5.70805150876808626177490879113E+00;
    x[36] = -5.49889066897390948452218926009E+00;
    x[37] = -5.29078548147717957643674180866E+00;
    x[38] = -5.08367616748933990673505368300E+00;
    x[39] = -4.87750603026481441216755173491E+00;
    x[40] = -4.67222117493263892214567470373E+00;
    x[41] = -4.46777025714858268344631831723E+00;
    x[42] = -4.26410425682551915674979043600E+00;
    x[43] = -4.06117627374927282427754765790E+00;
    x[44] = -3.85894134234428182659062673118E+00;
    x[45] = -3.65735626323530809623058740618E+00;
    x[46] = -3.45637944957173748220943445337E+00;
    x[47] = -3.25597078635065934665290567700E+00;
    x[48] = -3.05609150120268005595784091684E+00;
    x[49] = -2.85670404529740528265184910544E+00;
    x[50] = -2.65777198318948399631081621992E+00;
    x[51] = -2.45925989056573940193677619199E+00;
    x[52] = -2.26113325897306228028420817752E+00;
    x[53] = -2.06335840670856597768175136750E+00;
    x[54] = -1.86590239514059869664912407799E+00;
    x[55] = -1.66873294980372363048660121191E+00;
    x[56] = -1.47181838567448600067837560546E+00;
    x[57] = -1.27512753608915832143251082623E+00;
    x[58] = -1.07862968481090893047100757570E+00;
    x[59] = -0.882294500792981406000508343227E+00;
    x[60] = -0.686091975217334872045286432691E+00;
    x[61] = -0.489992360415458918089044385637E+00;
    x[62] = -0.293966110300295702813351867404E+00;
    x[63] = -0.0979838219558189543137713246862E+00;
    x[64] = 0.0979838219558189543137713246862E+00;
    x[65] = 0.293966110300295702813351867404E+00;
    x[66] = 0.489992360415458918089044385637E+00;
    x[67] = 0.686091975217334872045286432691E+00;
    x[68] = 0.882294500792981406000508343227E+00;
    x[69] = 1.07862968481090893047100757570E+00;
    x[70] = 1.27512753608915832143251082623E+00;
    x[71] = 1.47181838567448600067837560546E+00;
    x[72] = 1.66873294980372363048660121191E+00;
    x[73] = 1.86590239514059869664912407799E+00;
    x[74] = 2.06335840670856597768175136750E+00;
    x[75] = 2.26113325897306228028420817752E+00;
    x[76] = 2.45925989056573940193677619199E+00;
    x[77] = 2.65777198318948399631081621992E+00;
    x[78] = 2.85670404529740528265184910544E+00;
    x[79] = 3.05609150120268005595784091684E+00;
    x[80] = 3.25597078635065934665290567700E+00;
    x[81] = 3.45637944957173748220943445337E+00;
    x[82] = 3.65735626323530809623058740618E+00;
    x[83] = 3.85894134234428182659062673118E+00;
    x[84] = 4.06117627374927282427754765790E+00;
    x[85] = 4.26410425682551915674979043600E+00;
    x[86] = 4.46777025714858268344631831723E+00;
    x[87] = 4.67222117493263892214567470373E+00;
    x[88] = 4.87750603026481441216755173491E+00;
    x[89] = 5.08367616748933990673505368300E+00;
    x[90] = 5.29078548147717957643674180866E+00;
    x[91] = 5.49889066897390948452218926009E+00;
    x[92] = 5.70805150876808626177490879113E+00;
    x[93] = 5.91833117508581167511681743446E+00;
    x[94] = 6.12979658942216202462059597877E+00;
    x[95] = 6.34251881700177947172938955573E+00;
    x[96] = 6.55657351526448288962578894289E+00;
    x[97] = 6.77204144325592885820688621877E+00;
    x[98] = 6.98900904264477401185223744438E+00;
    x[99] = 7.20756910338733385441779947109E+00;
    x[100] = 7.42782152995230111565739552073E+00;
    x[101] = 7.64987422768100656113184995327E+00;
    x[102] = 7.87384413353543446678710891884E+00;
    x[103] = 8.09985842150789607545794348110E+00;
    x[104] = 8.32805592079014664500802003672E+00;
    x[105] = 8.55858879506450828896030380072E+00;
    x[106] = 8.79162454488868640635040291427E+00;
    x[107] = 9.02734841339478834482665573280E+00;
    x[108] = 9.26596630029617592185364037517E+00;
    x[109] = 9.50770832327905657695490182674E+00;
    x[110] = 9.75283321343916867454942614151E+00;
    x[111] = 10.0016337989301228460111363000E+00;
    x[112] = 10.2544439284709307170245436605E+00;
    x[113] = 10.5116473299148686173941369279E+00;
    x[114] = 10.7736891151614406713116609896E+00;
    x[115] = 11.0410909760196333842428986719E+00;
    x[116] = 11.3144716442899779172120028451E+00;
    x[117] = 11.5945750547414517467820845908E+00;
    x[118] = 11.8823101188783115808359168093E+00;
    x[119] = 12.1788086198312463132740693095E+00;
    x[120] = 12.4855125853494481606990566084E+00;
    x[121] = 12.8043120820671312950137141654E+00;
    x[122] = 13.1377747880276511010650586719E+00;
    x[123] = 13.4895564126231418263791177750E+00;
    x[124] = 13.8652069844762415897768433203E+00;
    x[125] = 14.2739813047878355625094356564E+00;
    x[126] = 14.7338424735892990556131447177E+00;
    x[127] = 15.2918197668827409717467886552E+00;

    w[0] = 1.79906598010928472082336338805E-102;
    w[1] = 2.60817240240911107924885148459E-95;
    w[2] = 1.40468977131508863479865725345E-89;
    w[3] = 1.2612494833385383033093221663E-84;
    w[4] = 3.4012300869366371268669286673E-80;
    w[5] = 3.75121586880472499656274624235E-76;
    w[6] = 2.04158579724398501580069247379E-72;
    w[7] = 6.21424416183031366240930730224E-69;
    w[8] = 1.15615516409637521334725409468E-65;
    w[9] = 1.40446725774048726044186592003E-62;
    w[10] = 1.17197850121298051738559888373E-59;
    w[11] = 6.9930729240519559879874741506E-57;
    w[12] = 3.08207738333929868710425541163E-54;
    w[13] = 1.03048625205569473422672330856E-51;
    w[14] = 2.67274375173606785452021989916E-49;
    w[15] = 5.48021702897879649820616283051E-47;
    w[16] = 9.02804013878664400917961564574E-45;
    w[17] = 1.21177953413059190735434940091E-42;
    w[18] = 1.34149748176436936696556841563E-40;
    w[19] = 1.2380855579763680376188381669E-38;
    w[20] = 9.6167080679675069775952182446E-37;
    w[21] = 6.33991352636648906076753997388E-35;
    w[22] = 3.57437889587942107216457034803E-33;
    w[23] = 1.73510302028206120881601688138E-31;
    w[24] = 7.29654500676840425381868704321E-30;
    w[25] = 2.67292362005807324017266437183E-28;
    w[26] = 8.5728304837693537445493254974E-27;
    w[27] = 2.41840345964766496960390574396E-25;
    w[28] = 6.02598403200645428864656987226E-24;
    w[29] = 1.33136785903358960440599429474E-22;
    w[30] = 2.61745758393481115586873166674E-21;
    w[31] = 4.59400767732972159221172605389E-20;
    w[32] = 7.22010731692829201964437734131E-19;
    w[33] = 1.01893323042329252403658204469E-17;
    w[34] = 1.2945481593393715343954569556E-16;
    w[35] = 1.4842238375138564829118955689E-15;
    w[36] = 1.53904973035354581424981070383E-14;
    w[37] = 1.44634732119041656320590928428E-13;
    w[38] = 1.23421448660055669081623604437E-12;
    w[39] = 9.58031650873585770862066358548E-12;
    w[40] = 6.77578048777455378630839649193E-11;
    w[41] = 4.37318665984840344563217253619E-10;
    w[42] = 2.57939722942639480114980569527E-9;
    w[43] = 1.39219071529351788119578816175E-8;
    w[44] = 6.88458112215435009064406266312E-8;
    w[45] = 3.12287298617890308197944991751E-7;
    w[46] = 1.30074700323819923351375586698E-6;
    w[47] = 4.97992453259098701134099270598E-6;
    w[48] = 0.0000175404858480939050383677619791E+00;
    w[49] = 0.0000568874376004024109270187885882E+00;
    w[50] = 0.000170014088262809409409897155763E+00;
    w[51] = 0.000468551537808411365479802126842E+00;
    w[52] = 0.00119156381445716723911680561041E+00;
    w[53] = 0.00279783940160578927319080368252E+00;
    w[54] = 0.00606886240692588762066801419927E+00;
    w[55] = 0.0121669188644693394910166328856E+00;
    w[56] = 0.0225543101678244224102498222492E+00;
    w[57] = 0.0386739548106369026550248867136E+00;
    w[58] = 0.061360721004490065664651069257E+00;
    w[59] = 0.090108678376448919548057439804E+00;
    w[60] = 0.122503273164135694618664605611E+00;
    w[61] = 0.154210435298354383363527713284E+00;
    w[62] = 0.179773083907799264988697956102E+00;
    w[63] = 0.194097611864087756977697028723E+00;
    w[64] = 0.194097611864087756977697028723E+00;
    w[65] = 0.179773083907799264988697956102E+00;
    w[66] = 0.154210435298354383363527713284E+00;
    w[67] = 0.122503273164135694618664605611E+00;
    w[68] = 0.090108678376448919548057439804E+00;
    w[69] = 0.061360721004490065664651069257E+00;
    w[70] = 0.0386739548106369026550248867136E+00;
    w[71] = 0.0225543101678244224102498222492E+00;
    w[72] = 0.0121669188644693394910166328856E+00;
    w[73] = 0.00606886240692588762066801419927E+00;
    w[74] = 0.00279783940160578927319080368252E+00;
    w[75] = 0.00119156381445716723911680561041E+00;
    w[76] = 0.000468551537808411365479802126842E+00;
    w[77] = 0.000170014088262809409409897155763E+00;
    w[78] = 0.0000568874376004024109270187885882E+00;
    w[79] = 0.0000175404858480939050383677619791E+00;
    w[80] = 4.97992453259098701134099270598E-06;
    w[81] = 1.30074700323819923351375586698E-06;
    w[82] = 3.12287298617890308197944991751E-07;
    w[83] = 6.88458112215435009064406266312E-08;
    w[84] = 1.39219071529351788119578816175E-08;
    w[85] = 2.57939722942639480114980569527E-09;
    w[86] = 4.37318665984840344563217253619E-10;
    w[87] = 6.77578048777455378630839649193E-11;
    w[88] = 9.58031650873585770862066358548E-12;
    w[89] = 1.23421448660055669081623604437E-12;
    w[90] = 1.44634732119041656320590928428E-13;
    w[91] = 1.53904973035354581424981070383E-14;
    w[92] = 1.4842238375138564829118955689E-15;
    w[93] = 1.2945481593393715343954569556E-16;
    w[94] = 1.01893323042329252403658204469E-17;
    w[95] = 7.22010731692829201964437734131E-19;
    w[96] = 4.59400767732972159221172605389E-20;
    w[97] = 2.61745758393481115586873166674E-21;
    w[98] = 1.33136785903358960440599429474E-22;
    w[99] = 6.02598403200645428864656987226E-24;
    w[100] = 2.41840345964766496960390574396E-25;
    w[101] = 8.5728304837693537445493254974E-27;
    w[102] = 2.67292362005807324017266437183E-28;
    w[103] = 7.29654500676840425381868704321E-30;
    w[104] = 1.73510302028206120881601688138E-31;
    w[105] = 3.57437889587942107216457034803E-33;
    w[106] = 6.33991352636648906076753997388E-35;
    w[107] = 9.6167080679675069775952182446E-37;
    w[108] = 1.2380855579763680376188381669E-38;
    w[109] = 1.34149748176436936696556841563E-40;
    w[110] = 1.21177953413059190735434940091E-42;
    w[111] = 9.02804013878664400917961564574E-45;
    w[112] = 5.48021702897879649820616283051E-47;
    w[113] = 2.67274375173606785452021989916E-49;
    w[114] = 1.03048625205569473422672330856E-51;
    w[115] = 3.08207738333929868710425541163E-54;
    w[116] = 6.9930729240519559879874741506E-57;
    w[117] = 1.17197850121298051738559888373E-59;
    w[118] = 1.40446725774048726044186592003E-62;
    w[119] = 1.15615516409637521334725409468E-65;
    w[120] = 6.21424416183031366240930730224E-69;
    w[121] = 2.04158579724398501580069247379E-72;
    w[122] = 3.75121586880472499656274624235E-76;
    w[123] = 3.4012300869366371268669286673E-80;
    w[124] = 1.2612494833385383033093221663E-84;
    w[125] = 1.40468977131508863479865725345E-89;
    w[126] = 2.60817240240911107924885148459E-95;
    w[127] = 1.79906598010928472082336338805E-102;
  }
  else if ( n == 129 )
  {
    x[0] = -15.3550496746831285549167746019E+00;
    x[1] = -14.7978308964903080628845608050E+00;
    x[2] = -14.3386115290089672811362217078E+00;
    x[3] = -13.9304208664791805435265533989E+00;
    x[4] = -13.5553179661308567022816946453E+00;
    x[5] = -13.2040593596741921982903144147E+00;
    x[6] = -12.8711017789036282758938040681E+00;
    x[7] = -12.5527939524445397878411009214E+00;
    x[8] = -12.2465713150240016840404064819E+00;
    x[9] = -11.9505460927691823148587203418E+00;
    x[10] = -11.6632780120689523111895976521E+00;
    x[11] = -11.3836366735119364919041401601E+00;
    x[12] = -11.1107142851416459382067369906E+00;
    x[13] = -10.8437678377155441232588867872E+00;
    x[14] = -10.5821793789735138638177686355E+00;
    x[15] = -10.3254278845735933555309803756E+00;
    x[16] = -10.0730688225840385168071109595E+00;
    x[17] = -9.82471897583315292981163227664E+00;
    x[18] = -9.58004495076523054718996925368E+00;
    x[19] = -9.33875432946329683313144773753E+00;
    x[20] = -9.10058875441877630705419698871E+00;
    x[21] = -8.86531845144460445006059884245E+00;
    x[22] = -8.63273783950843552405601767759E+00;
    x[23] = -8.40266197362535012572889592790E+00;
    x[24] = -8.17492363437398978801516910338E+00;
    x[25] = -7.94937092512605027069441566240E+00;
    x[26] = -7.72586527212128476858225880726E+00;
    x[27] = -7.50427974726352240601473698475E+00;
    x[28] = -7.28449765174015372725258835236E+00;
    x[29] = -7.06641131216039069912389691607E+00;
    x[30] = -6.84992105116014925339717178568E+00;
    x[31] = -6.63493430223598850302862096202E+00;
    x[32] = -6.42136484458510366813184121466E+00;
    x[33] = -6.20913213839958724139275362341E+00;
    x[34] = -5.99816074472179863235247556956E+00;
    x[35] = -5.78837981685588271189500366573E+00;
    x[36] = -5.57972265262736517721195076411E+00;
    x[37] = -5.37212629862206810544406569124E+00;
    x[38] = -5.16553119901808798749445925424E+00;
    x[39] = -4.95988088282680494441019859192E+00;
    x[40] = -4.75512168433945660698412105871E+00;
    x[41] = -4.55120249237992751786552441506E+00;
    x[42] = -4.34807452462720039287086617988E+00;
    x[43] = -4.14569112381985885943227755134E+00;
    x[44] = -3.94400757311174034591077592589E+00;
    x[45] = -3.74298092822942697020355662873E+00;
    x[46] = -3.54256986440235708810942555303E+00;
    x[47] = -3.34273453630584653564135188865E+00;
    x[48] = -3.14343644948499981855157130224E+00;
    x[49] = -2.94463834192045889859029200431E+00;
    x[50] = -2.74630407456093830886162508356E+00;
    x[51] = -2.54839852978722422234239962660E+00;
    x[52] = -2.35088751689161374336515868011E+00;
    x[53] = -2.15373768375879465132813086025E+00;
    x[54] = -1.95691643402151503173496953711E+00;
    x[55] = -1.76039184903921457004681976795E+00;
    x[56] = -1.56413261411186298656139617983E+00;
    x[57] = -1.36810794839605047973791359079E+00;
    x[58] = -1.17228753803712318264391216760E+00;
    x[59] = -0.976641472070867557126534700881E+00;
    x[60] = -0.781140180681760289238140546741E+00;
    x[61] = -0.585754375432805697119652981369E+00;
    x[62] = -0.390454991105046004780513867383E+00;
    x[63] = -0.195213128803407573801607754230E+00;
    x[64] = 0.0E+00;
    x[65] = 0.195213128803407573801607754230E+00;
    x[66] = 0.390454991105046004780513867383E+00;
    x[67] = 0.585754375432805697119652981369E+00;
    x[68] = 0.781140180681760289238140546741E+00;
    x[69] = 0.976641472070867557126534700881E+00;
    x[70] = 1.17228753803712318264391216760E+00;
    x[71] = 1.36810794839605047973791359079E+00;
    x[72] = 1.56413261411186298656139617983E+00;
    x[73] = 1.76039184903921457004681976795E+00;
    x[74] = 1.95691643402151503173496953711E+00;
    x[75] = 2.15373768375879465132813086025E+00;
    x[76] = 2.35088751689161374336515868011E+00;
    x[77] = 2.54839852978722422234239962660E+00;
    x[78] = 2.74630407456093830886162508356E+00;
    x[79] = 2.94463834192045889859029200431E+00;
    x[80] = 3.14343644948499981855157130224E+00;
    x[81] = 3.34273453630584653564135188865E+00;
    x[82] = 3.54256986440235708810942555303E+00;
    x[83] = 3.74298092822942697020355662873E+00;
    x[84] = 3.94400757311174034591077592589E+00;
    x[85] = 4.14569112381985885943227755134E+00;
    x[86] = 4.34807452462720039287086617988E+00;
    x[87] = 4.55120249237992751786552441506E+00;
    x[88] = 4.75512168433945660698412105871E+00;
    x[89] = 4.95988088282680494441019859192E+00;
    x[90] = 5.16553119901808798749445925424E+00;
    x[91] = 5.37212629862206810544406569124E+00;
    x[92] = 5.57972265262736517721195076411E+00;
    x[93] = 5.78837981685588271189500366573E+00;
    x[94] = 5.99816074472179863235247556956E+00;
    x[95] = 6.20913213839958724139275362341E+00;
    x[96] = 6.42136484458510366813184121466E+00;
    x[97] = 6.63493430223598850302862096202E+00;
    x[98] = 6.84992105116014925339717178568E+00;
    x[99] = 7.06641131216039069912389691607E+00;
    x[100] = 7.28449765174015372725258835236E+00;
    x[101] = 7.50427974726352240601473698475E+00;
    x[102] = 7.72586527212128476858225880726E+00;
    x[103] = 7.94937092512605027069441566240E+00;
    x[104] = 8.17492363437398978801516910338E+00;
    x[105] = 8.40266197362535012572889592790E+00;
    x[106] = 8.63273783950843552405601767759E+00;
    x[107] = 8.86531845144460445006059884245E+00;
    x[108] = 9.10058875441877630705419698871E+00;
    x[109] = 9.33875432946329683313144773753E+00;
    x[110] = 9.58004495076523054718996925368E+00;
    x[111] = 9.82471897583315292981163227664E+00;
    x[112] = 10.0730688225840385168071109595E+00;
    x[113] = 10.3254278845735933555309803756E+00;
    x[114] = 10.5821793789735138638177686355E+00;
    x[115] = 10.8437678377155441232588867872E+00;
    x[116] = 11.1107142851416459382067369906E+00;
    x[117] = 11.3836366735119364919041401601E+00;
    x[118] = 11.6632780120689523111895976521E+00;
    x[119] = 11.9505460927691823148587203418E+00;
    x[120] = 12.2465713150240016840404064819E+00;
    x[121] = 12.5527939524445397878411009214E+00;
    x[122] = 12.8711017789036282758938040681E+00;
    x[123] = 13.2040593596741921982903144147E+00;
    x[124] = 13.5553179661308567022816946453E+00;
    x[125] = 13.9304208664791805435265533989E+00;
    x[126] = 14.3386115290089672811362217078E+00;
    x[127] = 14.7978308964903080628845608050E+00;
    x[128] = 15.3550496746831285549167746019E+00;

    w[0] = 2.58755395082114927399064139631E-103;
    w[1] = 3.93601845908067608811461078697E-96;
    w[2] = 2.20725529577484588586177997021E-90;
    w[3] = 2.05563087297774646200941835216E-85;
    w[4] = 5.73584763407311509769038083955E-81;
    w[5] = 6.53456499014096713882711627986E-77;
    w[6] = 3.66903606454555600244832281797E-73;
    w[7] = 1.15105101975113879079427442365E-69;
    w[8] = 2.20553774145133363585421051568E-66;
    w[9] = 2.75763663311195533797446164671E-63;
    w[10] = 2.36731747071610805241477009401E-60;
    w[11] = 1.45257860403230704544333281907E-57;
    w[12] = 6.58119121529392093666305170751E-55;
    w[13] = 2.26137732951303228667152914802E-52;
    w[14] = 6.02643011329776195986432204924E-50;
    w[15] = 1.26938407638088455004457398255E-47;
    w[16] = 2.14791778799787733305388107076E-45;
    w[17] = 2.96092183537878053158423564486E-43;
    w[18] = 3.36616090532826422441501486485E-41;
    w[19] = 3.19014783528482307711547124192E-39;
    w[20] = 2.54439796712780366695746038013E-37;
    w[21] = 1.72239465322100711581154624691E-35;
    w[22] = 9.97105538735197785176257533162E-34;
    w[23] = 4.97009943352064894027841342072E-32;
    w[24] = 2.14620630607238052957787041268E-30;
    w[25] = 8.07377921555793000987040030256E-29;
    w[26] = 2.65936924028161851577004868287E-27;
    w[27] = 7.70515053183270746145645250031E-26;
    w[28] = 1.97204966381589933881729892459E-24;
    w[29] = 4.47579713475437089012921294273E-23;
    w[30] = 9.0403370335874459959906960673E-22;
    w[31] = 1.63036407035294103578729410788E-20;
    w[32] = 2.63320491826449443345354482912E-19;
    w[33] = 3.81944198027838553902522199764E-18;
    w[34] = 4.98833273307808866457667338365E-17;
    w[35] = 5.88022772755071728094452091844E-16;
    w[36] = 6.27023947714728862011748531319E-15;
    w[37] = 6.06072571359080078068964155295E-14;
    w[38] = 5.32049070753884105044682362639E-13;
    w[39] = 4.24955065877498808023415505556E-12;
    w[40] = 3.09330203932473692244204789801E-11;
    w[41] = 2.05524352987860630455845773203E-10;
    w[42] = 1.2482251818681431772545606389E-09;
    w[43] = 6.93896714453731562418029048785E-09;
    w[44] = 3.53518234605234028369262582274E-08;
    w[45] = 1.65252704577539544523562160076E-07;
    w[46] = 7.09535030601389014268064639021E-07;
    w[47] = 2.80106033677073567813925250808E-06;
    w[48] = 0.0000101764715414468349837840217278E+00;
    w[49] = 0.0000340541841724020078441933069804E+00;
    w[50] = 0.000105047486997647220847004754607E+00;
    w[51] = 0.000298922505941519029186629138321E+00;
    w[52] = 0.000785197220610268688197653195861E+00;
    w[53] = 0.00190506673927936544347937172051E+00;
    w[54] = 0.00427162074179231114560048384305E+00;
    w[55] = 0.00885609926394363549300290104701E+00;
    w[56] = 0.0169845117091580731620255503875E+00;
    w[57] = 0.0301436534848915408822025025241E+00;
    w[58] = 0.0495245901368945546436264039895E+00;
    w[59] = 0.0753454506416603410859342903275E+00;
    w[60] = 0.106172669789632918045101372957E+00;
    w[61] = 0.138604146980788427972651263139E+00;
    w[62] = 0.167654732143619067997522798882E+00;
    w[63] = 0.187923095463858179335367694936E+00;
    w[64] = 0.195208341719164170910088609838E+00;
    w[65] = 0.187923095463858179335367694936E+00;
    w[66] = 0.167654732143619067997522798882E+00;
    w[67] = 0.138604146980788427972651263139E+00;
    w[68] = 0.106172669789632918045101372957E+00;
    w[69] = 0.0753454506416603410859342903275E+00;
    w[70] = 0.0495245901368945546436264039895E+00;
    w[71] = 0.0301436534848915408822025025241E+00;
    w[72] = 0.0169845117091580731620255503875E+00;
    w[73] = 0.00885609926394363549300290104701E+00;
    w[74] = 0.00427162074179231114560048384305E+00;
    w[75] = 0.00190506673927936544347937172051E+00;
    w[76] = 0.000785197220610268688197653195861E+00;
    w[77] = 0.000298922505941519029186629138321E+00;
    w[78] = 0.000105047486997647220847004754607E+00;
    w[79] = 0.0000340541841724020078441933069804E+00;
    w[80] = 0.0000101764715414468349837840217278E+00;
    w[81] = 2.80106033677073567813925250808E-06;
    w[82] = 7.09535030601389014268064639021E-07;
    w[83] = 1.65252704577539544523562160076E-07;
    w[84] = 3.53518234605234028369262582274E-08;
    w[85] = 6.93896714453731562418029048785E-09;
    w[86] = 1.2482251818681431772545606389E-09;
    w[87] = 2.05524352987860630455845773203E-10;
    w[88] = 3.09330203932473692244204789801E-11;
    w[89] = 4.24955065877498808023415505556E-12;
    w[90] = 5.32049070753884105044682362639E-13;
    w[91] = 6.06072571359080078068964155295E-14;
    w[92] = 6.27023947714728862011748531319E-15;
    w[93] = 5.88022772755071728094452091844E-16;
    w[94] = 4.98833273307808866457667338365E-17;
    w[95] = 3.81944198027838553902522199764E-18;
    w[96] = 2.63320491826449443345354482912E-19;
    w[97] = 1.63036407035294103578729410788E-20;
    w[98] = 9.0403370335874459959906960673E-22;
    w[99] = 4.47579713475437089012921294273E-23;
    w[100] = 1.97204966381589933881729892459E-24;
    w[101] = 7.70515053183270746145645250031E-26;
    w[102] = 2.65936924028161851577004868287E-27;
    w[103] = 8.07377921555793000987040030256E-29;
    w[104] = 2.14620630607238052957787041268E-30;
    w[105] = 4.97009943352064894027841342072E-32;
    w[106] = 9.97105538735197785176257533162E-34;
    w[107] = 1.72239465322100711581154624691E-35;
    w[108] = 2.54439796712780366695746038013E-37;
    w[109] = 3.19014783528482307711547124192E-39;
    w[110] = 3.36616090532826422441501486485E-41;
    w[111] = 2.96092183537878053158423564486E-43;
    w[112] = 2.14791778799787733305388107076E-45;
    w[113] = 1.26938407638088455004457398255E-47;
    w[114] = 6.02643011329776195986432204924E-50;
    w[115] = 2.26137732951303228667152914802E-52;
    w[116] = 6.58119121529392093666305170751E-55;
    w[117] = 1.45257860403230704544333281907E-57;
    w[118] = 2.36731747071610805241477009401E-60;
    w[119] = 2.75763663311195533797446164671E-63;
    w[120] = 2.20553774145133363585421051568E-66;
    w[121] = 1.15105101975113879079427442365E-69;
    w[122] = 3.66903606454555600244832281797E-73;
    w[123] = 6.53456499014096713882711627986E-77;
    w[124] = 5.73584763407311509769038083955E-81;
    w[125] = 2.05563087297774646200941835216E-85;
    w[126] = 2.20725529577484588586177997021E-90;
    w[127] = 3.93601845908067608811461078697E-96;
    w[128] = 2.58755395082114927399064139631E-103;
  }
  else
  {
    cerr << "\n";
    cerr << "HERMITE_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1:20, 31/32/33, 63/64/65, 127/128/129.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void hermite_ss_compute ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_SS_COMPUTE computes a Gauss-Hermite quadrature rule.
//
//  Discussion:
//
//    The abscissas are the zeros of the N-th order Hermite polynomial.
//
//    The integral:
//
//      Integral ( -oo < X < +oo ) exp ( - X*X ) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the formula to be computed.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double cc;
  double dp2;
  int i;
  double p1;
  double pi = 3.141592653589793;
  double s;
  double x;

  cc = sqrt ( pi ) * gamma ( ( double ) ( order ) ) 
    / pow ( 2.0, order - 1 );

  s = pow ( 2.0 * ( double ) ( order ) + 1.0, 1.0 / 6.0 );

  for ( i = 0; i < ( order + 1 ) / 2; i++ )
  {
    if ( i == 0 )
    {
      x = s * s * s - 1.85575 / s;
    }
    else if ( i == 1 )
    {
      x = x - 1.14 * pow ( ( double ) ( order ), 0.426 ) / x;
    }
    else if ( i == 2 )
    {
      x = 1.86 * x - 0.86 * xtab[0];
    }
    else if ( i == 3 )
    {
      x = 1.91 * x - 0.91 * xtab[1];
    }
    else
    {
      x = 2.0 * x - xtab[i-2];
    }

    hermite_ss_root ( &x, order, &dp2, &p1 );

    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;

    xtab[order-i-1] = -x;
    weight[order-i-1] = weight[i];
  }
//
//  Reverse the order of the values.
//
  for ( i = 0; i < ( order / 2 ); i++ )
  {
    x               = xtab[i];
    xtab[i]         = xtab[order-1-i];
    xtab[order-1-i] = x;
  }

  return;
}
//****************************************************************************80

void hermite_ss_recur ( double *p2, double *dp2, double *p1, double x, int order )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_SS_RECUR finds the value and derivative of a Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of H(ORDER)(X).
//
//    Output, double *DP2, the value of H'(ORDER)(X).
//
//    Output, double *P1, the value of H(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
{
  int i;
  double dq0;
  double dq1;
  double dq2;
  double q0;
  double q1;
  double q2;

  q1 = 1.0;
  dq1 = 0.0;

  q2 = x;
  dq2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    q0 = q1;
    dq0 = dq1;

    q1 = q2;
    dq1 = dq2;

    q2  = x * q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * q0;
    dq2 = x * dq1 + q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * dq0;
  }

  *p2 = q2;
  *dp2 = dq2;
  *p1 = q1;

  return;
}
//****************************************************************************80

void hermite_ss_root ( double *x, int order, double *dp2, double *p1 )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_SS_ROOT improves an approximate root of a Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the Hermite polynomial.
//
//    Output, double *DP2, the value of H'(ORDER)(X).
//
//    Output, double *P1, the value of H(ORDER-1)(X).
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    hermite_ss_recur ( &p2, dp2, p1, *x, order );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }

  return;
}
//****************************************************************************80

int i4_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL2 computes the double factorial function.
//
//  Discussion:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    N!!
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial function.
//    If N is less than 1, I4_FACTORIAL2 is returned as 1.
//
//    Output, int I4_FACTORIAL2, the value of the double factorial of N.
//
{
  int value;

  if ( n < 1 )
  {
    return 1;
  }

  value = 1;

  while ( 1 < n )
  {
    value = value * n;
    n = n - 2;
  }

  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
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
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

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
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
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
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
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

void imtqlx ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentially) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( r8_abs ( e[m-1] ) <= prec * ( r8_abs ( d[m-1] ) + r8_abs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        cerr << "\n";
        cerr << "IMTQLX - Fatal error!\n";
        cerr << "  Iteration limit exceeded\n";
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r = sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + r8_abs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( r8_abs ( g ) <= r8_abs ( f ) )
        {
          c = g / f;
          r = sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r = sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
//****************************************************************************80

void jacobi_ek_compute ( int n, double alpha, double beta, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_EK_COMPUTE: Elhay-Kautsky method for Gauss-Jacobi quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) (1-X)**ALPHA * (1+X)**BETA * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, double ALPHA, BETA, the exponents of (1-X) and
//    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
//    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double abi;
  double *bj;
  int i;
  double i_r8;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = pow ( 2.0, alpha + beta + 1.0 )
    * r8_gamma ( alpha + 1.0 ) 
    * r8_gamma ( beta + 1.0 ) 
    / r8_gamma ( 2.0 + alpha + beta );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];

  x[0] = ( beta - alpha ) / ( 2.0 + alpha + beta );

  bj[0] = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
    / ( ( 3.0 + alpha + beta ) 
      * ( 2.0 + alpha + beta ) * ( 2.0 + alpha + beta ) );

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    abi = 2.0 * i_r8 + alpha + beta;
    x[i] = ( beta + alpha ) * ( beta - alpha ) / ( ( abi - 2.0 ) * abi );
    bj[i] = 4.0 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) 
      * ( i_r8 + alpha + beta ) 
      / ( ( abi - 1.0 ) * ( abi + 1.0 ) * abi * abi );
  }

  for ( i = 0; i < n; i++ )
  {
    bj[i] = sqrt ( bj[i] );
  }

  w[0] = sqrt ( zemu );

  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double jacobi_integral ( int expon, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_INTEGRAL evaluates the integral of a monomial with Jacobi weight.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 September 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Input, double ALPHA, the exponent of (1-X) in the weight factor.
//
//    Input, double BETA, the exponent of (1+X) in the weight factor.
//
//    Output, double JACOBI_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;
  double value2;

  c = ( double ) ( expon );

  if ( ( expon % 2 ) == 0 )
  {
    s = +1.0;
  }
  else
  {
    s = -1.0;
  }

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + beta + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  arg1 = - beta;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value2 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = gamma ( 1.0 + c ) * ( 
      s * gamma ( 1.0 + beta  ) * value1 
    / gamma ( 2.0 + beta  + c ) 
    +     gamma ( 1.0 + alpha ) * value2 
    / gamma ( 2.0 + alpha + c ) );

  return value;
}
//****************************************************************************80

void jacobi_ss_compute ( int order, double alpha, double beta, double xtab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SS_COMPUTE computes a Gauss-Jacobi quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) (1-X)**ALPHA * (1+X)**BETA * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//    Thanks to Xu Xiang of Fudan University for pointing out that
//    an earlier implementation of this routine was incorrect!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Input, double ALPHA, BETA, the exponents of (1-X) and
//    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
//    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double an;
  double *b;
  double bn;
  double *c;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double x;
//
//  Check ALPHA and BETA.
//
  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "JACOBI_SS_COMPUTE - Fatal error!\n";
    cerr << "  -1.0 < ALPHA is required.\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "JACOBI_SS_COMPUTE - Fatal error!\n";
    cerr << "  -1.0 < BETA is required.\n";
    exit ( 1 );
  }

  b = new double[order];
  c = new double[order];
//
//  Set the recursion coefficients.
//
  for ( i = 1; i <= order; i++ )
  {
    if ( alpha + beta == 0.0 || beta - alpha == 0.0 )
    {
      b[i-1] = 0.0;
    }
    else
    {
      b[i-1] = ( alpha + beta ) * ( beta - alpha ) / 
             ( ( alpha + beta + ( double ) ( 2 * i ) ) 
             * ( alpha + beta + ( double ) ( 2 * i - 2 ) ) );
    }

    if ( i == 1 )
    {
      c[i-1] = 0.0;
    }
    else
    {
      c[i-1] = 4.0 * ( double ) ( i - 1 ) 
         * ( alpha + ( double ) ( i - 1 ) ) 
          * ( beta + ( double ) ( i - 1 ) ) 
            * ( alpha + beta + ( double ) ( i - 1 ) ) / 
            ( ( alpha + beta + ( double ) ( 2 * i - 1 ) ) 
            * pow ( alpha + beta + ( double ) ( 2 * i - 2 ), 2 ) 
            * ( alpha + beta + ( double ) ( 2 * i - 3 ) ) );
    }
  }

  delta = exp ( lgamma ( alpha        + 1.0 ) 
              + lgamma (         beta + 1.0 ) 
              - lgamma ( alpha + beta + 2.0 ) );
  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + beta + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );
      bn = beta / ( double ) ( order );

      r1 = ( 1.0 + alpha ) 
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) ) 
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 1.48 * an + 0.96 * bn 
        + 0.452 * an * an + 0.83 * an * bn;

      x = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) / 
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) * 
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * beta * 
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x = x - r1 * r2 * r3 * ( 1.0 - x );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * beta / 
        ( ( 6.28 + beta ) * ( double ) ( order * order ) );

      x = x - r1 * r2 * r3 * ( xtab[0] - x );
    }
    else if ( i < order - 1 )
    {
      x = 3.0 * xtab[i-2] - 3.0 * xtab[i-3] + xtab[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta );

      r2 = 1.0 / ( 1.0 + 0.639 
        * ( ( double ) ( order ) - 4.0 ) 
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * 
        ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * beta ) / ( 1.67 + 0.28 * beta );

      r2 = 1.0 / 
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    jacobi_ss_root ( &x, order, alpha, beta, &dp2, &p1, b, c );

    xtab[i-1] = x;
    weight[i-1] = cc / ( dp2 * p1 );
  }
//
//  Reverse the order of the values.
//
  r8vec_reverse ( order, xtab );
  r8vec_reverse ( order, weight );

  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void jacobi_ss_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double beta, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SS_RECUR finds the value and derivative of a Jacobi polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of J(ORDER)(X).
//
//    Output, double *DP2, the value of J'(ORDER)(X).
//
//    Output, double *P1, the value of J(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial.
//
//    Input, double ALPHA, BETA, the exponents of (1-X) and
//    (1+X) in the quadrature rule.
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0 );
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = ( x - b[i-1] ) *  ( *p1 ) - c[i-1] * p0;
    *dp2 = ( x - b[i-1] ) * dp1 + ( *p1 ) - c[i-1] * dp0;
  }
  return;
}
//****************************************************************************80

void jacobi_ss_root ( double *x, int order, double alpha, double beta, 
  double *dp2, double *p1, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SS_ROOT improves an approximate root of a Jacobi polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial.
//
//    Input, double ALPHA, BETA, the exponents of (1-X) and
//    (1+X) in the quadrature rule.
//
//    Output, double *DP2, the value of J'(ORDER)(X).
//
//    Output, double *P1, the value of J(ORDER-1)(X).
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    jacobi_ss_recur ( &p2, dp2, p1, *x, order, alpha, beta, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
//****************************************************************************80

void kronrod_set ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    KRONROD_SET sets abscissas and weights for Gauss-Kronrod quadrature.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAb[I) )
//
//    A Kronrod rule is used in conjunction with a lower order
//    Gauss rule, and provides an efficient error estimation.
//
//    The error may be estimated as the difference in the two integral
//    approximations.
//
//    The efficiency comes about because the Kronrod uses the abscissas
//    of the Gauss rule, thus saving on the number of function evaluations
//    necessary.  If the Kronrod rule were replaced by a Gauss rule of
//    the same order, a higher precision integral estimate would be
//    made, but the function would have to be evaluated at many more
//    points.
//
//    The Gauss Kronrod pair of rules involves an ( ORDER + 1 ) / 2
//    point Gauss-Legendre rule and an ORDER point Kronrod rule.
//    Thus, the 15 point Kronrod rule should be paired with the
//    Gauss-Legendre 7 point rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenger, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK, A Subroutine Package for Automatic Integration,
//    Springer Verlag, 1983.
//
//  Parameters:
//
//    Input, int ORDER, the order, which may be
//    15, 21, 31 or 41, corresponding to Gauss-Legendre rules of
//    order 7, 10, 15 or 20.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 15 )
  {
    xtab[0] =  - 0.9914553711208126E+00;
    xtab[1] =  - 0.9491079123427585E+00;
    xtab[2] =  - 0.8648644233597691E+00;
    xtab[3] =  - 0.7415311855993944E+00;
    xtab[4] =  - 0.5860872354676911E+00;
    xtab[5] =  - 0.4058451513773972E+00;
    xtab[6] =  - 0.2077849550789850E+00;
    xtab[7] =    0.0E+00;
    xtab[8] =    0.2077849550789850E+00;
    xtab[9] =   0.4058451513773972E+00;
    xtab[10] =   0.5860872354676911E+00;
    xtab[11] =   0.7415311855993944E+00;
    xtab[12] =   0.8648644233597691E+00;
    xtab[13] =   0.9491079123427585E+00;
    xtab[14] =   0.9914553711208126E+00;

    weight[0] =  0.2293532201052922E-01;
    weight[1] =  0.6309209262997855E-01;
    weight[2] =  0.1047900103222502E+00;
    weight[3] =  0.1406532597155259E+00;
    weight[4] =  0.1690047266392679E+00;
    weight[5] =  0.1903505780647854E+00;
    weight[6] =  0.2044329400752989E+00;
    weight[7] =  0.2094821410847278E+00;
    weight[8] =  0.2044329400752989E+00;
    weight[9] = 0.1903505780647854E+00;
    weight[10] = 0.1690047266392679E+00;
    weight[11] = 0.1406532597155259E+00;
    weight[12] = 0.1047900103222502E+00;
    weight[13] = 0.6309209262997855E-01;
    weight[14] = 0.2293532201052922E-01;
  }
  else if ( order == 21 )
  {
    xtab[0] =  - 0.9956571630258081E+00;
    xtab[1] =  - 0.9739065285171717E+00;
    xtab[2] =  - 0.9301574913557082E+00;
    xtab[3] =  - 0.8650633666889845E+00;
    xtab[4] =  - 0.7808177265864169E+00;
    xtab[5] =  - 0.6794095682990244E+00;
    xtab[6] =  - 0.5627571346686047E+00;
    xtab[7] =  - 0.4333953941292472E+00;
    xtab[8] =  - 0.2943928627014602E+00;
    xtab[9] = - 0.1488743389816312E+00;
    xtab[10] =   0.0E+00;
    xtab[11] =   0.1488743389816312E+00;
    xtab[12] =   0.2943928627014602E+00;
    xtab[13] =   0.4333953941292472E+00;
    xtab[14] =   0.5627571346686047E+00;
    xtab[15] =   0.6794095682990244E+00;
    xtab[16] =   0.7808177265864169E+00;
    xtab[17] =   0.8650633666889845E+00;
    xtab[18] =   0.9301574913557082E+00;
    xtab[19] =   0.9739065285171717E+00;
    xtab[20] =   0.9956571630258081E+00;

    weight[0] =  0.1169463886737187E-01;
    weight[1] =  0.3255816230796473E-01;
    weight[2] =  0.5475589657435200E-01;
    weight[3] =  0.7503967481091995E-01;
    weight[4] =  0.9312545458369761E-01;
    weight[5] =  0.1093871588022976E+00;
    weight[6] =  0.1234919762620659E+00;
    weight[7] =  0.1347092173114733E+00;
    weight[8] =  0.1427759385770601E+00;
    weight[9] = 0.1477391049013385E+00;
    weight[10] = 0.1494455540029169E+00;
    weight[11] = 0.1477391049013385E+00;
    weight[12] = 0.1427759385770601E+00;
    weight[13] = 0.1347092173114733E+00;
    weight[14] = 0.1234919762620659E+00;
    weight[15] = 0.1093871588022976E+00;
    weight[16] = 0.9312545458369761E-01;
    weight[17] = 0.7503967481091995E-01;
    weight[18] = 0.5475589657435200E-01;
    weight[19] = 0.3255816230796473E-01;
    weight[20] = 0.1169463886737187E-01;
  }
  else if ( order == 31 )
  {
    xtab[0] =  - 0.9980022986933971E+00;
    xtab[1] =  - 0.9879925180204854E+00;
    xtab[2] =  - 0.9677390756791391E+00;
    xtab[3] =  - 0.9372733924007059E+00;
    xtab[4] =  - 0.8972645323440819E+00;
    xtab[5] =  - 0.8482065834104272E+00;
    xtab[6] =  - 0.7904185014424659E+00;
    xtab[7] =  - 0.7244177313601700E+00;
    xtab[8] =  - 0.6509967412974170E+00;
    xtab[9] = - 0.5709721726085388E+00;
    xtab[10] = - 0.4850818636402397E+00;
    xtab[11] = - 0.3941513470775634E+00;
    xtab[12] = - 0.2991800071531688E+00;
    xtab[13] = - 0.2011940939974345E+00;
    xtab[14] = - 0.1011420669187175E+00;
    xtab[15] =   0.0E+00;
    xtab[16] =   0.1011420669187175E+00;
    xtab[17] =   0.2011940939974345E+00;
    xtab[18] =   0.2991800071531688E+00;
    xtab[19] =   0.3941513470775634E+00;
    xtab[20] =   0.4850818636402397E+00;
    xtab[21] =   0.5709721726085388E+00;
    xtab[22] =   0.6509967412974170E+00;
    xtab[23] =   0.7244177313601700E+00;
    xtab[24] =   0.7904185014424659E+00;
    xtab[25] =   0.8482065834104272E+00;
    xtab[26] =   0.8972645323440819E+00;
    xtab[27] =   0.9372733924007059E+00;
    xtab[28] =   0.9677390756791391E+00;
    xtab[29] =   0.9879925180204854E+00;
    xtab[30] =   0.9980022986933971E+00;

    weight[0] =  0.5377479872923349E-02;
    weight[1] =  0.1500794732931612E-01;
    weight[2] =  0.2546084732671532E-01;
    weight[3] =  0.3534636079137585E-01;
    weight[4] =  0.4458975132476488E-01;
    weight[5] =  0.5348152469092809E-01;
    weight[6] =  0.6200956780067064E-01;
    weight[7] =  0.6985412131872826E-01;
    weight[8] =  0.7684968075772038E-01;
    weight[9] = 0.8308050282313302E-01;
    weight[10] = 0.8856444305621177E-01;
    weight[11] = 0.9312659817082532E-01;
    weight[12] = 0.9664272698362368E-01;
    weight[13] = 0.9917359872179196E-01;
    weight[14] = 0.1007698455238756E+00;
    weight[15] = 0.1013300070147915E+00;
    weight[16] = 0.1007698455238756E+00;
    weight[17] = 0.9917359872179196E-01;
    weight[18] = 0.9664272698362368E-01;
    weight[19] = 0.9312659817082532E-01;
    weight[20] = 0.8856444305621177E-01;
    weight[21] = 0.8308050282313302E-01;
    weight[22] = 0.7684968075772038E-01;
    weight[23] = 0.6985412131872826E-01;
    weight[24] = 0.6200956780067064E-01;
    weight[25] = 0.5348152469092809E-01;
    weight[26] = 0.4458975132476488E-01;
    weight[27] = 0.3534636079137585E-01;
    weight[28] = 0.2546084732671532E-01;
    weight[29] = 0.1500794732931612E-01;
    weight[30] = 0.5377479872923349E-02;
  }
  else if ( order == 41 )
  {
    xtab[0] =  - 0.9988590315882777E+00;
    xtab[1] =  - 0.9931285991850949E+00;
    xtab[2] =  - 0.9815078774502503E+00;
    xtab[3] =  - 0.9639719272779138E+00;
    xtab[4] =  - 0.9408226338317548E+00;
    xtab[5] =  - 0.9122344282513259E+00;
    xtab[6] =  - 0.8782768112522820E+00;
    xtab[7] =  - 0.8391169718222188E+00;
    xtab[8] =  - 0.7950414288375512E+00;
    xtab[9] = - 0.7463319064601508E+00;
    xtab[10] = - 0.6932376563347514E+00;
    xtab[11] = - 0.6360536807265150E+00;
    xtab[12] = - 0.5751404468197103E+00;
    xtab[13] = - 0.5108670019508271E+00;
    xtab[14] = - 0.4435931752387251E+00;
    xtab[15] = - 0.3737060887154196E+00;
    xtab[16] = - 0.3016278681149130E+00;
    xtab[17] = - 0.2277858511416451E+00;
    xtab[18] = - 0.1526054652409227E+00;
    xtab[19] = - 0.7652652113349733E-01;
    xtab[20] =   0.0E+00;
    xtab[21] =   0.7652652113349733E-01;
    xtab[22] =   0.1526054652409227E+00;
    xtab[23] =   0.2277858511416451E+00;
    xtab[24] =   0.3016278681149130E+00;
    xtab[25] =   0.3737060887154196E+00;
    xtab[26] =   0.4435931752387251E+00;
    xtab[27] =   0.5108670019508271E+00;
    xtab[28] =   0.5751404468197103E+00;
    xtab[29] =   0.6360536807265150E+00;
    xtab[30] =   0.6932376563347514E+00;
    xtab[31] =   0.7463319064601508E+00;
    xtab[32] =   0.7950414288375512E+00;
    xtab[33] =   0.8391169718222188E+00;
    xtab[34] =   0.8782768112522820E+00;
    xtab[35] =   0.9122344282513259E+00;
    xtab[36] =   0.9408226338317548E+00;
    xtab[37] =   0.9639719272779138E+00;
    xtab[38] =   0.9815078774502503E+00;
    xtab[39] =   0.9931285991850949E+00;
    xtab[40] =   0.9988590315882777E+00;

    weight[0] =  0.3073583718520532E-02;
    weight[1] =  0.8600269855642942E-02;
    weight[2] =  0.1462616925697125E-01;
    weight[3] =  0.2038837346126652E-01;
    weight[4] =  0.2588213360495116E-01;
    weight[5] =  0.3128730677703280E-01;
    weight[6] =  0.3660016975820080E-01;
    weight[7] =  0.4166887332797369E-01;
    weight[8] =  0.4643482186749767E-01;
    weight[9] = 0.5094457392372869E-01;
    weight[10] = 0.5519510534828599E-01;
    weight[11] = 0.5911140088063957E-01;
    weight[12] = 0.6265323755478117E-01;
    weight[13] = 0.6583459713361842E-01;
    weight[14] = 0.6864867292852162E-01;
    weight[15] = 0.7105442355344407E-01;
    weight[16] = 0.7303069033278667E-01;
    weight[17] = 0.7458287540049919E-01;
    weight[18] = 0.7570449768455667E-01;
    weight[19] = 0.7637786767208074E-01;
    weight[20] = 0.7660071191799966E-01;
    weight[21] = 0.7637786767208074E-01;
    weight[22] = 0.7570449768455667E-01;
    weight[23] = 0.7458287540049919E-01;
    weight[24] = 0.7303069033278667E-01;
    weight[25] = 0.7105442355344407E-01;
    weight[26] = 0.6864867292852162E-01;
    weight[27] = 0.6583459713361842E-01;
    weight[28] = 0.6265323755478117E-01;
    weight[29] = 0.5911140088063957E-01;
    weight[30] = 0.5519510534828599E-01;
    weight[31] = 0.5094457392372869E-01;
    weight[32] = 0.4643482186749767E-01;
    weight[33] = 0.4166887332797369E-01;
    weight[34] = 0.3660016975820080E-01;
    weight[35] = 0.3128730677703280E-01;
    weight[36] = 0.2588213360495116E-01;
    weight[37] = 0.2038837346126652E-01;
    weight[38] = 0.1462616925697125E-01;
    weight[39] = 0.8600269855642942E-02;
    weight[40] = 0.3073583718520532E-02;
  }
  else
  {
    cerr << "\n";
    cerr << "KRONROD_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are 15, 21, 31 or 41.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void laguerre_ek_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_EK_COMPUTE: Laguerre quadrature rule by the Elhay-Kautsky method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *bj;
  int i;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = 1.0;
//
//  Define the Jacobi matrix.
//
  bj = new double[n];

  for ( i = 0; i < n; i++ )
  {
    bj[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( 2 * i + 1 );
  }

  w[0] = sqrt ( zemu );

  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double laguerre_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
//
//  Discussion:
//
//    The integral:
//
//      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double EXACT, the value of the integral.
//
{
  double exact;

  exact = r8_factorial ( expon );

  return exact;
}
//****************************************************************************80

void laguerre_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_SET sets abscissas and weights for Laguerre quadrature.
//
//  Discussion:
//
//    The abscissas are the zeroes of the Laguerre polynomial L(N)(X).
//
//    The integral:
//
//      Integral ( 0 <= X < +oo ) exp ( -X ) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * f ( X(I) )
//
//    The integral:
//
//      Integral ( 0 <= X < +oo ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * exp ( X(I) ) * f ( X(I) )
//
//    Mathematica can numerically estimate the abscissas for the
//    n-th order polynomial to p digits of precision by the command:
//
//      NSolve [ LaguerreL[n,x] == 0, x, p ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20, 31/32/33, 63/64/65, 127/128/129.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[0] =  1.00000000000000000000000000000E+00;

    w[0] =  1.00000000000000000000000000000E+00;
  }
  else if ( n == 2 )
  {
    x[0] = 0.585786437626904951198311275790E+00;
    x[1] = 3.41421356237309504880168872421E+00;

    w[0] = 0.85355339059327376220042218105E+00;
    w[1] = 0.146446609406726237799577818948E+00;
  }
  else if ( n == 3 )
  {
    x[0] = 0.415774556783479083311533873128E+00;
    x[1] = 2.29428036027904171982205036136E+00;
    x[2] = 6.28994508293747919686641576551E+00;

    w[0] = 0.71109300992917301544959019114E+00;
    w[1] = 0.27851773356924084880144488846E+00;
    w[2] = 0.010389256501586135748964920401E+00;
  }
  else if ( n == 4 )
  {
    x[0] = 0.322547689619392311800361459104E+00;
    x[1] = 1.74576110115834657568681671252E+00;
    x[2] = 4.53662029692112798327928538496E+00;
    x[3] = 9.39507091230113312923353644342E+00;

    w[0] = 0.60315410434163360163596602382E+00;
    w[1] = 0.35741869243779968664149201746E+00;
    w[2] = 0.03888790851500538427243816816E+00;
    w[3] = 0.0005392947055613274501037905676E+00;
  }
  else if ( n == 5 )
  {
    x[0] = 0.263560319718140910203061943361E+00;
    x[1] = 1.41340305910651679221840798019E+00;
    x[2] = 3.59642577104072208122318658878E+00;
    x[3] = 7.08581000585883755692212418111E+00;
    x[4] = 12.6408008442757826594332193066E+00;

    w[0] = 0.52175561058280865247586092879E+00;
    w[1] = 0.3986668110831759274541333481E+00;
    w[2] = 0.0759424496817075953876533114E+00;
    w[3] = 0.00361175867992204845446126257E+00;
    w[4] = 0.00002336997238577622789114908455E+00;
  }
  else if ( n == 6 )
  {
    x[0] = 0.222846604179260689464354826787E+00;
    x[1] = 1.18893210167262303074315092194E+00;
    x[2] = 2.99273632605931407769132528451E+00;
    x[3] = 5.77514356910451050183983036943E+00;
    x[4] = 9.83746741838258991771554702994E+00;
    x[5] = 15.9828739806017017825457915674E+00;

    w[0] = 0.45896467394996359356828487771E+00;
    w[1] = 0.4170008307721209941133775662E+00;
    w[2] = 0.1133733820740449757387061851E+00;
    w[3] = 0.01039919745314907489891330285E+00;
    w[4] = 0.000261017202814932059479242860E+00;
    w[5] = 8.98547906429621238825292053E-07;
  }
  else if ( n == 7 )
  {
    x[0] = 0.193043676560362413838247885004E+00;
    x[1] = 1.02666489533919195034519944317E+00;
    x[2] = 2.56787674495074620690778622666E+00;
    x[3] = 4.90035308452648456810171437810E+00;
    x[4] = 8.18215344456286079108182755123E+00;
    x[5] = 12.7341802917978137580126424582E+00;
    x[6] = 19.3957278622625403117125820576E+00;

    w[0] = 0.40931895170127390213043288002E+00;
    w[1] = 0.4218312778617197799292810054E+00;
    w[2] = 0.1471263486575052783953741846E+00;
    w[3] = 0.0206335144687169398657056150E+00;
    w[4] = 0.00107401014328074552213195963E+00;
    w[5] = 0.0000158654643485642012687326223E+00;
    w[6] = 3.17031547899558056227132215E-08;
  }
  else if ( n == 8 )
  {
    x[0] = 0.170279632305100999788861856608E+00;
    x[1] = 0.903701776799379912186020223555E+00;
    x[2] = 2.25108662986613068930711836697E+00;
    x[3] = 4.26670017028765879364942182690E+00;
    x[4] = 7.04590540239346569727932548212E+00;
    x[5] = 10.7585160101809952240599567880E+00;
    x[6] = 15.7406786412780045780287611584E+00;
    x[7] = 22.8631317368892641057005342974E+00;

    w[0] = 0.36918858934163752992058283938E+00;
    w[1] = 0.4187867808143429560769785813E+00;
    w[2] = 0.175794986637171805699659867E+00;
    w[3] = 0.033343492261215651522132535E+00;
    w[4] = 0.0027945362352256725249389241E+00;
    w[5] = 0.00009076508773358213104238501E+00;
    w[6] = 8.4857467162725315448680183E-07;
    w[7] = 1.04800117487151038161508854E-09;
  }
  else if ( n == 9 )
  {
    x[0] = 0.152322227731808247428107073127E+00;
    x[1] = 0.807220022742255847741419210952E+00;
    x[2] = 2.00513515561934712298303324701E+00;
    x[3] = 3.78347397333123299167540609364E+00;
    x[4] = 6.20495677787661260697353521006E+00;
    x[5] = 9.37298525168757620180971073215E+00;
    x[6] = 13.4662369110920935710978818397E+00;
    x[7] = 18.8335977889916966141498992996E+00;
    x[8] = 26.3740718909273767961410072937E+00;

    w[0] = 0.336126421797962519673467717606E+00;
    w[1] = 0.411213980423984387309146942793E+00;
    w[2] = 0.199287525370885580860575607212E+00;
    w[3] = 0.0474605627656515992621163600479E+00;
    w[4] = 0.00559962661079458317700419900556E+00;
    w[5] = 0.000305249767093210566305412824291E+00;
    w[6] = 6.59212302607535239225572284875E-06;
    w[7] = 4.1107693303495484429024104033E-08;
    w[8] = 3.29087403035070757646681380323E-11;
  }
  else if ( n == 10 )
  {
    x[0] = 0.137793470540492430830772505653E+00;
    x[1] = 0.729454549503170498160373121676E+00;
    x[2] = 1.80834290174031604823292007575E+00;
    x[3] = 3.40143369785489951448253222141E+00;
    x[4] = 5.55249614006380363241755848687E+00;
    x[5] = 8.33015274676449670023876719727E+00;
    x[6] = 11.8437858379000655649185389191E+00;
    x[7] = 16.2792578313781020995326539358E+00;
    x[8] = 21.9965858119807619512770901956E+00;
    x[9] = 29.9206970122738915599087933408E+00;

    w[0] = 0.30844111576502014154747083468E+00;
    w[1] = 0.4011199291552735515157803099E+00;
    w[2] = 0.218068287611809421588648523E+00;
    w[3] = 0.062087456098677747392902129E+00;
    w[4] = 0.009501516975181100553839072E+00;
    w[5] = 0.0007530083885875387754559644E+00;
    w[6] = 0.00002825923349599565567422564E+00;
    w[7] = 4.249313984962686372586577E-07;
    w[8] = 1.839564823979630780921535E-09;
    w[9] = 9.911827219609008558377547E-13;
  }
  else if ( n == 11 )
  {
    x[0] = 0.125796442187967522675794577516E+00;
    x[1] = 0.665418255839227841678127839420E+00;
    x[2] = 1.64715054587216930958700321365E+00;
    x[3] = 3.09113814303525495330195934259E+00;
    x[4] = 5.02928440157983321236999508366E+00;
    x[5] = 7.50988786380661681941099714450E+00;
    x[6] = 10.6059509995469677805559216457E+00;
    x[7] = 14.4316137580641855353200450349E+00;
    x[8] = 19.1788574032146786478174853989E+00;
    x[9] = 25.2177093396775611040909447797E+00;
    x[10] = 33.4971928471755372731917259395E+00;

    w[0] = 0.28493321289420060505605102472E+00;
    w[1] = 0.3897208895278493779375535080E+00;
    w[2] = 0.232781831848991333940223796E+00;
    w[3] = 0.076564453546196686400854179E+00;
    w[4] = 0.014393282767350695091863919E+00;
    w[5] = 0.001518880846484873069847776E+00;
    w[6] = 0.0000851312243547192259720424E+00;
    w[7] = 2.29240387957450407857683E-06;
    w[8] = 2.48635370276779587373391E-08;
    w[9] = 7.71262693369132047028153E-11;
    w[10] = 2.883775868323623861597778E-14;
  }
  else if ( n == 12 )
  {
    x[0] = 0.115722117358020675267196428240E+00;
    x[1] = 0.611757484515130665391630053042E+00;
    x[2] = 1.51261026977641878678173792687E+00;
    x[3] = 2.83375133774350722862747177657E+00;
    x[4] = 4.59922763941834848460572922485E+00;
    x[5] = 6.84452545311517734775433041849E+00;
    x[6] = 9.62131684245686704391238234923E+00;
    x[7] = 13.0060549933063477203460524294E+00;
    x[8] = 17.1168551874622557281840528008E+00;
    x[9] = 22.1510903793970056699218950837E+00;
    x[10] = 28.4879672509840003125686072325E+00;
    x[11] = 37.0991210444669203366389142764E+00;

    w[0] = 0.26473137105544319034973889206E+00;
    w[1] = 0.3777592758731379820244905567E+00;
    w[2] = 0.244082011319877564254870818E+00;
    w[3] = 0.09044922221168093072750549E+00;
    w[4] = 0.02010238115463409652266129E+00;
    w[5] = 0.002663973541865315881054158E+00;
    w[6] = 0.000203231592662999392121433E+00;
    w[7] = 8.3650558568197987453363E-06;
    w[8] = 1.66849387654091026116990E-07;
    w[9] = 1.34239103051500414552392E-09;
    w[10] = 3.06160163503502078142408E-12;
    w[11] = 8.148077467426241682473119E-16;
  }
  else if ( n == 13 )
  {
    x[0] = 0.107142388472252310648493376977E+00;
    x[1] = 0.566131899040401853406036347177E+00;
    x[2] = 1.39856433645101971792750259921E+00;
    x[3] = 2.61659710840641129808364008472E+00;
    x[4] = 4.23884592901703327937303389926E+00;
    x[5] = 6.29225627114007378039376523025E+00;
    x[6] = 8.81500194118697804733348868036E+00;
    x[7] = 11.8614035888112425762212021880E+00;
    x[8] = 15.5107620377037527818478532958E+00;
    x[9] = 19.8846356638802283332036594634E+00;
    x[10] = 25.1852638646777580842970297823E+00;
    x[11] = 31.8003863019472683713663283526E+00;
    x[12] = 40.7230086692655795658979667001E+00;

    w[0] = 0.24718870842996262134624918596E+00;
    w[1] = 0.3656888229005219453067175309E+00;
    w[2] = 0.252562420057658502356824289E+00;
    w[3] = 0.10347075802418370511421863E+00;
    w[4] = 0.02643275441556161577815877E+00;
    w[5] = 0.00422039604025475276555209E+00;
    w[6] = 0.000411881770472734774892473E+00;
    w[7] = 0.0000235154739815532386882897E+00;
    w[8] = 7.3173116202490991040105E-07;
    w[9] = 1.10884162570398067979151E-08;
    w[10] = 6.7708266922058988406462E-11;
    w[11] = 1.15997995990507606094507E-13;
    w[12] = 2.245093203892758415991872E-17;
  }
  else if ( n == 14 )
  {
    x[0] = 0.0997475070325975745736829452514E+00;
    x[1] = 0.526857648851902896404583451502E+00;
    x[2] = 1.30062912125149648170842022116E+00;
    x[3] = 2.43080107873084463616999751038E+00;
    x[4] = 3.93210282229321888213134366778E+00;
    x[5] = 5.82553621830170841933899983898E+00;
    x[6] = 8.14024014156514503005978046052E+00;
    x[7] = 10.9164995073660188408130510904E+00;
    x[8] = 14.2108050111612886831059780825E+00;
    x[9] = 18.1048922202180984125546272083E+00;
    x[10] = 22.7233816282696248232280886985E+00;
    x[11] = 28.2729817232482056954158923218E+00;
    x[12] = 35.1494436605924265828643121364E+00;
    x[13] = 44.3660817111174230416312423666E+00;

    w[0] = 0.23181557714486497784077486110E+00;
    w[1] = 0.3537846915975431518023313013E+00;
    w[2] = 0.258734610245428085987320561E+00;
    w[3] = 0.11548289355692321008730499E+00;
    w[4] = 0.03319209215933736003874996E+00;
    w[5] = 0.00619286943700661021678786E+00;
    w[6] = 0.00073989037786738594242589E+00;
    w[7] = 0.000054907194668416983785733E+00;
    w[8] = 2.4095857640853774967578E-06;
    w[9] = 5.801543981676495180886E-08;
    w[10] = 6.819314692484974119616E-10;
    w[11] = 3.2212077518948479398089E-12;
    w[12] = 4.2213524405165873515980E-15;
    w[13] = 6.05237502228918880839871E-19;
  }
  else if ( n == 15 )
  {
    x[0] = 0.0933078120172818047629030383672E+00;
    x[1] = 0.492691740301883908960101791412E+00;
    x[2] = 1.21559541207094946372992716488E+00;
    x[3] = 2.26994952620374320247421741375E+00;
    x[4] = 3.66762272175143727724905959436E+00;
    x[5] = 5.42533662741355316534358132596E+00;
    x[6] = 7.56591622661306786049739555812E+00;
    x[7] = 10.1202285680191127347927394568E+00;
    x[8] = 13.1302824821757235640991204176E+00;
    x[9] = 16.6544077083299578225202408430E+00;
    x[10] = 20.7764788994487667729157175676E+00;
    x[11] = 25.6238942267287801445868285977E+00;
    x[12] = 31.4075191697539385152432196202E+00;
    x[13] = 38.5306833064860094162515167595E+00;
    x[14] = 48.0260855726857943465734308508E+00;

    w[0] = 0.21823488594008688985641323645E+00;
    w[1] = 0.3422101779228833296389489568E+00;
    w[2] = 0.263027577941680097414812275E+00;
    w[3] = 0.12642581810593053584303055E+00;
    w[4] = 0.04020686492100091484158548E+00;
    w[5] = 0.00856387780361183836391576E+00;
    w[6] = 0.00121243614721425207621921E+00;
    w[7] = 0.00011167439234425194199258E+00;
    w[8] = 6.459926762022900924653E-06;
    w[9] = 2.226316907096272630332E-07;
    w[10] = 4.227430384979365007351E-09;
    w[11] = 3.921897267041089290385E-11;
    w[12] = 1.4565152640731264063327E-13;
    w[13] = 1.4830270511133013354616E-16;
    w[14] = 1.60059490621113323104998E-20;
  }
  else if ( n == 16 )
  {
    x[0] = 0.0876494104789278403601980973401E+00;
    x[1] = 0.462696328915080831880838260664E+00;
    x[2] = 1.14105777483122685687794501811E+00;
    x[3] = 2.12928364509838061632615907066E+00;
    x[4] = 3.43708663389320664523510701675E+00;
    x[5] = 5.07801861454976791292305830814E+00;
    x[6] = 7.07033853504823413039598947080E+00;
    x[7] = 9.43831433639193878394724672911E+00;
    x[8] = 12.2142233688661587369391246088E+00;
    x[9] = 15.4415273687816170767647741622E+00;
    x[10] = 19.1801568567531348546631409497E+00;
    x[11] = 23.5159056939919085318231872752E+00;
    x[12] = 28.5787297428821403675206137099E+00;
    x[13] = 34.5833987022866258145276871778E+00;
    x[14] = 41.9404526476883326354722330252E+00;
    x[15] = 51.7011603395433183643426971197E+00;

    w[0] = 0.20615171495780099433427363674E+00;
    w[1] = 0.3310578549508841659929830987E+00;
    w[2] = 0.265795777644214152599502021E+00;
    w[3] = 0.13629693429637753997554751E+00;
    w[4] = 0.0473289286941252189780623E+00;
    w[5] = 0.0112999000803394532312490E+00;
    w[6] = 0.0018490709435263108642918E+00;
    w[7] = 0.00020427191530827846012602E+00;
    w[8] = 0.00001484458687398129877135E+00;
    w[9] = 6.828319330871199564396E-07;
    w[10] = 1.881024841079673213882E-08;
    w[11] = 2.862350242973881619631E-10;
    w[12] = 2.127079033224102967390E-12;
    w[13] = 6.297967002517867787174E-15;
    w[14] = 5.050473700035512820402E-18;
    w[15] = 4.1614623703728551904265E-22;
  }
  else if ( n == 17 )
  {
    x[0] = 0.0826382147089476690543986151980E+00;
    x[1] = 0.436150323558710436375959029847E+00;
    x[2] = 1.07517657751142857732980316755E+00;
    x[3] = 2.00519353164923224070293371933E+00;
    x[4] = 3.23425612404744376157380120696E+00;
    x[5] = 4.77351351370019726480932076262E+00;
    x[6] = 6.63782920536495266541643929703E+00;
    x[7] = 8.84668551116980005369470571184E+00;
    x[8] = 11.4255293193733525869726151469E+00;
    x[9] = 14.4078230374813180021982874959E+00;
    x[10] = 17.8382847307011409290658752412E+00;
    x[11] = 21.7782682577222653261749080522E+00;
    x[12] = 26.3153178112487997766149598369E+00;
    x[13] = 31.5817716804567331343908517497E+00;
    x[14] = 37.7960938374771007286092846663E+00;
    x[15] = 45.3757165339889661829258363215E+00;
    x[16] = 55.3897517898396106640900199790E+00;

    w[0] = 0.19533220525177083214592729770E+00;
    w[1] = 0.3203753572745402813366256320E+00;
    w[2] = 0.267329726357171097238809604E+00;
    w[3] = 0.14512985435875862540742645E+00;
    w[4] = 0.0544369432453384577793806E+00;
    w[5] = 0.0143572977660618672917767E+00;
    w[6] = 0.0026628247355727725684324E+00;
    w[7] = 0.0003436797271562999206118E+00;
    w[8] = 0.00003027551783782870109437E+00;
    w[9] = 1.768515053231676895381E-06;
    w[10] = 6.57627288681043332199E-08;
    w[11] = 1.469730932159546790344E-09;
    w[12] = 1.81691036255544979555E-11;
    w[13] = 1.095401388928687402976E-13;
    w[14] = 2.617373882223370421551E-16;
    w[15] = 1.6729356931461546908502E-19;
    w[16] = 1.06562631627404278815253E-23;
  }
  else if ( n == 18 )
  {
    x[0] = 0.0781691666697054712986747615334E+00;
    x[1] = 0.412490085259129291039101536536E+00;
    x[2] = 1.01652017962353968919093686187E+00;
    x[3] = 1.89488850996976091426727831954E+00;
    x[4] = 3.05435311320265975115241130719E+00;
    x[5] = 4.50420553888989282633795571455E+00;
    x[6] = 6.25672507394911145274209116326E+00;
    x[7] = 8.32782515660563002170470261564E+00;
    x[8] = 10.7379900477576093352179033397E+00;
    x[9] = 13.5136562075550898190863812108E+00;
    x[10] = 16.6893062819301059378183984163E+00;
    x[11] = 20.3107676262677428561313764553E+00;
    x[12] = 24.4406813592837027656442257980E+00;
    x[13] = 29.1682086625796161312980677805E+00;
    x[14] = 34.6279270656601721454012429438E+00;
    x[15] = 41.0418167728087581392948614284E+00;
    x[16] = 48.8339227160865227486586093290E+00;
    x[17] = 59.0905464359012507037157810181E+00;

    w[0] = 0.18558860314691880562333775228E+00;
    w[1] = 0.3101817663702252936495975957E+00;
    w[2] = 0.267866567148536354820854395E+00;
    w[3] = 0.15297974746807490655384308E+00;
    w[4] = 0.0614349178609616527076780E+00;
    w[5] = 0.0176872130807729312772600E+00;
    w[6] = 0.0036601797677599177980266E+00;
    w[7] = 0.0005406227870077353231284E+00;
    w[8] = 0.0000561696505121423113818E+00;
    w[9] = 4.01530788370115755859E-06;
    w[10] = 1.91466985667567497969E-07;
    w[11] = 5.8360952686315941292E-09;
    w[12] = 1.07171126695539012773E-10;
    w[13] = 1.08909871388883385562E-12;
    w[14] = 5.38666474837830887608E-15;
    w[15] = 1.049865978035703408779E-17;
    w[16] = 5.405398451631053643566E-21;
    w[17] = 2.6916532692010286270838E-25;
  }
  else if ( n == 19 )
  {
    x[0] = 0.0741587837572050877131369916024E+00;
    x[1] = 0.391268613319994607337648350299E+00;
    x[2] = 0.963957343997958058624878377130E+00;
    x[3] = 1.79617558206832812557725825252E+00;
    x[4] = 2.89365138187378399116494713237E+00;
    x[5] = 4.26421553962776647436040018167E+00;
    x[6] = 5.91814156164404855815360191408E+00;
    x[7] = 7.86861891533473373105668358176E+00;
    x[8] = 10.1324237168152659251627415800E+00;
    x[9] = 12.7308814638423980045092979656E+00;
    x[10] = 15.6912783398358885454136069861E+00;
    x[11] = 19.0489932098235501532136429732E+00;
    x[12] = 22.8508497608294829323930586693E+00;
    x[13] = 27.1606693274114488789963947149E+00;
    x[14] = 32.0691222518622423224362865906E+00;
    x[15] = 37.7129058012196494770647508283E+00;
    x[16] = 44.3173627958314961196067736013E+00;
    x[17] = 52.3129024574043831658644222420E+00;
    x[18] = 62.8024231535003758413504690673E+00;

    w[0] = 0.17676847491591250225103547981E+00;
    w[1] = 0.3004781436072543794821568077E+00;
    w[2] = 0.267599547038175030772695441E+00;
    w[3] = 0.15991337213558021678551215E+00;
    w[4] = 0.0682493799761491134552355E+00;
    w[5] = 0.0212393076065443249244062E+00;
    w[6] = 0.0048416273511483959672501E+00;
    w[7] = 0.0008049127473813667665946E+00;
    w[8] = 0.0000965247209315350170843E+00;
    w[9] = 8.20730525805103054409E-06;
    w[10] = 4.8305667247307725394E-07;
    w[11] = 1.90499136112328569994E-08;
    w[12] = 4.8166846309280615577E-10;
    w[13] = 7.3482588395511443768E-12;
    w[14] = 6.2022753875726163989E-14;
    w[15] = 2.54143084301542272372E-16;
    w[16] = 4.07886129682571235007E-19;
    w[17] = 1.707750187593837061004E-22;
    w[18] = 6.715064649908189959990E-27;
  }
  else if ( n == 20 )
  {
    x[0] = 0.0705398896919887533666890045842E+00;
    x[1] = 0.372126818001611443794241388761E+00;
    x[2] = 0.916582102483273564667716277074E+00;
    x[3] = 1.70730653102834388068768966741E+00;
    x[4] = 2.74919925530943212964503046049E+00;
    x[5] = 4.04892531385088692237495336913E+00;
    x[6] = 5.61517497086161651410453988565E+00;
    x[7] = 7.45901745367106330976886021837E+00;
    x[8] = 9.59439286958109677247367273428E+00;
    x[9] = 12.0388025469643163096234092989E+00;
    x[10] = 14.8142934426307399785126797100E+00;
    x[11] = 17.9488955205193760173657909926E+00;
    x[12] = 21.4787882402850109757351703696E+00;
    x[13] = 25.4517027931869055035186774846E+00;
    x[14] = 29.9325546317006120067136561352E+00;
    x[15] = 35.0134342404790000062849359067E+00;
    x[16] = 40.8330570567285710620295677078E+00;
    x[17] = 47.6199940473465021399416271529E+00;
    x[18] = 55.8107957500638988907507734445E+00;
    x[19] = 66.5244165256157538186403187915E+00;

    w[0] = 0.168746801851113862149223899689E+00;
    w[1] = 0.291254362006068281716795323812E+00;
    w[2] = 0.266686102867001288549520868998E+00;
    w[3] = 0.166002453269506840031469127816E+00;
    w[4] = 0.0748260646687923705400624639615E+00;
    w[5] = 0.0249644173092832210728227383234E+00;
    w[6] = 0.00620255084457223684744754785395E+00;
    w[7] = 0.00114496238647690824203955356969E+00;
    w[8] = 0.000155741773027811974779809513214E+00;
    w[9] = 0.0000154014408652249156893806714048E+00;
    w[10] = 1.08648636651798235147970004439E-06;
    w[11] = 5.33012090955671475092780244305E-08;
    w[12] = 1.7579811790505820035778763784E-09;
    w[13] = 3.72550240251232087262924585338E-11;
    w[14] = 4.76752925157819052449488071613E-13;
    w[15] = 3.37284424336243841236506064991E-15;
    w[16] = 1.15501433950039883096396247181E-17;
    w[17] = 1.53952214058234355346383319667E-20;
    w[18] = 5.28644272556915782880273587683E-24;
    w[19] = 1.65645661249902329590781908529E-28;
  }
  else if ( n == 31 )
  {
    x[0] = 0.0459019476211082907434960802752E+00;
    x[1] = 0.241980163824772048904089741517E+00;
    x[2] = 0.595253894222350737073301650054E+00;
    x[3] = 1.10668949953299871621113087898E+00;
    x[4] = 1.77759569287477272115937274827E+00;
    x[5] = 2.60970341525668065038933759253E+00;
    x[6] = 3.60519680234004426988058175542E+00;
    x[7] = 4.76674708447176113136291272711E+00;
    x[8] = 6.09755456718174092699254293285E+00;
    x[9] = 7.60140094923313742293601069429E+00;
    x[10] = 9.28271431347088941825366952977E+00;
    x[11] = 11.1466497556192913589938156296E+00;
    x[12] = 13.1991895762449985224649250286E+00;
    x[13] = 15.4472683155493100758093258918E+00;
    x[14] = 17.8989298266447576467257938178E+00;
    x[15] = 20.5635263367158221707430489688E+00;
    x[16] = 23.4519734820118585910502555759E+00;
    x[17] = 26.5770813521182604599758769865E+00;
    x[18] = 29.9539908723464455069519178400E+00;
    x[19] = 33.6007595329022027354103138858E+00;
    x[20] = 37.5391644073304408828879025580E+00;
    x[21] = 41.7958308701822199813479458533E+00;
    x[22] = 46.4038668064111231360292276044E+00;
    x[23] = 51.4053144767977551618614610884E+00;
    x[24] = 56.8549928687158436205119220557E+00;
    x[25] = 62.8268559087863214536775233048E+00;
    x[26] = 69.4252771910803456233222516564E+00;
    x[27] = 76.8070477638627328376099722855E+00;
    x[28] = 85.2303586075456691693870656070E+00;
    x[29] = 95.1889398915256299813086068540E+00;
    x[30] = 107.952243827578714750024401177E+00;

    w[0] = 0.112527895503725838208477280828E+00;
    w[1] = 0.21552760818089123795222505285E+00;
    w[2] = 0.238308251645696547319057880892E+00;
    w[3] = 0.195388309297902292499153033907E+00;
    w[4] = 0.126982832893061901436352729046E+00;
    w[5] = 0.0671861689238993006709294419935E+00;
    w[6] = 0.029303224993879487404888669312E+00;
    w[7] = 0.0105975699152957360895293803144E+00;
    w[8] = 0.0031851272582386980320974842433E+00;
    w[9] = 0.000795495483079403829220921490125E+00;
    w[10] = 0.000164800521266366873178629671164E+00;
    w[11] = 0.000028229237864310816393860971469E+00;
    w[12] = 3.98029025510085803871161749001E-06;
    w[13] = 4.59318398418010616737296945103E-07;
    w[14] = 4.30755451877311009301314574659E-08;
    w[15] = 3.25512499382715708551757492579E-09;
    w[16] = 1.96202466754105949962471515931E-10;
    w[17] = 9.31904990866175871295347164313E-12;
    w[18] = 3.43775418194116205203125978983E-13;
    w[19] = 9.67952471304467169974050357762E-15;
    w[20] = 2.03680661101152473980106242193E-16;
    w[21] = 3.12126872807135268317653586326E-18;
    w[22] = 3.37295817041610524533956783083E-20;
    w[23] = 2.46727963866166960110383632425E-22;
    w[24] = 1.15822019045256436348345645766E-24;
    w[25] = 3.2472922591425422434798022809E-27;
    w[26] = 4.91430173080574327408200762597E-30;
    w[27] = 3.45000711048083941322231359538E-33;
    w[28] = 8.76637101171620414729327607329E-37;
    w[29] = 5.03636439211614904112971723166E-41;
    w[30] = 1.99099845825314564824395490803E-46;
  }
  else if ( n == 32 )
  {
    x[0] = 0.0444893658332670184188501945244E+00;
    x[1] = 0.234526109519618537452909561302E+00;
    x[2] = 0.576884629301886426491552569378E+00;
    x[3] = 1.07244875381781763304091397718E+00;
    x[4] = 1.72240877644464544113093292797E+00;
    x[5] = 2.52833670642579488112419990556E+00;
    x[6] = 3.49221327302199448960880339977E+00;
    x[7] = 4.61645676974976738776205229617E+00;
    x[8] = 5.90395850417424394656152149158E+00;
    x[9] = 7.35812673318624111322198973719E+00;
    x[10] = 8.98294092421259610337824752677E+00;
    x[11] = 10.7830186325399720675019491381E+00;
    x[12] = 12.7636979867427251149690330822E+00;
    x[13] = 14.9311397555225573197969646873E+00;
    x[14] = 17.2924543367153147892357183836E+00;
    x[15] = 19.8558609403360547397899445841E+00;
    x[16] = 22.6308890131967744886757793394E+00;
    x[17] = 25.6286360224592477674761768768E+00;
    x[18] = 28.8621018163234747443426407115E+00;
    x[19] = 32.3466291539647370032321654237E+00;
    x[20] = 36.1004948057519738040171189479E+00;
    x[21] = 40.1457197715394415362093439289E+00;
    x[22] = 44.5092079957549379759066043775E+00;
    x[23] = 49.2243949873086391767222218066E+00;
    x[24] = 54.3337213333969073328671815512E+00;
    x[25] = 59.8925091621340181961304753247E+00;
    x[26] = 65.9753772879350527965630761193E+00;
    x[27] = 72.6876280906627086386753490878E+00;
    x[28] = 80.1874469779135230674916385687E+00;
    x[29] = 88.7353404178923986893554495243E+00;
    x[30] = 98.8295428682839725591844784095E+00;
    x[31] = 111.751398097937695213664716539E+00;

    w[0] = 0.109218341952384971136131337943E+00;
    w[1] = 0.210443107938813232936062071392E+00;
    w[2] = 0.235213229669848005394941106676E+00;
    w[3] = 0.195903335972881043413247901182E+00;
    w[4] = 0.129983786286071760607216822326E+00;
    w[5] = 0.0705786238657174415601643320433E+00;
    w[6] = 0.0317609125091750703058255211629E+00;
    w[7] = 0.0119182148348385570565446505756E+00;
    w[8] = 0.00373881629461152478966122847796E+00;
    w[9] = 0.000980803306614955132230630611308E+00;
    w[10] = 0.000214864918801364188023199483686E+00;
    w[11] = 0.0000392034196798794720432695682782E+00;
    w[12] = 5.93454161286863287835582893773E-06;
    w[13] = 7.4164045786675522190708220213E-07;
    w[14] = 7.60456787912078148111926545943E-08;
    w[15] = 6.35060222662580674242777108552E-09;
    w[16] = 4.28138297104092887881360582098E-10;
    w[17] = 2.30589949189133607927336809618E-11;
    w[18] = 9.79937928872709406333455225995E-13;
    w[19] = 3.23780165772926646231042646142E-14;
    w[20] = 8.17182344342071943320186059177E-16;
    w[21] = 1.54213383339382337217855949129E-17;
    w[22] = 2.11979229016361861204093474373E-19;
    w[23] = 2.05442967378804542665570987602E-21;
    w[24] = 1.34698258663739515580519340478E-23;
    w[25] = 5.66129413039735937112634432382E-26;
    w[26] = 1.41856054546303690595142933892E-28;
    w[27] = 1.91337549445422430937127829683E-31;
    w[28] = 1.19224876009822235654164532831E-34;
    w[29] = 2.67151121924013698599867893958E-38;
    w[30] = 1.33861694210625628271905701423E-42;
    w[31] = 4.51053619389897423222342830132E-48;
  }
  else if ( n == 33 )
  {
    x[0] = 0.0431611356173268921917334738206E+00;
    x[1] = 0.227517802803371123850290226913E+00;
    x[2] = 0.559616655851539887586282303916E+00;
    x[3] = 1.04026850775100205382209621927E+00;
    x[4] = 1.67055919607571519092562973257E+00;
    x[5] = 2.45192079589763054651073898192E+00;
    x[6] = 3.38615533758800483230187851832E+00;
    x[7] = 4.47545949839977145702059137905E+00;
    x[8] = 5.72245472027210352266790817933E+00;
    x[9] = 7.13022434440010801631414039534E+00;
    x[10] = 8.70235923062140624893696399459E+00;
    x[11] = 10.4430136502059824268455293839E+00;
    x[12] = 12.3569737593502859624441255236E+00;
    x[13] = 14.4497416815855402377145121178E+00;
    x[14] = 16.7276392186383223532615229942E+00;
    x[15] = 19.1979365872124466372283088222E+00;
    x[16] = 21.8690135249281898713512287043E+00;
    x[17] = 24.7505629061577956433730931987E+00;
    x[18] = 27.8538511114133567797747375537E+00;
    x[19] = 31.1920555455751298677734295989E+00;
    x[20] = 34.7807091535383377002292521853E+00;
    x[21] = 38.6382967177740302250360622751E+00;
    x[22] = 42.7870720782534794879639219927E+00;
    x[23] = 47.2542066029932658172690829767E+00;
    x[24] = 52.0734519015142202671640200482E+00;
    x[25] = 57.2876345410929400754514841078E+00;
    x[26] = 62.9525659469066302071906336861E+00;
    x[27] = 69.1435133801098924457366348147E+00;
    x[28] = 75.9666870142470623437939790250E+00;
    x[29] = 83.5816372232708807614192336050E+00;
    x[30] = 92.2511394441351012341481184391E+00;
    x[31] = 102.477844336823322575825984750E+00;
    x[32] = 115.554756448995807306876850793E+00;

    w[0] = 0.106097745553686759448980241986E+00;
    w[1] = 0.205582983661932603502389046086E+00;
    w[2] = 0.232126523496060850848680143719E+00;
    w[3] = 0.196207372769141916829837191484E+00;
    w[4] = 0.132744856705171098099698375677E+00;
    w[5] = 0.0738518038877138714048058524075E+00;
    w[6] = 0.0342232334108640270351258175641E+00;
    w[7] = 0.0132939751808086665861981468532E+00;
    w[8] = 0.00434094309504623645941229723703E+00;
    w[9] = 0.0011922509906686840510776728263E+00;
    w[10] = 0.000275158225582396584420253012954E+00;
    w[11] = 0.0000532433409922782412444424323192E+00;
    w[12] = 8.60957132646331618369236329569E-06;
    w[13] = 1.15837796102469195604266695945E-06;
    w[14] = 1.28985114856525884538052927779E-07;
    w[15] = 1.18096786980325241031580363325E-08;
    w[16] = 8.82276640967020246770192412272E-10;
    w[17] = 5.32961213410302701363055555216E-11;
    w[18] = 2.57550403748317439431393144398E-12;
    w[19] = 9.8314133225207825980561863437E-14;
    w[20] = 2.92041495556546845392792035892E-15;
    w[21] = 6.63077156752381637730149511056E-17;
    w[22] = 1.12609863704995018019882580368E-18;
    w[23] = 1.39311657122392009414616980902E-20;
    w[24] = 1.21481009891544297673141523063E-22;
    w[25] = 7.16158181142099840535743381278E-25;
    w[26] = 2.70320712488116872172720734689E-27;
    w[27] = 6.07192361286922243586897316027E-30;
    w[28] = 7.32134211132579407517616095945E-33;
    w[29] = 4.06131706145569795511645700604E-36;
    w[30] = 8.04952284545203726871355981553E-40;
    w[31] = 3.52902990360469937522008596417E-44;
    w[32] = 1.01716656299412569799194166119E-49;
  }
  else if ( n == 63 )
  {
    x[0] = 0.0227688937325761537859943302486E+00;
    x[1] = 0.119983252427278247157714164264E+00;
    x[2] = 0.294941854447701495774277385174E+00;
    x[3] = 0.547790878962377253638650737759E+00;
    x[4] = 0.878690611799319016738955670523E+00;
    x[5] = 1.28784643359197063023092077886E+00;
    x[6] = 1.77551238153885537639794632687E+00;
    x[7] = 2.34199255670859892560556283377E+00;
    x[8] = 2.98764232232464739399767310536E+00;
    x[9] = 3.71286959920180003462996374134E+00;
    x[10] = 4.51813633495035843911055685616E+00;
    x[11] = 5.40396017818259462869025997827E+00;
    x[12] = 6.37091637878653302203922508918E+00;
    x[13] = 7.41963993393117111548884931990E+00;
    x[14] = 8.55082800084033283125890487222E+00;
    x[15] = 9.76524259992453668070045929780E+00;
    x[16] = 11.0637136351406617362205504106E+00;
    x[17] = 12.4471422623564927497986875693E+00;
    x[18] = 13.9165046410578185629129670082E+00;
    x[19] = 15.4728561100362964247771436078E+00;
    x[20] = 17.1173358338635887531169003039E+00;
    x[21] = 18.8511719741548568508734837875E+00;
    x[22] = 20.6756874480565156603772656674E+00;
    x[23] = 22.5923063463115283812922777600E+00;
    x[24] = 24.6025610949726388837006427600E+00;
    x[25] = 26.7081004587373439697790879988E+00;
    x[26] = 28.9106985004513826401777181032E+00;
    x[27] = 31.2122646311759128854777738208E+00;
    x[28] = 33.6148549091011548365988428883E+00;
    x[29] = 36.1206847744848230563063287408E+00;
    x[30] = 38.7321434429335821456260416077E+00;
    x[31] = 41.4518102223187411911147261814E+00;
    x[32] = 44.2824730714792338393588571346E+00;
    x[33] = 47.2271497842956868989350952315E+00;
    x[34] = 50.2891122642406957617490218394E+00;
    x[35] = 53.4719144567886528083482806195E+00;
    x[36] = 56.7794246363420622130997810571E+00;
    x[37] = 60.2158629090198628864175501144E+00;
    x[38] = 63.7858450042359746317011396018E+00;
    x[39] = 67.4944337022938858303743256950E+00;
    x[40] = 71.3471996042952662866548033761E+00;
    x[41] = 75.3502934256532342542905047443E+00;
    x[42] = 79.5105326299863091495553913548E+00;
    x[43] = 83.8355060808722578433398176585E+00;
    x[44] = 88.3337015703543690861127663265E+00;
    x[45] = 93.0146627285585474053033990371E+00;
    x[46] = 97.8891841475781400433867276771E+00;
    x[47] = 102.969556907413816507839527468E+00;
    x[48] = 108.269881619615953922263509672E+00;
    x[49] = 113.806473502874627389344859559E+00;
    x[50] = 119.598395388304586669624529633E+00;
    x[51] = 125.668172558561194312911963033E+00;
    x[52] = 132.042772720911657465855905830E+00;
    x[53] = 138.754984181037890781675905675E+00;
    x[54] = 145.845413183135403582839942484E+00;
    x[55] = 153.365484594978636237108159627E+00;
    x[56] = 161.382151948137612435621726696E+00;
    x[57] = 169.985706006658394387951753012E+00;
    x[58] = 179.303662474015809102518278585E+00;
    x[59] = 189.527895965324754736687213330E+00;
    x[60] = 200.975211599246567416286718410E+00;
    x[61] = 214.253685366387886426980562964E+00;
    x[62] = 230.934657470897039712465629851E+00;

    w[0] = 0.0571186332138689798115872833905E+00;
    w[1] = 0.120674760906403952833199320363E+00;
    w[2] = 0.159250010965818737238705610965E+00;
    w[3] = 0.168751783275607992345961929636E+00;
    w[4] = 0.153666419776689566961937113101E+00;
    w[5] = 0.123687706147164816410866522619E+00;
    w[6] = 0.0892750988548486715452791500574E+00;
    w[7] = 0.0582584854461059449575718257252E+00;
    w[8] = 0.0345466575459925808747170858125E+00;
    w[9] = 0.0186756859857146567982865525912E+00;
    w[10] = 0.00922334490440935365284900752416E+00;
    w[11] = 0.00416712506848395927625826634702E+00;
    w[12] = 0.0017238120299900582715386728542E+00;
    w[13] = 0.00065320845029716311169340559359E+00;
    w[14] = 0.000226776446709095869524051732075E+00;
    w[15] = 0.0000721276741548106684107502702349E+00;
    w[16] = 0.0000210112611804664845988115368512E+00;
    w[17] = 5.60355008933572127491815360713E-06;
    w[18] = 1.3673642785604888017836641282E-06;
    w[19] = 3.050726393019581724073609719E-07;
    w[20] = 6.21800618393097635599817754E-08;
    w[21] = 1.1566529551931711260022449E-08;
    w[22] = 1.9614588267565478081534782E-09;
    w[23] = 3.028617119570941124433476E-10;
    w[24] = 4.252134453940068676901296E-11;
    w[25] = 5.42022205780738193346988E-12;
    w[26] = 6.2627306838597672554167E-13;
    w[27] = 6.5474443156573322992307E-14;
    w[28] = 6.18155758087291818463E-15;
    w[29] = 5.259272136350738140426E-16;
    w[30] = 4.023092009264648401539E-17;
    w[31] = 2.7600740511819536505E-18;
    w[32] = 1.69369467569682960533E-19;
    w[33] = 9.2689146872177087315E-21;
    w[34] = 4.509373906036563294E-22;
    w[35] = 1.9435162876132376574E-23;
    w[36] = 7.392627089516920704E-25;
    w[37] = 2.471436415443463262E-26;
    w[38] = 7.228864944674159766E-28;
    w[39] = 1.840761729261403936E-29;
    w[40] = 4.058349856684196011E-31;
    w[41] = 7.70004964164383681E-33;
    w[42] = 1.248850576499933433E-34;
    w[43] = 1.71850002267670107E-36;
    w[44] = 1.989637263667239694E-38;
    w[45] = 1.919967137880405827E-40;
    w[46] = 1.527858828552216692E-42;
    w[47] = 9.9054752688842143E-45;
    w[48] = 5.159752367302921188E-47;
    w[49] = 2.124984666408411125E-49;
    w[50] = 6.790385276685291059E-52;
    w[51] = 1.6466654148296177468E-54;
    w[52] = 2.9509065402691055027E-57;
    w[53] = 3.78384206475710519849E-60;
    w[54] = 3.33581300685424318782E-63;
    w[55] = 1.922346102227388098136E-66;
    w[56] = 6.781269696108301687278E-70;
    w[57] = 1.34047528024406046076205E-73;
    w[58] = 1.3109745101805029757648E-77;
    w[59] = 5.262486388140178738869458E-82;
    w[60] = 6.3780013856587414257760666E-87;
    w[61] = 1.299707894237292456634747392E-92;
    w[62] = 1.00085114969687540634437401684E-99;
  }
  else if ( n == 64 )
  {
    x[0] = 0.0224158741467052800228118848190E+00;
    x[1] = 0.118122512096770479797466436710E+00;
    x[2] = 0.290365744018036483999130066385E+00;
    x[3] = 0.539286221227979039318144947812E+00;
    x[4] = 0.865037004648113944619955074710E+00;
    x[5] = 1.26781404077524139811570887769E+00;
    x[6] = 1.74785962605943625282996395129E+00;
    x[7] = 2.30546373930750871854807054389E+00;
    x[8] = 2.94096515672525184067946815211E+00;
    x[9] = 3.65475265020729052703539791209E+00;
    x[10] = 4.44726634331309435674255098016E+00;
    x[11] = 5.31899925449639034352210985919E+00;
    x[12] = 6.27049904692365391291106464633E+00;
    x[13] = 7.30237000258739574722349840952E+00;
    x[14] = 8.41527523948302419449521859120E+00;
    x[15] = 9.60993919279610803576288204955E+00;
    x[16] = 10.8871503838863721425945504202E+00;
    x[17] = 12.2477645042443016181623692907E+00;
    x[18] = 13.6927078455475051527299325746E+00;
    x[19] = 15.2229811115247288480082687834E+00;
    x[20] = 16.8396636526487372105288380392E+00;
    x[21] = 18.5439181708591905236196259711E+00;
    x[22] = 20.3369959487302355011498158064E+00;
    x[23] = 22.2202426659508765399221543471E+00;
    x[24] = 24.1951048759332539898864438802E+00;
    x[25] = 26.2631372271184857851260239548E+00;
    x[26] = 28.4260105275010272994997715268E+00;
    x[27] = 30.6855207675259717710485823984E+00;
    x[28] = 33.0435992364378291255202106805E+00;
    x[29] = 35.5023238911412095869787785351E+00;
    x[30] = 38.0639321656464682603573179150E+00;
    x[31] = 40.7308354444586263657318695132E+00;
    x[32] = 43.5056354664215298527031849317E+00;
    x[33] = 46.3911429786161920736053999424E+00;
    x[34] = 49.3903990256246866792358008227E+00;
    x[35] = 52.5066993413463016501792769805E+00;
    x[36] = 55.7436224132783804633357112912E+00;
    x[37] = 59.1050619190171066088487420918E+00;
    x[38] = 62.5952644001513955960550179012E+00;
    x[39] = 66.2188732512475643822137626710E+00;
    x[40] = 69.9809803771468292285346579722E+00;
    x[41] = 73.8871872324829632109574031135E+00;
    x[42] = 77.9436774344631203136879758706E+00;
    x[43] = 82.1573037783193042951958683422E+00;
    x[44] = 86.5356933494565182102162783753E+00;
    x[45] = 91.0873756131330901456367153493E+00;
    x[46] = 95.8219400155207320947672154365E+00;
    x[47] = 100.750231969513979629259261451E+00;
    x[48] = 105.884599468799949356360427851E+00;
    x[49] = 111.239207524439582063484736638E+00;
    x[50] = 116.830445051306498463386669077E+00;
    x[51] = 122.677460268538576577419690565E+00;
    x[52] = 128.802878769237672512753623054E+00;
    x[53] = 135.233787949525827833980498879E+00;
    x[54] = 142.003121489931519025140038291E+00;
    x[55] = 149.151665900049388587293462932E+00;
    x[56] = 156.731075132671161233616960814E+00;
    x[57] = 164.808602655150522993190109025E+00;
    x[58] = 173.474946836424274522152844867E+00;
    x[59] = 182.858204691431463646342794510E+00;
    x[60] = 193.151136037072911479385527417E+00;
    x[61] = 204.672028485059455949064433343E+00;
    x[62] = 218.031851935328516332452384448E+00;
    x[63] = 234.809579171326164713055529725E+00;

    w[0] = 0.0562528423390298457410218545063E+00;
    w[1] = 0.119023987312426027814903505889E+00;
    w[2] = 0.157496403862144523820196434706E+00;
    w[3] = 0.167547050415773947880904411659E+00;
    w[4] = 0.153352855779236618085454792564E+00;
    w[5] = 0.124221053609329744512613782193E+00;
    w[6] = 0.0903423009864850577389741092016E+00;
    w[7] = 0.0594777557683550242122469974397E+00;
    w[8] = 0.0356275189040360718541657353369E+00;
    w[9] = 0.0194804104311664060433373802715E+00;
    w[10] = 0.00974359489938200224010796027138E+00;
    w[11] = 0.00446431036416627529236482419054E+00;
    w[12] = 0.00187535958132311482675012822252E+00;
    w[13] = 0.000722646981575005122719108706619E+00;
    w[14] = 0.00025548753283349670971444840218E+00;
    w[15] = 0.0000828714353439694217906322403988E+00;
    w[16] = 0.0000246568639678855874597337022172E+00;
    w[17] = 6.7267138788296685276125455363E-06;
    w[18] = 1.6817853699640888978212010865E-06;
    w[19] = 3.850812981546684414827886111E-07;
    w[20] = 8.06872804099049979041511092E-08;
    w[21] = 1.54572370675768882800370393E-08;
    w[22] = 2.7044801476174814099886074E-09;
    w[23] = 4.316775475427200912314164E-10;
    w[24] = 6.27775254176145220165296E-11;
    w[25] = 8.30631737628895806387893E-12;
    w[26] = 9.9840317872201640558973E-13;
    w[27] = 1.0883538871166626853261E-13;
    w[28] = 1.0740174034415901864828E-14;
    w[29] = 9.57573723157444210559E-16;
    w[30] = 7.69702802364858609886E-17;
    w[31] = 5.56488113745402536653E-18;
    w[32] = 3.6097564090104464983E-19;
    w[33] = 2.0950953695489462348E-20;
    w[34] = 1.0847933010975493612E-21;
    w[35] = 4.994699486363804116E-23;
    w[36] = 2.037836974598822311E-24;
    w[37] = 7.339537564278837039E-26;
    w[38] = 2.323783082198694261E-27;
    w[39] = 6.43823470690876242E-29;
    w[40] = 1.553121095788275271E-30;
    w[41] = 3.24425009201953731E-32;
    w[42] = 5.8323862678362015E-34;
    w[43] = 8.96325483310285406E-36;
    w[44] = 1.168703989550736241E-37;
    w[45] = 1.282055984359980381E-39;
    w[46] = 1.172094937405002292E-41;
    w[47] = 8.83533967232860498E-44;
    w[48] = 5.42495559030618659E-46;
    w[49] = 2.675542666678893829E-48;
    w[50] = 1.042917031411367078E-50;
    w[51] = 3.152902351957772624E-53;
    w[52] = 7.22954191064752234E-56;
    w[53] = 1.2242353012300822645E-58;
    w[54] = 1.48216850490191041178E-61;
    w[55] = 1.23251934881451880806E-64;
    w[56] = 6.69149900457126952681E-68;
    w[57] = 2.220465941850448995507E-71;
    w[58] = 4.1209460947388762499791E-75;
    w[59] = 3.77439906189648917041762E-79;
    w[60] = 1.414115052917619417463147E-83;
    w[61] = 1.5918330640413679178610318E-88;
    w[62] = 2.98948434886063430774131269E-94;
    w[63] = 2.0890635084369527708281542544E-101;
  }
  else if ( n == 65 )
  {
    x[0] = 0.0220736343882500875264595737093E+00;
    x[1] = 0.116318612213376151729622328698E+00;
    x[2] = 0.285929513070813951834551909860E+00;
    x[3] = 0.531041775784488438911664842933E+00;
    x[4] = 0.851801670809046586655299989571E+00;
    x[5] = 1.24839628201831707727195703687E+00;
    x[6] = 1.72105687981755734816623923359E+00;
    x[7] = 2.27006018491690256109784837001E+00;
    x[8] = 2.89572940756799471192249006458E+00;
    x[9] = 3.59843535756478540694788870796E+00;
    x[10] = 4.37859769744118464804519621154E+00;
    x[11] = 5.23668636694320050319815958793E+00;
    x[12] = 6.17322319658730631921440545744E+00;
    x[13] = 7.18878372629467202131889616235E+00;
    x[14] = 8.28399924588871831467573911308E+00;
    x[15] = 9.45955907610065448076092041454E+00;
    x[16] = 10.7162131111897048393914349801E+00;
    x[17] = 12.0547746472322707150580323735E+00;
    x[18] = 13.4761235235671154056092903386E+00;
    x[19] = 14.9812096088541520753799758172E+00;
    x[20] = 16.5710566677964491577962746244E+00;
    x[21] = 18.2467666498987672118974219097E+00;
    x[22] = 20.0095244478287060731090711532E+00;
    x[23] = 21.8606031801774165914390699117E+00;
    x[24] = 23.8013700618933285635853631176E+00;
    x[25] = 25.8332929356397220962712840077E+00;
    x[26] = 27.9579475491204073878989540109E+00;
    x[27] = 30.1770256774182582795717815112E+00;
    x[28] = 32.4923442060863420003362903722E+00;
    x[29] = 34.9058553107319841822666865568E+00;
    x[30] = 37.4196578929107315489849223258E+00;
    x[31] = 40.0360104612770088389674875188E+00;
    x[32] = 42.7573456823686148915804841094E+00;
    x[33] = 45.5862868687358177811114530444E+00;
    x[34] = 48.5256667254367437021219033875E+00;
    x[35] = 51.5785487419130140104165726852E+00;
    x[36] = 54.7482516984863756567192528045E+00;
    x[37] = 58.0383778598867354584946224426E+00;
    x[38] = 61.4528455586293273580430507845E+00;
    x[39] = 64.9959270371996369214889214525E+00;
    x[40] = 68.6722926314626183977917033408E+00;
    x[41] = 72.4870626544527259157609330997E+00;
    x[42] = 76.4458687019847672439750100021E+00;
    x[43] = 80.5549265807842719390299008467E+00;
    x[44] = 84.8211237010524840445780520727E+00;
    x[45] = 89.2521246438779093867045988317E+00;
    x[46] = 93.8564998060605244752886814683E+00;
    x[47] = 98.6438836854008198627513004124E+00;
    x[48] = 103.625171719648548339736534290E+00;
    x[49] = 108.812767977808831403825388295E+00;
    x[50] = 114.220900975902615729136419234E+00;
    x[51] = 119.866032356536657109576887729E+00;
    x[52] = 125.767394661236316522917056970E+00;
    x[53] = 131.947712599442095787498603776E+00;
    x[54] = 138.434191890473005752666966171E+00;
    x[55] = 145.259909991527292147157412636E+00;
    x[56] = 152.465831757831086161761248192E+00;
    x[57] = 160.103837850755827943698932303E+00;
    x[58] = 168.241478626055331550636477652E+00;
    x[59] = 176.969855979856624954098184957E+00;
    x[60] = 186.417642483351089954653964415E+00;
    x[61] = 196.778474440876949259615139650E+00;
    x[62] = 208.372107380940491095442543113E+00;
    x[63] = 221.812376576320945562507410640E+00;
    x[64] = 238.685811594674270102373709659E+00;

    w[0] = 0.0554129011565536469555170462551E+00;
    w[1] = 0.117417396564162014386769936525E+00;
    w[2] = 0.15577793159527363526750345236E+00;
    w[3] = 0.166347955884031811873697185854E+00;
    w[4] = 0.153013207065446887512281359026E+00;
    w[5] = 0.12471061536737329712442529837E+00;
    w[6] = 0.0913671486268474804579363734334E+00;
    w[7] = 0.0606693673532322224974724445698E+00;
    w[8] = 0.0366985078337756899608575633553E+00;
    w[9] = 0.0202884358923229233158787730215E+00;
    w[10] = 0.0102732022687699783894639990081E+00;
    w[11] = 0.00477128781081110626106879095602E+00;
    w[12] = 0.00203437021965744474885076853755E+00;
    w[13] = 0.00079674216708789273511811886929E+00;
    w[14] = 0.000286683812097562728436084278246E+00;
    w[15] = 0.0000947743836779584423250711498309E+00;
    w[16] = 0.0000287808776386491122522878225936E+00;
    w[17] = 8.025921887785674361327140519E-6;
    w[18] = 2.0542534210210521080625753887E-6;
    w[19] = 4.822974163227069271800228998E-7;
    w[20] = 1.037906330975825689484858253E-7;
    w[21] = 2.04556458252437904249066042E-8;
    w[22] = 3.6885834829334325209249956E-9;
    w[23] = 6.07893150020096839531226E-10;
    w[24] = 9.14524981574057744811517E-11;
    w[25] = 1.25427414380884616209197E-11;
    w[26] = 1.56600212784699337922365E-12;
    w[27] = 1.7771074191193507379193E-13;
    w[28] = 1.829851898040443279934E-14;
    w[29] = 1.70645929547049792703E-15;
    w[30] = 1.438412537484673440005E-16;
    w[31] = 1.09354404172789017674E-17;
    w[32] = 7.4805627565827998507E-19;
    w[33] = 4.5927490046001844919E-20;
    w[34] = 2.5237976161210732972E-21;
    w[35] = 1.2376030276622370811E-22;
    w[36] = 5.398129612170324038E-24;
    w[37] = 2.086925656272980914E-25;
    w[38] = 7.12363601216912351E-27;
    w[39] = 2.13797949495352507E-28;
    w[40] = 5.61585996656691837E-30;
    w[41] = 1.2845455347843963E-31;
    w[42] = 2.5444410680246895E-33;
    w[43] = 4.33793017530204526E-35;
    w[44] = 6.3222072485073684E-37;
    w[45] = 7.8174607001716233E-39;
    w[46] = 8.1320382378781478E-41;
    w[47] = 7.0491916339430245E-43;
    w[48] = 5.03749289878460876E-45;
    w[49] = 2.93162035252999008E-47;
    w[50] = 1.370002490040754814E-49;
    w[51] = 5.05827802154356236E-52;
    w[52] = 1.447822639997111427E-54;
    w[53] = 3.141466864703300069E-57;
    w[54] = 5.03056294905086669E-60;
    w[55] = 5.7547901892214244076E-63;
    w[56] = 4.5172924325696051375E-66;
    w[57] = 2.31224468889972999556E-69;
    w[58] = 7.22311071834335770012E-73;
    w[59] = 1.2595483022441380495635E-76;
    w[60] = 1.08123622593537779961256E-80;
    w[61] = 3.78397004144107162756714E-85;
    w[62] = 3.9595563374024690284872814E-90;
    w[63] = 6.85929105947869810086088696E-96;
    w[64] = 4.3543678721167358517710196959E-103;
  }
  else if ( n == 127 )
  {
    x[0] = 0.0113396352985186116918931696313E+00;
    x[1] = 0.0597497534357266202813482370574E+00;
    x[2] = 0.146850986907461676123882236874E+00;
    x[3] = 0.272675907358595531313780082789E+00;
    x[4] = 0.437246006441926655545770358699E+00;
    x[5] = 0.640586882225669295335764164000E+00;
    x[6] = 0.882729686390583644814876536500E+00;
    x[7] = 1.16371141601665376615605847010E+00;
    x[8] = 1.48357501528346138913135848610E+00;
    x[9] = 1.84236943516135653806863208099E+00;
    x[10] = 2.24014968395790242445133156565E+00;
    x[11] = 2.67697687801413036921678699612E+00;
    x[12] = 3.15291829570828255657715083088E+00;
    x[13] = 3.66804743603047525402263399265E+00;
    x[14] = 4.22244408233018884559778766674E+00;
    x[15] = 4.81619437158705024756655350873E+00;
    x[16] = 5.44939086945594167558621789084E+00;
    x[17] = 6.12213265129972541939445847632E+00;
    x[18] = 6.83452538941226681122379949733E+00;
    x[19] = 7.58668144663674721742059868368E+00;
    x[20] = 8.37871997659327252548421206595E+00;
    x[21] = 9.21076703074265587779225061024E+00;
    x[22] = 10.0829556725286438091664393536E+00;
    x[23] = 10.9954260988581254298031473588E+00;
    x[24] = 11.9483257691977259976106051279E+00;
    x[25] = 12.9418095425855310537233810982E+00;
    x[26] = 13.9760398228785065200144056687E+00;
    x[27] = 15.0511867125795236315747963654E+00;
    x[28] = 16.1674281756128529229773950518E+00;
    x[29] = 17.3249502094436734465611637126E+00;
    x[30] = 18.5239470269656885608117113093E+00;
    x[31] = 19.7646212486115041040716693869E+00;
    x[32] = 21.0471841051731836068770440201E+00;
    x[33] = 22.3718556518555428176481239181E+00;
    x[34] = 23.7388649941224971836523137887E+00;
    x[35] = 25.1484505259373682340772783856E+00;
    x[36] = 26.6008601810417496072533842798E+00;
    x[37] = 28.0963516979646192017539612921E+00;
    x[38] = 29.6351928995041789106102271386E+00;
    x[39] = 31.2176619874797591442144671526E+00;
    x[40] = 32.8440478536104304605229513413E+00;
    x[41] = 34.5146504074411491491056359474E+00;
    x[42] = 36.2297809223068040196153885089E+00;
    x[43] = 37.9897624003999564359687801403E+00;
    x[44] = 39.7949299580899617783964371417E+00;
    x[45] = 41.6456312327301807051539908975E+00;
    x[46] = 43.5422268122868595499508929938E+00;
    x[47] = 45.4850906892287911379961513367E+00;
    x[48] = 47.4746107402319647194687665991E+00;
    x[49] = 49.5111892333790877167288845844E+00;
    x[50] = 51.5952433646712444431827712669E+00;
    x[51] = 53.7272058258193167582881400691E+00;
    x[52] = 55.9075254054475533058306059917E+00;
    x[53] = 58.1366676260224391970775260257E+00;
    x[54] = 60.4151154190185902957071920538E+00;
    x[55] = 62.7433698410518097002071267427E+00;
    x[56] = 65.1219508339499963119560254171E+00;
    x[57] = 67.5513980319978863144118724431E+00;
    x[58] = 70.0322716198845845112298711920E+00;
    x[59] = 72.5651532452068490908886694168E+00;
    x[60] = 75.1506469897399352993543623251E+00;
    x[61] = 77.7893804040858160006474054621E+00;
    x[62] = 80.4820056107507292058039629268E+00;
    x[63] = 83.2292004811959148867961200190E+00;
    x[64] = 86.0316698929535829667982387326E+00;
    x[65] = 88.8901470735120510996525185443E+00;
    x[66] = 91.8053950383581779949712501705E+00;
    x[67] = 94.7782081313315832053870310348E+00;
    x[68] = 97.8094136763051164110541101154E+00;
    x[69] = 100.899873750172859403719397622E+00;
    x[70] = 104.050487088215989347040768450E+00;
    x[71] = 107.262191134146004284231164014E+00;
    x[72] = 110.535964248515005306027713513E+00;
    x[73] = 113.872828090758394853483761877E+00;
    x[74] = 117.273850191925177740954778864E+00;
    x[75] = 120.740146737188801061739780027E+00;
    x[76] = 124.272885579556983542595064469E+00;
    x[77] = 127.873289508859426450938417454E+00;
    x[78] = 131.542639803143669218093777421E+00;
    x[79] = 135.282280093118369701327381064E+00;
    x[80] = 139.093620574329700139644220870E+00;
    x[81] = 142.978142606436017768082277536E+00;
    x[82] = 146.937403744373665494410809691E+00;
    x[83] = 150.973043252521871274925114375E+00;
    x[84] = 155.086788160346125722296414206E+00;
    x[85] = 159.280459926632882354019569899E+00;
    x[86] = 163.555981789575711040159671821E+00;
    x[87] = 167.915386891943601342455471847E+00;
    x[88] = 172.360827284738125368381561917E+00;
    x[89] = 176.894583929601921763116749935E+00;
    x[90] = 181.519077840368130692275288340E+00;
    x[91] = 186.236882528281123738612025304E+00;
    x[92] = 191.050737944509291967908366108E+00;
    x[93] = 195.963566148798798378390025430E+00;
    x[94] = 200.978488976000251536964755261E+00;
    x[95] = 206.098848024688711121272830428E+00;
    x[96] = 211.328227356716552605723772570E+00;
    x[97] = 216.670479376582303234770894658E+00;
    x[98] = 222.129754459296872462673049638E+00;
    x[99] = 227.710535020722324190891324313E+00;
    x[100] = 233.417674882826024533677753226E+00;
    x[101] = 239.256444988303086200187496671E+00;
    x[102] = 245.232586778715671725312540190E+00;
    x[103] = 251.352374887181280300055009918E+00;
    x[104] = 257.622691237920614130761918823E+00;
    x[105] = 264.051113229082405517543772418E+00;
    x[106] = 270.646019457227967492991117186E+00;
    x[107] = 277.416717501636510717983882181E+00;
    x[108] = 284.373599742208703266744028731E+00;
    x[109] = 291.528335213464957195812820216E+00;
    x[110] = 298.894108370282486008788956154E+00;
    x[111] = 306.485919782626113204181124239E+00;
    x[112] = 314.320969864711774874000075076E+00;
    x[113] = 322.419155891286796833494403613E+00;
    x[114] = 330.803726638024056519338473349E+00;
    x[115] = 339.502161278324337477353675960E+00;
    x[116] = 348.547375594726973554807617874E+00;
    x[117] = 357.979420280298454540490074431E+00;
    x[118] = 367.847945200760045788583414229E+00;
    x[119] = 378.215906231355328183329791889E+00;
    x[120] = 389.165391412510041015794753252E+00;
    x[121] = 400.807293314517025899963612864E+00;
    x[122] = 413.298536817793844180082600819E+00;
    x[123] = 426.875791536636755382885090171E+00;
    x[124] = 441.930854853108414124603092718E+00;
    x[125] = 459.218046398884299819712673132E+00;
    x[126] = 480.693782633883738597042692293E+00;

    w[0] = 0.0287732466920001243557700103008E+00;
    w[1] = 0.0638174681751346493634809492651E+00;
    w[2] = 0.0919196697215705713898641946531E+00;
    w[3] = 0.11054167914413766381245463003E+00;
    w[4] = 0.118797716333758501883283294227E+00;
    w[5] = 0.117378185300526951488044516301E+00;
    w[6] = 0.108193059841805514883351455812E+00;
    w[7] = 0.0938270752904896280803772614011E+00;
    w[8] = 0.0769664509605888439958224859284E+00;
    w[9] = 0.0599349039129397143325707300635E+00;
    w[10] = 0.0444177420738890013717083162729E+00;
    w[11] = 0.0313850809662523209830093722151E+00;
    w[12] = 0.021172316041924506411370709025E+00;
    w[13] = 0.0136501453642305416521711855646E+00;
    w[14] = 0.00841728527105991722793666573854E+00;
    w[15] = 0.00496749900598827605159128586202E+00;
    w[16] = 0.00280699038950018846319619574464E+00;
    w[17] = 0.00151929510039419524604453410578E+00;
    w[18] = 0.000787890287517960840862172871405E+00;
    w[19] = 0.00039156751064868450584507324649E+00;
    w[20] = 0.000186524342688258605500935662601E+00;
    w[21] = 0.0000851731604155766219088098281602E+00;
    w[22] = 0.0000372856391978530377121453215777E+00;
    w[23] = 0.0000156484167917129939474478052968E+00;
    w[24] = 6.2964340695224829035692735525E-06;
    w[25] = 2.42889297113287245745413799382E-06;
    w[26] = 8.9824607890051007201922871545E-07;
    w[27] = 3.18441747407603537107429663281E-07;
    w[28] = 1.08212729055668392118618075427E-07;
    w[29] = 3.52450767506355360159027790853E-08;
    w[30] = 1.10012243657193474070638397617E-08;
    w[31] = 3.29040796167179321253293430033E-09;
    w[32] = 9.4289145237889976419772700773E-10;
    w[33] = 2.58825789046683181840501953093E-10;
    w[34] = 6.80474371033707626309422590176E-11;
    w[35] = 1.71313988051208378353995644756E-11;
    w[36] = 4.12917445240528654694439223049E-12;
    w[37] = 9.52641897188072732207076648735E-13;
    w[38] = 2.10326044324424259329629420475E-13;
    w[39] = 4.44271519387293528609404342858E-14;
    w[40] = 8.97605003628337033233198464055E-15;
    w[41] = 1.73415114077692870746279483468E-15;
    w[42] = 3.20280995489883566314943798352E-16;
    w[43] = 5.65313889507936820226607420952E-17;
    w[44] = 9.53296727990265912345880440259E-18;
    w[45] = 1.53534534773101425652885094376E-18;
    w[46] = 2.36089621794673656860578421322E-19;
    w[47] = 3.46487427944566113321938766532E-20;
    w[48] = 4.85152418970864613201269576635E-21;
    w[49] = 6.47862286335198134281373737907E-22;
    w[50] = 8.2476020965403242936448553126E-23;
    w[51] = 1.0005361880214719793491658283E-23;
    w[52] = 1.1561395116207304954233181264E-24;
    w[53] = 1.271934273116792265561213426E-25;
    w[54] = 1.331658471416537296734000416E-26;
    w[55] = 1.32612184546789440336461085E-27;
    w[56] = 1.25549954476439498072860741E-28;
    w[57] = 1.1294412178579462703240913E-29;
    w[58] = 9.649102027956211922850061E-31;
    w[59] = 7.82418467683020993967331E-32;
    w[60] = 6.01815035422196266582499E-33;
    w[61] = 4.38824827049617415515105E-34;
    w[62] = 3.0314137647517256304036E-35;
    w[63] = 1.9826016543944539545225E-36;
    w[64] = 1.2267623373665926559014E-37;
    w[65] = 7.176393169250888894381E-39;
    w[66] = 3.965937883383696358411E-40;
    w[67] = 2.068897055386804009958E-41;
    w[68] = 1.017958701797951724527E-42;
    w[69] = 4.72008277459863746257E-44;
    w[70] = 2.06068289855533748257E-45;
    w[71] = 8.4627575907305987246E-47;
    w[72] = 3.2661123687088798658E-48;
    w[73] = 1.1833939207883162381E-49;
    w[74] = 4.021120912389501381E-51;
    w[75] = 1.2799824394111125389E-52;
    w[76] = 3.81238777475488465E-54;
    w[77] = 1.061205754270115677E-55;
    w[78] = 2.757144694720040359E-57;
    w[79] = 6.67725442409284929E-59;
    w[80] = 1.505243838386823495E-60;
    w[81] = 3.15389868001137585E-62;
    w[82] = 6.13266142994831808E-64;
    w[83] = 1.10485100303248106E-65;
    w[84] = 1.84105635380913481E-67;
    w[85] = 2.83239265700528322E-69;
    w[86] = 4.01544098437636555E-71;
    w[87] = 5.23515302156837088E-73;
    w[88] = 6.2634476665005101E-75;
    w[89] = 6.861221053566653E-77;
    w[90] = 6.8651298840956019E-79;
    w[91] = 6.25813884337280849E-81;
    w[92] = 5.1833271237514904E-83;
    w[93] = 3.88936215719184435E-85;
    w[94] = 2.63577113794769328E-87;
    w[95] = 1.60788512939179797E-89;
    w[96] = 8.79780420709689396E-92;
    w[97] = 4.30134050774951099E-94;
    w[98] = 1.871343588134283853E-96;
    w[99] = 7.212574470806047168E-99;
    w[100] = 2.450874606217787438E-101;
    w[101] = 7.304209461947087578E-104;
    w[102] = 1.8983290818383463538E-106;
    w[103] = 4.2757400244246684123E-109;
    w[104] = 8.2894681420515755691E-112;
    w[105] = 1.37294322193244000131E-114;
    w[106] = 1.926546412640497322204E-117;
    w[107] = 2.269334450330135482614E-120;
    w[108] = 2.2209290603717355061909E-123;
    w[109] = 1.7851087685544512662857E-126;
    w[110] = 1.16309319903871644674312E-129;
    w[111] = 6.05244435846523922909528E-133;
    w[112] = 2.472956911506352864762838E-136;
    w[113] = 7.778906500648941036499721E-140;
    w[114] = 1.8409738662712607039570678E-143;
    w[115] = 3.1900921131079114970179072E-147;
    w[116] = 3.917948713917419973761766608E-151;
    w[117] = 3.2782158394188697053774429821E-155;
    w[118] = 1.77935907131388880628196401287E-159;
    w[119] = 5.88823534089326231574678353812E-164;
    w[120] = 1.09572365090711698777472032739E-168;
    w[121] = 1.02816211148670008982850769758E-173;
    w[122] = 4.1704725557697758145816510854E-179;
    w[123] = 5.8002877720316101774638319602E-185;
    w[124] = 1.88735077458255171061716191011E-191;
    w[125] = 6.91066018267309116827867059509E-199;
    w[126] = 4.35068132011058556283833133344E-208;
  }
  else if ( n == 128 )
  {
    x[0] = 0.0112513882636759629608518403162E+00;
    x[1] = 0.0592847412690264542879220089614E+00;
    x[2] = 0.145707966594312465141854059102E+00;
    x[3] = 0.270553178758665066190760897100E+00;
    x[4] = 0.433841407553836803056096580754E+00;
    x[5] = 0.635597665781621938340867677969E+00;
    x[6] = 0.875852384546520779155346013261E+00;
    x[7] = 1.15464170197439795008153355708E+00;
    x[8] = 1.47200756316673547554446633038E+00;
    x[9] = 1.82799777831235028535984528718E+00;
    x[10] = 2.22266607156190244817896452914E+00;
    x[11] = 2.65607212988348119522885329309E+00;
    x[12] = 3.12828165502791695498310369738E+00;
    x[13] = 3.63936641985240321221074522169E+00;
    x[14] = 4.18940432959404797478493079865E+00;
    x[15] = 4.77847948843487609165724239213E+00;
    x[16] = 5.40668227160049918527893820105E+00;
    x[17] = 6.07410940319653684309155844506E+00;
    x[18] = 6.78086403997562541104804929943E+00;
    x[19] = 7.52705586122937588585512279842E+00;
    x[20] = 8.31280116500777060337884381191E+00;
    x[21] = 9.13822297088039239600262641969E+00;
    x[22] = 10.0034511294682220892682761435E+00;
    x[23] = 10.9086224389908825478488613010E+00;
    x[24] = 11.8538807690918568038332644538E+00;
    x[25] = 12.8393771922232496874551935673E+00;
    x[26] = 13.8652701228920803029536799971E+00;
    x[27] = 14.9317254650919274737473553133E+00;
    x[28] = 16.0389167682670793509213783428E+00;
    x[29] = 17.1870253921812651027585044591E+00;
    x[30] = 18.3762406810896949333827523370E+00;
    x[31] = 19.6067601476416467279054690989E+00;
    x[32] = 20.8787896669713729158932014403E+00;
    x[33] = 22.1925436814678369066763923182E+00;
    x[34] = 23.5482454167489205609249730097E+00;
    x[35] = 24.9461271094034886279510396640E+00;
    x[36] = 26.3864302471052908976269305132E+00;
    x[37] = 27.8694058217463902295696818564E+00;
    x[38] = 29.3953145962849137215656288381E+00;
    x[39] = 30.9644273860527540023220861317E+00;
    x[40] = 32.5770253553237501781419486456E+00;
    x[41] = 34.2334003300022426604794108753E+00;
    x[42] = 35.9338551273561538107722924963E+00;
    x[43] = 37.6787039037883744300655582016E+00;
    x[44] = 39.4682725217157641271489033004E+00;
    x[45] = 41.3028989367070896380080417637E+00;
    x[46] = 43.1829336061203832438783635225E+00;
    x[47] = 45.1087399205772317441506148507E+00;
    x[48] = 47.0806946597172168560725351128E+00;
    x[49] = 49.0991884737910268021535852860E+00;
    x[50] = 51.1646263927766594446404335916E+00;
    x[51] = 53.2774283648407739161085367944E+00;
    x[52] = 55.4380298261178918683089638291E+00;
    x[53] = 57.6468823039452288144249811220E+00;
    x[54] = 59.9044540558720556965292635062E+00;
    x[55] = 62.2112307469614582456791552962E+00;
    x[56] = 64.5677161681212154290410515467E+00;
    x[57] = 66.9744329984415610548027156195E+00;
    x[58] = 69.4319236147834299557621097742E+00;
    x[59] = 71.9407509521543751573018481062E+00;
    x[60] = 74.5014994187340277930279831855E+00;
    x[61] = 77.1147758697705705283198924354E+00;
    x[62] = 79.7812106449685528544582124991E+00;
    x[63] = 82.5014586744314529140391768845E+00;
    x[64] = 85.2762006587153587377964042582E+00;
    x[65] = 88.1061443290995036940317393258E+00;
    x[66] = 90.9920257947926131560303030245E+00;
    x[67] = 93.9346109844796944955244642925E+00;
    x[68] = 96.9346971903819404516199551240E+00;
    x[69] = 99.9931147238642715216213000267E+00;
    x[70] = 103.110728692593987392319749158E+00;
    x[71] = 106.288440910345442668426129659E+00;
    x[72] = 109.527191951777550806618056918E+00;
    x[73] = 112.827963365904193877333487264E+00;
    x[74] = 116.191780063556780940235871708E+00;
    x[75] = 119.619712895932010462348887420E+00;
    x[76] = 123.112881443360190060911814509E+00;
    x[77] = 126.672457035760183662338694957E+00;
    x[78] = 130.299666028913462587217492864E+00;
    x[79] = 133.995793363747964343582120836E+00;
    x[80] = 137.762186439339380964302666578E+00;
    x[81] = 141.600259334393040305789642722E+00;
    x[82] = 145.511497416659592393640597008E+00;
    x[83] = 149.497462385177707088173175451E+00;
    x[84] = 153.559797796566440117982748261E+00;
    x[85] = 157.700235133978105059095336546E+00;
    x[86] = 161.920600485975634163753629031E+00;
    x[87] = 166.222821912768875092875739160E+00;
    x[88] = 170.608937589242234646550663310E+00;
    x[89] = 175.081104828414604880617405502E+00;
    x[90] = 179.641610105866994602634964639E+00;
    x[91] = 184.292880225846805341834505020E+00;
    x[92] = 189.037494793954109001292998345E+00;
    x[93] = 193.878200190472967540802875940E+00;
    x[94] = 198.817925273720602804745944994E+00;
    x[95] = 203.859799085769844571664824897E+00;
    x[96] = 209.007170885510867853387511181E+00;
    x[97] = 214.263632898788021280527758492E+00;
    x[98] = 219.633046255578174038387401024E+00;
    x[99] = 225.119570684209027756659796566E+00;
    x[100] = 230.727698658203619681658868680E+00;
    x[101] = 236.462294850177665966018904158E+00;
    x[102] = 242.328641949702698267864519866E+00;
    x[103] = 248.332494162357178892016601780E+00;
    x[104] = 254.480140044869131893543803358E+00;
    x[105] = 260.778476773579736538560064538E+00;
    x[106] = 267.235098528953836763992472029E+00;
    x[107] = 273.858402462693609793414602648E+00;
    x[108] = 280.657716776323492397504100977E+00;
    x[109] = 287.643456899219330638473677900E+00;
    x[110] = 294.827317787647739179806672104E+00;
    x[111] = 302.222513246449465380535981711E+00;
    x[112] = 309.844077326612663447772363643E+00;
    x[113] = 317.709248954906289495678052340E+00;
    x[114] = 325.837970121194949650401277931E+00;
    x[115] = 334.253542067654135375184450174E+00;
    x[116] = 342.983506273825316408508913329E+00;
    x[117] = 352.060853546526185083043426984E+00;
    x[118] = 361.525726392325047599066851839E+00;
    x[119] = 371.427889214327523078517984867E+00;
    x[120] = 381.830444119061080196207616882E+00;
    x[121] = 392.815671240808098809377819898E+00;
    x[122] = 404.494724750515074389666071660E+00;
    x[123] = 417.024902977989015820197277594E+00;
    x[124] = 430.643444166597381558323551668E+00;
    x[125] = 445.743096973927989652171720726E+00;
    x[126] = 463.080034109446258208013793406E+00;
    x[127] = 484.615543986443976044063131110E+00;

    w[0] = 0.0285518444532397286290731773612E+00;
    w[1] = 0.0633502117845051187797978127259E+00;
    w[2] = 0.0913083813661343144231616325903E+00;
    w[3] = 0.109913900410911746101013915833E+00;
    w[4] = 0.118274171034173698789809688874E+00;
    w[5] = 0.117045739000406721566458439207E+00;
    w[6] = 0.108089987545568415436473783125E+00;
    w[7] = 0.0939428886389285017878088436356E+00;
    w[8] = 0.0772536687978980532077800252359E+00;
    w[9] = 0.0603270562656615705389303086003E+00;
    w[10] = 0.0448473482471952140682424657998E+00;
    w[11] = 0.0317969479368768461739632484821E+00;
    w[12] = 0.0215301494537944439261107285438E+00;
    w[13] = 0.0139369517338463483277576885975E+00;
    w[14] = 0.00863158538020224714884473096489E+00;
    w[15] = 0.00511777701366922852873936722845E+00;
    w[16] = 0.00290634743648595585817980077219E+00;
    w[17] = 0.00158143294331667939723416013489E+00;
    w[18] = 0.000824738985098812150435438593253E+00;
    w[19] = 0.000412326088539694730970290830804E+00;
    w[20] = 0.000197649426442591498620529889783E+00;
    w[21] = 0.0000908515788782451508022826306153E+00;
    w[22] = 0.0000400484927835805298977887660442E+00;
    w[23] = 0.0000169307623980817855755102888475E+00;
    w[24] = 6.86452529111068208938636278412E-06;
    w[25] = 2.66921659814210266015872228584E-06;
    w[26] = 9.95364010286384477177483332196E-07;
    w[27] = 3.55943575300306543988020166563E-07;
    w[28] = 1.22053255194881194831205615734E-07;
    w[29] = 4.01279192093563506890167766024E-08;
    w[30] = 1.26481141474759786445650110908E-08;
    w[31] = 3.82148972942657229023411372003E-09;
    w[32] = 1.10664105922734169994480044024E-09;
    w[33] = 3.07100923709742319582290034639E-10;
    w[34] = 8.16549938415448956026437885004E-11;
    w[35] = 2.07985363278137784234189612116E-11;
    w[36] = 5.0739537708398704043296986402E-12;
    w[37] = 1.1853143771796112305093733131E-12;
    w[38] = 2.65092752372887358600565488195E-13;
    w[39] = 5.67463221575765876681065606161E-14;
    w[40] = 1.16237381470751589221529434901E-14;
    w[41] = 2.27776629270238637919733104451E-15;
    w[42] = 4.26883197029764927739172104126E-16;
    w[43] = 7.64928879936327510525948457803E-17;
    w[44] = 1.31013139198382464188082886821E-17;
    w[45] = 2.1441452341246636343706788692E-18;
    w[46] = 3.35194428720884780801470729044E-19;
    w[47] = 5.00373308645947376823179365121E-20;
    w[48] = 7.13003064195856212049702464626E-21;
    w[49] = 9.6945407403972664035320905829E-22;
    w[50] = 1.25728475563978459844059927432E-22;
    w[51] = 1.5546610955630634482202731199E-23;
    w[52] = 1.832109793253421778719084254E-24;
    w[53] = 2.056797978136734920722781372E-25;
    w[54] = 2.19866605262329119257657449E-26;
    w[55] = 2.23691600732428936729406222E-27;
    w[56] = 2.1649606446339054400256309E-28;
    w[57] = 1.9922276806187937873877251E-29;
    w[58] = 1.742153886325439585907653E-30;
    w[59] = 1.446949786106284637699605E-31;
    w[60] = 1.1407517061230822834189E-32;
    w[61] = 8.5318050978102090722116E-34;
    w[62] = 6.04970117793885843505E-35;
    w[63] = 4.0643432432648003017795E-36;
    w[64] = 2.585349374987909630703E-37;
    w[65] = 1.556028762522623447585E-38;
    w[66] = 8.85462584966333001103E-40;
    w[67] = 4.76045751736458068032E-41;
    w[68] = 2.416078510661232205E-42;
    w[69] = 1.15664705033873749321E-43;
    w[70] = 5.2185106194923759952E-45;
    w[71] = 2.2169743353361803305E-46;
    w[72] = 8.86010275661369606E-48;
    w[73] = 3.327811159201095553E-49;
    w[74] = 1.173490043078302544E-50;
    w[75] = 3.880967726420921431E-52;
    w[76] = 1.202426327933061418E-53;
    w[77] = 3.48602304410554638E-55;
    w[78] = 9.44554522159556681E-57;
    w[79] = 2.38888427455968395E-58;
    w[80] = 5.63188475075463052E-60;
    w[81] = 1.23592861191216019E-61;
    w[82] = 2.52100420237726743E-63;
    w[83] = 4.7722246219998052E-65;
    w[84] = 8.3700198919995783E-67;
    w[85] = 1.35782434112020985E-68;
    w[86] = 2.03368872715315416E-70;
    w[87] = 2.8068384806953538E-72;
    w[88] = 3.562567607062096E-74;
    w[89] = 4.1494527492937706E-76;
    w[90] = 4.4250079657663219E-78;
    w[91] = 4.3100842612898497E-80;
    w[92] = 3.8246610167617398E-82;
    w[93] = 3.08354784259879275E-84;
    w[94] = 2.25213982217062084E-86;
    w[95] = 1.48551474064504312E-88;
    w[96] = 8.8196354763726564E-91;
    w[97] = 4.69641782212598507E-93;
    w[98] = 2.23439382545477274E-95;
    w[99] = 9.45878703822074032E-98;
    w[100] = 3.546960831240672614E-100;
    w[101] = 1.17253213003488723E-102;
    w[102] = 3.399090555639915548E-105;
    w[103] = 8.591907200623898045E-108;
    w[104] = 1.8818913973535359647E-110;
    w[105] = 3.5473586323062565237E-113;
    w[106] = 5.7114822282836004745E-116;
    w[107] = 7.78947378804446095611E-119;
    w[108] = 8.91589869949126935148E-122;
    w[109] = 8.476856358868403207418E-125;
    w[110] = 6.617326935494900345408E-128;
    w[111] = 4.1862163574157095190077E-131;
    w[112] = 2.11438516898114207120093E-134;
    w[113] = 8.38216350136786953641675E-138;
    w[114] = 2.557202302197677884687798E-141;
    w[115] = 5.8667686421912043720461236E-145;
    w[116] = 9.8498610300648438019885689E-149;
    w[117] = 1.171383943342068942456857274E-152;
    w[118] = 9.483963265567383663821702301E-157;
    w[119] = 4.9770963811238028116653976343E-161;
    w[120] = 1.59089852775099765481980638695E-165;
    w[121] = 2.85630382911900292320607568044E-170;
    w[122] = 2.58225071969148999265031459122E-175;
    w[123] = 1.00735025005079740952983187255E-180;
    w[124] = 1.34425250044381631821772983363E-186;
    w[125] = 4.18296221403683473389726627221E-193;
    w[126] = 1.45716530772618631594481663188E-200;
    w[127] = 8.64059169046870867692891422354E-210;
  }
  else if ( n == 129 )
  {
    x[0] = 0.0111645041367687260935881187114E+00;
    x[1] = 0.0588269115255121725669144777376E+00;
    x[2] = 0.144582603939087375345544455104E+00;
    x[3] = 0.268463250498790809142537571727E+00;
    x[4] = 0.430489433028069665583513882755E+00;
    x[5] = 0.630685596971157529700818698614E+00;
    x[6] = 0.869081474989540465988995980646E+00;
    x[7] = 1.14571237269034129786358037349E+00;
    x[8] = 1.46061926689785560022252086358E+00;
    x[9] = 1.81384886225287260620048305182E+00;
    x[10] = 2.20545363849013952710373368048E+00;
    x[11] = 2.63549189753739459262727316570E+00;
    x[12] = 3.10402781353627480023526416641E+00;
    x[13] = 3.61113148701933289479734007535E+00;
    x[14] = 4.15687900382495881133416031205E+00;
    x[15] = 4.74135249908325871733484826319E+00;
    x[16] = 5.36464022650680264548807369539E+00;
    x[17] = 6.02683663318167548105631862177E+00;
    x[18] = 6.72804244004243132025609101021E+00;
    x[19] = 7.46836472821534963467632383543E+00;
    x[20] = 8.24791703142169723816558449856E+00;
    x[21] = 9.06681943464370270026626900050E+00;
    x[22] = 9.92519867926931734070041188408E+00;
    x[23] = 10.8231882749469306495297612192E+00;
    x[24] = 11.7609286183977310387181197615E+00;
    x[25] = 12.7385671194512351722694605084E+00;
    x[26] = 13.7562583345886271805335101149E+00;
    x[27] = 14.8141641082989854857712290559E+00;
    x[28] = 15.9124537225752979381236294324E+00;
    x[29] = 17.0513040549004685335351914932E+00;
    x[30] = 18.2308997450984136591429080617E+00;
    x[31] = 19.4514333714519620150207362048E+00;
    x[32] = 20.7131056365177555775299262985E+00;
    x[33] = 22.0161255630988608706489924430E+00;
    x[34] = 23.3607107008685190470514486749E+00;
    x[35] = 24.7470873441735867407744490421E+00;
    x[36] = 26.1754907615839641134296855243E+00;
    x[37] = 27.6461654377949106206644801830E+00;
    x[38] = 29.1593653285328756045576144321E+00;
    x[39] = 30.7153541291626095441732915451E+00;
    x[40] = 32.3144055577441922161871665693E+00;
    x[41] = 33.9568036533435689847296719094E+00;
    x[42] = 35.6428430904596160112634717165E+00;
    x[43] = 37.3728295104950910213327545948E+00;
    x[44] = 39.1470798712685466663878582421E+00;
    x[45] = 40.9659228156399190364649448364E+00;
    x[46] = 42.8296990604046437906422357564E+00;
    x[47] = 44.7387618067004519884778950119E+00;
    x[48] = 46.6934771732681867037686990052E+00;
    x[49] = 48.6942246540138734219622380649E+00;
    x[50] = 50.7413976014347803131042845818E+00;
    x[51] = 52.8354037375983340937979164025E+00;
    x[52] = 54.9766656945006500240481182310E+00;
    x[53] = 57.1656215857823649284158070179E+00;
    x[54] = 59.4027256119448606421881531943E+00;
    x[55] = 61.6884487013914405461822377003E+00;
    x[56] = 64.0232791898173852597210780437E+00;
    x[57] = 66.4077235406921080587914699631E+00;
    x[58] = 68.8423071098181639647332557636E+00;
    x[59] = 71.3275749572182499453797757024E+00;
    x[60] = 73.8640927098955268421782575110E+00;
    x[61] = 76.4524474793379566181613942983E+00;
    x[62] = 79.0932488379977030145472340597E+00;
    x[63] = 81.7871298593763443093790704140E+00;
    x[64] = 84.5347482267906647046323684124E+00;
    x[65] = 87.3367874163878117910865422310E+00;
    x[66] = 90.1939579605291478450570652459E+00;
    x[67] = 93.1069987982766656767050611186E+00;
    x[68] = 96.0766787204029972427806506124E+00;
    x[69] = 99.1037979171157474398207782757E+00;
    x[70] = 102.189189637550616978355114969E+00;
    x[71] = 105.333721971058838012388514189E+00;
    x[72] = 108.538299761408281119757506569E+00;
    x[73] = 111.803866666252185387269569516E+00;
    x[74] = 115.131407375615803792171876281E+00;
    x[75] = 118.521950004733905726449958829E+00;
    x[76] = 121.976568678369858697173472594E+00;
    x[77] = 125.496386325793836628100280130E+00;
    x[78] = 129.082577707933597477650969878E+00;
    x[79] = 132.736372700883616552797038522E+00;
    x[80] = 136.459059863023413416361147154E+00;
    x[81] = 140.251990316520692590651584246E+00;
    x[82] = 144.116581978059547282038191264E+00;
    x[83] = 148.054324178334554730971189024E+00;
    x[84] = 152.066782715303825347545842677E+00;
    x[85] = 156.155605392537787829354492826E+00;
    x[86] = 160.322528101405362530717313062E+00;
    x[87] = 164.569381514511906899139962575E+00;
    x[88] = 168.898098467996847713856358122E+00;
    x[89] = 173.310722122324145053009369479E+00;
    x[90] = 177.809415005439611927392788370E+00;
    x[91] = 182.396469059102149766157602559E+00;
    x[92] = 187.074316829415599837320402996E+00;
    x[93] = 191.845543966839763114071401295E+00;
    x[94] = 196.712903230183761859963922576E+00;
    x[95] = 201.679330224475912387142876872E+00;
    x[96] = 206.747961145685983577272619640E+00;
    x[97] = 211.922152858007833039218477193E+00;
    x[98] = 217.205505694330143211608701365E+00;
    x[99] = 222.601889450939732907762241579E+00;
    x[100] = 228.115473147766504188912042615E+00;
    x[101] = 233.750759251359480867215068547E+00;
    x[102] = 239.512623216994726852324857048E+00;
    x[103] = 245.406359409276117119061170837E+00;
    x[104] = 251.437734721509742305967783880E+00;
    x[105] = 257.613051552607945102309836371E+00;
    x[106] = 263.939222243647173894814944292E+00;
    x[107] = 270.423857663083214051346532320E+00;
    x[108] = 277.075373415313344577287378499E+00;
    x[109] = 283.903118212107869887400929941E+00;
    x[110] = 290.917530409009503510470042900E+00;
    x[111] = 298.130330747241946479391511151E+00;
    x[112] = 305.554762228622700877217556637E+00;
    x[113] = 313.205892212538716296350101328E+00;
    x[114] = 321.100997941634721519100026717E+00;
    x[115] = 329.260065894410473958350680155E+00;
    x[116] = 337.706449515634131989920236326E+00;
    x[117] = 346.467752279350659621376501841E+00;
    x[118] = 355.577039643063413893183224979E+00;
    x[119] = 365.074545471124791778391196263E+00;
    x[120] = 375.010148136708978052975802762E+00;
    x[121] = 385.447095254054417308720464640E+00;
    x[122] = 396.467858100744210106334636127E+00;
    x[123] = 408.183851152492844798297769341E+00;
    x[124] = 420.752744334742187498526476928E+00;
    x[125] = 434.412341688764625555428748148E+00;
    x[126] = 449.556338392256949417199002480E+00;
    x[127] = 466.942750921706688536121321308E+00;
    x[128] = 488.537715007400745716181291102E+00;

    w[0] = 0.0283338232816188129433412493366E+00;
    w[1] = 0.0628897352309939992519628028429E+00;
    w[2] = 0.0907050560197830441591715791845E+00;
    w[3] = 0.109292734964339745013347523543E+00;
    w[4] = 0.117753891824430328742552706746E+00;
    w[5] = 0.116712333575132760088854393741E+00;
    w[6] = 0.1079821092277907522768638822E+00;
    w[7] = 0.0940513886437790878162542877426E+00;
    w[8] = 0.0775328171368385256641246588694E+00;
    w[9] = 0.0607119801995722871258201910352E+00;
    w[10] = 0.0452716214541695196710137988047E+00;
    w[11] = 0.0322057586869443933590250840601E+00;
    w[12] = 0.0218870093879284288723521418152E+00;
    w[13] = 0.0142243242185532561642375502974E+00;
    w[14] = 0.00884734285745239479408590424342E+00;
    w[15] = 0.00526983370954167607842815218011E+00;
    w[16] = 0.00300740619275414763773247784756E+00;
    w[17] = 0.00164498171784021535901621253553E+00;
    w[18] = 0.000862641473273809069700952476134E+00;
    w[19] = 0.000433807488545501081264834235514E+00;
    w[20] = 0.000209234988721404556453070968853E+00;
    w[21] = 0.0000968044053231071525887634259114E+00;
    w[22] = 0.0000429650601010182583779356860953E+00;
    w[23] = 0.0000182943298240488545326843922155E+00;
    w[24] = 7.47320473307839845584026474317E-06;
    w[25] = 2.92876004890558731746712968433E-06;
    w[26] = 1.10111937532188602299646730309E-06;
    w[27] = 3.97133727854894494886436944708E-07;
    w[28] = 1.37391737873739072964678053016E-07;
    w[29] = 4.55898285044463980401770363171E-08;
    w[30] = 1.45082031554226827387170770004E-08;
    w[31] = 4.42736861865778798052557346184E-09;
    w[32] = 1.29540549841465072618582643105E-09;
    w[33] = 3.63353401896969889688016611161E-10;
    w[34] = 9.76889957112077658988662065766E-11;
    w[35] = 2.51697359198850123687093430919E-11;
    w[36] = 6.2136427115425329244941688694E-12;
    w[37] = 1.46947065273427272155102255087E-12;
    w[38] = 3.3283536396381226771786168693E-13;
    w[39] = 7.21860543546415622515782097245E-14;
    w[40] = 1.49874700296546634758941598894E-14;
    w[41] = 2.97813865190408297766537928957E-15;
    w[42] = 5.66223500996744709699260363288E-16;
    w[43] = 1.02976110977345161229212606736E-16;
    w[44] = 1.79087076765055918801501255712E-17;
    w[45] = 2.97741214327584722879794953728E-18;
    w[46] = 4.73066849378813640521244315218E-19;
    w[47] = 7.18076704552391091114386577815E-20;
    w[48] = 1.0409591754013912471892470954E-20;
    w[49] = 1.44063705945958837668569771815E-21;
    w[50] = 1.9027009013059586477368991424E-22;
    w[51] = 2.3972421860336028068385342016E-23;
    w[52] = 2.880065029076382866335001882E-24;
    w[53] = 3.2980570110683255202892323E-25;
    w[54] = 3.59822818119059018987046195E-26;
    w[55] = 3.7384843519427824153681456E-27;
    w[56] = 3.6971969670644497346136084E-28;
    w[57] = 3.478607942989822329014257E-29;
    w[58] = 3.112229078360896126467536E-30;
    w[59] = 2.64630166366922810478446E-31;
    w[60] = 2.13731385180863223984415E-32;
    w[61] = 1.63873356712820982018691E-33;
    w[62] = 1.1920639048111247727415E-34;
    w[63] = 8.221888191494076473793E-36;
    w[64] = 5.373327742595686629791E-37;
    w[65] = 3.325235584661609413228E-38;
    w[66] = 1.94717317556033610096E-39;
    w[67] = 1.07813497736466418105E-40;
    w[68] = 5.6402504582069233692E-42;
    w[69] = 2.785716667292756732E-43;
    w[70] = 1.2978694111929463222E-44;
    w[71] = 5.699117216622829387E-46;
    w[72] = 2.356556045713220169E-47;
    w[73] = 9.167179452095711245E-49;
    w[74] = 3.351643630271094859E-50;
    w[75] = 1.15053967361148792E-51;
    w[76] = 3.70428664291287775E-53;
    w[77] = 1.117334474142203311E-54;
    w[78] = 3.15377989811063792E-56;
    w[79] = 8.319920981942047E-58;
    w[80] = 2.04876111892933112E-59;
    w[81] = 4.7028955186049464E-61;
    w[82] = 1.00491633674668433E-62;
    w[83] = 1.9959187047623038E-64;
    w[84] = 3.6789923736675531E-66;
    w[85] = 6.2831482675040959E-68;
    w[86] = 9.925201342288209E-70;
    w[87] = 1.4475221077412768E-71;
    w[88] = 1.945364935931307E-73;
    w[89] = 2.4042822695448614E-75;
    w[90] = 2.7267496829701407E-77;
    w[91] = 2.8313374255297656E-79;
    w[92] = 2.6851895059223692E-81;
    w[93] = 2.3199549783717045E-83;
    w[94] = 1.821032672647817E-85;
    w[95] = 1.29486019972133753E-87;
    w[96] = 8.3146100960594316E-90;
    w[97] = 4.8053665090563748E-92;
    w[98] = 2.49071240066108676E-94;
    w[99] = 1.15335704284873844E-96;
    w[100] = 4.75169815023478164E-99;
    w[101] = 1.733951399870136754E-101;
    w[102] = 5.57731896834145892E-104;
    w[103] = 1.573010564351007982E-106;
    w[104] = 3.867845242632879313E-109;
    w[105] = 8.239883435606238718E-112;
    w[106] = 1.5104570697877326124E-114;
    w[107] = 2.3645657754433596259E-117;
    w[108] = 3.1349053289923477642E-120;
    w[109] = 3.48739145376585928069E-123;
    w[110] = 3.22170074744057989255E-126;
    w[111] = 2.443048415722317309221E-129;
    w[112] = 1.5008657805760609578501E-132;
    w[113] = 7.3592251345721592465131E-136;
    w[114] = 2.83121162238276011127992E-139;
    w[115] = 8.3785828758598937096069E-143;
    w[116] = 1.8637689328976254234922931E-146;
    w[117] = 3.0323700940390393081087066E-150;
    w[118] = 3.49260330326226204565809172E-154;
    w[119] = 2.736761201290944128360070077E-158;
    w[120] = 1.3888959774881077581342370711E-162;
    w[121] = 4.28912860126508716947322409477E-167;
    w[122] = 7.43133882324715291928018394696E-172;
    w[123] = 6.47421443374096511679045401121E-177;
    w[124] = 2.42953692988216878005255824922E-182;
    w[125] = 3.11143287762562176520260181694E-188;
    w[126] = 9.26127289624597363219192415542E-195;
    w[127] = 3.07023341560782650495387872798E-202;
    w[128] = 1.71530871887294016615286222244E-211;
  }
  else
  {
    cerr << "\n";
    cerr << "LAGUERRE_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1 to 20, 31/32/33, 63/64/65 or 127/128/129.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void laguerre_ss_compute ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_SS_COMPUTE computes a Gauss-Laguerre quadrature rule.
//
//  Discussion:
//
//    If the integral is:
//
//      Integral ( 0 <= X < +oo ) EXP ( - X ) * F(X) dX
//
//    then the quadrature rule is:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//    If the integral is:
//
//      Integral ( A <= X < +oo ) F(X) dX
//
//    then the quadrature rule is:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the quadrature rule to be computed.
//    ORDER must be at least 1.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double *b;
  double *c;
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double ratio;
  double x;

  b = new double[order];
  c = new double[order];
//
//  Set the recursion coefficients.
//
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( double ) ( 2 * i + 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i * i );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = prod;

  for ( i = 0; i < order; i++ )
  {
//
//  Compute an estimate for the root.
//
    if ( i == 0 )
    {
      x =  3.0 / ( 1.0 + 2.4 * ( double ) ( order ) );
    }
    else if ( i == 1 )
    {
      x = x + 15.0 / ( 1.0 + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) ) 
        / ( 1.9 * ( double ) ( i - 1 ) );

      x = x + r1 * ( x - xtab[i-2] );
    }
//
//  Use iteration to find the root.
//
    laguerre_ss_root ( &x, order, &dp2, &p1, b, c );
//
//  Set the abscissa and weight.
//
    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;
  }

  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void laguerre_ss_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_SS_RECUR evaluates a Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of L(ORDER)(X).
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial.
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
//****************************************************************************80

void laguerre_ss_root ( double *x, int order, double *dp2, double *p1, 
  double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_SS_ROOT improves a root of a Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    laguerre_ss_recur ( &p2, dp2, p1, *x, order, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

double laguerre_sum ( double func ( double x ), double a, int order, 
  double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_SUM carries out Laguerre quadrature over [ A, +oo ).
//
//  Discussion:
//
//    The simplest Laguerre integral is the
//    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
//    it is easy to modify the rule to approximate the integral from
//    A to +oo as well.
//
//    Another common Laguerre integral is the
//    integral from 0 to +oo of EXP(-X) * X**ALPHA * F(X).
//    This routine may be used to sum up the terms of the Laguerre
//    rule for such an integral as well.  However, if ALPHA is nonzero,
//    then there is no simple way to extend the rule to approximate the
//    integral from A to +oo.  The simplest procedures would be
//    to approximate the integral from 0 to A.
//
//    The integral:
//
//      Integral ( A <= X <= +oo ) EXP ( -X ) * F(X) dX
//    or
//      Integral ( 0 <= X <= +oo ) EXP ( -X ) * X^ALPHA * F(X) dX
//
//    The quadrature rule:
//
//      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) + A ) 
//    or
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, double FUNC ( double X ), the name of the function which
//    evaluates the integrand.  
//
//    Input, double A, the beginning of the integration interval.
//
//    Input, int ORDER, the order.
//
//    Input, double XTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights.
//
//    Output, double LAGUERRE_SUM, the approximate value of the integral.
//
{
  int i;
  double result;

  if ( order < 1 )
  {
    cerr << "\n";
    cerr << "LAGUERRE_SUM - Fatal error!\n";
    cerr << "  Nonpositive ORDER = " << order << "\n";
    exit ( 1 );
  }

  result = 0.0;
  for ( i = 0; i < order; i++ )
  {
    result = result + weight[i] * func ( xtab[i] + a );
  }
  result = exp ( - a ) * result;

  return result;
}
//****************************************************************************80

void legendre_dr_compute ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_DR_COMPUTE: Gauss-Legendre quadrature by Davis-Rabinowitz method.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2007
//
//  Author:
//
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be greater than 0.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    cerr << "\n";
    cerr << "LEGENDRE_DR_COMPUTE - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    exit ( 1 );
  }

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) ) 
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
//
//  Initial approximation H:
//
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
//
//  Refine H using one step of Newton's method:
//
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    weight[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
  }

  if ( ( order % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }
//
//  Shift the data up.
//
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
    weight[iback-1] = weight[iback-ncopy-1];
  }
//
//  Reflect values for the negative abscissas.
//
  for ( i = 1; i <= order - nmove; i++ )
  {
    xtab[i-1] = - xtab[order-i];
    weight[i-1] = weight[order-i];
  }

  return;
}
//****************************************************************************80

void legendre_ek_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_EK_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *bj;
  int i;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = 2.0;
//
//  Define the Jacobi matrix.
//
  bj = new double[n];

  for ( i = 0; i < n; i++ )
  {
    bj[i] = ( double ) ( ( i + 1 ) * ( i + 1 ) ) 
          / ( double ) ( 4 * ( i + 1 ) * ( i + 1 ) - 1 );
    bj[i] = sqrt ( bj[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  w[0] = sqrt ( zemu );

  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double legendre_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= +1 ) x^n dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Output, double LEGENDRE_INTEGRAL, the value of the exact integral.
//
{
  double exact;
//
//  Get the exact value of the integral.
//
  if ( ( expon % 2 ) == 0 )
  {
    exact = 2.0 / ( double ) ( expon + 1 );
  }
  else
  {
    exact = 0.0;
  }

  return exact;
}
//****************************************************************************80

void legendre_recur ( double *p2, double *dp2, double *p1, double x, int order )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_RECUR finds the value and derivative of a Legendre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of P(ORDER)(X).
//
//    Output, double *DP2, the value of P'(ORDER)(X).
//
//    Output, double *P1, the value of P(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x;
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = ( ( double ) ( 2 * i - 1 ) * x * ( *p1 ) 
          + ( double ) (   - i + 1 )     * p0 ) 
          / ( double ) (     i     );

    *dp2 = ( ( double ) ( 2 * i - 1 ) * ( ( *p1 ) + x * dp1 ) 
           - ( double ) (     i - 1 ) * dp0 ) 
           / ( double ) (     i     );
  }

  return;
}
//****************************************************************************80

void legendre_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    Quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//    The quadrature rule is exact for all polynomials through degree 2*N-1.
//
//    The abscissas are the zeroes of the Legendre polynomial P(N)(X).
//
//    Mathematica can compute the abscissas and weights of a Gauss-Legendre
//    rule of order N for the interval [A,B] with P digits of precision
//    by the commands:
//
//    Needs["NumericalDifferentialEquationAnalysis`"]
//    GaussianQuadratureWeights [ n, a, b, p ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798.
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 33 or 63/64/65, 127/128/129, 
//    255/256/257.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[0] = 0.000000000000000000000000000000;

    w[0] = 2.000000000000000000000000000000;
  }
  else if ( n == 2 )
  {
    x[0] = -0.577350269189625764509148780502;
    x[1] = 0.577350269189625764509148780502;

    w[0] = 1.000000000000000000000000000000;
    w[1] = 1.000000000000000000000000000000;
  }
  else if ( n == 3 )
  {
    x[0] = -0.774596669241483377035853079956;
    x[1] = 0.000000000000000000000000000000;
    x[2] = 0.774596669241483377035853079956;

    w[0] = 0.555555555555555555555555555556;
    w[1] = 0.888888888888888888888888888889;
    w[2] = 0.555555555555555555555555555556;
  }
  else if ( n == 4 )
  {
    x[0] = -0.861136311594052575223946488893;
    x[1] = -0.339981043584856264802665759103;
    x[2] = 0.339981043584856264802665759103;
    x[3] = 0.861136311594052575223946488893;

    w[0] = 0.347854845137453857373063949222;
    w[1] = 0.652145154862546142626936050778;
    w[2] = 0.652145154862546142626936050778;
    w[3] = 0.347854845137453857373063949222;
  }
  else if ( n == 5 )
  {
    x[0] = -0.906179845938663992797626878299;
    x[1] = -0.538469310105683091036314420700;
    x[2] = 0.000000000000000000000000000000;
    x[3] = 0.538469310105683091036314420700;
    x[4] = 0.906179845938663992797626878299;

    w[0] = 0.236926885056189087514264040720;
    w[1] = 0.478628670499366468041291514836;
    w[2] = 0.568888888888888888888888888889;
    w[3] = 0.478628670499366468041291514836;
    w[4] = 0.236926885056189087514264040720;
  }
  else if ( n == 6 )
  {
    x[0] = -0.932469514203152027812301554494;
    x[1] = -0.661209386466264513661399595020;
    x[2] = -0.238619186083196908630501721681;
    x[3] = 0.238619186083196908630501721681;
    x[4] = 0.661209386466264513661399595020;
    x[5] = 0.932469514203152027812301554494;

    w[0] = 0.171324492379170345040296142173;
    w[1] = 0.360761573048138607569833513838;
    w[2] = 0.467913934572691047389870343990;
    w[3] = 0.467913934572691047389870343990;
    w[4] = 0.360761573048138607569833513838;
    w[5] = 0.171324492379170345040296142173;
  }
  else if ( n == 7 )
  {
    x[0] = -0.949107912342758524526189684048;
    x[1] = -0.741531185599394439863864773281;
    x[2] = -0.405845151377397166906606412077;
    x[3] = 0.000000000000000000000000000000;
    x[4] = 0.405845151377397166906606412077;
    x[5] = 0.741531185599394439863864773281;
    x[6] = 0.949107912342758524526189684048;

    w[0] = 0.129484966168869693270611432679;
    w[1] = 0.279705391489276667901467771424;
    w[2] = 0.381830050505118944950369775489;
    w[3] = 0.417959183673469387755102040816;
    w[4] = 0.381830050505118944950369775489;
    w[5] = 0.279705391489276667901467771424;
    w[6] = 0.129484966168869693270611432679;
  }
  else if ( n == 8 )
  {
    x[0] = -0.960289856497536231683560868569;
    x[1] = -0.796666477413626739591553936476;
    x[2] = -0.525532409916328985817739049189;
    x[3] = -0.183434642495649804939476142360;
    x[4] = 0.183434642495649804939476142360;
    x[5] = 0.525532409916328985817739049189;
    x[6] = 0.796666477413626739591553936476;
    x[7] = 0.960289856497536231683560868569;

    w[0] = 0.101228536290376259152531354310;
    w[1] = 0.222381034453374470544355994426;
    w[2] = 0.313706645877887287337962201987;
    w[3] = 0.362683783378361982965150449277;
    w[4] = 0.362683783378361982965150449277;
    w[5] = 0.313706645877887287337962201987;
    w[6] = 0.222381034453374470544355994426;
    w[7] = 0.101228536290376259152531354310;
  }
  else if ( n == 9 )
  {
    x[0] = -0.968160239507626089835576203;
    x[1] = -0.836031107326635794299429788;
    x[2] = -0.613371432700590397308702039;
    x[3] = -0.324253423403808929038538015;
    x[4] = 0.000000000000000000000000000;
    x[5] = 0.324253423403808929038538015;
    x[6] = 0.613371432700590397308702039;
    x[7] = 0.836031107326635794299429788;
    x[8] = 0.968160239507626089835576203;

    w[0] = 0.081274388361574411971892158111;
    w[1] = 0.18064816069485740405847203124;
    w[2] = 0.26061069640293546231874286942;
    w[3] = 0.31234707704000284006863040658;
    w[4] = 0.33023935500125976316452506929;
    w[5] = 0.31234707704000284006863040658;
    w[6] = 0.26061069640293546231874286942;
    w[7] = 0.18064816069485740405847203124;
    w[8] = 0.081274388361574411971892158111;
  }
  else if ( n == 10 )
  {
    x[0] = -0.973906528517171720077964012;
    x[1] = -0.865063366688984510732096688;
    x[2] = -0.679409568299024406234327365;
    x[3] = -0.433395394129247190799265943;
    x[4] = -0.148874338981631210884826001;
    x[5] = 0.148874338981631210884826001;
    x[6] = 0.433395394129247190799265943;
    x[7] = 0.679409568299024406234327365;
    x[8] = 0.865063366688984510732096688;
    x[9] = 0.973906528517171720077964012;

    w[0] = 0.066671344308688137593568809893;
    w[1] = 0.14945134915058059314577633966;
    w[2] = 0.21908636251598204399553493423;
    w[3] = 0.26926671930999635509122692157;
    w[4] = 0.29552422471475287017389299465;
    w[5] = 0.29552422471475287017389299465;
    w[6] = 0.26926671930999635509122692157;
    w[7] = 0.21908636251598204399553493423;
    w[8] = 0.14945134915058059314577633966;
    w[9] = 0.066671344308688137593568809893;
  }
  else if ( n == 11 )
  {
    x[0] = -0.978228658146056992803938001;
    x[1] = -0.887062599768095299075157769;
    x[2] = -0.730152005574049324093416252;
    x[3] = -0.519096129206811815925725669;
    x[4] = -0.269543155952344972331531985;
    x[5] = 0.000000000000000000000000000;
    x[6] = 0.269543155952344972331531985;
    x[7] = 0.519096129206811815925725669;
    x[8] = 0.730152005574049324093416252;
    x[9] = 0.887062599768095299075157769;
    x[10] = 0.978228658146056992803938001;

    w[0] = 0.055668567116173666482753720443;
    w[1] = 0.12558036946490462463469429922;
    w[2] = 0.18629021092773425142609764143;
    w[3] = 0.23319376459199047991852370484;
    w[4] = 0.26280454451024666218068886989;
    w[5] = 0.27292508677790063071448352834;
    w[6] = 0.26280454451024666218068886989;
    w[7] = 0.23319376459199047991852370484;
    w[8] = 0.18629021092773425142609764143;
    w[9] = 0.12558036946490462463469429922;
    w[10] = 0.055668567116173666482753720443;
  }
  else if ( n == 12 )
  {
    x[0] = -0.981560634246719250690549090;
    x[1] = -0.904117256370474856678465866;
    x[2] = -0.769902674194304687036893833;
    x[3] = -0.587317954286617447296702419;
    x[4] = -0.367831498998180193752691537;
    x[5] = -0.125233408511468915472441369;
    x[6] = 0.125233408511468915472441369;
    x[7] = 0.367831498998180193752691537;
    x[8] = 0.587317954286617447296702419;
    x[9] = 0.769902674194304687036893833;
    x[10] = 0.904117256370474856678465866;
    x[11] = 0.981560634246719250690549090;

    w[0] = 0.047175336386511827194615961485;
    w[1] = 0.10693932599531843096025471819;
    w[2] = 0.16007832854334622633465252954;
    w[3] = 0.20316742672306592174906445581;
    w[4] = 0.23349253653835480876084989892;
    w[5] = 0.24914704581340278500056243604;
    w[6] = 0.24914704581340278500056243604;
    w[7] = 0.23349253653835480876084989892;
    w[8] = 0.20316742672306592174906445581;
    w[9] = 0.16007832854334622633465252954;
    w[10] = 0.10693932599531843096025471819;
    w[11] = 0.047175336386511827194615961485;
  }
  else if ( n == 13 )
  {
    x[0] = -0.984183054718588149472829449;
    x[1] = -0.917598399222977965206547837;
    x[2] = -0.801578090733309912794206490;
    x[3] = -0.642349339440340220643984607;
    x[4] = -0.448492751036446852877912852;
    x[5] = -0.230458315955134794065528121;
    x[6] = 0.000000000000000000000000000;
    x[7] = 0.230458315955134794065528121;
    x[8] = 0.448492751036446852877912852;
    x[9] = 0.642349339440340220643984607;
    x[10] = 0.80157809073330991279420649;
    x[11] = 0.91759839922297796520654784;
    x[12] = 0.98418305471858814947282945;

    w[0] = 0.040484004765315879520021592201;
    w[1] = 0.092121499837728447914421775954;
    w[2] = 0.13887351021978723846360177687;
    w[3] = 0.17814598076194573828004669200;
    w[4] = 0.20781604753688850231252321931;
    w[5] = 0.22628318026289723841209018604;
    w[6] = 0.23255155323087391019458951527;
    w[7] = 0.22628318026289723841209018604;
    w[8] = 0.20781604753688850231252321931;
    w[9] = 0.17814598076194573828004669200;
    w[10] = 0.13887351021978723846360177687;
    w[11] = 0.092121499837728447914421775954;
    w[12] = 0.040484004765315879520021592201;
  }
  else if ( n == 14 )
  {
    x[0] = -0.986283808696812338841597267;
    x[1] = -0.928434883663573517336391139;
    x[2] = -0.827201315069764993189794743;
    x[3] = -0.687292904811685470148019803;
    x[4] = -0.515248636358154091965290719;
    x[5] = -0.319112368927889760435671824;
    x[6] = -0.108054948707343662066244650;
    x[7] = 0.108054948707343662066244650;
    x[8] = 0.31911236892788976043567182;
    x[9] = 0.51524863635815409196529072;
    x[10] = 0.68729290481168547014801980;
    x[11] = 0.82720131506976499318979474;
    x[12] = 0.92843488366357351733639114;
    x[13] = 0.98628380869681233884159727;

    w[0] = 0.035119460331751863031832876138;
    w[1] = 0.08015808715976020980563327706;
    w[2] = 0.12151857068790318468941480907;
    w[3] = 0.15720316715819353456960193862;
    w[4] = 0.18553839747793781374171659013;
    w[5] = 0.20519846372129560396592406566;
    w[6] = 0.21526385346315779019587644332;
    w[7] = 0.21526385346315779019587644332;
    w[8] = 0.20519846372129560396592406566;
    w[9] = 0.18553839747793781374171659013;
    w[10] = 0.15720316715819353456960193862;
    w[11] = 0.12151857068790318468941480907;
    w[12] = 0.08015808715976020980563327706;
    w[13] = 0.035119460331751863031832876138;
  }
  else if ( n == 15 )
  {
    x[0] = -0.987992518020485428489565719;
    x[1] = -0.937273392400705904307758948;
    x[2] = -0.848206583410427216200648321;
    x[3] = -0.724417731360170047416186055;
    x[4] = -0.570972172608538847537226737;
    x[5] = -0.394151347077563369897207371;
    x[6] = -0.201194093997434522300628303;
    x[7] = 0.00000000000000000000000000;
    x[8] = 0.20119409399743452230062830;
    x[9] = 0.39415134707756336989720737;
    x[10] = 0.57097217260853884753722674;
    x[11] = 0.72441773136017004741618605;
    x[12] = 0.84820658341042721620064832;
    x[13] = 0.93727339240070590430775895;
    x[14] = 0.98799251802048542848956572;

    w[0] = 0.030753241996117268354628393577;
    w[1] = 0.070366047488108124709267416451;
    w[2] = 0.107159220467171935011869546686;
    w[3] = 0.13957067792615431444780479451;
    w[4] = 0.16626920581699393355320086048;
    w[5] = 0.18616100001556221102680056187;
    w[6] = 0.19843148532711157645611832644;
    w[7] = 0.20257824192556127288062019997;
    w[8] = 0.19843148532711157645611832644;
    w[9] = 0.18616100001556221102680056187;
    w[10] = 0.16626920581699393355320086048;
    w[11] = 0.13957067792615431444780479451;
    w[12] = 0.107159220467171935011869546686;
    w[13] = 0.070366047488108124709267416451;
    w[14] = 0.030753241996117268354628393577;
  }
  else if ( n == 16 )
  {
    x[0] = -0.989400934991649932596154173;
    x[1] = -0.944575023073232576077988416;
    x[2] = -0.865631202387831743880467898;
    x[3] = -0.755404408355003033895101195;
    x[4] = -0.617876244402643748446671764;
    x[5] = -0.458016777657227386342419443;
    x[6] = -0.281603550779258913230460501;
    x[7] = -0.09501250983763744018531934;
    x[8] = 0.09501250983763744018531934;
    x[9] = 0.28160355077925891323046050;
    x[10] = 0.45801677765722738634241944;
    x[11] = 0.61787624440264374844667176;
    x[12] = 0.75540440835500303389510119;
    x[13] = 0.86563120238783174388046790;
    x[14] = 0.94457502307323257607798842;
    x[15] = 0.98940093499164993259615417;

    w[0] = 0.027152459411754094851780572456;
    w[1] = 0.062253523938647892862843836994;
    w[2] = 0.09515851168249278480992510760;
    w[3] = 0.12462897125553387205247628219;
    w[4] = 0.14959598881657673208150173055;
    w[5] = 0.16915651939500253818931207903;
    w[6] = 0.18260341504492358886676366797;
    w[7] = 0.18945061045506849628539672321;
    w[8] = 0.18945061045506849628539672321;
    w[9] = 0.18260341504492358886676366797;
    w[10] = 0.16915651939500253818931207903;
    w[11] = 0.14959598881657673208150173055;
    w[12] = 0.12462897125553387205247628219;
    w[13] = 0.09515851168249278480992510760;
    w[14] = 0.062253523938647892862843836994;
    w[15] = 0.027152459411754094851780572456;
  }
  else if ( n == 17 )
  {
    x[0] = -0.990575475314417335675434020;
    x[1] = -0.950675521768767761222716958;
    x[2] = -0.880239153726985902122955694;
    x[3] = -0.781514003896801406925230056;
    x[4] = -0.657671159216690765850302217;
    x[5] = -0.512690537086476967886246569;
    x[6] = -0.35123176345387631529718552;
    x[7] = -0.17848418149584785585067749;
    x[8] = 0.00000000000000000000000000;
    x[9] = 0.17848418149584785585067749;
    x[10] = 0.35123176345387631529718552;
    x[11] = 0.51269053708647696788624657;
    x[12] = 0.65767115921669076585030222;
    x[13] = 0.78151400389680140692523006;
    x[14] = 0.88023915372698590212295569;
    x[15] = 0.95067552176876776122271696;
    x[16] = 0.99057547531441733567543402;

    w[0] = 0.024148302868547931960110026288;
    w[1] = 0.055459529373987201129440165359;
    w[2] = 0.085036148317179180883535370191;
    w[3] = 0.111883847193403971094788385626;
    w[4] = 0.13513636846852547328631998170;
    w[5] = 0.15404576107681028808143159480;
    w[6] = 0.16800410215645004450997066379;
    w[7] = 0.17656270536699264632527099011;
    w[8] = 0.17944647035620652545826564426;
    w[9] = 0.17656270536699264632527099011;
    w[10] = 0.16800410215645004450997066379;
    w[11] = 0.15404576107681028808143159480;
    w[12] = 0.13513636846852547328631998170;
    w[13] = 0.111883847193403971094788385626;
    w[14] = 0.085036148317179180883535370191;
    w[15] = 0.055459529373987201129440165359;
    w[16] = 0.024148302868547931960110026288;
  }
  else if ( n == 18 )
  {
    x[0] = -0.991565168420930946730016005;
    x[1] = -0.955823949571397755181195893;
    x[2] = -0.892602466497555739206060591;
    x[3] = -0.803704958972523115682417455;
    x[4] = -0.691687043060353207874891081;
    x[5] = -0.55977083107394753460787155;
    x[6] = -0.41175116146284264603593179;
    x[7] = -0.25188622569150550958897285;
    x[8] = -0.08477501304173530124226185;
    x[9] = 0.08477501304173530124226185;
    x[10] = 0.25188622569150550958897285;
    x[11] = 0.41175116146284264603593179;
    x[12] = 0.55977083107394753460787155;
    x[13] = 0.69168704306035320787489108;
    x[14] = 0.80370495897252311568241746;
    x[15] = 0.89260246649755573920606059;
    x[16] = 0.95582394957139775518119589;
    x[17] = 0.99156516842093094673001600;

    w[0] = 0.021616013526483310313342710266;
    w[1] = 0.049714548894969796453334946203;
    w[2] = 0.07642573025488905652912967762;
    w[3] = 0.10094204410628716556281398492;
    w[4] = 0.12255520671147846018451912680;
    w[5] = 0.14064291467065065120473130375;
    w[6] = 0.15468467512626524492541800384;
    w[7] = 0.16427648374583272298605377647;
    w[8] = 0.16914238296314359184065647013;
    w[9] = 0.16914238296314359184065647013;
    w[10] = 0.16427648374583272298605377647;
    w[11] = 0.15468467512626524492541800384;
    w[12] = 0.14064291467065065120473130375;
    w[13] = 0.12255520671147846018451912680;
    w[14] = 0.10094204410628716556281398492;
    w[15] = 0.07642573025488905652912967762;
    w[16] = 0.049714548894969796453334946203;
    w[17] = 0.021616013526483310313342710266;
  }
  else if ( n == 19 )
  {
    x[0] = -0.992406843843584403189017670;
    x[1] = -0.960208152134830030852778841;
    x[2] = -0.903155903614817901642660929;
    x[3] = -0.822714656537142824978922487;
    x[4] = -0.72096617733522937861709586;
    x[5] = -0.60054530466168102346963816;
    x[6] = -0.46457074137596094571726715;
    x[7] = -0.31656409996362983199011733;
    x[8] = -0.16035864564022537586809612;
    x[9] = 0.00000000000000000000000000;
    x[10] = 0.16035864564022537586809612;
    x[11] = 0.31656409996362983199011733;
    x[12] = 0.46457074137596094571726715;
    x[13] = 0.60054530466168102346963816;
    x[14] = 0.72096617733522937861709586;
    x[15] = 0.82271465653714282497892249;
    x[16] = 0.90315590361481790164266093;
    x[17] = 0.96020815213483003085277884;
    x[18] = 0.99240684384358440318901767;

    w[0] = 0.019461788229726477036312041464;
    w[1] = 0.044814226765699600332838157402;
    w[2] = 0.069044542737641226580708258006;
    w[3] = 0.091490021622449999464462094124;
    w[4] = 0.111566645547333994716023901682;
    w[5] = 0.12875396253933622767551578486;
    w[6] = 0.14260670217360661177574610944;
    w[7] = 0.15276604206585966677885540090;
    w[8] = 0.15896884339395434764995643946;
    w[9] = 0.16105444984878369597916362532;
    w[10] = 0.15896884339395434764995643946;
    w[11] = 0.15276604206585966677885540090;
    w[12] = 0.14260670217360661177574610944;
    w[13] = 0.12875396253933622767551578486;
    w[14] = 0.111566645547333994716023901682;
    w[15] = 0.091490021622449999464462094124;
    w[16] = 0.069044542737641226580708258006;
    w[17] = 0.044814226765699600332838157402;
    w[18] = 0.019461788229726477036312041464;
  }
  else if ( n == 20 )
  {
    x[0] = -0.993128599185094924786122388;
    x[1] = -0.963971927277913791267666131;
    x[2] = -0.912234428251325905867752441;
    x[3] = -0.83911697182221882339452906;
    x[4] = -0.74633190646015079261430507;
    x[5] = -0.63605368072651502545283670;
    x[6] = -0.51086700195082709800436405;
    x[7] = -0.37370608871541956067254818;
    x[8] = -0.22778585114164507808049620;
    x[9] = -0.07652652113349733375464041;
    x[10] = 0.07652652113349733375464041;
    x[11] = 0.22778585114164507808049620;
    x[12] = 0.37370608871541956067254818;
    x[13] = 0.51086700195082709800436405;
    x[14] = 0.63605368072651502545283670;
    x[15] = 0.74633190646015079261430507;
    x[16] = 0.83911697182221882339452906;
    x[17] = 0.91223442825132590586775244;
    x[18] = 0.96397192727791379126766613;
    x[19] = 0.99312859918509492478612239;

    w[0] = 0.017614007139152118311861962352;
    w[1] = 0.040601429800386941331039952275;
    w[2] = 0.062672048334109063569506535187;
    w[3] = 0.08327674157670474872475814322;
    w[4] = 0.10193011981724043503675013548;
    w[5] = 0.11819453196151841731237737771;
    w[6] = 0.13168863844917662689849449975;
    w[7] = 0.14209610931838205132929832507;
    w[8] = 0.14917298647260374678782873700;
    w[9] = 0.15275338713072585069808433195;
    w[10] = 0.15275338713072585069808433195;
    w[11] = 0.14917298647260374678782873700;
    w[12] = 0.14209610931838205132929832507;
    w[13] = 0.13168863844917662689849449975;
    w[14] = 0.11819453196151841731237737771;
    w[15] = 0.10193011981724043503675013548;
    w[16] = 0.08327674157670474872475814322;
    w[17] = 0.062672048334109063569506535187;
    w[18] = 0.040601429800386941331039952275;
    w[19] = 0.017614007139152118311861962352;
  }
  else if ( n == 21 )
  {
    x[ 0] =  -0.99375217062038950026024204;
    x[ 1] =  -0.96722683856630629431662221;
    x[ 2] =  -0.92009933415040082879018713;
    x[ 3] =  -0.85336336458331728364725064;
    x[ 4] =  -0.76843996347567790861587785;
    x[ 5] =  -0.66713880419741231930596667;
    x[ 6] =  -0.55161883588721980705901880;
    x[ 7] =  -0.42434212020743878357366889;
    x[ 8] =  -0.28802131680240109660079252;
    x[9] =  -0.14556185416089509093703098;
    x[10] =   0.00000000000000000000000000;
    x[11] =  +0.14556185416089509093703098;
    x[12] =  +0.28802131680240109660079252;
    x[13] =  +0.42434212020743878357366889;
    x[14] =  +0.55161883588721980705901880;
    x[15] =  +0.66713880419741231930596667;
    x[16] =  +0.76843996347567790861587785;
    x[17] =  +0.85336336458331728364725064;
    x[18] =  +0.92009933415040082879018713;
    x[19] =  +0.96722683856630629431662221;
    x[20] =  +0.99375217062038950026024204;

    w[ 0] =   0.016017228257774333324224616858;
    w[ 1] =   0.036953789770852493799950668299; 
    w[ 2] =   0.057134425426857208283635826472;
    w[ 3] =   0.076100113628379302017051653300;
    w[ 4] =   0.093444423456033861553289741114;
    w[ 5] =   0.108797299167148377663474578070;
    w[ 6] =   0.12183141605372853419536717713;
    w[ 7] =   0.13226893863333746178105257450;
    w[ 8] =   0.13988739479107315472213342387;
    w[9] =   0.14452440398997005906382716655;
    w[10] =   0.14608113364969042719198514768;
    w[11] =   0.14452440398997005906382716655; 
    w[12] =   0.13988739479107315472213342387; 
    w[13] =   0.13226893863333746178105257450;
    w[14] =   0.12183141605372853419536717713;
    w[15] =   0.108797299167148377663474578070;
    w[16] =   0.093444423456033861553289741114;
    w[17] =   0.076100113628379302017051653300;
    w[18] =   0.057134425426857208283635826472;
    w[19] =   0.036953789770852493799950668299;
    w[20] =   0.016017228257774333324224616858;
  }
  else if ( n == 22 )
  {
    x[0] = -0.99429458548239929207303142;
    x[1] = -0.97006049783542872712395099;
    x[2] = -0.92695677218717400052069294;
    x[3] = -0.86581257772030013653642564;
    x[4] = -0.78781680597920816200427796;
    x[5] = -0.69448726318668278005068984;
    x[6] = -0.58764040350691159295887693;
    x[7] = -0.46935583798675702640633071;
    x[8] = -0.34193582089208422515814742;
    x[9] = -0.20786042668822128547884653;
    x[10] = -0.06973927331972222121384180;
    x[11] = 0.06973927331972222121384180;
    x[12] = 0.20786042668822128547884653;
    x[13] = 0.34193582089208422515814742;
    x[14] = 0.46935583798675702640633071;
    x[15] = 0.58764040350691159295887693;
    x[16] = 0.69448726318668278005068984;
    x[17] = 0.78781680597920816200427796;
    x[18] = 0.86581257772030013653642564;
    x[19] = 0.92695677218717400052069294;
    x[20] = 0.97006049783542872712395099;
    x[21] = 0.99429458548239929207303142;
 
    w[0] = 0.014627995298272200684991098047;
    w[1] = 0.033774901584814154793302246866;
    w[2] = 0.052293335152683285940312051273;
    w[3] = 0.06979646842452048809496141893;
    w[4] = 0.08594160621706772741444368137;
    w[5] = 0.10041414444288096493207883783;
    w[6] = 0.11293229608053921839340060742;
    w[7] = 0.12325237681051242428556098615;
    w[8] = 0.13117350478706237073296499253;
    w[9] = 0.13654149834601517135257383123;
    w[10] = 0.13925187285563199337541024834;
    w[11] = 0.13925187285563199337541024834;
    w[12] = 0.13654149834601517135257383123;
    w[13] = 0.13117350478706237073296499253;
    w[14] = 0.12325237681051242428556098615;
    w[15] = 0.11293229608053921839340060742;
    w[16] = 0.10041414444288096493207883783;
    w[17] = 0.08594160621706772741444368137;
    w[18] = 0.06979646842452048809496141893;
    w[19] = 0.052293335152683285940312051273;
    w[20] = 0.033774901584814154793302246866;
    w[21] = 0.014627995298272200684991098047;
  }
  else if ( n == 23 )
  {
    x[0] = -0.99476933499755212352392572;
    x[1] = -0.97254247121811523195602408;
    x[2] = -0.93297108682601610234919699;
    x[3] = -0.87675235827044166737815689;
    x[4] = -0.80488840161883989215111841;
    x[5] = -0.71866136313195019446162448;
    x[6] = -0.61960987576364615638509731;
    x[7] = -0.50950147784600754968979305;
    x[8] = -0.39030103803029083142148887;
    x[9] = -0.26413568097034493053386954;
    x[10] = -0.13325682429846611093174268;
    x[11] = 0.00000000000000000000000000;
    x[12] = 0.13325682429846611093174268;
    x[13] = 0.26413568097034493053386954;
    x[14] = 0.39030103803029083142148887;
    x[15] = 0.50950147784600754968979305;
    x[16] = 0.61960987576364615638509731;
    x[17] = 0.71866136313195019446162448;
    x[18] = 0.80488840161883989215111841;
    x[19] = 0.87675235827044166737815689;
    x[20] = 0.93297108682601610234919699;
    x[21] = 0.97254247121811523195602408;
    x[22] = 0.99476933499755212352392572;

    w[0] = 0.013411859487141772081309493459;
    w[1] = 0.030988005856979444310694219642;
    w[2] = 0.048037671731084668571641071632;
    w[3] = 0.064232421408525852127169615159;
    w[4] = 0.079281411776718954922892524742;
    w[5] = 0.092915766060035147477018617370;
    w[6] = 0.104892091464541410074086185015;
    w[7] = 0.11499664022241136494164351293;
    w[8] = 0.12304908430672953046757840067;
    w[9] = 0.12890572218808214997859533940;
    w[10] = 0.13246203940469661737164246470;
    w[11] = 0.13365457218610617535145711055;
    w[12] = 0.13246203940469661737164246470;
    w[13] = 0.12890572218808214997859533940;
    w[14] = 0.12304908430672953046757840067;
    w[15] = 0.11499664022241136494164351293;
    w[16] = 0.104892091464541410074086185015;
    w[17] = 0.092915766060035147477018617370;
    w[18] = 0.079281411776718954922892524742;
    w[19] = 0.064232421408525852127169615159;
    w[20] = 0.048037671731084668571641071632;
    w[21] = 0.030988005856979444310694219642;
    w[22] = 0.013411859487141772081309493459;
  }
  else if ( n == 24 )
  {
    x[0] = -0.99518721999702136017999741;
    x[1] = -0.97472855597130949819839199;
    x[2] = -0.93827455200273275852364900;
    x[3] = -0.88641552700440103421315434;
    x[4] = -0.82000198597390292195394987;
    x[5] = -0.74012419157855436424382810;
    x[6] = -0.64809365193697556925249579;
    x[7] = -0.54542147138883953565837562;
    x[8] = -0.43379350762604513848708423;
    x[9] = -0.31504267969616337438679329;
    x[10] = -0.19111886747361630915863982;
    x[11] = -0.06405689286260562608504308;
    x[12] = 0.06405689286260562608504308;
    x[13] = 0.19111886747361630915863982;
    x[14] = 0.31504267969616337438679329;
    x[15] = 0.43379350762604513848708423;
    x[16] = 0.54542147138883953565837562;
    x[17] = 0.64809365193697556925249579;
    x[18] = 0.74012419157855436424382810;
    x[19] = 0.82000198597390292195394987;
    x[20] = 0.88641552700440103421315434;
    x[21] = 0.93827455200273275852364900;
    x[22] = 0.97472855597130949819839199;
    x[23] = 0.99518721999702136017999741;

    w[0] = 0.012341229799987199546805667070;
    w[1] = 0.028531388628933663181307815952;
    w[2] = 0.044277438817419806168602748211;
    w[3] = 0.059298584915436780746367758500;
    w[4] = 0.07334648141108030573403361525;
    w[5] = 0.08619016153195327591718520298;
    w[6] = 0.09761865210411388826988066446;
    w[7] = 0.10744427011596563478257734245;
    w[8] = 0.11550566805372560135334448391;
    w[9] = 0.12167047292780339120446315348;
    w[10] = 0.12583745634682829612137538251;
    w[11] = 0.12793819534675215697405616522;
    w[12] = 0.12793819534675215697405616522;
    w[13] = 0.12583745634682829612137538251;
    w[14] = 0.12167047292780339120446315348;
    w[15] = 0.11550566805372560135334448391;
    w[16] = 0.10744427011596563478257734245;
    w[17] = 0.09761865210411388826988066446;
    w[18] = 0.08619016153195327591718520298;
    w[19] = 0.07334648141108030573403361525;
    w[20] = 0.059298584915436780746367758500;
    w[21] = 0.044277438817419806168602748211;
    w[22] = 0.028531388628933663181307815952;
    w[23] = 0.012341229799987199546805667070;
  }
  else if ( n == 25 )
  {
    x[0] = -0.99555696979049809790878495;
    x[1] = -0.97666392145951751149831539;
    x[2] = -0.94297457122897433941401117;
    x[3] = -0.89499199787827536885104201;
    x[4] = -0.83344262876083400142102111;
    x[5] = -0.75925926303735763057728287;
    x[6] = -0.67356636847346836448512063;
    x[7] = -0.57766293024122296772368984;
    x[8] = -0.47300273144571496052218212;
    x[9] = -0.36117230580938783773582173;
    x[10] = -0.24386688372098843204519036;
    x[11] = -0.12286469261071039638735982;
    x[12] = 0.00000000000000000000000000;
    x[13] = 0.12286469261071039638735982;
    x[14] = 0.24386688372098843204519036;
    x[15] = 0.36117230580938783773582173;
    x[16] = 0.47300273144571496052218212;
    x[17] = 0.57766293024122296772368984;
    x[18] = 0.67356636847346836448512063;
    x[19] = 0.75925926303735763057728287;
    x[20] = 0.83344262876083400142102111;
    x[21] = 0.89499199787827536885104201;
    x[22] = 0.94297457122897433941401117;
    x[23] = 0.97666392145951751149831539;
    x[24] = 0.99555696979049809790878495;

    w[0] = 0.0113937985010262879479029641132;
    w[1] = 0.026354986615032137261901815295;
    w[2] = 0.040939156701306312655623487712;
    w[3] = 0.054904695975835191925936891541;
    w[4] = 0.068038333812356917207187185657;
    w[5] = 0.080140700335001018013234959669;
    w[6] = 0.091028261982963649811497220703;
    w[7] = 0.100535949067050644202206890393;
    w[8] = 0.108519624474263653116093957050;
    w[9] = 0.11485825914571164833932554587;
    w[10] = 0.11945576353578477222817812651;
    w[11] = 0.12224244299031004168895951895;
    w[12] = 0.12317605372671545120390287308;
    w[13] = 0.12224244299031004168895951895;
    w[14] = 0.11945576353578477222817812651;
    w[15] = 0.11485825914571164833932554587;
    w[16] = 0.108519624474263653116093957050;
    w[17] = 0.100535949067050644202206890393;
    w[18] = 0.091028261982963649811497220703;
    w[19] = 0.080140700335001018013234959669;
    w[20] = 0.068038333812356917207187185657;
    w[21] = 0.054904695975835191925936891541;
    w[22] = 0.040939156701306312655623487712;
    w[23] = 0.026354986615032137261901815295;
    w[24] = 0.0113937985010262879479029641132;
  }
  else if ( n == 26 )
  {
    x[0] = -0.99588570114561692900321696;
    x[1] = -0.97838544595647099110058035;
    x[2] = -0.94715906666171425013591528;
    x[3] = -0.90263786198430707421766560;
    x[4] = -0.84544594278849801879750706;
    x[5] = -0.77638594882067885619296725;
    x[6] = -0.69642726041995726486381391;
    x[7] = -0.60669229301761806323197875;
    x[8] = -0.50844071482450571769570306;
    x[9] = -0.40305175512348630648107738;
    x[10] = -0.29200483948595689514283538;
    x[11] = -0.17685882035689018396905775;
    x[12] = -0.05923009342931320709371858;
    x[13] = 0.05923009342931320709371858;
    x[14] = 0.17685882035689018396905775;
    x[15] = 0.29200483948595689514283538;
    x[16] = 0.40305175512348630648107738;
    x[17] = 0.50844071482450571769570306;
    x[18] = 0.60669229301761806323197875;
    x[19] = 0.69642726041995726486381391;
    x[20] = 0.77638594882067885619296725;
    x[21] = 0.84544594278849801879750706;
    x[22] = 0.90263786198430707421766560;
    x[23] = 0.94715906666171425013591528;
    x[24] = 0.97838544595647099110058035;
    x[25] = 0.99588570114561692900321696;

    w[0] = 0.010551372617343007155651187685;
    w[1] = 0.024417851092631908789615827520;
    w[2] = 0.037962383294362763950303141249;
    w[3] = 0.050975825297147811998319900724;
    w[4] = 0.063274046329574835539453689907;
    w[5] = 0.07468414976565974588707579610;
    w[6] = 0.08504589431348523921044776508;
    w[7] = 0.09421380035591414846366488307;
    w[8] = 0.10205916109442542323841407025;
    w[9] = 0.10847184052857659065657942673;
    w[10] = 0.11336181654631966654944071844;
    w[11] = 0.11666044348529658204466250754;
    w[12] = 0.11832141527926227651637108570;
    w[13] = 0.11832141527926227651637108570;
    w[14] = 0.11666044348529658204466250754;
    w[15] = 0.11336181654631966654944071844;
    w[16] = 0.10847184052857659065657942673;
    w[17] = 0.10205916109442542323841407025;
    w[18] = 0.09421380035591414846366488307;
    w[19] = 0.08504589431348523921044776508;
    w[20] = 0.07468414976565974588707579610;
    w[21] = 0.063274046329574835539453689907;
    w[22] = 0.050975825297147811998319900724;
    w[23] = 0.037962383294362763950303141249;
    w[24] = 0.024417851092631908789615827520;
    w[25] = 0.010551372617343007155651187685;
  }
  else if ( n == 27 )
  {
    x[0] = -0.99617926288898856693888721;
    x[1] = -0.97992347596150122285587336;
    x[2] = -0.95090055781470500685190803;
    x[3] = -0.90948232067749110430064502;
    x[4] = -0.85620790801829449030273722;
    x[5] = -0.79177163907050822714439734;
    x[6] = -0.71701347373942369929481621;
    x[7] = -0.63290797194649514092773464;
    x[8] = -0.54055156457945689490030094;
    x[9] = -0.44114825175002688058597416;
    x[10] = -0.33599390363850889973031903;
    x[11] = -0.22645936543953685885723911;
    x[12] = -0.11397258560952996693289498;
    x[13] = 0.00000000000000000000000000;
    x[14] = 0.11397258560952996693289498;
    x[15] = 0.22645936543953685885723911;
    x[16] = 0.33599390363850889973031903;
    x[17] = 0.44114825175002688058597416;
    x[18] = 0.54055156457945689490030094;
    x[19] = 0.63290797194649514092773464;
    x[20] = 0.71701347373942369929481621;
    x[21] = 0.79177163907050822714439734;
    x[22] = 0.85620790801829449030273722;
    x[23] = 0.90948232067749110430064502;
    x[24] = 0.95090055781470500685190803;
    x[25] = 0.97992347596150122285587336;
    x[26] = 0.99617926288898856693888721;

    w[0] = 0.0097989960512943602611500550912;
    w[1] = 0.022686231596180623196034206447;
    w[2] = 0.035297053757419711022578289305;
    w[3] = 0.047449412520615062704096710114;
    w[4] = 0.058983536859833599110300833720;
    w[5] = 0.069748823766245592984322888357;
    w[6] = 0.079604867773057771263074959010;
    w[7] = 0.088423158543756950194322802854;
    w[8] = 0.096088727370028507565652646558;
    w[9] = 0.102501637817745798671247711533;
    w[10] = 0.107578285788533187212162984427;
    w[11] = 0.111252488356845192672163096043;
    w[12] = 0.113476346108965148620369948092;
    w[13] = 0.11422086737895698904504573690;
    w[14] = 0.113476346108965148620369948092;
    w[15] = 0.111252488356845192672163096043;
    w[16] = 0.107578285788533187212162984427;
    w[17] = 0.102501637817745798671247711533;
    w[18] = 0.096088727370028507565652646558;
    w[19] = 0.088423158543756950194322802854;
    w[20] = 0.079604867773057771263074959010;
    w[21] = 0.069748823766245592984322888357;
    w[22] = 0.058983536859833599110300833720;
    w[23] = 0.047449412520615062704096710114;
    w[24] = 0.035297053757419711022578289305;
    w[25] = 0.022686231596180623196034206447;
    w[26] = 0.0097989960512943602611500550912;
  }
  else if ( n == 28 )
  {
    x[0] = -0.99644249757395444995043639;
    x[1] = -0.98130316537087275369455995;
    x[2] = -0.95425928062893819725410184;
    x[3] = -0.91563302639213207386968942;
    x[4] = -0.86589252257439504894225457;
    x[5] = -0.80564137091717917144788596;
    x[6] = -0.73561087801363177202814451;
    x[7] = -0.65665109403886496121989818;
    x[8] = -0.56972047181140171930800328;
    x[9] = -0.47587422495511826103441185;
    x[10] = -0.37625151608907871022135721;
    x[11] = -0.27206162763517807767682636;
    x[12] = -0.16456928213338077128147178;
    x[13] = -0.05507928988403427042651653;
    x[14] = 0.05507928988403427042651653;
    x[15] = 0.16456928213338077128147178;
    x[16] = 0.27206162763517807767682636;
    x[17] = 0.37625151608907871022135721;
    x[18] = 0.47587422495511826103441185;
    x[19] = 0.56972047181140171930800328;
    x[20] = 0.65665109403886496121989818;
    x[21] = 0.73561087801363177202814451;
    x[22] = 0.80564137091717917144788596;
    x[23] = 0.86589252257439504894225457;
    x[24] = 0.91563302639213207386968942;
    x[25] = 0.95425928062893819725410184;
    x[26] = 0.98130316537087275369455995;
    x[27] = 0.99644249757395444995043639;

    w[0] = 0.009124282593094517738816153923;
    w[1] = 0.021132112592771259751500380993;
    w[2] = 0.032901427782304379977630819171;
    w[3] = 0.044272934759004227839587877653;
    w[4] = 0.055107345675716745431482918227;
    w[5] = 0.06527292396699959579339756678;
    w[6] = 0.07464621423456877902393188717;
    w[7] = 0.08311341722890121839039649824;
    w[8] = 0.09057174439303284094218603134;
    w[9] = 0.09693065799792991585048900610;
    w[10] = 0.10211296757806076981421663851;
    w[11] = 0.10605576592284641791041643700;
    w[12] = 0.10871119225829413525357151930;
    w[13] = 0.11004701301647519628237626560;
    w[14] = 0.11004701301647519628237626560;
    w[15] = 0.10871119225829413525357151930;
    w[16] = 0.10605576592284641791041643700;
    w[17] = 0.10211296757806076981421663851;
    w[18] = 0.09693065799792991585048900610;
    w[19] = 0.09057174439303284094218603134;
    w[20] = 0.08311341722890121839039649824;
    w[21] = 0.07464621423456877902393188717;
    w[22] = 0.06527292396699959579339756678;
    w[23] = 0.055107345675716745431482918227;
    w[24] = 0.044272934759004227839587877653;
    w[25] = 0.032901427782304379977630819171;
    w[26] = 0.021132112592771259751500380993;
    w[27] = 0.009124282593094517738816153923;
  }
  else if ( n == 29 )
  {
    x[0] = -0.99667944226059658616319153;
    x[1] = -0.98254550526141317487092602;
    x[2] = -0.95728559577808772579820804;
    x[3] = -0.92118023295305878509375344;
    x[4] = -0.87463780492010279041779342;
    x[5] = -0.81818548761525244498957221;
    x[6] = -0.75246285173447713391261008;
    x[7] = -0.67821453760268651515618501;
    x[8] = -0.59628179713822782037958621;
    x[9] = -0.50759295512422764210262792;
    x[10] = -0.41315288817400866389070659;
    x[11] = -0.31403163786763993494819592;
    x[12] = -0.21135228616600107450637573;
    x[13] = -0.10627823013267923017098239;
    x[14] = 0.00000000000000000000000000;
    x[15] = 0.10627823013267923017098239;
    x[16] = 0.21135228616600107450637573;
    x[17] = 0.31403163786763993494819592;
    x[18] = 0.41315288817400866389070659;
    x[19] = 0.50759295512422764210262792;
    x[20] = 0.59628179713822782037958621;
    x[21] = 0.67821453760268651515618501;
    x[22] = 0.75246285173447713391261008;
    x[23] = 0.81818548761525244498957221;
    x[24] = 0.87463780492010279041779342;
    x[25] = 0.92118023295305878509375344;
    x[26] = 0.95728559577808772579820804;
    x[27] = 0.98254550526141317487092602;
    x[28] = 0.99667944226059658616319153;

    w[0] = 0.0085169038787464096542638133022;
    w[1] = 0.019732085056122705983859801640;
    w[2] = 0.030740492202093622644408525375;
    w[3] = 0.041402062518682836104830010114;
    w[4] = 0.051594826902497923912594381180;
    w[5] = 0.061203090657079138542109848024;
    w[6] = 0.070117933255051278569581486949;
    w[7] = 0.078238327135763783828144888660;
    w[8] = 0.085472257366172527545344849297;
    w[9] = 0.091737757139258763347966411077;
    w[10] = 0.096963834094408606301900074883;
    w[11] = 0.101091273759914966121820546907;
    w[12] = 0.104073310077729373913328471285;
    w[13] = 0.105876155097320941406591327852;
    w[14] = 0.10647938171831424424651112691;
    w[15] = 0.105876155097320941406591327852;
    w[16] = 0.104073310077729373913328471285;
    w[17] = 0.101091273759914966121820546907;
    w[18] = 0.096963834094408606301900074883;
    w[19] = 0.091737757139258763347966411077;
    w[20] = 0.085472257366172527545344849297;
    w[21] = 0.078238327135763783828144888660;
    w[22] = 0.070117933255051278569581486949;
    w[23] = 0.061203090657079138542109848024;
    w[24] = 0.051594826902497923912594381180;
    w[25] = 0.041402062518682836104830010114;
    w[26] = 0.030740492202093622644408525375;
    w[27] = 0.019732085056122705983859801640;
    w[28] = 0.0085169038787464096542638133022;
  }
  else if ( n == 30 )
  {
    x[0] = -0.99689348407464954027163005;
    x[1] = -0.98366812327974720997003258;
    x[2] = -0.96002186496830751221687103;
    x[3] = -0.92620004742927432587932428;
    x[4] = -0.88256053579205268154311646;
    x[5] = -0.82956576238276839744289812;
    x[6] = -0.76777743210482619491797734;
    x[7] = -0.69785049479331579693229239;
    x[8] = -0.62052618298924286114047756;
    x[9] = -0.53662414814201989926416979;
    x[10] = -0.44703376953808917678060990;
    x[11] = -0.35270472553087811347103721;
    x[12] = -0.25463692616788984643980513;
    x[13] = -0.15386991360858354696379467;
    x[14] = -0.05147184255531769583302521;
    x[15] = 0.05147184255531769583302521;
    x[16] = 0.15386991360858354696379467;
    x[17] = 0.25463692616788984643980513;
    x[18] = 0.35270472553087811347103721;
    x[19] = 0.44703376953808917678060990;
    x[20] = 0.53662414814201989926416979;
    x[21] = 0.62052618298924286114047756;
    x[22] = 0.69785049479331579693229239;
    x[23] = 0.76777743210482619491797734;
    x[24] = 0.82956576238276839744289812;
    x[25] = 0.88256053579205268154311646;
    x[26] = 0.92620004742927432587932428;
    x[27] = 0.96002186496830751221687103;
    x[28] = 0.98366812327974720997003258;
    x[29] = 0.99689348407464954027163005;

    w[0] = 0.007968192496166605615465883475;
    w[1] = 0.018466468311090959142302131912;
    w[2] = 0.028784707883323369349719179611;
    w[3] = 0.038799192569627049596801936446;
    w[4] = 0.048402672830594052902938140423;
    w[5] = 0.057493156217619066481721689402;
    w[6] = 0.06597422988218049512812851512;
    w[7] = 0.07375597473770520626824385002;
    w[8] = 0.08075589522942021535469493846;
    w[9] = 0.08689978720108297980238753072;
    w[10] = 0.09212252223778612871763270709;
    w[11] = 0.09636873717464425963946862635;
    w[12] = 0.09959342058679526706278028210;
    w[13] = 0.10176238974840550459642895217;
    w[14] = 0.10285265289355884034128563671;
    w[15] = 0.10285265289355884034128563671;
    w[16] = 0.10176238974840550459642895217;
    w[17] = 0.09959342058679526706278028210;
    w[18] = 0.09636873717464425963946862635;
    w[19] = 0.09212252223778612871763270709;
    w[20] = 0.08689978720108297980238753072;
    w[21] = 0.08075589522942021535469493846;
    w[22] = 0.07375597473770520626824385002;
    w[23] = 0.06597422988218049512812851512;
    w[24] = 0.057493156217619066481721689402;
    w[25] = 0.048402672830594052902938140423;
    w[26] = 0.038799192569627049596801936446;
    w[27] = 0.028784707883323369349719179611;
    w[28] = 0.018466468311090959142302131912;
    w[29] = 0.007968192496166605615465883475;
  }
  else if ( n == 31 )
  {
    x[0] = -0.99708748181947707405562655;
    x[1] = -0.98468590966515248400246517;
    x[2] = -0.96250392509294966178905240;
    x[3] = -0.93075699789664816495694576;
    x[4] = -0.88976002994827104337419201;
    x[5] = -0.83992032014626734008690454;
    x[6] = -0.78173314841662494040636002;
    x[7] = -0.71577678458685328390597087;
    x[8] = -0.64270672292426034618441820;
    x[9] = -0.56324916140714926272094492;
    x[10] = -0.47819378204490248044059404;
    x[11] = -0.38838590160823294306135146;
    x[12] = -0.29471806998170161661790390;
    x[13] = -0.19812119933557062877241300;
    x[14] = -0.09955531215234152032517479;
    x[15] = 0.00000000000000000000000000;
    x[16] = 0.09955531215234152032517479;
    x[17] = 0.19812119933557062877241300;
    x[18] = 0.29471806998170161661790390;
    x[19] = 0.38838590160823294306135146;
    x[20] = 0.47819378204490248044059404;
    x[21] = 0.56324916140714926272094492;
    x[22] = 0.64270672292426034618441820;
    x[23] = 0.71577678458685328390597087;
    x[24] = 0.78173314841662494040636002;
    x[25] = 0.83992032014626734008690454;
    x[26] = 0.88976002994827104337419201;
    x[27] = 0.93075699789664816495694576;
    x[28] = 0.96250392509294966178905240;
    x[29] = 0.98468590966515248400246517;
    x[30] = 0.99708748181947707405562655;

    w[0] = 0.0074708315792487758586968750322;
    w[1] = 0.017318620790310582463157996087;
    w[2] = 0.027009019184979421800608708092;
    w[3] = 0.036432273912385464024392010468;
    w[4] = 0.045493707527201102902315857895;
    w[5] = 0.054103082424916853711666259087;
    w[6] = 0.062174786561028426910343543687;
    w[7] = 0.069628583235410366167756126255;
    w[8] = 0.076390386598776616426357674901;
    w[9] = 0.082392991761589263903823367432;
    w[10] = 0.087576740608477876126198069695;
    w[11] = 0.091890113893641478215362871607;
    w[12] = 0.095290242912319512807204197488;
    w[13] = 0.097743335386328725093474010979;
    w[14] = 0.099225011226672307874875514429;
    w[15] = 0.09972054479342645142753383373;
    w[16] = 0.099225011226672307874875514429;
    w[17] = 0.097743335386328725093474010979;
    w[18] = 0.095290242912319512807204197488;
    w[19] = 0.091890113893641478215362871607;
    w[20] = 0.087576740608477876126198069695;
    w[21] = 0.082392991761589263903823367432;
    w[22] = 0.076390386598776616426357674901;
    w[23] = 0.069628583235410366167756126255;
    w[24] = 0.062174786561028426910343543687;
    w[25] = 0.054103082424916853711666259087;
    w[26] = 0.045493707527201102902315857895;
    w[27] = 0.036432273912385464024392010468;
    w[28] = 0.027009019184979421800608708092;
    w[29] = 0.017318620790310582463157996087;
    w[30] = 0.0074708315792487758586968750322;
  }
  else if ( n == 32 )
  {
    x[0] = -0.99726386184948156354498113;
    x[1] = -0.98561151154526833540017504;
    x[2] = -0.96476225558750643077381193;
    x[3] = -0.93490607593773968917091913;
    x[4] = -0.89632115576605212396530724;
    x[5] = -0.84936761373256997013369300;
    x[6] = -0.79448379596794240696309730;
    x[7] = -0.73218211874028968038742667;
    x[8] = -0.66304426693021520097511517;
    x[9] = -0.58771575724076232904074548;
    x[10] = -0.50689990893222939002374747;
    x[11] = -0.42135127613063534536411944;
    x[12] = -0.33186860228212764977991681;
    x[13] = -0.23928736225213707454460321;
    x[14] = -0.14447196158279649348518637;
    x[15] = -0.04830766568773831623481257;
    x[16] = 0.04830766568773831623481257;
    x[17] = 0.14447196158279649348518637;
    x[18] = 0.23928736225213707454460321;
    x[19] = 0.33186860228212764977991681;
    x[20] = 0.42135127613063534536411944;
    x[21] = 0.50689990893222939002374747;
    x[22] = 0.58771575724076232904074548;
    x[23] = 0.66304426693021520097511517;
    x[24] = 0.73218211874028968038742667;
    x[25] = 0.79448379596794240696309730;
    x[26] = 0.84936761373256997013369300;
    x[27] = 0.89632115576605212396530724;
    x[28] = 0.93490607593773968917091913;
    x[29] = 0.96476225558750643077381193;
    x[30] = 0.98561151154526833540017504;
    x[31] = 0.99726386184948156354498113;

    w[0] = 0.007018610009470096600407063739;
    w[1] = 0.016274394730905670605170562206;
    w[2] = 0.025392065309262059455752589789;
    w[3] = 0.034273862913021433102687732252;
    w[4] = 0.042835898022226680656878646606;
    w[5] = 0.050998059262376176196163244690;
    w[6] = 0.058684093478535547145283637300;
    w[7] = 0.06582222277636184683765006371;
    w[8] = 0.07234579410884850622539935648;
    w[9] = 0.07819389578707030647174091883;
    w[10] = 0.08331192422694675522219907460;
    w[11] = 0.08765209300440381114277146275;
    w[12] = 0.09117387869576388471286857711;
    w[13] = 0.09384439908080456563918023767;
    w[14] = 0.09563872007927485941908200220;
    w[15] = 0.09654008851472780056676483006;
    w[16] = 0.09654008851472780056676483006;
    w[17] = 0.09563872007927485941908200220;
    w[18] = 0.09384439908080456563918023767;
    w[19] = 0.09117387869576388471286857711;
    w[20] = 0.08765209300440381114277146275;
    w[21] = 0.08331192422694675522219907460;
    w[22] = 0.07819389578707030647174091883;
    w[23] = 0.07234579410884850622539935648;
    w[24] = 0.06582222277636184683765006371;
    w[25] = 0.058684093478535547145283637300;
    w[26] = 0.050998059262376176196163244690;
    w[27] = 0.042835898022226680656878646606;
    w[28] = 0.034273862913021433102687732252;
    w[29] = 0.025392065309262059455752589789;
    w[30] = 0.016274394730905670605170562206;
    w[31] = 0.007018610009470096600407063739;
  }
  else if ( n == 33 )
  {
    x[0] = -0.99742469424645521726616802;
    x[1] = -0.98645572623064248811037570;
    x[2] = -0.96682290968999276892837771;
    x[3] = -0.93869437261116835035583512;
    x[4] = -0.90231676774343358304053133;
    x[5] = -0.85800965267650406464306148;
    x[6] = -0.80616235627416658979620087;
    x[7] = -0.74723049644956215785905512;
    x[8] = -0.68173195996974278626821595;
    x[9] = -0.61024234583637902730728751;
    x[10] = -0.53338990478634764354889426;
    x[11] = -0.45185001727245069572599328;
    x[12] = -0.36633925774807334107022062;
    x[13] = -0.27760909715249702940324807;
    x[14] = -0.18643929882799157233579876;
    x[15] = -0.09363106585473338567074292;
    x[16] = 0.00000000000000000000000000;
    x[17] = 0.09363106585473338567074292;
    x[18] = 0.18643929882799157233579876;
    x[19] = 0.27760909715249702940324807;
    x[20] = 0.36633925774807334107022062;
    x[21] = 0.45185001727245069572599328;
    x[22] = 0.53338990478634764354889426;
    x[23] = 0.61024234583637902730728751;
    x[24] = 0.68173195996974278626821595;
    x[25] = 0.74723049644956215785905512;
    x[26] = 0.80616235627416658979620087;
    x[27] = 0.85800965267650406464306148;
    x[28] = 0.90231676774343358304053133;
    x[29] = 0.93869437261116835035583512;
    x[30] = 0.96682290968999276892837771;
    x[31] = 0.98645572623064248811037570;
    x[32] = 0.99742469424645521726616802;

    w[0] = 0.0066062278475873780586492352085;
    w[1] = 0.015321701512934676127945768534;
    w[2] = 0.023915548101749480350533257529;
    w[3] = 0.032300358632328953281561447250;
    w[4] = 0.040401541331669591563409790527;
    w[5] = 0.048147742818711695670146880138;
    w[6] = 0.055470846631663561284944495439;
    w[7] = 0.062306482530317480031627725771;
    w[8] = 0.068594572818656712805955073015;
    w[9] = 0.074279854843954149342472175919;
    w[10] = 0.079312364794886738363908384942;
    w[11] = 0.083647876067038707613928014518;
    w[12] = 0.087248287618844337607281670945;
    w[13] = 0.090081958660638577239743705500;
    w[14] = 0.092123986643316846213240977717;
    w[15] = 0.093356426065596116160999126274;
    w[16] = 0.09376844616020999656730454155;
    w[17] = 0.093356426065596116160999126274;
    w[18] = 0.092123986643316846213240977717;
    w[19] = 0.090081958660638577239743705500;
    w[20] = 0.087248287618844337607281670945;
    w[21] = 0.083647876067038707613928014518;
    w[22] = 0.079312364794886738363908384942;
    w[23] = 0.074279854843954149342472175919;
    w[24] = 0.068594572818656712805955073015;
    w[25] = 0.062306482530317480031627725771;
    w[26] = 0.055470846631663561284944495439;
    w[27] = 0.048147742818711695670146880138;
    w[28] = 0.040401541331669591563409790527;
    w[29] = 0.032300358632328953281561447250;
    w[30] = 0.023915548101749480350533257529;
    w[31] = 0.015321701512934676127945768534;
    w[32] = 0.0066062278475873780586492352085;
  }
  else if ( n == 63 )
  {
    x[0] = -0.99928298402912378037893614;
    x[1] = -0.99622401277797010860219336;
    x[2] = -0.99072854689218946681089467;
    x[3] = -0.98280881059372723486251141;
    x[4] = -0.97248403469757002280196068;
    x[5] = -0.95977944975894192707035417;
    x[6] = -0.94472613404100980296637532;
    x[7] = -0.92736092062184320544703138;
    x[8] = -0.90772630277853155803695313;
    x[9] = -0.88587032850785342629029846;
    x[10] = -0.86184648236412371953961184;
    x[11] = -0.83571355431950284347180777;
    x[12] = -0.80753549577345676005146599;
    x[13] = -0.7773812629903723355633302;
    x[14] = -0.7453246483178474178293217;
    x[15] = -0.7114440995848458078514315;
    x[16] = -0.6758225281149860901311033;
    x[17] = -0.6385471058213653850003070;
    x[18] = -0.5997090518776252357390089;
    x[19] = -0.5594034094862850132676978;
    x[20] = -0.5177288132900332481244776;
    x[21] = -0.4747872479948043999222123;
    x[22] = -0.4306837987951116006620889;
    x[23] = -0.3855263942122478924776150;
    x[24] = -0.3394255419745844024688344;
    x[25] = -0.2924940585862514400361572;
    x[26] = -0.2448467932459533627484046;
    x[27] = -0.1966003467915066845576275;
    x[28] = -0.1478727863578719685698391;
    x[29] = -0.0987833564469452795297037;
    x[30] = -0.0494521871161596272342338;
    x[31] = 0.0000000000000000000000000;
    x[32] = 0.0494521871161596272342338;
    x[33] = 0.0987833564469452795297037;
    x[34] = 0.1478727863578719685698391;
    x[35] = 0.1966003467915066845576275;
    x[36] = 0.2448467932459533627484046;
    x[37] = 0.2924940585862514400361572;
    x[38] = 0.3394255419745844024688344;
    x[39] = 0.3855263942122478924776150;
    x[40] = 0.4306837987951116006620889;
    x[41] = 0.4747872479948043999222123;
    x[42] = 0.5177288132900332481244776;
    x[43] = 0.5594034094862850132676978;
    x[44] = 0.5997090518776252357390089;
    x[45] = 0.6385471058213653850003070;
    x[46] = 0.6758225281149860901311033;
    x[47] = 0.7114440995848458078514315;
    x[48] = 0.7453246483178474178293217;
    x[49] = 0.7773812629903723355633302;
    x[50] = 0.8075354957734567600514660;
    x[51] = 0.8357135543195028434718078;
    x[52] = 0.8618464823641237195396118;
    x[53] = 0.8858703285078534262902985;
    x[54] = 0.9077263027785315580369531;
    x[55] = 0.9273609206218432054470314;
    x[56] = 0.9447261340410098029663753;
    x[57] = 0.9597794497589419270703542;
    x[58] = 0.9724840346975700228019607;
    x[59] = 0.9828088105937272348625114;
    x[60] = 0.9907285468921894668108947;
    x[61] = 0.9962240127779701086021934;
    x[62] = 0.9992829840291237803789361;

    w[0] = 0.0018398745955770841170924455540;
    w[1] = 0.0042785083468637618660784110826;
    w[2] = 0.0067102917659601362519069307298;
    w[3] = 0.0091259686763266563540586454218;
    w[4] = 0.011519376076880041750750606149;
    w[5] = 0.013884612616115610824866086368;
    w[6] = 0.016215878410338338882283672975;
    w[7] = 0.018507464460161270409260545805;
    w[8] = 0.020753761258039090775341953421;
    w[9] = 0.022949271004889933148942319562;
    w[10] = 0.025088620553344986618630138068;
    w[11] = 0.027166574359097933225189839439;
    w[12] = 0.029178047208280526945551502154;
    w[13] = 0.031118116622219817508215988557;
    w[14] = 0.032982034883779341765683179672;
    w[15] = 0.034765240645355877697180504643;
    w[16] = 0.036463370085457289630452409788;
    w[17] = 0.038072267584349556763638324928;
    w[18] = 0.039587995891544093984807928149;
    w[19] = 0.041006845759666398635110037009;
    w[20] = 0.042325345020815822982505485403;
    w[21] = 0.043540267083027590798964315704;
    w[22] = 0.044648638825941395370332669517;
    w[23] = 0.045647747876292608685885992609;
    w[24] = 0.046535149245383696510395418747;
    w[25] = 0.047308671312268919080604988339;
    w[26] = 0.047966421137995131411052756195;
    w[27] = 0.048506789097883847864090099146;
    w[28] = 0.048928452820511989944709361549;
    w[29] = 0.049230380423747560785043116988;
    w[30] = 0.049411833039918178967039646117;
    w[31] = 0.04947236662393102088866936042;
    w[32] = 0.049411833039918178967039646117;
    w[33] = 0.049230380423747560785043116988;
    w[34] = 0.048928452820511989944709361549;
    w[35] = 0.048506789097883847864090099146;
    w[36] = 0.047966421137995131411052756195;
    w[37] = 0.047308671312268919080604988339;
    w[38] = 0.046535149245383696510395418747;
    w[39] = 0.045647747876292608685885992609;
    w[40] = 0.044648638825941395370332669517;
    w[41] = 0.043540267083027590798964315704;
    w[42] = 0.042325345020815822982505485403;
    w[43] = 0.041006845759666398635110037009;
    w[44] = 0.039587995891544093984807928149;
    w[45] = 0.038072267584349556763638324928;
    w[46] = 0.036463370085457289630452409788;
    w[47] = 0.034765240645355877697180504643;
    w[48] = 0.032982034883779341765683179672;
    w[49] = 0.031118116622219817508215988557;
    w[50] = 0.029178047208280526945551502154;
    w[51] = 0.027166574359097933225189839439;
    w[52] = 0.025088620553344986618630138068;
    w[53] = 0.022949271004889933148942319562;
    w[54] = 0.020753761258039090775341953421;
    w[55] = 0.018507464460161270409260545805;
    w[56] = 0.016215878410338338882283672975;
    w[57] = 0.013884612616115610824866086368;
    w[58] = 0.011519376076880041750750606149;
    w[59] = 0.0091259686763266563540586454218;
    w[60] = 0.0067102917659601362519069307298;
    w[61] = 0.0042785083468637618660784110826;
    w[62] = 0.0018398745955770841170924455540;
  }
  else if ( n == 64 )
  {
    x[0] = -0.99930504173577213945690562;
    x[1] = -0.99634011677195527934692450;
    x[2] = -0.99101337147674432073938238;
    x[3] = -0.98333625388462595693129930;
    x[4] = -0.97332682778991096374185351;
    x[5] = -0.96100879965205371891861412;
    x[6] = -0.94641137485840281606248149;
    x[7] = -0.92956917213193957582149015;
    x[8] = -0.91052213707850280575638067;
    x[9] = -0.88931544599511410585340404;
    x[10] = -0.86599939815409281976078339;
    x[11] = -0.8406292962525803627516915;
    x[12] = -0.8132653151227975597419233;
    x[13] = -0.7839723589433414076102205;
    x[14] = -0.7528199072605318966118638;
    x[15] = -0.7198818501716108268489402;
    x[16] = -0.6852363130542332425635584;
    x[17] = -0.6489654712546573398577612;
    x[18] = -0.6111553551723932502488530;
    x[19] = -0.5718956462026340342838781;
    x[20] = -0.5312794640198945456580139;
    x[21] = -0.4894031457070529574785263;
    x[22] = -0.4463660172534640879849477;
    x[23] = -0.4022701579639916036957668;
    x[24] = -0.3572201583376681159504426;
    x[25] = -0.3113228719902109561575127;
    x[26] = -0.2646871622087674163739642;
    x[27] = -0.2174236437400070841496487;
    x[28] = -0.1696444204239928180373136;
    x[29] = -0.1214628192961205544703765;
    x[30] = -0.0729931217877990394495429;
    x[31] = -0.0243502926634244325089558;
    x[32] = 0.0243502926634244325089558;
    x[33] = 0.0729931217877990394495429;
    x[34] = 0.1214628192961205544703765;
    x[35] = 0.1696444204239928180373136;
    x[36] = 0.2174236437400070841496487;
    x[37] = 0.2646871622087674163739642;
    x[38] = 0.3113228719902109561575127;
    x[39] = 0.3572201583376681159504426;
    x[40] = 0.4022701579639916036957668;
    x[41] = 0.4463660172534640879849477;
    x[42] = 0.4894031457070529574785263;
    x[43] = 0.5312794640198945456580139;
    x[44] = 0.5718956462026340342838781;
    x[45] = 0.6111553551723932502488530;
    x[46] = 0.6489654712546573398577612;
    x[47] = 0.6852363130542332425635584;
    x[48] = 0.7198818501716108268489402;
    x[49] = 0.7528199072605318966118638;
    x[50] = 0.7839723589433414076102205;
    x[51] = 0.8132653151227975597419233;
    x[52] = 0.8406292962525803627516915;
    x[53] = 0.8659993981540928197607834;
    x[54] = 0.8893154459951141058534040;
    x[55] = 0.9105221370785028057563807;
    x[56] = 0.9295691721319395758214902;
    x[57] = 0.9464113748584028160624815;
    x[58] = 0.9610087996520537189186141;
    x[59] = 0.9733268277899109637418535;
    x[60] = 0.9833362538846259569312993;
    x[61] = 0.9910133714767443207393824;
    x[62] = 0.9963401167719552793469245;
    x[63] = 0.9993050417357721394569056;

    w[0] = 0.0017832807216964329472960791450;
    w[1] = 0.0041470332605624676352875357286;
    w[2] = 0.006504457968978362856117360400;
    w[3] = 0.008846759826363947723030914660;
    w[4] = 0.011168139460131128818590493019;
    w[5] = 0.013463047896718642598060766686;
    w[6] = 0.015726030476024719321965995298;
    w[7] = 0.017951715775697343085045302001;
    w[8] = 0.020134823153530209372340316729;
    w[9] = 0.022270173808383254159298330384;
    w[10] = 0.024352702568710873338177550409;
    w[11] = 0.026377469715054658671691792625;
    w[12] = 0.028339672614259483227511305200;
    w[13] = 0.030234657072402478867974059820;
    w[14] = 0.032057928354851553585467504348;
    w[15] = 0.033805161837141609391565482111;
    w[16] = 0.035472213256882383810693146715;
    w[17] = 0.037055128540240046040415101810;
    w[18] = 0.038550153178615629128962496947;
    w[19] = 0.039953741132720341386656926128;
    w[20] = 0.041262563242623528610156297474;
    w[21] = 0.042473515123653589007339767909;
    w[22] = 0.043583724529323453376827860974;
    w[23] = 0.044590558163756563060134710031;
    w[24] = 0.045491627927418144479770996971;
    w[25] = 0.046284796581314417295953249232;
    w[26] = 0.046968182816210017325326285755;
    w[27] = 0.047540165714830308662282206944;
    w[28] = 0.04799938859645830772812617987;
    w[29] = 0.04834476223480295716976952716;
    w[30] = 0.04857546744150342693479906678;
    w[31] = 0.04869095700913972038336539073;
    w[32] = 0.04869095700913972038336539073;
    w[33] = 0.04857546744150342693479906678;
    w[34] = 0.04834476223480295716976952716;
    w[35] = 0.04799938859645830772812617987;
    w[36] = 0.047540165714830308662282206944;
    w[37] = 0.046968182816210017325326285755;
    w[38] = 0.046284796581314417295953249232;
    w[39] = 0.045491627927418144479770996971;
    w[40] = 0.044590558163756563060134710031;
    w[41] = 0.043583724529323453376827860974;
    w[42] = 0.042473515123653589007339767909;
    w[43] = 0.041262563242623528610156297474;
    w[44] = 0.039953741132720341386656926128;
    w[45] = 0.038550153178615629128962496947;
    w[46] = 0.037055128540240046040415101810;
    w[47] = 0.035472213256882383810693146715;
    w[48] = 0.033805161837141609391565482111;
    w[49] = 0.032057928354851553585467504348;
    w[50] = 0.030234657072402478867974059820;
    w[51] = 0.028339672614259483227511305200;
    w[52] = 0.026377469715054658671691792625;
    w[53] = 0.024352702568710873338177550409;
    w[54] = 0.022270173808383254159298330384;
    w[55] = 0.020134823153530209372340316729;
    w[56] = 0.017951715775697343085045302001;
    w[57] = 0.015726030476024719321965995298;
    w[58] = 0.013463047896718642598060766686;
    w[59] = 0.011168139460131128818590493019;
    w[60] = 0.008846759826363947723030914660;
    w[61] = 0.006504457968978362856117360400;
    w[62] = 0.0041470332605624676352875357286;
    w[63] = 0.0017832807216964329472960791450;
  }
  else if ( n == 65 )
  {
    x[0] = -0.99932609707541287726569361;
    x[1] = -0.99645094806184916305579494;
    x[2] = -0.99128527617680166872182118;
    x[3] = -0.98383981218703494137763778;
    x[4] = -0.97413153983355116907496789;
    x[5] = -0.96218275471805523771198375;
    x[6] = -0.94802092816840750637376974;
    x[7] = -0.93167862822874933796567699;
    x[8] = -0.91319344054284626173654692;
    x[9] = -0.89260788050473893142328554;
    x[10] = -0.8699692949264070361941320;
    x[11] = -0.8453297528999302839424500;
    x[12] = -0.8187459259226514534339191;
    x[13] = -0.7902789574921218430473804;
    x[14] = -0.7599943224419997868739828;
    x[15] = -0.7279616763294246790119737;
    x[16] = -0.6942546952139916335526225;
    x[17] = -0.6589509061936251330409408;
    x[18] = -0.6221315090854002415825996;
    x[19] = -0.5838811896604873133271545;
    x[20] = -0.5442879248622271385455725;
    x[21] = -0.5034427804550068823410431;
    x[22] = -0.4614397015691450576978341;
    x[23] = -0.4183752966234090092641990;
    x[24] = -0.3743486151220660120087939;
    x[25] = -0.3294609198374864076452867;
    x[26] = -0.2838154539022487306176554;
    x[27] = -0.2375172033464168065707124;
    x[28] = -0.1906726556261427697749124;
    x[29] = -0.1433895546989751711312496;
    x[30] = -0.0957766532091975056522186;
    x[31] = -0.0479434623531718575225298;
    x[32] = 0.0000000000000000000000000;
    x[33] = 0.0479434623531718575225298;
    x[34] = 0.0957766532091975056522186;
    x[35] = 0.1433895546989751711312496;
    x[36] = 0.1906726556261427697749124;
    x[37] = 0.2375172033464168065707124;
    x[38] = 0.2838154539022487306176554;
    x[39] = 0.3294609198374864076452867;
    x[40] = 0.3743486151220660120087939;
    x[41] = 0.4183752966234090092641990;
    x[42] = 0.4614397015691450576978341;
    x[43] = 0.5034427804550068823410431;
    x[44] = 0.5442879248622271385455725;
    x[45] = 0.5838811896604873133271545;
    x[46] = 0.6221315090854002415825996;
    x[47] = 0.6589509061936251330409408;
    x[48] = 0.6942546952139916335526225;
    x[49] = 0.7279616763294246790119737;
    x[50] = 0.7599943224419997868739828;
    x[51] = 0.7902789574921218430473804;
    x[52] = 0.8187459259226514534339191;
    x[53] = 0.8453297528999302839424500;
    x[54] = 0.8699692949264070361941320;
    x[55] = 0.8926078805047389314232855;
    x[56] = 0.9131934405428462617365469;
    x[57] = 0.9316786282287493379656770;
    x[58] = 0.9480209281684075063737697;
    x[59] = 0.9621827547180552377119837;
    x[60] = 0.9741315398335511690749679;
    x[61] = 0.9838398121870349413776378;
    x[62] = 0.9912852761768016687218212;
    x[63] = 0.9964509480618491630557949;
    x[64] = 0.9993260970754128772656936;

    w[0] = 0.0017292582513002508983395851463;
    w[1] = 0.0040215241720037363470786599528;
    w[2] = 0.0063079425789717545501888719039;
    w[3] = 0.0085801482668814598936358121592;
    w[4] = 0.0108326787895979686215140551272;
    w[5] = 0.013060311639994846336168342922;
    w[6] = 0.015257912146448310349265388145;
    w[7] = 0.017420421997670248495365759969;
    w[8] = 0.019542865836750062826837429313;
    w[9] = 0.021620361284934062841654274667;
    w[10] = 0.023648129691287236698780978994;
    w[11] = 0.025621506938037758214084978694;
    w[12] = 0.027535954088450343942499722327;
    w[13] = 0.029387067789310668062644859210;
    w[14] = 0.031170590380189142464431845777;
    w[15] = 0.032882419676368574984049638008;
    w[16] = 0.034518618398549058625221276859;
    w[17] = 0.036075423225565273932166270524;
    w[18] = 0.037549253448257709809772223198;
    w[19] = 0.038936719204051197616673806364;
    w[20] = 0.040234629273005533815446337743;
    w[21] = 0.041439998417240293022686299233;
    w[22] = 0.042550054246755802719217150803;
    w[23] = 0.043562243595800486532284821661;
    w[24] = 0.044474238395082974427323504000;
    w[25] = 0.045283941026300230657128240574;
    w[26] = 0.045989489146651696963893390818;
    w[27] = 0.046589259972233498302255136790;
    w[28] = 0.047081874010454522246006808290;
    w[29] = 0.047466198232885503152644458740;
    w[30] = 0.047741348681240621559038972227;
    w[31] = 0.047906692500495862031347289176;
    w[32] = 0.04796184939446661812070762137;
    w[33] = 0.047906692500495862031347289176;
    w[34] = 0.047741348681240621559038972227;
    w[35] = 0.047466198232885503152644458740;
    w[36] = 0.047081874010454522246006808290;
    w[37] = 0.046589259972233498302255136790;
    w[38] = 0.045989489146651696963893390818;
    w[39] = 0.045283941026300230657128240574;
    w[40] = 0.044474238395082974427323504000;
    w[41] = 0.043562243595800486532284821661;
    w[42] = 0.042550054246755802719217150803;
    w[43] = 0.041439998417240293022686299233;
    w[44] = 0.040234629273005533815446337743;
    w[45] = 0.038936719204051197616673806364;
    w[46] = 0.037549253448257709809772223198;
    w[47] = 0.036075423225565273932166270524;
    w[48] = 0.034518618398549058625221276859;
    w[49] = 0.032882419676368574984049638008;
    w[50] = 0.031170590380189142464431845777;
    w[51] = 0.029387067789310668062644859210;
    w[52] = 0.027535954088450343942499722327;
    w[53] = 0.025621506938037758214084978694;
    w[54] = 0.023648129691287236698780978994;
    w[55] = 0.021620361284934062841654274667;
    w[56] = 0.019542865836750062826837429313;
    w[57] = 0.017420421997670248495365759969;
    w[58] = 0.015257912146448310349265388145;
    w[59] = 0.013060311639994846336168342922;
    w[60] = 0.0108326787895979686215140551272;
    w[61] = 0.0085801482668814598936358121592;
    w[62] = 0.0063079425789717545501888719039;
    w[63] = 0.0040215241720037363470786599528;
    w[64] = 0.0017292582513002508983395851463;
  }
  else if ( n == 127 )
  {
    x[0] = -0.9998221304153061462673512;
    x[1] = -0.9990629343553118951383159;
    x[2] = -0.9976975661898046210744170;
    x[3] = -0.9957265513520272266354334;
    x[4] = -0.9931510492545171473611308;
    x[5] = -0.9899726145914841576077867;
    x[6] = -0.9861931740169316667104383;
    x[7] = -0.9818150208038141100334631;
    x[8] = -0.9768408123430703268174439;
    x[9] = -0.9712735681615291922889469;
    x[10] = -0.9651166679452921210908251;
    x[11] = -0.9583738494252387711491029;
    x[12] = -0.9510492060778803105479076;
    x[13] = -0.9431471846248148273454496;
    x[14] = -0.9346725823247379685736349;
    x[15] = -0.9256305440562338491274647;
    x[16] = -0.9160265591914658093130886;
    x[17] = -0.9058664582618213828024613;
    x[18] = -0.8951564094170837089690438;
    x[19] = -0.8839029146800265699452579;
    x[20] = -0.8721128059985607114196375;
    x[21] = -0.8597932410977408098120313;
    x[22] = -0.8469516991340975984533393;
    x[23] = -0.8335959761548995143795572;
    x[24] = -0.8197341803650786741551191;
    x[25] = -0.8053747272046802146665608;
    x[26] = -0.7905263342398137999454500;
    x[27] = -0.7751980158702023824449628;
    x[28] = -0.7593990778565366715566637;
    x[29] = -0.7431391116709545129205669;
    x[30] = -0.7264279886740726855356929;
    x[31] = -0.7092758541221045609994446;
    x[32] = -0.6916931210077006701564414;
    x[33] = -0.6736904637382504853466825;
    x[34] = -0.6552788116554826302767651;
    x[35] = -0.6364693424002972413476082;
    x[36] = -0.6172734751268582838576392;
    x[37] = -0.5977028635700652293844120;
    x[38] = -0.5777693889706125800032517;
    x[39] = -0.5574851528619322329218619;
    x[40] = -0.5368624697233975674581664;
    x[41] = -0.5159138595042493572772773;
    x[42] = -0.4946520400227821173949402;
    x[43] = -0.4730899192454052416450999;
    x[44] = -0.4512405874502662273318986;
    x[45] = -0.4291173092801933762625441;
    x[46] = -0.4067335156897825634086729;
    x[47] = -0.3841027957915169357790778;
    x[48] = -0.3612388886058697060709248;
    x[49] = -0.3381556747203985013760003;
    x[50] = -0.3148671678628949814860148;
    x[51] = -0.2913875063937056207945188;
    x[52] = -0.2677309447223886208883435;
    x[53] = -0.2439118446539178579707132;
    x[54] = -0.2199446666696875424545234;
    x[55] = -0.1958439611486108515042816;
    x[56] = -0.1716243595336421650083449;
    x[57] = -0.1473005654490856693893293;
    x[58] = -0.1228873457740829717260337;
    x[59] = -0.0983995216776989707510918;
    x[60] = -0.0738519596210485452734404;
    x[61] = -0.0492595623319266303153793;
    x[62] = -0.0246372597574209446148971;
    x[63] = 0.0000000000000000000000000;
    x[64] = 0.0246372597574209446148971;
    x[65] = 0.0492595623319266303153793;
    x[66] = 0.0738519596210485452734404;
    x[67] = 0.0983995216776989707510918;
    x[68] = 0.1228873457740829717260337;
    x[69] = 0.1473005654490856693893293;
    x[70] = 0.1716243595336421650083449;
    x[71] = 0.1958439611486108515042816;
    x[72] = 0.2199446666696875424545234;
    x[73] = 0.2439118446539178579707132;
    x[74] = 0.2677309447223886208883435;
    x[75] = 0.2913875063937056207945188;
    x[76] = 0.3148671678628949814860148;
    x[77] = 0.3381556747203985013760003;
    x[78] = 0.3612388886058697060709248;
    x[79] = 0.3841027957915169357790778;
    x[80] = 0.4067335156897825634086729;
    x[81] = 0.4291173092801933762625441;
    x[82] = 0.4512405874502662273318986;
    x[83] = 0.4730899192454052416450999;
    x[84] = 0.4946520400227821173949402;
    x[85] = 0.5159138595042493572772773;
    x[86] = 0.5368624697233975674581664;
    x[87] = 0.5574851528619322329218619;
    x[88] = 0.5777693889706125800032517;
    x[89] = 0.5977028635700652293844120;
    x[90] = 0.6172734751268582838576392;
    x[91] = 0.6364693424002972413476082;
    x[92] = 0.6552788116554826302767651;
    x[93] = 0.6736904637382504853466825;
    x[94] = 0.6916931210077006701564414;
    x[95] = 0.7092758541221045609994446;
    x[96] = 0.7264279886740726855356929;
    x[97] = 0.7431391116709545129205669;
    x[98] = 0.7593990778565366715566637;
    x[99] = 0.7751980158702023824449628;
    x[100] = 0.7905263342398137999454500;
    x[101] = 0.8053747272046802146665608;
    x[102] = 0.8197341803650786741551191;
    x[103] = 0.8335959761548995143795572;
    x[104] = 0.8469516991340975984533393;
    x[105] = 0.8597932410977408098120313;
    x[106] = 0.8721128059985607114196375;
    x[107] = 0.8839029146800265699452579;
    x[108] = 0.8951564094170837089690438;
    x[109] = 0.9058664582618213828024613;
    x[110] = 0.9160265591914658093130886;
    x[111] = 0.9256305440562338491274647;
    x[112] = 0.9346725823247379685736349;
    x[113] = 0.9431471846248148273454496;
    x[114] = 0.9510492060778803105479076;
    x[115] = 0.9583738494252387711491029;
    x[116] = 0.965116667945292121090825;
    x[117] = 0.971273568161529192288947;
    x[118] = 0.976840812343070326817444;
    x[119] = 0.981815020803814110033463;
    x[120] = 0.986193174016931666710438;
    x[121] = 0.989972614591484157607787;
    x[122] = 0.993151049254517147361131;
    x[123] = 0.995726551352027226635433;
    x[124] = 0.997697566189804621074417;
    x[125] = 0.999062934355311895138316;
    x[126] = 0.999822130415306146267351;

    w[0] = 0.00045645726109586662791936519265;
    w[1] = 0.00106227668695384869596523598532;
    w[2] = 0.0016683488125171936761028862915;
    w[3] = 0.0022734860707492547802810840776;
    w[4] = 0.0028772587656289004082883197514;
    w[5] = 0.0034792893810051465908910894100;
    w[6] = 0.0040792095178254605327114733457;
    w[7] = 0.0046766539777779034772638165663;
    w[8] = 0.0052712596565634400891303815906;
    w[9] = 0.0058626653903523901033648343751;
    w[10] = 0.0064505120486899171845442463869;
    w[11] = 0.0070344427036681608755685893033;
    w[12] = 0.0076141028256526859356393930849;
    w[13] = 0.0081891404887415730817235884719;
    w[14] = 0.0087592065795403145773316804234;
    w[15] = 0.0093239550065309714787536985834;
    w[16] = 0.0098830429087554914716648010900;
    w[17] = 0.0104361308631410052256731719977;
    w[18] = 0.0109828830900689757887996573761;
    w[19] = 0.011522967656921087154811609735;
    w[20] = 0.012056056679400848183529562145;
    w[21] = 0.012581826520465013101514365424;
    w[22] = 0.013099957986718627426172681913;
    w[23] = 0.013610136522139249906034237534;
    w[24] = 0.014112052399003395774044161634;
    w[25] = 0.014605400905893418351737288079;
    w[26] = 0.015089882532666922992635733981;
    w[27] = 0.015565203152273955098532590263;
    w[28] = 0.016031074199309941802254151843;
    w[29] = 0.016487212845194879399346060358;
    w[30] = 0.016933342169871654545878815295;
    w[31] = 0.017369191329918731922164721250;
    w[32] = 0.017794495722974774231027912900;
    w[33] = 0.018208997148375106468721469154;
    w[34] = 0.018612443963902310429440419899;
    w[35] = 0.019004591238555646611148901045;
    w[36] = 0.019385200901246454628112623489;
    w[37] = 0.019754041885329183081815217323;
    w[38] = 0.020110890268880247225644623956;
    w[39] = 0.020455529410639508279497065713;
    w[40] = 0.020787750081531811812652137291;
    w[41] = 0.021107350591688713643523847922;
    w[42] = 0.021414136912893259295449693234;
    w[43] = 0.021707922796373466052301324695;
    w[44] = 0.021988529885872983756478409759;
    w[45] = 0.022255787825930280235631416460;
    w[46] = 0.022509534365300608085694429903;
    w[47] = 0.022749615455457959852242553241;
    w[48] = 0.022975885344117206754377437839;
    w[49] = 0.023188206663719640249922582982;
    w[50] = 0.023386450514828194170722043497;
    w[51] = 0.023570496544381716050033676844;
    w[52] = 0.023740233018760777777714726703;
    w[53] = 0.023895556891620665983864481754;
    w[54] = 0.024036373866450369675132086026;
    w[55] = 0.024162598453819584716522917711;
    w[56] = 0.024274154023278979833195063937;
    w[57] = 0.024370972849882214952813561907;
    w[58] = 0.024452996155301467956140198472;
    w[59] = 0.024520174143511508275183033290;
    w[60] = 0.024572466031020653286354137335;
    w[61] = 0.024609840071630254092545634003;
    w[62] = 0.024632273575707679066033370218;
    w[63] = 0.02463975292396109441957941748;
    w[64] = 0.024632273575707679066033370218;
    w[65] = 0.024609840071630254092545634003;
    w[66] = 0.024572466031020653286354137335;
    w[67] = 0.024520174143511508275183033290;
    w[68] = 0.024452996155301467956140198472;
    w[69] = 0.024370972849882214952813561907;
    w[70] = 0.024274154023278979833195063937;
    w[71] = 0.024162598453819584716522917711;
    w[72] = 0.024036373866450369675132086026;
    w[73] = 0.023895556891620665983864481754;
    w[74] = 0.023740233018760777777714726703;
    w[75] = 0.023570496544381716050033676844;
    w[76] = 0.023386450514828194170722043497;
    w[77] = 0.023188206663719640249922582982;
    w[78] = 0.022975885344117206754377437839;
    w[79] = 0.022749615455457959852242553241;
    w[80] = 0.022509534365300608085694429903;
    w[81] = 0.022255787825930280235631416460;
    w[82] = 0.021988529885872983756478409759;
    w[83] = 0.021707922796373466052301324695;
    w[84] = 0.021414136912893259295449693234;
    w[85] = 0.021107350591688713643523847922;
    w[86] = 0.020787750081531811812652137291;
    w[87] = 0.020455529410639508279497065713;
    w[88] = 0.020110890268880247225644623956;
    w[89] = 0.019754041885329183081815217323;
    w[90] = 0.019385200901246454628112623489;
    w[91] = 0.019004591238555646611148901045;
    w[92] = 0.018612443963902310429440419899;
    w[93] = 0.018208997148375106468721469154;
    w[94] = 0.017794495722974774231027912900;
    w[95] = 0.017369191329918731922164721250;
    w[96] = 0.016933342169871654545878815295;
    w[97] = 0.016487212845194879399346060358;
    w[98] = 0.016031074199309941802254151843;
    w[99] = 0.015565203152273955098532590263;
    w[100] = 0.015089882532666922992635733981;
    w[101] = 0.014605400905893418351737288079;
    w[102] = 0.014112052399003395774044161634;
    w[103] = 0.013610136522139249906034237534;
    w[104] = 0.013099957986718627426172681913;
    w[105] = 0.012581826520465013101514365424;
    w[106] = 0.012056056679400848183529562145;
    w[107] = 0.011522967656921087154811609735;
    w[108] = 0.0109828830900689757887996573761;
    w[109] = 0.0104361308631410052256731719977;
    w[110] = 0.0098830429087554914716648010900;
    w[111] = 0.0093239550065309714787536985834;
    w[112] = 0.0087592065795403145773316804234;
    w[113] = 0.0081891404887415730817235884719;
    w[114] = 0.0076141028256526859356393930849;
    w[115] = 0.0070344427036681608755685893033;
    w[116] = 0.0064505120486899171845442463869;
    w[117] = 0.0058626653903523901033648343751;
    w[118] = 0.0052712596565634400891303815906;
    w[119] = 0.0046766539777779034772638165663;
    w[120] = 0.0040792095178254605327114733457;
    w[121] = 0.0034792893810051465908910894100;
    w[122] = 0.0028772587656289004082883197514;
    w[123] = 0.0022734860707492547802810840776;
    w[124] = 0.0016683488125171936761028862915;
    w[125] = 0.00106227668695384869596523598532;
    w[126] = 0.00045645726109586662791936519265;
  }
  else if ( n == 128 )
  {
    x[0] = -0.9998248879471319144736081;
    x[1] = -0.9990774599773758950119878;
    x[2] = -0.9977332486255140198821574;
    x[3] = -0.9957927585349811868641612;
    x[4] = -0.9932571129002129353034372;
    x[5] = -0.9901278184917343833379303;
    x[6] = -0.9864067427245862088712355;
    x[7] = -0.9820961084357185360247656;
    x[8] = -0.9771984914639073871653744;
    x[9] = -0.9717168187471365809043384;
    x[10] = -0.9656543664319652686458290;
    x[11] = -0.9590147578536999280989185;
    x[12] = -0.9518019613412643862177963;
    x[13] = -0.9440202878302201821211114;
    x[14] = -0.9356743882779163757831268;
    x[15] = -0.9267692508789478433346245;
    x[16] = -0.9173101980809605370364836;
    x[17] = -0.9073028834017568139214859;
    x[18] = -0.8967532880491581843864474;
    x[19] = -0.8856677173453972174082924;
    x[20] = -0.8740527969580317986954180;
    x[21] = -0.8619154689395484605906323;
    x[22] = -0.8492629875779689691636001;
    x[23] = -0.8361029150609068471168753;
    x[24] = -0.8224431169556438424645942;
    x[25] = -0.8082917575079136601196422;
    x[26] = -0.7936572947621932902433329;
    x[27] = -0.7785484755064119668504941;
    x[28] = -0.7629743300440947227797691;
    x[29] = -0.7469441667970619811698824;
    x[30] = -0.7304675667419088064717369;
    x[31] = -0.7135543776835874133438599;
    x[32] = -0.6962147083695143323850866;
    x[33] = -0.6784589224477192593677557;
    x[34] = -0.6602976322726460521059468;
    x[35] = -0.6417416925623075571535249;
    x[36] = -0.6228021939105849107615396;
    x[37] = -0.6034904561585486242035732;
    x[38] = -0.5838180216287630895500389;
    x[39] = -0.5637966482266180839144308;
    x[40] = -0.5434383024128103634441936;
    x[41] = -0.5227551520511754784539479;
    x[42] = -0.5017595591361444642896063;
    x[43] = -0.4804640724041720258582757;
    x[44] = -0.4588814198335521954490891;
    x[45] = -0.4370245010371041629370429;
    x[46] = -0.4149063795522750154922739;
    x[47] = -0.3925402750332674427356482;
    x[48] = -0.3699395553498590266165917;
    x[49] = -0.3471177285976355084261628;
    x[50] = -0.3240884350244133751832523;
    x[51] = -0.3008654388776772026671541;
    x[52] = -0.2774626201779044028062316;
    x[53] = -0.2538939664226943208556180;
    x[54] = -0.2301735642266599864109866;
    x[55] = -0.2063155909020792171540580;
    x[56] = -0.1823343059853371824103826;
    x[57] = -0.1582440427142249339974755;
    x[58] = -0.1340591994611877851175753;
    x[59] = -0.1097942311276437466729747;
    x[60] = -0.0854636405045154986364980;
    x[61] = -0.0610819696041395681037870;
    x[62] = -0.0366637909687334933302153;
    x[63] = -0.0122236989606157641980521;
    x[64] = 0.0122236989606157641980521;
    x[65] = 0.0366637909687334933302153;
    x[66] = 0.0610819696041395681037870;
    x[67] = 0.0854636405045154986364980;
    x[68] = 0.1097942311276437466729747;
    x[69] = 0.1340591994611877851175753;
    x[70] = 0.1582440427142249339974755;
    x[71] = 0.1823343059853371824103826;
    x[72] = 0.2063155909020792171540580;
    x[73] = 0.2301735642266599864109866;
    x[74] = 0.2538939664226943208556180;
    x[75] = 0.2774626201779044028062316;
    x[76] = 0.3008654388776772026671541;
    x[77] = 0.3240884350244133751832523;
    x[78] = 0.3471177285976355084261628;
    x[79] = 0.3699395553498590266165917;
    x[80] = 0.3925402750332674427356482;
    x[81] = 0.4149063795522750154922739;
    x[82] = 0.4370245010371041629370429;
    x[83] = 0.4588814198335521954490891;
    x[84] = 0.4804640724041720258582757;
    x[85] = 0.5017595591361444642896063;
    x[86] = 0.5227551520511754784539479;
    x[87] = 0.5434383024128103634441936;
    x[88] = 0.5637966482266180839144308;
    x[89] = 0.5838180216287630895500389;
    x[90] = 0.6034904561585486242035732;
    x[91] = 0.6228021939105849107615396;
    x[92] = 0.6417416925623075571535249;
    x[93] = 0.6602976322726460521059468;
    x[94] = 0.6784589224477192593677557;
    x[95] = 0.6962147083695143323850866;
    x[96] = 0.7135543776835874133438599;
    x[97] = 0.7304675667419088064717369;
    x[98] = 0.7469441667970619811698824;
    x[99] = 0.7629743300440947227797691;
    x[100] = 0.7785484755064119668504941;
    x[101] = 0.7936572947621932902433329;
    x[102] = 0.8082917575079136601196422;
    x[103] = 0.8224431169556438424645942;
    x[104] = 0.8361029150609068471168753;
    x[105] = 0.8492629875779689691636001;
    x[106] = 0.8619154689395484605906323;
    x[107] = 0.8740527969580317986954180;
    x[108] = 0.8856677173453972174082924;
    x[109] = 0.8967532880491581843864474;
    x[110] = 0.9073028834017568139214859;
    x[111] = 0.9173101980809605370364836;
    x[112] = 0.926769250878947843334625;
    x[113] = 0.935674388277916375783127;
    x[114] = 0.944020287830220182121111;
    x[115] = 0.951801961341264386217796;
    x[116] = 0.959014757853699928098919;
    x[117] = 0.965654366431965268645829;
    x[118] = 0.971716818747136580904338;
    x[119] = 0.977198491463907387165374;
    x[120] = 0.982096108435718536024766;
    x[121] = 0.986406742724586208871236;
    x[122] = 0.990127818491734383337930;
    x[123] = 0.993257112900212935303437;
    x[124] = 0.995792758534981186864161;
    x[125] = 0.997733248625514019882157;
    x[126] = 0.999077459977375895011988;
    x[127] = 0.999824887947131914473608;

    w[0] = 0.00044938096029209037639429223999;
    w[1] = 0.0010458126793403487793128516001;
    w[2] = 0.0016425030186690295387908755948;
    w[3] = 0.0022382884309626187436220542727;
    w[4] = 0.0028327514714579910952857346468;
    w[5] = 0.0034255260409102157743377846601;
    w[6] = 0.0040162549837386423131943434863;
    w[7] = 0.0046045842567029551182905419803;
    w[8] = 0.0051901618326763302050707671348;
    w[9] = 0.0057726375428656985893346176261;
    w[10] = 0.006351663161707188787214327826;
    w[11] = 0.006926892566898813563426670360;
    w[12] = 0.007497981925634728687671962688;
    w[13] = 0.008064589890486057972928598698;
    w[14] = 0.008626377798616749704978843782;
    w[15] = 0.009183009871660874334478743688;
    w[16] = 0.009734153415006805863548266094;
    w[17] = 0.010279479015832157133215340326;
    w[18] = 0.010818660739503076247659646277;
    w[19] = 0.011351376324080416693281668453;
    w[20] = 0.011877307372740279575891106926;
    w[21] = 0.012396139543950922968821728197;
    w[22] = 0.012907562739267347220442834004;
    w[23] = 0.013411271288616332314488951616;
    w[24] = 0.013906964132951985244288007396;
    w[25] = 0.014394345004166846176823892009;
    w[26] = 0.014873122602147314252385498520;
    w[27] = 0.015343010768865144085990853741;
    w[28] = 0.015803728659399346858965631687;
    w[29] = 0.016255000909785187051657456477;
    w[30] = 0.016696557801589204589091507954;
    w[31] = 0.017128135423111376830680987619;
    w[32] = 0.017549475827117704648706925634;
    w[33] = 0.017960327185008685940196927525;
    w[34] = 0.018360443937331343221289290991;
    w[35] = 0.018749586940544708650919548474;
    w[36] = 0.019127523609950945486518531668;
    w[37] = 0.019494028058706602823021918681;
    w[38] = 0.019848881232830862219944413265;
    w[39] = 0.020191871042130041180673158406;
    w[40] = 0.020522792486960069432284967788;
    w[41] = 0.020841447780751149113583948423;
    w[42] = 0.021147646468221348537019535180;
    w[43] = 0.021441205539208460137111853878;
    w[44] = 0.021721949538052075375260957768;
    w[45] = 0.021989710668460491434122106599;
    w[46] = 0.022244328893799765104629133607;
    w[47] = 0.022485652032744966871824603941;
    w[48] = 0.022713535850236461309712635923;
    w[49] = 0.022927844143686846920410987209;
    w[50] = 0.023128448824387027879297902403;
    w[51] = 0.023315229994062760122415671273;
    w[52] = 0.023488076016535913153025273282;
    w[53] = 0.023646883584447615143651392303;
    w[54] = 0.023791557781003400638780709885;
    w[55] = 0.023922012136703455672450408817;
    w[56] = 0.024038168681024052637587316820;
    w[57] = 0.024139957989019284997716653890;
    w[58] = 0.024227319222815248120093308442;
    w[59] = 0.024300200167971865323442606364;
    w[60] = 0.024358557264690625853268520246;
    w[61] = 0.024402355633849582093297989694;
    w[62] = 0.02443156909785004505484856143;
    w[63] = 0.02444618019626251821132585261;
    w[64] = 0.02444618019626251821132585261;
    w[65] = 0.02443156909785004505484856143;
    w[66] = 0.024402355633849582093297989694;
    w[67] = 0.024358557264690625853268520246;
    w[68] = 0.024300200167971865323442606364;
    w[69] = 0.024227319222815248120093308442;
    w[70] = 0.024139957989019284997716653890;
    w[71] = 0.024038168681024052637587316820;
    w[72] = 0.023922012136703455672450408817;
    w[73] = 0.023791557781003400638780709885;
    w[74] = 0.023646883584447615143651392303;
    w[75] = 0.023488076016535913153025273282;
    w[76] = 0.023315229994062760122415671273;
    w[77] = 0.023128448824387027879297902403;
    w[78] = 0.022927844143686846920410987209;
    w[79] = 0.022713535850236461309712635923;
    w[80] = 0.022485652032744966871824603941;
    w[81] = 0.022244328893799765104629133607;
    w[82] = 0.021989710668460491434122106599;
    w[83] = 0.021721949538052075375260957768;
    w[84] = 0.021441205539208460137111853878;
    w[85] = 0.021147646468221348537019535180;
    w[86] = 0.020841447780751149113583948423;
    w[87] = 0.020522792486960069432284967788;
    w[88] = 0.020191871042130041180673158406;
    w[89] = 0.019848881232830862219944413265;
    w[90] = 0.019494028058706602823021918681;
    w[91] = 0.019127523609950945486518531668;
    w[92] = 0.018749586940544708650919548474;
    w[93] = 0.018360443937331343221289290991;
    w[94] = 0.017960327185008685940196927525;
    w[95] = 0.017549475827117704648706925634;
    w[96] = 0.017128135423111376830680987619;
    w[97] = 0.016696557801589204589091507954;
    w[98] = 0.016255000909785187051657456477;
    w[99] = 0.015803728659399346858965631687;
    w[100] = 0.015343010768865144085990853741;
    w[101] = 0.014873122602147314252385498520;
    w[102] = 0.014394345004166846176823892009;
    w[103] = 0.013906964132951985244288007396;
    w[104] = 0.013411271288616332314488951616;
    w[105] = 0.012907562739267347220442834004;
    w[106] = 0.012396139543950922968821728197;
    w[107] = 0.011877307372740279575891106926;
    w[108] = 0.011351376324080416693281668453;
    w[109] = 0.010818660739503076247659646277;
    w[110] = 0.010279479015832157133215340326;
    w[111] = 0.009734153415006805863548266094;
    w[112] = 0.009183009871660874334478743688;
    w[113] = 0.008626377798616749704978843782;
    w[114] = 0.008064589890486057972928598698;
    w[115] = 0.007497981925634728687671962688;
    w[116] = 0.006926892566898813563426670360;
    w[117] = 0.006351663161707188787214327826;
    w[118] = 0.0057726375428656985893346176261;
    w[119] = 0.0051901618326763302050707671348;
    w[120] = 0.0046045842567029551182905419803;
    w[121] = 0.0040162549837386423131943434863;
    w[122] = 0.0034255260409102157743377846601;
    w[123] = 0.0028327514714579910952857346468;
    w[124] = 0.0022382884309626187436220542727;
    w[125] = 0.0016425030186690295387908755948;
    w[126] = 0.0010458126793403487793128516001;
    w[127] = 0.00044938096029209037639429223999;
  }
  else if ( n == 129 )
  {
    x[0] = -0.9998275818477487191077441;
    x[1] = -0.9990916504696409986514389;
    x[2] = -0.9977681080525852721429460;
    x[3] = -0.9958574393142831982149111;
    x[4] = -0.9933607326210712814854011;
    x[5] = -0.9902794486488178389207689;
    x[6] = -0.9866153978313475022005761;
    x[7] = -0.9823707352517413115507418;
    x[8] = -0.9775479582993672474447814;
    x[9] = -0.9721499048427034297274163;
    x[10] = -0.9661797514202097197778763;
    x[11] = -0.9596410113101918904168119;
    x[12] = -0.9525375324342090471027732;
    x[13] = -0.9448734950776734726784764;
    x[14] = -0.9366534094216514605284616;
    x[15] = -0.9278821128840036204317296;
    x[16] = -0.9185647672698286252225115;
    x[17] = -0.9087068557320696331245539;
    x[18] = -0.8983141795436338850435985;
    x[19] = -0.8873928546826803665034968;
    x[20] = -0.8759493082329433892035217;
    x[21] = -0.8639902746011257878940216;
    x[22] = -0.8515227915535356930243826;
    x[23] = -0.8385541960742664442975407;
    x[24] = -0.8250921200473358809210133;
    x[25] = -0.8111444857653120742087717;
    x[26] = -0.7967195012670592680339606;
    x[27] = -0.7818256555073413245387500;
    x[28] = -0.7664717133611208816717785;
    x[29] = -0.7506667104654910227632368;
    x[30] = -0.7344199479022727047791516;
    x[31] = -0.7177409867244055767721220;
    x[32] = -0.7006396423293521790044710;
    x[33] = -0.6831259786828258512462248;
    x[34] = -0.6652103023962409818802202;
    x[35] = -0.6469031566613704719753373;
    x[36] = -0.6282153150457794374886895;
    x[37] = -0.6091577751526861909563306;
    x[38] = -0.5897417521489813916767844;
    x[39] = -0.5699786721652138894754096;
    x[40] = -0.5498801655714271702189358;
    x[41] = -0.5294580601328034000099406;
    x[42] = -0.5087243740491428186199463;
    x[43] = -0.4876913088822746111853066;
    x[44] = -0.4663712423755613514331869;
    x[45] = -0.4447767211697226217818454;
    x[46] = -0.4229204534192644388475065;
    x[47] = -0.4008153013138596117693121;
    x[48] = -0.3784742735090801012801265;
    x[49] = -0.3559105174709357969672656;
    x[50] = -0.3331373117387248575049982;
    x[51] = -0.3101680581107488341147318;
    x[52] = -0.2870162737574911929568755;
    x[53] = -0.2636955832669005409666949;
    x[54] = -0.2402197106264598167721148;
    x[55] = -0.2166024711467599103221439;
    x[56] = -0.1928577633313305998663880;
    x[57] = -0.1689995606975133227390302;
    x[58] = -0.1450419035531891084328306;
    x[59] = -0.1209988907342009817690539;
    x[60] = -0.0968846713073332753086909;
    x[61] = -0.0727134362437305599118207;
    x[62] = -0.0484994100676562986191764;
    x[63] = -0.0242568424855058415749954;
    x[64] = 0.0000000000000000000000000;
    x[65] = 0.0242568424855058415749954;
    x[66] = 0.0484994100676562986191764;
    x[67] = 0.0727134362437305599118207;
    x[68] = 0.0968846713073332753086909;
    x[69] = 0.1209988907342009817690539;
    x[70] = 0.1450419035531891084328306;
    x[71] = 0.1689995606975133227390302;
    x[72] = 0.1928577633313305998663880;
    x[73] = 0.2166024711467599103221439;
    x[74] = 0.2402197106264598167721148;
    x[75] = 0.2636955832669005409666949;
    x[76] = 0.2870162737574911929568755;
    x[77] = 0.3101680581107488341147318;
    x[78] = 0.3331373117387248575049982;
    x[79] = 0.3559105174709357969672656;
    x[80] = 0.3784742735090801012801265;
    x[81] = 0.4008153013138596117693121;
    x[82] = 0.4229204534192644388475065;
    x[83] = 0.4447767211697226217818454;
    x[84] = 0.4663712423755613514331869;
    x[85] = 0.4876913088822746111853066;
    x[86] = 0.5087243740491428186199463;
    x[87] = 0.5294580601328034000099406;
    x[88] = 0.5498801655714271702189358;
    x[89] = 0.5699786721652138894754096;
    x[90] = 0.5897417521489813916767844;
    x[91] = 0.6091577751526861909563306;
    x[92] = 0.6282153150457794374886895;
    x[93] = 0.6469031566613704719753373;
    x[94] = 0.6652103023962409818802202;
    x[95] = 0.6831259786828258512462248;
    x[96] = 0.7006396423293521790044710;
    x[97] = 0.7177409867244055767721220;
    x[98] = 0.7344199479022727047791516;
    x[99] = 0.7506667104654910227632368;
    x[100] = 0.7664717133611208816717785;
    x[101] = 0.7818256555073413245387500;
    x[102] = 0.7967195012670592680339606;
    x[103] = 0.8111444857653120742087717;
    x[104] = 0.8250921200473358809210133;
    x[105] = 0.8385541960742664442975407;
    x[106] = 0.8515227915535356930243826;
    x[107] = 0.8639902746011257878940216;
    x[108] = 0.875949308232943389203522;
    x[109] = 0.887392854682680366503497;
    x[110] = 0.898314179543633885043599;
    x[111] = 0.908706855732069633124554;
    x[112] = 0.918564767269828625222511;
    x[113] = 0.927882112884003620431730;
    x[114] = 0.936653409421651460528462;
    x[115] = 0.944873495077673472678476;
    x[116] = 0.952537532434209047102773;
    x[117] = 0.959641011310191890416812;
    x[118] = 0.966179751420209719777876;
    x[119] = 0.972149904842703429727416;
    x[120] = 0.977547958299367247444781;
    x[121] = 0.982370735251741311550742;
    x[122] = 0.986615397831347502200576;
    x[123] = 0.990279448648817838920769;
    x[124] = 0.993360732621071281485401;
    x[125] = 0.995857439314283198214911;
    x[126] = 0.997768108052585272142946;
    x[127] = 0.999091650469640998651439;
    x[128] = 0.999827581847748719107744;

    w[0] = 0.00044246794182939296923668005717;
    w[1] = 0.00102972844619622394463273519315;
    w[2] = 0.0016172530556785534682413679271;
    w[3] = 0.0022039015180966937075786419741;
    w[4] = 0.0027892681877797554940944677057;
    w[5] = 0.0033729979506246246117755709288;
    w[6] = 0.0039547444682113562172392974765;
    w[7] = 0.0045341644298525434513226874954;
    w[8] = 0.0051109164669246267289761565766;
    w[9] = 0.0056846609912469045788016012203;
    w[10] = 0.0062550602724461408889348709586;
    w[11] = 0.0068217785893519121070498527769;
    w[12] = 0.0073844824072454014447165055698;
    w[13] = 0.0079428405646668029041114107832;
    w[14] = 0.0084965244635723279730542832506;
    w[15] = 0.0090452082602137316404219313819;
    w[16] = 0.0095885690555104190787301294510;
    w[17] = 0.0101262870842733548093160774580;
    w[18] = 0.0106580459029055185304204093001;
    w[19] = 0.0111835325753305049735380697538;
    w[20] = 0.011702437856964778185746436834;
    w[21] = 0.012214456376582979416221105914;
    w[22] = 0.012719286815944623465099036330;
    w[23] = 0.013216632087061724231482387345;
    w[24] = 0.013706199506993971244060563234;
    w[25] = 0.014187700970062900419317230938;
    w[26] = 0.014660853117380060971041027493;
    w[27] = 0.015125377503587024690403432771;
    w[28] = 0.015581000760707523415881287558;
    w[29] = 0.016027454759014214436403950465;
    w[30] = 0.016464476764814667467169189640;
    w[31] = 0.016891809595063204177526208819;
    w[32] = 0.017309201768707240731293596444;
    w[33] = 0.017716407654678809269702031810;
    w[34] = 0.018113187616443980503999783812;
    w[35] = 0.018499308153024985727791918518;
    w[36] = 0.018874542036411948181617592169;
    w[37] = 0.019238668445283284085199492202;
    w[38] = 0.019591473094956024580283987216;
    w[39] = 0.019932748363489542089706675388;
    w[40] = 0.020262293413868438317104423081;
    w[41] = 0.020579914312192665948185517085;
    w[42] = 0.020885424141805311409990024684;
    w[43] = 0.021178643113290860912881038703;
    w[44] = 0.021459398670279205389981598196;
    w[45] = 0.021727525590993110687305178710;
    w[46] = 0.021982866085479386179554968899;
    w[47] = 0.022225269888466526554736910919;
    w[48] = 0.022454594347794176432066564511;
    w[49] = 0.022670704508362374313093970958;
    w[50] = 0.022873473191551169638592083492;
    w[51] = 0.023062781070063872924670495006;
    w[52] = 0.023238516738149892544490435771;
    w[53] = 0.023400576777165831146714346635;
    w[54] = 0.023548865816436258377269094263;
    w[55] = 0.023683296589378342897341543485;
    w[56] = 0.023803789984857314051325299744;
    w[57] = 0.023910275093742530302367230296;
    w[58] = 0.024002689250636756075547029720;
    w[59] = 0.024080978070754089272959634041;
    w[60] = 0.024145095481924836783843156014;
    w[61] = 0.024195003751708503129818111597;
    w[62] = 0.024230673509598936275508460625;
    w[63] = 0.024252083764308562906498864071;
    w[64] = 0.02425922191612154143202867472;
    w[65] = 0.024252083764308562906498864071;
    w[66] = 0.024230673509598936275508460625;
    w[67] = 0.024195003751708503129818111597;
    w[68] = 0.024145095481924836783843156014;
    w[69] = 0.024080978070754089272959634041;
    w[70] = 0.024002689250636756075547029720;
    w[71] = 0.023910275093742530302367230296;
    w[72] = 0.023803789984857314051325299744;
    w[73] = 0.023683296589378342897341543485;
    w[74] = 0.023548865816436258377269094263;
    w[75] = 0.023400576777165831146714346635;
    w[76] = 0.023238516738149892544490435771;
    w[77] = 0.023062781070063872924670495006;
    w[78] = 0.022873473191551169638592083492;
    w[79] = 0.022670704508362374313093970958;
    w[80] = 0.022454594347794176432066564511;
    w[81] = 0.022225269888466526554736910919;
    w[82] = 0.021982866085479386179554968899;
    w[83] = 0.021727525590993110687305178710;
    w[84] = 0.021459398670279205389981598196;
    w[85] = 0.021178643113290860912881038703;
    w[86] = 0.020885424141805311409990024684;
    w[87] = 0.020579914312192665948185517085;
    w[88] = 0.020262293413868438317104423081;
    w[89] = 0.019932748363489542089706675388;
    w[90] = 0.019591473094956024580283987216;
    w[91] = 0.019238668445283284085199492202;
    w[92] = 0.018874542036411948181617592169;
    w[93] = 0.018499308153024985727791918518;
    w[94] = 0.018113187616443980503999783812;
    w[95] = 0.017716407654678809269702031810;
    w[96] = 0.017309201768707240731293596444;
    w[97] = 0.016891809595063204177526208819;
    w[98] = 0.016464476764814667467169189640;
    w[99] = 0.016027454759014214436403950465;
    w[100] = 0.015581000760707523415881287558;
    w[101] = 0.015125377503587024690403432771;
    w[102] = 0.014660853117380060971041027493;
    w[103] = 0.014187700970062900419317230938;
    w[104] = 0.013706199506993971244060563234;
    w[105] = 0.013216632087061724231482387345;
    w[106] = 0.012719286815944623465099036330;
    w[107] = 0.012214456376582979416221105914;
    w[108] = 0.011702437856964778185746436834;
    w[109] = 0.0111835325753305049735380697538;
    w[110] = 0.0106580459029055185304204093001;
    w[111] = 0.0101262870842733548093160774580;
    w[112] = 0.0095885690555104190787301294510;
    w[113] = 0.0090452082602137316404219313819;
    w[114] = 0.0084965244635723279730542832506;
    w[115] = 0.0079428405646668029041114107832;
    w[116] = 0.0073844824072454014447165055698;
    w[117] = 0.0068217785893519121070498527769;
    w[118] = 0.0062550602724461408889348709586;
    w[119] = 0.0056846609912469045788016012203;
    w[120] = 0.0051109164669246267289761565766;
    w[121] = 0.0045341644298525434513226874954;
    w[122] = 0.0039547444682113562172392974765;
    w[123] = 0.0033729979506246246117755709288;
    w[124] = 0.0027892681877797554940944677057;
    w[125] = 0.0022039015180966937075786419741;
    w[126] = 0.0016172530556785534682413679271;
    w[127] = 0.00102972844619622394463273519315;
    w[128] = 0.00044246794182939296923668005717;
  }
  else if ( n == 255 )
  {
    x[0] = -0.999955705317563751730191;
    x[1] = -0.999766621312000569367063;
    x[2] = -0.999426474680169959344386;
    x[3] = -0.998935241284654635142155;
    x[4] = -0.998292986136967889228248;
    x[5] = -0.997499804126615814044844;
    x[6] = -0.996555814435198617028738;
    x[7] = -0.995461159480026294089975;
    x[8] = -0.994216004616630164799381;
    x[9] = -0.992820538021989138984811;
    x[10] = -0.991274970630385567164523;
    x[11] = -0.989579536085920123498574;
    x[12] = -0.987734490699732356281248;
    x[13] = -0.985740113407419277752900;
    x[14] = -0.983596705724776358640192;
    x[15] = -0.981304591701017185126565;
    x[16] = -0.978864117869068155239121;
    x[17] = -0.976275653192735980815246;
    x[18] = -0.973539589010643617645393;
    x[19] = -0.970656338976880365477697;
    x[20] = -0.967626338998338798105523;
    x[21] = -0.964450047168726298761719;
    x[22] = -0.961127943699247839572910;
    x[23] = -0.957660530845962076295490;
    x[24] = -0.954048332833816317950921;
    x[25] = -0.950291895777368285733522;
    x[26] = -0.946391787598204251752103;
    x[27] = -0.942348597939064408301480;
    x[28] = -0.938162938074687317626793;
    x[29] = -0.933835440819386124349338;
    x[30] = -0.929366760431369935739045;
    x[31] = -0.924757572513824425220425;
    x[32] = -0.920008573912766315142721;
    x[33] = -0.915120482611686961035103;
    x[34] = -0.910094037623000801254172;
    x[35] = -0.904929998876314959753358;
    x[36] = -0.899629147103536800144342;
    x[37] = -0.894192283720836729335637;
    x[38] = -0.888620230707484040924981;
    x[39] = -0.882913830481574073645470;
    x[40] = -0.877073945772665439532627;
    x[41] = -0.871101459491346550796200;
    x[42] = -0.864997274595751144137121;
    x[43] = -0.858762313955042966785823;
    x[44] = -0.852397520209890250084237;
    x[45] = -0.845903855629951054143931;
    x[46] = -0.839282301968391021084600;
    x[47] = -0.832533860313455524647230;
    x[48] = -0.825659550937118650611534;
    x[49] = -0.818660413140831885432406;
    x[50] = -0.811537505098395829833580;
    x[51] = -0.804291903695978689734633;
    x[52] = -0.796924704369305728807154;
    x[53] = -0.789437020938044295117764;
    x[54] = -0.781829985437409458675147;
    x[55] = -0.774104747947015717207115;
    x[56] = -0.766262476417000644100858;
    x[57] = -0.758304356491446765092016;
    x[58] = -0.750231591329128358931528;
    x[59] = -0.742045401421610281838045;
    x[60] = -0.733747024408726316001889;
    x[61] = -0.725337714891464938687812;
    x[62] = -0.716818744242290800531501;
    x[63] = -0.708191400412930589382399;
    x[64] = -0.699456987739652339456557;
    x[65] = -0.690616826746067624571761;
    x[66] = -0.681672253943486448787259;
    x[67] = -0.672624621628855017806731;
    x[68] = -0.663475297680306939970658;
    x[69] = -0.654225665350358766508700;
    x[70] = -0.644877123056781136890077;
    x[71] = -0.635431084171177146547142;
    x[72] = -0.625888976805299900901619;
    x[73] = -0.616252243595141561442344;
    x[74] = -0.606522341482826526536576;
    x[75] = -0.596700741496341721653202;
    x[76] = -0.586788928527137300685706;
    x[77] = -0.576788401105631382036211;
    x[78] = -0.566700671174652760010815;
    x[79] = -0.556527263860855843833077;
    x[80] = -0.546269717244142383159817;
    x[81] = -0.535929582125124840335150;
    x[82] = -0.525508421790666565699453;
    x[83] = -0.515007811777534223035005;
    x[84] = -0.504429339634198197635551;
    x[85] = -0.493774604680816999489812;
    x[86] = -0.483045217767441948626854;
    x[87] = -0.472242801030478698742627;
    x[88] = -0.461368987647442418771401;
    x[89] = -0.450425421590043710043279;
    x[90] = -0.439413757375642589040685;
    x[91] = -0.428335659817108112494341;
    x[92] = -0.417192803771121462605751;
    x[93] = -0.405986873884960545511889;
    x[94] = -0.394719564341804385683361;
    x[95] = -0.383392578604595822734854;
    x[96] = -0.372007629158501235092510;
    x[97] = -0.360566437252006227074021;
    x[98] = -0.349070732636686422161576;
    x[99] = -0.337522253305692705554261;
    x[100] = -0.325922745230990453444769;
    x[101] = -0.314273962099392474845918;
    x[102] = -0.302577665047425574167140;
    x[103] = -0.290835622395070819082047;
    x[104] = -0.279049609378417768508970;
    x[105] = -0.267221407881273079721012;
    x[106] = -0.255352806165764071686080;
    x[107] = -0.243445598601977973686482;
    x[108] = -0.231501585396677734059116;
    x[109] = -0.219522572321135403508985;
    x[110] = -0.207510370438124240859625;
    x[111] = -0.195466795828110816293869;
    x[112] = -0.183393669314688508087976;
    x[113] = -0.171292816189293903533225;
    x[114] = -0.159166065935247723154292;
    x[115] = -0.147015251951161989456661;
    x[116] = -0.134842211273755257250625;
    x[117] = -0.122648784300117812092492;
    x[118] = -0.110436814509468826540991;
    x[119] = -0.098208148184447540736015;
    x[120] = -0.085964634131980604256000;
    x[121] = -0.073708123403767780288977;
    x[122] = -0.061440469016428270850728;
    x[123] = -0.049163525671349973093019;
    x[124] = -0.036879149474284021657652;
    x[125] = -0.024589197654727010541405;
    x[126] = -0.012295528285133320036860;
    x[127] = 0.000000000000000000000000;
    x[128] = 0.012295528285133320036860;
    x[129] = 0.024589197654727010541405;
    x[130] = 0.036879149474284021657652;
    x[131] = 0.049163525671349973093019;
    x[132] = 0.061440469016428270850728;
    x[133] = 0.073708123403767780288977;
    x[134] = 0.085964634131980604256000;
    x[135] = 0.098208148184447540736015;
    x[136] = 0.110436814509468826540991;
    x[137] = 0.122648784300117812092492;
    x[138] = 0.134842211273755257250625;
    x[139] = 0.147015251951161989456661;
    x[140] = 0.159166065935247723154292;
    x[141] = 0.171292816189293903533225;
    x[142] = 0.183393669314688508087976;
    x[143] = 0.195466795828110816293869;
    x[144] = 0.207510370438124240859625;
    x[145] = 0.219522572321135403508985;
    x[146] = 0.231501585396677734059116;
    x[147] = 0.243445598601977973686482;
    x[148] = 0.255352806165764071686080;
    x[149] = 0.267221407881273079721012;
    x[150] = 0.279049609378417768508970;
    x[151] = 0.290835622395070819082047;
    x[152] = 0.302577665047425574167140;
    x[153] = 0.314273962099392474845918;
    x[154] = 0.325922745230990453444769;
    x[155] = 0.337522253305692705554261;
    x[156] = 0.349070732636686422161576;
    x[157] = 0.360566437252006227074021;
    x[158] = 0.372007629158501235092510;
    x[159] = 0.383392578604595822734854;
    x[160] = 0.394719564341804385683361;
    x[161] = 0.405986873884960545511889;
    x[162] = 0.417192803771121462605751;
    x[163] = 0.428335659817108112494341;
    x[164] = 0.439413757375642589040685;
    x[165] = 0.450425421590043710043279;
    x[166] = 0.461368987647442418771401;
    x[167] = 0.472242801030478698742627;
    x[168] = 0.483045217767441948626854;
    x[169] = 0.493774604680816999489812;
    x[170] = 0.504429339634198197635551;
    x[171] = 0.515007811777534223035005;
    x[172] = 0.525508421790666565699453;
    x[173] = 0.535929582125124840335150;
    x[174] = 0.546269717244142383159817;
    x[175] = 0.556527263860855843833077;
    x[176] = 0.566700671174652760010815;
    x[177] = 0.576788401105631382036211;
    x[178] = 0.586788928527137300685706;
    x[179] = 0.596700741496341721653202;
    x[180] = 0.606522341482826526536576;
    x[181] = 0.616252243595141561442344;
    x[182] = 0.625888976805299900901619;
    x[183] = 0.635431084171177146547142;
    x[184] = 0.644877123056781136890077;
    x[185] = 0.654225665350358766508700;
    x[186] = 0.663475297680306939970658;
    x[187] = 0.672624621628855017806731;
    x[188] = 0.681672253943486448787259;
    x[189] = 0.690616826746067624571761;
    x[190] = 0.699456987739652339456557;
    x[191] = 0.708191400412930589382399;
    x[192] = 0.716818744242290800531501;
    x[193] = 0.725337714891464938687812;
    x[194] = 0.733747024408726316001889;
    x[195] = 0.742045401421610281838045;
    x[196] = 0.750231591329128358931528;
    x[197] = 0.758304356491446765092016;
    x[198] = 0.766262476417000644100858;
    x[199] = 0.774104747947015717207115;
    x[200] = 0.781829985437409458675147;
    x[201] = 0.789437020938044295117764;
    x[202] = 0.796924704369305728807154;
    x[203] = 0.804291903695978689734633;
    x[204] = 0.811537505098395829833580;
    x[205] = 0.818660413140831885432406;
    x[206] = 0.825659550937118650611534;
    x[207] = 0.832533860313455524647230;
    x[208] = 0.839282301968391021084600;
    x[209] = 0.845903855629951054143931;
    x[210] = 0.852397520209890250084237;
    x[211] = 0.858762313955042966785823;
    x[212] = 0.864997274595751144137121;
    x[213] = 0.871101459491346550796200;
    x[214] = 0.877073945772665439532627;
    x[215] = 0.882913830481574073645470;
    x[216] = 0.888620230707484040924981;
    x[217] = 0.894192283720836729335637;
    x[218] = 0.899629147103536800144342;
    x[219] = 0.904929998876314959753358;
    x[220] = 0.910094037623000801254172;
    x[221] = 0.915120482611686961035103;
    x[222] = 0.920008573912766315142721;
    x[223] = 0.924757572513824425220425;
    x[224] = 0.929366760431369935739045;
    x[225] = 0.933835440819386124349338;
    x[226] = 0.938162938074687317626793;
    x[227] = 0.942348597939064408301480;
    x[228] = 0.946391787598204251752103;
    x[229] = 0.950291895777368285733522;
    x[230] = 0.954048332833816317950921;
    x[231] = 0.957660530845962076295490;
    x[232] = 0.961127943699247839572910;
    x[233] = 0.964450047168726298761719;
    x[234] = 0.967626338998338798105523;
    x[235] = 0.970656338976880365477697;
    x[236] = 0.973539589010643617645393;
    x[237] = 0.976275653192735980815246;
    x[238] = 0.978864117869068155239121;
    x[239] = 0.981304591701017185126565;
    x[240] = 0.983596705724776358640192;
    x[241] = 0.985740113407419277752900;
    x[242] = 0.987734490699732356281248;
    x[243] = 0.989579536085920123498574;
    x[244] = 0.991274970630385567164523;
    x[245] = 0.992820538021989138984811;
    x[246] = 0.994216004616630164799381;
    x[247] = 0.995461159480026294089975;
    x[248] = 0.996555814435198617028738;
    x[249] = 0.997499804126615814044844;
    x[250] = 0.998292986136967889228248;
    x[251] = 0.998935241284654635142155;
    x[252] = 0.999426474680169959344386;
    x[253] = 0.999766621312000569367063;
    x[254] = 0.999955705317563751730191;

    w[0] = 0.00011367361999142272115645954414;
    w[1] = 0.00026459387119083065532790838855;
    w[2] = 0.00041569762526823913616284210066;
    w[3] = 0.00056675794564824918946626058353;
    w[4] = 0.00071773647800611087798371518325;
    w[5] = 0.00086860766611945667949717690640;
    w[6] = 0.00101934797642732530281229369360;
    w[7] = 0.0011699343729388079886897709773;
    w[8] = 0.0013203439900221692090523602144;
    w[9] = 0.0014705540427783843160097204304;
    w[10] = 0.0016205417990415653896921100325;
    w[11] = 0.0017702845706603213070421243905;
    w[12] = 0.0019197597117132050055085980675;
    w[13] = 0.0020689446195015801533643667413;
    w[14] = 0.0022178167367540171700373764020;
    w[15] = 0.0023663535543962867157201855305;
    w[16] = 0.0025145326145997073931298921370;
    w[17] = 0.0026623315139717112732749157331;
    w[18] = 0.0028097279068204407457332299361;
    w[19] = 0.0029566995084575002760043344138;
    w[20] = 0.0031032240985191112621977893133;
    w[21] = 0.0032492795242943133198690930777;
    w[22] = 0.0033948437040533928255056951665;
    w[23] = 0.0035398946303722552150296713510;
    w[24] = 0.0036844103734499176530742235517;
    w[25] = 0.0038283690844171626400743524999;
    w[26] = 0.0039717489986349171988699773906;
    w[27] = 0.0041145284389812475901826468094;
    w[28] = 0.0042566858191260658425395494472;
    w[29] = 0.0043981996467927779838546384780;
    w[30] = 0.0045390485270061921259394035112;
    w[31] = 0.0046792111653260640506279893190;
    w[32] = 0.0048186663710656988918572043815;
    w[33] = 0.0049573930604950563104281084148;
    w[34] = 0.0050953702600278273039420404117;
    w[35] = 0.0052325771093919661294970523234;
    w[36] = 0.0053689928647831724787741258653;
    w[37] = 0.0055045969020008281904902120813;
    w[38] = 0.0056393687195659001929970994675;
    w[39] = 0.0057732879418203275712033691864;
    w[40] = 0.0059063343220074160130475409466;
    w[41] = 0.0060384877453327676663371666884;
    w[42] = 0.0061697282320052788060812561217;
    w[43] = 0.0063000359402577418025981070425;
    w[44] = 0.0064293911693465917826140832500;
    w[45] = 0.0065577743625303421548456356354;
    w[46] = 0.0066851661100262568757892743568;
    w[47] = 0.0068115471519448109954345674817;
    w[48] = 0.0069368983812014946719507501243;
    w[49] = 0.0070612008464055194979848418291;
    w[50] = 0.0071844357547249896530757997058;
    w[51] = 0.0073065844747281040972736443146;
    w[52] = 0.0074276285391999597581348419714;
    w[53] = 0.0075475496479345294426435656724;
    w[54] = 0.0076663296705013920315933272426;
    w[55] = 0.0077839506489867963897419914623;
    w[56] = 0.0079003948007086443529587296692;
    w[57] = 0.0080156445209049821352946484008;
    w[58] = 0.0081296823853955935356080649925;
    w[59] = 0.0082424911532162924158504385939;
    w[60] = 0.0083540537692255160718568405530;
    w[61] = 0.0084643533666828253227353760036;
    w[62] = 0.0085733732697989214067758505840;
    w[63] = 0.0086810969962567940901133439612;
    w[64] = 0.0087875082597036197689825483144;
    w[65] = 0.0088925909722130327769834298578;
    w[66] = 0.0089963292467173975949700110383;
    w[67] = 0.0090987073994097142025303711406;
    w[68] = 0.0091997099521147934060534414075;
    w[69] = 0.0092993216346293436285393234867;
    w[70] = 0.0093975273870306153500305317074;
    w[71] = 0.0094943123619532541442165010292;
    w[72] = 0.0095896619268340180657610209655;
    w[73] = 0.0096835616661240200035669970076;
    w[74] = 0.0097759973834681605268499842249;
    w[75] = 0.0098669551038514217128483481814;
    w[76] = 0.0099564210757116974565448593910;
    w[77] = 0.0100443817730188408231888789497;
    w[78] = 0.0101308238973196141129538950955;
    w[79] = 0.0102157343797482324629939488415;
    w[80] = 0.0102991003830021970147153502911;
    w[81] = 0.0103809093032831189224876935085;
    w[82] = 0.0104611487722022407735015844669;
    w[83] = 0.0105398066586503673262517188088;
    w[84] = 0.0106168710706319228563864391054;
    w[85] = 0.0106923303570628578226139809571;
    w[86] = 0.0107661731095321330311788312990;
    w[87] = 0.0108383881640265149842990798832;
    w[88] = 0.0109089646026184216450603134401;
    w[89] = 0.0109778917551165634377595759712;
    w[90] = 0.0110451592006791299277436662993;
    w[91] = 0.0111107567693892782875426356195;
    w[92] = 0.0111746745437926853557086684962;
    w[93] = 0.0112369028603969308303734810332;
    w[94] = 0.0112974323111324849102690558722;
    w[95] = 0.0113562537447750795009464486204;
    w[96] = 0.011413358268329247942299599697;
    w[97] = 0.011468737248372824084374355981;
    w[98] = 0.011522382312362197440930930031;
    w[99] = 0.011574285349898127083439539046;
    w[100] = 0.011624438513951922901227922331;
    w[101] = 0.011672834222051808845465154244;
    w[102] = 0.011719465157429288794653489478;
    w[103] = 0.011764324270125341726399410909;
    w[104] = 0.011807404778056278953532930501;
    w[105] = 0.011848700168039102281222824051;
    w[106] = 0.011888204196776208064673282076;
    w[107] = 0.011925910891799288293359117699;
    w[108] = 0.011961814552372285996633285380;
    w[109] = 0.011995909750353268455989686823;
    w[110] = 0.012028191331015087920350431142;
    w[111] = 0.012058654413824705751531083631;
    w[112] = 0.012087294393181062176578184854;
    w[113] = 0.012114106939111380091025793650;
    w[114] = 0.012139087997925797641334635250;
    w[115] = 0.012162233792830230614908682534;
    w[116] = 0.012183540824497371981177306326;
    w[117] = 0.012203005871595742256331865516;
    w[118] = 0.012220625991276710706457005806;
    w[119] = 0.012236398519619413758040249691;
    w[120] = 0.012250321072033503350218104906;
    w[121] = 0.012262391543619664338660618398;
    w[122] = 0.012272608109487846445745237751;
    w[123] = 0.012280969225033162644659793962;
    w[124] = 0.012287473626169412265336919908;
    w[125] = 0.012292120329520193516690694701;
    w[126] = 0.012294908632567576531532225710;
    w[127] = 0.01229583811375831445681490730;
    w[128] = 0.012294908632567576531532225710;
    w[129] = 0.012292120329520193516690694701;
    w[130] = 0.012287473626169412265336919908;
    w[131] = 0.012280969225033162644659793962;
    w[132] = 0.012272608109487846445745237751;
    w[133] = 0.012262391543619664338660618398;
    w[134] = 0.012250321072033503350218104906;
    w[135] = 0.012236398519619413758040249691;
    w[136] = 0.012220625991276710706457005806;
    w[137] = 0.012203005871595742256331865516;
    w[138] = 0.012183540824497371981177306326;
    w[139] = 0.012162233792830230614908682534;
    w[140] = 0.012139087997925797641334635250;
    w[141] = 0.012114106939111380091025793650;
    w[142] = 0.012087294393181062176578184854;
    w[143] = 0.012058654413824705751531083631;
    w[144] = 0.012028191331015087920350431142;
    w[145] = 0.011995909750353268455989686823;
    w[146] = 0.011961814552372285996633285380;
    w[147] = 0.011925910891799288293359117699;
    w[148] = 0.011888204196776208064673282076;
    w[149] = 0.011848700168039102281222824051;
    w[150] = 0.011807404778056278953532930501;
    w[151] = 0.011764324270125341726399410909;
    w[152] = 0.011719465157429288794653489478;
    w[153] = 0.011672834222051808845465154244;
    w[154] = 0.011624438513951922901227922331;
    w[155] = 0.011574285349898127083439539046;
    w[156] = 0.011522382312362197440930930031;
    w[157] = 0.011468737248372824084374355981;
    w[158] = 0.011413358268329247942299599697;
    w[159] = 0.0113562537447750795009464486204;
    w[160] = 0.0112974323111324849102690558722;
    w[161] = 0.0112369028603969308303734810332;
    w[162] = 0.0111746745437926853557086684962;
    w[163] = 0.0111107567693892782875426356195;
    w[164] = 0.0110451592006791299277436662993;
    w[165] = 0.0109778917551165634377595759712;
    w[166] = 0.0109089646026184216450603134401;
    w[167] = 0.0108383881640265149842990798832;
    w[168] = 0.0107661731095321330311788312990;
    w[169] = 0.0106923303570628578226139809571;
    w[170] = 0.0106168710706319228563864391054;
    w[171] = 0.0105398066586503673262517188088;
    w[172] = 0.0104611487722022407735015844669;
    w[173] = 0.0103809093032831189224876935085;
    w[174] = 0.0102991003830021970147153502911;
    w[175] = 0.0102157343797482324629939488415;
    w[176] = 0.0101308238973196141129538950955;
    w[177] = 0.0100443817730188408231888789497;
    w[178] = 0.0099564210757116974565448593910;
    w[179] = 0.0098669551038514217128483481814;
    w[180] = 0.0097759973834681605268499842249;
    w[181] = 0.0096835616661240200035669970076;
    w[182] = 0.0095896619268340180657610209655;
    w[183] = 0.0094943123619532541442165010292;
    w[184] = 0.0093975273870306153500305317074;
    w[185] = 0.0092993216346293436285393234867;
    w[186] = 0.0091997099521147934060534414075;
    w[187] = 0.0090987073994097142025303711406;
    w[188] = 0.0089963292467173975949700110383;
    w[189] = 0.0088925909722130327769834298578;
    w[190] = 0.0087875082597036197689825483144;
    w[191] = 0.0086810969962567940901133439612;
    w[192] = 0.0085733732697989214067758505840;
    w[193] = 0.0084643533666828253227353760036;
    w[194] = 0.0083540537692255160718568405530;
    w[195] = 0.0082424911532162924158504385939;
    w[196] = 0.0081296823853955935356080649925;
    w[197] = 0.0080156445209049821352946484008;
    w[198] = 0.0079003948007086443529587296692;
    w[199] = 0.0077839506489867963897419914623;
    w[200] = 0.0076663296705013920315933272426;
    w[201] = 0.0075475496479345294426435656724;
    w[202] = 0.0074276285391999597581348419714;
    w[203] = 0.0073065844747281040972736443146;
    w[204] = 0.0071844357547249896530757997058;
    w[205] = 0.0070612008464055194979848418291;
    w[206] = 0.0069368983812014946719507501243;
    w[207] = 0.0068115471519448109954345674817;
    w[208] = 0.0066851661100262568757892743568;
    w[209] = 0.0065577743625303421548456356354;
    w[210] = 0.0064293911693465917826140832500;
    w[211] = 0.0063000359402577418025981070425;
    w[212] = 0.0061697282320052788060812561217;
    w[213] = 0.0060384877453327676663371666884;
    w[214] = 0.0059063343220074160130475409466;
    w[215] = 0.0057732879418203275712033691864;
    w[216] = 0.0056393687195659001929970994675;
    w[217] = 0.0055045969020008281904902120813;
    w[218] = 0.0053689928647831724787741258653;
    w[219] = 0.0052325771093919661294970523234;
    w[220] = 0.0050953702600278273039420404117;
    w[221] = 0.0049573930604950563104281084148;
    w[222] = 0.0048186663710656988918572043815;
    w[223] = 0.0046792111653260640506279893190;
    w[224] = 0.0045390485270061921259394035112;
    w[225] = 0.0043981996467927779838546384780;
    w[226] = 0.0042566858191260658425395494472;
    w[227] = 0.0041145284389812475901826468094;
    w[228] = 0.0039717489986349171988699773906;
    w[229] = 0.0038283690844171626400743524999;
    w[230] = 0.0036844103734499176530742235517;
    w[231] = 0.0035398946303722552150296713510;
    w[232] = 0.0033948437040533928255056951665;
    w[233] = 0.0032492795242943133198690930777;
    w[234] = 0.0031032240985191112621977893133;
    w[235] = 0.0029566995084575002760043344138;
    w[236] = 0.0028097279068204407457332299361;
    w[237] = 0.0026623315139717112732749157331;
    w[238] = 0.0025145326145997073931298921370;
    w[239] = 0.0023663535543962867157201855305;
    w[240] = 0.0022178167367540171700373764020;
    w[241] = 0.0020689446195015801533643667413;
    w[242] = 0.0019197597117132050055085980675;
    w[243] = 0.0017702845706603213070421243905;
    w[244] = 0.0016205417990415653896921100325;
    w[245] = 0.0014705540427783843160097204304;
    w[246] = 0.0013203439900221692090523602144;
    w[247] = 0.0011699343729388079886897709773;
    w[248] = 0.00101934797642732530281229369360;
    w[249] = 0.00086860766611945667949717690640;
    w[250] = 0.00071773647800611087798371518325;
    w[251] = 0.00056675794564824918946626058353;
    w[252] = 0.00041569762526823913616284210066;
    w[253] = 0.00026459387119083065532790838855;
    w[254] = 0.00011367361999142272115645954414;
  }
  else if ( n == 256 )
  {
    x[0] = -0.999956050018992230734801;
    x[1] = -0.999768437409263186104879;
    x[2] = -0.999430937466261408240854;
    x[3] = -0.998943525843408856555026;
    x[4] = -0.998306266473006444055500;
    x[5] = -0.997519252756720827563409;
    x[6] = -0.996582602023381540430504;
    x[7] = -0.995496454481096356592647;
    x[8] = -0.994260972922409664962878;
    x[9] = -0.992876342608822117143534;
    x[10] = -0.991342771207583086922189;
    x[11] = -0.989660488745065218319244;
    x[12] = -0.987829747564860608916488;
    x[13] = -0.985850822286125956479245;
    x[14] = -0.983724009760315496166686;
    x[15] = -0.981449629025464405769303;
    x[16] = -0.979028021257622038824238;
    x[17] = -0.976459549719234155621011;
    x[18] = -0.973744599704370405266079;
    x[19] = -0.970883578480743029320923;
    x[20] = -0.967876915228489454909004;
    x[21] = -0.964725060975706430932612;
    x[22] = -0.961428488530732144006407;
    x[23] = -0.957987692411178129365790;
    x[24] = -0.954403188769716241764448;
    x[25] = -0.950675515316628276363852;
    x[26] = -0.946805231239127481372052;
    x[27] = -0.942792917117462443183076;
    x[28] = -0.938639174837814804981926;
    x[29] = -0.934344627502003094292477;
    x[30] = -0.929909919334005641180246;
    x[31] = -0.925335715583316202872730;
    x[32] = -0.920622702425146495505047;
    x[33] = -0.915771586857490384526670;
    x[34] = -0.910783096595065011890907;
    x[35] = -0.905657979960144647082682;
    x[36] = -0.900397005770303544771620;
    x[37] = -0.895000963223084577441223;
    x[38] = -0.889470661777610888828677;
    x[39] = -0.883806931033158284859826;
    x[40] = -0.878010620604706543986435;
    x[41] = -0.872082599995488289130046;
    x[42] = -0.866023758466554519297515;
    x[43] = -0.859835004903376350696173;
    x[44] = -0.853517267679502965073036;
    x[45] = -0.847071494517296207187072;
    x[46] = -0.840498652345762713895068;
    x[47] = -0.833799727155504894348444;
    x[48] = -0.826975723850812514289093;
    x[49] = -0.820027666098917067403478;
    x[50] = -0.812956596176431543136410;
    x[51] = -0.805763574812998623257389;
    x[52] = -0.798449681032170758782543;
    x[53] = -0.791016011989545994546707;
    x[54] = -0.783463682808183820750670;
    x[55] = -0.775793826411325739132053;
    x[56] = -0.768007593352445635975891;
    x[57] = -0.760106151642655454941907;
    x[58] = -0.752090686575492059587530;
    x[59] = -0.743962400549111568455683;
    x[60] = -0.735722512885917834620373;
    x[61] = -0.727372259649652126586894;
    x[62] = -0.718912893459971448372640;
    x[63] = -0.710345683304543313394566;
    x[64] = -0.701671914348685159406084;
    x[65] = -0.692892887742576960105342;
    x[66] = -0.684009920426075953124877;
    x[67] = -0.675024344931162763855919;
    x[68] = -0.665937509182048559906408;
    x[69] = -0.656750776292973221887500;
    x[70] = -0.647465524363724862617016;
    x[71] = -0.638083146272911368668689;
    x[72] = -0.628605049469014975432210;
    x[73] = -0.619032655759261219430968;
    x[74] = -0.609367401096333939522311;
    x[75] = -0.599610735362968321730388;
    x[76] = -0.589764122154454300785786;
    x[77] = -0.579829038559082944921832;
    x[78] = -0.569806974936568759057668;
    x[79] = -0.559699434694481145136907;
    x[80] = -0.549507934062718557042427;
    x[81] = -0.539234001866059181127936;
    x[82] = -0.528879179294822261951476;
    x[83] = -0.518445019673674476221662;
    x[84] = -0.507933088228616036231925;
    x[85] = -0.497344961852181477119512;
    x[86] = -0.486682228866890350103621;
    x[87] = -0.475946488786983306390738;
    x[88] = -0.465139352078479313645570;
    x[89] = -0.454262439917589998774455;
    x[90] = -0.443317383947527357216926;
    x[91] = -0.432305826033741309953441;
    x[92] = -0.421229418017623824976812;
    x[93] = -0.410089821468716550006434;
    x[94] = -0.398888707435459127713463;
    x[95] = -0.387627756194515583637985;
    x[96] = -0.376308656998716390283056;
    x[97] = -0.364933107823654018533465;
    x[98] = -0.353502815112969989537790;
    x[99] = -0.342019493522371636480730;
    x[100] = -0.330484865662416976229187;
    x[101] = -0.318900661840106275631683;
    x[102] = -0.307268619799319076258610;
    x[103] = -0.295590484460135614563787;
    x[104] = -0.283868007657081741799766;
    x[105] = -0.272102947876336609505245;
    x[106] = -0.260297069991942541978561;
    x[107] = -0.248452145001056666833243;
    x[108] = -0.236569949758284018477508;
    x[109] = -0.224652266709131967147878;
    x[110] = -0.212700883622625957937040;
    x[111] = -0.200717593323126670068001;
    x[112] = -0.188704193421388826461504;
    x[113] = -0.176662486044901997403722;
    x[114] = -0.164594277567553849829285;
    x[115] = -0.152501378338656395374607;
    x[116] = -0.140385602411375885913025;
    x[117] = -0.128248767270607094742050;
    x[118] = -0.116092693560332804940735;
    x[119] = -0.103919204810509403639197;
    x[120] = -0.091730127163519552031146;
    x[121] = -0.079527289100232965903227;
    x[122] = -0.067312521165716400242290;
    x[123] = -0.055087655694633984104561;
    x[124] = -0.042854526536379098381242;
    x[125] = -0.030614968779979029366279;
    x[126] = -0.018370818478813665117926;
    x[127] = -0.006123912375189529501170;
    x[128] = 0.006123912375189529501170;
    x[129] = 0.018370818478813665117926;
    x[130] = 0.030614968779979029366279;
    x[131] = 0.042854526536379098381242;
    x[132] = 0.055087655694633984104561;
    x[133] = 0.067312521165716400242290;
    x[134] = 0.079527289100232965903227;
    x[135] = 0.091730127163519552031146;
    x[136] = 0.103919204810509403639197;
    x[137] = 0.116092693560332804940735;
    x[138] = 0.128248767270607094742050;
    x[139] = 0.140385602411375885913025;
    x[140] = 0.152501378338656395374607;
    x[141] = 0.164594277567553849829285;
    x[142] = 0.176662486044901997403722;
    x[143] = 0.188704193421388826461504;
    x[144] = 0.200717593323126670068001;
    x[145] = 0.212700883622625957937040;
    x[146] = 0.224652266709131967147878;
    x[147] = 0.236569949758284018477508;
    x[148] = 0.248452145001056666833243;
    x[149] = 0.260297069991942541978561;
    x[150] = 0.272102947876336609505245;
    x[151] = 0.283868007657081741799766;
    x[152] = 0.295590484460135614563787;
    x[153] = 0.307268619799319076258610;
    x[154] = 0.318900661840106275631683;
    x[155] = 0.330484865662416976229187;
    x[156] = 0.342019493522371636480730;
    x[157] = 0.353502815112969989537790;
    x[158] = 0.364933107823654018533465;
    x[159] = 0.376308656998716390283056;
    x[160] = 0.387627756194515583637985;
    x[161] = 0.398888707435459127713463;
    x[162] = 0.410089821468716550006434;
    x[163] = 0.421229418017623824976812;
    x[164] = 0.432305826033741309953441;
    x[165] = 0.443317383947527357216926;
    x[166] = 0.454262439917589998774455;
    x[167] = 0.465139352078479313645570;
    x[168] = 0.475946488786983306390738;
    x[169] = 0.486682228866890350103621;
    x[170] = 0.497344961852181477119512;
    x[171] = 0.507933088228616036231925;
    x[172] = 0.518445019673674476221662;
    x[173] = 0.528879179294822261951476;
    x[174] = 0.539234001866059181127936;
    x[175] = 0.549507934062718557042427;
    x[176] = 0.559699434694481145136907;
    x[177] = 0.569806974936568759057668;
    x[178] = 0.579829038559082944921832;
    x[179] = 0.589764122154454300785786;
    x[180] = 0.599610735362968321730388;
    x[181] = 0.609367401096333939522311;
    x[182] = 0.619032655759261219430968;
    x[183] = 0.628605049469014975432210;
    x[184] = 0.638083146272911368668689;
    x[185] = 0.647465524363724862617016;
    x[186] = 0.656750776292973221887500;
    x[187] = 0.665937509182048559906408;
    x[188] = 0.675024344931162763855919;
    x[189] = 0.684009920426075953124877;
    x[190] = 0.692892887742576960105342;
    x[191] = 0.701671914348685159406084;
    x[192] = 0.710345683304543313394566;
    x[193] = 0.718912893459971448372640;
    x[194] = 0.727372259649652126586894;
    x[195] = 0.735722512885917834620373;
    x[196] = 0.743962400549111568455683;
    x[197] = 0.752090686575492059587530;
    x[198] = 0.760106151642655454941907;
    x[199] = 0.768007593352445635975891;
    x[200] = 0.775793826411325739132053;
    x[201] = 0.783463682808183820750670;
    x[202] = 0.791016011989545994546707;
    x[203] = 0.798449681032170758782543;
    x[204] = 0.805763574812998623257389;
    x[205] = 0.812956596176431543136410;
    x[206] = 0.820027666098917067403478;
    x[207] = 0.826975723850812514289093;
    x[208] = 0.833799727155504894348444;
    x[209] = 0.840498652345762713895068;
    x[210] = 0.847071494517296207187072;
    x[211] = 0.853517267679502965073036;
    x[212] = 0.859835004903376350696173;
    x[213] = 0.866023758466554519297515;
    x[214] = 0.872082599995488289130046;
    x[215] = 0.878010620604706543986435;
    x[216] = 0.883806931033158284859826;
    x[217] = 0.889470661777610888828677;
    x[218] = 0.895000963223084577441223;
    x[219] = 0.900397005770303544771620;
    x[220] = 0.905657979960144647082682;
    x[221] = 0.910783096595065011890907;
    x[222] = 0.915771586857490384526670;
    x[223] = 0.920622702425146495505047;
    x[224] = 0.925335715583316202872730;
    x[225] = 0.929909919334005641180246;
    x[226] = 0.934344627502003094292477;
    x[227] = 0.938639174837814804981926;
    x[228] = 0.942792917117462443183076;
    x[229] = 0.946805231239127481372052;
    x[230] = 0.950675515316628276363852;
    x[231] = 0.954403188769716241764448;
    x[232] = 0.957987692411178129365790;
    x[233] = 0.961428488530732144006407;
    x[234] = 0.964725060975706430932612;
    x[235] = 0.967876915228489454909004;
    x[236] = 0.970883578480743029320923;
    x[237] = 0.973744599704370405266079;
    x[238] = 0.976459549719234155621011;
    x[239] = 0.979028021257622038824238;
    x[240] = 0.981449629025464405769303;
    x[241] = 0.983724009760315496166686;
    x[242] = 0.985850822286125956479245;
    x[243] = 0.987829747564860608916488;
    x[244] = 0.989660488745065218319244;
    x[245] = 0.991342771207583086922189;
    x[246] = 0.992876342608822117143534;
    x[247] = 0.994260972922409664962878;
    x[248] = 0.995496454481096356592647;
    x[249] = 0.996582602023381540430504;
    x[250] = 0.997519252756720827563409;
    x[251] = 0.998306266473006444055500;
    x[252] = 0.998943525843408856555026;
    x[253] = 0.999430937466261408240854;
    x[254] = 0.999768437409263186104879;
    x[255] = 0.999956050018992230734801;

    w[0] = 0.00011278901782227217551253887725;
    w[1] = 0.00026253494429644590628745756250;
    w[2] = 0.00041246325442617632843218583774;
    w[3] = 0.00056234895403140980281523674759;
    w[4] = 0.0007121541634733206669089891511;
    w[5] = 0.0008618537014200890378140934163;
    w[6] = 0.0010114243932084404526058128414;
    w[7] = 0.0011608435575677247239705981135;
    w[8] = 0.0013100886819025044578316804271;
    w[9] = 0.0014591373333107332010883864996;
    w[10] = 0.0016079671307493272424499395690;
    w[11] = 0.0017565557363307299936069145295;
    w[12] = 0.0019048808534997184044191411746;
    w[13] = 0.0020529202279661431745487818492;
    w[14] = 0.0022006516498399104996848834189;
    w[15] = 0.0023480529563273120170064609087;
    w[16] = 0.0024951020347037068508395354372;
    w[17] = 0.0026417768254274905641208292516;
    w[18] = 0.0027880553253277068805747610763;
    w[19] = 0.0029339155908297166460123254142;
    w[20] = 0.0030793357411993375832053528316;
    w[21] = 0.0032242939617941981570107134269;
    w[22] = 0.0033687685073155510120191062489;
    w[23] = 0.0035127377050563073309710549844;
    w[24] = 0.0036561799581425021693892413052;
    w[25] = 0.0037990737487662579981170192082;
    w[26] = 0.0039413976414088336277290349840;
    w[27] = 0.0040831302860526684085997759212;
    w[28] = 0.0042242504213815362723565049060;
    w[29] = 0.0043647368779680566815684200621;
    w[30] = 0.0045045685814478970686417923159;
    w[31] = 0.0046437245556800603139790923525;
    w[32] = 0.0047821839258926913729317340448;
    w[33] = 0.0049199259218138656695587765655;
    w[34] = 0.0050569298807868423875578160762;
    w[35] = 0.0051931752508692809303287536296;
    w[36] = 0.0053286415939159303170811114788;
    w[37] = 0.0054633085886443102775705318566;
    w[38] = 0.0055971560336829100775514452572;
    w[39] = 0.005730163850601437177384417555;
    w[40] = 0.005862312086922653060661598801;
    w[41] = 0.005993580919115338221127696870;
    w[42] = 0.006123950655567932542389081187;
    w[43] = 0.006253401739542401272063645975;
    w[44] = 0.006381914752107880570375164275;
    w[45] = 0.006509470415053660267809899951;
    w[46] = 0.006636049593781065044590038355;
    w[47] = 0.006761633300173798780927861108;
    w[48] = 0.006886202695446320346713323775;
    w[49] = 0.007009739092969822621234436194;
    w[50] = 0.007132223961075390071672422986;
    w[51] = 0.007253638925833913783829137214;
    w[52] = 0.007373965773812346437572440695;
    w[53] = 0.007493186454805883358599761133;
    w[54] = 0.007611283084545659461618719618;
    w[55] = 0.007728237947381555631110194958;
    w[56] = 0.007844033498939711866810316151;
    w[57] = 0.007958652368754348353613161227;
    w[58] = 0.008072077362873499500946974804;
    w[59] = 0.008184291466438269935619761004;
    w[60] = 0.008295277846235225425171412553;
    w[61] = 0.008405019853221535756180301698;
    w[62] = 0.008513501025022490693838354790;
    w[63] = 0.008620705088401014305368838410;
    w[64] = 0.008726615961698807140336632217;
    w[65] = 0.008831217757248750025318272685;
    w[66] = 0.008934494783758207548408417085;
    w[67] = 0.009036431548662873680227775572;
    w[68] = 0.009137012760450806402000472219;
    w[69] = 0.009236223330956302687378716714;
    w[70] = 0.009334048377623269712466014486;
    w[71] = 0.009430473225737752747352764482;
    w[72] = 0.009525483410629284811829685754;
    w[73] = 0.009619064679840727857162164401;
    w[74] = 0.009711202995266279964249670496;
    w[75] = 0.009801884535257327825498800250;
    w[76] = 0.009891095696695828602630683809;
    w[77] = 0.009978823097034910124733949495;
    w[78] = 0.010065053576306383309460978930;
    w[79] = 0.010149774199094865654634066042;
    w[80] = 0.010232972256478219656954857160;
    w[81] = 0.010314635267934015068260713997;
    w[82] = 0.010394750983211728997101725205;
    w[83] = 0.010473307384170403003569566927;
    w[84] = 0.010550292686581481517533575536;
    w[85] = 0.010625695341896561133961681801;
    w[86] = 0.010699504038979785603048200583;
    w[87] = 0.010771707705804626636653631927;
    w[88] = 0.010842295511114795995293477058;
    w[89] = 0.010911256866049039700796847788;
    w[90] = 0.010978581425729570637988203448;
    w[91] = 0.011044259090813901263517571044;
    w[92] = 0.011108280009009843630460815451;
    w[93] = 0.011170634576553449462710881938;
    w[94] = 0.011231313439649668572656802083;
    w[95] = 0.011290307495875509508367594121;
    w[96] = 0.011347607895545491941625714297;
    w[97] = 0.011403206043039185964847059552;
    w[98] = 0.011457093598090639152334392298;
    w[99] = 0.011509262477039497958586392439;
    w[100] = 0.011559704854043635772668656950;
    w[101] = 0.011608413162253105722084706677;
    w[102] = 0.011655380094945242121298939730;
    w[103] = 0.011700598606620740288189823359;
    w[104] = 0.011744061914060550305376732759;
    w[105] = 0.011785763497343426181690117627;
    w[106] = 0.011825697100823977771160737958;
    w[107] = 0.011863856734071078731904572908;
    w[108] = 0.011900236672766489754287204237;
    w[109] = 0.011934831459563562255873201696;
    w[110] = 0.011967635904905893729007282670;
    w[111] = 0.011998645087805811934536710071;
    w[112] = 0.012027854356582571161267533498;
    w[113] = 0.012055259329560149814347085327;
    w[114] = 0.012080855895724544655975183976;
    w[115] = 0.012104640215340463097757829736;
    w[116] = 0.012126608720527321034718492205;
    w[117] = 0.012146758115794459815559837664;
    w[118] = 0.012165085378535502061307291839;
    w[119] = 0.012181587759481772174047585032;
    w[120] = 0.012196262783114713518180974196;
    w[121] = 0.012209108248037240407514094371;
    w[122] = 0.012220122227303969191708737227;
    w[123] = 0.012229303068710278904146266083;
    w[124] = 0.012236649395040158109242574767;
    w[125] = 0.012242160104272800769728083260;
    w[126] = 0.012245834369747920142463857550;
    w[127] = 0.01224767164028975590407032649;
    w[128] = 0.01224767164028975590407032649;
    w[129] = 0.012245834369747920142463857550;
    w[130] = 0.012242160104272800769728083260;
    w[131] = 0.012236649395040158109242574767;
    w[132] = 0.012229303068710278904146266083;
    w[133] = 0.012220122227303969191708737227;
    w[134] = 0.012209108248037240407514094371;
    w[135] = 0.012196262783114713518180974196;
    w[136] = 0.012181587759481772174047585032;
    w[137] = 0.012165085378535502061307291839;
    w[138] = 0.012146758115794459815559837664;
    w[139] = 0.012126608720527321034718492205;
    w[140] = 0.012104640215340463097757829736;
    w[141] = 0.012080855895724544655975183976;
    w[142] = 0.012055259329560149814347085327;
    w[143] = 0.012027854356582571161267533498;
    w[144] = 0.011998645087805811934536710071;
    w[145] = 0.011967635904905893729007282670;
    w[146] = 0.011934831459563562255873201696;
    w[147] = 0.011900236672766489754287204237;
    w[148] = 0.011863856734071078731904572908;
    w[149] = 0.011825697100823977771160737958;
    w[150] = 0.011785763497343426181690117627;
    w[151] = 0.011744061914060550305376732759;
    w[152] = 0.011700598606620740288189823359;
    w[153] = 0.011655380094945242121298939730;
    w[154] = 0.011608413162253105722084706677;
    w[155] = 0.011559704854043635772668656950;
    w[156] = 0.011509262477039497958586392439;
    w[157] = 0.011457093598090639152334392298;
    w[158] = 0.011403206043039185964847059552;
    w[159] = 0.011347607895545491941625714297;
    w[160] = 0.011290307495875509508367594121;
    w[161] = 0.011231313439649668572656802083;
    w[162] = 0.011170634576553449462710881938;
    w[163] = 0.011108280009009843630460815451;
    w[164] = 0.011044259090813901263517571044;
    w[165] = 0.010978581425729570637988203448;
    w[166] = 0.010911256866049039700796847788;
    w[167] = 0.010842295511114795995293477058;
    w[168] = 0.010771707705804626636653631927;
    w[169] = 0.010699504038979785603048200583;
    w[170] = 0.010625695341896561133961681801;
    w[171] = 0.010550292686581481517533575536;
    w[172] = 0.010473307384170403003569566927;
    w[173] = 0.010394750983211728997101725205;
    w[174] = 0.010314635267934015068260713997;
    w[175] = 0.010232972256478219656954857160;
    w[176] = 0.010149774199094865654634066042;
    w[177] = 0.010065053576306383309460978930;
    w[178] = 0.009978823097034910124733949495;
    w[179] = 0.009891095696695828602630683809;
    w[180] = 0.009801884535257327825498800250;
    w[181] = 0.009711202995266279964249670496;
    w[182] = 0.009619064679840727857162164401;
    w[183] = 0.009525483410629284811829685754;
    w[184] = 0.009430473225737752747352764482;
    w[185] = 0.009334048377623269712466014486;
    w[186] = 0.009236223330956302687378716714;
    w[187] = 0.009137012760450806402000472219;
    w[188] = 0.009036431548662873680227775572;
    w[189] = 0.008934494783758207548408417085;
    w[190] = 0.008831217757248750025318272685;
    w[191] = 0.008726615961698807140336632217;
    w[192] = 0.008620705088401014305368838410;
    w[193] = 0.008513501025022490693838354790;
    w[194] = 0.008405019853221535756180301698;
    w[195] = 0.008295277846235225425171412553;
    w[196] = 0.008184291466438269935619761004;
    w[197] = 0.008072077362873499500946974804;
    w[198] = 0.007958652368754348353613161227;
    w[199] = 0.007844033498939711866810316151;
    w[200] = 0.007728237947381555631110194958;
    w[201] = 0.007611283084545659461618719618;
    w[202] = 0.007493186454805883358599761133;
    w[203] = 0.007373965773812346437572440695;
    w[204] = 0.007253638925833913783829137214;
    w[205] = 0.007132223961075390071672422986;
    w[206] = 0.007009739092969822621234436194;
    w[207] = 0.006886202695446320346713323775;
    w[208] = 0.006761633300173798780927861108;
    w[209] = 0.006636049593781065044590038355;
    w[210] = 0.006509470415053660267809899951;
    w[211] = 0.006381914752107880570375164275;
    w[212] = 0.006253401739542401272063645975;
    w[213] = 0.006123950655567932542389081187;
    w[214] = 0.005993580919115338221127696870;
    w[215] = 0.005862312086922653060661598801;
    w[216] = 0.005730163850601437177384417555;
    w[217] = 0.0055971560336829100775514452572;
    w[218] = 0.0054633085886443102775705318566;
    w[219] = 0.0053286415939159303170811114788;
    w[220] = 0.0051931752508692809303287536296;
    w[221] = 0.0050569298807868423875578160762;
    w[222] = 0.0049199259218138656695587765655;
    w[223] = 0.0047821839258926913729317340448;
    w[224] = 0.0046437245556800603139790923525;
    w[225] = 0.0045045685814478970686417923159;
    w[226] = 0.0043647368779680566815684200621;
    w[227] = 0.0042242504213815362723565049060;
    w[228] = 0.0040831302860526684085997759212;
    w[229] = 0.0039413976414088336277290349840;
    w[230] = 0.0037990737487662579981170192082;
    w[231] = 0.0036561799581425021693892413052;
    w[232] = 0.0035127377050563073309710549844;
    w[233] = 0.0033687685073155510120191062489;
    w[234] = 0.0032242939617941981570107134269;
    w[235] = 0.0030793357411993375832053528316;
    w[236] = 0.0029339155908297166460123254142;
    w[237] = 0.0027880553253277068805747610763;
    w[238] = 0.0026417768254274905641208292516;
    w[239] = 0.0024951020347037068508395354372;
    w[240] = 0.0023480529563273120170064609087;
    w[241] = 0.0022006516498399104996848834189;
    w[242] = 0.0020529202279661431745487818492;
    w[243] = 0.0019048808534997184044191411746;
    w[244] = 0.0017565557363307299936069145295;
    w[245] = 0.0016079671307493272424499395690;
    w[246] = 0.0014591373333107332010883864996;
    w[247] = 0.0013100886819025044578316804271;
    w[248] = 0.0011608435575677247239705981135;
    w[249] = 0.0010114243932084404526058128414;
    w[250] = 0.0008618537014200890378140934163;
    w[251] = 0.0007121541634733206669089891511;
    w[252] = 0.00056234895403140980281523674759;
    w[253] = 0.00041246325442617632843218583774;
    w[254] = 0.00026253494429644590628745756250;
    w[255] = 0.00011278901782227217551253887725;
  }
  else if ( n == 257 )
  {
    x[0] = -0.999956390712330402472857;
    x[1] = -0.999770232390338019056053;
    x[2] = -0.999435348366365078441838;
    x[3] = -0.998951714093223210129834;
    x[4] = -0.998319392445383847808766;
    x[5] = -0.997538475365520218731818;
    x[6] = -0.996609078365487004512326;
    x[7] = -0.995531339486830143483750;
    x[8] = -0.994305419008553630362377;
    x[9] = -0.992931499332908653172844;
    x[10] = -0.991409784923101705201254;
    x[11] = -0.989740502257507526030375;
    x[12] = -0.987923899788618253106809;
    x[13] = -0.985960247902290665366669;
    x[14] = -0.983849838875444644048531;
    x[15] = -0.981592986831381877693095;
    x[16] = -0.979190027692327124191591;
    x[17] = -0.976641319128992592610888;
    x[18] = -0.973947240507062326750976;
    x[19] = -0.971108192830542793021113;
    x[20] = -0.968124598681952354372943;
    x[21] = -0.964996902159337170373447;
    x[22] = -0.961725568810109767190665;
    x[23] = -0.958311085561711847074814;
    x[24] = -0.954753960649106318830855;
    x[25] = -0.951054723539105826691801;
    x[26] = -0.947213924851546682950881;
    x[27] = -0.943232136277318328151464;
    x[28] = -0.939109950493259404355123;
    x[29] = -0.934847981073932324370129;
    x[30] = -0.930446862400288909805510;
    x[31] = -0.925907249565240289235888;
    x[32] = -0.921229818276144817520964;
    x[33] = -0.916415264754228313295468;
    x[34] = -0.911464305630951423630955;
    x[35] = -0.906377677841339419411308;
    x[36] = -0.901156138514290206476301;
    x[37] = -0.895800464859876809085345;
    x[38] = -0.890311454053661045810287;
    x[39] = -0.884689923118035575018750;
    x[40] = -0.878936708800611938658765;
    x[41] = -0.873052667449672679799858;
    x[42] = -0.867038674886706051812473;
    x[43] = -0.860895626276042275514686;
    x[44] = -0.854624435991610735314055;
    x[45] = -0.848226037480837936478636;
    x[46] = -0.841701383125706473284556;
    x[47] = -0.835051444100995681967937;
    x[48] = -0.828277210229725073186687;
    x[49] = -0.821379689835822056081139;
    x[50] = -0.814359909594035880004229;
    x[51] = -0.807218914377120130552073;
    x[52] = -0.799957767100306523636066;
    x[53] = -0.792577548563093144962574;
    x[54] = -0.785079357288370682385816;
    x[55] = -0.777464309358910595129671;
    x[56] = -0.769733538251239556788216;
    x[57] = -0.761888194666924898264210;
    x[58] = -0.753929446361296162339238;
    x[59] = -0.745858477969628263337895;
    x[60] = -0.737676490830812123299244;
    x[61] = -0.729384702808539030149808;
    x[62] = -0.720984348110025333531072;
    x[63] = -0.712476677102304460118510;
    x[64] = -0.703862956126113592426171;
    x[65] = -0.695144467307402713168813;
    x[66] = -0.686322508366494071200553;
    x[67] = -0.677398392424920474813593;
    x[68] = -0.668373447809971163711735;
    x[69] = -0.659249017856974352220492;
    x[70] = -0.650026460709345873208532;
    x[71] = -0.640707149116433684724434;
    x[72] = -0.631292470229188329449219;
    x[73] = -0.621783825393689760680446;
    x[74] = -0.612182629942561267650033;
    x[75] = -0.602490312984301547488097;
    x[76] = -0.592708317190566281032495;
    x[77] = -0.582838098581430874902446;
    x[78] = -0.572881126308666332759406;
    x[79] = -0.562838882437060514424546;
    x[80] = -0.552712861723817332466074;
    x[81] = -0.542504571396066721967792;
    x[82] = -0.532215530926518500400434;
    x[83] = -0.521847271807293510797499;
    x[84] = -0.511401337321965712746629;
    x[85] = -0.500879282315849152005553;
    x[86] = -0.490282672964564000798817;
    x[87] = -0.479613086540916117008992;
    x[88] = -0.468872111180124821505728;
    x[89] = -0.458061345643433838720630;
    x[90] = -0.447182399080140586238810;
    x[91] = -0.436236890788079234603398;
    x[92] = -0.425226449972593188682213;
    x[93] = -0.414152715504032866791986;
    x[94] = -0.403017335673814873281489;
    x[95] = -0.391821967949078874408131;
    x[96] = -0.380568278725978696070941;
    x[97] = -0.369257943081644365255611;
    x[98] = -0.357892644524852014873858;
    x[99] = -0.346474074745438764010632;
    x[100] = -0.335003933362499872399782;
    x[101] = -0.323483927671405649204085;
    x[102] = -0.311915772389675771851948;
    x[103] = -0.300301189401748840754520;
    x[104] = -0.288641907502685160168097;
    x[105] = -0.276939662140840894253032;
    x[106] = -0.265196195159551900488370;
    x[107] = -0.253413254537865690008131;
    x[108] = -0.241592594130360106108882;
    x[109] = -0.229735973406087448117604;
    x[110] = -0.217845157186682897983880;
    x[111] = -0.205921915383676231351599;
    x[112] = -0.193968022735045913454182;
    x[113] = -0.181985258541054792946197;
    x[114] = -0.169975406399406713716337;
    x[115] = -0.157940253939763465806087;
    x[116] = -0.145881592557661591770148;
    x[117] = -0.133801217147868654144405;
    x[118] = -0.121700925837218653121859;
    x[119] = -0.109582519716966361063898;
    x[120] = -0.097447802574700412082119;
    x[121] = -0.085298580625855050603929;
    x[122] = -0.073136662244860502573600;
    x[123] = -0.060963857695971986730406;
    x[124] = -0.048781978863817431238958;
    x[125] = -0.036592838983704002816750;
    x[126] = -0.024398252371723591403953;
    x[127] = -0.012200034154697423345412;
    x[128] = 0.000000000000000000000000;
    x[129] = 0.012200034154697423345412;
    x[130] = 0.024398252371723591403953;
    x[131] = 0.036592838983704002816750;
    x[132] = 0.048781978863817431238958;
    x[133] = 0.060963857695971986730406;
    x[134] = 0.073136662244860502573600;
    x[135] = 0.085298580625855050603929;
    x[136] = 0.097447802574700412082119;
    x[137] = 0.109582519716966361063898;
    x[138] = 0.121700925837218653121859;
    x[139] = 0.133801217147868654144405;
    x[140] = 0.145881592557661591770148;
    x[141] = 0.157940253939763465806087;
    x[142] = 0.169975406399406713716337;
    x[143] = 0.181985258541054792946197;
    x[144] = 0.193968022735045913454182;
    x[145] = 0.205921915383676231351599;
    x[146] = 0.217845157186682897983880;
    x[147] = 0.229735973406087448117604;
    x[148] = 0.241592594130360106108882;
    x[149] = 0.253413254537865690008131;
    x[150] = 0.265196195159551900488370;
    x[151] = 0.276939662140840894253032;
    x[152] = 0.288641907502685160168097;
    x[153] = 0.300301189401748840754520;
    x[154] = 0.311915772389675771851948;
    x[155] = 0.323483927671405649204085;
    x[156] = 0.335003933362499872399782;
    x[157] = 0.346474074745438764010632;
    x[158] = 0.357892644524852014873858;
    x[159] = 0.369257943081644365255611;
    x[160] = 0.380568278725978696070941;
    x[161] = 0.391821967949078874408131;
    x[162] = 0.403017335673814873281489;
    x[163] = 0.414152715504032866791986;
    x[164] = 0.425226449972593188682213;
    x[165] = 0.436236890788079234603398;
    x[166] = 0.447182399080140586238810;
    x[167] = 0.458061345643433838720630;
    x[168] = 0.468872111180124821505728;
    x[169] = 0.479613086540916117008992;
    x[170] = 0.490282672964564000798817;
    x[171] = 0.500879282315849152005553;
    x[172] = 0.511401337321965712746629;
    x[173] = 0.521847271807293510797499;
    x[174] = 0.532215530926518500400434;
    x[175] = 0.542504571396066721967792;
    x[176] = 0.552712861723817332466074;
    x[177] = 0.562838882437060514424546;
    x[178] = 0.572881126308666332759406;
    x[179] = 0.582838098581430874902446;
    x[180] = 0.592708317190566281032495;
    x[181] = 0.602490312984301547488097;
    x[182] = 0.612182629942561267650033;
    x[183] = 0.621783825393689760680446;
    x[184] = 0.631292470229188329449219;
    x[185] = 0.640707149116433684724434;
    x[186] = 0.650026460709345873208532;
    x[187] = 0.659249017856974352220492;
    x[188] = 0.668373447809971163711735;
    x[189] = 0.677398392424920474813593;
    x[190] = 0.686322508366494071200553;
    x[191] = 0.695144467307402713168813;
    x[192] = 0.703862956126113592426171;
    x[193] = 0.712476677102304460118510;
    x[194] = 0.720984348110025333531072;
    x[195] = 0.729384702808539030149808;
    x[196] = 0.737676490830812123299244;
    x[197] = 0.745858477969628263337895;
    x[198] = 0.753929446361296162339238;
    x[199] = 0.761888194666924898264210;
    x[200] = 0.769733538251239556788216;
    x[201] = 0.777464309358910595129671;
    x[202] = 0.785079357288370682385816;
    x[203] = 0.792577548563093144962574;
    x[204] = 0.799957767100306523636066;
    x[205] = 0.807218914377120130552073;
    x[206] = 0.814359909594035880004229;
    x[207] = 0.821379689835822056081139;
    x[208] = 0.828277210229725073186687;
    x[209] = 0.835051444100995681967937;
    x[210] = 0.841701383125706473284556;
    x[211] = 0.848226037480837936478636;
    x[212] = 0.854624435991610735314055;
    x[213] = 0.860895626276042275514686;
    x[214] = 0.867038674886706051812473;
    x[215] = 0.873052667449672679799858;
    x[216] = 0.878936708800611938658765;
    x[217] = 0.884689923118035575018750;
    x[218] = 0.890311454053661045810287;
    x[219] = 0.895800464859876809085345;
    x[220] = 0.901156138514290206476301;
    x[221] = 0.906377677841339419411308;
    x[222] = 0.911464305630951423630955;
    x[223] = 0.916415264754228313295468;
    x[224] = 0.921229818276144817520964;
    x[225] = 0.925907249565240289235888;
    x[226] = 0.930446862400288909805510;
    x[227] = 0.934847981073932324370129;
    x[228] = 0.939109950493259404355123;
    x[229] = 0.943232136277318328151464;
    x[230] = 0.947213924851546682950881;
    x[231] = 0.951054723539105826691801;
    x[232] = 0.954753960649106318830855;
    x[233] = 0.958311085561711847074814;
    x[234] = 0.961725568810109767190665;
    x[235] = 0.964996902159337170373447;
    x[236] = 0.968124598681952354372943;
    x[237] = 0.971108192830542793021113;
    x[238] = 0.973947240507062326750976;
    x[239] = 0.976641319128992592610888;
    x[240] = 0.979190027692327124191591;
    x[241] = 0.981592986831381877693095;
    x[242] = 0.983849838875444644048531;
    x[243] = 0.985960247902290665366669;
    x[244] = 0.987923899788618253106809;
    x[245] = 0.989740502257507526030375;
    x[246] = 0.991409784923101705201254;
    x[247] = 0.992931499332908653172844;
    x[248] = 0.994305419008553630362377;
    x[249] = 0.995531339486830143483750;
    x[250] = 0.996609078365487004512326;
    x[251] = 0.997538475365520218731818;
    x[252] = 0.998319392445383847808766;
    x[253] = 0.998951714093223210129834;
    x[254] = 0.999435348366365078441838;
    x[255] = 0.999770232390338019056053;
    x[256] = 0.999956390712330402472857;

    w[0] = 0.00011191470145601756450862287886;
    w[1] = 0.00026049995580176964436806680831;
    w[2] = 0.00040926648283531339591138751432;
    w[3] = 0.00055799120546880640169677292533;
    w[4] = 0.00070663671051592291949335494247;
    w[5] = 0.00085517818446696565626595950963;
    w[6] = 0.00100359280467969441299468763292;
    w[7] = 0.0011518582377826677880963146741;
    w[8] = 0.0012999523174235227389668643832;
    w[9] = 0.0014478529559255120065233994722;
    w[10] = 0.0015955381166175133369701690235;
    w[11] = 0.0017429858051468299509941139300;
    w[12] = 0.0018901740676190104269878470891;
    w[13] = 0.0020370809914723626741694800322;
    w[14] = 0.0021836847075455253317921866057;
    w[15] = 0.0023299633927021828561308282641;
    w[16] = 0.0024758952727301488651840215879;
    w[17] = 0.0026214586253808109266552781372;
    w[18] = 0.0027666317834818283552560256501;
    w[19] = 0.0029113931380877846359302447381;
    w[20] = 0.0030557211416493711130936102459;
    w[21] = 0.0031995943111899437356540290142;
    w[22] = 0.0033429912314827618499065991316;
    w[23] = 0.0034858905582247143702551557840;
    w[24] = 0.0036282710212037760873102463983;
    w[25] = 0.0037701114274582873548537007645;
    w[26] = 0.0039113906644266662571543468015;
    w[27] = 0.0040520877030864825223229951262;
    w[28] = 0.0041921816010820254766367595011;
    w[29] = 0.0043316515058396297504806208252;
    w[30] = 0.0044704766576701092218388764046;
    w[31] = 0.0046086363928577081326523656522;
    w[32] = 0.0047461101467350184936945641585;
    w[33] = 0.0048828774567433411142588306018;
    w[34] = 0.0050189179654779878773297516544;
    w[35] = 0.0051542114237180378340642003713;
    w[36] = 0.0052887376934400710240953933529;
    w[37] = 0.0054224767508154127788846727083;
    w[38] = 0.0055554086891904284012033890901;
    w[39] = 0.0056875137220494140577838938236;
    w[40] = 0.0058187721859596348346566361185;
    w[41] = 0.0059491645434980654366600347567;
    w[42] = 0.0060786713861593931405204596709;
    w[43] = 0.0062072734372448464599330978665;
    w[44] = 0.0063349515547314166407936938524;
    w[45] = 0.0064616867341210426397202932350;
    w[46] = 0.0065874601112693336961737372300;
    w[47] = 0.0067122529651934070221351960200;
    w[48] = 0.0068360467208584215286561508406;
    w[49] = 0.0069588229519423919043121805236;
    w[50] = 0.0070805633835788707705149901066;
    w[51] = 0.0072012498950770900730828552207;
    w[52] = 0.0073208645226191563361371026044;
    w[53] = 0.0074393894619338979090297315972;
    w[54] = 0.0075568070709469658838993300454;
    w[55] = 0.0076730998724067939537782250476;
    w[56] = 0.0077882505564860261212726654404;
    w[57] = 0.0079022419833580248574070864277;
    w[58] = 0.0080150571857480760504667455353;
    w[59] = 0.0081266793714589108764118189068;
    w[60] = 0.0082370919258701685661946145361;
    w[61] = 0.0083462784144114279413811886655;
    w[62] = 0.0084542225850084395379670551258;
    w[63] = 0.0085609083705021941391459209280;
    w[64] = 0.0086663198910404675908861979240;
    w[65] = 0.0087704414564414858792445834744;
    w[66] = 0.0088732575685293586050755892934;
    w[67] = 0.0089747529234409331997949023068;
    w[68] = 0.0090749124139037264846862498962;
    w[69] = 0.0091737211314845944854270065178;
    w[70] = 0.0092711643688088057725325917169;
    w[71] = 0.0093672276217491880067391857021;
    w[72] = 0.0094618965915850218253881576301;
    w[73] = 0.0095551571871303607110514249099;
    w[74] = 0.0096469955268314600363329731559;
    w[75] = 0.0097373979408330030783691793250;
    w[76] = 0.0098263509730128164423854701706;
    w[77] = 0.0099138413829847720250916955489;
    w[78] = 0.0099998561480695773850435626986;
    w[79] = 0.0100843824652331611676814627839;
    w[80] = 0.0101674077529923650568895461852;
    w[81] = 0.0102489196532876585918958554047;
    w[82] = 0.0103289060333225980974485876288;
    w[83] = 0.0104073549873697559257355517893;
    w[84] = 0.0104842548385428511997370260353;
    w[85] = 0.0105595941405348182788823332058;
    w[86] = 0.0106333616793215542382761147904;
    w[87] = 0.0107055464748310917616231511294;
    w[88] = 0.0107761377825779489945556541150;
    w[89] = 0.0108451250952624130885928632830;
    w[90] = 0.0109124981443345193856719616965;
    w[91] = 0.0109782469015224934483083029166;
    w[92] = 0.0110423615803254284301924654946;
    w[93] = 0.0111048326374699756056269264803;
    w[94] = 0.0111656507743308312328559850485;
    w[95] = 0.0112248069383148083152535688671;
    w[96] = 0.0112822923242082872447042603128;
    w[97] = 0.0113380983754878447625379269120;
    w[98] = 0.011392216785593866154247619654;
    w[99] = 0.011444639499166951104119199270;
    w[100] = 0.011495358713246929174010288914;
    w[101] = 0.011544366878434306436012137033;
    w[102] = 0.011591656700013970380783131035;
    w[103] = 0.011637221139040985841125311445;
    w[104] = 0.011681053413388320313049670635;
    w[105] = 0.011723146998756342723302879656;
    w[106] = 0.011763495629643945382264331878;
    w[107] = 0.011802093300281144573421477037;
    w[108] = 0.011838934265523020964443424791;
    w[109] = 0.011874013041704866779344562066;
    w[110] = 0.011907324407458412445505183140;
    w[111] = 0.011938863404489011222535627643;
    w[112] = 0.011968625338313666131272065445;
    w[113] = 0.011996605778959789329711050159;
    w[114] = 0.012022800561624589927558893338;
    w[115] = 0.012047205787294992091420946532;
    w[116] = 0.012069817823327991167612855626;
    w[117] = 0.012090633303991361438266420912;
    w[118] = 0.012109649130964635027950450318;
    w[119] = 0.012126862473800277391553601370;
    w[120] = 0.012142270770344990738801546574;
    w[121] = 0.012155871727121082685623083829;
    w[122] = 0.012167663319667843366755737416;
    w[123] = 0.012177643792842880196606249581;
    w[124] = 0.012185811661083365425569178819;
    w[125] = 0.012192165708627157605870499188;
    w[126] = 0.012196704989693764053654538465;
    w[127] = 0.012199428828625117371582840212;
    w[128] = 0.01220033681998614507777289232;
    w[129] = 0.012199428828625117371582840212;
    w[130] = 0.012196704989693764053654538465;
    w[131] = 0.012192165708627157605870499188;
    w[132] = 0.012185811661083365425569178819;
    w[133] = 0.012177643792842880196606249581;
    w[134] = 0.012167663319667843366755737416;
    w[135] = 0.012155871727121082685623083829;
    w[136] = 0.012142270770344990738801546574;
    w[137] = 0.012126862473800277391553601370;
    w[138] = 0.012109649130964635027950450318;
    w[139] = 0.012090633303991361438266420912;
    w[140] = 0.012069817823327991167612855626;
    w[141] = 0.012047205787294992091420946532;
    w[142] = 0.012022800561624589927558893338;
    w[143] = 0.011996605778959789329711050159;
    w[144] = 0.011968625338313666131272065445;
    w[145] = 0.011938863404489011222535627643;
    w[146] = 0.011907324407458412445505183140;
    w[147] = 0.011874013041704866779344562066;
    w[148] = 0.011838934265523020964443424791;
    w[149] = 0.011802093300281144573421477037;
    w[150] = 0.011763495629643945382264331878;
    w[151] = 0.011723146998756342723302879656;
    w[152] = 0.011681053413388320313049670635;
    w[153] = 0.011637221139040985841125311445;
    w[154] = 0.011591656700013970380783131035;
    w[155] = 0.011544366878434306436012137033;
    w[156] = 0.011495358713246929174010288914;
    w[157] = 0.011444639499166951104119199270;
    w[158] = 0.011392216785593866154247619654;
    w[159] = 0.0113380983754878447625379269120;
    w[160] = 0.0112822923242082872447042603128;
    w[161] = 0.0112248069383148083152535688671;
    w[162] = 0.0111656507743308312328559850485;
    w[163] = 0.0111048326374699756056269264803;
    w[164] = 0.0110423615803254284301924654946;
    w[165] = 0.0109782469015224934483083029166;
    w[166] = 0.0109124981443345193856719616965;
    w[167] = 0.0108451250952624130885928632830;
    w[168] = 0.0107761377825779489945556541150;
    w[169] = 0.0107055464748310917616231511294;
    w[170] = 0.0106333616793215542382761147904;
    w[171] = 0.0105595941405348182788823332058;
    w[172] = 0.0104842548385428511997370260353;
    w[173] = 0.0104073549873697559257355517893;
    w[174] = 0.0103289060333225980974485876288;
    w[175] = 0.0102489196532876585918958554047;
    w[176] = 0.0101674077529923650568895461852;
    w[177] = 0.0100843824652331611676814627839;
    w[178] = 0.0099998561480695773850435626986;
    w[179] = 0.0099138413829847720250916955489;
    w[180] = 0.0098263509730128164423854701706;
    w[181] = 0.0097373979408330030783691793250;
    w[182] = 0.0096469955268314600363329731559;
    w[183] = 0.0095551571871303607110514249099;
    w[184] = 0.0094618965915850218253881576301;
    w[185] = 0.0093672276217491880067391857021;
    w[186] = 0.0092711643688088057725325917169;
    w[187] = 0.0091737211314845944854270065178;
    w[188] = 0.0090749124139037264846862498962;
    w[189] = 0.0089747529234409331997949023068;
    w[190] = 0.0088732575685293586050755892934;
    w[191] = 0.0087704414564414858792445834744;
    w[192] = 0.0086663198910404675908861979240;
    w[193] = 0.0085609083705021941391459209280;
    w[194] = 0.0084542225850084395379670551258;
    w[195] = 0.0083462784144114279413811886655;
    w[196] = 0.0082370919258701685661946145361;
    w[197] = 0.0081266793714589108764118189068;
    w[198] = 0.0080150571857480760504667455353;
    w[199] = 0.0079022419833580248574070864277;
    w[200] = 0.0077882505564860261212726654404;
    w[201] = 0.0076730998724067939537782250476;
    w[202] = 0.0075568070709469658838993300454;
    w[203] = 0.0074393894619338979090297315972;
    w[204] = 0.0073208645226191563361371026044;
    w[205] = 0.0072012498950770900730828552207;
    w[206] = 0.0070805633835788707705149901066;
    w[207] = 0.0069588229519423919043121805236;
    w[208] = 0.0068360467208584215286561508406;
    w[209] = 0.0067122529651934070221351960200;
    w[210] = 0.0065874601112693336961737372300;
    w[211] = 0.0064616867341210426397202932350;
    w[212] = 0.0063349515547314166407936938524;
    w[213] = 0.0062072734372448464599330978665;
    w[214] = 0.0060786713861593931405204596709;
    w[215] = 0.0059491645434980654366600347567;
    w[216] = 0.0058187721859596348346566361185;
    w[217] = 0.0056875137220494140577838938236;
    w[218] = 0.0055554086891904284012033890901;
    w[219] = 0.0054224767508154127788846727083;
    w[220] = 0.0052887376934400710240953933529;
    w[221] = 0.0051542114237180378340642003713;
    w[222] = 0.0050189179654779878773297516544;
    w[223] = 0.0048828774567433411142588306018;
    w[224] = 0.0047461101467350184936945641585;
    w[225] = 0.0046086363928577081326523656522;
    w[226] = 0.0044704766576701092218388764046;
    w[227] = 0.0043316515058396297504806208252;
    w[228] = 0.0041921816010820254766367595011;
    w[229] = 0.0040520877030864825223229951262;
    w[230] = 0.0039113906644266662571543468015;
    w[231] = 0.0037701114274582873548537007645;
    w[232] = 0.0036282710212037760873102463983;
    w[233] = 0.0034858905582247143702551557840;
    w[234] = 0.0033429912314827618499065991316;
    w[235] = 0.0031995943111899437356540290142;
    w[236] = 0.0030557211416493711130936102459;
    w[237] = 0.0029113931380877846359302447381;
    w[238] = 0.0027666317834818283552560256501;
    w[239] = 0.0026214586253808109266552781372;
    w[240] = 0.0024758952727301488651840215879;
    w[241] = 0.0023299633927021828561308282641;
    w[242] = 0.0021836847075455253317921866057;
    w[243] = 0.0020370809914723626741694800322;
    w[244] = 0.0018901740676190104269878470891;
    w[245] = 0.0017429858051468299509941139300;
    w[246] = 0.0015955381166175133369701690235;
    w[247] = 0.0014478529559255120065233994722;
    w[248] = 0.0012999523174235227389668643832;
    w[249] = 0.0011518582377826677880963146741;
    w[250] = 0.00100359280467969441299468763292;
    w[251] = 0.00085517818446696565626595950963;
    w[252] = 0.00070663671051592291949335494247;
    w[253] = 0.00055799120546880640169677292533;
    w[254] = 0.00040926648283531339591138751432;
    w[255] = 0.00026049995580176964436806680831;
    w[256] = 0.00011191470145601756450862287886;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1:33, 63/64/65, 127/128/129, 255/256/257\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_cos ( int n, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_COS: Gauss-Legendre rules for COS(X)*F(X) on [-PI/2,PI/2].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -PI/2 <= X <= PI/2 ) COS(X) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) WEIGHT[I] * F ( XTAB[I] )
//
//  Discussion:
//
//    The same rule can be used to approximate
//
//      Integral ( 0 <= X <= PI ) SIN(X) * F(X) dX
//
//    as
//
//      Sum ( 1 <= I <= N ) WEIGHT[I] * F ( XTAB[I] + PI/2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993, QA299.3E93, page 310.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1, 2, 4, 8 or 16.
//
//    Output, double XTAB[N], the abscissas.
//
//    Output, double WEIGHT[N], the weights.
//
{
  if ( n == 1 )
  {
    xtab[0] = 0.0E+00;

    weight[0] = 2.0E+00;
  }
  else if ( n == 2 )
  {
    xtab[0] = - 0.68366739008990304094E+00;
    xtab[1] =   0.68366739008990304094E+00;

    weight[0] = 1.0E+00;
    weight[1] = 1.0E+00;
  }
  else if ( n == 4 )
  {
    xtab[0] = - 1.1906765638948557415E+00;
    xtab[1] = - 0.43928746686001514756E+00;
    xtab[2] =   0.43928746686001514756E+00;
    xtab[3] =   1.1906765638948557415E+00;

    weight[0] = 0.22407061812762016065E+00;
    weight[1] = 0.77592938187237983935E+00;
    weight[2] = 0.77592938187237983935E+00;
    weight[3] = 0.22407061812762016065E+00;
  }
  else if ( n == 8 )
  {
    xtab[0] = - 1.4414905401823575701E+00;
    xtab[1] = - 1.1537256454567275850E+00;
    xtab[2] = - 0.74346864787549244989E+00;
    xtab[3] = - 0.25649650741623123020E+00;
    xtab[4] =   0.25649650741623123020E+00;
    xtab[5] =   0.74346864787549244989E+00;
    xtab[6] =   1.1537256454567275850E+00;
    xtab[7] =   1.4414905401823575701E+00;

    weight[0] = 0.027535633513767011149E+00;
    weight[1] = 0.14420409203022750950E+00;
    weight[2] = 0.33626447785280459621E+00;
    weight[3] = 0.49199579660320088314E+00;
    weight[4] = 0.49199579660320088314E+00;
    weight[5] = 0.33626447785280459621E+00;
    weight[6] = 0.14420409203022750950E+00;
    weight[7] = 0.027535633513767011149E+00;
  }
  else if ( n == 16 )
  {
    xtab[ 0] = - 1.5327507132362304779E+00;
    xtab[ 1] = - 1.4446014873666514608E+00;
    xtab[ 2] = - 1.3097818904452936698E+00;
    xtab[ 3] = - 1.1330068786005003695E+00;
    xtab[ 4] = - 0.92027786206637096497E+00;
    xtab[ 5] = - 0.67861108097560545347E+00;
    xtab[ 6] = - 0.41577197673418943962E+00;
    xtab[ 7] = - 0.14003444424696773778E+00;
    xtab[ 8] =   0.14003444424696773778E+00;
    xtab[9] =   0.41577197673418943962E+00;
    xtab[10] =   0.67861108097560545347E+00;
    xtab[11] =   0.92027786206637096497E+00;
    xtab[12] =   1.1330068786005003695E+00;
    xtab[13] =   1.3097818904452936698E+00;
    xtab[14] =   1.4446014873666514608E+00;
    xtab[15] =   1.5327507132362304779E+00;

    weight[ 0] = 0.0024194677567615628193E+00;
    weight[ 1] = 0.014115268156854008264E+00;
    weight[ 2] = 0.040437893946503669410E+00;
    weight[ 3] = 0.083026647573217742131E+00;
    weight[ 4] = 0.13834195526951273359E+00;
    weight[ 5] = 0.19741148870253455567E+00;
    weight[ 6] = 0.24763632094635522403E+00;
    weight[ 7] = 0.27661095764826050408E+00;
    weight[ 8] = 0.27661095764826050408E+00;
    weight[9] = 0.24763632094635522403E+00;
    weight[10] = 0.19741148870253455567E+00;
    weight[11] = 0.13834195526951273359E+00;
    weight[12] = 0.083026647573217742131E+00;
    weight[13] = 0.040437893946503669410E+00;
    weight[14] = 0.014115268156854008264E+00;
    weight[15] = 0.0024194677567615628193E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_COS - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_cos2 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_COS2: Gauss-Legendre rules for COS(X)*F(X) on [0,PI/2].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= PI/2 ) COS(X) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
//
//    The integral:
//
//      Integral ( 0 <= X <= PI/2 ) SIN(X) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( PI/2 - XTAb[I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993, QA299.3E93, page 311.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 2, 4, 8 or 16.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 2 )
  {
    xtab[0] = 0.26587388056307823382E+00;
    xtab[1] = 1.0351526093171315182E+00;

    weight[0] = 0.60362553280827113087E+00;
    weight[1] = 0.39637446719172886913E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] = 0.095669389196858636773E+00;
    xtab[1] = 0.45240902327067096554E+00;
    xtab[2] = 0.93185057672024082424E+00;
    xtab[3] = 1.3564439599666466230E+00;

    weight[ 0] = 0.23783071419515504517E+00;
    weight[ 1] = 0.40265695523581253512E+00;
    weight[ 2] = 0.28681737948564715225E+00;
    weight[ 3] = 0.072694951083385267446E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] = 0.029023729768913933432E+00;
    xtab[1] = 0.14828524404581819442E+00;
    xtab[2] = 0.34531111151664787488E+00;
    xtab[3] = 0.59447696797658360178E+00;
    xtab[4] = 0.86538380686123504827E+00;
    xtab[5] = 1.1263076093187456632E+00;
    xtab[6] = 1.3470150460281258016E+00;
    xtab[7] = 1.5015603622059195568E+00;

    weight[ 0] = 0.073908998095117384985E+00;
    weight[ 1] = 0.16002993702338006099E+00;
    weight[ 2] = 0.21444434341803549108E+00;
    weight[ 3] = 0.21979581268851903339E+00;
    weight[ 4] = 0.17581164478209568886E+00;
    weight[ 5] = 0.10560448025308322171E+00;
    weight[ 6] = 0.042485497299217201089E+00;
    weight[ 7] = 0.0079192864405519178899E+00;
  }
  else if ( order == 16 )
  {
    xtab[ 0] = 0.0080145034906295973494E+00;
    xtab[ 1] = 0.041893031354246254797E+00;
    xtab[ 2] = 0.10149954486757579459E+00;
    xtab[ 3] = 0.18463185923836617507E+00;
    xtab[ 4] = 0.28826388487760574589E+00;
    xtab[ 5] = 0.40870579076464794191E+00;
    xtab[ 6] = 0.54176054986913847463E+00;
    xtab[ 7] = 0.68287636658719416893E+00;
    xtab[ 8] = 0.82729287620416833520E+00;
    xtab[9] = 0.97018212594829367065E+00;
    xtab[10] = 1.1067865150286247873E+00;
    xtab[11] = 1.2325555697227748824E+00;
    xtab[12] = 1.3432821921580721861E+00;
    xtab[13] = 1.4352370549295032923E+00;
    xtab[14] = 1.5052970876794669248E+00;
    xtab[15] = 1.5510586944086135769E+00;

    weight[ 0] = 0.020528714977215248902E+00;
    weight[ 1] = 0.046990919853597958123E+00;
    weight[ 2] = 0.071441021312218541698E+00;
    weight[ 3] = 0.092350338329243052271E+00;
    weight[ 4] = 0.10804928026816236935E+00;
    weight[ 5] = 0.11698241243306261791E+00;
    weight[ 6] = 0.11812395361762037649E+00;
    weight[ 7] = 0.11137584940420091049E+00;
    weight[ 8] = 0.097778236145946543110E+00;
    weight[9] = 0.079418758985944482077E+00;
    weight[10] = 0.059039620053768691402E+00;
    weight[11] = 0.039458876783728165671E+00;
    weight[12] = 0.022987785677206847531E+00;
    weight[13] = 0.011010405600421536861E+00;
    weight[14] = 0.0038123928030499915653E+00;
    weight[15] = 0.00065143375461266656171E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_COS2 - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_log ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_LOG sets a Gauss-Legendre rule for - LOG(X) * F(X) on [0,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) - LOG(X) * F(X) dX
// 
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Gwynne Evans,
//    Practical Numerical Integration,
//    Wiley, 1993, QA299.3E93, page 309.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 through 8, or 16.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0] = 0.25E+00;

    weight[0] = 1.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[0] = 0.112008806166976182957205488948E+00;
    xtab[1] = 0.602276908118738102757080225338E+00;

    weight[0] = 0.718539319030384440665510200891E+00;
    weight[1] = 0.281460680969615559334489799109E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] = 0.0638907930873254049961166031363E+00;
    xtab[1] = 0.368997063715618765546197645857E+00;
    xtab[2] = 0.766880303938941455423682659817E+00;

    weight[0] = 0.513404552232363325129300497567E+00;
    weight[1] = 0.391980041201487554806287180966E+00;
    weight[2] = 0.0946154065661491200644123214672E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] = 0.0414484801993832208033213101564E+00;
    xtab[1] = 0.245274914320602251939675759523E+00;
    xtab[2] = 0.556165453560275837180184354376E+00;
    xtab[3] = 0.848982394532985174647849188085E+00;

    weight[0] = 0.383464068145135124850046522343E+00;
    weight[1] = 0.386875317774762627336008234554E+00;
    weight[2] = 0.190435126950142415361360014547E+00;
    weight[3] = 0.0392254871299598324525852285552E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] = 0.0291344721519720533037267621154E+00;
    xtab[1] = 0.173977213320897628701139710829E+00;
    xtab[2] = 0.411702520284902043174931924646E+00;
    xtab[3] = 0.677314174582820380701802667998E+00;
    xtab[4] = 0.894771361031008283638886204455E+00;

    weight[0] = 0.297893471782894457272257877888E+00;
    weight[1] = 0.349776226513224180375071870307E+00;
    weight[2] = 0.234488290044052418886906857943E+00;
    weight[3] = 0.0989304595166331469761807114404E+00;
    weight[4] = 0.0189115521431957964895826824218E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] = 0.0216340058441169489956958558537E+00;
    xtab[1] = 0.129583391154950796131158505009E+00;
    xtab[2] = 0.314020449914765508798248188420E+00;
    xtab[3] = 0.538657217351802144548941893993E+00;
    xtab[4] = 0.756915337377402852164544156139E+00;
    xtab[5] = 0.922668851372120237333873231507E+00;

    weight[0] = 0.238763662578547569722268303330E+00;
    weight[1] = 0.308286573273946792969383109211E+00;
    weight[2] = 0.245317426563210385984932540188E+00;
    weight[3] = 0.142008756566476685421345576030E+00;
    weight[4] = 0.0554546223248862900151353549662E+00;
    weight[5] = 0.0101689586929322758869351162755E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] = 0.0167193554082585159416673609320E+00;
    xtab[1] = 0.100185677915675121586885031757E+00;
    xtab[2] = 0.246294246207930599046668547239E+00;
    xtab[3] = 0.433463493257033105832882482601E+00;
    xtab[4] = 0.632350988047766088461805812245E+00;
    xtab[5] = 0.811118626740105576526226796782E+00;
    xtab[6] = 0.940848166743347721760134113379E+00;

    weight[0] = 0.196169389425248207525427377585E+00;
    weight[1] = 0.270302644247272982145271719533E+00;
    weight[2] = 0.239681873007690948308072785252E+00;
    weight[3] = 0.165775774810432906560869687736E+00;
    weight[4] = 0.0889432271376579644357238403458E+00;
    weight[5] = 0.0331943043565710670254494111034E+00;
    weight[6] = 0.00593278701512592399918517844468E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] = 0.0133202441608924650122526725243E+00;
    xtab[1] = 0.0797504290138949384098277291424E+00;
    xtab[2] = 0.197871029326188053794476159516E+00;
    xtab[3] = 0.354153994351909419671463603538E+00;
    xtab[4] = 0.529458575234917277706149699996E+00;
    xtab[5] = 0.701814529939099963837152670310E+00;
    xtab[6] = 0.849379320441106676048309202301E+00;
    xtab[7] = 0.953326450056359788767379678514E+00;

    weight[0] = 0.164416604728002886831472568326E+00;
    weight[1] = 0.237525610023306020501348561960E+00;
    weight[2] = 0.226841984431919126368780402936E+00;
    weight[3] = 0.175754079006070244988056212006E+00;
    weight[4] = 0.112924030246759051855000442086E+00;
    weight[5] = 0.0578722107177820723985279672940E+00;
    weight[6] = 0.0209790737421329780434615241150E+00;
    weight[7] = 0.00368640710402761901335232127647E+00;
  }
  else if ( order == 16 )
  {
    xtab[ 0] = 0.00389783448711591592405360527037E+00;
    xtab[ 1] = 0.0230289456168732398203176309848E+00;
    xtab[ 2] = 0.0582803983062404123483532298394E+00;
    xtab[ 3] = 0.108678365091054036487713613051E+00;
    xtab[ 4] = 0.172609454909843937760843776232E+00;
    xtab[ 5] = 0.247937054470578495147671753047E+00;
    xtab[ 6] = 0.332094549129917155984755859320E+00;
    xtab[ 7] = 0.422183910581948600115088366745E+00;
    xtab[ 8] = 0.515082473381462603476277704052E+00;
    xtab[9] = 0.607556120447728724086384921709E+00;
    xtab[10] = 0.696375653228214061156318166581E+00;
    xtab[11] = 0.778432565873265405203868167732E+00;
    xtab[12] = 0.850850269715391083233822761319E+00;
    xtab[13] = 0.911086857222271905418818994060E+00;
    xtab[14] = 0.957025571703542157591520509383E+00;
    xtab[15] = 0.987047800247984476758697436516E+00;

    weight[ 0] = 0.0607917100435912328511733871235E+00;
    weight[ 1] = 0.102915677517582144387691736210E+00;
    weight[ 2] = 0.122355662046009193557547513197E+00;
    weight[ 3] = 0.127569246937015988717042209239E+00;
    weight[ 4] = 0.123013574600070915423123365137E+00;
    weight[ 5] = 0.111847244855485722621848903429E+00;
    weight[ 6] = 0.0965963851521243412529681650802E+00;
    weight[ 7] = 0.0793566643514731387824416770520E+00;
    weight[ 8] = 0.0618504945819652070951360853113E+00;
    weight[9] = 0.0454352465077266686288299526509E+00;
    weight[10] = 0.0310989747515818064092528627927E+00;
    weight[11] = 0.0194597659273608420780860268669E+00;
    weight[12] = 0.0107762549632055256455393162159E+00;
    weight[13] = 0.00497254289008764171250524951646E+00;
    weight[14] = 0.00167820111005119451503546419059E+00;
    weight[15] = 0.000282353764668436321778085987413E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_LOG - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_sqrtx_01 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_SQRTX_01 sets Gauss-Legendre rules for SQRT(X)*F(X) on [0,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) SQRT ( X ) * F(X) dX =
//      Integral ( 0 <= Y <= 1 ) 2 * Y*Y * F(Y*Y) dY.
//      (using Y = SQRT(X) )
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  int i;
  int order2;
  double *xtab2;
  double *weight2;

  order2 = 2 * order + 1;

  xtab2 = new double[order2];
  weight2 = new double[order2];

  legendre_set ( order2, xtab2, weight2 );

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = pow ( xtab2[order+1+i], 2 );
  }
  for ( i = 0; i < order; i++ )
  {
    weight[i] = 2.0 * weight2[order+1+i] * xtab[i];
  }

  delete [] xtab2;
  delete [] weight2;

  return;
}
//****************************************************************************80

void legendre_set_sqrtx2_01 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_SQRTX2_01: Gauss-Legendre rules for F(X)/SQRT(X) on [0,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) F(X) / SQRT ( X ) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  int i;
  int order2;
  double *xtab2;
  double *weight2;

  order2 = 2 * order + 1;

  xtab2 = new double[order2];
  weight2 = new double[order2];

  legendre_set ( order2, xtab2, weight2 );

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = pow ( xtab2[order+1+i], 2 );
  }
  for ( i = 0; i < order; i++ )
  {
    weight[i] = 2.0 * weight2[order+1+i];
  }

  delete [] xtab2;
  delete [] weight2;

  return;
}
//****************************************************************************80

void legendre_set_x0_01 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_X0_01 sets a Gauss-Legendre rule for F(X) on [0,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 8.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0] =   0.5E+00;

    weight[0] = 1.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[0] = 0.2113248654E+00;
    xtab[1] = 0.7886751346E+00;

    weight[0] = 0.5E+00;
    weight[1] = 0.5E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] = 0.1127016654E+00;
    xtab[1] = 0.5000000000E+00;
    xtab[2] = 0.8872983346E+00;

    weight[0] = 5.0E+00 / 18.0E+00;
    weight[1] = 8.0E+00 / 18.0E+00;
    weight[2] = 5.0E+00 / 18.0E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] = 0.0694318442E+00;
    xtab[1] = 0.3300094782E+00;
    xtab[2] = 0.6699905218E+00;
    xtab[3] = 0.9305681558E+00;

    weight[0] = 0.1739274226E+00;
    weight[1] = 0.3260725774E+00;
    weight[2] = 0.3260725774E+00;
    weight[3] = 0.1739274226E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] = 0.0469100770E+00;
    xtab[1] = 0.2307653449E+00;
    xtab[2] = 0.5000000000E+00;
    xtab[3] = 0.7692346551E+00;
    xtab[4] = 0.9530899230E+00;

    weight[0] = 0.1184634425E+00;
    weight[1] = 0.2393143352E+00;
    weight[2] = 0.2844444444E+00;
    weight[3] = 0.2393143352E+00;
    weight[4] = 0.1184634425E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] = 0.0337652429E+00;
    xtab[1] = 0.1693953068E+00;
    xtab[2] = 0.3806904070E+00;
    xtab[3] = 0.6193095930E+00;
    xtab[4] = 0.8306046932E+00;
    xtab[5] = 0.9662347571E+00;

    weight[0] = 0.0856622462E+00;
    weight[1] = 0.1803807865E+00;
    weight[2] = 0.2339569673E+00;
    weight[3] = 0.2339569673E+00;
    weight[4] = 0.1803807865E+00;
    weight[5] = 0.0856622462E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] = 0.0254460438E+00;
    xtab[1] = 0.1292344072E+00;
    xtab[2] = 0.2970774243E+00;
    xtab[3] = 0.5000000000E+00;
    xtab[4] = 0.7029225757E+00;
    xtab[5] = 0.8707655928E+00;
    xtab[6] = 0.9745539562E+00;

    weight[0] = 0.0647424831E+00;
    weight[1] = 0.1398526957E+00;
    weight[2] = 0.1909150253E+00;
    weight[3] = 0.2089795918E+00;
    weight[4] = 0.1909150253E+00;
    weight[5] = 0.1398526957E+00;
    weight[6] = 0.0647424831E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] = 0.0198550718E+00;
    xtab[1] = 0.1016667613E+00;
    xtab[2] = 0.2372337950E+00;
    xtab[3] = 0.4082826788E+00;
    xtab[4] = 0.5917173212E+00;
    xtab[5] = 0.7627662050E+00;
    xtab[6] = 0.8983332387E+00;
    xtab[7] = 0.9801449282E+00;

    weight[0] = 0.0506142681E+00;
    weight[1] = 0.1111905172E+00;
    weight[2] = 0.1568533229E+00;
    weight[3] = 0.1813418917E+00;
    weight[4] = 0.1813418917E+00;
    weight[5] = 0.1568533229E+00;
    weight[6] = 0.1111905172E+00;
    weight[7] = 0.0506142681E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_X0_01 - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_x1 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_X1 sets a Gauss-Legendre rule for ( 1 + X ) * F(X) on [-1,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) ( 1 + X ) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 9.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0] =  0.333333333333333333333333333333E+00;

    weight[0] = 2.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[0] = -0.289897948556635619639456814941E+00;
    xtab[1] =  0.689897948556635619639456814941E+00;

    weight[0] =  0.727834473024091322422523991699E+00;
    weight[1] =  1.27216552697590867757747600830E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] = -0.575318923521694112050483779752E+00;
    xtab[1] =  0.181066271118530578270147495862E+00;
    xtab[2] =  0.822824080974592105208907712461E+00;

    weight[0] =  0.279307919605816490135525088716E+00;
    weight[1] =  0.916964425438344986775682378225E+00;
    weight[2] =  0.803727654955838523088792533058E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] = -0.720480271312438895695825837750E+00;
    xtab[1] = -0.167180864737833640113395337326E+00;
    xtab[2] =  0.446313972723752344639908004629E+00;
    xtab[3] =  0.885791607770964635613757614892E+00;

    weight[0] =  0.124723883800032328695500588386E+00;
    weight[1] =  0.519390190432929763305824811559E+00;
    weight[2] =  0.813858272041085443165617903743E+00;
    weight[3] =  0.542027653725952464833056696312E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] = -0.802929828402347147753002204224E+00;
    xtab[1] = -0.390928546707272189029229647442E+00;
    xtab[2] =  0.124050379505227711989974959990E+00;
    xtab[3] =  0.603973164252783654928415726409E+00;
    xtab[4] =  0.920380285897062515318386619813E+00;

    weight[0] =  0.0629916580867691047411692662740E+00;
    weight[1] =  0.295635480290466681402532877367E+00;
    weight[2] =  0.585547948338679234792151477424E+00;
    weight[3] =  0.668698552377478261966702492391E+00;
    weight[4] =  0.387126360906606717097443886545E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] = -0.853891342639482229703747931639E+00;
    xtab[1] = -0.538467724060109001833766720231E+00;
    xtab[2] = -0.117343037543100264162786683611E+00;
    xtab[3] =  0.326030619437691401805894055838E+00;
    xtab[4] =  0.703842800663031416300046295008E+00;
    xtab[5] =  0.941367145680430216055899446174E+00;

    weight[0] =  0.0349532072544381270240692132496E+00;
    weight[1] =  0.175820662202035902032706497222E+00;
    weight[2] =  0.394644603562621056482338042193E+00;
    weight[3] =  0.563170215152795712476307356284E+00;
    weight[4] =  0.542169988926074467362761586552E+00;
    weight[5] =  0.289241322902034734621817304499E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] = -0.887474878926155707068695617935E+00;
    xtab[1] = -0.639518616526215270024840114382E+00;
    xtab[2] = -0.294750565773660725252184459658E+00;
    xtab[3] =  0.0943072526611107660028971153047E+00;
    xtab[4] =  0.468420354430821063046421216613E+00;
    xtab[5] =  0.770641893678191536180719525865E+00;
    xtab[6] =  0.955041227122575003782349000858E+00;

    weight[0] =  0.0208574488112296163587654972151E+00;
    weight[1] =  0.109633426887493901777324193433E+00;
    weight[2] =  0.265538785861965879934591955055E+00;
    weight[3] =  0.428500262783494679963649011999E+00;
    weight[4] =  0.509563589198353307674937943100E+00;
    weight[5] =  0.442037032763498409684482945478E+00;
    weight[6] =  0.223869453693964204606248453720E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] = -0.910732089420060298533757956283E+00;
    xtab[1] = -0.711267485915708857029562959544E+00;
    xtab[2] = -0.426350485711138962102627520502E+00;
    xtab[3] = -0.0903733696068532980645444599064E+00;
    xtab[4] =  0.256135670833455395138292079035E+00;
    xtab[5] =  0.571383041208738483284917464837E+00;
    xtab[6] =  0.817352784200412087992517083851E+00;
    xtab[7] =  0.964440169705273096373589797925E+00;

    weight[0] =  0.0131807657689951954189692640444E+00;
    weight[1] =  0.0713716106239448335742111888042E+00;
    weight[2] =  0.181757278018795592332221684383E+00;
    weight[3] =  0.316798397969276640481632757440E+00;
    weight[4] =  0.424189437743720042818124385645E+00;
    weight[5] =  0.450023197883549464687088394417E+00;
    weight[6] =  0.364476094545494505382889847132E+00;
    weight[7] =  0.178203217446223725304862478136E+00;
  }
  else if ( order == 9 )
  {
    xtab[0] = -0.927484374233581078117671398464E+00;
    xtab[1] = -0.763842042420002599615429776011E+00;
    xtab[2] = -0.525646030370079229365386614293E+00;
    xtab[3] = -0.236234469390588049278459503207E+00;
    xtab[4] =  0.0760591978379781302337137826389E+00;
    xtab[5] =  0.380664840144724365880759065541E+00;
    xtab[6] =  0.647766687674009436273648507855E+00;
    xtab[7] =  0.851225220581607910728163628088E+00;
    xtab[8] =  0.971175180702246902734346518378E+00;

    weight[0] =  0.00872338834309252349019620448007E+00;
    weight[1] =  0.0482400171391415162069086091476E+00;
    weight[2] =  0.127219285964216005046760427743E+00;
    weight[3] =  0.233604781180660442262926091607E+00;
    weight[4] =  0.337433287379681397577000079834E+00;
    weight[5] =  0.401235236773473158616600898930E+00;
    weight[6] =  0.394134968689382820640692081477E+00;
    weight[7] =  0.304297020437232650320317215016E+00;
    weight[8] =  0.145112014093119485838598391765E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_X1 - Fatal error!\n";
    cerr << "  Illegal input value of ORDER = " << order << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_x1_01 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_X1_01 sets a Gauss-Legendre rule for X * F(X) on [0,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) X * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 8.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0] =   0.6666666667E+00;

    weight[0] = 0.5000000000E+00;
  }
  else if ( order == 2 )
  {
    xtab[0] = 0.3550510257E+00;
    xtab[1] = 0.8449489743E+00;

    weight[0] = 0.1819586183E+00;
    weight[1] = 0.3180413817E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] = 0.2123405382E+00;
    xtab[1] = 0.5905331356E+00;
    xtab[2] = 0.9114120405E+00;

    weight[0] = 0.0698269799E+00;
    weight[1] = 0.2292411064E+00;
    weight[2] = 0.2009319137E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] = 0.1397598643E+00;
    xtab[1] = 0.4164095676E+00;
    xtab[2] = 0.7231569864E+00;
    xtab[3] = 0.9428958039E+00;

    weight[0] = 0.0311809710E+00;
    weight[1] = 0.1298475476E+00;
    weight[2] = 0.2034645680E+00;
    weight[3] = 0.1355069134E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] = 0.0985350858E+00;
    xtab[1] = 0.3045357266E+00;
    xtab[2] = 0.5620251898E+00;
    xtab[3] = 0.8019865821E+00;
    xtab[4] = 0.9601901429E+00;

    weight[0] = 0.0157479145E+00;
    weight[1] = 0.0739088701E+00;
    weight[2] = 0.1463869871E+00;
    weight[3] = 0.1671746381E+00;
    weight[4] = 0.0967815902E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] = 0.0730543287E+00;
    xtab[1] = 0.2307661380E+00;
    xtab[2] = 0.4413284812E+00;
    xtab[3] = 0.6630153097E+00;
    xtab[4] = 0.8519214003E+00;
    xtab[5] = 0.9706835728E+00;

    weight[0] = 0.0087383108E+00;
    weight[1] = 0.0439551656E+00;
    weight[2] = 0.0986611509E+00;
    weight[3] = 0.1407925538E+00;
    weight[4] = 0.1355424972E+00;
    weight[5] = 0.0723103307E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] = 0.0562625605E+00;
    xtab[1] = 0.1802406917E+00;
    xtab[2] = 0.3526247171E+00;
    xtab[3] = 0.5471536263E+00;
    xtab[4] = 0.7342101772E+00;
    xtab[5] = 0.8853209468E+00;
    xtab[6] = 0.9775206136E+00;

    weight[0] = 0.0052143622E+00;
    weight[1] = 0.0274083567E+00;
    weight[2] = 0.0663846965E+00;
    weight[3] = 0.1071250657E+00;
    weight[4] = 0.1273908973E+00;
    weight[5] = 0.1105092582E+00;
    weight[6] = 0.0559673634E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] = 0.0446339553E+00;
    xtab[1] = 0.1443662570E+00;
    xtab[2] = 0.2868247571E+00;
    xtab[3] = 0.4548133152E+00;
    xtab[4] = 0.6280678354E+00;
    xtab[5] = 0.7856915206E+00;
    xtab[6] = 0.9086763921E+00;
    xtab[7] = 0.9822200849E+00;

    weight[0] = 0.0032951914E+00;
    weight[1] = 0.0178429027E+00;
    weight[2] = 0.0454393195E+00;
    weight[3] = 0.0791995995E+00;
    weight[4] = 0.1060473594E+00;
    weight[5] = 0.1125057995E+00;
    weight[6] = 0.0911190236E+00;
    weight[7] = 0.0445508044E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_X1_01 - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void legendre_set_x2 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_X2 sets Gauss-Legendre rules for ( 1 + X )^2*F(X) on [-1,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) ( 1 + X ) * ( 1 + X ) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 9.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0] =  0.5E+00;

    weight[0] =  2.66666666666666666666666666666E+00;
  }
  else if ( order == 2 )
  {
    xtab[0] = -0.0883036880224505775998524725910E+00;
    xtab[1] =  0.754970354689117244266519139258E+00;

    weight[0] =  0.806287056638603444666851075928E+00;
    weight[1] =  1.86037961002806322199981559074E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] = -0.410004419776996766244796955168E+00;
    xtab[1] =  0.305992467923296230556472913192E+00;
    xtab[2] =  0.854011951853700535688324041976E+00;

    weight[0] =  0.239605624068645584091811926047E+00;
    weight[1] =  1.16997015407892817602809616291E+00;
    weight[2] =  1.25709088851909290654675857771E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] = -0.591702835793545726606755921586E+00;
    xtab[1] = -0.0340945902087350046811467387661E+00;
    xtab[2] =  0.522798524896275389882037174551E+00;
    xtab[3] =  0.902998901106005341405865485802E+00;

    weight[0] =  0.0828179259993445222751812523731E+00;
    weight[1] =  0.549071097383384602539010760334E+00;
    weight[2] =  1.14767031839371367238662411421E+00;
    weight[3] =  0.887107324890223869465850539752E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] = -0.702108425894032836232448374820E+00;
    xtab[1] = -0.268666945261773544694327777841E+00;
    xtab[2] =  0.220227225868961343518209179230E+00;
    xtab[3] =  0.653039358456608553790815164028E+00;
    xtab[4] =  0.930842120163569816951085142737E+00;

    weight[0] =  0.0329106016247920636689299329544E+00;
    weight[1] =  0.256444805783695354037991444453E+00;
    weight[2] =  0.713601289772720001490035944563E+00;
    weight[3] =  1.00959169519929190423066348132E+00;
    weight[4] =  0.654118274286167343239045863379E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] = -0.773611232355123732602532012021E+00;
    xtab[1] = -0.431362254623427837535325249187E+00;
    xtab[2] = -0.0180728263295041680220798103354E+00;
    xtab[3] =  0.395126163954217534500188844163E+00;
    xtab[4] =  0.736872116684029732026178298518E+00;
    xtab[5] =  0.948190889812665614490712786006E+00;

    weight[0] =  0.0146486064549543818622276447204E+00;
    weight[1] =  0.125762377479560410622810097040E+00;
    weight[2] =  0.410316569036929681761034600615E+00;
    weight[3] =  0.756617493988329628546336413760E+00;
    weight[4] =  0.859011997894245060846045458784E+00;
    weight[5] =  0.500309621812647503028212451747E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] = -0.822366333126005527278634734418E+00;
    xtab[1] = -0.547034493182875002223997992852E+00;
    xtab[2] = -0.200043026557985860387937545780E+00;
    xtab[3] =  0.171995710805880507163425502299E+00;
    xtab[4] =  0.518891747903884926692601716998E+00;
    xtab[5] =  0.793821941703901970495546427988E+00;
    xtab[6] =  0.959734452453198985538996625765E+00;

    weight[0] =  0.00714150426951365443207221475404E+00;
    weight[1] =  0.0653034050584375560578544725498E+00;
    weight[2] =  0.235377690316228918725962815880E+00;
    weight[3] =  0.505171029671130381676271523850E+00;
    weight[4] =  0.733870426238362032891332767175E+00;
    weight[5] =  0.725590596901489156295739839779E+00;
    weight[6] =  0.394212014211504966587433032679E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] = -0.857017929919813794402037235698E+00;
    xtab[1] = -0.631543407166567521509503573952E+00;
    xtab[2] = -0.339104543648722903660229021109E+00;
    xtab[3] = -0.0111941563689783438801237300122E+00;
    xtab[4] =  0.316696017045595559454075475675E+00;
    xtab[5] =  0.609049663022520165351466780939E+00;
    xtab[6] =  0.834198765028697794599267293239E+00;
    xtab[7] =  0.967804480896157932935972899807E+00;

    weight[0] =  0.00374814227227757804631954025851E+00;
    weight[1] =  0.0357961737041152639660521680263E+00;
    weight[2] =  0.137974910241879862433949246199E+00;
    weight[3] =  0.326515411108352185491692769217E+00;
    weight[4] =  0.547577467373226177976217604887E+00;
    weight[5] =  0.682278153375510121675529810121E+00;
    weight[6] =  0.614544746137780998436053880546E+00;
    weight[7] =  0.318231662453524478640851647411E+00;
  }
  else if ( order == 9 )
  {
    xtab[0] = -0.882491728426548422828684254270E+00;
    xtab[1] = -0.694873684026474640346360850039E+00;
    xtab[2] = -0.446537143480670863635920316400E+00;
    xtab[3] = -0.159388112702326252531544826624E+00;
    xtab[4] =  0.141092709224374414981503995427E+00;
    xtab[5] =  0.428217823321559204544020866175E+00;
    xtab[6] =  0.676480966471850715860378175342E+00;
    xtab[7] =  0.863830940812464825046988286026E+00;
    xtab[8] =  0.973668228805771018909618924364E+00;

    weight[0] =  0.00209009877215570354392734918986E+00;
    weight[1] =  0.0205951891648697848186537272448E+00;
    weight[2] =  0.0832489326348178964194106978875E+00;
    weight[3] =  0.210746247220398685903797568021E+00;
    weight[4] =  0.388325022916052063676224499399E+00;
    weight[5] =  0.554275165518437673725822282791E+00;
    weight[6] =  0.621388553284444032628761363828E+00;
    weight[7] =  0.523916296267173054255512857631E+00;
    weight[8] =  0.262081160888317771694556320674E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_X2 - Fatal error!\n";
    cerr << "  Illegal input value of ORDER = " << order << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_x2_01 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_X2_01 sets a Gauss-Legendre rule for X*X * F(X) on [0,1].
//
//  Discussion:
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) X * X * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 8.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0] =   0.75E+00;

    weight[0] = 1.0E+00 / 3.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[0] = 0.4558481560E+00;
    xtab[1] = 0.8774851773E+00;

    weight[0] = 0.1007858821E+00;
    weight[1] = 0.2325474513E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] = 0.2949977901E+00;
    xtab[1] = 0.6529962340E+00;
    xtab[2] = 0.9270059759E+00;

    weight[0] = 0.0299507030E+00;
    weight[1] = 0.1462462693E+00;
    weight[2] = 0.1571363611E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] = 0.2041485821E+00;
    xtab[1] = 0.4829527049E+00;
    xtab[2] = 0.7613992624E+00;
    xtab[3] = 0.9514994506E+00;

    weight[0] = 0.0103522408E+00;
    weight[1] = 0.0686338872E+00;
    weight[2] = 0.1434587898E+00;
    weight[3] = 0.1108884156E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] = 0.1489457871E+00;
    xtab[1] = 0.3656665274E+00;
    xtab[2] = 0.6101136129E+00;
    xtab[3] = 0.8265196792E+00;
    xtab[4] = 0.9654210601E+00;

    weight[0] = 0.0041138252E+00;
    weight[1] = 0.0320556007E+00;
    weight[2] = 0.0892001612E+00;
    weight[3] = 0.1261989619E+00;
    weight[4] = 0.0817647843E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] = 0.1131943838E+00;
    xtab[1] = 0.2843188727E+00;
    xtab[2] = 0.4909635868E+00;
    xtab[3] = 0.6975630820E+00;
    xtab[4] = 0.8684360583E+00;
    xtab[5] = 0.9740954449E+00;

    weight[0] = 0.0018310758E+00;
    weight[1] = 0.0157202972E+00;
    weight[2] = 0.0512895711E+00;
    weight[3] = 0.0945771867E+00;
    weight[4] = 0.1073764997E+00;
    weight[5] = 0.0625387027E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] = 0.0888168334E+00;
    xtab[1] = 0.2264827534E+00;
    xtab[2] = 0.3999784867E+00;
    xtab[3] = 0.5859978554E+00;
    xtab[4] = 0.7594458740E+00;
    xtab[5] = 0.8969109709E+00;
    xtab[6] = 0.9798672262E+00;

    weight[0] = 0.0008926880E+00;
    weight[1] = 0.0081629256E+00;
    weight[2] = 0.0294222113E+00;
    weight[3] = 0.0631463787E+00;
    weight[4] = 0.0917338033E+00;
    weight[5] = 0.0906988246E+00;
    weight[6] = 0.0492765018E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] = 0.0714910350E+00;
    xtab[1] = 0.1842282964E+00;
    xtab[2] = 0.3304477282E+00;
    xtab[3] = 0.4944029218E+00;
    xtab[4] = 0.6583480085E+00;
    xtab[5] = 0.8045248315E+00;
    xtab[6] = 0.9170993825E+00;
    xtab[7] = 0.9839022404E+00;

    weight[0] = 0.0004685178E+00;
    weight[1] = 0.0044745217E+00;
    weight[2] = 0.0172468638E+00;
    weight[3] = 0.0408144264E+00;
    weight[4] = 0.0684471834E+00;
    weight[5] = 0.0852847692E+00;
    weight[6] = 0.0768180933E+00;
    weight[7] = 0.0397789578E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_X2_01 - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void lobatto_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_COMPUTE computes a Lobatto quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) )
//
//    The quadrature rule will integrate exactly all polynomials up to
//    X**(2*N-3).
//
//    The Lobatto rule is distinguished by the fact that both endpoints
//    (-1 and 1) are always abscissas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2007
//
//  Author:
//
//    Original MATLAB version by Greg von Winckel.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
//    Spectral Methods in Fluid Dynamics,
//    Springer, 1993,
//    ISNB13: 978-3540522058,
//    LC: QA377.S676.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.  N must be at least 2.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  int j;
  double *p;
  double pi = 3.141592653589793;
  double test;
  double tolerance;
  double *xold;

  if ( n < 2 )
  {
    cerr << "\n";
    cerr << "LOBATTO_COMPUTE - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << n << "\n";
    cerr << "  ORDER must be at least 2.\n";
    exit ( 1 );
  }
  tolerance = 100.0 * r8_epsilon ( );
//
//  Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto nodes.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( pi * ( double ) ( i ) / ( double ) ( n - 1 ) );
  }

  xold = new double[n];
  p = new double[n*n];

  do
  {
    for ( i = 0; i < n; i++ )
    {
      xold[i] = x[i];
    }
    for ( i = 0; i < n; i++ )
    {
      p[i+0*n] = 1.0;
    }
    for ( i = 0; i < n; i++ )
    {
      p[i+1*n] = x[i];
    }

    for ( j = 2; j <= n-1; j++ )
    {
      for ( i = 0; i < n; i++)
      {
        p[i+j*n] = ( ( double ) ( 2 * j - 1 ) * x[i] * p[i+(j-1)*n]     
                   + ( double ) (   - j + 1 ) *        p[i+(j-2)*n] ) 
                   / ( double ) (     j     );
      }
    }

    for ( i = 0; i < n; i++ )
    {
      x[i] = xold[i] - ( x[i] * p[i+(n-1)*n] - p[i+(n-2)*n] ) 
           / ( ( double ) ( n ) * p[i+(n-1)*n] );
    }

    test = 0.0;
    for ( i = 0; i < n; i++ )
    {
      test = r8_max ( test, r8_abs ( x[i] - xold[i] ) );
    }

  } while ( tolerance < test );

  r8vec_reverse ( n, x );

  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( ( double ) ( ( n - 1 ) * n ) * pow ( p[i+(n-1)*n], 2 ) );
  }

  delete [] p;
  delete [] xold;

  return;
}
//****************************************************************************80

void lobatto_set ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_SET sets abscissas and weights for Lobatto quadrature.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
//
//    The quadrature rule will integrate exactly all polynomials up to
//    X**(2*ORDER-3).
//
//    The Lobatto rule is distinguished by the fact that both endpoints
//    (-1 and 1) are always abscissas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 2 and 20.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 2 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =    1.0E+00;

    weight[0] =  1.0E+00;
    weight[1] =  1.0E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =    0.0E+00;
    xtab[2] =    1.0E+00;

    weight[0] =  1.0 / 3.0E+00;
    weight[1] =  4.0 / 3.0E+00;
    weight[2] =  1.0 / 3.0E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.447213595499957939281834733746E+00;
    xtab[2] =    0.447213595499957939281834733746E+00;
    xtab[3] =    1.0E+00;

    weight[0] =  1.0E+00 / 6.0E+00;
    weight[1] =  5.0E+00 / 6.0E+00;
    weight[2] =  5.0E+00 / 6.0E+00;
    weight[3] =  1.0E+00 / 6.0E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.654653670707977143798292456247E+00;
    xtab[2] =    0.0E+00;
    xtab[3] =    0.654653670707977143798292456247E+00;
    xtab[4] =    1.0E+00;

    weight[0] =  9.0E+00 / 90.0E+00;
    weight[1] = 49.0E+00 / 90.0E+00;
    weight[2] = 64.0E+00 / 90.0E+00;
    weight[3] = 49.0E+00 / 90.0E+00;
    weight[4] =  9.0E+00 / 90.0E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.765055323929464692851002973959E+00;
    xtab[2] =  - 0.285231516480645096314150994041E+00;
    xtab[3] =    0.285231516480645096314150994041E+00;
    xtab[4] =    0.765055323929464692851002973959E+00;
    xtab[5] =    1.0E+00;

    weight[0] =  0.066666666666666666666666666667E+00;
    weight[1] =  0.378474956297846980316612808212E+00;
    weight[2] =  0.554858377035486353016720525121E+00;
    weight[3] =  0.554858377035486353016720525121E+00;
    weight[4] =  0.378474956297846980316612808212E+00;
    weight[5] =  0.066666666666666666666666666667E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.830223896278566929872032213967E+00;
    xtab[2] =  - 0.468848793470714213803771881909E+00;
    xtab[3] =    0.0E+00;
    xtab[4] =    0.468848793470714213803771881909E+00;
    xtab[5] =    0.830223896278566929872032213967E+00;
    xtab[6] =    1.0E+00;

    weight[0] =  0.476190476190476190476190476190E-01;
    weight[1] =  0.276826047361565948010700406290E+00;
    weight[2] =  0.431745381209862623417871022281E+00;
    weight[3] =  0.487619047619047619047619047619E+00;
    weight[4] =  0.431745381209862623417871022281E+00;
    weight[5] =  0.276826047361565948010700406290E+00;
    weight[6] =  0.476190476190476190476190476190E-01;
  }
  else if ( order == 8 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.871740148509606615337445761221E+00;
    xtab[2] =  - 0.591700181433142302144510731398E+00;
    xtab[3] =  - 0.209299217902478868768657260345E+00;
    xtab[4] =    0.209299217902478868768657260345E+00;
    xtab[5] =    0.591700181433142302144510731398E+00;
    xtab[6] =    0.871740148509606615337445761221E+00;
    xtab[7] =    1.0E+00;

    weight[0] =  0.357142857142857142857142857143E-01;
    weight[1] =  0.210704227143506039382991065776E+00;
    weight[2] =  0.341122692483504364764240677108E+00;
    weight[3] =  0.412458794658703881567052971402E+00;
    weight[4] =  0.412458794658703881567052971402E+00;
    weight[5] =  0.341122692483504364764240677108E+00;
    weight[6] =  0.210704227143506039382991065776E+00;
    weight[7] =  0.357142857142857142857142857143E-01;
  }
  else if ( order == 9 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.899757995411460157312345244418E+00;
    xtab[2] =  - 0.677186279510737753445885427091E+00;
    xtab[3] =  - 0.363117463826178158710752068709E+00;
    xtab[4] =    0.0E+00;
    xtab[5] =    0.363117463826178158710752068709E+00;
    xtab[6] =    0.677186279510737753445885427091E+00;
    xtab[7] =    0.899757995411460157312345244418E+00;
    xtab[8] =    1.0E+00;

    weight[0] =  0.277777777777777777777777777778E-01;
    weight[1] =  0.165495361560805525046339720029E+00;
    weight[2] =  0.274538712500161735280705618579E+00;
    weight[3] =  0.346428510973046345115131532140E+00;
    weight[4] =  0.371519274376417233560090702948E+00;
    weight[5] =  0.346428510973046345115131532140E+00;
    weight[6] =  0.274538712500161735280705618579E+00;
    weight[7] =  0.165495361560805525046339720029E+00;
    weight[8] =  0.277777777777777777777777777778E-01;
  }
  else if ( order == 10 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.919533908166458813828932660822E+00;
    xtab[2] =  - 0.738773865105505075003106174860E+00;
    xtab[3] =  - 0.477924949810444495661175092731E+00;
    xtab[4] =  - 0.165278957666387024626219765958E+00;
    xtab[5] =    0.165278957666387024626219765958E+00;
    xtab[6] =    0.477924949810444495661175092731E+00;
    xtab[7] =    0.738773865105505075003106174860E+00;
    xtab[8] =    0.919533908166458813828932660822E+00;
    xtab[9] =   1.0E+00;

    weight[0] =  0.222222222222222222222222222222E-01;
    weight[1] =  0.133305990851070111126227170755E+00;
    weight[2] =  0.224889342063126452119457821731E+00;
    weight[3] =  0.292042683679683757875582257374E+00;
    weight[4] =  0.327539761183897456656510527917E+00;
    weight[5] =  0.327539761183897456656510527917E+00;
    weight[6] =  0.292042683679683757875582257374E+00;
    weight[7] =  0.224889342063126452119457821731E+00;
    weight[8] =  0.133305990851070111126227170755E+00;
    weight[9] = 0.222222222222222222222222222222E-01;
  }
  else if ( order == 11 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.934001430408059134332274136099E+00;
    xtab[2] =  - 0.784483473663144418622417816108E+00;
    xtab[3] =  - 0.565235326996205006470963969478E+00;
    xtab[4] =  - 0.295758135586939391431911515559E+00;
    xtab[5] =    0.0E+00;
    xtab[6] =    0.295758135586939391431911515559E+00;
    xtab[7] =    0.565235326996205006470963969478E+00;
    xtab[8] =    0.784483473663144418622417816108E+00;
    xtab[9] =   0.934001430408059134332274136099E+00;
    xtab[10] =   1.0E+00;

    weight[0] =  0.181818181818181818181818181818E-01;
    weight[1] =  0.109612273266994864461403449580E+00;
    weight[2] =  0.187169881780305204108141521899E+00;
    weight[3] =  0.248048104264028314040084866422E+00;
    weight[4] =  0.286879124779008088679222403332E+00;
    weight[5] =  0.300217595455690693785931881170E+00;
    weight[6] =  0.286879124779008088679222403332E+00;
    weight[7] =  0.248048104264028314040084866422E+00;
    weight[8] =  0.187169881780305204108141521899E+00;
    weight[9] = 0.109612273266994864461403449580E+00;
    weight[10] = 0.181818181818181818181818181818E-01;
  }
  else if ( order == 12 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.944899272222882223407580138303E+00;
    xtab[2] =  - 0.819279321644006678348641581717E+00;
    xtab[3] =  - 0.632876153031869677662404854444E+00;
    xtab[4] =  - 0.399530940965348932264349791567E+00;
    xtab[5] =  - 0.136552932854927554864061855740E+00;
    xtab[6] =    0.136552932854927554864061855740E+00;
    xtab[7] =    0.399530940965348932264349791567E+00;
    xtab[8] =    0.632876153031869677662404854444E+00;
    xtab[9] =   0.819279321644006678348641581717E+00;
    xtab[10] =   0.944899272222882223407580138303E+00;
    xtab[11] =   1.0E+00;

    weight[0] =  0.151515151515151515151515151515E-01;
    weight[1] =  0.916845174131961306683425941341E-01;
    weight[2] =  0.157974705564370115164671062700E+00;
    weight[3] =  0.212508417761021145358302077367E+00;
    weight[4] =  0.251275603199201280293244412148E+00;
    weight[5] =  0.271405240910696177000288338500E+00;
    weight[6] =  0.271405240910696177000288338500E+00;
    weight[7] =  0.251275603199201280293244412148E+00;
    weight[8] =  0.212508417761021145358302077367E+00;
    weight[9] = 0.157974705564370115164671062700E+00;
    weight[10] = 0.916845174131961306683425941341E-01;
    weight[11] = 0.151515151515151515151515151515E-01;
  }
  else if ( order == 13 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.953309846642163911896905464755E+00;
    xtab[2] =  - 0.846347564651872316865925607099E+00;
    xtab[3] =  - 0.686188469081757426072759039566E+00;
    xtab[4] =  - 0.482909821091336201746937233637E+00;
    xtab[5] =  - 0.249286930106239992568673700374E+00;
    xtab[6] =    0.0E+00;
    xtab[7] =    0.249286930106239992568673700374E+00;
    xtab[8] =    0.482909821091336201746937233637E+00;
    xtab[9] =   0.686188469081757426072759039566E+00;
    xtab[10] =   0.846347564651872316865925607099E+00;
    xtab[11] =   0.953309846642163911896905464755E+00;
    xtab[12] =   1.0E+00;

    weight[0] =  0.128205128205128205128205128205E-01;
    weight[1] =  0.778016867468189277935889883331E-01;
    weight[2] =  0.134981926689608349119914762589E+00;
    weight[3] =  0.183646865203550092007494258747E+00;
    weight[4] =  0.220767793566110086085534008379E+00;
    weight[5] =  0.244015790306676356458578148360E+00;
    weight[6] =  0.251930849333446736044138641541E+00;
    weight[7] =  0.244015790306676356458578148360E+00;
    weight[8] =  0.220767793566110086085534008379E+00;
    weight[9] = 0.183646865203550092007494258747E+00;
    weight[10] = 0.134981926689608349119914762589E+00;
    weight[11] = 0.778016867468189277935889883331E-01;
    weight[12] = 0.128205128205128205128205128205E-01;
  }
  else if ( order == 14 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.959935045267260901355100162015E+00;
    xtab[2] =  - 0.867801053830347251000220202908E+00;
    xtab[3] =  - 0.728868599091326140584672400521E+00;
    xtab[4] =  - 0.550639402928647055316622705859E+00;
    xtab[5] =  - 0.342724013342712845043903403642E+00;
    xtab[6] =  - 0.116331868883703867658776709736E+00;
    xtab[7] =    0.116331868883703867658776709736E+00;
    xtab[8] =    0.342724013342712845043903403642E+00;
    xtab[9] =   0.550639402928647055316622705859E+00;
    xtab[10] =   0.728868599091326140584672400521E+00;
    xtab[11] =   0.867801053830347251000220202908E+00;
    xtab[12] =   0.959935045267260901355100162015E+00;
    xtab[13] =   1.0E+00;

    weight[0] =  0.109890109890109890109890109890E-01;
    weight[1] =  0.668372844976812846340706607461E-01;
    weight[2] =  0.116586655898711651540996670655E+00;
    weight[3] =  0.160021851762952142412820997988E+00;
    weight[4] =  0.194826149373416118640331778376E+00;
    weight[5] =  0.219126253009770754871162523954E+00;
    weight[6] =  0.231612794468457058889628357293E+00;
    weight[7] =  0.231612794468457058889628357293E+00;
    weight[8] =  0.219126253009770754871162523954E+00;
    weight[9] = 0.194826149373416118640331778376E+00;
    weight[10] = 0.160021851762952142412820997988E+00;
    weight[11] = 0.116586655898711651540996670655E+00;
    weight[12] = 0.668372844976812846340706607461E-01;
    weight[13] = 0.109890109890109890109890109890E-01;
  }
  else if ( order == 15 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.965245926503838572795851392070E+00;
    xtab[2] =  - 0.885082044222976298825401631482E+00;
    xtab[3] =  - 0.763519689951815200704118475976E+00;
    xtab[4] =  - 0.606253205469845711123529938637E+00;
    xtab[5] =  - 0.420638054713672480921896938739E+00;
    xtab[6] =  - 0.215353955363794238225679446273E+00;
    xtab[7] =    0.0E+00;
    xtab[8] =    0.215353955363794238225679446273E+00;
    xtab[9] =   0.420638054713672480921896938739E+00;
    xtab[10] =   0.606253205469845711123529938637E+00;
    xtab[11] =   0.763519689951815200704118475976E+00;
    xtab[12] =   0.885082044222976298825401631482E+00;
    xtab[13] =   0.965245926503838572795851392070E+00;
    xtab[14] =   1.0E+00;

    weight[0] =  0.952380952380952380952380952381E-02;
    weight[1] =  0.580298930286012490968805840253E-01;
    weight[2] =  0.101660070325718067603666170789E+00;
    weight[3] =  0.140511699802428109460446805644E+00;
    weight[4] =  0.172789647253600949052077099408E+00;
    weight[5] =  0.196987235964613356092500346507E+00;
    weight[6] =  0.211973585926820920127430076977E+00;
    weight[7] =  0.217048116348815649514950214251E+00;
    weight[8] =  0.211973585926820920127430076977E+00;
    weight[9] = 0.196987235964613356092500346507E+00;
    weight[10] = 0.172789647253600949052077099408E+00;
    weight[11] = 0.140511699802428109460446805644E+00;
    weight[12] = 0.101660070325718067603666170789E+00;
    weight[13] = 0.580298930286012490968805840253E-01;
    weight[14] = 0.952380952380952380952380952381E-02;
  }
  else if ( order == 16 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.969568046270217932952242738367E+00;
    xtab[2] =  - 0.899200533093472092994628261520E+00;
    xtab[3] =  - 0.792008291861815063931088270963E+00;
    xtab[4] =  - 0.652388702882493089467883219641E+00;
    xtab[5] =  - 0.486059421887137611781890785847E+00;
    xtab[6] =  - 0.299830468900763208098353454722E+00;
    xtab[7] =  - 0.101326273521949447843033005046E+00;
    xtab[8] =    0.101326273521949447843033005046E+00;
    xtab[9] =   0.299830468900763208098353454722E+00;
    xtab[10] =   0.486059421887137611781890785847E+00;
    xtab[11] =   0.652388702882493089467883219641E+00;
    xtab[12] =   0.792008291861815063931088270963E+00;
    xtab[13] =   0.899200533093472092994628261520E+00;
    xtab[14] =   0.969568046270217932952242738367E+00;
    xtab[15] =   1.0E+00;

    weight[0] =  0.833333333333333333333333333333E-02;
    weight[1] =  0.508503610059199054032449195655E-01;
    weight[2] =  0.893936973259308009910520801661E-01;
    weight[3] =  0.124255382132514098349536332657E+00;
    weight[4] =  0.154026980807164280815644940485E+00;
    weight[5] =  0.177491913391704125301075669528E+00;
    weight[6] =  0.193690023825203584316913598854E+00;
    weight[7] =  0.201958308178229871489199125411E+00;
    weight[8] =  0.201958308178229871489199125411E+00;
    weight[9] = 0.193690023825203584316913598854E+00;
    weight[10] = 0.177491913391704125301075669528E+00;
    weight[11] = 0.154026980807164280815644940485E+00;
    weight[12] = 0.124255382132514098349536332657E+00;
    weight[13] = 0.893936973259308009910520801661E-01;
    weight[14] = 0.508503610059199054032449195655E-01;
    weight[15] = 0.833333333333333333333333333333E-02;
  }
  else if ( order == 17 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.973132176631418314156979501874E+00;
    xtab[2] =  - 0.910879995915573595623802506398E+00;
    xtab[3] =  - 0.815696251221770307106750553238E+00;
    xtab[4] =  - 0.691028980627684705394919357372E+00;
    xtab[5] =  - 0.541385399330101539123733407504E+00;
    xtab[6] =  - 0.372174433565477041907234680735E+00;
    xtab[7] =  - 0.189511973518317388304263014753E+00;
    xtab[8] =    0.0E+00;
    xtab[9] =   0.189511973518317388304263014753E+00;
    xtab[10] =   0.372174433565477041907234680735E+00;
    xtab[11] =   0.541385399330101539123733407504E+00;
    xtab[12] =   0.691028980627684705394919357372E+00;
    xtab[13] =   0.815696251221770307106750553238E+00;
    xtab[14] =   0.910879995915573595623802506398E+00;
    xtab[15] =   0.973132176631418314156979501874E+00;
    xtab[16] =   1.0E+00;

    weight[0] =  0.735294117647058823529411764706E-02;
    weight[1] =  0.449219405432542096474009546232E-01;
    weight[2] =  0.791982705036871191902644299528E-01;
    weight[3] =  0.110592909007028161375772705220E+00;
    weight[4] =  0.137987746201926559056201574954E+00;
    weight[5] =  0.160394661997621539516328365865E+00;
    weight[6] =  0.177004253515657870436945745363E+00;
    weight[7] =  0.187216339677619235892088482861E+00;
    weight[8] =  0.190661874753469433299407247028E+00;
    weight[9] = 0.187216339677619235892088482861E+00;
    weight[10] = 0.177004253515657870436945745363E+00;
    weight[11] = 0.160394661997621539516328365865E+00;
    weight[12] = 0.137987746201926559056201574954E+00;
    weight[13] = 0.110592909007028161375772705220E+00;
    weight[14] = 0.791982705036871191902644299528E-01;
    weight[15] = 0.449219405432542096474009546232E-01;
    weight[16] = 0.735294117647058823529411764706E-02;
  }
  else if ( order == 18 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.976105557412198542864518924342E+00;
    xtab[2] =  - 0.920649185347533873837854625431E+00;
    xtab[3] =  - 0.835593535218090213713646362328E+00;
    xtab[4] =  - 0.723679329283242681306210365302E+00;
    xtab[5] =  - 0.588504834318661761173535893194E+00;
    xtab[6] =  - 0.434415036912123975342287136741E+00;
    xtab[7] =  - 0.266362652878280984167665332026E+00;
    xtab[8] =  - 0.897490934846521110226450100886E-01;
    xtab[9] =   0.897490934846521110226450100886E-01;
    xtab[10] =   0.266362652878280984167665332026E+00;
    xtab[11] =   0.434415036912123975342287136741E+00;
    xtab[12] =   0.588504834318661761173535893194E+00;
    xtab[13] =   0.723679329283242681306210365302E+00;
    xtab[14] =   0.835593535218090213713646362328E+00;
    xtab[15] =   0.920649185347533873837854625431E+00;
    xtab[16] =   0.976105557412198542864518924342E+00;
    xtab[17] =   1.0E+00;

    weight[0] =  0.653594771241830065359477124183E-02;
    weight[1] =  0.399706288109140661375991764101E-01;
    weight[2] =  0.706371668856336649992229601678E-01;
    weight[3] =  0.990162717175028023944236053187E-01;
    weight[4] =  0.124210533132967100263396358897E+00;
    weight[5] =  0.145411961573802267983003210494E+00;
    weight[6] =  0.161939517237602489264326706700E+00;
    weight[7] =  0.173262109489456226010614403827E+00;
    weight[8] =  0.179015863439703082293818806944E+00;
    weight[9] = 0.179015863439703082293818806944E+00;
    weight[10] = 0.173262109489456226010614403827E+00;
    weight[11] = 0.161939517237602489264326706700E+00;
    weight[12] = 0.145411961573802267983003210494E+00;
    weight[13] = 0.124210533132967100263396358897E+00;
    weight[14] = 0.990162717175028023944236053187E-01;
    weight[15] = 0.706371668856336649992229601678E-01;
    weight[16] = 0.399706288109140661375991764101E-01;
    weight[17] = 0.653594771241830065359477124183E-02;
  }
  else if ( order == 19 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.978611766222080095152634063110E+00;
    xtab[2] =  - 0.928901528152586243717940258797E+00;
    xtab[3] =  - 0.852460577796646093085955970041E+00;
    xtab[4] =  - 0.751494202552613014163637489634E+00;
    xtab[5] =  - 0.628908137265220497766832306229E+00;
    xtab[6] =  - 0.488229285680713502777909637625E+00;
    xtab[7] =  - 0.333504847824498610298500103845E+00;
    xtab[8] =  - 0.169186023409281571375154153445E+00;
    xtab[9] =   0.0E+00;
    xtab[10] =   0.169186023409281571375154153445E+00;
    xtab[11] =   0.333504847824498610298500103845E+00;
    xtab[12] =   0.488229285680713502777909637625E+00;
    xtab[13] =   0.628908137265220497766832306229E+00;
    xtab[14] =   0.751494202552613014163637489634E+00;
    xtab[15] =   0.852460577796646093085955970041E+00;
    xtab[16] =   0.928901528152586243717940258797E+00;
    xtab[17] =   0.978611766222080095152634063110E+00;
    xtab[18] =   1.0E+00;

    weight[0] =  0.584795321637426900584795321637E-02;
    weight[1] =  0.357933651861764771154255690351E-01;
    weight[2] =  0.633818917626297368516956904183E-01;
    weight[3] =  0.891317570992070844480087905562E-01;
    weight[4] =  0.112315341477305044070910015464E+00;
    weight[5] =  0.132267280448750776926046733910E+00;
    weight[6] =  0.148413942595938885009680643668E+00;
    weight[7] =  0.160290924044061241979910968184E+00;
    weight[8] =  0.167556584527142867270137277740E+00;
    weight[9] = 0.170001919284827234644672715617E+00;
    weight[10] = 0.167556584527142867270137277740E+00;
    weight[11] = 0.160290924044061241979910968184E+00;
    weight[12] = 0.148413942595938885009680643668E+00;
    weight[13] = 0.132267280448750776926046733910E+00;
    weight[14] = 0.112315341477305044070910015464E+00;
    weight[15] = 0.891317570992070844480087905562E-01;
    weight[16] = 0.633818917626297368516956904183E-01;
    weight[17] = 0.357933651861764771154255690351E-01;
    weight[18] = 0.584795321637426900584795321637E-02;
  }
  else if ( order == 20 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.980743704893914171925446438584E+00;
    xtab[2] =  - 0.935934498812665435716181584931E+00;
    xtab[3] =  - 0.866877978089950141309847214616E+00;
    xtab[4] =  - 0.775368260952055870414317527595E+00;
    xtab[5] =  - 0.663776402290311289846403322971E+00;
    xtab[6] =  - 0.534992864031886261648135961829E+00;
    xtab[7] =  - 0.392353183713909299386474703816E+00;
    xtab[8] =  - 0.239551705922986495182401356927E+00;
    xtab[9] = - 0.805459372388218379759445181596E-01;
    xtab[10] =   0.805459372388218379759445181596E-01;
    xtab[11] =   0.239551705922986495182401356927E+00;
    xtab[12] =   0.392353183713909299386474703816E+00;
    xtab[13] =   0.534992864031886261648135961829E+00;
    xtab[14] =   0.663776402290311289846403322971E+00;
    xtab[15] =   0.775368260952055870414317527595E+00;
    xtab[16] =   0.866877978089950141309847214616E+00;
    xtab[17] =   0.935934498812665435716181584931E+00;
    xtab[18] =   0.980743704893914171925446438584E+00;
    xtab[19] =   1.0E+00;

    weight[0] =  0.526315789473684210526315789474E-02;
    weight[1] =  0.322371231884889414916050281173E-01;
    weight[2] =  0.571818021275668260047536271732E-01;
    weight[3] =  0.806317639961196031447768461137E-01;
    weight[4] =  0.101991499699450815683781205733E+00;
    weight[5] =  0.120709227628674725099429705002E+00;
    weight[6] =  0.136300482358724184489780792989E+00;
    weight[7] =  0.148361554070916825814713013734E+00;
    weight[8] =  0.156580102647475487158169896794E+00;
    weight[9] = 0.160743286387845749007726726449E+00;
    weight[10] = 0.160743286387845749007726726449E+00;
    weight[11] = 0.156580102647475487158169896794E+00;
    weight[12] = 0.148361554070916825814713013734E+00;
    weight[13] = 0.136300482358724184489780792989E+00;
    weight[14] = 0.120709227628674725099429705002E+00;
    weight[15] = 0.101991499699450815683781205733E+00;
    weight[16] = 0.806317639961196031447768461137E-01;
    weight[17] = 0.571818021275668260047536271732E-01;
    weight[18] = 0.322371231884889414916050281173E-01;
    weight[19] = 0.526315789473684210526315789474E-02;
  }
  else
  {
    cerr << "\n";
    cerr << "LOBATTO_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are between 2 and 20.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void moulton_set ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    MOULTON_SET sets weights for Adams-Moulton quadrature.
//
//  Discussion:
//
//    Adams-Moulton quadrature formulas are normally used in solving
//    ordinary differential equations, and are not suitable for general
//    quadrature computations.  However, an Adams-Moulton formula is
//    equivalent to approximating the integral of F(Y(X)) between X(M)
//    and X(M+1), using an implicit formula that relies on known values
//    of F(Y(X)) at X(M-N+1) through X(M), plus the unknown value at X(M+1).
//
//    Suppose the unknown function is denoted by Y(X), with derivative F(Y(X)),
//    and that approximate values of the function are known at a series of
//    X values, which we write as X(1), X(2), ..., X(M).  We write the value
//    Y(X(1)) as Y(1) and so on.
//
//    Then the solution of the ODE Y' = F(X,Y) at the next point X(M+1) is
//    computed by:
//
//      Y(M+0] = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
//             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+2-I)) 
//             approximately.
//
//    Note that this formula is implicit, since the unknown value Y(M+1)
//    appears on the right hand side.  Hence, in ODE applications, this
//    equation must be solved via a nonlinear equation solver.  For
//    quadrature problems, where the function to be integrated is known
//    beforehand, this is not a problem, and the calculation is explicit.
//
//    In the documentation that follows, we replace F(Y(X)) by F(X).
//
//
//    The Adams-Moulton formulas require equally spaced data.
//
//    Here is how the formula is applied in the case with non-unit spacing:
//
//      Integral ( A <= X <= A+H ) F(X) dX =
//      H * Sum ( 1 <= I <= ORDER ) weight[I) * F ( A - (I-2)*H ),
//      approximately.
//
//    The integral:
//
//      Integral ( 0 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) weight[I) * F ( 2 - I )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Jean Lapidus, John Seinfeld,
//    Numerical Solution of Ordinary Differential Equations,
//    Academic Press, 1971.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int ORDER, the order.  ORDER must be
//    between 1 and 10 or 12, 14, 16, 18 or 20.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//    weight[1) is the weight at X = 1, weight[2) the weight at X = 0, 
//    and so on.  The weights are rational.  The weights are not
//    symmetric, and some weights may be negative.  
//    They should sum to 1.
//
{
  double d;
  int i;

  if ( order == 1 )
  {
    weight[0] =  1.0;
  }
  else if ( order == 2 )
  {
    d = 2.0;

    weight[0] =  1.0 / d;
    weight[1] =  1.0 / d;
  }
  else if ( order == 3 )
  {
    d = 12.0;

    weight[0] =    5.0 / d;
    weight[1] =    8.0 / d;
    weight[2] =  - 1.0 / d;
  }
  else if ( order == 4 )
  {
    d = 24.0;

    weight[0] =    9.0 / d;
    weight[1] =   19.0 / d;
    weight[2] =  - 5.0 / d;
    weight[3] =    1.0 / d;
  }
  else if ( order == 5 )
  {
    d = 720.0;

    weight[0] =    251.0 / d;
    weight[1] =    646.0 / d;
    weight[2] =  - 264.0 / d;
    weight[3] =    106.0 / d;
    weight[4] =   - 19.0 / d;
  }
  else if ( order == 6 )
  {
    d = 1440.0;

    weight[0] =    475.0 / d;
    weight[1] =   1427.0 / d;
    weight[2] =  - 798.0 / d;
    weight[3] =    482.0 / d;
    weight[4] =  - 173.0 / d;
    weight[5] =     27.0 / d;
  }
  else if ( order == 7 )
  {
    d = 60480.0;

    weight[0] =    19087.0 / d;
    weight[1] =    65112.0 / d;
    weight[2] =  - 46461.0 / d;
    weight[3] =    37504.0 / d;
    weight[4] =  - 20211.0 / d;
    weight[5] =     6312.0 / d;
    weight[6] =    - 863.0 / d;
  }
  else if ( order == 8 )
  {
    d = 120960.0;

    weight[0] =    36799.0 / d;
    weight[1] =   139849.0 / d;
    weight[2] = - 121797.0 / d;
    weight[3] =   123133.0 / d;
    weight[4] =  - 88547.0 / d;
    weight[5] =    41499.0 / d;
    weight[6] =  - 11351.0 / d;
    weight[7] =     1375.0 / d;
  }
  else if ( order == 9 )
  {
    d = 3628800.0;

    weight[0] =   1070017.0 / d;
    weight[1] =   4467094.0 / d;
    weight[2] = - 4604594.0 / d;
    weight[3] =   5595358.0 / d;
    weight[4] = - 5033120.0 / d;
    weight[5] =   3146338.0 / d;
    weight[6] = - 1291214.0 / d;
    weight[7] =    312874.0 / d;
    weight[8] =   - 33953.0 / d;
  }
  else if ( order == 10 )
  {
    d = 7257600.0;

    weight[0] =    2082753.0 / d;
    weight[1] =    9449717.0 / d;
    weight[2] = - 11271304.0 / d;
    weight[3] =   16002320.0 / d;
    weight[4] = - 17283646.0 / d;
    weight[5] =   13510082.0 / d;
    weight[6] =  - 7394032.0 / d;
    weight[7] =    2687864.0 / d;
    weight[8] =   - 583435.0 / d;
    weight[9] =     57281.0 / d;
  }
  else if ( order == 12 )
  {
    d = 958003200.0E+00;

    weight[0] =    262747265.0E+00 / d;
    weight[1] =   1374799219.0E+00 / d;
    weight[2] =  -2092490673.0E+00 / d;
    weight[3] =   3828828885.0E+00 / d;
    weight[4] =  -5519460582.0E+00 / d;
    weight[5] =   6043521486.0E+00 / d;
    weight[6] =  -4963166514.0E+00 / d;
    weight[7] =   3007739418.0E+00 / d;
    weight[8] =  -1305971115.0E+00 / d;
    weight[9] =   384709327.0E+00 / d;
    weight[10] =   -68928781.0E+00 / d;
    weight[11] =     5675265.0E+00 / d;
  }
  else if ( order == 14 )
  {
    d = 5230697472000.0E+00;

    weight[0] =    1382741929621.0E+00 / d;
    weight[1] =    8153167962181.0E+00 / d;
    weight[2] =  -15141235084110.0E+00 / d;
    weight[3] =   33928990133618.0E+00 / d;
    weight[4] =  -61188680131285.0E+00 / d;
    weight[5] =   86180228689563.0E+00 / d;
    weight[6] =  -94393338653892.0E+00 / d;
    weight[7] =   80101021029180.0E+00 / d;
    weight[8] =  -52177910882661.0E+00 / d;
    weight[9] =  25620259777835.0E+00 / d;
    weight[10] =  -9181635605134.0E+00 / d;
    weight[11] =   2268078814386.0E+00 / d;
    weight[12] =   -345457086395.0E+00 / d;
    weight[13] =     24466579093.0E+00 / d;
  }
  else if ( order == 16 )
  {
    d = 62768369664000.0E+00;

    weight[0] =     16088129229375.0E+00 / d;
    weight[1] =    105145058757073.0E+00 / d;
    weight[2] =   -230992163723849.0E+00 / d;
    weight[3] =    612744541065337.0E+00 / d;
    weight[4] =  -1326978663058069.0E+00 / d;
    weight[5] =   2285168598349733.0E+00 / d;
    weight[6] =  -3129453071993581.0E+00 / d;
    weight[7] =   3414941728852893.0E+00 / d;
    weight[8] =  -2966365730265699.0E+00 / d;
    weight[9] =  2039345879546643.0E+00 / d;
    weight[10] = -1096355235402331.0E+00 / d;
    weight[11] =   451403108933483.0E+00 / d;
    weight[12] =  -137515713789319.0E+00 / d;
    weight[13] =    29219384284087.0E+00 / d;
    weight[14] =    -3867689367599.0E+00 / d;
    weight[15] =      240208245823.0E+00 / d;
  }
  else if ( order == 18 )
  {
    d = 64023737057280000.0E+00;
       
    weight[0] =      15980174332775873.0E+00 / d;
    weight[1] =     114329243705491117.0E+00 / d;
    weight[2] =    -290470969929371220.0E+00 / d;
    weight[3] =     890337710266029860.0E+00 / d;
    weight[4] =   -2250854333681641520.0E+00 / d;
    weight[5] =    4582441343348851896.0E+00 / d;
    weight[6] =   -7532171919277411636.0E+00 / d;
    weight[7] =   10047287575124288740.0E+00 / d;
    weight[8] =  -10910555637627652470.0E+00 / d;
    weight[9] =   9644799218032932490.0E+00 / d;
    weight[10] =  -6913858539337636636.0E+00 / d;
    weight[11] =   3985516155854664396.0E+00 / d;
    weight[12] =  -1821304040326216520.0E+00 / d;
    weight[13] =    645008976643217360.0E+00 / d;
    weight[14] =   -170761422500096220.0E+00 / d;
    weight[15] =     31816981024600492.0E+00 / d;
    weight[16] =     -3722582669836627.0E+00 / d;
    weight[17] =       205804074290625.0E+00 / d;
  }
  else if ( order == 20 )
  {
    d = 102181884343418880000.0E+00;      

    weight[0] =       24919383499187492303.0E+00 / d;
    weight[1] =      193280569173472261637.0E+00 / d;
    weight[2] =     -558160720115629395555.0E+00 / d;
    weight[3] =     1941395668950986461335.0E+00 / d;
    weight[4] =    -5612131802364455926260.0E+00 / d;
    weight[5] =    13187185898439270330756.0E+00 / d;
    weight[6] =   -25293146116627869170796.0E+00 / d;
    weight[7] =    39878419226784442421820.0E+00 / d;
    weight[8] =   -51970649453670274135470.0E+00 / d;
    weight[9] =   56154678684618739939910.0E+00 / d;
    weight[10] =  -50320851025594566473146.0E+00 / d;
    weight[11] =   37297227252822858381906.0E+00 / d;
    weight[12] =  -22726350407538133839300.0E+00 / d;
    weight[13] =   11268210124987992327060.0E+00 / d;
    weight[14] =   -4474886658024166985340.0E+00 / d;
    weight[15] =    1389665263296211699212.0E+00 / d;
    weight[16] =    -325187970422032795497.0E+00 / d;
    weight[17] =      53935307402575440285.0E+00 / d;
    weight[18] =      -5652892248087175675.0E+00 / d;
    weight[19] =        281550972898020815.0E+00 / d;
  }
  else
  {
    cerr << "\n";
    cerr << "MOULTON_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are 1 through 10, 12, 14, 16, 18, or 20.\n";
    exit ( 1 );
  }
  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( double ) ( 1 - i );
  }
  return;
}
//****************************************************************************80

void nc_compute ( int n, double x_min, double x_max, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NC_COMPUTE computes a Newton-Cotes quadrature rule.
//
//  Discussion:
//
//    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule 
//    estimates
//
//      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
//
//    using N abscissas X and weights W:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
//
//    For the CLOSED rule, the abscissas include the end points.
//    For the OPEN rule, the abscissas do not include the end points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, double X_MIN, X_MAX, the endpoints of the interval.
//
//    Input, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *d;
  int i;
  int j;
  int k;
  double yvala;
  double yvalb;

  d = new double[n];

  for ( i = 1; i <= n; i++ )
  {
//
//  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
//  and zero at the other nodes.
//
    for ( j = 1; j <= n; j++ )
    {
      d[j-1] = 0.0;
    }
    d[i-1] = 1.0;

    for ( j = 2; j <= n; j++ )
    {
      for ( k = j; k <= n; k++ )
      {
        d[n+j-k-1] = ( d[n+j-k-2] - d[n+j-k-1] ) / ( x[n-k] - x[n+j-k-1] );
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( k = 1; k <= n - j; k++ )
      {
        d[n-k-1] = d[n-k-1] - x[n-k-j] * d[n-k];
      }
    }
//
//  Evaluate the antiderivative of the polynomial at the left and
//  right endpoints.
//
    yvala = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvala = yvala * x_min + d[j] / ( double ) ( j + 1 );
    }
    yvala = yvala * x_min;

    yvalb = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvalb = yvalb * x_max + d[j] / ( double ) ( j + 1 );
    }
    yvalb = yvalb * x_max;

    w[i-1] = yvalb - yvala;
  }

  delete [] d;

  return;
}
//****************************************************************************80

void ncc_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_COMPUTE computes a Newton-Cotes closed quadrature rule.
//
//  Discussion:
//
//    For the interval [-1.0,+1.0], the Newton-Cotes closed quadrature rule
//    estimates
//
//      Integral ( -1 <= X <= +1 ) F(X) dX
//
//    using N abscissas X and weights W:
//
//      sum ( 1 <= I <= N ) W[I] * F ( X[I] ).
//
//    For the CLOSED rule, the abscissas are equally spaced and include 
//    the end points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  ncc_compute_points ( n, x );

  ncc_compute_weights ( n, w );

  return;
}
//****************************************************************************80

void ncc_compute_points ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_COMPUTE_POINTS: points of a Newton-Cotes Closed quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;

  if ( n == 1 )
  {
    x[0] = ( x_max + x_min ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * x_min   
             + ( double ) (     i     ) * x_max ) 
             / ( double ) ( n     - 1 );
    }
  }
  return;
}
//****************************************************************************80

void ncc_compute_weights ( int n, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_COMPUTE_WEIGHTS: weights of a Newton-Cotes Closed quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double W[N], the weights.
//
{
  int i;
  double *x;
  double x_max = 1.0;
  double x_min = -1.0;

  if ( n == 1 )
  {
    w[0] = x_max - x_min;
  }
  else
  {
    x = new double[n];

    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * x_min   
             + ( double ) (     i     ) * x_max ) 
             / ( double ) ( n     - 1 );
    }
    nc_compute ( n, x_min, x_max, x, w );

    delete [] x;
  }
  return;
}
//****************************************************************************80

void ncc_set ( int order, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_SET sets abscissas and weights for closed Newton-Cotes quadrature.
//
//  Discussion:
//
//    The closed Newton-Cotes rules use equally spaced abscissas, and
//    hence may be used with tabulated function data.
//
//    The rules are called "closed" because they include the endpoints.
//    As a favor, we include an order 1 rule, the midpoint rule, even 
//    though this does not satisfy the requirement that the endpoints
//    be included//
//
//    The higher order rules involve negative weights.  These can produce
//    loss of accuracy due to the subtraction of large, nearly equal quantities.
//
//    ORDER = 1 is the midpoint rule (and is not really an NCC rule//)
//    ORDER = 2 is the trapezoidal rule.
//    ORDER = 3 is Simpson's rule.
//    ORDER = 4 is Simpson's 3/8 rule.
//    ORDER = 5 is Bode's rule.
//
//    The Kopal reference for ORDER = 12 lists
//      W(6) = 15494566.0D+00 / 43545600.0D+00
//    but this results in a set of coeffients that don't add up to 2.
//    The correct value is
//      W(6) = 15493566.0D+00 / 43545600.0.
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//    In Mathematica, the closed Newton-Cotes weights
//    can be computed by:
//
//      Needs["NumericalDifferentialEquationAnalysis`"]
//      NewtonCotesWeights [ order, -1, 1, QuadratureType -> Closed ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Johnson,
//    Quarterly Journal of Mathematics,
//    Volume 46, Number 52, 1915.
//
//    Zdenek Kopal,
//    Numerical Analysis,
//    John Wiley, 1955,
//    LC: QA297.K6.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 21.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double W[ORDER], the weights.
//
{
  if ( order == 1 )
  {
//
//  2
//
    x[0] = 0.00000000000000000000;

    w[0] = 2.00000000000000000000;
  }
  else if ( order == 2 )
  {
//
//  1
//  1
//
    x[0] = -1.00000000000000000000;
    x[1] =  1.00000000000000000000;

    w[0] = 1.00000000000000000000;
    w[1] = 1.00000000000000000000;
  }
  else if ( order == 3 )
  {
//
//  1 / 3
//  4 / 3
//  1 / 3
//
    x[0] = -1.00000000000000000000;
    x[1] =  0.00000000000000000000;
    x[2] =  1.00000000000000000000;

    w[0] = 0.33333333333333333333;
    w[1] = 1.33333333333333333333;
    w[2] = 0.33333333333333333333;
  }
  else if ( order == 4 )
  {
//
//  1 / 4
//  3 / 4
//  3 / 4
//  1 / 4
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.33333333333333333333;
    x[2] =  0.33333333333333333333;
    x[3] =  1.00000000000000000000;

    w[0] = 0.25000000000000000000;
    w[1] = 0.75000000000000000000;
    w[2] = 0.75000000000000000000;
    w[3] = 0.25000000000000000000;
  }
  else if ( order == 5 )
  {
//
//   7 / 45
//  32 / 45
//  12 / 45
//  32 / 45
//   7 / 45
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.50000000000000000000;
    x[2] =  0.00000000000000000000;
    x[3] =  0.50000000000000000000;
    x[4] =  1.00000000000000000000;

    w[0] = 0.15555555555555555556;
    w[1] = 0.71111111111111111111;
    w[2] = 0.26666666666666666667;
    w[3] = 0.71111111111111111111;
    w[4] = 0.15555555555555555556;
  }
  else if ( order == 6 )
  {
//
//  19 / 144
//  75 / 144
//  50 / 144
//  50 / 144
//  75 / 144
//  19 / 144
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.60000000000000000000;
    x[2] = -0.20000000000000000000;
    x[3] =  0.20000000000000000000;
    x[4] =  0.60000000000000000000;
    x[5] =  1.00000000000000000000;

    w[0] = 0.13194444444444444444;
    w[1] = 0.52083333333333333333;
    w[2] = 0.34722222222222222222;
    w[3] = 0.34722222222222222222;
    w[4] = 0.52083333333333333333;
    w[5] = 0.13194444444444444444;
  }
  else if ( order == 7 )
  {
//
//   41 / 420
//  216 / 420
//   27 / 420
//  272 / 420
//   27 / 420
//  216 / 420
//   41 / 420
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.66666666666666666667;
    x[2] = -0.33333333333333333333;
    x[3] =  0.00000000000000000000;
    x[4] =  0.33333333333333333333;
    x[5] =  0.66666666666666666667;
    x[6] =  1.00000000000000000000;

    w[0] = 0.097619047619047619048;
    w[1] = 0.51428571428571428571;
    w[2] = 0.064285714285714285714;
    w[3] = 0.64761904761904761905;
    w[4] = 0.064285714285714285714;
    w[5] = 0.51428571428571428571;
    w[6] = 0.097619047619047619048;
  }
  else if ( order == 8 )
  {
//
//   751 / 8640
//  3577 / 8640
//  1323 / 8640
//  2989 / 8640
//  2989 / 8640
//  1323 / 8640
//  3577 / 8640
//   751 / 8640
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.71428571428571428571;
    x[2] = -0.42857142857142857143;
    x[3] = -0.14285714285714285714;
    x[4] =  0.14285714285714285714;
    x[5] =  0.42857142857142857143;
    x[6] =  0.71428571428571428571;
    x[7] =  1.00000000000000000000;

    w[0] = 0.086921296296296296296;
    w[1] = 0.41400462962962962963;
    w[2] = 0.15312500000000000000;
    w[3] = 0.34594907407407407407;
    w[4] = 0.34594907407407407407;
    w[5] = 0.15312500000000000000;
    w[6] = 0.41400462962962962963;
    w[7] = 0.086921296296296296296;
  }
  else if ( order == 9 )
  {
//
//    989 / 14175
//   5888 / 14175
//   -928 / 14175
//  10496 / 14175
//  -4540 / 14175
//  10496 / 14175
//   -928 / 14175
//   5888 / 14175
//    989 / 14175
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.75000000000000000000;
    x[2] = -0.50000000000000000000;
    x[3] = -0.25000000000000000000;
    x[4] =  0.00000000000000000000;
    x[5] =  0.25000000000000000000;
    x[6] =  0.50000000000000000000;
    x[7] =  0.75000000000000000000;
    x[8] =  1.00000000000000000000;

    w[0] =  0.069770723104056437390;
    w[1] =  0.41537918871252204586;
    w[2] = -0.065467372134038800705;
    w[3] =  0.74045855379188712522;
    w[4] = -0.32028218694885361552;
    w[5] =  0.74045855379188712522;
    w[6] = -0.065467372134038800705;
    w[7] =  0.41537918871252204586;
    w[8] =  0.069770723104056437390;
  }
  else if ( order == 10 )
  {
//
//   2857 / 44800
//  15741 / 44800
//   1080 / 44800
//  19344 / 44800
//   5778 / 44800
//   5778 / 44800
//  19344 / 44800
//   1080 / 44800
//  15741 / 44800
//   2857 / 44800
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.77777777777777777778;
    x[2] = -0.55555555555555555556;
    x[3] = -0.33333333333333333333;
    x[4] = -0.11111111111111111111;
    x[5] =  0.11111111111111111111;
    x[6] =  0.33333333333333333333;
    x[7] =  0.55555555555555555556;
    x[8] =  0.77777777777777777778;
    x[9] =  1.00000000000000000000;

    w[0] =  0.063772321428571428571;
    w[1] =  0.35136160714285714286;
    w[2] =  0.024107142857142857143;
    w[3] =  0.43178571428571428571;
    w[4] =  0.12897321428571428571;
    w[5] =  0.12897321428571428571;
    w[6] =  0.43178571428571428571;
    w[7] =  0.024107142857142857143;
    w[8] =  0.35136160714285714286;
    w[9] =  0.063772321428571428571;
  }
  else if ( order == 11 )
  {
//
//     16067 / 299376
//    106300 / 299376
//   - 48525 / 299376
//    272400 / 299376
//  - 260550 / 299376
//    427368 / 299376
//  - 260550 / 299376
//    272400 / 299376
//   - 48525 / 299376
//    106300 / 299376
//     16067 / 299376
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.80000000000000000000;
    x[2] = -0.60000000000000000000;
    x[3] = -0.40000000000000000000;
    x[4] = -0.20000000000000000000;
    x[5] =  0.00000000000000000000;
    x[6] =  0.20000000000000000000;
    x[7] =  0.40000000000000000000;
    x[8] =  0.60000000000000000000;
    x[9] =  0.80000000000000000000;
    x[10] = 1.00000000000000000000;

    w[0] =  0.053668296723852279408;
    w[1] =  0.35507188284966062744;
    w[2] = -0.16208714125380792047;
    w[3] =  0.90989257655924322591;
    w[4] = -0.87031024531024531025;
    w[5] =  1.4275292608625941959;
    w[6] = -0.87031024531024531025;
    w[7] =  0.90989257655924322591;
    w[8] = -0.16208714125380792047;
    w[9] =  0.35507188284966062744;
    w[10] = 0.053668296723852279408;
  }
  else if ( order == 12 )
  {
//
//     2171465 / 43545600
//    13486539 / 43545600
//   - 3237113 / 43545600
//    25226685 / 43545600
//   - 9595542 / 43545600
//    15493566 / 43545600
//    15493566 / 43545600
//   - 9595542 / 43545600
//    25226685 / 43545600
//   - 3237113 / 43545600
//    13486539 / 43545600
//     2171465 / 43545600
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.81818181818181818182;
    x[2] = -0.63636363636363636364;
    x[3] = -0.45454545454545454545;
    x[4] = -0.27272727272727272727;
    x[5] = -0.090909090909090909091;
    x[6] =  0.090909090909090909091;
    x[7] =  0.27272727272727272727;
    x[8] =  0.45454545454545454545;
    x[9] =  0.63636363636363636364;
    x[10] = 0.81818181818181818182;
    x[11] = 1.00000000000000000000;

    w[0] =   0.049866461823927101705;
    w[1] =   0.30971071704144620811;
    w[2] =  -0.074338463587595532040;
    w[3] =   0.57931650958994708995;
    w[4] =  -0.22035617835097001764;
    w[5] =   0.35580095348324514991;
    w[6] =   0.35580095348324514991;
    w[7] =  -0.22035617835097001764;
    w[8] =   0.57931650958994708995;
    w[9] =  -0.074338463587595532040;
    w[10] =  0.30971071704144620811;
    w[11] =  0.049866461823927101705;
  }
  else if ( order == 13 )
  {
//
//      1364651 / 31531500
//      9903168 / 31531500
//    - 7587864 / 31531500
//     35725120 / 31531500
//   - 51491295 / 31531500
//     87516288 / 31531500
//   - 87797136 / 31531500
//     87516288 / 31531500
//   - 51491295 / 31531500
//     35725120 / 31531500
//    - 7587864 / 31531500
//      9903168 / 31531500
//      1364651 / 31531500
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.83333333333333333333;
    x[2] = -0.66666666666666666667;
    x[3] = -0.50000000000000000000;
    x[4] = -0.33333333333333333333;
    x[5] = -0.16666666666666666667;
    x[6] =  0.00000000000000000000;
    x[7] =  0.16666666666666666667;
    x[8] =  0.33333333333333333333;
    x[9] =  0.50000000000000000000;
    x[10] = 0.66666666666666666667;
    x[11] = 0.83333333333333333333;
    x[12] = 1.00000000000000000000;

    w[0] =   0.043278974993260707546;
    w[1] =   0.31407221350078492936;
    w[2] =  -0.24064392750107035821;
    w[3] =   1.1329977958549387121;
    w[4] =  -1.6330112744398458684;
    w[5] =   2.7755193378050520908;
    w[6] =  -2.7844262404262404262;
    w[7] =   2.7755193378050520908;
    w[8] =  -1.6330112744398458684;
    w[9] =   1.1329977958549387121;
    w[10] = -0.24064392750107035821;
    w[11] =  0.31407221350078492936;
    w[12] =  0.043278974993260707546;
  }
  else if ( order == 14 )
  {
//
//      6137698213 / 150885504000
//     42194238652 / 150885504000
//   - 23361540993 / 150885504000
//    116778274403 / 150885504000
//  - 113219777650 / 150885504000
//    154424590209 / 150885504000
//   - 32067978834 / 150885504000
//   - 32067978834 / 150885504000
//    154424590209 / 150885504000
//  - 113219777650 / 150885504000
//    116778274403 / 150885504000
//   - 23361540993 / 150885504000
//     42194238652 / 150885504000
//      6137698213 / 150885504000
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.84615384615384615385;
    x[2] = -0.69230769230769230769;
    x[3] = -0.53846153846153846154;
    x[4] = -0.38461538461538461538;
    x[5] = -0.23076923076923076923;
    x[6] = -0.076923076923076923077;
    x[7] =  0.076923076923076923077;
    x[8] =  0.23076923076923076923;
    x[9] =  0.38461538461538461538;
    x[10] = 0.53846153846153846154;
    x[11] = 0.69230769230769230769;
    x[12] = 0.84615384615384615385;
    x[13] = 1.00000000000000000000;

    w[0] =   0.040669438210247155353;
    w[1] =   0.27975217053157074652;
    w[2] =  -0.15542374057682837445;
    w[3] =   0.77579230848776566369;
    w[4] =  -0.75384763266423526013;
    w[5] =   1.0273523591123107492;
    w[6] =  -0.21429490310083068020;
    w[7] =  -0.21429490310083068020;
    w[8] =   1.0273523591123107492;
    w[9] =  -0.75384763266423526013;
    w[10] =  0.77579230848776566369;
    w[11] = -0.15542374057682837445;
    w[12] =  0.27975217053157074652;
    w[13] =  0.040669438210247155353;
  }
  else if ( order == 15 )
  {
//
//       90241897 / 2501928000
//      710986864 / 2501928000
//    - 770720657 / 2501928000
//     3501442784 / 2501928000
//   - 6625093363 / 2501928000
//    12630121616 / 2501928000
//  - 16802270373 / 2501928000
//    19534438464 / 2501928000
//  - 16802270373 / 2501928000
//    12630121616 / 2501928000
//   - 6625093363 / 2501928000
//     3501442784 / 2501928000
//    - 770720657 / 2501928000
//      710986864 / 2501928000
//       90241897 / 2501928000
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.85714285714285714286;
    x[2] = -0.71428571428571428571;
    x[3] = -0.57142857142857142857;
    x[4] = -0.42857142857142857143;
    x[5] = -0.28571428571428571429;
    x[6] = -0.14285714285714285714;
    x[7] =  0.00000000000000000000;
    x[8] =  0.14285714285714285714;
    x[9] =  0.28571428571428571429;
    x[10] = 0.42857142857142857143;
    x[11] = 0.57142857142857142857;
    x[12] = 0.71428571428571428571;
    x[13] = 0.85714285714285714286;
    x[14] = 1.00000000000000000000;

    w[0] =   0.036068942431596752584;
    w[1] =   0.28417558938546592868;
    w[2] =  -0.30805069410470645039;
    w[3] =   1.3994978208805369299;
    w[4] =  -2.6479952112930507992;
    w[5] =   5.0481555088715582543;
    w[6] =  -6.7157289790113864188;
    w[7] =   7.8077540456799716059;
    w[8] =  -6.7157289790113864188;
    w[9] =   5.0481555088715582543;
    w[10] = -2.6479952112930507992;
    w[11] =  1.3994978208805369299;
    w[12] = -0.30805069410470645039;
    w[13] =  0.28417558938546592868;
    w[14] =  0.036068942431596752584;
  }
  else if ( order == 16 )
  {
//
//     105930069 / 3099672576
//     796661595 / 3099672576
//   - 698808195 / 3099672576
//    3143332755 / 3099672576
//  - 4688522055 / 3099672576
//    7385654007 / 3099672576
//  - 6000998415 / 3099672576
//    3056422815 / 3099672576
//    3056422815 / 3099672576
//  - 6000998415 / 3099672576
//    7385654007 / 3099672576
//  - 4688522055 / 3099672576
//    3143332755 / 3099672576
//   - 698808195 / 3099672576
//     796661595 / 3099672576
//     105930069 / 3099672576
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.86666666666666666667;
    x[2] = -0.73333333333333333333;
    x[3] = -0.60000000000000000000;
    x[4] = -0.46666666666666666667;
    x[5] = -0.33333333333333333333;
    x[6] = -0.20000000000000000000;
    x[7] = -0.066666666666666666667;
    x[8] =  0.066666666666666666667;
    x[9] =  0.20000000000000000000;
    x[10] = 0.33333333333333333333;
    x[11] = 0.46666666666666666667;
    x[12] = 0.60000000000000000000;
    x[13] = 0.73333333333333333333;
    x[14] = 0.86666666666666666667;
    x[15] = 1.00000000000000000000;

    w[0] =   0.034174599543251887002;
    w[1] =   0.25701475735481036820;
    w[2] =  -0.22544581011901045383;
    w[3] =   1.0140854164204471124;
    w[4] =  -1.5125862296882804695;
    w[5] =   2.3827206990135980091;
    w[6] =  -1.9360104229924960952;
    w[7] =   0.98604699046767964179;
    w[8] =   0.98604699046767964179;
    w[9] =  -1.9360104229924960952;
    w[10] =  2.3827206990135980091;
    w[11] = -1.5125862296882804695;
    w[12] =  1.0140854164204471124;
    w[13] = -0.22544581011901045383;
    w[14] =  0.25701475735481036820;
    w[15] =  0.034174599543251887002;
  }
  else if ( order == 17 )
  {
//
//       15043611773 / 488462349375
//      127626606592 / 488462349375
//    - 179731134720 / 488462349375
//      832211855360 / 488462349375
//   - 1929498607520 / 488462349375
//     4177588893696 / 488462349375
//   - 6806534407936 / 488462349375
//     9368875018240 / 488462349375
//  - 10234238972220 / 488462349375
//     9368875018240 / 488462349375
//   - 6806534407936 / 488462349375
//     4177588893696 / 488462349375
//   - 1929498607520 / 488462349375
//      832211855360 / 488462349375
//    - 179731134720 / 488462349375
//      127626606592 / 488462349375
//       15043611773 / 488462349375
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.87500000000000000000;
    x[2] = -0.75000000000000000000;
    x[3] = -0.62500000000000000000;
    x[4] = -0.50000000000000000000;
    x[5] = -0.37500000000000000000;
    x[6] = -0.25000000000000000000;
    x[7] = -0.12500000000000000000;
    x[8] =  0.00000000000000000000;
    x[9] =  0.12500000000000000000;
    x[10] = 0.25000000000000000000;
    x[11] = 0.37500000000000000000;
    x[12] = 0.50000000000000000000;
    x[13] = 0.62500000000000000000;
    x[14] = 0.75000000000000000000;
    x[15] = 0.87500000000000000000;
    x[16] = 1.00000000000000000000;

    w[0]  =   0.030797894233299012495;
    w[1]  =   0.26128238288028031086;
    w[2]  =  -0.36795289329867605622;
    w[3]  =   1.7037379778090086905;
    w[4]  =  -3.9501480717783930427;
    w[5]  =   8.5525299934402953388;
    w[6]  = -13.934614237197880038;
    w[7]  =  19.180342211078732848;
    w[8]  = -20.951950514333334128;
    w[9] =   19.180342211078732848;
    w[10] = -13.934614237197880038;
    w[11] =   8.5525299934402953388;
    w[12] =  -3.9501480717783930427;
    w[13] =   1.7037379778090086905;
    w[14] =  -0.36795289329867605622;
    w[15] =   0.26128238288028031086;
    w[16] =   0.030797894233299012495;
  }
  else if ( order == 18 )
  {
//
//       55294720874657 / 1883051089920000
//      450185515446285 / 1883051089920000
//    - 542023437008852 / 1883051089920000
//     2428636525764260 / 1883051089920000
//   - 4768916800123440 / 1883051089920000
//     8855416648684984 / 1883051089920000
//  - 10905371859796660 / 1883051089920000
//    10069615750132836 / 1883051089920000
//   - 3759785974054070 / 1883051089920000
//   - 3759785974054070 / 1883051089920000
//    10069615750132836 / 1883051089920000
//  - 10905371859796660 / 1883051089920000
//     8855416648684984 / 1883051089920000
//   - 4768916800123440 / 1883051089920000
//     2428636525764260 / 1883051089920000
//    - 542023437008852 / 1883051089920000
//      450185515446285 / 1883051089920000
//       55294720874657 / 1883051089920000
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.88235294117647058824;
    x[2] = -0.76470588235294117647;
    x[3] = -0.64705882352941176471;
    x[4] = -0.52941176470588235294;
    x[5] = -0.41176470588235294118;
    x[6] = -0.29411764705882352941;
    x[7] = -0.17647058823529411765;
    x[8] = -0.058823529411764705882;
    x[9] =  0.058823529411764705882;
    x[10] = 0.17647058823529411765;
    x[11] = 0.29411764705882352941;
    x[12] = 0.41176470588235294118;
    x[13] = 0.52941176470588235294;
    x[14] = 0.64705882352941176471;
    x[15] = 0.76470588235294117647;
    x[16] = 0.88235294117647058824;
    x[17] = 1.00000000000000000000;

    w[0] =   0.029364429446790078519;
    w[1] =   0.23907238516051669677;
    w[2] =  -0.28784319231183443641;
    w[3] =   1.2897348026109258587;
    w[4] =  -2.5325477495812627261;
    w[5] =   4.7026959045817496499;
    w[6] =  -5.7913308450170443690;
    w[7] =   5.3475000248456540826;
    w[8] =  -1.9966457597354948350;
    w[9] =  -1.9966457597354948350;
    w[10] =  5.3475000248456540826;
    w[11] = -5.7913308450170443690;
    w[12] =  4.7026959045817496499;
    w[13] = -2.5325477495812627261;
    w[14] =  1.2897348026109258587;
    w[15] = -0.28784319231183443641;
    w[16] =  0.23907238516051669677;
    w[17] =  0.029364429446790078519;
  }
  else if ( order == 19 )
  {
//
//       203732352169 / 7604556960000
//      1848730221900 / 7604556960000
//    - 3212744374395 / 7604556960000
//     15529830312096 / 7604556960000
//   - 42368630685840 / 7604556960000
//    103680563465808 / 7604556960000
//  - 198648429867720 / 7604556960000
//    319035784479840 / 7604556960000
//  - 419127951114198 / 7604556960000
//    461327344340680 / 7604556960000
//  - 419127951114198 / 7604556960000
//    319035784479840 / 7604556960000
//  - 198648429867720 / 7604556960000
//    103680563465808 / 7604556960000
//   - 42368630685840 / 7604556960000
//     15529830312096 / 7604556960000
//    - 3212744374395 / 7604556960000
//      1848730221900 / 7604556960000
//       203732352169 / 7604556960000
//
    x[0] = -1.00000000000000000000;
    x[1] = -0.88888888888888888889;
    x[2] = -0.77777777777777777778;
    x[3] = -0.66666666666666666667;
    x[4] = -0.55555555555555555556;
    x[5] = -0.44444444444444444444;
    x[6] = -0.33333333333333333333;
    x[7] = -0.22222222222222222222;
    x[8] = -0.11111111111111111111;
    x[9] =  0.00000000000000000000;
    x[10] = 0.11111111111111111111;
    x[11] = 0.22222222222222222222;
    x[12] = 0.33333333333333333333;
    x[13] = 0.44444444444444444444;
    x[14] = 0.55555555555555555556;
    x[15] = 0.66666666666666666667;
    x[16] = 0.77777777777777777778;
    x[17] = 0.88888888888888888889;
    x[18] = 1.00000000000000000000;

    w[0] =    0.026790824664820447344;
    w[1] =    0.24310820888374278151;
    w[2] =   -0.42247620621346493274;
    w[3] =    2.0421742376029227612;
    w[4] =   -5.5714791681749728126;
    w[5] =   13.634004454324976218;
    w[6] =  -26.122288374274995239;
    w[7] =   41.953237533490708445;
    w[8] =  -55.115367445968607749;
    w[9] =   60.664591871329740161;
    w[10] = -55.115367445968607749;
    w[11] =  41.953237533490708445;
    w[12] = -26.122288374274995239;
    w[13] =  13.634004454324976218;
    w[14] =  -5.5714791681749728126;
    w[15] =   2.0421742376029227612;
    w[16] =  -0.42247620621346493274;
    w[17] =   0.24310820888374278151;
    w[18] =   0.026790824664820447344;
  }
  else if ( order == 20 )
  {
//
//       69028763155644023 / 2688996956405760000
//      603652082270808125 / 2688996956405760000
//    - 926840515700222955 / 2688996956405760000
//     4301581538450500095 / 2688996956405760000
//  - 10343692234243192788 / 2688996956405760000
//    22336420328479961316 / 2688996956405760000
//  - 35331888421114781580 / 2688996956405760000
//    43920768370565135580 / 2688996956405760000
//  - 37088370261379851390 / 2688996956405760000
//    15148337305921759574 / 2688996956405760000
//    15148337305921759574 / 2688996956405760000
//  - 37088370261379851390 / 2688996956405760000
//    43920768370565135580 / 2688996956405760000
//  - 35331888421114781580 / 2688996956405760000
//    22336420328479961316 / 2688996956405760000
//  - 10343692234243192788 / 2688996956405760000
//     4301581538450500095 / 2688996956405760000
//    - 926840515700222955 / 2688996956405760000
//      603652082270808125 / 2688996956405760000
//       69028763155644023 / 2688996956405760000
//
    x[0] =  -1.00000000000000000000;
    x[1] =  -0.89473684210526315789;
    x[2] =  -0.78947368421052631579;
    x[3] =  -0.68421052631578947368;
    x[4] =  -0.57894736842105263158;
    x[5] =  -0.47368421052631578947;
    x[6] =  -0.36842105263157894737;
    x[7] =  -0.26315789473684210526;
    x[8] =  -0.15789473684210526316;
    x[9] =  -0.052631578947368421053;
    x[10] =  0.052631578947368421053;
    x[11] =  0.15789473684210526316;
    x[12] =  0.26315789473684210526;
    x[13] =  0.36842105263157894737;
    x[14] =  0.47368421052631578947;
    x[15] =  0.57894736842105263158;
    x[16] =  0.68421052631578947368;
    x[17] =  0.78947368421052631579;
    x[18] =  0.89473684210526315789;
    x[19] =  1.00000000000000000000;

    w[0] =    0.025670822345560078100;
    w[1] =    0.22448968595251886556;
    w[2] =   -0.34467890099030890987;
    w[3] =    1.5996974366978074270;
    w[4] =   -3.8466730910952978835;
    w[5] =    8.3065993344729824120;
    w[6] =  -13.139430424771119113;
    w[7] =   16.333513604742678295;
    w[8] =  -13.792641220001198577;
    w[9] =    5.6334527526463774045;
    w[10] =   5.6334527526463774045;
    w[11] = -13.792641220001198577;
    w[12] =  16.333513604742678295;
    w[13] = -13.139430424771119113;
    w[14] =   8.3065993344729824120;
    w[15] =  -3.8466730910952978835;
    w[16] =   1.5996974366978074270;
    w[17] =  -0.34467890099030890987;
    w[18] =   0.22448968595251886556;
    w[19] =   0.025670822345560078100;
  }
  else if ( order == 21 )
  {
    x[0] =  -1.00000000000000000000;
    x[1] =  -0.90000000000000000000;
    x[2] =  -0.80000000000000000000;
    x[3] =  -0.70000000000000000000;
    x[4] =  -0.60000000000000000000;
    x[5] =  -0.50000000000000000000;
    x[6] =  -0.40000000000000000000;
    x[7] =  -0.30000000000000000000;
    x[8] =  -0.20000000000000000000;
    x[9] =  -0.10000000000000000000;
    x[10] =  0.00000000000000000000;
    x[11] =  0.10000000000000000000;
    x[12] =  0.20000000000000000000;
    x[13] =  0.30000000000000000000;
    x[14] =  0.40000000000000000000;
    x[15] =  0.50000000000000000000;
    x[16] =  0.60000000000000000000;
    x[17] =  0.70000000000000000000;
    x[18] =  0.80000000000000000000;
    x[19] =  0.90000000000000000000;
    x[20] =  1.00000000000000000000;

    w[0] =     0.023650546498063206389;
    w[1] =     0.22827543528921394997;
    w[2] =    -0.47295674102285392846;
    w[3] =     2.4123737869637513288;
    w[4] =    -7.5420634534306609355;
    w[5] =    20.673596439879602287;
    w[6] =   -45.417631687959024596;
    w[7] =    83.656114844387109207;
    w[8] =  -128.15055898030800930;
    w[9] =   165.59456694494570344;
    w[10] = -180.01073427048578932;
    w[11] =  165.59456694494570344;
    w[12] = -128.15055898030800930;
    w[13] =   83.656114844387109207;
    w[14] =  -45.417631687959024596;
    w[15] =   20.673596439879602287;
    w[16] =   -7.5420634534306609355;
    w[17] =    2.4123737869637513288;
    w[18] =   -0.47295674102285392846;
    w[19] =    0.22827543528921394997;
    w[20] =    0.023650546498063206389;
  }
  else
  {
    cerr << "\n";
    cerr << "NCC_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are 1 through 21.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void nco_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_COMPUTE computes a Newton-Cotes Open quadrature rule.
//
//  Discussion:
//
//    For the interval [-1,+1], the Newton-Cotes Open quadrature rule
//    estimates
//
//      Integral ( -1 <= X <= +1 ) F(X) dX
//
//    using N equally spaced abscissas X and weights W:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
//
//    For the OPEN rule, the abscissas do not include the end points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  nco_compute_points ( n, x );

  nco_compute_weights ( n, w );

  return;
}
//****************************************************************************80

void nco_compute_points ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_COMPUTE_POINTS: points for a Newton-Cotes Open quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * x_min   
           + ( double ) (   + i + 1 ) * x_max ) 
           / ( double ) ( n     + 1 );
  }

  return;
}
//****************************************************************************80

void nco_compute_weights ( int n, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_COMPUTE_WEIGHTS: weights for a Newton-Cotes Open quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double W[N], the weights.
//
{
  int i;
  double *x;
  double x_max = 1.0;
  double x_min = -1.0;

  x = new double[n];

  nco_compute_points ( n, x );

  nc_compute ( n, x_min, x_max, x, w );

  delete [] x;

  return;
}
//****************************************************************************80

void nco_set ( int n, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_SET sets abscissas and weights for open Newton-Cotes quadrature.
//
//  Discussion:
//
//    The open Newton-Cotes rules use equally spaced abscissas, and
//    hence may be used with equally spaced data.
//
//    The rules are called "open" because they do not include the interval
//    endpoints.
//
//    Most of the rules involve negative weights.  These can produce loss
//    of accuracy due to the subtraction of large, nearly equal quantities.
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 7, and 9.
//
//    Output, double XTAB[N], the abscissas.
//
//    Output, double WEIGHT[N], the weights.
//
{
  double d;
  int i;

  if ( n == 1 )
  {
    weight[0] = 2.0;
  }
  else if ( n == 2 )
  {
    weight[0] = 1.0;
    weight[1] = 1.0;
  }
  else if ( n == 3 )
  {
    d = 3.0;

    weight[0] =   4.0 / d;
    weight[1] = - 2.0 / d;
    weight[2] =   4.0 / d;
  }
  else if ( n == 4 )
  {
    d = 12.0;

    weight[0] = 11.0 / d;
    weight[1] =  1.0 / d;
    weight[2] =  1.0 / d;
    weight[3] = 11.0 / d;
  }
  else if ( n == 5 )
  {
    d = 10.0;

    weight[0] =   11.0 / d;
    weight[1] = - 14.0 / d;
    weight[2] =   26.0 / d;
    weight[3] = - 14.0 / d;
    weight[4] =   11.0 / d;
  }
  else if ( n == 6 )
  {
    d = 1440.0;

    weight[0] =  1222.0 / d;
    weight[1] = - 906.0 / d;
    weight[2] =  1124.0 / d;
    weight[3] =  1124.0 / d;
    weight[4] = - 906.0 / d;
    weight[5] =  1222.0 / d;
  }
  else if ( n == 7 )
  {
    d = 945.0;

    weight[0] =    920.0 / d;
    weight[1] = - 1908.0 / d;
    weight[2] =   4392.0 / d;
    weight[3] = - 4918.0 / d;
    weight[4] =   4392.0 / d;
    weight[5] = - 1908.0 / d;
    weight[6] =    920.0 / d;
  }
  else if ( n == 9 )
  {
    d = 4536.0;

    weight[0] =    4045.0 / d;
    weight[1] = - 11690.0 / d;
    weight[2] =   33340.0 / d;
    weight[3] = - 55070.0 / d;
    weight[4] =   67822.0 / d;
    weight[5] = - 55070.0 / d;
    weight[6] =   33340.0 / d;
    weight[7] = - 11690.0 / d;
    weight[8] =    4045.0 / d;
  }
  else
  {
    cerr << "\n";
    cerr << "NCO_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1 to 7, and 9.\n";
    exit ( 1 );
  }
//
//  Set the abscissas.
//
  for ( i = 0; i < n; i++ )
  {
    xtab[i] = ( double ) ( 2 * i - n + 1 ) 
            / ( double ) ( n + 1 );
  }

  return;
}
//****************************************************************************80

void ncoh_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCOH_COMPUTE computes a Newton-Cotes "open half" quadrature rule.
//
//  Discussion:
//
//    The input value N is used to define N equal subintervals of [-1,+1].
//    The I-th abscissa is the center of the I-th subinterval.
//
//    The integral:
//
//      Integral ( X_MIN <= X <= X_MAX ) F(X) dx
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  const double x_max = 1.0;
  const double x_min = -1.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( 2 * n - 2 * i - 1 ) * x_min   
           + ( double ) (         2 * i + 1 ) * x_max ) 
           / ( double ) ( 2 * n             );
  }

  nc_compute ( n, x_min, x_max, x, w );

  return;
}
//****************************************************************************80

void ncoh_compute_points ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCOH_COMPUTE_POINTS computes points for a Newton-Cotes "open half" quadrature rule.
//
//  Discussion:
//
//    The input value N is used to define N equal subintervals of [-1,+1].
//    The I-th abscissa is the center of the I-th subinterval.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
{
  int i;
  const double x_max = 1.0;
  const double x_min = -1.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( 2 * n - 2 * i - 1 ) * x_min   
           + ( double ) (         2 * i + 1 ) * x_max ) 
           / ( double ) ( 2 * n             );
  }

  return;
}
//****************************************************************************80

void ncoh_compute_weights ( int n, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCOH_COMPUTE_WEIGHTS computes weights for a Newton-Cotes "open half" quadrature rule.
//
//  Discussion:
//
//    The input value N is used to define N equal subintervals of [-1,+1].
//    The I-th abscissa is the center of the I-th subinterval.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double W[N], the weights.
//
{
  int i;
  double *x;
  const double x_max = 1.0;
  const double x_min = -1.0;

  x = new double[n];

  ncoh_compute_points ( n, x );

  nc_compute ( n, x_min, x_max, x, w );

  delete [] x;

  return;
}
//****************************************************************************80

void ncoh_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCOH_SET sets abscissas and weights for Newton-Cotes "open half" rules.
//
//  Discussion:
//
//    The open Newton-Cotes rules use equally spaced abscissas, and
//    hence may be used with equally spaced data.
//
//    The rules are called "open" because the abscissas do not include
//    the interval endpoints.
//
//    For this uncommon type of open Newton-Cotes rule, the abscissas for
//    rule N are found by dividing the interval into N equal subintervals,
//    and using the midpoint of each subinterval as the abscissa.
//
//    Most of the rules involve negative weights.  These can produce loss
//    of accuracy due to the subtraction of large, nearly equal quantities.
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//    In Mathematica, the Newton-Cotes "open half" weights and abscissas
//    can be computed by the commands:
//
//      Needs["NumericalDifferentialEquationAnalysis`"]
//      NewtonCotesWeights [ n, -1, 1, QuadratureType -> Open ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N <= 17.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double a;
  double b;
  double d;
  int i;

  if ( n == 1 )
  {
    w[0] = 2.0E+00;
  }
  else if ( n == 2 )
  {
    w[0] = 1.0E+00;
    w[1] = 1.0E+00;
  }
  else if ( n == 3 )
  {
    d = 4.0E+00;

    w[0] =   3.0E+00 / d;
    w[1] =   2.0E+00 / d;
    w[2] =   3.0E+00 / d;
  }
  else if ( n == 4 )
  {
    d = 24.0E+00;

    w[0] = 13.0E+00 / d;
    w[1] = 11.0E+00 / d;
    w[2] = 11.0E+00 / d;
    w[3] = 13.0E+00 / d;
  }
  else if ( n == 5 )
  {
    d = 576.0E+00;

    w[0] =  275.0E+00 / d;
    w[1] =  100.0E+00 / d;
    w[2] =  402.0E+00 / d;
    w[3] =  100.0E+00 / d;
    w[4] =  275.0E+00 / d;
  }
  else if ( n == 6 )
  {
    d = 640.0E+00;

    w[0] =   247.0E+00 / d;
    w[1] =   139.0E+00 / d;
    w[2] =   254.0E+00 / d;
    w[3] =   254.0E+00 / d;
    w[4] =   139.0E+00 / d;
    w[5] =   247.0E+00 / d;
  }
  else if ( n == 7 )
  {
    d = 138240.0E+00;

    w[0] =   49490.0E+00 / d;
    w[1] =    1764.0E+00 / d;
    w[2] =  112014.0E+00 / d;
    w[3] =  -50056.0E+00 / d;
    w[4] =  112014.0E+00 / d;
    w[5] =    1764.0E+00 / d;
    w[6] =   49490.0E+00 / d;
  }
  else if ( n == 8 )
  {
    d = 967680.0E+00;

    w[0] =  295627.0E+00 / d;
    w[1] =   71329.0E+00 / d;
    w[2] =  471771.0E+00 / d;
    w[3] =  128953.0E+00 / d;
    w[4] =  128953.0E+00 / d;
    w[5] =  471771.0E+00 / d;
    w[6] =   71329.0E+00 / d;
    w[7] =  295627.0E+00 / d;
  }
  else if ( n == 9 )
  {
    d = 2867200.0E+00;

    w[0] =    832221.0E+00 / d;
    w[1] =   -260808.0E+00 / d;
    w[2] =   2903148.0E+00 / d;
    w[3] =  -3227256.0E+00 / d;
    w[4] =   5239790.0E+00 / d;
    w[5] =  -3227256.0E+00 / d;
    w[6] =   2903148.0E+00 / d;
    w[7] =   -260808.0E+00 / d;
    w[8] =    832221.0E+00 / d;
  }
  else if ( n == 10 )
  {
    d = 18579456.0E+00;

    w[0] =    4751285.0E+00 / d;
    w[1] =    -492755.0E+00 / d;
    w[2] =   12269956.0E+00 / d;
    w[3] =   -6274220.0E+00 / d;
    w[4] =    8325190.0E+00 / d;
    w[5] =    8325190.0E+00 / d;
    w[6] =   -6274220.0E+00 / d;
    w[7] =   12269956.0E+00 / d;
    w[8] =    -492755.0E+00 / d;
    w[9] =   4751285.0E+00 / d;
  }
  else if ( n == 11 )
  {
    w[ 0] = + 0.246271364278193434009406231628;
    w[ 1] = - 0.167027133984260177836566725456;
    w[ 2] = + 1.27129728179339588844797178131;
    w[ 3] = - 2.19004533609595458553791887125;
    w[ 4] = + 3.91748917836991567460317460317;
    w[ 5] = - 4.15597070872258046737213403880;
    w[ 6] = + 3.91748917836991567460317460317;
    w[ 7] = - 2.19004533609595458553791887125;
    w[ 8] = + 1.27129728179339588844797178131;
    w[9] = - 0.167027133984260177836566725456;
    w[10] = + 0.246271364278193434009406231628;
  }
  else if ( n == 12 )
  {
    w[ 0] = + 0.221603210581329477813852813853;
    w[ 1] = - 0.103156166902352205086580086580;
    w[ 2] = + 0.889254983348763866341991341991;
    w[ 3] = - 1.08160728355506797889610389610;
    w[ 4] = + 1.49180546087620062229437229437;
    w[ 5] = - 0.417900204348873782467532467532;
    w[ 6] = - 0.417900204348873782467532467532;
    w[ 7] = + 1.49180546087620062229437229437;
    w[ 8] = - 1.08160728355506797889610389610; 
    w[9] = + 0.889254983348763866341991341991;
    w[10] = - 0.103156166902352205086580086580;
    w[11] = + 0.221603210581329477813852813853;
  }
  else if ( n == 13 )
  {
    w[0] = 0.215232356419153566228270676022;
    w[1] = -0.227154289276070155983970468097;
    w[2] = 1.57154640756958579127322926926;
    w[3] = -3.60188931740556785445074962271;
    w[4] = 7.51615534838963020202557032914;
    w[5] = -10.7785343238762239297023523214;
    w[6] = 12.6092876363589847612200042756;
    w[7] = -10.7785343238762239297023523214;
    w[8] = 7.51615534838963020202557032914;
    w[9] = -3.60188931740556785445074962271;
    w[10] = 1.57154640756958579127322926926;
    w[11] = -0.227154289276070155983970468097;
    w[12] = 0.215232356419153566228270676022;
  }
  else if ( n == 14 )
  {
    w[0] = 0.196600731862944474955289480752;
    w[1] = -0.165179242362168469504443173425;
    w[2] = 1.16085790162743923526801130968;
    w[3] = -2.14582271238684154514413484321;
    w[4] = 3.66584923423684682693019643251;
    w[5] = -3.34045051168652382743365816282;
    w[6] = 1.62814459870830330492873895652;
    w[7] = 1.62814459870830330492873895652;
    w[8] = -3.34045051168652382743365816282;
    w[9] = 3.66584923423684682693019643251;
    w[10] = -2.14582271238684154514413484321;
    w[11] = 1.16085790162743923526801130968;
    w[12] = -0.165179242362168469504443173425;
    w[13] = 0.196600731862944474955289480752;
  }
  else if ( n == 15 )
  {
    w[0] = 0.192053656112251156523443074782;
    w[1] = -0.277042941258250537556131864168;
    w[2] = 1.90509434600895135399617123947;
    w[3] = -5.39701622989083452471078029114;
    w[4] = 13.1085281753546466152623727959;
    w[5] = -23.3466945206436771323681898459;
    w[6] = 33.4478422682091702199516443378;
    w[7] = -37.2655295077845143021970588935;
    w[8] = 33.4478422682091702199516443378;
    w[9] = -23.3466945206436771323681898459;
    w[10] = 13.1085281753546466152623727959;
    w[11] = -5.39701622989083452471078029114;
    w[12] = 1.90509434600895135399617123947;
    w[13] = -0.277042941258250537556131864168;
    w[14] = 0.192053656112251156523443074782;
  }
  else if ( n == 16 )
  {
    w[0] = 0.177408479879589716830780293564;
    w[1] = -0.217359399771056183974705378131;
    w[2] = 1.46740967914595726066296033468;
    w[3] = -3.56820982596198712548876407280;
    w[4] = 7.42429624597447227175662173974;
    w[5] = -10.1614344802943189309930887295;
    w[6] = 9.74825566269696996625640284529;
    w[7] = -3.87036636166962697505020703289;
    w[8] = -3.87036636166962697505020703289;
    w[9] = 9.74825566269696996625640284529;
    w[10] = -10.1614344802943189309930887295;
    w[11] = 7.42429624597447227175662173974;
    w[12] = -3.56820982596198712548876407280;
    w[13] = 1.46740967914595726066296033468;
    w[14] = -0.217359399771056183974705378131;
    w[15] = 0.177408479879589716830780293564;
  }
  else if ( n == 17 )
  {
    w[0] = 0.174021728363203743659784786159;
    w[1] = -0.319844636797863622878597303396;
    w[2] = 2.26685253417620917889510925819;
    w[3] = -7.60565246092744264614795072379;
    w[4] = 21.2205863313196208783745036601;
    w[5] = -44.9326914054546061828308816595;
    w[6] = 76.6598740687724224896458863733;
    w[7] = -104.621704713086021393464459433;
    w[8] = 116.317117107268955109493210084;
    w[9] = -104.621704713086021393464459433;
    w[10] = 76.6598740687724224896458863733;
    w[11] = -44.9326914054546061828308816595;
    w[12] = 21.2205863313196208783745036601;
    w[13] = -7.60565246092744264614795072379;
    w[14] = 2.26685253417620917889510925819;
    w[15] = -0.319844636797863622878597303396;
    w[16] = 0.174021728363203743659784786159;
  }
  else
  {
    cerr << "\n";
    cerr << "NCOH_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1 to 17.\n";
    exit ( 1 );
  }
//
//  Set the abscissas.
//
  a = -1.0E+00;
  b = +1.0E+00;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( 2 * n - 2 * i - 1 ) * a   
           + ( double ) (         2 * i + 1 ) * b ) 
           / ( double ) ( 2 * n                   );
  }

  return;
}
//****************************************************************************80

void patterson_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    PATTERSON_SET sets abscissas and weights for Gauss-Patterson quadrature.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W[I) * F ( X[I) )
//
//    The zeroth rule, of order 1, is the standard Gauss-Legendre rule.
//
//    The first rule, of order 3, is the standard Gauss-Legendre rule.
//
//    The second rule, of order 7, includes the abscissas of the previous
//    rule.
//
//    Rules are available of orders 1, 3, 7, 15, 31, 63, 127, 255, and 511.
//
//    These rules constitute a nested family.  The rules can integrate exactly
//    any polynomial of degree 1, 5, 11, 23, 47, 95, 191, 383 or 767, 
//    respectively.
//
//    The data for N = 511 was supplied by Dirk Laurie, and is derived
//    from a NAG Library function d01arf.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Prem Kythe, Michael Schaeferkotter,
//    Handbook of Computational Methods for Integration,
//    Chapman and Hall, 2004,
//    ISBN: 1-58488-428-2,
//    LC: QA299.3.K98.
//
//    NAG Library Documentation,
//    D01ARF,
//    The Numerical Algorithms Group.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 7, 15, 31, 63, 127, 255 or 511.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[0] =   0.0;

    w[0] = 2.0;
  }
  else if ( n == 3 )
  {
    x[0] = -0.77459666924148337704;
    x[1] =  0.0;
    x[2] =  0.77459666924148337704;

    w[0] = 0.555555555555555555556;
    w[1] = 0.888888888888888888889;
    w[2] = 0.555555555555555555556;
  }
  else if ( n == 7 )
  {
    x[0] = -0.96049126870802028342;
    x[1] = -0.77459666924148337704;
    x[2] = -0.43424374934680255800;
    x[3] =  0.0;
    x[4] =  0.43424374934680255800;
    x[5] =  0.77459666924148337704;
    x[6] =  0.96049126870802028342;

    w[0] = 0.104656226026467265194;
    w[1] = 0.268488089868333440729;
    w[2] = 0.401397414775962222905;
    w[3] = 0.450916538658474142345;
    w[4] = 0.401397414775962222905;
    w[5] = 0.268488089868333440729;
    w[6] = 0.104656226026467265194;
  }
  else if ( n == 15 )
  {
    x[ 0] = -0.99383196321275502221;
    x[ 1] = -0.96049126870802028342;
    x[ 2] = -0.88845923287225699889;
    x[ 3] = -0.77459666924148337704;
    x[ 4] = -0.62110294673722640294;
    x[ 5] = -0.43424374934680255800;
    x[ 6] = -0.22338668642896688163;
    x[ 7] =  0.0;
    x[ 8] =  0.22338668642896688163;
    x[ 9] =  0.43424374934680255800;
    x[10] =  0.62110294673722640294;
    x[11] =  0.77459666924148337704;
    x[12] =  0.88845923287225699889;
    x[13] =  0.96049126870802028342;
    x[14] =  0.99383196321275502221;

    w[ 0] = 0.0170017196299402603390;
    w[ 1] = 0.0516032829970797396969;
    w[ 2] = 0.0929271953151245376859;
    w[ 3] = 0.134415255243784220360;
    w[ 4] = 0.171511909136391380787;
    w[ 5] = 0.200628529376989021034;
    w[ 6] = 0.219156858401587496404;
    w[ 7] = 0.225510499798206687386;
    w[ 8] = 0.219156858401587496404;
    w[ 9] = 0.200628529376989021034;
    w[10] = 0.171511909136391380787;
    w[11] = 0.134415255243784220360;
    w[12] = 0.0929271953151245376859;
    w[13] = 0.0516032829970797396969;
    w[14] = 0.0170017196299402603390;
  }
  else if ( n == 31 )
  {
    x[ 0] = -0.99909812496766759766;
    x[ 1] = -0.99383196321275502221;
    x[ 2] = -0.98153114955374010687;
    x[ 3] = -0.96049126870802028342;
    x[ 4] = -0.92965485742974005667;
    x[ 5] = -0.88845923287225699889;
    x[ 6] = -0.83672593816886873550;
    x[ 7] = -0.77459666924148337704;
    x[ 8] = -0.70249620649152707861;
    x[ 9] = -0.62110294673722640294;
    x[10] = -0.53131974364437562397;
    x[11] = -0.43424374934680255800;
    x[12] = -0.33113539325797683309;
    x[13] = -0.22338668642896688163;
    x[14] = -0.11248894313318662575;
    x[15] =  0.0;
    x[16] =  0.11248894313318662575;
    x[17] =  0.22338668642896688163;
    x[18] =  0.33113539325797683309;
    x[19] =  0.43424374934680255800;
    x[20] =  0.53131974364437562397;
    x[21] =  0.62110294673722640294;
    x[22] =  0.70249620649152707861;
    x[23] =  0.77459666924148337704;
    x[24] =  0.83672593816886873550;
    x[25] =  0.88845923287225699889;
    x[26] =  0.92965485742974005667;
    x[27] =  0.96049126870802028342;
    x[28] =  0.98153114955374010687;
    x[29] =  0.99383196321275502221;
    x[30] =  0.99909812496766759766;

    w[ 0] = 0.00254478079156187441540;
    w[ 1] = 0.00843456573932110624631;
    w[ 2] = 0.0164460498543878109338;
    w[ 3] = 0.0258075980961766535646;
    w[ 4] = 0.0359571033071293220968;
    w[ 5] = 0.0464628932617579865414;
    w[ 6] = 0.0569795094941233574122;
    w[ 7] = 0.0672077542959907035404;
    w[ 8] = 0.0768796204990035310427;
    w[ 9] = 0.0857559200499903511542;
    w[10] = 0.0936271099812644736167;
    w[11] = 0.100314278611795578771;
    w[12] = 0.105669893580234809744;
    w[13] = 0.109578421055924638237;
    w[14] = 0.111956873020953456880;
    w[15] = 0.112755256720768691607;
    w[16] = 0.111956873020953456880;
    w[17] = 0.109578421055924638237;
    w[18] = 0.105669893580234809744;
    w[19] = 0.100314278611795578771;
    w[20] = 0.0936271099812644736167;
    w[21] = 0.0857559200499903511542;
    w[22] = 0.0768796204990035310427;
    w[23] = 0.0672077542959907035404;
    w[24] = 0.0569795094941233574122;
    w[25] = 0.0464628932617579865414;
    w[26] = 0.0359571033071293220968;
    w[27] = 0.0258075980961766535646;
    w[28] = 0.0164460498543878109338;
    w[29] = 0.00843456573932110624631;
    w[30] = 0.00254478079156187441540;
  }
  else if ( n == 63 )
  {
    x[ 0] = -0.99987288812035761194;
    x[ 1] = -0.99909812496766759766;
    x[ 2] = -0.99720625937222195908;
    x[ 3] = -0.99383196321275502221;
    x[ 4] = -0.98868475754742947994;
    x[ 5] = -0.98153114955374010687;
    x[ 6] = -0.97218287474858179658;
    x[ 7] = -0.96049126870802028342;
    x[ 8] = -0.94634285837340290515;
    x[ 9] = -0.92965485742974005667;
    x[10] = -0.91037115695700429250;
    x[11] = -0.88845923287225699889;
    x[12] = -0.86390793819369047715;
    x[13] = -0.83672593816886873550;
    x[14] = -0.80694053195021761186;
    x[15] = -0.77459666924148337704;
    x[16] = -0.73975604435269475868;
    x[17] = -0.70249620649152707861;
    x[18] = -0.66290966002478059546;
    x[19] = -0.62110294673722640294;
    x[20] = -0.57719571005204581484;
    x[21] = -0.53131974364437562397;
    x[22] = -0.48361802694584102756;
    x[23] = -0.43424374934680255800;
    x[24] = -0.38335932419873034692;
    x[25] = -0.33113539325797683309;
    x[26] = -0.27774982202182431507;
    x[27] = -0.22338668642896688163;
    x[28] = -0.16823525155220746498;
    x[29] = -0.11248894313318662575;
    x[30] = -0.056344313046592789972;
    x[31] =  0.0;
    x[32] =  0.056344313046592789972;
    x[33] =  0.11248894313318662575;
    x[34] =  0.16823525155220746498;
    x[35] =  0.22338668642896688163;
    x[36] =  0.27774982202182431507;
    x[37] =  0.33113539325797683309;
    x[38] =  0.38335932419873034692;
    x[39] =  0.43424374934680255800;
    x[40] =  0.48361802694584102756;
    x[41] =  0.53131974364437562397;
    x[42] =  0.57719571005204581484;
    x[43] =  0.62110294673722640294;
    x[44] =  0.66290966002478059546;
    x[45] =  0.70249620649152707861;
    x[46] =  0.73975604435269475868;
    x[47] =  0.77459666924148337704;
    x[48] =  0.80694053195021761186;
    x[49] =  0.83672593816886873550;
    x[50] =  0.86390793819369047715;
    x[51] =  0.88845923287225699889;
    x[52] =  0.91037115695700429250;
    x[53] =  0.92965485742974005667;
    x[54] =  0.94634285837340290515;
    x[55] =  0.96049126870802028342;
    x[56] =  0.97218287474858179658;
    x[57] =  0.98153114955374010687;
    x[58] =  0.98868475754742947994;
    x[59] =  0.99383196321275502221;
    x[60] =  0.99720625937222195908;
    x[61] =  0.99909812496766759766;
    x[62] =  0.99987288812035761194;

    w[ 0] = 0.000363221481845530659694;
    w[ 1] = 0.00126515655623006801137;
    w[ 2] = 0.00257904979468568827243;
    w[ 3] = 0.00421763044155885483908;
    w[ 4] = 0.00611550682211724633968;
    w[ 5] = 0.00822300795723592966926;
    w[ 6] = 0.0104982469096213218983;
    w[ 7] = 0.0129038001003512656260;
    w[ 8] = 0.0154067504665594978021;
    w[ 9] = 0.0179785515681282703329;
    w[10] = 0.0205942339159127111492;
    w[11] = 0.0232314466399102694433;
    w[12] = 0.0258696793272147469108;
    w[13] = 0.0284897547458335486125;
    w[14] = 0.0310735511116879648799;
    w[15] = 0.0336038771482077305417;
    w[16] = 0.0360644327807825726401;
    w[17] = 0.0384398102494555320386;
    w[18] = 0.0407155101169443189339;
    w[19] = 0.0428779600250077344929;
    w[20] = 0.0449145316536321974143;
    w[21] = 0.0468135549906280124026;
    w[22] = 0.0485643304066731987159;
    w[23] = 0.0501571393058995374137;
    w[24] = 0.0515832539520484587768;
    w[25] = 0.0528349467901165198621;
    w[26] = 0.0539054993352660639269;
    w[27] = 0.0547892105279628650322;
    w[28] = 0.0554814043565593639878;
    w[29] = 0.0559784365104763194076;
    w[30] = 0.0562776998312543012726;
    w[31] = 0.0563776283603847173877;
    w[32] = 0.0562776998312543012726;
    w[33] = 0.0559784365104763194076;
    w[34] = 0.0554814043565593639878;
    w[35] = 0.0547892105279628650322;
    w[36] = 0.0539054993352660639269;
    w[37] = 0.0528349467901165198621;
    w[38] = 0.0515832539520484587768;
    w[39] = 0.0501571393058995374137;
    w[40] = 0.0485643304066731987159;
    w[41] = 0.0468135549906280124026;
    w[42] = 0.0449145316536321974143;
    w[43] = 0.0428779600250077344929;
    w[44] = 0.0407155101169443189339;
    w[45] = 0.0384398102494555320386;
    w[46] = 0.0360644327807825726401;
    w[47] = 0.0336038771482077305417;
    w[48] = 0.0310735511116879648799;
    w[49] = 0.0284897547458335486125;
    w[50] = 0.0258696793272147469108;
    w[51] = 0.0232314466399102694433;
    w[52] = 0.0205942339159127111492;
    w[53] = 0.0179785515681282703329;
    w[54] = 0.0154067504665594978021;
    w[55] = 0.0129038001003512656260;
    w[56] = 0.0104982469096213218983;
    w[57] = 0.00822300795723592966926;
    w[58] = 0.00611550682211724633968;
    w[59] = 0.00421763044155885483908;
    w[60] = 0.00257904979468568827243;
    w[61] = 0.00126515655623006801137;
    w[62] = 0.000363221481845530659694;
  }
  else if ( n == 127 )
  {
    x[  0] = -0.99998243035489159858;
    x[  1] = -0.99987288812035761194;
    x[  2] = -0.99959879967191068325;
    x[  3] = -0.99909812496766759766;
    x[  4] = -0.99831663531840739253;
    x[  5] = -0.99720625937222195908;
    x[  6] = -0.99572410469840718851;
    x[  7] = -0.99383196321275502221;
    x[  8] = -0.99149572117810613240;
    x[  9] = -0.98868475754742947994;
    x[ 10] = -0.98537149959852037111;
    x[ 11] = -0.98153114955374010687;
    x[ 12] = -0.97714151463970571416;
    x[ 13] = -0.97218287474858179658;
    x[ 14] = -0.96663785155841656709;
    x[ 15] = -0.96049126870802028342;
    x[ 16] = -0.95373000642576113641;
    x[ 17] = -0.94634285837340290515;
    x[ 18] = -0.93832039777959288365;
    x[ 19] = -0.92965485742974005667;
    x[ 20] = -0.92034002547001242073;
    x[ 21] = -0.91037115695700429250;
    x[ 22] = -0.89974489977694003664;
    x[ 23] = -0.88845923287225699889;
    x[ 24] = -0.87651341448470526974;
    x[ 25] = -0.86390793819369047715;
    x[ 26] = -0.85064449476835027976;
    x[ 27] = -0.83672593816886873550;
    x[ 28] = -0.82215625436498040737;
    x[ 29] = -0.80694053195021761186;
    x[ 30] = -0.79108493379984836143;
    x[ 31] = -0.77459666924148337704;
    x[ 32] = -0.75748396638051363793;
    x[ 33] = -0.73975604435269475868;
    x[ 34] = -0.72142308537009891548;
    x[ 35] = -0.70249620649152707861;
    x[ 36] = -0.68298743109107922809;
    x[ 37] = -0.66290966002478059546;
    x[ 38] = -0.64227664250975951377;
    x[ 39] = -0.62110294673722640294;
    x[ 40] = -0.59940393024224289297;
    x[ 41] = -0.57719571005204581484;
    x[ 42] = -0.55449513263193254887;
    x[ 43] = -0.53131974364437562397;
    x[ 44] = -0.50768775753371660215;
    x[ 45] = -0.48361802694584102756;
    x[ 46] = -0.45913001198983233287;
    x[ 47] = -0.43424374934680255800;
    x[ 48] = -0.40897982122988867241;
    x[ 49] = -0.38335932419873034692;
    x[ 50] = -0.35740383783153215238;
    x[ 51] = -0.33113539325797683309;
    x[ 52] = -0.30457644155671404334;
    x[ 53] = -0.27774982202182431507;
    x[ 54] = -0.25067873030348317661;
    x[ 55] = -0.22338668642896688163;
    x[ 56] = -0.19589750271110015392;
    x[ 57] = -0.16823525155220746498;
    x[ 58] = -0.14042423315256017459;
    x[ 59] = -0.11248894313318662575;
    x[ 60] = -0.084454040083710883710;
    x[ 61] = -0.056344313046592789972;
    x[ 62] = -0.028184648949745694339;
    x[ 63] =  0.0;
    x[ 64] =  0.028184648949745694339;
    x[ 65] =  0.056344313046592789972;
    x[ 66] =  0.084454040083710883710;
    x[ 67] =  0.11248894313318662575;
    x[ 68] =  0.14042423315256017459;
    x[ 69] =  0.16823525155220746498;
    x[ 70] =  0.19589750271110015392;
    x[ 71] =  0.22338668642896688163;
    x[ 72] =  0.25067873030348317661;
    x[ 73] =  0.27774982202182431507;
    x[ 74] =  0.30457644155671404334;
    x[ 75] =  0.33113539325797683309;
    x[ 76] =  0.35740383783153215238;
    x[ 77] =  0.38335932419873034692;
    x[ 78] =  0.40897982122988867241;
    x[ 79] =  0.43424374934680255800;
    x[ 80] =  0.45913001198983233287;
    x[ 81] =  0.48361802694584102756;
    x[ 82] =  0.50768775753371660215;
    x[ 83] =  0.53131974364437562397;
    x[ 84] =  0.55449513263193254887;
    x[ 85] =  0.57719571005204581484;
    x[ 86] =  0.59940393024224289297;
    x[ 87] =  0.62110294673722640294;
    x[ 88] =  0.64227664250975951377;
    x[ 89] =  0.66290966002478059546;
    x[ 90] =  0.68298743109107922809;
    x[ 91] =  0.70249620649152707861;
    x[ 92] =  0.72142308537009891548;
    x[ 93] =  0.73975604435269475868;
    x[ 94] =  0.75748396638051363793;
    x[ 95] =  0.77459666924148337704;
    x[ 96] =  0.79108493379984836143;
    x[ 97] =  0.80694053195021761186;
    x[ 98] =  0.82215625436498040737;
    x[ 99] =  0.83672593816886873550;
    x[100] =  0.85064449476835027976;
    x[101] =  0.86390793819369047715;
    x[102] =  0.87651341448470526974;
    x[103] =  0.88845923287225699889;
    x[104] =  0.89974489977694003664;
    x[105] =  0.91037115695700429250;
    x[106] =  0.92034002547001242073;
    x[107] =  0.92965485742974005667;
    x[108] =  0.93832039777959288365;
    x[109] =  0.94634285837340290515;
    x[110] =  0.95373000642576113641;
    x[111] =  0.96049126870802028342;
    x[112] =  0.96663785155841656709;
    x[113] =  0.97218287474858179658;
    x[114] =  0.97714151463970571416;
    x[115] =  0.98153114955374010687;
    x[116] =  0.98537149959852037111;
    x[117] =  0.98868475754742947994;
    x[118] =  0.99149572117810613240;
    x[119] =  0.99383196321275502221;
    x[120] =  0.99572410469840718851;
    x[121] =  0.99720625937222195908;
    x[122] =  0.99831663531840739253;
    x[123] =  0.99909812496766759766;
    x[124] =  0.99959879967191068325;
    x[125] =  0.99987288812035761194;
    x[126] =  0.99998243035489159858;

    w[  0] = 0.0000505360952078625176247;
    w[  1] = 0.000180739564445388357820;
    w[  2] = 0.000377746646326984660274;
    w[  3] = 0.000632607319362633544219;
    w[  4] = 0.000938369848542381500794;
    w[  5] = 0.00128952408261041739210;
    w[  6] = 0.00168114286542146990631;
    w[  7] = 0.00210881524572663287933;
    w[  8] = 0.00256876494379402037313;
    w[  9] = 0.00305775341017553113613;
    w[ 10] = 0.00357289278351729964938;
    w[ 11] = 0.00411150397865469304717;
    w[ 12] = 0.00467105037211432174741;
    w[ 13] = 0.00524912345480885912513;
    w[ 14] = 0.00584344987583563950756;
    w[ 15] = 0.00645190005017573692280;
    w[ 16] = 0.00707248999543355546805;
    w[ 17] = 0.00770337523327974184817;
    w[ 18] = 0.00834283875396815770558;
    w[ 19] = 0.00898927578406413572328;
    w[ 20] = 0.00964117772970253669530;
    w[ 21] = 0.0102971169579563555237;
    w[ 22] = 0.0109557333878379016480;
    w[ 23] = 0.0116157233199551347270;
    w[ 24] = 0.0122758305600827700870;
    w[ 25] = 0.0129348396636073734547;
    w[ 26] = 0.0135915710097655467896;
    w[ 27] = 0.0142448773729167743063;
    w[ 28] = 0.0148936416648151820348;
    w[ 29] = 0.0155367755558439824399;
    w[ 30] = 0.0161732187295777199419;
    w[ 31] = 0.0168019385741038652709;
    w[ 32] = 0.0174219301594641737472;
    w[ 33] = 0.0180322163903912863201;
    w[ 34] = 0.0186318482561387901863;
    w[ 35] = 0.0192199051247277660193;
    w[ 36] = 0.0197954950480974994880;
    w[ 37] = 0.0203577550584721594669;
    w[ 38] = 0.0209058514458120238522;
    w[ 39] = 0.0214389800125038672465;
    w[ 40] = 0.0219563663053178249393;
    w[ 41] = 0.0224572658268160987071;
    w[ 42] = 0.0229409642293877487608;
    w[ 43] = 0.0234067774953140062013;
    w[ 44] = 0.0238540521060385400804;
    w[ 45] = 0.0242821652033365993580;
    w[ 46] = 0.0246905247444876769091;
    w[ 47] = 0.0250785696529497687068;
    w[ 48] = 0.0254457699654647658126;
    w[ 49] = 0.0257916269760242293884;
    w[ 50] = 0.0261156733767060976805;
    w[ 51] = 0.0264174733950582599310;
    w[ 52] = 0.0266966229274503599062;
    w[ 53] = 0.0269527496676330319634;
    w[ 54] = 0.0271855132296247918192;
    w[ 55] = 0.0273946052639814325161;
    w[ 56] = 0.0275797495664818730349;
    w[ 57] = 0.0277407021782796819939;
    w[ 58] = 0.0278772514766137016085;
    w[ 59] = 0.0279892182552381597038;
    w[ 60] = 0.0280764557938172466068;
    w[ 61] = 0.0281388499156271506363;
    w[ 62] = 0.0281763190330166021307;
    w[ 63] = 0.0281888141801923586938;
    w[ 64] = 0.0281763190330166021307;
    w[ 65] = 0.0281388499156271506363;
    w[ 66] = 0.0280764557938172466068;
    w[ 67] = 0.0279892182552381597038;
    w[ 68] = 0.0278772514766137016085;
    w[ 69] = 0.0277407021782796819939;
    w[ 70] = 0.0275797495664818730349;
    w[ 71] = 0.0273946052639814325161;
    w[ 72] = 0.0271855132296247918192;
    w[ 73] = 0.0269527496676330319634;
    w[ 74] = 0.0266966229274503599062;
    w[ 75] = 0.0264174733950582599310;
    w[ 76] = 0.0261156733767060976805;
    w[ 77] = 0.0257916269760242293884;
    w[ 78] = 0.0254457699654647658126;
    w[ 79] = 0.0250785696529497687068;
    w[ 80] = 0.0246905247444876769091;
    w[ 81] = 0.0242821652033365993580;
    w[ 82] = 0.0238540521060385400804;
    w[ 83] = 0.0234067774953140062013;
    w[ 84] = 0.0229409642293877487608;
    w[ 85] = 0.0224572658268160987071;
    w[ 86] = 0.0219563663053178249393;
    w[ 87] = 0.0214389800125038672465;
    w[ 88] = 0.0209058514458120238522;
    w[ 89] = 0.0203577550584721594669;
    w[ 90] = 0.0197954950480974994880;
    w[ 91] = 0.0192199051247277660193;
    w[ 92] = 0.0186318482561387901863;
    w[ 93] = 0.0180322163903912863201;
    w[ 94] = 0.0174219301594641737472;
    w[ 95] = 0.0168019385741038652709;
    w[ 96] = 0.0161732187295777199419;
    w[ 97] = 0.0155367755558439824399;
    w[ 98] = 0.0148936416648151820348;
    w[ 99] = 0.0142448773729167743063;
    w[100] = 0.0135915710097655467896;
    w[101] = 0.0129348396636073734547;
    w[102] = 0.0122758305600827700870;
    w[103] = 0.0116157233199551347270;
    w[104] = 0.0109557333878379016480;
    w[105] = 0.0102971169579563555237;
    w[106] = 0.00964117772970253669530;
    w[107] = 0.00898927578406413572328;
    w[108] = 0.00834283875396815770558;
    w[109] = 0.00770337523327974184817;
    w[110] = 0.00707248999543355546805;
    w[111] = 0.00645190005017573692280;
    w[112] = 0.00584344987583563950756;
    w[113] = 0.00524912345480885912513;
    w[114] = 0.00467105037211432174741;
    w[115] = 0.00411150397865469304717;
    w[116] = 0.00357289278351729964938;
    w[117] = 0.00305775341017553113613;
    w[118] = 0.00256876494379402037313;
    w[119] = 0.00210881524572663287933;
    w[120] = 0.00168114286542146990631;
    w[121] = 0.00128952408261041739210;
    w[122] = 0.000938369848542381500794;
    w[123] = 0.000632607319362633544219;
    w[124] = 0.000377746646326984660274;
    w[125] = 0.000180739564445388357820;
    w[126] = 0.0000505360952078625176247;
  }
  else if ( n == 255 )
  {
    x[  0] = -0.99999759637974846462;
    x[  1] = -0.99998243035489159858;
    x[  2] = -0.99994399620705437576;
    x[  3] = -0.99987288812035761194;
    x[  4] = -0.99976049092443204733;
    x[  5] = -0.99959879967191068325;
    x[  6] = -0.99938033802502358193;
    x[  7] = -0.99909812496766759766;
    x[  8] = -0.99874561446809511470;
    x[  9] = -0.99831663531840739253;
    x[ 10] = -0.99780535449595727456;
    x[ 11] = -0.99720625937222195908;
    x[ 12] = -0.99651414591489027385;
    x[ 13] = -0.99572410469840718851;
    x[ 14] = -0.99483150280062100052;
    x[ 15] = -0.99383196321275502221;
    x[ 16] = -0.99272134428278861533;
    x[ 17] = -0.99149572117810613240;
    x[ 18] = -0.99015137040077015918;
    x[ 19] = -0.98868475754742947994;
    x[ 20] = -0.98709252795403406719;
    x[ 21] = -0.98537149959852037111;
    x[ 22] = -0.98351865757863272876;
    x[ 23] = -0.98153114955374010687;
    x[ 24] = -0.97940628167086268381;
    x[ 25] = -0.97714151463970571416;
    x[ 26] = -0.97473445975240266776;
    x[ 27] = -0.97218287474858179658;
    x[ 28] = -0.96948465950245923177;
    x[ 29] = -0.96663785155841656709;
    x[ 30] = -0.96364062156981213252;
    x[ 31] = -0.96049126870802028342;
    x[ 32] = -0.95718821610986096274;
    x[ 33] = -0.95373000642576113641;
    x[ 34] = -0.95011529752129487656;
    x[ 35] = -0.94634285837340290515;
    x[ 36] = -0.94241156519108305981;
    x[ 37] = -0.93832039777959288365;
    x[ 38] = -0.93406843615772578800;
    x[ 39] = -0.92965485742974005667;
    x[ 40] = -0.92507893290707565236;
    x[ 41] = -0.92034002547001242073;
    x[ 42] = -0.91543758715576504064;
    x[ 43] = -0.91037115695700429250;
    x[ 44] = -0.90514035881326159519;
    x[ 45] = -0.89974489977694003664;
    x[ 46] = -0.89418456833555902286;
    x[ 47] = -0.88845923287225699889;
    x[ 48] = -0.88256884024734190684;
    x[ 49] = -0.87651341448470526974;
    x[ 50] = -0.87029305554811390585;
    x[ 51] = -0.86390793819369047715;
    x[ 52] = -0.85735831088623215653;
    x[ 53] = -0.85064449476835027976;
    x[ 54] = -0.84376688267270860104;
    x[ 55] = -0.83672593816886873550;
    x[ 56] = -0.82952219463740140018;
    x[ 57] = -0.82215625436498040737;
    x[ 58] = -0.81462878765513741344;
    x[ 59] = -0.80694053195021761186;
    x[ 60] = -0.79909229096084140180;
    x[ 61] = -0.79108493379984836143;
    x[ 62] = -0.78291939411828301639;
    x[ 63] = -0.77459666924148337704;
    x[ 64] = -0.76611781930376009072;
    x[ 65] = -0.75748396638051363793;
    x[ 66] = -0.74869629361693660282;
    x[ 67] = -0.73975604435269475868;
    x[ 68] = -0.73066452124218126133;
    x[ 69] = -0.72142308537009891548;
    x[ 70] = -0.71203315536225203459;
    x[ 71] = -0.70249620649152707861;
    x[ 72] = -0.69281376977911470289;
    x[ 73] = -0.68298743109107922809;
    x[ 74] = -0.67301883023041847920;
    x[ 75] = -0.66290966002478059546;
    x[ 76] = -0.65266166541001749610;
    x[ 77] = -0.64227664250975951377;
    x[ 78] = -0.63175643771119423041;
    x[ 79] = -0.62110294673722640294;
    x[ 80] = -0.61031811371518640016;
    x[ 81] = -0.59940393024224289297;
    x[ 82] = -0.58836243444766254143;
    x[ 83] = -0.57719571005204581484;
    x[ 84] = -0.56590588542365442262;
    x[ 85] = -0.55449513263193254887;
    x[ 86] = -0.54296566649831149049;
    x[ 87] = -0.53131974364437562397;
    x[ 88] = -0.51955966153745702199;
    x[ 89] = -0.50768775753371660215;
    x[ 90] = -0.49570640791876146017;
    x[ 91] = -0.48361802694584102756;
    x[ 92] = -0.47142506587165887693;
    x[ 93] = -0.45913001198983233287;
    x[ 94] = -0.44673538766202847374;
    x[ 95] = -0.43424374934680255800;
    x[ 96] = -0.42165768662616330006;
    x[ 97] = -0.40897982122988867241;
    x[ 98] = -0.39621280605761593918;
    x[ 99] = -0.38335932419873034692;
    x[100] = -0.37042208795007823014;
    x[101] = -0.35740383783153215238;
    x[102] = -0.34430734159943802278;
    x[103] = -0.33113539325797683309;
    x[104] = -0.31789081206847668318;
    x[105] = -0.30457644155671404334;
    x[106] = -0.29119514851824668196;
    x[107] = -0.27774982202182431507;
    x[108] = -0.26424337241092676194;
    x[109] = -0.25067873030348317661;
    x[110] = -0.23705884558982972721;
    x[111] = -0.22338668642896688163;
    x[112] = -0.20966523824318119477;
    x[113] = -0.19589750271110015392;
    x[114] = -0.18208649675925219825;
    x[115] = -0.16823525155220746498;
    x[116] = -0.15434681148137810869;
    x[117] = -0.14042423315256017459;
    x[118] = -0.12647058437230196685;
    x[119] = -0.11248894313318662575;
    x[120] = -0.098482396598119202090;
    x[121] = -0.084454040083710883710;
    x[122] = -0.070406976042855179063;
    x[123] = -0.056344313046592789972;
    x[124] = -0.042269164765363603212;
    x[125] = -0.028184648949745694339;
    x[126] = -0.014093886410782462614;
    x[127] =  0.0;
    x[128] =  0.014093886410782462614;
    x[129] =  0.028184648949745694339;
    x[130] =  0.042269164765363603212;
    x[131] =  0.056344313046592789972;
    x[132] =  0.070406976042855179063;
    x[133] =  0.084454040083710883710;
    x[134] =  0.098482396598119202090;
    x[135] =  0.11248894313318662575;
    x[136] =  0.12647058437230196685;
    x[137] =  0.14042423315256017459;
    x[138] =  0.15434681148137810869;
    x[139] =  0.16823525155220746498;
    x[140] =  0.18208649675925219825;
    x[141] =  0.19589750271110015392;
    x[142] =  0.20966523824318119477;
    x[143] =  0.22338668642896688163;
    x[144] =  0.23705884558982972721;
    x[145] =  0.25067873030348317661;
    x[146] =  0.26424337241092676194;
    x[147] =  0.27774982202182431507;
    x[148] =  0.29119514851824668196;
    x[149] =  0.30457644155671404334;
    x[150] =  0.31789081206847668318;
    x[151] =  0.33113539325797683309;
    x[152] =  0.34430734159943802278;
    x[153] =  0.35740383783153215238;
    x[154] =  0.37042208795007823014;
    x[155] =  0.38335932419873034692;
    x[156] =  0.39621280605761593918;
    x[157] =  0.40897982122988867241;
    x[158] =  0.42165768662616330006;
    x[159] =  0.43424374934680255800;
    x[160] =  0.44673538766202847374;
    x[161] =  0.45913001198983233287;
    x[162] =  0.47142506587165887693;
    x[163] =  0.48361802694584102756;
    x[164] =  0.49570640791876146017;
    x[165] =  0.50768775753371660215;
    x[166] =  0.51955966153745702199;
    x[167] =  0.53131974364437562397;
    x[168] =  0.54296566649831149049;
    x[169] =  0.55449513263193254887;
    x[170] =  0.56590588542365442262;
    x[171] =  0.57719571005204581484;
    x[172] =  0.58836243444766254143;
    x[173] =  0.59940393024224289297;
    x[174] =  0.61031811371518640016;
    x[175] =  0.62110294673722640294;
    x[176] =  0.63175643771119423041;
    x[177] =  0.64227664250975951377;
    x[178] =  0.65266166541001749610;
    x[179] =  0.66290966002478059546;
    x[180] =  0.67301883023041847920;
    x[181] =  0.68298743109107922809;
    x[182] =  0.69281376977911470289;
    x[183] =  0.70249620649152707861;
    x[184] =  0.71203315536225203459;
    x[185] =  0.72142308537009891548;
    x[186] =  0.73066452124218126133;
    x[187] =  0.73975604435269475868;
    x[188] =  0.74869629361693660282;
    x[189] =  0.75748396638051363793;
    x[190] =  0.76611781930376009072;
    x[191] =  0.77459666924148337704;
    x[192] =  0.78291939411828301639;
    x[193] =  0.79108493379984836143;
    x[194] =  0.79909229096084140180;
    x[195] =  0.80694053195021761186;
    x[196] =  0.81462878765513741344;
    x[197] =  0.82215625436498040737;
    x[198] =  0.82952219463740140018;
    x[199] =  0.83672593816886873550;
    x[200] =  0.84376688267270860104;
    x[201] =  0.85064449476835027976;
    x[202] =  0.85735831088623215653;
    x[203] =  0.86390793819369047715;
    x[204] =  0.87029305554811390585;
    x[205] =  0.87651341448470526974;
    x[206] =  0.88256884024734190684;
    x[207] =  0.88845923287225699889;
    x[208] =  0.89418456833555902286;
    x[209] =  0.89974489977694003664;
    x[210] =  0.90514035881326159519;
    x[211] =  0.91037115695700429250;
    x[212] =  0.91543758715576504064;
    x[213] =  0.92034002547001242073;
    x[214] =  0.92507893290707565236;
    x[215] =  0.92965485742974005667;
    x[216] =  0.93406843615772578800;
    x[217] =  0.93832039777959288365;
    x[218] =  0.94241156519108305981;
    x[219] =  0.94634285837340290515;
    x[220] =  0.95011529752129487656;
    x[221] =  0.95373000642576113641;
    x[222] =  0.95718821610986096274;
    x[223] =  0.96049126870802028342;
    x[224] =  0.96364062156981213252;
    x[225] =  0.96663785155841656709;
    x[226] =  0.96948465950245923177;
    x[227] =  0.97218287474858179658;
    x[228] =  0.97473445975240266776;
    x[229] =  0.97714151463970571416;
    x[230] =  0.97940628167086268381;
    x[231] =  0.98153114955374010687;
    x[232] =  0.98351865757863272876;
    x[233] =  0.98537149959852037111;
    x[234] =  0.98709252795403406719;
    x[235] =  0.98868475754742947994;
    x[236] =  0.99015137040077015918;
    x[237] =  0.99149572117810613240;
    x[238] =  0.99272134428278861533;
    x[239] =  0.99383196321275502221;
    x[240] =  0.99483150280062100052;
    x[241] =  0.99572410469840718851;
    x[242] =  0.99651414591489027385;
    x[243] =  0.99720625937222195908;
    x[244] =  0.99780535449595727456;
    x[245] =  0.99831663531840739253;
    x[246] =  0.99874561446809511470;
    x[247] =  0.99909812496766759766;
    x[248] =  0.99938033802502358193;
    x[249] =  0.99959879967191068325;
    x[250] =  0.99976049092443204733;
    x[251] =  0.99987288812035761194;
    x[252] =  0.99994399620705437576;
    x[253] =  0.99998243035489159858;
    x[254] =  0.99999759637974846462;

    w[  0] = 0.69379364324108267170E-05;
    w[  1] = 0.25157870384280661489E-04;
    w[  2] = 0.53275293669780613125E-04;
    w[  3] = 0.90372734658751149261E-04;
    w[  4] = 0.13575491094922871973E-03;
    w[  5] = 0.18887326450650491366E-03;
    w[  6] = 0.24921240048299729402E-03;
    w[  7] = 0.31630366082226447689E-03;
    w[  8] = 0.38974528447328229322E-03;
    w[  9] = 0.46918492424785040975E-03;
    w[ 10] = 0.55429531493037471492E-03;
    w[ 11] = 0.64476204130572477933E-03;
    w[ 12] = 0.74028280424450333046E-03;
    w[ 13] = 0.84057143271072246365E-03;
    w[ 14] = 0.94536151685852538246E-03;
    w[ 15] = 0.10544076228633167722E-02;
    w[ 16] = 0.11674841174299594077E-02;
    w[ 17] = 0.12843824718970101768E-02;
    w[ 18] = 0.14049079956551446427E-02;
    w[ 19] = 0.15288767050877655684E-02;
    w[ 20] = 0.16561127281544526052E-02;
    w[ 21] = 0.17864463917586498247E-02;
    w[ 22] = 0.19197129710138724125E-02;
    w[ 23] = 0.20557519893273465236E-02;
    w[ 24] = 0.21944069253638388388E-02;
    w[ 25] = 0.23355251860571608737E-02;
    w[ 26] = 0.24789582266575679307E-02;
    w[ 27] = 0.26245617274044295626E-02;
    w[ 28] = 0.27721957645934509940E-02;
    w[ 29] = 0.29217249379178197538E-02;
    w[ 30] = 0.30730184347025783234E-02;
    w[ 31] = 0.32259500250878684614E-02;
    w[ 32] = 0.33803979910869203823E-02;
    w[ 33] = 0.35362449977167777340E-02;
    w[ 34] = 0.36933779170256508183E-02;
    w[ 35] = 0.38516876166398709241E-02;
    w[ 36] = 0.40110687240750233989E-02;
    w[ 37] = 0.41714193769840788528E-02;
    w[ 38] = 0.43326409680929828545E-02;
    w[ 39] = 0.44946378920320678616E-02;
    w[ 40] = 0.46573172997568547773E-02;
    w[ 41] = 0.48205888648512683476E-02;
    w[ 42] = 0.49843645647655386012E-02;
    w[ 43] = 0.51485584789781777618E-02;
    w[ 44] = 0.53130866051870565663E-02;
    w[ 45] = 0.54778666939189508240E-02;
    w[ 46] = 0.56428181013844441585E-02;
    w[ 47] = 0.58078616599775673635E-02;
    w[ 48] = 0.59729195655081658049E-02;
    w[ 49] = 0.61379152800413850435E-02;
    w[ 50] = 0.63027734490857587172E-02;
    w[ 51] = 0.64674198318036867274E-02;
    w[ 52] = 0.66317812429018878941E-02;
    w[ 53] = 0.67957855048827733948E-02;
    w[ 54] = 0.69593614093904229394E-02;
    w[ 55] = 0.71224386864583871532E-02;
    w[ 56] = 0.72849479805538070639E-02;
    w[ 57] = 0.74468208324075910174E-02;
    w[ 58] = 0.76079896657190565832E-02;
    w[ 59] = 0.77683877779219912200E-02;
    w[ 60] = 0.79279493342948491103E-02;
    w[ 61] = 0.80866093647888599710E-02;
    w[ 62] = 0.82443037630328680306E-02;
    w[ 63] = 0.84009692870519326354E-02;
    w[ 64] = 0.85565435613076896192E-02;
    w[ 65] = 0.87109650797320868736E-02;
    w[ 66] = 0.88641732094824942641E-02;
    w[ 67] = 0.90161081951956431600E-02;
    w[ 68] = 0.91667111635607884067E-02;
    w[ 69] = 0.93159241280693950932E-02;
    w[ 70] = 0.94636899938300652943E-02;
    w[ 71] = 0.96099525623638830097E-02;
    w[ 72] = 0.97546565363174114611E-02;
    w[ 73] = 0.98977475240487497440E-02;
    w[ 74] = 0.10039172044056840798E-01;
    w[ 75] = 0.10178877529236079733E-01;
    w[ 76] = 0.10316812330947621682E-01;
    w[ 77] = 0.10452925722906011926E-01;
    w[ 78] = 0.10587167904885197931E-01;
    w[ 79] = 0.10719490006251933623E-01;
    w[ 80] = 0.10849844089337314099E-01;
    w[ 81] = 0.10978183152658912470E-01;
    w[ 82] = 0.11104461134006926537E-01;
    w[ 83] = 0.11228632913408049354E-01;
    w[ 84] = 0.11350654315980596602E-01;
    w[ 85] = 0.11470482114693874380E-01;
    w[ 86] = 0.11588074033043952568E-01;
    w[ 87] = 0.11703388747657003101E-01;
    w[ 88] = 0.11816385890830235763E-01;
    w[ 89] = 0.11927026053019270040E-01;
    w[ 90] = 0.12035270785279562630E-01;
    w[ 91] = 0.12141082601668299679E-01;
    w[ 92] = 0.12244424981611985899E-01;
    w[ 93] = 0.12345262372243838455E-01;
    w[ 94] = 0.12443560190714035263E-01;
    w[ 95] = 0.12539284826474884353E-01;
    w[ 96] = 0.12632403643542078765E-01;
    w[ 97] = 0.12722884982732382906E-01;
    w[ 98] = 0.12810698163877361967E-01;
    w[ 99] = 0.12895813488012114694E-01;
    w[100] = 0.12978202239537399286E-01;
    w[101] = 0.13057836688353048840E-01;
    w[102] = 0.13134690091960152836E-01;
    w[103] = 0.13208736697529129966E-01;
    w[104] = 0.13279951743930530650E-01;
    w[105] = 0.13348311463725179953E-01;
    w[106] = 0.13413793085110098513E-01;
    w[107] = 0.13476374833816515982E-01;
    w[108] = 0.13536035934956213614E-01;
    w[109] = 0.13592756614812395910E-01;
    w[110] = 0.13646518102571291428E-01;
    w[111] = 0.13697302631990716258E-01;
    w[112] = 0.13745093443001896632E-01;
    w[113] = 0.13789874783240936517E-01;
    w[114] = 0.13831631909506428676E-01;
    w[115] = 0.13870351089139840997E-01;
    w[116] = 0.13906019601325461264E-01;
    w[117] = 0.13938625738306850804E-01;
    w[118] = 0.13968158806516938516E-01;
    w[119] = 0.13994609127619079852E-01;
    w[120] = 0.14017968039456608810E-01;
    w[121] = 0.14038227896908623303E-01;
    w[122] = 0.14055382072649964277E-01;
    w[123] = 0.14069424957813575318E-01;
    w[124] = 0.14080351962553661325E-01;
    w[125] = 0.14088159516508301065E-01;
    w[126] = 0.14092845069160408355E-01;
    w[127] = 0.14094407090096179347E-01;
    w[128] = 0.14092845069160408355E-01;
    w[129] = 0.14088159516508301065E-01;
    w[130] = 0.14080351962553661325E-01;
    w[131] = 0.14069424957813575318E-01;
    w[132] = 0.14055382072649964277E-01;
    w[133] = 0.14038227896908623303E-01;
    w[134] = 0.14017968039456608810E-01;
    w[135] = 0.13994609127619079852E-01;
    w[136] = 0.13968158806516938516E-01;
    w[137] = 0.13938625738306850804E-01;
    w[138] = 0.13906019601325461264E-01;
    w[139] = 0.13870351089139840997E-01;
    w[140] = 0.13831631909506428676E-01;
    w[141] = 0.13789874783240936517E-01;
    w[142] = 0.13745093443001896632E-01;
    w[143] = 0.13697302631990716258E-01;
    w[144] = 0.13646518102571291428E-01;
    w[145] = 0.13592756614812395910E-01;
    w[146] = 0.13536035934956213614E-01;
    w[147] = 0.13476374833816515982E-01;
    w[148] = 0.13413793085110098513E-01;
    w[149] = 0.13348311463725179953E-01;
    w[150] = 0.13279951743930530650E-01;
    w[151] = 0.13208736697529129966E-01;
    w[152] = 0.13134690091960152836E-01;
    w[153] = 0.13057836688353048840E-01;
    w[154] = 0.12978202239537399286E-01;
    w[155] = 0.12895813488012114694E-01;
    w[156] = 0.12810698163877361967E-01;
    w[157] = 0.12722884982732382906E-01;
    w[158] = 0.12632403643542078765E-01;
    w[159] = 0.12539284826474884353E-01;
    w[160] = 0.12443560190714035263E-01;
    w[161] = 0.12345262372243838455E-01;
    w[162] = 0.12244424981611985899E-01;
    w[163] = 0.12141082601668299679E-01;
    w[164] = 0.12035270785279562630E-01;
    w[165] = 0.11927026053019270040E-01;
    w[166] = 0.11816385890830235763E-01;
    w[167] = 0.11703388747657003101E-01;
    w[168] = 0.11588074033043952568E-01;
    w[169] = 0.11470482114693874380E-01;
    w[170] = 0.11350654315980596602E-01;
    w[171] = 0.11228632913408049354E-01;
    w[172] = 0.11104461134006926537E-01;
    w[173] = 0.10978183152658912470E-01;
    w[174] = 0.10849844089337314099E-01;
    w[175] = 0.10719490006251933623E-01;
    w[176] = 0.10587167904885197931E-01;
    w[177] = 0.10452925722906011926E-01;
    w[178] = 0.10316812330947621682E-01;
    w[179] = 0.10178877529236079733E-01;
    w[180] = 0.10039172044056840798E-01;
    w[181] = 0.98977475240487497440E-02;
    w[182] = 0.97546565363174114611E-02;
    w[183] = 0.96099525623638830097E-02;
    w[184] = 0.94636899938300652943E-02;
    w[185] = 0.93159241280693950932E-02;
    w[186] = 0.91667111635607884067E-02;
    w[187] = 0.90161081951956431600E-02;
    w[188] = 0.88641732094824942641E-02;
    w[189] = 0.87109650797320868736E-02;
    w[190] = 0.85565435613076896192E-02;
    w[191] = 0.84009692870519326354E-02;
    w[192] = 0.82443037630328680306E-02;
    w[193] = 0.80866093647888599710E-02;
    w[194] = 0.79279493342948491103E-02;
    w[195] = 0.77683877779219912200E-02;
    w[196] = 0.76079896657190565832E-02;
    w[197] = 0.74468208324075910174E-02;
    w[198] = 0.72849479805538070639E-02;
    w[199] = 0.71224386864583871532E-02;
    w[200] = 0.69593614093904229394E-02;
    w[201] = 0.67957855048827733948E-02;
    w[202] = 0.66317812429018878941E-02;
    w[203] = 0.64674198318036867274E-02;
    w[204] = 0.63027734490857587172E-02;
    w[205] = 0.61379152800413850435E-02;
    w[206] = 0.59729195655081658049E-02;
    w[207] = 0.58078616599775673635E-02;
    w[208] = 0.56428181013844441585E-02;
    w[209] = 0.54778666939189508240E-02;
    w[210] = 0.53130866051870565663E-02;
    w[211] = 0.51485584789781777618E-02;
    w[212] = 0.49843645647655386012E-02;
    w[213] = 0.48205888648512683476E-02;
    w[214] = 0.46573172997568547773E-02;
    w[215] = 0.44946378920320678616E-02;
    w[216] = 0.43326409680929828545E-02;
    w[217] = 0.41714193769840788528E-02;
    w[218] = 0.40110687240750233989E-02;
    w[219] = 0.38516876166398709241E-02;
    w[220] = 0.36933779170256508183E-02;
    w[221] = 0.35362449977167777340E-02;
    w[222] = 0.33803979910869203823E-02;
    w[223] = 0.32259500250878684614E-02;
    w[224] = 0.30730184347025783234E-02;
    w[225] = 0.29217249379178197538E-02;
    w[226] = 0.27721957645934509940E-02;
    w[227] = 0.26245617274044295626E-02;
    w[228] = 0.24789582266575679307E-02;
    w[229] = 0.23355251860571608737E-02;
    w[230] = 0.21944069253638388388E-02;
    w[231] = 0.20557519893273465236E-02;
    w[232] = 0.19197129710138724125E-02;
    w[233] = 0.17864463917586498247E-02;
    w[234] = 0.16561127281544526052E-02;
    w[235] = 0.15288767050877655684E-02;
    w[236] = 0.14049079956551446427E-02;
    w[237] = 0.12843824718970101768E-02;
    w[238] = 0.11674841174299594077E-02;
    w[239] = 0.10544076228633167722E-02;
    w[240] = 0.94536151685852538246E-03;
    w[241] = 0.84057143271072246365E-03;
    w[242] = 0.74028280424450333046E-03;
    w[243] = 0.64476204130572477933E-03;
    w[244] = 0.55429531493037471492E-03;
    w[245] = 0.46918492424785040975E-03;
    w[246] = 0.38974528447328229322E-03;
    w[247] = 0.31630366082226447689E-03;
    w[248] = 0.24921240048299729402E-03;
    w[249] = 0.18887326450650491366E-03;
    w[250] = 0.13575491094922871973E-03;
    w[251] = 0.90372734658751149261E-04;
    w[252] = 0.53275293669780613125E-04;
    w[253] = 0.25157870384280661489E-04;
    w[254] = 0.69379364324108267170E-05;
  }
  else if ( n == 511 )
  {
    x[  0] = -0.999999672956734384381;
    x[  1] = -0.999997596379748464620;
    x[  2] = -0.999992298136257588028;
    x[  3] = -0.999982430354891598580;
    x[  4] = -0.999966730098486276883;
    x[  5] = -0.999943996207054375764;
    x[  6] = -0.999913081144678282800;
    x[  7] = -0.999872888120357611938;
    x[  8] = -0.999822363679787739196;
    x[  9] = -0.999760490924432047330;
    x[ 10] = -0.999686286448317731776;
    x[ 11] = -0.999598799671910683252;
    x[ 12] = -0.999497112467187190535;
    x[ 13] = -0.999380338025023581928;
    x[ 14] = -0.999247618943342473599;
    x[ 15] = -0.999098124967667597662;
    x[ 16] = -0.998931050830810562236;
    x[ 17] = -0.998745614468095114704;
    x[ 18] = -0.998541055697167906027;
    x[ 19] = -0.998316635318407392531;
    x[ 20] = -0.998071634524930323302;
    x[ 21] = -0.997805354495957274562;
    x[ 22] = -0.997517116063472399965;
    x[ 23] = -0.997206259372221959076;
    x[ 24] = -0.996872143485260161299;
    x[ 25] = -0.996514145914890273849;
    x[ 26] = -0.996131662079315037786;
    x[ 27] = -0.995724104698407188509;
    x[ 28] = -0.995290903148810302261;
    x[ 29] = -0.994831502800621000519;
    x[ 30] = -0.994345364356723405931;
    x[ 31] = -0.993831963212755022209;
    x[ 32] = -0.993290788851684966211;
    x[ 33] = -0.992721344282788615328;
    x[ 34] = -0.992123145530863117683;
    x[ 35] = -0.991495721178106132399;
    x[ 36] = -0.990838611958294243677;
    x[ 37] = -0.990151370400770159181;
    x[ 38] = -0.989433560520240838716;
    x[ 39] = -0.988684757547429479939;
    x[ 40] = -0.987904547695124280467;
    x[ 41] = -0.987092527954034067190;
    x[ 42] = -0.986248305913007552681;
    x[ 43] = -0.985371499598520371114;
    x[ 44] = -0.984461737328814534596;
    x[ 45] = -0.983518657578632728762;
    x[ 46] = -0.982541908851080604251;
    x[ 47] = -0.981531149553740106867;
    x[ 48] = -0.980486047876721339416;
    x[ 49] = -0.979406281670862683806;
    x[ 50] = -0.978291538324758539526;
    x[ 51] = -0.977141514639705714156;
    x[ 52] = -0.975955916702011753129;
    x[ 53] = -0.974734459752402667761;
    x[ 54] = -0.973476868052506926773;
    x[ 55] = -0.972182874748581796578;
    x[ 56] = -0.970852221732792443256;
    x[ 57] = -0.969484659502459231771;
    x[ 58] = -0.968079947017759947964;
    x[ 59] = -0.966637851558416567092;
    x[ 60] = -0.965158148579915665979;
    x[ 61] = -0.963640621569812132521;
    x[ 62] = -0.962085061904651475741;
    x[ 63] = -0.960491268708020283423;
    x[ 64] = -0.958859048710200221356;
    x[ 65] = -0.957188216109860962736;
    x[ 66] = -0.955478592438183697574;
    x[ 67] = -0.953730006425761136415;
    x[ 68] = -0.951942293872573589498;
    x[ 69] = -0.950115297521294876558;
    x[ 70] = -0.948248866934137357063;
    x[ 71] = -0.946342858373402905148;
    x[ 72] = -0.944397134685866648591;
    x[ 73] = -0.942411565191083059813;
    x[ 74] = -0.940386025573669721370;
    x[ 75] = -0.938320397779592883655;
    x[ 76] = -0.936214569916450806625;
    x[ 77] = -0.934068436157725787999;
    x[ 78] = -0.931881896650953639345;
    x[ 79] = -0.929654857429740056670;
    x[ 80] = -0.927387230329536696843;
    x[ 81] = -0.925078932907075652364;
    x[ 82] = -0.922729888363349241523;
    x[ 83] = -0.920340025470012420730;
    x[ 84] = -0.917909278499077501636;
    x[ 85] = -0.915437587155765040644;
    x[ 86] = -0.912924896514370590080;
    x[ 87] = -0.910371156957004292498;
    x[ 88] = -0.907776324115058903624;
    x[ 89] = -0.905140358813261595189;
    x[ 90] = -0.902463227016165675048;
    x[ 91] = -0.899744899776940036639;
    x[ 92] = -0.896985353188316590376;
    x[ 93] = -0.894184568335559022859;
    x[ 94] = -0.891342531251319871666;
    x[ 95] = -0.888459232872256998890;
    x[ 96] = -0.885534668997285008926;
    x[ 97] = -0.882568840247341906842;
    x[ 98] = -0.879561752026556262568;
    x[ 99] = -0.876513414484705269742;
    x[100] = -0.873423842480859310192;
    x[101] = -0.870293055548113905851;
    x[102] = -0.867121077859315215614;
    x[103] = -0.863907938193690477146;
    x[104] = -0.860653669904299969802;
    x[105] = -0.857358310886232156525;
    x[106] = -0.854021903545468625813;
    x[107] = -0.850644494768350279758;
    x[108] = -0.847226135891580884381;
    x[109] = -0.843766882672708601038;
    x[110] = -0.840266795261030442350;
    x[111] = -0.836725938168868735503;
    x[112] = -0.833144380243172624728;
    x[113] = -0.829522194637401400178;
    x[114] = -0.825859458783650001088;
    x[115] = -0.822156254364980407373;
    x[116] = -0.818412667287925807395;
    x[117] = -0.814628787655137413436;
    x[118] = -0.810804709738146594361;
    x[119] = -0.806940531950217611856;
    x[120] = -0.803036356819268687782;
    x[121] = -0.799092290960841401800;
    x[122] = -0.795108445051100526780;
    x[123] = -0.791084933799848361435;
    x[124] = -0.787021875923539422170;
    x[125] = -0.782919394118283016385;
    x[126] = -0.778777615032822744702;
    x[127] = -0.774596669241483377036;
    x[128] = -0.770376691217076824278;
    x[129] = -0.766117819303760090717;
    x[130] = -0.761820195689839149173;
    x[131] = -0.757483966380513637926;
    x[132] = -0.753109281170558142523;
    x[133] = -0.748696293616936602823;
    x[134] = -0.744245161011347082309;
    x[135] = -0.739756044352694758677;
    x[136] = -0.735229108319491547663;
    x[137] = -0.730664521242181261329;
    x[138] = -0.726062455075389632685;
    x[139] = -0.721423085370098915485;
    x[140] = -0.716746591245747095767;
    x[141] = -0.712033155362252034587;
    x[142] = -0.707282963891961103412;
    x[143] = -0.702496206491527078610;
    x[144] = -0.697673076273711232906;
    x[145] = -0.692813769779114702895;
    x[146] = -0.687918486947839325756;
    x[147] = -0.682987431091079228087;
    x[148] = -0.678020808862644517838;
    x[149] = -0.673018830230418479199;
    x[150] = -0.667981708447749702165;
    x[151] = -0.662909660024780595461;
    x[152] = -0.657802904699713735422;
    x[153] = -0.652661665410017496101;
    x[154] = -0.647486168263572388782;
    x[155] = -0.642276642509759513774;
    x[156] = -0.637033320510492495071;
    x[157] = -0.631756437711194230414;
    x[158] = -0.626446232611719746542;
    x[159] = -0.621102946737226402941;
    x[160] = -0.615726824608992638014;
    x[161] = -0.610318113715186400156;
    x[162] = -0.604877064481584353319;
    x[163] = -0.599403930242242892974;
    x[164] = -0.593898967210121954393;
    x[165] = -0.588362434447662541434;
    x[166] = -0.582794593837318850840;
    x[167] = -0.577195710052045814844;
    x[168] = -0.571566050525742833992;
    x[169] = -0.565905885423654422623;
    x[170] = -0.560215487612728441818;
    x[171] = -0.554495132631932548866;
    x[172] = -0.548745098662529448608;
    x[173] = -0.542965666498311490492;
    x[174] = -0.537157119515795115982;
    x[175] = -0.531319743644375623972;
    x[176] = -0.525453827336442687395;
    x[177] = -0.519559661537457021993;
    x[178] = -0.513637539655988578507;
    x[179] = -0.507687757533716602155;
    x[180] = -0.501710613415391878251;
    x[181] = -0.495706407918761460170;
    x[182] = -0.489675444004456155436;
    x[183] = -0.483618026945841027562;
    x[184] = -0.477534464298829155284;
    x[185] = -0.471425065871658876934;
    x[186] = -0.465290143694634735858;
    x[187] = -0.459130011989832332874;
    x[188] = -0.452944987140767283784;
    x[189] = -0.446735387662028473742;
    x[190] = -0.440501534168875795783;
    x[191] = -0.434243749346802558002;
    x[192] = -0.427962357921062742583;
    x[193] = -0.421657686626163300056;
    x[194] = -0.415330064175321663764;
    x[195] = -0.408979821229888672409;
    x[196] = -0.402607290368737092671;
    x[197] = -0.396212806057615939183;
    x[198] = -0.389796704618470795479;
    x[199] = -0.383359324198730346916;
    x[200] = -0.376901004740559344802;
    x[201] = -0.370422087950078230138;
    x[202] = -0.363922917266549655269;
    x[203] = -0.357403837831532152376;
    x[204] = -0.350865196458001209011;
    x[205] = -0.344307341599438022777;
    x[206] = -0.337730623318886219621;
    x[207] = -0.331135393257976833093;
    x[208] = -0.324522004605921855207;
    x[209] = -0.317890812068476683182;
    x[210] = -0.311242171836871800300;
    x[211] = -0.304576441556714043335;
    x[212] = -0.297893980296857823437;
    x[213] = -0.291195148518246681964;
    x[214] = -0.284480308042725577496;
    x[215] = -0.277749822021824315065;
    x[216] = -0.271004054905512543536;
    x[217] = -0.264243372410926761945;
    x[218] = -0.257468141491069790481;
    x[219] = -0.250678730303483176613;
    x[220] = -0.243875508178893021593;
    x[221] = -0.237058845589829727213;
    x[222] = -0.230229114119222177156;
    x[223] = -0.223386686428966881628;
    x[224] = -0.216531936228472628081;
    x[225] = -0.209665238243181194766;
    x[226] = -0.202786968183064697557;
    x[227] = -0.195897502711100153915;
    x[228] = -0.188997219411721861059;
    x[229] = -0.182086496759252198246;
    x[230] = -0.175165714086311475707;
    x[231] = -0.168235251552207464982;
    x[232] = -0.161295490111305257361;
    x[233] = -0.154346811481378108692;
    x[234] = -0.147389598111939940054;
    x[235] = -0.140424233152560174594;
    x[236] = -0.133451100421161601344;
    x[237] = -0.126470584372301966851;
    x[238] = -0.119483070065440005133;
    x[239] = -0.112488943133186625746;
    x[240] = -0.105488589749541988533;
    x[241] = -0.984823965981192020903E-01;
    x[242] = -0.914707508403553909095E-01;
    x[243] = -0.844540400837108837102E-01;
    x[244] = -0.774326523498572825675E-01;
    x[245] = -0.704069760428551790633E-01;
    x[246] = -0.633773999173222898797E-01;
    x[247] = -0.563443130465927899720E-01;
    x[248] = -0.493081047908686267156E-01;
    x[249] = -0.422691647653636032124E-01;
    x[250] = -0.352278828084410232603E-01;
    x[251] = -0.281846489497456943394E-01;
    x[252] = -0.211398533783310883350E-01;
    x[253] = -0.140938864107824626142E-01;
    x[254] = -0.704713845933674648514E-02;
    x[255] = +0.000000000000000000000;
    x[256] = +0.704713845933674648514E-02;
    x[257] = +0.140938864107824626142E-01;
    x[258] = +0.211398533783310883350E-01;
    x[259] = +0.281846489497456943394E-01;
    x[260] = +0.352278828084410232603E-01;
    x[261] = +0.422691647653636032124E-01;
    x[262] = +0.493081047908686267156E-01;
    x[263] = +0.563443130465927899720E-01;
    x[264] = +0.633773999173222898797E-01;
    x[265] = +0.704069760428551790633E-01;
    x[266] = +0.774326523498572825675E-01;
    x[267] = +0.844540400837108837102E-01;
    x[268] = +0.914707508403553909095E-01;
    x[269] = +0.984823965981192020903E-01;
    x[270] = +0.105488589749541988533;
    x[271] = +0.112488943133186625746;
    x[272] = +0.119483070065440005133;
    x[273] = +0.126470584372301966851;
    x[274] = +0.133451100421161601344;
    x[275] = +0.140424233152560174594;
    x[276] = +0.147389598111939940054;
    x[277] = +0.154346811481378108692;
    x[278] = +0.161295490111305257361;
    x[279] = +0.168235251552207464982;
    x[280] = +0.175165714086311475707;
    x[281] = +0.182086496759252198246;
    x[282] = +0.188997219411721861059;
    x[283] = +0.195897502711100153915;
    x[284] = +0.202786968183064697557;
    x[285] = +0.209665238243181194766;
    x[286] = +0.216531936228472628081;
    x[287] = +0.223386686428966881628;
    x[288] = +0.230229114119222177156;
    x[289] = +0.237058845589829727213;
    x[290] = +0.243875508178893021593;
    x[291] = +0.250678730303483176613;
    x[292] = +0.257468141491069790481;
    x[293] = +0.264243372410926761945;
    x[294] = +0.271004054905512543536;
    x[295] = +0.277749822021824315065;
    x[296] = +0.284480308042725577496;
    x[297] = +0.291195148518246681964;
    x[298] = +0.297893980296857823437;
    x[299] = +0.304576441556714043335;
    x[300] = +0.311242171836871800300;
    x[301] = +0.317890812068476683182;
    x[302] = +0.324522004605921855207;
    x[303] = +0.331135393257976833093;
    x[304] = +0.337730623318886219621;
    x[305] = +0.344307341599438022777;
    x[306] = +0.350865196458001209011;
    x[307] = +0.357403837831532152376;
    x[308] = +0.363922917266549655269;
    x[309] = +0.370422087950078230138;
    x[310] = +0.376901004740559344802;
    x[311] = +0.383359324198730346916;
    x[312] = +0.389796704618470795479;
    x[313] = +0.396212806057615939183;
    x[314] = +0.402607290368737092671;
    x[315] = +0.408979821229888672409;
    x[316] = +0.415330064175321663764;
    x[317] = +0.421657686626163300056;
    x[318] = +0.427962357921062742583;
    x[319] = +0.434243749346802558002;
    x[320] = +0.440501534168875795783;
    x[321] = +0.446735387662028473742;
    x[322] = +0.452944987140767283784;
    x[323] = +0.459130011989832332874;
    x[324] = +0.465290143694634735858;
    x[325] = +0.471425065871658876934;
    x[326] = +0.477534464298829155284;
    x[327] = +0.483618026945841027562;
    x[328] = +0.489675444004456155436;
    x[329] = +0.495706407918761460170;
    x[330] = +0.501710613415391878251;
    x[331] = +0.507687757533716602155;
    x[332] = +0.513637539655988578507;
    x[333] = +0.519559661537457021993;
    x[334] = +0.525453827336442687395;
    x[335] = +0.531319743644375623972;
    x[336] = +0.537157119515795115982;
    x[337] = +0.542965666498311490492;
    x[338] = +0.548745098662529448608;
    x[339] = +0.554495132631932548866;
    x[340] = +0.560215487612728441818;
    x[341] = +0.565905885423654422623;
    x[342] = +0.571566050525742833992;
    x[343] = +0.577195710052045814844;
    x[344] = +0.582794593837318850840;
    x[345] = +0.588362434447662541434;
    x[346] = +0.593898967210121954393;
    x[347] = +0.599403930242242892974;
    x[348] = +0.604877064481584353319;
    x[349] = +0.610318113715186400156;
    x[350] = +0.615726824608992638014;
    x[351] = +0.621102946737226402941;
    x[352] = +0.626446232611719746542;
    x[353] = +0.631756437711194230414;
    x[354] = +0.637033320510492495071;
    x[355] = +0.642276642509759513774;
    x[356] = +0.647486168263572388782;
    x[357] = +0.652661665410017496101;
    x[358] = +0.657802904699713735422;
    x[359] = +0.662909660024780595461;
    x[360] = +0.667981708447749702165;
    x[361] = +0.673018830230418479199;
    x[362] = +0.678020808862644517838;
    x[363] = +0.682987431091079228087;
    x[364] = +0.687918486947839325756;
    x[365] = +0.692813769779114702895;
    x[366] = +0.697673076273711232906;
    x[367] = +0.702496206491527078610;
    x[368] = +0.707282963891961103412;
    x[369] = +0.712033155362252034587;
    x[370] = +0.716746591245747095767;
    x[371] = +0.721423085370098915485;
    x[372] = +0.726062455075389632685;
    x[373] = +0.730664521242181261329;
    x[374] = +0.735229108319491547663;
    x[375] = +0.739756044352694758677;
    x[376] = +0.744245161011347082309;
    x[377] = +0.748696293616936602823;
    x[378] = +0.753109281170558142523;
    x[379] = +0.757483966380513637926;
    x[380] = +0.761820195689839149173;
    x[381] = +0.766117819303760090717;
    x[382] = +0.770376691217076824278;
    x[383] = +0.774596669241483377036;
    x[384] = +0.778777615032822744702;
    x[385] = +0.782919394118283016385;
    x[386] = +0.787021875923539422170;
    x[387] = +0.791084933799848361435;
    x[388] = +0.795108445051100526780;
    x[389] = +0.799092290960841401800;
    x[390] = +0.803036356819268687782;
    x[391] = +0.806940531950217611856;
    x[392] = +0.810804709738146594361;
    x[393] = +0.814628787655137413436;
    x[394] = +0.818412667287925807395;
    x[395] = +0.822156254364980407373;
    x[396] = +0.825859458783650001088;
    x[397] = +0.829522194637401400178;
    x[398] = +0.833144380243172624728;
    x[399] = +0.836725938168868735503;
    x[400] = +0.840266795261030442350;
    x[401] = +0.843766882672708601038;
    x[402] = +0.847226135891580884381;
    x[403] = +0.850644494768350279758;
    x[404] = +0.854021903545468625813;
    x[405] = +0.857358310886232156525;
    x[406] = +0.860653669904299969802;
    x[407] = +0.863907938193690477146;
    x[408] = +0.867121077859315215614;
    x[409] = +0.870293055548113905851;
    x[410] = +0.873423842480859310192;
    x[411] = +0.876513414484705269742;
    x[412] = +0.879561752026556262568;
    x[413] = +0.882568840247341906842;
    x[414] = +0.885534668997285008926;
    x[415] = +0.888459232872256998890;
    x[416] = +0.891342531251319871666;
    x[417] = +0.894184568335559022859;
    x[418] = +0.896985353188316590376;
    x[419] = +0.899744899776940036639;
    x[420] = +0.902463227016165675048;
    x[421] = +0.905140358813261595189;
    x[422] = +0.907776324115058903624;
    x[423] = +0.910371156957004292498;
    x[424] = +0.912924896514370590080;
    x[425] = +0.915437587155765040644;
    x[426] = +0.917909278499077501636;
    x[427] = +0.920340025470012420730;
    x[428] = +0.922729888363349241523;
    x[429] = +0.925078932907075652364;
    x[430] = +0.927387230329536696843;
    x[431] = +0.929654857429740056670;
    x[432] = +0.931881896650953639345;
    x[433] = +0.934068436157725787999;
    x[434] = +0.936214569916450806625;
    x[435] = +0.938320397779592883655;
    x[436] = +0.940386025573669721370;
    x[437] = +0.942411565191083059813;
    x[438] = +0.944397134685866648591;
    x[439] = +0.946342858373402905148;
    x[440] = +0.948248866934137357063;
    x[441] = +0.950115297521294876558;
    x[442] = +0.951942293872573589498;
    x[443] = +0.953730006425761136415;
    x[444] = +0.955478592438183697574;
    x[445] = +0.957188216109860962736;
    x[446] = +0.958859048710200221356;
    x[447] = +0.960491268708020283423;
    x[448] = +0.962085061904651475741;
    x[449] = +0.963640621569812132521;
    x[450] = +0.965158148579915665979;
    x[451] = +0.966637851558416567092;
    x[452] = +0.968079947017759947964;
    x[453] = +0.969484659502459231771;
    x[454] = +0.970852221732792443256;
    x[455] = +0.972182874748581796578;
    x[456] = +0.973476868052506926773;
    x[457] = +0.974734459752402667761;
    x[458] = +0.975955916702011753129;
    x[459] = +0.977141514639705714156;
    x[460] = +0.978291538324758539526;
    x[461] = +0.979406281670862683806;
    x[462] = +0.980486047876721339416;
    x[463] = +0.981531149553740106867;
    x[464] = +0.982541908851080604251;
    x[465] = +0.983518657578632728762;
    x[466] = +0.984461737328814534596;
    x[467] = +0.985371499598520371114;
    x[468] = +0.986248305913007552681;
    x[469] = +0.987092527954034067190;
    x[470] = +0.987904547695124280467;
    x[471] = +0.988684757547429479939;
    x[472] = +0.989433560520240838716;
    x[473] = +0.990151370400770159181;
    x[474] = +0.990838611958294243677;
    x[475] = +0.991495721178106132399;
    x[476] = +0.992123145530863117683;
    x[477] = +0.992721344282788615328;
    x[478] = +0.993290788851684966211;
    x[479] = +0.993831963212755022209;
    x[480] = +0.994345364356723405931;
    x[481] = +0.994831502800621000519;
    x[482] = +0.995290903148810302261;
    x[483] = +0.995724104698407188509;
    x[484] = +0.996131662079315037786;
    x[485] = +0.996514145914890273849;
    x[486] = +0.996872143485260161299;
    x[487] = +0.997206259372221959076;
    x[488] = +0.997517116063472399965;
    x[489] = +0.997805354495957274562;
    x[490] = +0.998071634524930323302;
    x[491] = +0.998316635318407392531;
    x[492] = +0.998541055697167906027;
    x[493] = +0.998745614468095114704;
    x[494] = +0.998931050830810562236;
    x[495] = +0.999098124967667597662;
    x[496] = +0.999247618943342473599;
    x[497] = +0.999380338025023581928;
    x[498] = +0.999497112467187190535;
    x[499] = +0.999598799671910683252;
    x[500] = +0.999686286448317731776;
    x[501] = +0.999760490924432047330;
    x[502] = +0.999822363679787739196;
    x[503] = +0.999872888120357611938;
    x[504] = +0.999913081144678282800;
    x[505] = +0.999943996207054375764;
    x[506] = +0.999966730098486276883;
    x[507] = +0.999982430354891598580;
    x[508] = +0.999992298136257588028;
    x[509] = +0.999997596379748464620;
    x[510] = +0.999999672956734384381;

    w[  0] = 0.945715933950007048827E-06;
    w[  1] = 0.345456507169149134898E-05;
    w[  2] = 0.736624069102321668857E-05;
    w[  3] = 0.125792781889592743525E-04;
    w[  4] = 0.190213681905875816679E-04;
    w[  5] = 0.266376412339000901358E-04;
    w[  6] = 0.353751372055189588628E-04;
    w[  7] = 0.451863674126296143105E-04;
    w[  8] = 0.560319507856164252140E-04;
    w[  9] = 0.678774554733972416227E-04;
    w[ 10] = 0.806899228014035293851E-04;
    w[ 11] = 0.944366322532705527066E-04;
    w[ 12] = 0.109085545645741522051E-03;
    w[ 13] = 0.124606200241498368482E-03;
    w[ 14] = 0.140970302204104791413E-03;
    w[ 15] = 0.158151830411132242924E-03;
    w[ 16] = 0.176126765545083195474E-03;
    w[ 17] = 0.194872642236641146532E-03;
    w[ 18] = 0.214368090034216937149E-03;
    w[ 19] = 0.234592462123925204879E-03;
    w[ 20] = 0.255525589595236862014E-03;
    w[ 21] = 0.277147657465187357459E-03;
    w[ 22] = 0.299439176850911730874E-03;
    w[ 23] = 0.322381020652862389664E-03;
    w[ 24] = 0.345954492129903871350E-03;
    w[ 25] = 0.370141402122251665232E-03;
    w[ 26] = 0.394924138246873704434E-03;
    w[ 27] = 0.420285716355361231823E-03;
    w[ 28] = 0.446209810101403247488E-03;
    w[ 29] = 0.472680758429262691232E-03;
    w[ 30] = 0.499683553312800484519E-03;
    w[ 31] = 0.527203811431658386125E-03;
    w[ 32] = 0.555227733977307579715E-03;
    w[ 33] = 0.583742058714979703847E-03;
    w[ 34] = 0.612734008012225209294E-03;
    w[ 35] = 0.642191235948505088403E-03;
    w[ 36] = 0.672101776960108194646E-03;
    w[ 37] = 0.702453997827572321358E-03;
    w[ 38] = 0.733236554224767912055E-03;
    w[ 39] = 0.764438352543882784191E-03;
    w[ 40] = 0.796048517297550871506E-03;
    w[ 41] = 0.828056364077226302608E-03;
    w[ 42] = 0.860451377808527848128E-03;
    w[ 43] = 0.893223195879324912340E-03;
    w[ 44] = 0.926361595613111283368E-03;
    w[ 45] = 0.959856485506936206261E-03;
    w[ 46] = 0.993697899638760857945E-03;
    w[ 47] = 0.102787599466367326179E-02;
    w[ 48] = 0.106238104885340071375E-02;
    w[ 49] = 0.109720346268191941940E-02;
    w[ 50] = 0.113233376051597664917E-02;
    w[ 51] = 0.116776259302858043685E-02;
    w[ 52] = 0.120348074001265964881E-02;
    w[ 53] = 0.123947911332878396534E-02;
    w[ 54] = 0.127574875977346947345E-02;
    w[ 55] = 0.131228086370221478128E-02;
    w[ 56] = 0.134906674928353113127E-02;
    w[ 57] = 0.138609788229672549700E-02;
    w[ 58] = 0.142336587141720519900E-02;
    w[ 59] = 0.146086246895890987689E-02;
    w[ 60] = 0.149857957106456636214E-02;
    w[ 61] = 0.153650921735128916170E-02;
    w[ 62] = 0.157464359003212166189E-02;
    w[ 63] = 0.161297501254393423070E-02;
    w[ 64] = 0.165149594771914570655E-02;
    w[ 65] = 0.169019899554346019117E-02;
    w[ 66] = 0.172907689054461607168E-02;
    w[ 67] = 0.176812249885838886701E-02;
    w[ 68] = 0.180732881501808930079E-02;
    w[ 69] = 0.184668895851282540913E-02;
    w[ 70] = 0.188619617015808475394E-02;
    w[ 71] = 0.192584380831993546204E-02;
    w[ 72] = 0.196562534503150547732E-02;
    w[ 73] = 0.200553436203751169944E-02;
    w[ 74] = 0.204556454679958293446E-02;
    w[ 75] = 0.208570968849203942640E-02;
    w[ 76] = 0.212596367401472533045E-02;
    w[ 77] = 0.216632048404649142727E-02;
    w[ 78] = 0.220677418916003329194E-02;
    w[ 79] = 0.224731894601603393082E-02;
    w[ 80] = 0.228794899365195972378E-02;
    w[ 81] = 0.232865864987842738864E-02;
    w[ 82] = 0.236944230779380495146E-02;
    w[ 83] = 0.241029443242563417382E-02;
    w[ 84] = 0.245120955750556483923E-02;
    w[ 85] = 0.249218228238276930060E-02;
    w[ 86] = 0.253320726907925325750E-02;
    w[ 87] = 0.257427923948908888092E-02;
    w[ 88] = 0.261539297272236109225E-02;
    w[ 89] = 0.265654330259352828314E-02;
    w[ 90] = 0.269772511525294586667E-02;
    w[ 91] = 0.273893334695947541201E-02;
    w[ 92] = 0.278016298199139435045E-02;
    w[ 93] = 0.282140905069222207923E-02;
    w[ 94] = 0.286266662764757868253E-02;
    w[ 95] = 0.290393082998878368175E-02;
    w[ 96] = 0.294519681581857582284E-02;
    w[ 97] = 0.298645978275408290247E-02;
    w[ 98] = 0.302771496658198544480E-02;
    w[ 99] = 0.306895764002069252174E-02;
    w[100] = 0.311018311158427546158E-02;
    w[101] = 0.315138672454287935858E-02;
    w[102] = 0.319256385597434736790E-02;
    w[103] = 0.323370991590184336368E-02;
    w[104] = 0.327482034651233969564E-02;
    w[105] = 0.331589062145094394706E-02;
    w[106] = 0.335691624518616761342E-02;
    w[107] = 0.339789275244138669739E-02;
    w[108] = 0.343881570768790591876E-02;
    w[109] = 0.347968070469521146972E-02;
    w[110] = 0.352048336613417922682E-02;
    w[111] = 0.356121934322919357659E-02;
    w[112] = 0.360188431545532431869E-02;
    w[113] = 0.364247399027690353194E-02;
    w[114] = 0.368298410292403911967E-02;
    w[115] = 0.372341041620379550870E-02;
    w[116] = 0.376374872034296338241E-02;
    w[117] = 0.380399483285952829161E-02;
    w[118] = 0.384414459846013158917E-02;
    w[119] = 0.388419388896099560998E-02;
    w[120] = 0.392413860322995774660E-02;
    w[121] = 0.396397466714742455513E-02;
    w[122] = 0.400369803358421688562E-02;
    w[123] = 0.404330468239442998549E-02;
    w[124] = 0.408279062042157838350E-02;
    w[125] = 0.412215188151643401528E-02;
    w[126] = 0.416138452656509745764E-02;
    w[127] = 0.420048464352596631772E-02;
    w[128] = 0.423944834747438184434E-02;
    w[129] = 0.427827178065384480959E-02;
    w[130] = 0.431695111253279479928E-02;
    w[131] = 0.435548253986604343679E-02;
    w[132] = 0.439386228676004195260E-02;
    w[133] = 0.443208660474124713206E-02;
    w[134] = 0.447015177282692726900E-02;
    w[135] = 0.450805409759782158001E-02;
    w[136] = 0.454578991327213285488E-02;
    w[137] = 0.458335558178039420335E-02;
    w[138] = 0.462074749284080687482E-02;
    w[139] = 0.465796206403469754658E-02;
    w[140] = 0.469499574088179046532E-02;
    w[141] = 0.473184499691503264714E-02;
    w[142] = 0.476850633375474925263E-02;
    w[143] = 0.480497628118194150483E-02;
    w[144] = 0.484125139721057135214E-02;
    w[145] = 0.487732826815870573054E-02;
    w[146] = 0.491320350871841897367E-02;
    w[147] = 0.494887376202437487201E-02;
    w[148] = 0.498433569972103029914E-02;
    w[149] = 0.501958602202842039909E-02;
    w[150] = 0.505462145780650125058E-02;
    w[151] = 0.508943876461803986674E-02;
    w[152] = 0.512403472879005351831E-02;
    w[153] = 0.515840616547381084096E-02;
    w[154] = 0.519254991870341614863E-02;
    w[155] = 0.522646286145300596306E-02;
    w[156] = 0.526014189569259311205E-02;
    w[157] = 0.529358395244259896547E-02;
    w[158] = 0.532678599182711857974E-02;
    w[159] = 0.535974500312596681161E-02;
    w[160] = 0.539245800482555593606E-02;
    w[161] = 0.542492204466865704951E-02;
    w[162] = 0.545713419970309863995E-02;
    w[163] = 0.548909157632945623482E-02;
    w[164] = 0.552079131034778706457E-02;
    w[165] = 0.555223056700346326850E-02;
    w[166] = 0.558340654103215637610E-02;
    w[167] = 0.561431645670402467678E-02;
    w[168] = 0.564495756786715368885E-02;
    w[169] = 0.567532715799029830087E-02;
    w[170] = 0.570542254020497332312E-02;
    w[171] = 0.573524105734693719020E-02;
    w[172] = 0.576478008199711142954E-02;
    w[173] = 0.579403701652197628421E-02;
    w[174] = 0.582300929311348057702E-02;
    w[175] = 0.585169437382850155033E-02;
    w[176] = 0.588008975062788803205E-02;
    w[177] = 0.590819294541511788161E-02;
    w[178] = 0.593600151007459827614E-02;
    w[179] = 0.596351302650963502011E-02;
    w[180] = 0.599072510668009471472E-02;
    w[181] = 0.601763539263978131522E-02;
    w[182] = 0.604424155657354634589E-02;
    w[183] = 0.607054130083414983949E-02;
    w[184] = 0.609653235797888692923E-02;
    w[185] = 0.612221249080599294931E-02;
    w[186] = 0.614757949239083790214E-02;
    w[187] = 0.617263118612191922727E-02;
    w[188] = 0.619736542573665996342E-02;
    w[189] = 0.622178009535701763157E-02;
    w[190] = 0.624587310952490748541E-02;
    w[191] = 0.626964241323744217671E-02;
    w[192] = 0.629308598198198836688E-02;
    w[193] = 0.631620182177103938227E-02;
    w[194] = 0.633898796917690165912E-02;
    w[195] = 0.636144249136619145314E-02;
    w[196] = 0.638356348613413709795E-02;
    w[197] = 0.640534908193868098342E-02;
    w[198] = 0.642679743793437438922E-02;
    w[199] = 0.644790674400605734710E-02;
    w[200] = 0.646867522080231481688E-02;
    w[201] = 0.648910111976869964292E-02;
    w[202] = 0.650918272318071200827E-02;
    w[203] = 0.652891834417652442012E-02;
    w[204] = 0.654830632678944064054E-02;
    w[205] = 0.656734504598007641819E-02;
    w[206] = 0.658603290766824937794E-02;
    w[207] = 0.660436834876456498276E-02;
    w[208] = 0.662234983720168509457E-02;
    w[209] = 0.663997587196526532519E-02;
    w[210] = 0.665724498312454708217E-02;
    w[211] = 0.667415573186258997654E-02;
    w[212] = 0.669070671050613006584E-02;
    w[213] = 0.670689654255504925648E-02;
    w[214] = 0.672272388271144108036E-02;
    w[215] = 0.673818741690825799086E-02;
    w[216] = 0.675328586233752529078E-02;
    w[217] = 0.676801796747810680683E-02;
    w[218] = 0.678238251212300746082E-02;
    w[219] = 0.679637830740619795480E-02;
    w[220] = 0.681000419582894688374E-02;
    w[221] = 0.682325905128564571420E-02;
    w[222] = 0.683614177908911221841E-02;
    w[223] = 0.684865131599535812903E-02;
    w[224] = 0.686078663022780697951E-02;
    w[225] = 0.687254672150094831613E-02;
    w[226] = 0.688393062104341470995E-02;
    w[227] = 0.689493739162046825872E-02;
    w[228] = 0.690556612755588354803E-02;
    w[229] = 0.691581595475321433825E-02;
    w[230] = 0.692568603071643155621E-02;
    w[231] = 0.693517554456992049848E-02;
    w[232] = 0.694428371707782549438E-02;
    w[233] = 0.695300980066273063177E-02;
    w[234] = 0.696135307942366551493E-02;
    w[235] = 0.696931286915342540213E-02;
    w[236] = 0.697688851735519545845E-02;
    w[237] = 0.698407940325846925786E-02;
    w[238] = 0.699088493783425207545E-02;
    w[239] = 0.699730456380953992594E-02;
    w[240] = 0.700333775568106572820E-02;
    w[241] = 0.700898401972830440494E-02;
    w[242] = 0.701424289402572916425E-02;
    w[243] = 0.701911394845431165171E-02;
    w[244] = 0.702359678471225911031E-02;
    w[245] = 0.702769103632498213858E-02;
    w[246] = 0.703139636865428709508E-02;
    w[247] = 0.703471247890678765907E-02;
    w[248] = 0.703763909614153052319E-02;
    w[249] = 0.704017598127683066242E-02;
    w[250] = 0.704232292709631209597E-02;
    w[251] = 0.704407975825415053266E-02;
    w[252] = 0.704544633127951476780E-02;
    w[253] = 0.704642253458020417748E-02;
    w[254] = 0.704700828844548013730E-02;
    w[255] = 0.704720354504808967346E-02;
    w[256] = 0.704700828844548013730E-02;
    w[257] = 0.704642253458020417748E-02;
    w[258] = 0.704544633127951476780E-02;
    w[259] = 0.704407975825415053266E-02;
    w[260] = 0.704232292709631209597E-02;
    w[261] = 0.704017598127683066242E-02;
    w[262] = 0.703763909614153052319E-02;
    w[263] = 0.703471247890678765907E-02;
    w[264] = 0.703139636865428709508E-02;
    w[265] = 0.702769103632498213858E-02;
    w[266] = 0.702359678471225911031E-02;
    w[267] = 0.701911394845431165171E-02;
    w[268] = 0.701424289402572916425E-02;
    w[269] = 0.700898401972830440494E-02;
    w[270] = 0.700333775568106572820E-02;
    w[271] = 0.699730456380953992594E-02;
    w[272] = 0.699088493783425207545E-02;
    w[273] = 0.698407940325846925786E-02;
    w[274] = 0.697688851735519545845E-02;
    w[275] = 0.696931286915342540213E-02;
    w[276] = 0.696135307942366551493E-02;
    w[277] = 0.695300980066273063177E-02;
    w[278] = 0.694428371707782549438E-02;
    w[279] = 0.693517554456992049848E-02;
    w[280] = 0.692568603071643155621E-02;
    w[281] = 0.691581595475321433825E-02;
    w[282] = 0.690556612755588354803E-02;
    w[283] = 0.689493739162046825872E-02;
    w[284] = 0.688393062104341470995E-02;
    w[285] = 0.687254672150094831613E-02;
    w[286] = 0.686078663022780697951E-02;
    w[287] = 0.684865131599535812903E-02;
    w[288] = 0.683614177908911221841E-02;
    w[289] = 0.682325905128564571420E-02;
    w[290] = 0.681000419582894688374E-02;
    w[291] = 0.679637830740619795480E-02;
    w[292] = 0.678238251212300746082E-02;
    w[293] = 0.676801796747810680683E-02;
    w[294] = 0.675328586233752529078E-02;
    w[295] = 0.673818741690825799086E-02;
    w[296] = 0.672272388271144108036E-02;
    w[297] = 0.670689654255504925648E-02;
    w[298] = 0.669070671050613006584E-02;
    w[299] = 0.667415573186258997654E-02;
    w[300] = 0.665724498312454708217E-02;
    w[301] = 0.663997587196526532519E-02;
    w[302] = 0.662234983720168509457E-02;
    w[303] = 0.660436834876456498276E-02;
    w[304] = 0.658603290766824937794E-02;
    w[305] = 0.656734504598007641819E-02;
    w[306] = 0.654830632678944064054E-02;
    w[307] = 0.652891834417652442012E-02;
    w[308] = 0.650918272318071200827E-02;
    w[309] = 0.648910111976869964292E-02;
    w[310] = 0.646867522080231481688E-02;
    w[311] = 0.644790674400605734710E-02;
    w[312] = 0.642679743793437438922E-02;
    w[313] = 0.640534908193868098342E-02;
    w[314] = 0.638356348613413709795E-02;
    w[315] = 0.636144249136619145314E-02;
    w[316] = 0.633898796917690165912E-02;
    w[317] = 0.631620182177103938227E-02;
    w[318] = 0.629308598198198836688E-02;
    w[319] = 0.626964241323744217671E-02;
    w[320] = 0.624587310952490748541E-02;
    w[321] = 0.622178009535701763157E-02;
    w[322] = 0.619736542573665996342E-02;
    w[323] = 0.617263118612191922727E-02;
    w[324] = 0.614757949239083790214E-02;
    w[325] = 0.612221249080599294931E-02;
    w[326] = 0.609653235797888692923E-02;
    w[327] = 0.607054130083414983949E-02;
    w[328] = 0.604424155657354634589E-02;
    w[329] = 0.601763539263978131522E-02;
    w[330] = 0.599072510668009471472E-02;
    w[331] = 0.596351302650963502011E-02;
    w[332] = 0.593600151007459827614E-02;
    w[333] = 0.590819294541511788161E-02;
    w[334] = 0.588008975062788803205E-02;
    w[335] = 0.585169437382850155033E-02;
    w[336] = 0.582300929311348057702E-02;
    w[337] = 0.579403701652197628421E-02;
    w[338] = 0.576478008199711142954E-02;
    w[339] = 0.573524105734693719020E-02;
    w[340] = 0.570542254020497332312E-02;
    w[341] = 0.567532715799029830087E-02;
    w[342] = 0.564495756786715368885E-02;
    w[343] = 0.561431645670402467678E-02;
    w[344] = 0.558340654103215637610E-02;
    w[345] = 0.555223056700346326850E-02;
    w[346] = 0.552079131034778706457E-02;
    w[347] = 0.548909157632945623482E-02;
    w[348] = 0.545713419970309863995E-02;
    w[349] = 0.542492204466865704951E-02;
    w[350] = 0.539245800482555593606E-02;
    w[351] = 0.535974500312596681161E-02;
    w[352] = 0.532678599182711857974E-02;
    w[353] = 0.529358395244259896547E-02;
    w[354] = 0.526014189569259311205E-02;
    w[355] = 0.522646286145300596306E-02;
    w[356] = 0.519254991870341614863E-02;
    w[357] = 0.515840616547381084096E-02;
    w[358] = 0.512403472879005351831E-02;
    w[359] = 0.508943876461803986674E-02;
    w[360] = 0.505462145780650125058E-02;
    w[361] = 0.501958602202842039909E-02;
    w[362] = 0.498433569972103029914E-02;
    w[363] = 0.494887376202437487201E-02;
    w[364] = 0.491320350871841897367E-02;
    w[365] = 0.487732826815870573054E-02;
    w[366] = 0.484125139721057135214E-02;
    w[367] = 0.480497628118194150483E-02;
    w[368] = 0.476850633375474925263E-02;
    w[369] = 0.473184499691503264714E-02;
    w[370] = 0.469499574088179046532E-02;
    w[371] = 0.465796206403469754658E-02;
    w[372] = 0.462074749284080687482E-02;
    w[373] = 0.458335558178039420335E-02;
    w[374] = 0.454578991327213285488E-02;
    w[375] = 0.450805409759782158001E-02;
    w[376] = 0.447015177282692726900E-02;
    w[377] = 0.443208660474124713206E-02;
    w[378] = 0.439386228676004195260E-02;
    w[379] = 0.435548253986604343679E-02;
    w[380] = 0.431695111253279479928E-02;
    w[381] = 0.427827178065384480959E-02;
    w[382] = 0.423944834747438184434E-02;
    w[383] = 0.420048464352596631772E-02;
    w[384] = 0.416138452656509745764E-02;
    w[385] = 0.412215188151643401528E-02;
    w[386] = 0.408279062042157838350E-02;
    w[387] = 0.404330468239442998549E-02;
    w[388] = 0.400369803358421688562E-02;
    w[389] = 0.396397466714742455513E-02;
    w[390] = 0.392413860322995774660E-02;
    w[391] = 0.388419388896099560998E-02;
    w[392] = 0.384414459846013158917E-02;
    w[393] = 0.380399483285952829161E-02;
    w[394] = 0.376374872034296338241E-02;
    w[395] = 0.372341041620379550870E-02;
    w[396] = 0.368298410292403911967E-02;
    w[397] = 0.364247399027690353194E-02;
    w[398] = 0.360188431545532431869E-02;
    w[399] = 0.356121934322919357659E-02;
    w[400] = 0.352048336613417922682E-02;
    w[401] = 0.347968070469521146972E-02;
    w[402] = 0.343881570768790591876E-02;
    w[403] = 0.339789275244138669739E-02;
    w[404] = 0.335691624518616761342E-02;
    w[405] = 0.331589062145094394706E-02;
    w[406] = 0.327482034651233969564E-02;
    w[407] = 0.323370991590184336368E-02;
    w[408] = 0.319256385597434736790E-02;
    w[409] = 0.315138672454287935858E-02;
    w[410] = 0.311018311158427546158E-02;
    w[411] = 0.306895764002069252174E-02;
    w[412] = 0.302771496658198544480E-02;
    w[413] = 0.298645978275408290247E-02;
    w[414] = 0.294519681581857582284E-02;
    w[415] = 0.290393082998878368175E-02;
    w[416] = 0.286266662764757868253E-02;
    w[417] = 0.282140905069222207923E-02;
    w[418] = 0.278016298199139435045E-02;
    w[419] = 0.273893334695947541201E-02;
    w[420] = 0.269772511525294586667E-02;
    w[421] = 0.265654330259352828314E-02;
    w[422] = 0.261539297272236109225E-02;
    w[423] = 0.257427923948908888092E-02;
    w[424] = 0.253320726907925325750E-02;
    w[425] = 0.249218228238276930060E-02;
    w[426] = 0.245120955750556483923E-02;
    w[427] = 0.241029443242563417382E-02;
    w[428] = 0.236944230779380495146E-02;
    w[429] = 0.232865864987842738864E-02;
    w[430] = 0.228794899365195972378E-02;
    w[431] = 0.224731894601603393082E-02;
    w[432] = 0.220677418916003329194E-02;
    w[433] = 0.216632048404649142727E-02;
    w[434] = 0.212596367401472533045E-02;
    w[435] = 0.208570968849203942640E-02;
    w[436] = 0.204556454679958293446E-02;
    w[437] = 0.200553436203751169944E-02;
    w[438] = 0.196562534503150547732E-02;
    w[439] = 0.192584380831993546204E-02;
    w[440] = 0.188619617015808475394E-02;
    w[441] = 0.184668895851282540913E-02;
    w[442] = 0.180732881501808930079E-02;
    w[443] = 0.176812249885838886701E-02;
    w[444] = 0.172907689054461607168E-02;
    w[445] = 0.169019899554346019117E-02;
    w[446] = 0.165149594771914570655E-02;
    w[447] = 0.161297501254393423070E-02;
    w[448] = 0.157464359003212166189E-02;
    w[449] = 0.153650921735128916170E-02;
    w[450] = 0.149857957106456636214E-02;
    w[451] = 0.146086246895890987689E-02;
    w[452] = 0.142336587141720519900E-02;
    w[453] = 0.138609788229672549700E-02;
    w[454] = 0.134906674928353113127E-02;
    w[455] = 0.131228086370221478128E-02;
    w[456] = 0.127574875977346947345E-02;
    w[457] = 0.123947911332878396534E-02;
    w[458] = 0.120348074001265964881E-02;
    w[459] = 0.116776259302858043685E-02;
    w[460] = 0.113233376051597664917E-02;
    w[461] = 0.109720346268191941940E-02;
    w[462] = 0.106238104885340071375E-02;
    w[463] = 0.102787599466367326179E-02;
    w[464] = 0.993697899638760857945E-03;
    w[465] = 0.959856485506936206261E-03;
    w[466] = 0.926361595613111283368E-03;
    w[467] = 0.893223195879324912340E-03;
    w[468] = 0.860451377808527848128E-03;
    w[469] = 0.828056364077226302608E-03;
    w[470] = 0.796048517297550871506E-03;
    w[471] = 0.764438352543882784191E-03;
    w[472] = 0.733236554224767912055E-03;
    w[473] = 0.702453997827572321358E-03;
    w[474] = 0.672101776960108194646E-03;
    w[475] = 0.642191235948505088403E-03;
    w[476] = 0.612734008012225209294E-03;
    w[477] = 0.583742058714979703847E-03;
    w[478] = 0.555227733977307579715E-03;
    w[479] = 0.527203811431658386125E-03;
    w[480] = 0.499683553312800484519E-03;
    w[481] = 0.472680758429262691232E-03;
    w[482] = 0.446209810101403247488E-03;
    w[483] = 0.420285716355361231823E-03;
    w[484] = 0.394924138246873704434E-03;
    w[485] = 0.370141402122251665232E-03;
    w[486] = 0.345954492129903871350E-03;
    w[487] = 0.322381020652862389664E-03;
    w[488] = 0.299439176850911730874E-03;
    w[489] = 0.277147657465187357459E-03;
    w[490] = 0.255525589595236862014E-03;
    w[491] = 0.234592462123925204879E-03;
    w[492] = 0.214368090034216937149E-03;
    w[493] = 0.194872642236641146532E-03;
    w[494] = 0.176126765545083195474E-03;
    w[495] = 0.158151830411132242924E-03;
    w[496] = 0.140970302204104791413E-03;
    w[497] = 0.124606200241498368482E-03;
    w[498] = 0.109085545645741522051E-03;
    w[499] = 0.944366322532705527066E-04;
    w[500] = 0.806899228014035293851E-04;
    w[501] = 0.678774554733972416227E-04;
    w[502] = 0.560319507856164252140E-04;
    w[503] = 0.451863674126296143105E-04;
    w[504] = 0.353751372055189588628E-04;
    w[505] = 0.266376412339000901358E-04;
    w[506] = 0.190213681905875816679E-04;
    w[507] = 0.125792781889592743525E-04;
    w[508] = 0.736624069102321668857E-05;
    w[509] = 0.345456507169149134898E-05;
    w[510] = 0.945715933950007048827E-06;
  }
  else
  {
    cerr << "\n";
    cerr << "PATTERSON_SET - Fatal error!\n";
    cerr << "  Illegal input value of N.\n";
    cerr << "  N must be 1, 3, 7, 15, 31, 63, 127, 255 or 511.\n";
    exit ( 1 );
  }
  return;
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

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the property 
//    that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the double precision round-off unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Formula:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL2 computes the double factorial function.
//
//  Formula:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    N!!
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial 
//    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
//
//    Output, double R8_FACTORIAL2, the value of N!!.
//
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
//
//  Discussion:
//
//    The program uses rational functions that theoretically approximate
//    LOG(GAMMA(X)) to at least 18 significant decimal digits.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody and Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
//    Mathematics of Computation,
//    Volume 21, 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    Hart, Et. Al.,
//    Computer Approximations,
//    Wiley and sons, New York, 1968.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.  X must be positive.
//
//    Output, double R8_GAMMA_LOG, the logarithm of the Gamma function of X.
//    If X <= 0.0, or if overflow would occur, the program returns the
//    value XINF, the largest representable floating point number.
//
//  Local Parameters:
//
//  BETA   - radix for the floating-point representation.
//
//  MAXEXP - the smallest positive power of BETA that overflows.
//
//  XBIG   - largest argument for which LN(GAMMA(X)) is representable
//           in the machine, i.e., the solution to the equation
//             LN(GAMMA(XBIG)) = BETA**MAXEXP.
//
//  FRTBIG - Rough estimate of the fourth root of XBIG
//
//
//  Approximate values for some important machines are:
//
//                            BETA      MAXEXP         XBIG
//
//  CRAY-1        (S.P.)        2        8191       9.62E+2461
//  Cyber 180/855
//    under NOS   (S.P.)        2        1070       1.72E+319
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)        2         128       4.08E+36
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)        2        1024       2.55D+305
//  IBM 3033      (D.P.)       16          63       4.29D+73
//  VAX D-Format  (D.P.)        2         127       2.05D+36
//  VAX G-Format  (D.P.)        2        1023       1.28D+305
//
//
//                           FRTBIG
//
//  CRAY-1        (S.P.)   3.13E+615
//  Cyber 180/855
//    under NOS   (S.P.)   6.44E+79
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)   1.42E+9
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)   2.25D+76
//  IBM 3033      (D.P.)   2.56D+18
//  VAX D-Format  (D.P.)   1.20D+9
//  VAX G-Format  (D.P.)   1.89D+76
//
{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04, 
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double corr;
  double d1 = - 5.772156649015328605195174E-01;
  double d2 =   4.227843350984671393993777E-01;
  double d4 =   1.791759469228055000094023;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = {
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+010,
    1.702665737765398868392998E+011,
    4.926125793377430887588120E+011, 
    5.606251856223951465078242E+011 };
  double pnt68 = 0.6796875;
  double q1[8] = {
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = {
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = {
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+010, 
    1.016803586272438228077304E+011, 
    3.417476345507377132798597E+011, 
    4.463158187419713286462081E+011 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
//
//  Return immediately if the argument is out of range.
//
  if ( x <= 0.0 || xbig < x )
  {
    return r8_huge ( );
  }

  if ( x <= r8_epsilon ( ) )
  {
    res = - log ( x );
  }
  else if ( x <= 1.5 )
  {
    if ( x < pnt68 )
    {
      corr = - log ( x );
      xm1 = x;
    }
    else
    {
      corr = 0.0;
      xm1 = ( x - 0.5 ) - 0.5;
    }

    if ( x <= 0.5 || pnt68 <= x )
    {
      xden = 1.0;
      xnum = 0.0;

      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm1 + p1[i];
        xden = xden * xm1 + q1[i];
      }

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
    }
    else
    {
      xm2 = ( x - 0.5 ) - 0.5;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );

    }
  }
  else if ( x <= 4.0 )
  {
    xm2 = x - 2.0;
    xden = 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm2 + p2[i];
      xden = xden * xm2 + q2[i];
    }

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
  }
  else if ( x <= 12.0 )
  {
    xm4 = x - 4.0;
    xden = - 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm4 + p4[i];
      xden = xden * xm4 + q4[i];
    }

    res = d4 + xm4 * ( xnum / xden );
  }
  else
  {
    res = 0.0;

    if ( x <= frtbig )
    {

      res = c[6];
      xsq = x * x;

      for ( i = 0; i < 6; i++ )
      {
        res = res / xsq + c[i];
      }

    }

    res = res / x;
    corr = log ( x );
    res = res + sqrtpi - 0.5 * corr;
    res = res + x * ( corr - 1.0 );

  }

  return res;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8_hyper_2f1 ( double a, double b, double c, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HYPER_2F1 evaluates the hypergeometric function 2F1(A,B,C,X).
//
//  Discussion:
//
//    A minor bug was corrected.  The HW variable, used in several places as
//    the "old" value of a quantity being iteratively improved, was not
//    being initialized.  JVB, 11 February 2008.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
//    C++ version by John Burkardt.
//
//    The original FORTRAN77 version of this routine is copyrighted by
//    Shanjie Zhang and Jianming Jin.  However, they give permission to
//    incorporate this routine into a user program provided that the copyright
//    is acknowledged.
//
//  Reference:
//
//    Shanjie Zhang, Jianming Jin,
//    Computation of Special Functions,
//    Wiley, 1996,
//    ISBN: 0-471-11963-6,
//    LC: QA351.C45
//
//  Parameters:
//
//    Input, double A, B, C, X, the arguments of the function.
//    C must not be equal to a nonpositive integer.
//    X < 1.
//
//    Output, double R8_HYPER_2F1, the value of the function.
//
{
  double a0;
  double aa;
  double bb;
  double c0;
  double c1;
  double el = 0.5772156649015329;
  double eps;
  double f0;
  double f1;
  double g0;
  double g1;
  double g2;
  double g3;
  double ga;
  double gabc;
  double gam;
  double gb;
  double gbm;
  double gc;
  double gca;
  double gcab;
  double gcb;
  double gm;
  double hf;
  double hw;
  int j;
  int k;
  bool l0;
  bool l1;
  bool l2;
  bool l3;
  bool l4;
  bool l5;
  int m;
  int nm;
  double pa;
  double pb;
  double pi = 3.141592653589793;
  double r;
  double r0;
  double r1;
  double rm;
  double rp;
  double sm;
  double sp;
  double sp0;
  double x1;

  l0 = ( c == ( int ) ( c ) ) && ( c < 0.0 );
  l1 = ( 1.0 - x < 1.0E-15 ) && ( c - a - b <= 0.0 );
  l2 = ( a == ( int ) ( a ) ) && ( a < 0.0 );
  l3 = ( b == ( int ) ( b ) ) && ( b < 0.0 );
  l4 = ( c - a == ( int ) ( c - a ) ) && ( c - a <= 0.0 );
  l5 = ( c - b == ( int ) ( c - b ) ) && ( c - b <= 0.0 );

  if ( l0 )
  {
    cerr << "\n";
    cerr << "R8_HYPER_2F1 - Fatal error!\n";
    cerr << "  The hypergeometric series is divergent.\n";
    cerr << "  C is integral and negative.\n";
    cerr << "  C = " << c << "\n";
    exit ( 1 );
  }

  if ( l1 )
  {
    cerr << "\n";
    cerr << "R8_HYPER_2F1 - Fatal error!\n";
    cerr << "  The hypergeometric series is divergent.\n";
    cerr << "  1 - X < 0, C - A - B <= 0\n";
    cerr << "  A = " << a << "\n";
    cerr << "  B = " << b << "\n";
    cerr << "  C = " << c << "\n";
    cerr << "  X = " << x << "\n";
    exit ( 1 );
  }

  if ( 0.95 < x )
  {
    eps = 1.0E-08;
  }
  else
  {
    eps = 1.0E-15;
  }

  if ( x == 0.0 || a == 0.0 || b == 0.0 )
  {
    hf = 1.0;
    return hf;
  }
  else if ( 1.0 - x == eps && 0.0 < c - a - b )
  {
    gc = gamma ( c );
    gcab = gamma ( c - a - b );
    gca = gamma ( c - a );
    gcb = gamma ( c - b );
    hf = gc * gcab / ( gca * gcb );
    return hf;
  }
  else if ( 1.0 + x <= eps && r8_abs ( c - a + b - 1.0 ) <= eps )
  {
    g0 = sqrt ( pi ) * pow ( 2.0, - a );
    g1 = gamma ( c );
    g2 = gamma ( 1.0 + a / 2.0 - b );
    g3 = gamma ( 0.5 + 0.5 * a );
    hf = g0 * g1 / ( g2 * g3 );
    return hf;
  }
  else if ( l2 || l3 )
  {
    if ( l2 )
    {
      nm = ( int ) ( r8_abs ( a ) );
    }

    if ( l3 )
    {
      nm = ( int ) ( r8_abs ( b ) );
    }

    hf = 1.0;
    r = 1.0;

    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }

    return hf;
  }
  else if ( l4 || l5 )
  {
    if ( l4 )
    {
      nm = ( int ) ( r8_abs ( c - a ) );
    }

    if ( l5 )
    {
      nm = ( int ) ( r8_abs ( c - b ) );
    }

    hf = 1.0;
    r  = 1.0;
    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }
    hf = pow ( 1.0 - x, c - a - b ) * hf;
    return hf;
  }

  aa = a;
  bb = b;
  x1 = x;

  if ( x < 0.0 )
  {
    x = x / ( x - 1.0 );
    if ( a < c && b < a && 0.0 < b )
    {
      a = bb;
      b = aa;
    }
    b = c - b;
  }

  if ( 0.75 <= x )
  {
    gm = 0.0;

    if ( r8_abs ( c - a - b - ( int ) ( c - a - b ) ) < 1.0E-15 )
    {
      m = int ( c - a - b );
      ga = gamma ( a );
      gb = gamma ( b );
      gc = gamma ( c );
      gam = gamma ( a + m );
      gbm = gamma ( b + m );

      pa = r8_psi ( a );
      pb = r8_psi ( b );

      if ( m != 0 )
      {
        gm = 1.0;
      }

      for ( j = 1; j <= abs ( m ) - 1; j++ )
      {
        gm = gm * j;
      }

      rm = 1.0;
      for ( j = 1; j <= abs ( m ); j++ )
      {
        rm = rm * j;
      }

      f0 = 1.0;
      r0 = 1.0;;
      r1 = 1.0;
      sp0 = 0.0;;
      sp = 0.0;

      if ( 0 <= m )
      {
        c0 = gm * gc / ( gam * gbm );
        c1 = - gc * pow ( x - 1.0, m ) / ( ga * gb * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( a + k - 1.0 ) + 1.0 / ( b + k - 1.0 ) 
          - 1.0 / ( double ) ( k );
        }

        f1 = pa + pb + sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + ( 1.0 - a ) 
              / ( ( j + k ) * ( a + j + k - 1.0 ) ) 
              + 1.0 / ( b + j + k - 1.0 );
          }

          rp = pa + pb + 2.0 * el + sp + sm + log ( 1.0 - x );

          r1 = r1 * ( a + m + k - 1.0 ) * ( b + m + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }
          hw = f1;
        }
        hf = f0 * c0 + f1 * c1;
      }
      else if ( m < 0 )
      {
        m = - m;
        c0 = gm * gc / ( ga * gb * pow ( 1.0 - x, m ) );
        c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a - m + k - 1.0 ) * ( b - m + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( double ) ( k );
        }

        f1 = pa + pb - sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) 
            / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + 1.0 / ( double ) ( j + k );
          }

          rp = pa + pb + 2.0 * el + sp - sm + log ( 1.0 - x );

          r1 = r1 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }

          hw = f1;
        }

        hf = f0 * c0 + f1 * c1;
      }
    }
    else
    {
      ga = gamma ( a );
      gb = gamma ( b );
      gc = gamma ( c );
      gca = gamma ( c - a );
      gcb = gamma ( c - b );
      gcab = gamma ( c - a - b );
      gabc = gamma ( a + b - c );
      c0 = gc * gcab / ( gca * gcb );
      c1 = gc * gabc / ( ga * gb ) * pow ( 1.0 - x, c - a - b );
      hf = 0.0;
      hw = hf;
      r0 = c0;
      r1 = c1;

      for ( k = 1; k <= 250; k++ )
      {
        r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
          / ( k * ( a + b - c + k ) ) * ( 1.0 - x );

        r1 = r1 * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
          / ( k * ( c - a - b + k ) ) * ( 1.0 - x );

        hf = hf + r0 + r1;

        if ( r8_abs ( hf - hw ) < r8_abs ( hf ) * eps )
        {
          break;
        }
        hw = hf;
      }
      hf = hf + c0 + c1;
    }
  }
  else
  {
    a0 = 1.0;

    if ( a < c && c < 2.0 * a && b < c && c < 2.0 * b )
    {
      a0 = pow ( 1.0 - x, c - a - b );
      a = c - a;
      b = c - b;
    }

    hf = 1.0;
    hw = hf;
    r = 1.0;

    for ( k = 1; k <= 250; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;

      hf = hf + r;

      if ( r8_abs ( hf - hw ) <= r8_abs ( hf ) * eps )
      {
        break;
      }

      hw = hf;
    }
    hf = a0 * hf;
  }

  if ( x1 < 0.0 )
  {
    x = x1;
    c0 = 1.0 / pow ( 1.0 - x, aa );
    hf = c0 * hf;
  }

  a = aa;
  b = bb;

  if ( 120 < k )
  {
    cout << "\n";
    cout << "R8_HYPER_2F1 - Warning!\n";
    cout << "  A large number of iterations were needed.\n";
    cout << "  The accuracy of the results should be checked.\n";
  }

  return hf;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_psi ( double xx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI evaluates the function Psi(X).
//
//  Discussion:
//
//    This routine evaluates the logarithmic derivative of the
//    Gamma function,
//
//      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
//             = d/dX LN ( GAMMA(X) )
//
//    for real X, where either
//
//      - XMAX1 < X < - XMIN, and X is not a negative integer,
//
//    or
//
//      XMIN < X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody, Anthony Strecok, Henry Thacher,
//    Chebyshev Approximations for the Psi Function,
//    Mathematics of Computation,
//    Volume 27, Number 121, January 1973, pages 123-127.
//
//  Parameters:
//
//    Input, double XX, the argument of the function.
//
//    Output, double R8_PSI, the value of the function.
//
{
  double aug;
  double den;
  double four = 4.0;
  double fourth = 0.25;
  double half = 0.5;
  int i;
  int n;
  int nq;
  double one = 1.0;
  double p1[9] = { 
   4.5104681245762934160E-03, 
   5.4932855833000385356, 
   3.7646693175929276856E+02, 
   7.9525490849151998065E+03, 
   7.1451595818951933210E+04, 
   3.0655976301987365674E+05, 
   6.3606997788964458797E+05, 
   5.8041312783537569993E+05, 
   1.6585695029761022321E+05 };
  double p2[7] = { 
  -2.7103228277757834192, 
  -1.5166271776896121383E+01, 
  -1.9784554148719218667E+01, 
  -8.8100958828312219821, 
  -1.4479614616899842986, 
  -7.3689600332394549911E-02, 
  -6.5135387732718171306E-21 };
  double piov4 = 0.78539816339744830962;
  double q1[8] = { 
   9.6141654774222358525E+01, 
   2.6287715790581193330E+03, 
   2.9862497022250277920E+04, 
   1.6206566091533671639E+05, 
   4.3487880712768329037E+05, 
   5.4256384537269993733E+05, 
   2.4242185002017985252E+05, 
   6.4155223783576225996E-08 };
  double q2[6] = { 
   4.4992760373789365846E+01, 
   2.0240955312679931159E+02, 
   2.4736979003315290057E+02, 
   1.0742543875702278326E+02, 
   1.7463965060678569906E+01, 
   8.8427520398873480342E-01 };
  double sgn;
  double three = 3.0;
  double upper;
  double value;
  double w;
  double x;
  double x01 = 187.0;
  double x01d = 128.0;
  double x02 = 6.9464496836234126266E-04;
  double xinf = 1.70E+38;
  double xlarge = 2.04E+15;
  double xmax1 = 3.60E+16;
  double xmin1 = 5.89E-39;
  double xsmall = 2.05E-09;
  double z;
  double zero = 0.0;

  x = xx;
  w = r8_abs ( x );
  aug = zero;
//
//  Check for valid arguments, then branch to appropriate algorithm.
//
  if ( xmax1 <= - x || w < xmin1 )
  {
    if ( zero < x )
    {
      value = - xinf;
    }
    else
    {
      value = xinf;
    }
    return value;
  }

  if ( x < half )
  {
//
//  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
//  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
//
    if ( w <= xsmall )
    {
      aug = - one / x;
    }
//
//  Argument reduction for cotangent.
//
    else
    {
      if ( x < zero )
      {
        sgn = piov4;
      }
      else
      {
        sgn = - piov4;
      }

      w = w - ( double ) ( ( int ) ( w ) );
      nq = int ( w * four );
      w = four * ( w - ( double ) ( nq ) * fourth );
//
//  W is now related to the fractional part of 4.0 * X.
//  Adjust argument to correspond to values in the first
//  quadrant and determine the sign.
//
      n = nq / 2;

      if ( n + n != nq )
      {
        w = one - w;
      }

      z = piov4 * w;

      if ( ( n % 2 ) != 0 )
      {
        sgn = - sgn;
      }
//
//  Determine the final value for  -pi * cotan(pi*x).
//
      n = ( nq + 1 ) / 2;
      if ( ( n % 2 ) == 0 )
      {
//
//  Check for singularity.
//
        if ( z == zero )
        {
          if ( zero < x )
          {
            value = -xinf;
          }
          else
          {
            value = xinf;
          }
          return value;
        }
        aug = sgn * ( four / tan ( z ) );
      }
      else
      {
        aug = sgn * ( four * tan ( z ) );
      }
    }
    x = one - x;
  }
//
//  0.5 <= X <= 3.0.
//
  if ( x <= three )
  {
    den = x;
    upper = p1[0] * x;
    for ( i = 1; i <= 7; i++ )
    {
      den = ( den + q1[i-1] ) * x;
      upper = ( upper + p1[i]) * x;
    }
    den = ( upper + p1[8] ) / ( den + q1[7] );
    x = ( x - x01 / x01d ) - x02;
    value = den * x + aug;
    return value;
  }
//
//  3.0 < X.
//
  if ( x < xlarge )
  {
    w = one / ( x * x );
    den = w;
    upper = p2[0] * w;
    for ( i = 1; i <= 5; i++ )
    {
      den = ( den + q2[i-1] ) * w;
      upper = ( upper + p2[i] ) * w;
    }
    aug = ( upper + p2[6] ) / ( den + q2[5] ) - half / x + aug;
  }

  value = aug + log ( x );

  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  } 
  else
  {
    value = 1.0;
  }
  return value;
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

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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

void r8vec_reverse ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_REVERSE reverses the elements of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of double precision values.
//
//    Input:
//
//      N = 5,
//      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 )
//
//    Output:
//
//      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, double A[N], the array to be reversed.
//
{
  int i;
  double temp;

  for ( i = 0; i < n/2; i++ )
  {
    temp     = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = temp;
  }
  return;
}
//****************************************************************************80

void radau_compute ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RADAU_COMPUTE computes a Radau quadrature rule.
//
//  Discussion:
//
//    The Radau rule is distinguished by the fact that the left endpoint
//    (-1) is always an abscissa.
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//    The quadrature rule will integrate exactly all polynomials up to
//    X**(2*NORDER-2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2007
//
//  Author:
//
//    Original MATLAB version by Greg von Winckel.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
//    Spectral Methods in Fluid Dynamics,
//    Springer, 1993,
//    ISNB13: 978-3540522058,
//    LC: QA377.S676.
//
//    Francis Hildebrand,
//    Section 8.11,
//    Introduction to Numerical Analysis,
//    Dover, 1987,
//    ISBN13: 978-0486653631,
//    LC: QA300.H5.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.  N must be at least 1.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  int iterate;
  int iterate_max = 25;
  int j;
  double *p;
  double pi = 3.141592653589793;
  double temp;
  double test;
  double tolerance;
  double *xold;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "RADAU_COMPUTE - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << n << "\n";
    cerr << "  ORDER must be at least 1.\n";
    exit ( 1 );
  }

  tolerance = 100.0 * r8_epsilon ( );
//
//  Initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = - cos ( 2.0 * pi * ( double ) (         i ) 
                            / ( double ) ( 2 * n - 1 ) );
  }
  xold = new double[n];
  p = new double[n*(n+1)];
  iterate = 0;

  do
  {
    for ( i = 0; i < n; i++ )
    {
      xold[i] = x[i];
    }

    temp = 1.0;
    for ( j = 0; j < n + 1; j++ )
    {
      p[0+j*n] = temp;
      temp = -temp;
    }

    for ( i = 1; i < n; i++ )
    {
      p[i+0*n] = 1.0;
    }
    for ( i = 1; i < n; i++ )
    {
      p[i+1*n] = x[i];
    }

    for ( j = 2; j <= n; j++ )
    {
      for ( i = 1; i < n; i++ )
      {
        p[i+j*n] = ( ( double ) ( 2 * j - 1 ) * x[i] * p[i+(j-1)*n]     
                   + ( double ) (   - j + 1 ) *        p[i+(j-2)*n] ) 
                   / ( double ) (     j     );
      }
    }
    for ( i = 1; i < n; i++ )
    {
      x[i] = xold[i] - ( ( 1.0 - xold[i] ) / ( double ) ( n ) ) 
        * ( p[i+(n-1)*n] + p[i+n*n] ) / ( p[i+(n-1)*n] - p[i+n*n] );
    }
    test = 0.0;
    for ( i = 0; i < n; i++ )
    {
      test = r8_max ( test, r8_abs ( x[i] - xold[i] ) );
    }
    iterate = iterate + 1;
  } while ( tolerance < test && iterate < iterate_max );

  w[0] = 2.0 / ( double ) ( n * n );
  for ( i = 1; i < n; i++ )
  {
    w[i] = ( 1.0 - x[i] ) / pow ( ( double ) ( n ) * p[i+(n-1)*n], 2 );
  }
  delete [] xold;
  delete [] p;

  return;
}
//****************************************************************************80

void radau_set ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    RADAU_SET sets abscissas and weights for Radau quadrature.
//
//  Discussion:
//
//    The Radau rule is distinguished by the fact that the left endpoint
//    (-1) is always an abscissa.
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
//
//    The quadrature rule will integrate exactly all polynomials up to
//    X**(2*ORDER-2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be between 1 and 15.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  if ( order == 1 )
  {
    xtab[0] =   - 1.0E+00;
    weight[0] =   2.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =    1.0E+00 / 3.0E+00;

    weight[0] =  0.5E+00;
    weight[1] =  1.5E+00;
  }
  else if ( order == 3 )
  {
    xtab[0] =   - 1.0E+00;
    xtab[1] =   - 0.289897948556635619639456814941E+00;
    xtab[2] =     0.689897948556635619639456814941E+00;

    weight[0] =  0.222222222222222222222222222222E+00;
    weight[1] =  0.102497165237684322767762689304E+01;
    weight[2] =  0.752806125400934550100150884739E+00;
  }
  else if ( order == 4 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.575318923521694112050483779752E+00;
    xtab[2] =    0.181066271118530578270147495862E+00;
    xtab[3] =    0.822824080974592105208907712461E+00;

    weight[0] =  0.125E+00;
    weight[1] =  0.657688639960119487888578442146E+00;
    weight[2] =  0.776386937686343761560464613780E+00;
    weight[3] =  0.440924422353536750550956944074E+00;
  }
  else if ( order == 5 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.720480271312438895695825837750E+00;
    xtab[2] =  - 0.167180864737833640113395337326E+00;
    xtab[3] =    0.446313972723752344639908004629E+00;
    xtab[4] =    0.885791607770964635613757614892E+00;

    weight[0] =  0.08E+00;
    weight[1] =  0.446207802167141488805120436457E+00;
    weight[2] =  0.623653045951482508163709823153E+00;
    weight[3] =  0.562712030298924120384345300681E+00;
    weight[4] =  0.287427121582451882646824439708E+00;
  }
  else if ( order == 6 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.802929828402347147753002204224E+00;
    xtab[2] =  - 0.390928546707272189029229647442E+00;
    xtab[3] =    0.124050379505227711989974959990E+00;
    xtab[4] =    0.603973164252783654928415726409E+00;
    xtab[5] =    0.920380285897062515318386619813E+00;

    weight[0] =  0.555555555555555555555555555556E-01;
    weight[1] =  0.319640753220510966545779983796E+00;
    weight[2] =  0.485387188468969916159827915587E+00;
    weight[3] =  0.520926783189574982570229406570E+00;
    weight[4] =  0.416901334311907738959406382743E+00;
    weight[5] =  0.201588385253480840209200755749E+00;
  }
  else if ( order == 7 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.853891342639482229703747931639E+00;
    xtab[2] =  - 0.538467724060109001833766720231E+00;
    xtab[3] =  - 0.117343037543100264162786683611E+00;
    xtab[4] =    0.326030619437691401805894055838E+00;
    xtab[5] =    0.703842800663031416300046295008E+00;
    xtab[6] =    0.941367145680430216055899446174E+00;

    weight[0] =  0.408163265306122448979591836735E-01;
    weight[1] =  0.239227489225312405787077480770E+00;
    weight[2] =  0.380949873644231153805938347876E+00;
    weight[3] =  0.447109829014566469499348953642E+00;
    weight[4] =  0.424703779005955608398308039150E+00;
    weight[5] =  0.318204231467301481744870434470E+00;
    weight[6] =  0.148988471112020635866497560418E+00;
  }
  else if ( order == 8 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.887474878926155707068695617935E+00;
    xtab[2] =  - 0.639518616526215270024840114382E+00;
    xtab[3] =  - 0.294750565773660725252184459658E+00;
    xtab[4] =    0.943072526611107660028971153047E-01;
    xtab[5] =    0.468420354430821063046421216613E+00;
    xtab[6] =    0.770641893678191536180719525865E+00;
    xtab[7] =    0.955041227122575003782349000858E+00;

    weight[0] =  0.03125E+00;
    weight[1] =  0.185358154802979278540728972699E+00;
    weight[2] =  0.304130620646785128975743291400E+00;
    weight[3] =  0.376517545389118556572129261442E+00;
    weight[4] =  0.391572167452493593082499534004E+00;
    weight[5] =  0.347014795634501280228675918422E+00;
    weight[6] =  0.249647901329864963257869293513E+00;
    weight[7] =  0.114508814744257199342353728520E+00;
  }
  else if ( order == 9 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.910732089420060298533757956283E+00;
    xtab[2] =  - 0.711267485915708857029562959544E+00;
    xtab[3] =  - 0.426350485711138962102627520502E+00;
    xtab[4] =  - 0.903733696068532980645444599064E-01;
    xtab[5] =    0.256135670833455395138292079035E+00;
    xtab[6] =    0.571383041208738483284917464837E+00;
    xtab[7] =    0.817352784200412087992517083851E+00;
    xtab[8] =    0.964440169705273096373589797925E+00;

    weight[0] =  0.246913580246913580246913580247E-01;
    weight[1] =  0.147654019046315385819588499802E+00;
    weight[2] =  0.247189378204593052361239794969E+00;
    weight[3] =  0.316843775670437978338000849642E+00;
    weight[4] =  0.348273002772966594071991031186E+00;
    weight[5] =  0.337693966975929585803724239792E+00;
    weight[6] =  0.286386696357231171146705637752E+00;
    weight[7] =  0.200553298024551957421165090417E+00;
    weight[8] =  0.907145049232829170128934984159E-01;
  }
  else if ( order == 10 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.927484374233581078117671398464E+00;
    xtab[2] =  - 0.763842042420002599615429776011E+00;
    xtab[3] =  - 0.525646030370079229365386614293E+00;
    xtab[4] =  - 0.236234469390588049278459503207E+00;
    xtab[5] =    0.760591978379781302337137826389E-01;
    xtab[6] =    0.380664840144724365880759065541E+00;
    xtab[7] =    0.647766687674009436273648507855E+00;
    xtab[8] =    0.851225220581607910728163628088E+00;
    xtab[9] =   0.971175180702246902734346518378E+00;

    weight[0] =  0.02E+00;
    weight[1] =  0.120296670557481631517310522702E+00;
    weight[2] =  0.204270131879000675555788672223E+00;
    weight[3] =  0.268194837841178696058554475262E+00;
    weight[4] =  0.305859287724422621016275475401E+00;
    weight[5] =  0.313582457226938376695902847302E+00;
    weight[6] =  0.290610164832918311146863077963E+00;
    weight[7] =  0.239193431714379713376571966160E+00;
    weight[8] =  0.164376012736921475701681668908E+00;
    weight[9] = 0.736170054867584989310512940790E-01;
  }
  else if ( order == 11 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.939941935677027005913871284731E+00;
    xtab[2] =  - 0.803421975580293540697597956820E+00;
    xtab[3] =  - 0.601957842073797690275892603234E+00;
    xtab[4] =  - 0.351888923353330214714301017870E+00;
    xtab[5] =  - 0.734775314313212657461903554238E-01;
    xtab[6] =    0.210720306228426314076095789845E+00;
    xtab[7] =    0.477680647983087519467896683890E+00;
    xtab[8] =    0.705777100713859519144801128840E+00;
    xtab[9] =   0.876535856245703748954741265611E+00;
    xtab[10] =   0.976164773135168806180508826082E+00;

    weight[0] =  0.165289256198347107438016528926E-01;
    weight[1] =  0.998460819079680638957534695802E-01;
    weight[2] =  0.171317619206659836486712649042E+00;
    weight[3] =  0.228866123848976624401683231126E+00;
    weight[4] =  0.267867086189684177806638163355E+00;
    weight[5] =  0.285165563941007337460004408915E+00;
    weight[6] =  0.279361333103383045188962195720E+00;
    weight[7] =  0.250925377697128394649140267633E+00;
    weight[8] =  0.202163108540024418349931754266E+00;
    weight[9] = 0.137033682133202256310153880580E+00;
    weight[10] = 0.609250978121311347072183268883E-01;
  }
  else if ( order == 12 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.949452759204959300493337627077E+00;
    xtab[2] =  - 0.833916773105189706586269254036E+00;
    xtab[3] =  - 0.661649799245637148061133087811E+00;
    xtab[4] =  - 0.444406569781935851126642615609E+00;
    xtab[5] =  - 0.196994559534278366455441427346E+00;
    xtab[6] =    0.637247738208319158337792384845E-01;
    xtab[7] =    0.319983684170669623532789532206E+00;
    xtab[8] =    0.554318785912324288984337093085E+00;
    xtab[9] =   0.750761549711113852529400825472E+00;
    xtab[10] =   0.895929097745638894832914608454E+00;
    xtab[11] =   0.979963439076639188313950540264E+00;

    weight[0] =  0.138888888888888888888888888888E-01;
    weight[1] =  0.841721349386809762415796536813E-01;
    weight[2] =  0.145563668853995128522547654706E+00;
    weight[3] =  0.196998534826089634656049637969E+00;
    weight[4] =  0.235003115144985839348633985940E+00;
    weight[5] =  0.256991338152707776127974253598E+00;
    weight[6] =  0.261465660552133103438074715743E+00;
    weight[7] =  0.248121560804009959403073107079E+00;
    weight[8] =  0.217868879026192438848747482023E+00;
    weight[9] = 0.172770639313308564306065766966E+00;
    weight[10] = 0.115907480291738392750341908272E+00;
    weight[11] = 0.512480992072692974680229451351E-01;
  }
  else if ( order == 13 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.956875873668299278183813833834E+00;
    xtab[2] =  - 0.857884202528822035697620310269E+00;
    xtab[3] =  - 0.709105087529871761580423832811E+00;
    xtab[4] =  - 0.519197779050454107485205148087E+00;
    xtab[5] =  - 0.299201300554509985532583446686E+00;
    xtab[6] =  - 0.619016986256353412578604857936E-01;
    xtab[7] =    0.178909837597084635021931298881E+00;
    xtab[8] =    0.409238231474839556754166331248E+00;
    xtab[9] =   0.615697890940291918017885487543E+00;
    xtab[10] =   0.786291018233046684731786459135E+00;
    xtab[11] =   0.911107073689184553949066402429E+00;
    xtab[12] =   0.982921890023145161262671078244E+00;

    weight[0] =  0.118343195266272189349112426036E-01;
    weight[1] =  0.719024162924955289397537405641E-01;
    weight[2] =  0.125103834331152358133769287976E+00;
    weight[3] =  0.171003460470616642463758674512E+00;
    weight[4] =  0.206960611455877074631132560829E+00;
    weight[5] =  0.230888862886995434012203758668E+00;
    weight[6] =  0.241398342287691148630866924129E+00;
    weight[7] =  0.237878547660712031342685189180E+00;
    weight[8] =  0.220534229288451464691077164199E+00;
    weight[9] = 0.190373715559631732254759820746E+00;
    weight[10] = 0.149150950090000205151491864242E+00;
    weight[11] = 0.992678068818470859847363877478E-01;
    weight[12] = 0.437029032679020748288533846051E-01;
  }
  else if ( order == 14 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.962779269978024297120561244319E+00;
    xtab[2] =  - 0.877048918201462024795266773531E+00;
    xtab[3] =  - 0.747389642613378838735429134263E+00;
    xtab[4] =  - 0.580314056546874971105726664999E+00;
    xtab[5] =  - 0.384202003439203313794083903375E+00;
    xtab[6] =  - 0.168887928042680911008441695622E+00;
    xtab[7] =    0.548312279917645496498107146428E-01;
    xtab[8] =    0.275737205435522399182637403545E+00;
    xtab[9] =   0.482752918588474966820418534355E+00;
    xtab[10] =   0.665497977216884537008955042481E+00;
    xtab[11] =   0.814809550601994729434217249123E+00;
    xtab[12] =   0.923203722520643299246334950272E+00;
    xtab[13] =   0.985270697947821356698617003172E+00;

    weight[0] =  0.102040816326530612244897959184E-01;
    weight[1] =  0.621220169077714601661329164668E-01;
    weight[2] =  0.108607722744362826826720935229E+00;
    weight[3] =  0.149620539353121355950520836946E+00;
    weight[4] =  0.183127002125729654123867302103E+00;
    weight[5] =  0.207449763335175672668082886489E+00;
    weight[6] =  0.221369811499570948931671683021E+00;
    weight[7] =  0.224189348002707794238414632220E+00;
    weight[8] =  0.215767100604618851381187446115E+00;
    weight[9] = 0.196525518452982430324613091930E+00;
    weight[10] = 0.167429727891086278990102277038E+00;
    weight[11] = 0.129939668737342347807425737146E+00;
    weight[12] = 0.859405354429804030893077310866E-01;
    weight[13] = 0.377071632698969142774627282919E-01;
  }
  else if ( order == 15 )
  {
    xtab[0] =  - 1.0E+00;
    xtab[1] =  - 0.967550468197200476562456018282E+00;
    xtab[2] =  - 0.892605400120550767066811886849E+00;
    xtab[3] =  - 0.778685617639031079381743321893E+00;
    xtab[4] =  - 0.630779478886949283946148437224E+00;
    xtab[5] =  - 0.455352905778529370872053455981E+00;
    xtab[6] =  - 0.260073376740807915768961188263E+00;
    xtab[7] =  - 0.534757226797460641074538896258E-01;
    xtab[8] =    0.155410685384859484319182024964E+00;
    xtab[9] =   0.357456512022127651195319205174E+00;
    xtab[10] =   0.543831458701484016930711802760E+00;
    xtab[11] =   0.706390264637572540152679669478E+00;
    xtab[12] =   0.838029000636089631215097384520E+00;
    xtab[13] =   0.932997190935973719928072142859E+00;
    xtab[14] =   0.987166478414363086378359071811E+00;

    weight[0] =  0.888888888888888888888888888889E-02;
    weight[1] =  0.542027800486444943382142368018E-01;
    weight[2] =  0.951295994604808992038477266346E-01;
    weight[3] =  0.131875462504951632186262157944E+00;
    weight[4] =  0.162854477303832629448732245828E+00;
    weight[5] =  0.186715145839450908083795103799E+00;
    weight[6] =  0.202415187030618429872703310435E+00;
    weight[7] =  0.209268608147694581430889790306E+00;
    weight[8] =  0.206975960249553755479027321787E+00;
    weight[9] = 0.195637503045116116473556617575E+00;
    weight[10] = 0.175748872642447685670310440476E+00;
    weight[11] = 0.148179527003467253924682058743E+00;
    weight[12] = 0.114135203489752753013075582569E+00;
    weight[13] = 0.751083927605064397329716653914E-01;
    weight[14] = 0.328643915845935322530428528231E-01;
  }
  else
  {
    cerr << "\n";
    cerr << "RADAU_SET - Fatal error!\n";
    cerr << "  Illegal value of ORDER = " << order << "\n";
    cerr << "  Legal values are 1 to 15.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void rule_adjust ( double a, double b, double c, double d, int order, 
  double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
//
//  Discussion:
//
//    Most quadrature rules are defined on a special interval, like
//    [-1,1] or [0,1].  To integrate over an interval, the abscissas
//    and weights must be adjusted.  This can be done on the fly,
//    or by calling this routine.
//
//    If the weight function W(X) is not 1, then the weight vector W will
//    require further adjustment by the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the definition interval.
//
//    Input, double C, D, the endpoints of the integration interval.
//
//    Input, int ORDER, the number of abscissas and weights.
//
//    Input/output, double X[ORDER], W[ORDER], the abscissas
//    and weights.
//
{
  int i;

  for ( i = 0; i < order; i++ )
  {
    x[i] = ( ( b - x[i]     ) * c   
           + (     x[i] - a ) * d ) 
           / ( b        - a );
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = ( ( d - c ) / ( b - a ) ) * w[i];
  }

  return;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal.
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ )
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length )
  {
    for ( i = nchar; i < s1_length; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return false;
      }
    }
  }
  else if ( nchar < s2_length )
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' )
      {
        return false;
      }
    }
  }

  return true;
}
//****************************************************************************80

double sum_sub ( double func ( double x ), double a, double b, int nsub, 
  int order, double xlo, double xhi, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SUM_SUB carries out a composite quadrature rule.
//
//  Discussion:
//
//    SUM_SUB assumes the original rule was written for [XLO,XHI].
//
//    The integral:
//
//      Integral ( A <= X <= B ) F(X) dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double X ), the name of the function which
//    evaluates the integrand.
//
//    Input, double A, B, the lower and upper limits of integration.
//
//    Input, int NSUB, the number of equal subintervals into
//    which the finite interval (A,B) is to be subdivided for
//    higher accuracy.  NSUB must be at least 1.
//
//    Input, int ORDER, the order.
//    ORDER must be at least 1.
//
//    Input, double XLO, XHI, the left and right endpoints of
//    the interval over which the quadrature rule was defined.
//
//    Input, double XTAB[ORDER], the abscissas of a quadrature
//    rule for the interval [XLO,XHI].
//
//    Input, double WEIGHT[ORDER], the weights of the
//    quadrature rule.
//
//    Output, double SUM_SUB, the approximate value of the integral.
//
{
  double a_sub;
  double b_sub;
  int i;
  int j;
  double quad_sub;
  double result;
  double result_sub;
  double x;
  double volume;
  double volume_sub;

  if ( a == b )
  {
    result = 0.0;
    return result;
  }

  if ( order < 1 )
  {
    cerr << "\n";
    cerr << "SUM_SUB - Fatal error!\n";
    cerr << "  Nonpositive value of ORDER = " << order << "\n";
    exit ( 1 );
  }

  if ( nsub < 1 )
  {
    cerr << "\n";
    cerr << "SUM_SUB - Fatal error!\n";
    cerr << "  Nonpositive value of NSUB = " << nsub << "\n";
    exit ( 1 );
  }

  if ( xlo == xhi )
  {
    cerr << "\n";
    cerr << "SUM_SUB - Fatal error!\n";
    cerr << "  XLO = XHI.\n";
    exit ( 1 );
  }

  volume = 0.0;
  result = 0.0;

  for ( j = 1; j <= nsub; j++ )
  {
    a_sub = ( ( double ) ( nsub - j + 1 ) * a   
            + ( double ) (        j - 1 ) * b )
            / ( double ) ( nsub         );

    b_sub = ( ( double ) ( nsub - j )     * a   
            + ( double ) (        j )     * b ) 
            / ( double ) ( nsub     );

    quad_sub = 0.0;
    for ( i = 0; i < order; i++ )
    {
      x = ( ( xhi - xtab[i]       ) * a_sub   
          + (       xtab[i] - xlo ) * b_sub ) 
          / ( xhi           - xlo );
      quad_sub = quad_sub + weight[i] * func ( x );
    }

    volume_sub = ( b - a ) / ( ( xhi - xlo ) * ( double ) ( nsub ) );
    result_sub = quad_sub * volume_sub;

    volume = volume + volume_sub;
    result = result + result_sub;
  }

  return result;
}
//****************************************************************************80

double summer ( double func ( double x ), int order, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    SUMMER carries out a quadrature rule over a single interval.
//
//  Formula:
//
//    RESULT = sum ( 1 <= I <= ORDER ) W(I) * FUNC ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double X ), the name of the function which
//    evaluates the integrand.
//
//    Input, int ORDER, the order.
//
//    Input, double X[ORDER], the abscissas.
//
//    Input, double W[ORDER], the weights.
//
//    Output, double SUMMER, the approximate value of the integral.
//
{
  int i;
  double result;

  result = 0.0;
  for ( i = 0; i < order; i++ )
  {
    result = result + w[i] * func ( x[i] );
  }

  return result;
}
//****************************************************************************80

void summer_gk ( double func ( double x ), int orderg, double weightg[],
  double *resultg, int orderk, double xtabk[], double weightk[], 
  double *resultk )

//****************************************************************************80
//
//  Purpose:
//
//    SUMMER_GK carries out Gauss-Kronrod quadrature over a single interval.
//
//  Discussion:
//
//    The abscissas for the Gauss-Legendre rule of order ORDERG are
//    not required, since they are assumed to be the even-indexed
//    entries of the corresponding Kronrod rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x ), the name of the function which
//    evaluates the integrand.
//
//    Input, int ORDERG, the order of the Gauss-Legendre rule.
//
//    Input, double WEIGHTG[ORDERG], the weights of the
//    Gauss-Legendre rule.
//
//    Output, double *RESULTG, the approximate value of the
//    integral, based on the Gauss-Legendre rule.
//
//    Input, int ORDERK, the order of the Kronrod rule.  ORDERK
//    must equal 2 * ORDERG + 1.
//
//    Input, double XTABK[ORDERK], the abscissas of the Kronrod rule.
//
//    Input, double WEIGHTK[ORDERK], the weights of the Kronrod rule.
//
//    Output, double *RESULTK, the approximate value of the integral,
//    based on the Kronrod rule.
//
{
  double fk;
  int i;

  if ( orderk != 2 * orderg + 1 )
  {
    cerr << "\n";
    cerr << "SUMMER_GK - Fatal error!\n";
    cerr << "  ORDERK must equal 2 * ORDERG + 1.\n";
    cerr << "  The input value was ORDERG = " << orderg << "\n";
    cerr << "  The input value was ORDERK = " << orderk << "\n";
    exit ( 1 );
  }

  *resultg = 0.0;
  *resultk = 0.0;

  for ( i = 0; i < orderk; i++ )
  {
    fk = func ( xtabk[i] );

    *resultk = *resultk + weightk[i] * fk;

    if ( ( i % 2 ) == 1 )
    {
      *resultg = *resultg + weightg[(i-1)/2] * fk;
    }
  }
  return;
}
//****************************************************************************80

void sum_sub_gk ( double func ( double x ), double a, double b, int nsub, 
  int orderg, double weightg[], double *resultg, int orderk, double xtabk[], 
  double weightk[], double *resultk, double *error )

//****************************************************************************80
//
//  Purpose:
//
//    SUM_SUB_GK carries out a composite Gauss-Kronrod rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( A <= X <= B ) F(X) dX
//
//    The quadrature rule:
//
//      H = ( B - A ) / NSUB
//      XMID(J) = A + 0.5 * H * ( 2 * J - 1 )
//
//      Sum ( 1 <= J <= NSUB )
//        Sum ( 1 <= I <= ORDERK )
//          WEIGHTK(I) * F ( XMID(J) + 0.5 * H * XTABK(I) )
//
//    The Gauss-Legendre weights should be computed by LEGCOM or LEGSET.
//    The Kronrod abscissas and weights should be computed by KRONSET.
//
//    The orders of the Gauss-Legendre and Kronrod rules must satisfy
//    ORDERK = 2 * ORDERG + 1.
//
//    The Kronrod rule uses the abscissas of the Gauss-Legendre rule,
//    plus more points, resulting in an efficient and higher order estimate.
//
//    The difference between the Gauss-Legendre and Kronrod estimates
//    is taken as an estimate of the error in the approximation to the
//    integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double X ), the name of the function which
//    evaluates the integrand.  
//
//    Input, double A, B, the lower and upper limits of integration.
//
//    Input, int NSUB, the number of equal subintervals into
//    which the finite interval (A,B) is to be subdivided for
//    higher accuracy.  NSUB must be at least 1.
//
//    Input, int ORDERG, the order of the Gauss-Legendre rule.
//    ORDERG must be at least 1.
//
//    Input, double WEIGHTG[ORDERG], the weights of the
//    Gauss-Legendre rule.
//
//    Output, double *RESULTG, the approximate value of the
//    integral based on the Gauss-Legendre rule.
//
//    Input, int ORDERK, the order of the Kronrod rule.
//    ORDERK must be at least 1.
//
//    Input, double XTABK[ORDERK], the abscissas of the
//    Kronrod rule.
//
//    Input, double WEIGHTK[ORDERK], the weights of the
//    Kronrod rule.
//
//    Output, double *RESULTK, the approximate value of the
//    integral based on the Kronrod rule.
//
//    Output, double *ERROR, an estimate of the approximation
//    error.  This is computed by taking the sum of the absolute values of
//    the differences between the Gauss-Legendre and Kronrod rules
//    over each subinterval.  This is usually a good estimate of
//    the error in the value RESULTG.  The error in the Kronrod
//    estimate RESULTK is usually much smaller.
//
{
  double fk;
  double h;
  int i;
  int j;
  double partg;
  double partk;
  double x;
  double xmid;

  *resultg = 0.0;
  *resultk = 0.0;
  *error = 0.0;

  if ( a == b )
  {
    return;
  }

  if ( orderk != 2 * orderg + 1 )
  {
    cerr << "\n";
    cerr << "SUM_SUB_GK - Fatal error!\n";
    cerr << "  ORDERK must equal 2 * ORDERG + 1.\n";
    cerr << "  The input value was ORDERG = " << orderg << "\n";
    cerr << "  The input value was ORDERK = " << orderk << "\n";
    exit ( 1 );
  }

  h = ( b - a ) / ( double ) ( nsub );

  for ( j = 0; j < nsub; j++ )
  {
    xmid = a + 0.5 * h * ( double ) ( 2 * j + 1 );

    partg = 0.0;
    partk = 0.0;

    for ( i = 0; i < orderk; i++ )
    {
      x = xmid + 0.5 * h * xtabk[i];
      fk = func ( x );
      partk = partk + 0.5 * h * weightk[i] * fk;

      if ( ( i % 2 ) == 1 )
      {
        partg = partg + 0.5 * h * weightg[(i-1)/2] * fk;
      }
    }
    *resultg = *resultg + partg;
    *resultk = *resultk + partk;
    *error = *error + r8_abs ( partk - partg );
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
