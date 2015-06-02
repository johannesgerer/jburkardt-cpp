# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "fem1d_lagrange.hpp"

//****************************************************************************80

void fem1d_lagrange_stiffness ( int x_num, double x[], int q_num, 
  double f ( double x ), double a[], double m[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_LAGRANGE_STIFFNESS evaluates the Lagrange polynomial stiffness matrix.
//
//  Discussion:
//
//    The finite element method is to be applied over a given interval that
//    has been meshed with X_NUM points X.
//
//    The finite element basis functions are to be the X_NUM Lagrange
//    basis polynomials L(i)(X), such that
//      L(i)(X(j)) = delta(i,j).
//
//    The following items are computed:
//    * A, the stiffness matrix, with A(I,J) = integral L'(i)(x) L'(j)(x)
//    * M, the mass matrix, with M(I,J) = integral L(i)(x) L(j)(x)
//    * B, the load matrix, with B(I) = integral L(i)(x) F(x)
//
//    The integrals are approximated by quadrature.
//
//    Boundary conditions are not handled by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X[X_NUM], the coordinates of the nodes.
//
//    Input, int Q_NUM, the number of quadrature points to use.
//
//    Input, double F ( double X ), the right hand side function.
//
//    Output, double A[X_NUM*X_NUM], the stiffness matrix.
//
//    Output, double M[X_NUM*X_NUM], the mass matrix.
//
//    Output, double B[X_NUM], the right hand side vector.
//
{
  double *l;
  double li;
  double lj;
  double *lp;
  double lpi;
  double lpj;
  int q_i;
  double *q_w;
  double *q_x;
  int x_i;
  int x_j;
//
//  Get the quadrature rule for [-1,+1].
//
  q_x = new double[q_num];
  q_w = new double[q_num];

  legendre_set ( q_num, q_x, q_w );
//
//  Adjust the quadrature rule to the interval [ x(1), x(x_num) }
//
  for ( q_i = 0; q_i < q_num; q_i++ )
  {
    q_x[q_i] =  ( ( 1.0 - q_x[q_i] ) * x[0] 
                + ( 1.0 + q_x[q_i] ) * x[x_num-1] ) 
                /   2.0;

    q_w[q_i] = q_w[q_i] * ( x[x_num-1] - x[0] ) / 2.0;
  }
//
//  Evaluate all the Lagrange basis polynomials at all the quadrature points.
//
  l = lagrange_value ( x_num, x, q_num, q_x );
//
//  Evaluate all the Lagrange basis polynomial derivatives at all 
//  the quadrature points.
//
  lp = lagrange_derivative ( x_num, x, q_num, q_x );
//
//  Assemble the matrix and right hand side.
//
  for ( x_j = 0; x_j < x_num; x_j++ )
  {
    for ( x_i = 0; x_i < x_num; x_i++ )
    {
      a[x_i+x_j*x_num] = 0.0;
      m[x_i+x_j*x_num] = 0.0;
    }
    b[x_j] = 0.0;
  }

  for ( x_i = 0; x_i < x_num; x_i++ )
  {
    for ( q_i = 0; q_i < q_num; q_i++ )
    {
      li = l[q_i+x_i*q_num];
      lpi = lp[q_i+x_i*q_num];
      for ( x_j = 0; x_j < x_num; x_j++ )
      {
        lj = l[q_i+x_j*q_num];
        lpj = lp[q_i+x_j*q_num];
        a[x_i+x_j*x_num] = a[x_i+x_j*x_num] + q_w[q_i] * lpi * lpj;
        m[x_i+x_j*x_num] = m[x_i+x_j*x_num] + q_w[q_i] * li * lj;
      }
      b[x_i] = b[x_i] + q_w[q_i] * li * f ( q_x[q_i] );
    }
  }

  delete [] l;
  delete [] lp;
  delete [] q_w;
  delete [] q_x;

  return;
}
//****************************************************************************80

double *lagrange_derivative ( int nd, double xd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_DERIVATIVE evaluates the Lagrange basis derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//    ND must be at least 1.
//
//    Input, double XD[ND], the data points.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double LAGRANGE_DERIVATIVE[NI*ND], the values
//    of the Lagrange basis derivatives at the interpolation points.
//
{
  int i;
  int j;
  int j1;
  int j2;
  double *lpi;
  double p;

  lpi = new double[ni*nd];

  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < ni; i++ )
    {
      lpi[i+j*ni] = 0.0;
    }
  }

  for ( i = 0; i < ni; i++ )
  {
    for ( j = 0; j < nd; j++ )
    {
      for ( j1 = 0; j1 < nd; j1++ )
      {
        if ( j1 != j )
        {
          p = 1.0;
          for ( j2 = 0; j2 < nd; j2++ )
          {
            if ( j2 == j1 )
            {
              p = p / ( xd[j] - xd[j2] );
            }
            else if ( j2 != j )
            {
              p = p * ( xi[i] - xd[j2] ) / ( xd[j] - xd[j2] );
            }
          }
          lpi[i+j*ni] = lpi[i+j*ni] + p;
        }
      }
    }
  }

  return lpi;
}
//****************************************************************************80

double *lagrange_value ( int nd, double xd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_VALUE evaluates the Lagrange basis polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//    ND must be at least 1.
//
//    Input, double XD[ND], the data points.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double LAGRANGE_BASIS[NI*ND], the values
//    of the Lagrange basis polynomials at the interpolation points.
//
{
  int i;
  int j;
  int k;
  double *li;
//
//  Evaluate the polynomial.
//
  li = new double[ni*nd];

  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < ni; i++ )
    {
      li[i+j*ni] = 1.0;
    }
  }

  for ( i = 0; i < nd; i++ )
  {
    for ( j = 0; j < nd; j++ )
    {
      if ( j != i )
      {
        for ( k = 0; k < ni; k++ )
        {
          li[k+i*ni] = li[k+i*ni] * ( xi[k] - xd[j] ) / ( xd[i] - xd[j] );
        }
      }
    }
  }

  return li;
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
//    The quadrature rule is exact for polynomials through degree 2*N-1.
//  
//    The abscissas are the zeroes of the Legendre polynomial P(ORDER)(X).
//  
//    Mathematica can compute the abscissas and weights of a Gauss-Legendre
//    rule of order N for the interval [A,B] with P digits of precision
//    by the commands:
//
//    Needs["NumericalDifferentialEquationAnalysis`"]
//    GaussianQuadratureWeights [n, a, b, p ]
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
//    N must be between 1 and 16.
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
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET - Fatal error!\n";
    cerr << "  Illegal value of N = %d\n";
    cerr << "  Legal values are 1 to 16\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double *r8mat_fs_new ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_FS_NEW factors and solves a system with one right hand side.
//
//  Discussion:
//
//    This routine differs from R8MAT_FSS_NEW in two ways:
//    * only one right hand side is allowed;
//    * the input matrix A is not modified.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double X[N], the solution of the linear system.
//
{
  double *a2;
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  a2 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = fabs ( a2[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < fabs ( a2[i-1+(jcol-1)*n] ) )
      {
        piv = fabs ( a2[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_FS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                  = a2[jcol-1+(j-1)*n];
        a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
        a2[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
//
//  Scale the pivot row.
//
    t = a2[jcol-1+(jcol-1)*n];
    a2[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a2[jcol-1+(j-1)*n] = a2[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a2[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a2[i-1+(jcol-1)*n];
        a2[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a2[i-1+(j-1)*n] = a2[i-1+(j-1)*n] + t * a2[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a2[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  delete [] a2;

  return x;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
//
//    In other words, the interval is divided into N-1 even subintervals,
//    and the endpoints of intervals are used as the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
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
