# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( void );
void adapt ( int ndim, double a[], double b[], int *minpts, int maxpts, 
  double functn ( int indx, int ndim, double z[], double alpha[], 
  double beta[] ), double rel_tol, int itest, double alpha[], double beta[],
  int lenwrk, double wrkstr[], double *relerr, double *finest, int *ifail );
double genz_function ( int indx, int ndim, double z[], double alpha[], 
  double beta[] );
double genz_integral ( int indx, int ndim, double a[], double b[], 
  double alpha[], double beta[] );
char *genz_name ( int indx );
double genz_phi ( double z );
double genz_random ( int *seed );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
int i4vec_sum ( int n, int a[] );
void multst ( int nsamp, int tstlim, int tstfns[], int tstmax, double difclt[], 
  double expnts[], int ndiml, int ndims[], char *sbname, 
  void subrtn ( int ndim, double a[], double b[], int *minpts, int maxpts, 
    double functn ( int indx, int ndim, double z[], double alpha[], 
      double beta[] ), 
    double rel_tol, int itest, double alpha[], double beta[], int lenwrk, 
    double wrkstr[], double *errest, double *finest, int *ifail ), 
  double rel_tol, int maxpts );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8vec_dot ( int n, double a1[], double a2[] );
void r8vec_median ( int n, double r[], double rmed[3] );
double r8vec_product ( int n, double a[] );
double r8vec_sum ( int n, double a[] );
void timestamp ( void );
void tuple_next ( int m1, int m2, int n, int *rank, int x[] );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TESTPACK.
//
//  Discussion:
//
//    TESTPACK is a collection of several items, including six test
//    integrand functions, an early version of ADAPT, a multidimensional
//    quadrature program, and MULTST, a routine that tests quadrature programs
//    on the test integrands.  These have all been combined to make
//    an executable program that demonstrates the testing process.
//
//  Modified:
//
//    21 March 2007
//
//  Author:
//
//    Alan Genz
//
//    C++ version by John Burkardt
//
{
# define NDIML 5
# define TSTLIM 6
# define TSTMAX 6

  double difclt[TSTMAX] = { 110.0, 600.0, 600.0, 100.0, 150.0, 100.0 };
  double expnts[TSTMAX] = { 1.5, 2.0, 2.0, 1.0, 2.0, 2.0 };
  int i;
  int maxpts = 10000;
  int ndims[NDIML] = { 2, 3, 4, 6, 8 };
  int nsamp = 20;
  double rel_tol = 1.0E-06;
  char *sbname = "ADAPT";
  int tstfns[TSTLIM] = { 1, 2, 3, 4, 5, 6 };

  timestamp ( );
  cout << "\n";
  cout << "TESTPACK\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Call MULTST, which can test a routine that\n";
  cout << "  is designed to estimate multidimensional\n";
  cout << "  integrals, by numerical quadrature.\n";
  cout << "\n";
  cout << "  The routine to be tested here is called ADAPT.\n";
  cout << "\n";
  cout << "  The test integrands are Genz's standard set.\n";
  cout << "\n";
  cout << "  MULTST, ADAPT and the test integrands were\n";
  cout << "  written in FORTRAN77 by Alan Genz.\n";

  multst ( nsamp, TSTLIM, tstfns, TSTMAX, difclt, 
    expnts, NDIML, ndims, sbname, adapt, rel_tol, maxpts );

  cout << "\n";
  cout << "TESTPACK\n";
  cout << "  Normal end of execution\n";
  cout << "\n";
  timestamp ( );

  return 0;
# undef NDIML
# undef TSTLIM
# undef TSTMAX 
}
//****************************************************************************80

void adapt ( int ndim, double a[], double b[], int *minpts, int maxpts, 
  double functn ( int indx, int ndim, double z[], double alpha[], 
  double beta[] ), double rel_tol, int itest, double alpha[], double beta[],
  int lenwrk, double wrkstr[], double *relerr, double *finest, int *ifail )

//****************************************************************************80
//
//  Purpose:
//
//    ADAPT carries out adaptive multidimensional quadrature.
//
//  Modified:
//
//    22 March 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan Genz
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the number of variables.
//    2 <= NDIM.
//
//    Input, double A[NDIM], the lower limits of integration.
//
//    Input, double B[NDIM], the upper limits of integration.
//
//    Input/output, int *MINPTS, the minimum number of function evaluations
//    to be allowed,  MINPTS must not exceed MAXPTS.  If MINPTS < 0 then the
//    routine assumes a previous call has been made with the same integrand
//    and continues that calculation.
//
//    Input, int MAXPTS, the maximum number of function
//    evaluations allowed, which must be at least RULCLS, where
//    RULCLS = 2**NDIM + 2 * NDIM**2 + 2 * NDIM + 1, when NDIM <= 15 and
//    RULCLS = ( NDIM * ( 14 - NDIM * ( 6 - 4 * NDIM ) ) ) / 3 + 1,
//    when 15 < NDIM.
//    for NDIM  =  2   3   4   5   6   7   8   9   10   11   12
//    RULCLS   =  17  33  57  93 149 241 401 693 1245 2313 4409
//    A suggested starting value for MAXPTS is 100*RULCLS.  If
//    this is not large enough for the required accuracy, then
//    MAXPTS and LENWRK should be increased accordingly.
//
//    Input, external, double FUNCTN, the user-defined function
//    to be integrated.  It must have the form
//      double functn ( int indx, ind ntim, double z[], double alpha[],
//        double beta[] )
//    where
//      INDX is the index of the test function,
//      NDIM is the spatial dimension,
//      Z is the evaluation point,
//      ALPHA is a parameter vector,
//      BETA is a parameter vector.
//
//    Input, double REL_TOL, the user's requested relative accuracy.
//
//    Input, int ITEST, the index of the test.
//
//    Input, double ALPHA[NDIM], BETA[NDIM], parameters
//    associated with the integrand function.
//
//    Input, int LENWRK, the length of the array WRKSTR.
//    The routine needs (2*NDIM+3)*(1+MAXPTS/RULCLS)/2 for LENWRK if
//    MAXPTS function calls are used.
//
//    Input/output, double WRKSTR[LENWRK].  This array does not
//    need to be set or inspected by the user.  However, the output value of
//    WKRSTR from one call may be needed by the program on a followup call
//    if the input value of MINPTS < 0, which signals that another calculation
//    is requested for the same integrand.
//
//    Output, double *RELERR, the estimated relative accuracy
//    of the integral estimate.
//
//    Output, double *FINEST, the estimated value of integral.
//
//    Output, int *IFAIL
//    * 0, for normal exit, when estimated relative error RELERR is less
//    than REL_TOL, and with MAXPTS or less function calls made.
//    * 1, if MAXPTS was too small for ADAPT to obtain the required relative
//    error REL_TOL.  In this case ADAPT returns a value of FINEST with
//    estimated relative error RELERR.
//    * 2, if LENWRK was too small for MAXPTS function calls.  In
//    this case ADAPT returns a value of FINEST with estimated error
//    RELERR using the working storage available, but RELERR is likely to
//    be greater than REL_TOL.
//    * 3, if NDIM < 2 or MAXPTS < MINPTS or MAXPTS < RULCLS.
//
{
  double *center;
  double df1;
  double df2;
  double dif;
  double difmax;
  int divaxn;
  int divaxo;
  int divflg;
  double f1;
  double f2;
  double f3;
  double f4;
  int funcls;
  int i;
  int index1;
  int index2;
  int j;
  int k;
  int l;
  double lambda2;
  double lambda4;
  double lambda5;
  int m;
  int n;
  double ratio;
  double rgncmp;
  double rgnerr;
  int rgnstr = 0;
  double rgnval;
  double rgnvol;
  int rulcls;
  int sbrgns = 0;
  int sbtmpp;
  int subrgn;
  int subtmp;
  double sum1;
  double sum2;
  double sum3;
  double sum4;
  double sum5;
  double weit1;
  double weit2;
  double weit3;
  double weit4;
  double weit5;
  double weitp1;
  double weitp2;
  double weitp3;
  double weitp4;
  double *width;
  double *widthl;
  double *z;

  *ifail = 3;
  *relerr = 1.0;
  funcls = 0;

  if ( ndim < 2 )
  {
    *minpts = 0;
    wrkstr[lenwrk-2] = sbrgns;
    *relerr = 1.0;
    *finest = 0.0;
    *ifail = 3;
    return;
  }

  if ( maxpts < *minpts )
  {
    *minpts = 0;
    wrkstr[lenwrk-2] = sbrgns;
    *relerr = 1.0;
    *finest = 0.0;
    *ifail = 3;
    return;
  }

  if ( ndim <= 15 )
  {
    rulcls = i4_power ( 2, ndim ) + 2 * ndim * ndim + 2 * ndim + 1;
  }
  else if ( 15 < ndim )
  {
    rulcls = 1 + ( ndim * ( 12 + ( ndim - 1 ) 
      * ( 6 + ( ndim - 2 ) * 4 ) ) ) / 3;
  }

  if ( maxpts < rulcls )
  {
    *relerr = 1.0;
    *finest = 0.0;
    *ifail = 3;
    return;
  }
//
//  Initialization.
//
  rgnstr = 2 * ndim + 3;
  divaxo = 0;

  center = new double[ndim];
  width = new double[ndim];
  widthl = new double[ndim];
  z = new double[ndim];
//
//  Basic rule initialization.
//
  lambda5 = 9.0 / 19.0;

  if ( ndim <= 15 )
  {
    lambda4 = 9.0 / 10.0;
    lambda2 = 9.0 / 70.0;
    weit5 = 1.0 / pow ( 3.0 * lambda5, 3 ) / pow ( 2.0, ndim );
  }
  else
  {
    ratio = ( double ) ( ndim - 2 ) / 9.0;

    lambda4 = ( 1.0 / 5.0 - ratio ) / ( 1.0 / 3.0 - ratio / lambda5 );

    ratio = ( 1.0 - lambda4 / lambda5 ) 
      * ( double ) ( ndim - 1 ) * ratio / 6.0;

    lambda2 = ( 1.0 / 7.0 - lambda4 / 5.0 - ratio ) 
      / ( 1.0 / 5.0 - lambda4 / 3.0 - ratio / lambda5 );

    weit5 = 1.0 / pow ( 6.0 * lambda5, 3 );
  }

  weit4 = ( 1.0 / 15.0 - lambda5 / 9.0 ) 
    / ( 4.0 * ( lambda4 - lambda5 ) * lambda4 * lambda4 );

  weit3 = ( 1.0 / 7.0 - ( lambda5 + lambda2 ) / 5.0 
    + lambda5 * lambda2 / 3.0 ) / ( 2.0 * lambda4 
    * ( lambda4 - lambda5 ) * ( lambda4 - lambda2 ) ) 
    - 2.0 * ( double ) ( ndim - 1 ) * weit4;

  weit2 = ( 1.0 / 7.0 - ( lambda5 + lambda4 ) / 5.0 
    + lambda5 * lambda4 / 3.0 ) / ( 2.0 * lambda2 
    * ( lambda2 - lambda5 ) * ( lambda2 - lambda4 ) );

  if ( ndim <= 15 )
  {
    weit1 = 1.0 - 2.0 * ( double ) ( ndim )
      * ( weit2 + weit3 + ( double ) ( ndim - 1 ) * weit4 ) 
      - pow ( 2.0, ndim ) * weit5;
  }
  else
  {
    weit1 = 1.0 - 2.0 * ( double ) ndim 
      * ( weit2 + weit3 + ( double ) ( ndim - 1 ) * 
      ( weit4 + 2.0 * ( double ) ( ndim - 2 ) * weit5 / 3.0 ) );
  }

  weitp4 = 1.0 / pow ( 6.0 * lambda4, 2 );

  weitp3 = ( 1.0 / 5.0 - lambda2 / 3.0 ) / 
    ( 2.0 * lambda4 * ( lambda4 - lambda2 ) ) 
    - 2.0 * ( double ) ( ndim - 1 ) * weitp4;

  weitp2 = ( 1.0 / 5.0 - lambda4 / 3.0 ) 
    / ( 2.0 * lambda2 * ( lambda2 - lambda4 ) );

  weitp1 = 1.0 - 2.0 * ( double ) ( ndim ) * 
    ( weitp2 + weitp3 + ( double ) ( ndim - 1 ) * weitp4 );

  ratio = lambda2 / lambda4;

  lambda5 = sqrt ( lambda5 );
  lambda4 = sqrt ( lambda4 );
  lambda2 = sqrt ( lambda2 );
//
//  End basic rule initialization.
//
  if ( *minpts < 0 )
  {
    sbrgns = ( int ) wrkstr[lenwrk-2];
    divflg = 0;
    subrgn = rgnstr;
    wrkstr[lenwrk-1] = wrkstr[lenwrk-1] - wrkstr[subrgn-1];
    *finest = *finest - wrkstr[subrgn-2];
    divaxo = ( int ) wrkstr[subrgn-3];

    for ( j = 1; j <= ndim; j++ )
    {
      subtmp = subrgn - 2 * ( j + 1 );
      center[j-1] = wrkstr[subtmp];
      width[j-1] = wrkstr[subtmp-1];
    }
    width[divaxo-1] = width[divaxo-1] / 2.0;
    center[divaxo-1] = center[divaxo-1] - width[divaxo-1];
  }
  else
  {
    for ( j = 0; j < ndim; j++ )
    {
      width[j] = ( b[j] - a[j] ) / 2.0;
    }
    for ( j = 0; j < ndim; j++ )
    {
      center[j] = a[j] + width[j];
    }

    *finest = 0.0;
    wrkstr[lenwrk-1] = 0.0;
    divflg = 1;
    subrgn = rgnstr;
    sbrgns = rgnstr;
  }
//
//  Begin basic rule.
//
  for ( ; ; )
  {
    rgnvol = pow ( 2.0, ndim ) * r8vec_product ( ndim, width );

    for ( j = 0; j < ndim; j++ )
    {
      z[j] = center[j];
    }

    sum1 = functn ( itest, ndim, z, alpha, beta );
//
//  Compute symmetric sums of functn(lambda2,0,0,...,0) and
//  functn(lambda4,0,0,...,0), and maximum fourth difference.
//
    difmax = -1.0;
    sum2 = 0.0;
    sum3 = 0.0;

    for ( j = 0; j < ndim; j++ )
    {
      z[j] = center[j] - lambda2 * width[j];
      f1 = functn ( itest, ndim, z, alpha, beta );
      z[j] = center[j] + lambda2 * width[j];
      f2 = functn ( itest, ndim, z, alpha, beta );
      widthl[j] = lambda4 * width[j];
      z[j] = center[j] - widthl[j];
      f3 = functn ( itest, ndim, z, alpha, beta );
      z[j] = center[j] + widthl[j];
      f4 = functn ( itest, ndim, z, alpha, beta );
      sum2 = sum2 + f1 + f2;
      sum3 = sum3 + f3 + f4;
      df1 = f1 + f2 - 2.0 * sum1;
      df2 = f3 + f4 - 2.0 * sum1;
      dif = r8_abs ( df1 - ratio * df2 );

      if ( difmax < dif )
      {
        difmax = dif;
        divaxn = j + 1;
      }
      z[j] = center[j];
    }

    if ( sum1 == sum1 + difmax / 8.0 )
    {
      divaxn = ( divaxo % ndim ) + 1;
    }
//
//  Compute symmetric sum of functn(lambda4,lambda4,0,0,...,0).
//
    sum4 = 0.0;

    for ( j = 2; j <= ndim; j++ )
    {
      for ( k = j; k <= ndim; k++ )
      {
        for ( l = 1; l <= 2; l++ )
        {
          widthl[j-2] = -widthl[j-2];
          z[j-2] = center[j-2] + widthl[j-2];
          for ( m = 1; m <= 2; m++ )
          {
            widthl[k-1] = -widthl[k-1];
            z[k-1] = center[k-1] + widthl[k-1];
            sum4 = sum4 + functn ( itest, ndim, z, alpha, beta );
          }
        }
        z[k-1] = center[k-1];
      }
      z[j-2] = center[j-2];
    }
//
//  If NDIM < 16 compute symmetric sum of functn(lambda5,lambda5,...,lambda5).
//
    if ( ndim <= 15 )
    {
      sum5 = 0.0;

      for ( j = 0; j < ndim; j++ )
      {
        widthl[j] = -lambda5 * width[j];
      }
      for ( j = 0; j < ndim; j++ )
      {
        z[j] = center[j] + widthl[j];
      }

      for ( ; ; )
      {
        sum5 = sum5 + functn ( itest, ndim, z, alpha, beta );

        j = ndim;

        for ( ; ; )
        {
          widthl[j-1] = - widthl[j-1];
          z[j-1] = center[j-1] + widthl[j-1];

          if ( 0.0 <= widthl[j-1] )
          {
            break;
          }
          j = j - 1;

          if ( j < 1 )
          {
            break;
          }
        }

        if ( j < 1 )
        {
          break;
        }
      }
    }
//
//  If 15 < NDIM, compute symmetric sum of
//  FUNCTN(lambda5,lambda5,lambda5,0,0,...,0).
//
    else
    {
      sum5 = 0.0;


      for ( j = 0; j < ndim; j++ )
      {
        widthl[j] = lambda5 * width[j];
      }

      for ( i = 3; i <= ndim; i++ )
      {
        for ( j = i; j <= ndim; j++ )
        {
          for ( k = j; k <= ndim; k++ )
          {
            for ( l = 1; l <= 2; l++ )
            {
              widthl[i-3] = -widthl[i-3];
              z[i-3] = center[i-3] + widthl[i-3];
              for ( m = 1; m <= 2; m++ )
              {
                widthl[j-2] = -widthl[j-2];
                z[j-2] = center[j-2] + widthl[j-2];
                for ( n = 1; n <= 2; n++ )
                {
                  widthl[k-1] = -widthl[k-1];
                  z[k-1] = center[k-1] + widthl[k-1];
                  sum5 = sum5 + functn ( itest, ndim, z, alpha, beta );
                }
              }
            }
            z[k-1] = center[k-1];
          }
          z[j-2] = center[j-2];
        }
        z[i-3] = center[i-3];
      }
    }
//
//  Compute fifth and seventh degree rules and error.
//
    rgncmp = rgnvol * ( weitp1 * sum1 
                      + weitp2 * sum2 
                      + weitp3 * sum3 
                      + weitp4 * sum4 );

    rgnval = rgnvol * ( weit1 * sum1 
                      + weit2 * sum2 
                      + weit3 * sum3 
                      + weit4 * sum4 
                      + weit5 * sum5 );

    rgnerr = r8_abs ( rgnval - rgncmp );
//
//  End basic rule.
//
    *finest = *finest + rgnval;
    wrkstr[lenwrk-1] = wrkstr[lenwrk-1] + rgnerr;
    funcls = funcls + rulcls;
//
//  Place results of basic rule into partially ordered list
//  according to subregion error.
//
//  When DIVFLG = 0, start at the top of the list and move down the
//  list tree to find the correct position for the results from the
//  first half of the recently divided subregion.
//
    if ( divflg != 1 )
    {
      for ( ; ; )
      {
        subtmp = 2 * subrgn;
        if ( sbrgns < subtmp )
        {
          break;
        }
        if ( subtmp != sbrgns )
        {
          sbtmpp = subtmp + rgnstr;
          if ( wrkstr[subtmp-1] < wrkstr[sbtmpp-1] )
          {
            subtmp = sbtmpp;
          }
        }
        if ( wrkstr[subtmp-1] <= rgnerr )
        {
          break;
        }
        for ( k = 1; k <= rgnstr; k++ )
        {
          wrkstr[subrgn-k] = wrkstr[subtmp-k];
        }
        subrgn = subtmp;
      }
    }
//
//  When DIVFLG = 1 start at bottom right branch and move up list
//  tree to find correct position for results from second half of
//  recently divided subregion.
//
    else
    {
      for ( ; ; )
      {
        subtmp = ( subrgn / ( 2 * rgnstr ) ) * rgnstr;

        if ( subtmp < rgnstr )
        {
          break;
        }
        if ( rgnerr <= wrkstr[subtmp-1] )
        {
          break;
        }
        for ( k = 1; k <= rgnstr; k++ )
        {
          index1 = subrgn - k + 1;
          index2 = subtmp - k + 1;
          wrkstr[index1-1] = wrkstr[index2-1];
        }
        subrgn = subtmp;
      }
    }
//
//  Store results of basic rule in correct position in list.
//
    wrkstr[subrgn-1] = rgnerr;
    wrkstr[subrgn-2] = rgnval;
    wrkstr[subrgn-3] = divaxn;

    for ( j = 1; j <= ndim; j++ )
    {
      subtmp = subrgn - 2 * ( j + 1 );
      wrkstr[subtmp] = center[j-1];
      wrkstr[subtmp-1] = width[j-1];
    }
//
//  When DIVFLG = 0 prepare for second application of basic rule.
//
    if ( divflg != 1 )
    {
      center[divaxo-1] = center[divaxo-1] + 2.0 * width[divaxo-1];
      sbrgns = sbrgns + rgnstr;
      subrgn = sbrgns;
      divflg = 1;
      continue;
    }
//
//  End ordering and storage of basic rule results.
//  Make checks for possible termination of routine.
//
    *relerr = 1.0;

    if ( wrkstr[lenwrk-1] <= 0.0 )
    {
      wrkstr[lenwrk-1] = 0.0;
    }

    if ( r8_abs ( *finest ) != 0.0 )
    {
      *relerr = wrkstr[lenwrk-1] / r8_abs ( *finest );
    }

    if ( 1.0 < *relerr )
    {
      *relerr = 1.0;
    }

    if ( lenwrk < sbrgns + rgnstr + 2 )
    {
      *ifail = 2;
    }

    if ( maxpts < funcls + 2 * rulcls )
    {
      *ifail = 1;
    }

    if ( *relerr < rel_tol && *minpts <= funcls )
    {
      *ifail = 0;
    }

    if ( *ifail < 3 )
    {
      *minpts = funcls;
      wrkstr[lenwrk-2] = sbrgns;
      break;
    }
//
//  Prepare to use basic rule on each half of subregion with largest
//  error.
//
    divflg = 0;
    subrgn = rgnstr;
    wrkstr[lenwrk-1] = wrkstr[lenwrk-1] - wrkstr[subrgn-1];
    *finest = *finest - wrkstr[subrgn-2];
    divaxo = ( int ) wrkstr[subrgn-3];

    for ( j = 1; j <= ndim; j++ )
    {
      subtmp = subrgn - 2 * ( j + 1 );
      center[j-1] = wrkstr[subtmp];
      width[j-1] = wrkstr[subtmp-1];
    }

    width[divaxo-1] = width[divaxo-1] / 2.0;
    center[divaxo-1] = center[divaxo-1] - width[divaxo-1];
  }

  delete [] center;
  delete [] width;
  delete [] widthl;
  delete [] z;

  return;
}
//****************************************************************************80

double genz_function ( int indx, int ndim, double z[], double alpha[], 
  double beta[] )

//****************************************************************************80
//
//  Purpose:
//
//    GENZ_FUNCTION evaluates one of the test integrand functions.
//
//  Modified:
//
//    26 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan Genz
//    C++ version by John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration:
//    Recent Developments, Software and Applications,
//    edited by Patrick Keast, Graeme Fairweather,
//    D Reidel, 1987, pages 337-340,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int INDX, the index of the test function.
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, double Z[NDIM], the point at which the integrand 
//    is to be evaluated.
//
//    Input, double ALPHA[NDIM], BETA[NDIM], parameters 
//    associated with the integrand function.
//
//    Output, double GENZ_FUNCTION, the value of the test function.
//
{
  int j;
  const double pi = 3.14159265358979323844;
  bool test;
  double total;
  double value;

  value = 0.0;
//
//  Oscillatory.
//
  if ( indx == 1 )
  {
    total = 2.0 * pi * beta[0] + r8vec_sum ( ndim, z );
    value = cos ( total );
  }
//
//  Product Peak.
//
  else if ( indx == 2 )
  {
    total = 1.0;
    for ( j = 0; j < ndim; j++ )
    {
      total = total * (
        1.0 / pow ( alpha[j], 2) + pow ( z[j] - beta[j], 2 ) );
    }
    value = 1.0 / total;
  }
//
//  Corner Peak.
//
  else if ( indx == 3 )
  {
//
//  For this case, the BETA's are used to randomly select
//  a corner for the peak.
//
    total = 1.0;
    for ( j = 0; j < ndim; j++ )
    {
      if ( beta[j] < 0.5 )
      {
        total = total + z[j];
      }
      else
      {
        total = total + alpha[j] - z[j];
      }
    }
    value = 1.0 / pow ( total, ndim + 1 );
  }
//
//  Gaussian.
//  C++ math library complains about things like exp ( -700 )!
//
  else if ( indx == 4 )
  {
    total = 0.0;
    for ( j = 0; j < ndim; j++ )
    {
      total = total + pow ( alpha[j] * ( z[j] - beta[j] ), 2 );
    }
    total = r8_min ( total, 100.0 );
    value = exp ( - total );
  }
//
//  C0 Function.
//
  else if ( indx == 5 )
  {
    total = 0.0;
    for ( j = 0; j < ndim; j++ )
    {
      total = total + alpha[j] * r8_abs ( z[j] - beta[j] );
    }
    value = exp ( - total );
  }
//
//  Discontinuous.
//
  else if ( indx == 6 )
  {
    test = false;

    for ( j = 0; j < ndim; j++ )
    {
      if ( beta[j] < z[j] )
      {
        test = true;
        break;
      }
    }
    if ( test )
    {
      value = 0.0;
    }
    else
    {
      total = r8vec_dot ( ndim, alpha, z );
      value = exp ( total );
    }
  }
  return value;
}
//****************************************************************************80

double genz_integral ( int indx, int ndim, double a[], double b[], 
  double alpha[], double beta[] )

//****************************************************************************80
//
//  Purpose:
//
//    GENZ_INTEGRAL computes the exact integrals of the test functions.
//
//  Modified:
//
//    26 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan Genz
//    C++ version by John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration:
//    Recent Developments, Software and Applications,
//    edited by Patrick Keast, Graeme Fairweather,
//    D Reidel, 1987, pages 337-340,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int INDX, the index of the test.
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, double A[NDIM], B[NDIM], the lower and upper limits
//    of integration.
//
//    Input, double ALPHA[NDIM], BETA[NDIM], parameters 
//    associated with the integrand function.
//
//    Output, double GENZ_INTEGRAL, the exact value of the integral.
//
{
  double ab;
  int *ic;
  int isum;
  int j;
  const double pi = 3.14159265358979323844;
  int rank;
  double s;
  double sgndm;
  double total;
  double value;
//
//  Oscillatory.
//
  if ( indx == 1 )
  {
    value = 0.0;
//
//  Generate all sequences of NDIM 0's and 1's.
//
    rank = 0;
    ic = new int[ndim];

    for ( ; ; )
    {
      tuple_next ( 0, 1, ndim, &rank, ic );

      if ( rank == 0 )
      {
        break;
      }

      total = 2.0 * pi * beta[0];
      for ( j = 0; j < ndim; j++ )
      {
        if ( ic[j] != 1 )
        {
          total = total + alpha[j];
        }
      }

      isum = i4vec_sum ( ndim, ic );

      s = 1 + 2 * ( ( isum / 2 ) * 2 - isum );

      if ( ( ndim % 2 ) == 0 )
      {
        value = value + s * cos ( total );
      }
      else
      {
        value = value + s * sin ( total );
      }
    }
    delete [] ic;

    if ( 1 < ( ndim % 4 ) )
    {
      value = - value;
    }
  }
//
//  Product Peak.
//
  else if ( indx == 2 )
  {
    value = 1.0;

    for ( j = 0; j < ndim; j++ )
    {
      value = value * alpha[j] * ( 
          atan ( ( 1.0 - beta[j] ) * alpha[j] ) 
        + atan (       + beta[j]   * alpha[j] ) );
    }
  }
//
//  Corner Peak.
//
  else if ( indx == 3 )
  {
    value = 0.0;

    sgndm = 1.0;
    for ( j = 1; j <= ndim; j++ )
    {
      sgndm = - sgndm / ( double ) ( j );
    }

    rank = 0;
    ic = new int[ndim];

    for ( ; ; )
    {
      tuple_next ( 0, 1, ndim, &rank, ic );

      if ( rank == 0 )
      {
        break;
      }

      total = 1.0;

      for ( j = 0; j < ndim; j++ )
      {
        if ( ic[j] != 1 )
        {
          total = total + alpha[j];
        }
      }

      isum = i4vec_sum ( ndim, ic );

      s = 1 + 2 * ( ( isum / 2 ) * 2 - isum );
      value = value + ( double ) s / total;

    }

    delete [] ic;

    value = value * sgndm;
  }
//
//  Gaussian.
//
  else if ( indx == 4 )
  {
    value = 1.0;

    ab = sqrt ( 2.0 );
    for ( j = 0; j < ndim; j++ )
    {
      value = value * ( sqrt ( pi ) / alpha[j] ) * 
        (   genz_phi ( ( 1.0 - beta[j] ) * ab * alpha[j] ) 
          - genz_phi (       - beta[j]   * ab * alpha[j] ) );
    }
  }
//
//  C0 Function.
//
  else if ( indx == 5 )
  {
    value = 1.0;
    for ( j = 0; j < ndim; j++ )
    {
      ab = alpha[j] * beta[j];
      value = value * 
        ( 2.0 - exp ( - ab ) - exp ( ab - alpha[j] ) ) / alpha[j];
    }
  }
//
//  Discontinuous.
//
  else if ( indx == 6 )
  {
    value = 1.0;
    for ( j = 0; j < ndim; j++ )
    {
      value = value * ( exp ( alpha[j] * beta[j] ) - 1.0 ) / alpha[j];
    }
  }

  return value;
}
//****************************************************************************80

char *genz_name ( int indx )

//****************************************************************************80
//
//  Purpose:
//
//    GENZ_NAME returns the name of a Genz test integrand.
//
//  Modified:
//
//    26 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int INDX, the index of the test integrand.
//
//    Output, char *GENZ_NAME, the name of the test integrand.
//
{
  char *name;

  name = new char[14];

  if ( indx == 1 )
  {
    strcpy ( name, "Oscillatory  " );
  }
  else if ( indx == 2 )
  {
    strcpy ( name, "Product Peak " );
  }
  else if ( indx == 3 )
  {
    strcpy ( name, "Corner Peak  " );
  }
  else if ( indx == 4 )
  {
    strcpy ( name, "Gaussian     " );
  }
  else if ( indx == 5 )
  {
    strcpy ( name, "C0 Function  " );
  }
  else if ( indx == 6 )
  {
    strcpy ( name, "Discontinuous" );
  }
  else
  {
    cout << "\n";
    cout << "  GENZ_NAME - Fatal error!\n";
    cout << "  1 <= INDX <= 6 is required.\n";
    exit ( 1 );
  }
  return name;
}
//****************************************************************************80

double genz_phi ( double z )

//****************************************************************************80
//
//  Purpose:
//
//    GENZ_PHI estimates the normal cumulative density function.
//
//  Discussion:
//
//    The approximation is accurate to 1.0E-07.
//
//    This routine is based upon algorithm 5666 for the error function,
//    from Hart et al.
//
//  Modified:
//
//    20 March 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan Miller
//    C++ version by John Burkardt
//
//  Reference:
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
//    Input, double Z, a value which can be regarded as the distance,
//    in standard deviations, from the mean.
//
//    Output, double GENZ_PHI, the integral of the normal PDF from negative
//    infinity to Z.
//
//  Local parameters:
//
//    Local, double ROOTPI, despite the name, is actually the 
//    square root of TWO * pi.
//
{
  double expntl;
  double p;
  const double p0 = 220.2068679123761;
  const double p1 = 221.2135961699311;
  const double p2 = 112.0792914978709;
  const double p3 = 33.91286607838300;
  const double p4 = 6.373962203531650;
  const double p5 = 0.7003830644436881;
  const double p6 = 0.03526249659989109;
  const double q0 = 440.4137358247522;
  const double q1 = 793.8265125199484;
  const double q2 = 637.3336333788311;
  const double q3 = 296.5642487796737;
  const double q4 = 86.78073220294608;
  const double q5 = 16.06417757920695;
  const double q6 = 1.755667163182642;
  const double q7 = 0.08838834764831844;
  const double rootpi = 2.506628274631001;
  double zabs;

  zabs = r8_abs ( z );
//
//  12 < |Z|.
//
  if ( 12.0 < zabs )
  {
    p = 0.0;
  }
  else
  {
//
//  |Z| <= 12
//
    expntl = exp ( - zabs * zabs / 2.0 );
//
//  |Z| < 7
//
    if ( zabs < 7.0 )
    {
      p = expntl * (((((( 
                  p6 
         * zabs + p5 ) 
         * zabs + p4 ) 
         * zabs + p3 ) 
         * zabs + p2 ) 
         * zabs + p1 ) 
         * zabs + p0 ) / ((((((( 
                  q7 
         * zabs + q6 ) 
         * zabs + q5 ) 
         * zabs + q4 ) 
         * zabs + q3 ) 
         * zabs + q2 ) 
         * zabs + q1 ) 
         * zabs + q0 );
    }
//
//  CUTOFF <= |Z|
//
    else
    {
      p = expntl / ( 
        zabs + 1.0 / (
        zabs + 2.0 / ( 
        zabs + 3.0 / ( 
        zabs + 4.0 / ( 
        zabs + 0.65 ))))) / rootpi;
    }
  }

  if ( 0.0 < z )
  {
    p = 1.0 - p;
  }

  return p;
}
//****************************************************************************80

double genz_random ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    GENZ_RANDOM is a portable random number generator
//
//  Modified:
//
//    27 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Linus Schrage
//    C++ version by John Burkardt
//
//  Reference:
//
//    Linus Schrage,
//    A More Portable Fortran Random Number Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 2, June 1979, pages 132-138.
//
//  Parameters:
//
//    Input, integer/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, double GENZ_RANDOM, a pseudorandom value.
//
{
  const int a = 16807;
  const int b15 = 32768;
  const int b16 = 65536;
  int fhi;
  int k;
  int leftlo;
  const int p = 2147483647;
  double value;
  int xalo;
  int xhi;

  xhi = *seed / b16;
  xalo = ( *seed - xhi * b16 ) * a;
  leftlo = xalo / b16;
  fhi = xhi * a + leftlo;
  k = fhi / b15;

  *seed = ( 
            ( 
              ( xalo - leftlo * b16 ) - p 
            ) 
          + ( fhi - k * b15 ) * b16 
          ) + k;

  if ( *seed < 0 )
  {
    *seed = *seed + p;
  }

  value = ( double ) ( *seed ) / ( double ) ( p );

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

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
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

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

void multst ( int nsamp, int tstlim, int tstfns[], int tstmax, double difclt[], 
  double expnts[], int ndiml, int ndims[], char *sbname, 
  void subrtn ( int ndim, double a[], double b[], int *minpts, int maxpts, 
    double functn ( int indx, int ndim, double z[], double alpha[], 
      double beta[] ), 
    double rel_tol, int itest, double alpha[], double beta[], int lenwrk, 
    double wrkstr[], double *errest, double *finest, int *ifail ), 
  double rel_tol, int maxpts )

//****************************************************************************80
//
//  Purpose:
//
//    MULTST tests a multidimensional integration routine.
//
//  Discussion:
//
//    The routine uses the Genz test integrand functions, with
//    the user selecting the particular subset of test integrands,
//    the set of difficulty factors, and the spatial dimensions.
//
//  Modified:
//
//    26 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan Genz
//    C++ version by John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration:
//    Recent Developments, Software and Applications,
//    edited by Patrick Keast, Graeme Fairweather,
//    D Reidel, 1987, pages 337-340,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int NSAMP, the number of samples.
//    1 <= NSAMP.
//
//    Input, int TSTLIM, the number of test integrands.
//
//    Input, int TSTFNS[TSTLIM], the indices of the test integrands.
//    Each index is between 1 and 6.
//
//    Input, int TSTMAX, the number of difficulty levels to be tried.
//
//    Input, double DIFCLT[TSTMAX], difficulty levels.
//
//    Input, double EXPNTS[TSTMAX], the difficulty exponents.
//
//    Input, int NDIML, the number of sets of variable sizes.
//
//    Input, int NDIMS[NDIML], the number of variables for the integrals
//    in each test.
//
//    Input, char *SBNAME, the name of the integration
//    subroutine to be tested.
//
//    Input, external SUBRTN, the integration subroutine to be tested.
//
//    Input, double REL_TOL, the relative error tolerance.
//
//    Input, int MAXPTS, the maximum number of integrand calls
//    for all tests.
//
{
# define MXTSFN 6

  double *a;
  double *alpha;
  double *b;
  double *beta;
  double callsa[MXTSFN*MXTSFN];
  double callsb[MXTSFN*MXTSFN];
  double concof;
  double dfact;
  double dfclt;
  int digits;
  double errest;
  double errlog;
  double ersacb[MXTSFN*MXTSFN];
  double ersact[MXTSFN*MXTSFN];
  double ersdsb[MXTSFN*MXTSFN];
  double ersdsc[MXTSFN*MXTSFN];
  double ersesb[MXTSFN*MXTSFN];
  double ersest[MXTSFN*MXTSFN];
  double ersrel[MXTSFN*MXTSFN];
  double estlog;
  double exn;
  double expons[MXTSFN];
  double finest;
  int i;
  int idfclt[MXTSFN];
  int ifail;
  int ifails;
  int it;
  int itest;
  int j;
  int k;
  int lenwrk;
  double medacb[MXTSFN];
  double medacb_med[3];
  double *medact;
  double medact_med[3];
  double medcla[MXTSFN];
  double medcla_med[3];
  double medclb[MXTSFN];
  double medclb_med[3];
  double *medcls;
  double medcls_med[3];
  double meddsb[MXTSFN];
  double meddsb_med[3];
  double *meddsc;
  double meddsc_med[3];
  double medesb[MXTSFN];
  double medesb_med[3];
  double *medest;
  double medest_med[3];
  double medrel;
  double *medrll;
  double medrll_med[3];
  int minpts;
  int n;
  char *name;
  int nconf;
  int ndim;
  int ndimv;
  double qality;
  double *qallty;
  double qallty_med[3];
  double qualty[MXTSFN*MXTSFN];
  int rcalsa;
  int rcalsb;
  double relerr;
  int rulcls;
  int seed;
  double small;
  double tactrb[MXTSFN];
  double tactrb_med[3];
  double tactrs[MXTSFN];
  double tactrs_med[3];
  double tcalsa[MXTSFN];
  double tcalsa_med[3];
  double tcalsb[MXTSFN];
  double tcalsb_med[3];
  double terdsb[MXTSFN];
  double terdsb_med[3];
  double terdsc[MXTSFN];
  double terdsc_med[3];
  double testrb[MXTSFN];
  double testrb_med[3];
  double testrs[MXTSFN];
  double testrs_med[3];
  double tqualt[MXTSFN];
  double tqualt_med[3];
  double total;
  double trelib[MXTSFN];
  double trelib_med[3];
  double value;
  double *wrkstr;

  medact = new double[nsamp];
  medcls = new double[nsamp];
  meddsc = new double[nsamp];
  medest = new double[nsamp];
  medrll = new double[nsamp];
  qallty = new double[nsamp];
//
//  Initialize and compute confidence coefficient.
//
  concof = 0.0;
  nconf = i4_max ( 1, ( 2 * nsamp ) / 5 - 2 );

  for ( i = 1; i <= nconf; i++ )
  {
    concof = 1.0 + ( double ) ( nsamp - nconf + i ) * concof 
      / ( double ) ( nconf - i + 1 );
  }

  concof = 1.0 - concof / ( double ) ( i4_power ( 2, nsamp - 1 ) );

  seed = 123456;

  small = r8_epsilon ( );

  for ( i = 0; i < tstlim; i++ )
  {
    idfclt[i] = ( int ) difclt[tstfns[i]-1];
  }
  for ( i = 0; i < tstlim; i++ )
  {
    expons[i] = expnts[tstfns[i]-1];
  }
//
//  Begin main loop for different numbers of variables.
//
  for ( ndimv = 0; ndimv < ndiml; ndimv++ )
  {
    ndim = ndims[ndimv];

    a = new double[ndim];
    alpha = new double[ndim];
    b = new double[ndim];
    beta = new double[ndim];

    if ( ndim <= 15 )
    {
      rulcls = i4_power ( 2, ndim ) + 2 * i4_power ( ndim, 2 ) + 2 * ndim + 1;
    }
    else
    {
      rulcls = ( ndim * ( 14 - ndim * ( 6 - 4 * ndim ) ) ) / 3 + 1;
    }

    lenwrk = ( 2 * ndim + 3 ) * ( 1 + maxpts / rulcls ) / 2;
    wrkstr = new double[lenwrk];

    if ( ( ndimv % 6 ) == 0 )
    {
      cout << "\n";
      cout << "  Test results with " << nsamp << " samples per test.\n";
      cout << "\n";
      cout << "  Difficulty levels";
      for ( j = 0; j < tstlim; j++ )
      {
        cout << setw(6) << idfclt[j];
      }
      cout << "\n";
      cout << "          Exponents";
      for ( j = 0; j < tstlim; j++ )
      {
        cout << setw(6) << expons[j];
      }
      cout << "\n";

      digits = ( int ) ( -log10 ( rel_tol ) );

      cout << "\n";
      cout << "  Requested digits = " << digits
           << " Maximum values = " << maxpts << "\n";
      cout << "\n";
      cout << sbname << " tests, variable results with confidence "
           << concof << "\n";
      cout << "\n";
      cout << " Vari-  Integrand     Correct digits   Relia-  Wrong"
           << "   Integrand   Quality Total\n";
      cout << " ables              Estimated   Actual bility Digits"
           << "    Values             Fails\n";
      cout << "\n";
    }
//
//  Begin loop for different test integrands.
//
    for ( it = 0; it < tstlim; it++ )
    {
      itest = tstfns[it];
      exn = expnts[itest-1];
      dfclt = difclt[itest-1];

      for ( j = 0; j < ndim; j++ )
      {
        a[j] = 0.0;
      }
      for ( j = 0; j < ndim; j++ )
      {
        b[j] = 1.0;
      }
      ifails = 0;
      medrel = 0;
//
//  Begin loop for different samples.
//
      for ( k = 0; k < nsamp; k++ )
      {
        ifail = 1;
//
//  Choose the integrand function parameters at random.
//
        for ( n = 0; n < ndim; n++ )
        {
          alpha[n] = genz_random ( &seed );
          beta[n] = genz_random ( &seed );
        }
//
//  Modify ALPHA to account for difficulty parameter.
//
        total = r8vec_sum ( ndim, alpha );
        dfact = total * pow ( ndim, exn ) / dfclt;
        for ( j = 0; j < ndim; j++ )
        {
          alpha[j] = alpha[j] / dfact;
        }
//
//  For tests 1 and 3, we modify the value of B.
//
        if ( itest == 1 || itest == 3 )
        {
          for ( j = 0; j < ndim; j++ )
          {
            b[j] = alpha[j];
          }
        }
//
//  For test 6, we modify the value of BETA.
//
        if ( itest == 6 )
        {
          for ( n = 2; n < ndim; n++ )
          {
            beta[n] = 1.0;
          }
        }
//
//  Get the exact value of the integral.
//
        value = genz_integral ( itest, ndim, a, b, alpha, beta );
//
//  Call the integration subroutine.
//
        minpts = 4 * i4_power ( 2, ndim );

        subrtn ( ndim, a, b, &minpts, maxpts, genz_function, rel_tol, 
          itest, alpha, beta, lenwrk, wrkstr, &errest, &finest, &ifail );

        relerr = r8_abs ( ( finest - value ) / value );
        ifails = ifails + i4_min ( ifail, 1 );
        relerr = r8_max ( r8_min ( 1.0, relerr ), small );
        errlog = r8_max ( 0.0, -log10 ( relerr ) );
        errest = r8_max ( r8_min ( 1.0, errest ), small );
        estlog = r8_max ( 0.0, -log10 ( errest ) );
        meddsc[k] = r8_max ( 0.0, estlog - errlog );
        medest[k] = estlog;
        medact[k] = errlog;
        medcls[k] = minpts;

        if ( relerr <= errest )
        {
          medrel = medrel + 1;
        }
      }
//
//  End loop for different samples and compute medians.
//
      r8vec_median ( nsamp, medest, medest_med );
      r8vec_median ( nsamp, medact, medact_med );
      r8vec_median ( nsamp, medcls, medcls_med );
      r8vec_median ( nsamp, meddsc, meddsc_med );

      medrel = medrel / ( double ) ( nsamp );

      trelib[it] = medrel;

      tactrs[it] = medact_med[1];
      testrs[it] = medest_med[1];
      terdsc[it] = meddsc_med[1];
      tcalsa[it] = medcls_med[1];

      tcalsb[it] = medcls_med[2];
      tactrb[it] = medact_med[2];
      testrb[it] = medest_med[2];
      terdsb[it] = meddsc_med[2];

      ersrel[itest-1+ndimv*MXTSFN] = medrel;

      ersest[itest-1+ndimv*MXTSFN] = medest_med[1];
      ersact[itest-1+ndimv*MXTSFN] = medact_med[1];
      ersdsc[itest-1+ndimv*MXTSFN] = meddsc_med[1];

      ersesb[itest-1+ndimv*MXTSFN] = medest_med[2];
      ersacb[itest-1+ndimv*MXTSFN] = medact_med[2];
      ersdsb[itest-1+ndimv*MXTSFN] = meddsc_med[2];

      callsa[itest-1+ndimv*MXTSFN] = medcls_med[1];

      callsb[itest-1+ndimv*MXTSFN] = medcls_med[2];

      qality = 0.0;

      if ( medcls_med[0] != 0.0 )
      {
        qality = ( medact_med[0] + 1.0 ) * 
          ( medest_med[0] + 1.0 - meddsc_med[0] ) / log ( medcls_med[0] );
      }

      tqualt[it] = qality;
      qualty[itest-1+ndimv*MXTSFN] = qality;
      rcalsa = ( int ) medcls_med[1];
      rcalsb = ( int ) medcls_med[2];
      name = genz_name ( itest );

      cout << setw(4) << ndim
           << "  " << setw(14) << name
           << setprecision(2) << setw(4) << medest_med[1]
           << setprecision(2) << setw(5) << medest_med[2]
           << setprecision(2) << setw(5) << medact_med[1]
           << setprecision(2) << setw(5) << medact_med[2]
           << setprecision(3) << setw(5) << medrel
           << setprecision(2) << setw(4) << meddsc_med[1]
           << setprecision(2) << setw(4) << meddsc_med[2]
           << setw(7) << rcalsa
           << setw(8) << rcalsb
           << setprecision(3) << setw(6) << qality
           << setw(5) << ifails << "\n";
      delete [] name;
    }
//
//  End loop for different test integrands.
//
    r8vec_median ( tstlim, tactrs, tactrs_med );
    r8vec_median ( tstlim, trelib, trelib_med );
    r8vec_median ( tstlim, testrs, testrs_med );
    r8vec_median ( tstlim, terdsc, terdsc_med );
    r8vec_median ( tstlim, tactrb, tactrb_med );
    r8vec_median ( tstlim, testrb, testrb_med );
    r8vec_median ( tstlim, terdsb, terdsb_med );
    r8vec_median ( tstlim, tqualt, tqualt_med );
    r8vec_median ( tstlim, tcalsa, tcalsa_med );
    r8vec_median ( tstlim, tcalsb, tcalsb_med );

    rcalsa = ( int ) tcalsa_med[0];
    rcalsb = ( int ) tcalsb_med[0];

    cout << setw(4) << ndim
         << "  Medians       "
         << setprecision(2) << setw(4) << testrs_med[0]
         << setprecision(2) << setw(5) << testrb_med[0]
         << setprecision(2) << setw(5) << testrs_med[0]
         << setprecision(2) << setw(5) << tactrb_med[0]
         << setprecision(3) << setw(5) << trelib_med[0]
         << setprecision(2) << setw(4) << terdsc_med[0]
         << setprecision(2) << setw(4) << terdsb_med[0]
         << setw(7) << rcalsa
         << setw(8) << rcalsb
         << setprecision(3) << setw(6) << tqualt_med[0] << "\n";

    cout << "\n";

    delete [] a;
    delete [] alpha;
    delete [] b;
    delete [] beta;
    delete [] wrkstr;
  }
//
//  End loop for different numbers of variables.
//
  if ( 1 < ndiml )
  {
    cout << "\n";
    cout << "      " << sbname << " Test integrand medians for variables";
    for ( j = 0; j < ndiml; j++ )
    {
      cout << setw(3) << ndims[j];
    }
    cout << "\n";

    cout << "\n";
    cout << "        Integrand     Correct digits   Relia-  Wrong"
         << "   Integrand   Quality\n";
    cout << "                    Estimated   Actual bility digits"
         << "     Values\n";
    cout << "\n";

    for ( it = 0; it < tstlim; it++ )
    {
      itest = tstfns[it];

      for ( j = 0; j < ndiml; j++ )
      {
        medact[j] = ersact[itest-1+j*MXTSFN];
        medest[j] = ersest[itest-1+j*MXTSFN];
        meddsc[j] = ersdsc[itest-1+j*MXTSFN];
        medacb[j] = ersacb[itest-1+j*MXTSFN];
        medesb[j] = ersesb[itest-1+j*MXTSFN];
        meddsb[j] = ersdsb[itest-1+j*MXTSFN];
        medrll[j] = ersrel[itest-1+j*MXTSFN];
        qallty[j] = qualty[itest-1+j*MXTSFN];
        medcla[j] = callsa[itest-1+j*MXTSFN];
        medclb[j] = callsb[itest-1+j*MXTSFN];
      }

      r8vec_median ( ndiml, medrll, medrll_med );
      r8vec_median ( ndiml, medact, medact_med );
      r8vec_median ( ndiml, medest, medest_med );
      r8vec_median ( ndiml, meddsc, meddsc_med );
      r8vec_median ( ndiml, medacb, medacb_med );
      r8vec_median ( ndiml, medesb, medesb_med );
      r8vec_median ( ndiml, meddsb, meddsb_med );
      r8vec_median ( ndiml, qallty, qallty_med );
      r8vec_median ( ndiml, medcla, medcla_med );
      r8vec_median ( ndiml, medclb, medclb_med );

      rcalsa = ( int ) medcla_med[0];
      rcalsb = ( int ) medclb_med[0];
      name = genz_name ( itest );

      cout << "      " << setw(14) << name
           << setprecision(2) << setw(4) << medest_med[0]
           << setprecision(2) << setw(5) << medesb_med[0]
           << setprecision(2) << setw(5) << medact_med[0]
           << setprecision(2) << setw(5) << medacb_med[0]
           << setprecision(3) << setw(5) << medrll_med[0]
           << setprecision(2) << setw(4) << meddsc_med[0]
           << setprecision(2) << setw(4) << meddsb_med[0]
           << setw(7) << rcalsa
           << setw(8) << rcalsb
           << setprecision(3) << setw(6) << qallty_med[0]
           << setw(5) << ifails << "\n";

      delete [] name;

      tactrs[it] = medact_med[0];
      testrs[it] = medest_med[0];
      terdsc[it] = meddsc_med[0];
      tactrb[it] = medacb_med[0];
      testrb[it] = medesb_med[0];
      terdsb[it] = meddsb_med[0];
      tcalsa[it] = medcla_med[0];
      tcalsb[it] = medclb_med[0];
      trelib[it] = medrll_med[0];
      tqualt[it] = qallty_med[0];
    }

    r8vec_median ( tstlim, tactrs, tactrs_med );
    r8vec_median ( tstlim, testrs, testrs_med );
    r8vec_median ( tstlim, terdsc, terdsc_med );
    r8vec_median ( tstlim, tactrb, tactrb_med );
    r8vec_median ( tstlim, testrb, testrb_med );
    r8vec_median ( tstlim, terdsb, terdsb_med );
    r8vec_median ( tstlim, trelib, trelib_med );
    r8vec_median ( tstlim, tqualt, tqualt_med );
    r8vec_median ( tstlim, tcalsa, tcalsa_med );
    r8vec_median ( tstlim, tcalsb, tcalsb_med );

    rcalsa = ( int ) tcalsa_med[0];
    rcalsb = ( int ) tcalsb_med[0];

    cout << "      Global medians"
         << setprecision(2) << setw(4) << testrs_med[0]
         << setprecision(2) << setw(5) << testrb_med[0]
         << setprecision(2) << setw(5) << tactrs_med[0]
         << setprecision(2) << setw(5) << tactrb_med[0]
         << setprecision(3) << setw(5) << trelib_med[0]
         << setprecision(2) << setw(4) << terdsc_med[0]
         << setprecision(2) << setw(4) << terdsb_med[0]
         << setw(7) << rcalsa
         << setw(8) << rcalsb
         << setprecision(3) << setw(6) << tqualt_med[0]
         << setw(5) << ifails << "\n";

    cout << "\n";
  }

  delete [] medact;
  delete [] medcls;
  delete [] meddsc;
  delete [] medest;
  delete [] medrll;
  delete [] qallty;

  return;
# undef MXTSFN
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
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
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
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
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
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

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
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

void r8vec_median ( int n, double r[], double rmed[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEDIAN computes the median of an R8VEC.
//
//  Modified:
//
//    20 March 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan Genz
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the array.
//
//    Input, double R[N], the array to be examined.
//
//    Output, double RMED[3].  RMED[0] contains the median,
//    RMED[1] and RMED[2] specify the confidence interval.
//
{
  int j;
  int k;
  int kmax;
  int nconf;
  int nd;
  double rmax;

  for ( j = 0; j < n; j++ )
  {
    kmax = j;

    for ( k = j + 1; k < n; k++ )
    {
      if ( r[kmax] < r[k] )
      {
        kmax = k;
      }
    }
    rmax = r[kmax];
    r[kmax] = r[j];
    r[j] = rmax;
  }

  nd = n / 2;

  if ( ( n % 2 ) == 0 )
  {
    rmed[0] = ( r[nd-1] + r[nd] ) / 2.0;
  }
  else
  {
    rmed[0] = r[nd];
  }

  nconf = i4_max ( 1, ( 2 * n ) / 5 - 2 );

  rmed[1] = r[n-nconf];
  rmed[2] = r[nconf-1];

  return;
}
//****************************************************************************80

double r8vec_product ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRODUCT returns the product of the entries of an R8VEC.
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
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_PRODUCT, the product of the vector.
//
{
  int i;
  double product;

  product = 1.0;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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

void tuple_next ( int m1, int m2, int n, int *rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT computes the next element of a tuple space.
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between M1 and M2.  The elements are produced one at a time.
//    The first element is
//      (M1,M1,...,M1),
//    the second element is
//      (M1,M1,...,M1+1),
//    and the last element is
//      (M2,M2,...,M2)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2, M1 = 1, M2 = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      Rank   X
//    ----  ---    -----  ---
//    0     * *    1      1 1
//    1     1 1    2      1 2
//    2     1 2    3      1 3
//    3     1 3    4      2 1
//    4     2 1    5      2 2
//    5     2 2    6      2 3
//    6     2 3    7      3 1
//    7     3 1    8      3 2
//    8     3 2    9      3 3
//    9     3 3    0      0 0
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M1, M2, the minimum and maximum entries.
//
//    Input, int N, the number of components.
//
//    Input/output, int *RANK, counts the elements.
//    On first call, set RANK to 0.  Thereafter, the output value of RANK
//    will indicate the order of the element returned.  When there are no
//    more elements, RANK will be returned as 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
{
  int i;
  int j;

  if ( m2 < m1 )
  {
    *rank = 0;
    return;
  }

  if ( *rank <= 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = m1;
    }
    *rank = 1;
  }
  else
  {
    *rank = *rank + 1;
    i = n - 1;

    for ( ; ; )
    {

      if ( x[i] < m2 )
      {
        x[i] = x[i] + 1;
        break;
      }

      x[i] = m1;

      if ( i == 0 )
      {
        *rank = 0;
        for ( j = 0; j < n; j++ )
        {
          x[j] = m1;
        }
        break;
      }
      i = i - 1;
    }
  }
  return;
}

