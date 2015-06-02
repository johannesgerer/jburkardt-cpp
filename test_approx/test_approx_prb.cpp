# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "test_approx.hpp"
# include "spline.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_APPROX_PRB.
//
//  Discussion:
//
//    TEST_APPROX_PRB tests the TEST_APPROX library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp (  );
  cout << "\n";
  cout << "TEST_APPROX_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_APPROX library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_APPROX_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 shows how P00_TITLE can be called.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;
  string title;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Demonstrate some of the bookkeeping routines.\n";
  cout << "  P00_PROB_NUM returns the number of problems.\n";
  cout << "  P00_TITLE returns the problem title.\n";
  cout << "  P00_LIMIT returns the problem limits.\n";

  prob_num = p00_prob_num ( );

  cout << "\n";
  cout << "  Number of problems = " << prob_num << "\n";
  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );
    cout << "  " << setw(2) << prob
         << "  \"" << title << "\"\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 shows how P00_STORY can be called.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  P00_STORY prints the problem \"story\".\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    cout << "\n";
    cout << "  Problem " << prob << "\n";
    p00_story ( prob );
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses polynomial interpolation on data vector problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *diftab;
  int i;
  int j;
  int jhi;
  char mark;
  int max_tab = 12;
  int ntab;
  int data_num;
  int prob;
  int prob_num;
  string title;
  double x;
  double *xdata;
  double yapprox;
  double *ydata;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Polynomial interpolation to a vector of data.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );

    cout << "\n";
    cout << "  Problem " << prob << "\n";
    cout << "  " << title << "\n";

    data_num = p00_data_num ( prob );

    cout << "  DATA_NUM = " << data_num << "\n";

    if ( max_tab < data_num )
    {
      cout << "\n";
      cout << "  Skipped problem " << prob << "\n";
      cout << "  Too big.\n";
    }
    else
    {
      xdata = new double[data_num];
      ydata = new double[data_num];
      diftab = new double[data_num];

      p00_dat ( prob, data_num, xdata, ydata );

      ntab = data_num;

      cout << "\n";
      cout << "  Interpolating polynomial order = " << ntab << "\n";
      cout << "\n";
//
//  Construct the interpolating polynomial via finite differences.
//
      data_to_dif ( ntab, xdata, ydata, diftab );
//
//  Print out the approximation, including midpoints of the intervals.
//
      for ( i = 1; i <= ntab; i++ )
      {
        if ( i < ntab )
        {
          jhi = 2;
        }
        else
        {
          jhi = 1;
        }

        for ( j = 1; j <= jhi; j++ )
        {
          if ( i < ntab )
          {
            x = ( ( double ) ( jhi - j + 1 ) * xdata[i-1]   
                + ( double ) (       j - 1 ) * xdata[i] ) 
                / ( double ) ( jhi         );
          }
          else
          {
            x = xdata[ntab-1];
          }

          if ( j == 1 )
          {
            mark = '*';
          }
          else
          {
            mark = ' ';
          }

          yapprox = dif_val ( ntab, xdata, diftab, x );

          cout << "  " << mark
               << "  " << setw(14) << x
               << "  " << setw(14) << yapprox << "\n";
        }
      }
      delete [] diftab;
      delete [] xdata;
      delete [] ydata;
    }
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 uses linear spline interpolation on all problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int imax;
  char mark;
  int data_num;
  int prob;
  int prob_num;
  string title;
  double *xdata;
  double xval;
  double *ydata;
  double ypval;
  double yval;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Linear spline interpolation.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );

    data_num = p00_data_num ( prob );
    xdata = new double[data_num];
    ydata = new double[data_num];

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];

    cout << "\n";
    cout << "  Problem " << prob << "\n";
    cout << "  " << title << "\n";
    cout << "\n";
    cout << "       X          Y          Y'\n";
    cout << "\n";
//
//  Evaluate the interpolation function.
//
    imax = 2 * data_num - 1;

    for ( i = 1; i <= imax; i++ )
    {
      xval = ( ( double ) ( imax - i     ) * a   
             + ( double ) (        i - 1 ) * b ) 
             / ( double ) ( imax     - 1 );

      spline_linear_val ( data_num, xdata, ydata, xval, &yval, &ypval );

      if ( ( i % 2 ) == 1 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout << "  " << mark
           << "  " << setw(14) << xval
           << "  " << setw(14) << yval
           << "  " << setw(14) << ypval << "\n";
    }
    delete [] xdata;
    delete [] ydata;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 uses Overhauser spline interpolation on all problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int j;
  int jhi;
  int jmax;
  char mark;
  int data_num;
  int num_dim = 1;
  int prob;
  int prob_num;
  string title;
  double *xdata;
  double xval;
  double *ydata;
  double yval;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Overhauser spline interpolation.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );

    data_num = p00_data_num ( prob );

    xdata = new double[data_num];
    ydata = new double[data_num];

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];

    cout << "\n";
    cout << "  Problem " << prob << "\n";
    cout << "  " << title << "\n";
    cout << "\n";
    cout << "  X   Y\n";
    cout << "\n";
//
//  Evaluate the interpolation function.
//
    for ( i = 1; i < data_num; i++ )
    {
      jmax = 3;

      if ( i == data_num - 1 )
      {
        jhi = jmax;
      }
      else
      {
        jhi = jmax - 1;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        xval = ( ( double ) ( jmax - j     ) * xdata[i-1]     
               + ( double ) (        j - 1 ) * xdata[i] ) 
               / ( double ) ( jmax     - 1 );

        spline_overhauser_val ( num_dim, data_num, xdata, ydata, xval, &yval );

        if ( j == 1 || j == 3 )
        {
          mark = '*';
        }
        else
        {
          mark = ' ';
        }
        cout << "  " << mark
             << "  " << setw(14) << xval
             << "  " << setw(14) << yval << "\n";
      }

    }
    delete [] xdata;
    delete [] ydata;
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 uses cubic spline interpolation on all problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int jmax;
  char mark;
  int data_num;
  int prob;
  int prob_num;
  string title;
  double *xdata;
  double xval;
  double ybcbeg;
  double ybcend;
  double *ydata;
  double *ypp;
  double yppval;
  double ypval;
  double yval;

  ibcbeg = 0;
  ibcend = 0;
  ybcbeg = 0.0;
  ybcend = 0.0;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Cubic spline interpolation.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );

    data_num = p00_data_num ( prob );

    xdata = new double[data_num];
    ydata = new double[data_num];
    ypp = new double[data_num];

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];
//
//  Set up the interpolation function.
//
    ypp = spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, 
      ibcend, ybcend );

    cout << "\n";
    cout << "  Problem " << prob << "\n";
    cout << "  " << title << "\n";
    cout << "\n";
    cout << "    X   Y\n";
    cout << "\n";
//
//  Evaluate the interpolation function.
//
    for ( i = 1; i < data_num; i++ )
    {
      jmax = 3;

      if ( i == data_num - 1 )
      {
        jhi = jmax;
      }
      else
      {
        jhi = jmax - 1;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        xval = ( ( double ) ( jmax - j     ) * xdata[i-1]
               + ( double ) (        j - 1 ) * xdata[i] ) 
               / ( double ) ( jmax     - 1 );

        yval = spline_cubic_val ( data_num, xdata, ydata, ypp, xval, &ypval, &yppval );

        if ( j == 1 || j == 3 )
        {
          mark = '*';
        }
        else
        {
          mark = ' ';
        }
        cout << "  " << mark
             << "  " << setw(14) << xval
             << "  " << setw(14) << yval << "\n";
      }
    }
    delete [] xdata;
    delete [] ydata;
    delete [] ypp;
  }

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 plots an Overhauser spline interpolant for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  string approx_filename = "test07_approx.txt";
  string data_filename = "test07_data.txt";
  int i;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int num_dim = 1;
  int plot;
  int prob;
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Plot an Overhauser spline interpolant for problem 7.\n";
//
//  Get the problem data.
//
  prob = 7;

  data_num = p00_data_num ( prob );

  xdata = new double[data_num];
  ydata = new double[data_num];

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  cout << "\n";
  cout << "  Data values stored in \"" << data_filename << "\".\n";
//
//  Evaluate the approximating function.
//
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = new double[nplot];
  yplot = new double[nplot];

  plot = 0;

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1] 
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      spline_overhauser_val ( num_dim, data_num, xdata, ydata, xval, &yval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  cout << "  Approximant values stored in \"" << approx_filename << "\".\n"; 

  delete [] xdata;
  delete [] xplot;
  delete [] ydata;
  delete [] yplot;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 plots a cubic spline interpolant for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  string approx_filename = "test08_approx.txt";
  string data_filename = "test08_data.txt";
  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int plot;
  int prob;
  double *xdata;
  double *xplot;
  double xval;
  double ybcbeg;
  double ybcend;
  double *ydata;
  double *yplot;
  double *ypp;
  double yppval;
  double ypval;
  double yval;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Plot a cubic spline interpolant for problem 7.\n";

  prob = 7;
//
//  Get the data.
//
  data_num = p00_data_num ( prob );

  xdata = new double[data_num];
  ydata = new double[data_num];
  ypp = new double[data_num];

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  cout << "\n";
  cout << "  Data values stored in \"" << data_filename << "\".\n";
//
//  Set up the interpolation function.
//
  ibcbeg = 0;
  ibcend = 0;
  ybcbeg = 0.0;
  ybcend = 0.0;

  ypp = spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, ibcend, ybcend );
//
//  Evaluate the interpolation function.
//
  plot = 0;
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = new double[nplot];
  yplot = new double[nplot];

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }
    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1]
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      yval = spline_cubic_val ( data_num, xdata, ydata, ypp, xval, &ypval, 
        &yppval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  cout << "  Approximant values stored in \"" << approx_filename << "\".\n"; 

  delete [] xdata;
  delete [] xplot;
  delete [] ydata;
  delete [] yplot;
  delete [] ypp;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 uses B spline approximation on all problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int j;
  int jhi;
  int jmax;
  char mark;
  int data_num;
  int prob;
  int prob_num;
  string title;
  double *xdata;
  double xval;
  double *ydata;
  double yval;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  B spline approximation.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );

    data_num = p00_data_num ( prob );

    xdata = new double[data_num];
    ydata = new double[data_num];

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];

    cout << "\n";
    cout << "  Problem " << prob << "\n";
    cout << "  " << title << "\n";
    cout << "\n";
    cout << "       X        Y\n";
    cout << "\n";
//
//  Evaluate the interpolation function.
//
    for ( i = 1; i < data_num; i++ )
    {
      jmax = 3;

      if ( i == data_num - 1 )
      {
        jhi = jmax;
      }
      else
      {
        jhi = jmax - 1;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        xval = ( ( double ) ( jmax - j     ) * xdata[i-1]
               + ( double ) (        j - 1 ) * xdata[i] ) 
               / ( double ) ( jmax     - 1 );

        yval = spline_b_val ( data_num, xdata, ydata, xval );

        if ( j == 1 || j == 3 )
        {
          mark = '*';
        }
        else
        {
          mark = ' ';
        }
        cout << "  " << mark
             << "  " << setw(14) << xval
             << "  " << setw(14) << yval << "\n";
      }
    }
    delete [] xdata;
    delete [] ydata;
  }

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 plots a B spline approximant for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  string approx_filename = "test10_approx.txt";
  string data_filename = "test10_data.txt";
  int i;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int plot;
  int prob;
  string title;
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  Plot a B spline approximant for problem 7\n";

  prob = 7;

  title = p00_title ( prob );
//
//  Get the data.
//
  data_num = p00_data_num ( prob );

  xdata = new double[data_num];
  ydata = new double[data_num];

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  cout << "\n";
  cout << "  Data values stored in \"" << data_filename << "\".\n";
//
//  Evaluate the approximation function.
//
  plot = 0;
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = new double[nplot];
  yplot = new double[nplot];

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }
    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1]  
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      yval = spline_b_val ( data_num, xdata, ydata, xval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  cout << "  Approximant values stored in \"" << approx_filename << "\".\n"; 

  delete [] xdata;
  delete [] xplot;
  delete [] ydata;
  delete [] yplot;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 plots a beta spline approximant for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  string approx_filename = "test11_approx.txt";
  double beta1;
  double beta2;
  string data_filename = "test11_data.txt";
  int i;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int plot;
  int prob;
  string title;
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  beta1 = 100.0;
  beta2 = 0.0;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  Plot a beta spline approximant for problem 7\n";
  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";

  prob = 7;

  title = p00_title ( prob );
//
//  Get the data.
//
  data_num = p00_data_num ( prob );

  xdata = new double[data_num];
  ydata = new double[data_num];

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  cout << "\n";
  cout << "  Data values stored in \"" << data_filename << "\".\n";
//
//  Evaluate the interpolation function.
//
  plot = 0;
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = new double[nplot];
  yplot = new double[nplot];

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1] 
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      yval = spline_beta_val ( beta1, beta2, data_num, xdata, ydata, xval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  cout << "  Approximant values stored in \"" << approx_filename << "\".\n"; 

  delete [] xdata;
  delete [] xplot;
  delete [] ydata;
  delete [] yplot;

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 plots a Bernstein spline approximant for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  string approx_filename = "test12_approx.txt";
  double b;
  string data_filename = "test12_data.txt";
  int i;
  int data_num;
  int nplot = 101;
  int plot;
  int prob;
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  Plot a Bernstein approximant for problem 5.\n";
  cout << "  Note that the Bernstein approximant requires equally\n";
  cout << "  spaced data!\n";

  prob = 5;
//
//  Get the data.
//
  data_num = p00_data_num ( prob );

  xdata = new double[data_num];
  ydata = new double[data_num];

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  cout << "\n";
  cout << "  Data values stored in \"" << data_filename << "\".\n";
//
//  Evaluate the approximant function.
//
  xplot = new double[nplot];
  yplot = new double[nplot];

  a = xdata[0];
  b = xdata[data_num-1];

  for ( plot = 1; plot <= nplot; plot++ )
  {
    xval = ( ( double ) ( nplot - plot     ) * a     
           + ( double ) (         plot - 1 ) * b ) 
           / ( double ) ( nplot        - 1 );

    yval = bpab_approx ( data_num - 1, a, b, ydata, xval );

    xplot[plot-1] = xval;
    yplot[plot-1] = yval;
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  cout << "  Approximant values stored in \"" << approx_filename << "\".\n"; 

  delete [] xdata;
  delete [] xplot;
  delete [] ydata;
  delete [] yplot;

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 plots a cubic spline interpolant for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2011
//
//  Author:
//
//    John Burkardt
//
{
# define NPLOT 101

  string approx_filename = "test13_approx.txt";
  string data_filename = "test13_data.txt";
  int ibcbeg;
  int ibcend;
  int j;
  int data_num;
  int nplot = NPLOT;
  int prob;
  double *xdata;
  double xplot[NPLOT];
  double xval;
  double ybcbeg;
  double ybcend;
  double *ydata;
  double yplot[NPLOT];
  double *ypp;
  double yppval;
  double ypval;
  double yval;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  Plot a cubic spline interpolant for problem 5\n";

  prob = 5;

  data_num = p00_data_num ( prob );

  xdata = new double[data_num];
  ydata = new double[data_num];
  ypp = new double[data_num];

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  cout << "\n";
  cout << "  Data values stored in \"" << data_filename << "\".\n";
//
//  Set up the interpolation function.
//
  ibcbeg = 0;
  ibcend = 0;
  ybcbeg = 0.0;
  ybcend = 0.0;

  ypp = spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, ibcend, ybcend );
//
//  Evaluate the interpolation function.
//
  for ( j = 1; j <= nplot; j++ )
  {
    xval = ( ( double ) ( nplot - j     ) * xdata[0]
           + ( double ) (         j - 1 ) * xdata[data_num-1] ) 
           / ( double ) ( nplot     - 1 );

    yval = spline_cubic_val ( data_num, xdata, ydata, ypp, xval, &ypval, 
      &yppval );

    xplot[j-1] = xval;
    yplot[j-1] = yval;
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  cout << "  Approximant values stored in \"" << approx_filename << "\".\n"; 

  delete [] xdata;
  delete [] ydata;
  delete [] ypp;

  return;
# undef NPLOT
}
