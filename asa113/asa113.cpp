# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa113.hpp"

//****************************************************************************80

void swap ( double varval[], int klass[], int clsize[], int in, int ik, int iv, 
  double *critvl, int *ntrans, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    SWAP interchanges objects between different classes to improve a criterion.
//
//  Discussion:
//
//    This routine is given a classification of objects, including the
//    number of objects in each class, and the current value of some criterion
//    which is desired to be minimized.
//
//    The routine calculates the change in criterion for all possible swaps,
//    that is, operations in which two objects in different classes exchange 
//    places. Each swap that would result in a lowering of the criterion is 
//    executed, and the related quantities are updated.
//
//    When no more advantageous swaps can be found, the routine returns.
//
//    The routine relies on a user-supplied routine, CRSWAP, to report the
//    expected change in the criterion for a given swap, and to carry
//    out that transfer if requested.
//
//    The variables CLASS and CRITVL have been added to the argument list
//    of CRSWAP.
//
//    Also, the order of the two classes "L" and "M" was interchanged in
//    the call to CRSWAP.  The original order was counterintuitive.
//
//    Sinced CLASS is a reserved keyword in C++, the variable originally
//    named "CLASS" has been unoriginally renamed "KLASS".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Colin Banfield, LC Bassill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Colin Banfield, LC Bassill,
//    Algorithm AS 113:
//    A transfer for non-hierarchichal classification,
//    Applied Statistics,
//    Volume 26, Number 2, 1977, pages 206-210.
//
//  Parameters:
//
//    Input, double VARVAL(IN,IV), the data values.  There are 
//    IN objects, each having spatial dimension IV.
//
//    Input/output, int KLASS(IN), the classification of 
//    each object.
//
//    Input/output, int CLSIZE(IK), the number of objects 
//    in each class.
//
//    Input, int IN, the number of objects.
//
//    Input, int IK, the number of classes.
//
//    Input, int IV, the number of spatial dimensions, 
//    or variates, of the objects.
//
//    Input/output, double *CRITVL, the current value of the criterion.
//
//    Output, int *NTRANS, the number of transfers executed.
//
//    Output, int *IFAULT, error indicator.
//    0, no error detected.
//    1, the number of classes was less than 2.
//    2, the number of objects was less than the number of classes.
//
{
  double eps = 1.0E-38;
  int i;
  int icount;
  double inc;
  int iswitch;
  int it;
  int itop;
  int j;
  int k;
  int l;
  int m;

  if ( ik <= 1 )
  {
    *ifault = 1;
    return;
  }

  if ( in <= ik )
  {
    *ifault = 2;
    return;
  }

  *ifault = 0;
  icount = 0;
  *ntrans = 0;
  itop = ( in * ( in - 1 ) ) / 2;

  i = 1;

  for ( ; ; )
  {
    i = i + 1;

    if ( itop <= icount )
    {
      break;
    }

    if ( in < i )
    {
      i = 1;
      continue;
    }

    l = klass[i-1];
    k = l;
    it = i - 1;
//
//  Test the swap of object I from class M to L, 
//  and object J from class L to M.
//
    for ( j = 1; j <= it; j++ )
    {
      icount = icount + 1;
      m = klass[j-1];

      if ( l != j )
      {
        if ( clsize[l-1] != 1 || clsize[m-1] != 1 )
        {
          iswitch = 1;
          inc = crswap ( varval, klass, clsize, in, ik, iv, critvl, 
            i, j, l, m, iswitch );

          if ( inc < - eps )
          {
            *critvl = *critvl + inc;
            icount = 0;

            iswitch = 2;
            crswap ( varval, klass, clsize, in, ik, iv, critvl, 
              i, j, l, m, iswitch );

            *ntrans = *ntrans + 1;
            klass[i-1] = m;
            klass[j-1] = l;
            l = m;
          }
        }
      }
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

void trnsfr ( double varval[], int klass[], int clsize[], int in, int ik, 
  int iv, double *critvl, int *ntrans, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    TRNSFR transfers objects between classes to improve a criterion.
//
//  Discussion:
//
//    This routine is given a classification of objects, including the
//    number of objects in each class, and the current value of some criterion
//    which is desired to be minimized.
//
//    The routine calculates the change in criterion for all possible transfers
//    of any object from its current class to a different class.  Each transfer
//    that would result in a lowering of the criterion is executed, and the
//    related quantities are updated.
//
//    When no more advantageous transfers can be found, the routine returns.
//
//    The routine relies on a user-supplied routine, CRTRAN, to report the
//    expected change in the criterion for a given transfer, and to carry
//    out that transfer if requested.
//
//    The variables CLASS and CRITVL have been added to the argument list
//    of CRTRAN.
//
//    Also, the order of the two classes "L" and "M" was interchanged in
//    the call to CRTRAN.  The original order was counterintuitive.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Colin Banfield, LC Bassill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Colin Banfield, LC Bassill,
//    Algorithm AS 113:
//    A transfer for non-hierarchichal classification,
//    Applied Statistics,
//    Volume 26, Number 2, 1977, pages 206-210.
//
//  Parameters:
//
//    Input, double VARVAL(IN,IV), the data values.  There are IN 
//    objects, each having spatial dimension IV.
//
//    Input/output, int KLASS(IN), the classification of 
//    each object.
//
//    Input/output, int CLSIZE(IK), the number of objects in 
//    each class.
//
//    Input, int IN, the number of objects.
//
//    Input, int IK, the number of classes.
//
//    Input, int IV, the number of spatial dimensions, or 
//    variates, of the objects.
//
//    Input/output, double *CRITVL, the current value of the criterion.
//
//    Output, int *NTRANS, the number of transfers executed.
//
//    Output, int *IFAULT, error indicator.
//    0, no error detected.
//    1, the number of classes was less than 2.
//    2, the number of objects was less than the number of classes.
//
{
  double eps = 1.0E-38;
  int i;
  int icount;
  double inc;
  double inco;
  int iswitch;
  int l;
  int lo;
  int m;

  if ( ik <= 1 )
  {
    *ifault = 1;
    return;
  }

  if ( in <= ik )
  {
    *ifault = 2;
    return;
  }

  *ifault = 0;
  *ntrans = 0;
  i = 0;
  icount = 0;

  for ( ; ; )
  {
    i = i + 1;

    if ( in <= icount )
    {
      break;
    }

    if ( in < i )
    {
      i = 0;
      icount = 0;
      continue;
    }

    m = klass[i-1];
    if ( clsize[m-1] <= 1 )
    {
      icount = icount + 1;
      continue;
    }

    inco = - eps;
    lo = m;
//
//  Test the transfer of object I from class M to class L.
//
    for ( l = 1; l <= ik; l++ )
    {
      if ( l != m )
      {
        iswitch = 1;
        inc = crtran ( varval, klass, clsize, in, ik, iv, critvl, 
          i, m, l, iswitch );
//
//  Remember the values of L and INC.
//
        if ( inc < inco )
        {
          lo = l;
          inco = inc;
        }
      }
    }

    icount = icount + 1;
//
//  Execute the transfer of object I from class M to class LO.
//
    if ( lo != m )
    {
      l = lo;
      *critvl = *critvl + inco;
      icount = 0;

      iswitch = 2;
      crtran ( varval, klass, clsize, in, ik, iv, critvl, 
        i, m, l, iswitch );

      *ntrans = *ntrans + 1;
      klass[i-1] = l;
      clsize[l-1] = clsize[l-1] + 1;
      clsize[m-1] = clsize[m-1] - 1;
    }
  }
  return;
}
