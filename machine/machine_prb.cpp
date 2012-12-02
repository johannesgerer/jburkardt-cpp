# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "machine.hpp"

int main ( );
void d1mach_prb ( );
void i1mach_prb ( );
void r1mach_prb ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MACHINE_PRB.
//
//  Discussion:
//
//    MACHINE_PRB runs the MACHINE tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "MACHINE_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the MACHINE library.\n";

  d1mach_prb ( );
  i1mach_prb ( );
  r1mach_prb ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "MACHINE_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void d1mach_prb ( )

//****************************************************************************80
//
//  Purpose:
//
//    D1MACH_PRB reports the constants returned by D1MACH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "D1MACH_PRB\n";
  cout << "  D1MACH reports the value of constants associated\n";
  cout << "  with real double precision computer arithmetic.\n";

  cout << "\n";
  cout << "  Assume that double precision numbers are stored\n";
  cout << "  with a mantissa of T digits in base B, with an\n";
  cout << "  exponent whose value must lie between EMIN and EMAX.\n";

  cout << "\n";
  cout << "  For input arguments of 1 <= I <= 5,\n";
  cout << "  D1MACH will return the following values:\n";

  cout << "\n";
  cout << "  D1MACH(1) = B^(EMIN-1), the smallest positive magnitude.\n";
  cout << setw(26) << setprecision(16) << d1mach(1) << "\n";

  cout << "\n";
  cout << "  D1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.\n";
  cout << setw(26) << setprecision(16) << d1mach(2) << "\n";

  cout << "\n";
  cout << "  D1MACH(3) = B^(-T), the smallest relative spacing.\n";
  cout << setw(26) << setprecision(16) << d1mach(3) << "\n";

  cout << "\n";
  cout << "  D1MACH(4) = B^(1-T), the largest relative spacing.\n";
  cout << setw(26) << setprecision(16) << d1mach(4) << "\n";

  cout << "\n";
  cout << "  D1MACH(5) = log10(B).\n";
  cout << setw(26) << setprecision(16) << d1mach(5) << "\n";

  return;
}
//****************************************************************************80

void i1mach_prb ( )

//****************************************************************************80
//
//  Purpose:
//
//    I1MACH_PRB reports the constants returned by I1MACH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "I1MACH_PRB\n";
  cout << "  I1MACH reports the value of constants associated\n";
  cout << "  with integer computer arithmetic.\n";

  cout << "\n";
  cout << "  Numbers associated with input/output units:\n";

  cout << "\n";
  cout << "  I1MACH(1) = the standard input unit.\n";
  cout << i1mach(1) << "\n";

  cout << "\n";
  cout << "  I1MACH(2) = the standard output unit.\n";
  cout << i1mach(2) << "\n";

  cout << "\n";
  cout << "  I1MACH(3) = the standard punch unit.\n";
  cout << i1mach(3) << "\n";

  cout << "\n";
  cout << "  I1MACH(4) = the standard error message unit.\n";
  cout << i1mach(4) << "\n";

  cout << "\n";
  cout << "  Numbers associated with words:\n";

  cout << "\n";
  cout << "  I1MACH(5) = the number of bits per integer.\n";
  cout << i1mach(5) << "\n";

  cout << "\n";
  cout << "  I1MACH(6) = the number of characters per integer.\n";
  cout << i1mach(6) << "\n";

  cout << "\n";
  cout << "  Numbers associated with integer values:\n";

  cout << "\n";
  cout << "  Assume integers are represented in the S digit \n";
  cout << "  base A form:\n";
  cout << "\n";
  cout << "    Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))\n";
  cout << "\n";
  cout << "  where the digits X satisfy 0 <= X(1:S-1) < A.\n";

  cout << "\n";
  cout << "  I1MACH(7) = A, the base.\n";
  cout << i1mach(7) << "\n";

  cout << "\n";
  cout << "  I1MACH(8) = S, the number of base A digits.\n";
  cout << i1mach(8) << "\n";

  cout << "\n";
  cout << "  I1MACH(9) = A^S-1, the largest integer.\n";
  cout << i1mach(9) << "\n";

  cout << "\n";
  cout << "  Numbers associated with floating point values:\n";
  cout << "\n";
  cout << "  Assume floating point numbers are represented \n";
  cout << "  in the T digit base B form:\n";
  cout << "\n";
  cout << "    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B^T) )\n";
  cout << "\n";
  cout << "  where\n";
  cout << "\n";
  cout << "    0 <= X(1:T) < B,\n";
  cout << "    0 < X(1) (unless the value being represented is 0),\n";
  cout << "    EMIN <= E <= EMAX.\n";

  cout << "\n";
  cout << "  I1MACH(10) = B, the base.\n";
  cout << i1mach(10) << "\n";

  cout << "\n";
  cout << "  Numbers associated with single precision values:\n";
  cout << "\n";
  cout << "  I1MACH(11) = T, the number of base B digits.\n";
  cout << i1mach(11) << "\n";

  cout << "\n";
  cout << "  I1MACH(12) = EMIN, the smallest exponent E.\n";
  cout << i1mach(12) << "\n";

  cout << "\n";
  cout << "  I1MACH(13) = EMAX, the largest exponent E.\n";
  cout << i1mach(13) << "\n";

  cout << "\n";
  cout << "  Numbers associated with double precision values:\n";
  cout << "\n";
  cout << "  I1MACH(14) = T, the number of base B digits.\n";
  cout << i1mach(14) << "\n";

  cout << "\n";
  cout << "  I1MACH(15) = EMIN, the smallest exponent E.\n";
  cout << i1mach(15) << "\n";

  cout << "\n";
  cout << "  I1MACH(16) = EMAX, the largest exponent E.\n";
  cout << i1mach(16) << "\n";

  return;
}
//****************************************************************************80

void r1mach_prb ( )

//****************************************************************************80
//
//  Purpose:
//
//    R1MACH_PRB reports the constants returned by R1MACH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "R1MACH_PRB\n";
  cout << "  R1MACH reports the value of constants associated\n";
  cout << "  with real single precision computer arithmetic.\n";

  cout << "\n";
  cout << "  Assume that single precision numbers are stored \n";
  cout << "  with a mantissa of T digits in base B, with an \n";
  cout << "  exponent whose value must lie between EMIN and EMAX.\n";

  cout << "\n";
  cout << "  For input arguments of 1 <= I <= 5,\n";
  cout << "  R1MACH will return the following values:\n";

  cout << "\n";
  cout << "  R1MACH(1) = B^(EMIN-1), the smallest positive magnitude.\n";
  cout << setw(26) << setprecision(16) << r1mach(1) << "\n";

  cout << "\n";
  cout << "  R1MACH(2) = B^EMAX*(1-B**(-T)), the largest magnitude.\n";
  cout << setw(26) << setprecision(16) << r1mach(2) << "\n";

  cout << "\n";
  cout << "  R1MACH(3) = B^(-T), the smallest relative spacing.\n";
  cout << setw(26) << setprecision(16) << r1mach(3) << "\n";

  cout << "\n";
  cout << "  R1MACH(4) = B^(1-T), the largest relative spacing.\n";
  cout << setw(26) << setprecision(16) << r1mach(4) << "\n";

  cout << "\n";
  cout << "  R1MACH(5) = log10(B).\n";
  cout << setw(26) << setprecision(16) << r1mach(5) << "\n";

  return;
}
