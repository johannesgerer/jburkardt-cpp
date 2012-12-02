# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa172.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA172_PRB.
//
//  Discussion:
//
//    ASA172_PRB calls the ASA172 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 July 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA172_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA172 library.\n";

  test01 ( );
  test02 ( );

  cout << "\n";
  cout << "ASA172_PRB:\n";
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
//    TEST01 compares indices computed by a triple loop.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 July 2008
//
//  Author:
//
//    John Burkardt
//
{
# define KDIM 3

  int i;
  int i1;
  int i2;
  int i3;
  int ifault;
  int iprod[KDIM];
  int ivec[KDIM];
  int j;
  int jsub;
  int kdim = KDIM;
  int n;
  int nr[KDIM] = { 3, 2, 4 };
  bool qfor;
  bool qind;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  SIMDO can convert between compressed and\n";
  cout << "  vector indices representing a nested loop.\n";
  cout << "\n";
  cout << "  Here, we set QFOR = FALSE, meaning we do\n";
  cout << "  NOT want to convert from FORTRAN ordering\n";
  cout << "  to lexical ordering.\n";
  cout << "\n";
  cout << "  Here, we actually carry out a triple loop\n";
  cout << "  list the indices, and then compare.\n";

  qfor = false;
//
//  If QFOR is FALSE, then the definition of IPROD is reversed...
//
  iprod[0] = nr[kdim-1];
  for ( i = 1; i < kdim; i++ )
  {
    iprod[i] = iprod[i-1] * nr[kdim-1-i];
  }

  n = iprod[kdim-1];
//
//  Carry out the nested loops, and use JSUB to count each iteration.
//  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
//
  jsub = 0;

  cout << "\n";
  cout << "  #1: Generate JSUB by counting as we DO the loops:\n";
  cout << "\n";
  cout << "  DO I1 = 1, N1\n";
  cout << "    DO I2 = 1, N2\n";
  cout << "      DO I3 = 1, N3\n";
  cout << "\n";
  cout << "      JSUB            I1        I2        I3\n";
  cout << "\n";
  for ( i1 = 1; i1 <= nr[0]; i1++ )
  {
    ivec[0] = i1;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i3 = 1; i3 <= nr[2]; i3++ )
      {
        ivec[2] = i3;
        jsub = jsub + 1;
        cout << "  " << setw(8) << jsub << "    "
             << "  " << setw(8) << i1
             << "  " << setw(8) << i2
             << "  " << setw(8) << i3 << "\n";
      }
    }
  }
//
//  Now for each value of JSUB, retrieve the corresponding index subscript.
//  In order to use the QFOR = .FALSE. switch, I apparently have to reverse
//  the sense of the NR vector//
//
  qind = true;

  cout << "\n";
  cout << "  #2: Loop on JSUB, retrieve loop indices\n";
  cout << "      QIND = TRUE J ->I(J)\n";
  cout << "      QFOR = FALSE\n";
  cout << "\n";
  cout << "      JSUB            I1        I2        I3\n";
  cout << "\n";

  for ( j = 1; j <= n; j++ )
  {
    jsub = j;
    ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
    cout << "  " << setw(8) << jsub << "    "
         << "  " << setw(8) << ivec[0]
         << "  " << setw(8) << ivec[1]
         << "  " << setw(8) << ivec[2] << "\n";
  }
//
//  Carry out the nested loops, and DO NOT compute JSUB.
//  Have SIMDO determine JSUB.
//
  qind = false;

  cout << "\n";
  cout << "  #3: For any set of loop indices, retrieve JSUB\n";
  cout << "      QIND = FALSE I(J) -> J\n";
  cout << "      QFOR = FALSE\n";
  cout << "\n";
  cout << "      JSUB            I1        I2        I3\n";
  cout << "\n";
  for ( i1 = 1; i1 <= nr[0]; i1++ )
  {
    ivec[0] = i1;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i3 = 1; i3 <= nr[2]; i3++ )
      {
        ivec[2] = i3;
        ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
        cout << "  " << setw(8) << jsub << "    "
             << "  " << setw(8) << i1
             << "  " << setw(8) << i2
             << "  " << setw(8) << i3 << "\n";
      }
    }
  }
  return;
# undef KDIM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 compares indices computed by a triple loop.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 July 2008
//
//  Author:
//
//    John Burkardt
//
{
# define KDIM 3

  int i;
  int i1;
  int i2;
  int i3;
  int ifault;
  int iprod[KDIM];
  int ivec[KDIM];
  int j;
  int jsub;
  int kdim = KDIM;
  int n;
  int nr[KDIM] = { 3, 2, 4 };
  bool qfor;
  bool qind;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  SIMDO can convert between compressed and\n";
  cout << "  vector indices representing a nested loop.\n";
  cout << "\n";
  cout << "  Here, we set QFOR = TRUE, meaning we DO\n";
  cout << "  want to convert from the FORTRAN \n";
  cout << "  ordering to lexical convention.\n";
  cout << "\n";
  cout << "  Here, we actually carry out a triple loop\n";
  cout << "  list the indices, and then compare.\n";

  qfor = true;

  iprod[0] = nr[0];
  for ( i = 1; i < kdim; i++ )
  {
    iprod[i] = iprod[i-1] * nr[i];
  }

  n = iprod[kdim-1];
//
//  Carry out the nested loops, and use JSUB to count each iteration.
//  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
//
  jsub = 0;

  cout << "\n";
  cout << "  #1: Generate JSUB by counting as we do the loops.\n";
  cout << "\n";
  cout << "  DO I3 = 1, N3\n";
  cout << "    DO I2 = 1, N2\n";
  cout << "      DO I1 = 1, N1\n";
  cout << "\n";
  cout << "      JSUB            I1        I2        I3\n";
  cout << "\n";
  for ( i3 = 1; i3 <= nr[2]; i3++ )
  {
    ivec[2] = i3;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i1 = 1; i1 <= nr[0]; i1++ )
      {
        ivec[0] = i1;
        jsub = jsub + 1;
        cout << "  " << setw(8) << jsub << "    "
             << "  " << setw(8) << i1
             << "  " << setw(8) << i2
             << "  " << setw(8) << i3 << "\n";
      }
    }
  }
//
//  Reverse the order, so that the loop indices are generated in lexical order.
//
  qind = true;

  cout << "\n";
  cout << "  #2: Setting QFOR false means loop indices\n";
  cout << "  are generated in lexical order.\n";
  cout << "      QIND = TRUE J -> I(J)\n";
  cout << "      QFOR = TRUE\n";
  cout << "\n";
  cout << "      JSUB            I1        I2        I3\n";
  cout << "\n";

  for ( j = 1; j <= n; j++ )
  {
    jsub = j;
    ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
    cout << "  " << setw(8) << jsub << "    "
         << "  " << setw(8) << ivec[0]
         << "  " << setw(8) << ivec[1]
         << "  " << setw(8) << ivec[2] << "\n";
  }
//
//  Carry out the nested loops, and DO NOT compute JSUB.
//  Have SIMDO determine JSUB.
//
  qind = false;

  cout << "\n";
  cout << "  #3: For any set of loop indices, retrieve JSUB\n";
  cout << "      QIND = FALSE I(J) -> J\n";
  cout << "      QFOR = TRUE\n";
  cout << "\n";
  cout << "      JSUB            I1        I2        I3\n";
  cout << "\n";
  for ( i3 = 1; i3 <= nr[2]; i3++ )
  {
    ivec[2] = i3;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i1 = 1; i1 <= nr[0]; i1++ )
      {
        ivec[0] = i1;
        ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
        cout << "  " << setw(8) << jsub << "    "
             << "  " << setw(8) << i1
             << "  " << setw(8) << i2
             << "  " << setw(8) << i3 << "\n";
      }
    }
  }
  return;
}
