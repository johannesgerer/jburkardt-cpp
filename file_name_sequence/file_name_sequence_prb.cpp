# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <sstream>

using namespace std;

# include "file_name_sequence.hpp"

int main ( );
void test02 ( string prefix, string suffix, int first, int last );
void test03 ( string prefix, string suffix, int first, int last );
void test04 ( string filename, int filename_num );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FILE_NAME_SEQUENCE_PRB.
//
//  Discussion:
//
//    FILE_NAME_SEQUENCE_PRB tests the FILE_NAME_SEQUENCE library.
//
//    There are situations such as animations or parallel processing in which
//    it is necessary to generate a sequence of file names which include
//    an embedded index that increases.  A simple example might be
//
//      "fred0.txt", "fred1.txt", "fred2.txt"
//
//    A side issue arises when the number of files is large enough that the
//    number of digits in the index will vary.  Thus, if we are going to have
//    15 files, do we want to number them as
//
//      "fred00.txt" through "fred14.txt"
//
//    which means, for one thing, that they will alphabetize properly, or
//    will we be satisfied with
//
//      "fred0.txt" through "fred14.txt" ?
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "frodo_01345_lives.txt";

  timestamp ( );
  cout << "\n";
  cout << "FILE_NAME_SEQUENCE_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the FILE_NAME_SEQUENCE library.\n";

  test02 ( "fred", ".txt", 0, 12 );
  test03 ( "frid", ".txt", 99, 105 );
  test04 ( filename, 10 );
//
//  Terminate.
//
  cout << "\n";
  cout << "FILE_NAME_SEQUENCE_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test02 ( string prefix, string suffix, int first, int last )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses I4_TO_S().
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int i;
  string number;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  FILENAME(I) = PREFIX + i4_to_s ( I ) + SUFFIX )\n";
  cout << "  PREFIX = \"" << prefix << "\"\n";
  cout << "  SUFFIX = \"" << suffix << "\"\n";
  cout << "  " << first << " <= I <= " << last << "\n";
  cout << "  Numbers do NOT include leading zeros.\n";
  cout << "\n";

  for ( i = first; i <= last; i++ )
  {
    number = i4_to_s ( i );
    filename = prefix + number + suffix;
    cout << "  " << setw(4) << i
         << "  \"" << filename << "\".\n";
  }

  return;
}
//****************************************************************************80

void test03 ( string prefix, string suffix, int first, int last )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses I4_TO_S0().
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int digits = 3;
  int i;
  string number;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  FILENAME(I) = PREFIX + i4_to_s0 ( I, 3 ) + SUFFIX )\n";
  cout << "  PREFIX = \"" << prefix << "\"\n";
  cout << "  SUFFIX = \"" << suffix << "\"\n";
  cout << "  " << first << " <= I <= " << last << "\n";
  cout << "  Numbers DO include leading zeros.\n";
  cout << "\n";

  for ( i = first; i <= last; i++ )
  {
    number = i4_to_s0 ( i, digits );
    filename = prefix + number + suffix;
    cout << "  " << setw(4) << i
         << "  \"" << filename << "\".\n";
  }
  return;
}
//****************************************************************************80

void test04 ( string filename, int filename_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 uses FILENAME_INC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  FILENAME(I+1) = FILENAME_INC ( FILENAME(I) )\n";
  cout << "  First FILENAME = \"" << filename << "\"\n";
  cout << "  Number of filenames = " << filename_num << "\n";
  cout << "  Numbers may include leading zeros.\n";
  cout << "\n";

  for ( i = 1; i <= filename_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  \"" << filename << "\".\n";
    filename_inc ( &filename );
  }

  return;
}
