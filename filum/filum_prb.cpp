# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "filum.hpp"

int main ( );
void test02 ( );
void test03 ( );
void test06 ( );
void test08 ( );
void test085 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test165 ( );
void test17 ( );
void test22 ( );
void test225 ( );
void test24 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FILUM_PRB.
//
//  Discussion:
//
//    FILUM_PRB tests the FILUM library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FILUM_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the FILUM library.\n";

  test02 ( );
  test03 ( );
  test06 ( );
  test08 ( );
  test085 ( );

  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test165 ( );
  test17 ( );

  test22 ( );
  test225 ( );
  test24 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FILUM_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests FILE_CHAR_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "filum_prb_test.txt";

  cout << "\n";
  cout << "TEST02\n";
  cout << "  FILE_CHAR_COUNT counts the characters in a file.\n";
  cout << "\n";
  cout << "  Examining file \"" << filename << "\".\n";
  cout << "\n";
  cout << "  Number of characters: " 
       << file_char_count ( filename ) << ".\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests FILE_COLUMN_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2009
//
//  Author:
//
//    John Burkardt 
//
{
  int column_num;
  string filename = "filum_prb_4by5.txt";

  cout << "\n";
  cout << "TEST03\n";
  cout << "  FILE_COLUMN_COUNT counts the columns in a file.\n";
  cout << "\n";
  cout << "  It is assumed that the file contains a number of lines,\n";
  cout << "  with each line containing the same number of words.\n";
  cout << "  The task is to determine the number of words in a line,\n";
  cout << "  that is, the number of \"columns\" of text.\n";

  cout << "\n";
  cout << "  Examining the file:\"" << filename << "\".\n";

  column_num = file_column_count ( filename );

  cout << "\n";
  cout << "  Number of columns: " << column_num << "\n";

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests FILE_EXIST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  string filename;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  FILE_EXIST reports whether a file 'exists'.\n";
  cout << "\n";
  cout << "  Exist?   File_name\n";
  cout << "\n";

  filename = "filum_prb.C";

  cout << "       "
       << file_exist ( filename ) << "  "
       << filename << "\n";

  filename = "filum.C";
  cout << "       "
       << file_exist ( filename ) << "  "
       << filename << "\n";

  filename = "raisin.txt";
  cout << "       "
       << file_exist ( filename ) << "  "
       << filename << "\n";

  filename = "make.money.fast";
  cout << "       "
       << file_exist ( filename ) << "  "
       << filename << "\n";

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests FILE_LINE_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "filum_prb_test.txt";
  string line;
//char line[80];
  int line_index;
  int line_num;
  int seed;
  int test;
  int test_num = 5;

  seed = 123456789;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  FILE_LINE_UNIFORM selects one line of a file,\n";
  cout << "  uniformly at random, reading the file one time.\n";
  cout << "\n";
  cout << "  Examining the file: \"" << filename << "\".\n";

  cout << "\n";
  cout << "      Line  Text:\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    line = file_line_uniform ( filename, &seed, &line_index, &line_num );
    cout << "  " << setw(8) << line_index
         << "  \"" << line << "\".\n";
  }

  cout << "\n";
  cout << "  Total number of lines in file is " << line_num << "\n";

  return;
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests FILE_LINE_WIDTH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "filum_prb_test.txt";

  cout << "\n";
  cout << "TEST085\n";
  cout << "  FILE_LINE_WIDTH counts the longest line in a file.\n";
  cout << "\n";
  cout << "  Examining file \"" << filename << "\".\n";
  cout << "\n";
  cout << "  Longest line is: " 
       << file_line_width ( filename ) << " characters.\n";

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests FILENAME_DEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  string filename;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  FILENAME_DEC decrements a string\n";
  cout << "\n";
  cout << "     Input             Output\n";

  for ( i = 0; i < 4; i++ )
  {
    if ( i == 0 )
    {
      filename = "file_02_99.dat";
    }
    else if ( i == 1 ) 
    {
      filename = "file076.dat";
    }
    else if ( i == 2 ) 
    {
      filename = "3cat3.dat  ";
    }
    else if ( i == 3 ) 
    {
      filename = "fred03.txt ";
    }
    cout << "\n";
    for ( j = 1; j <= 4; j++ )
    {
      cout << "  " << setw(11) << filename << "  ";

      filename_dec ( &filename );

      cout << "  " << setw(11) << filename << "\n";

    }
  }

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests FILENAME_EXT_GET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int itest;
  int j;
  int ntest = 5;
  string text;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  FILENAME_EXT_GET finds a file extension.\n";
  cout << "\n";
  cout << "  FILENAME     Begin    End\n";
  cout << "\n";

  for (  itest = 0; itest < ntest; itest++ )
  {
    if ( itest == 0 ) 
    {
      text = "bob.for";
    }
    else if ( itest == 1 )
    {
      text = "N.B.C.D";
    }
    else if ( itest == 2 )
    {
      text = "Naomi.";
    }
    else if ( itest == 3 )
    {
      text = "Arthur";
    }
    else if ( itest == 4 )
    {
      text = ".amos";
    }

    filename_ext_get ( text, &i, &j );

    cout << "  " << setw(10) << text 
         << "  " << setw(6)  << i
         << "  " << setw(6)  << j      << "\n";
  }

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests FILENAME_EXT_SWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  string ext;
  string filename;
  string filename2;
  int itest;
  int ntest = 5;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  FILENAME_EXT_SWAP changes a file extension.\n";
  cout << "\n";
  cout << "  FILENAME   EXT   FILENAME_EXT_SWAP\n";
  cout << "\n";

  for ( itest = 0; itest < ntest; itest++ )
  {
    if ( itest == 0 ) 
    {
      filename = "bob.for";
      ext = "obj";
    }
    else if ( itest == 1 )
    {
      filename = "bob.bob.bob";
      ext = "txt";
    }
    else if ( itest == 2 )
    {
      filename = "bob.";
      ext = "yak";
    }
    else if ( itest == 3 )
    {
      filename = "bob";
      ext = "ps";
    }
    else if ( itest == 4 )
    {
      filename = ".oops";
      ext ="eek";
    }

    cout << "  " << setw(12) << filename  
         << "  " << setw(3)  << ext
         << "  ";

    filename2 = filename_ext_swap ( filename, ext );

    cout << "  " << setw(12) << filename2 << "\n";
  }

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests FILENAME_INC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  string filename;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  FILENAME_INC increments a string\n";
  cout << "\n";
  cout << "     Input             Output\n";

  for ( i = 0; i < 4; i++ )
  {
    if ( i == 0 )
    {
      filename = "file???.dat";
    }
    else if ( i == 1 ) 
    {
      filename = "file072.dat";
    }
    else if ( i == 2 ) 
    {
      filename = "2cat9.dat  ";
    }
    else if ( i == 3 ) 
    {
      filename = "fred98.txt ";
    }
    cout << "\n";
    for ( j = 1; j <= 4; j++ )
    {
      cout << "  " << setw(11) << filename << "  ";

      filename_inc ( &filename );

      cout << "  " << setw(11) << filename << "\n";

      if ( s_len_trim ( filename ) <= 0 )
      {
        cout << "  (File name not incrementable.  Quit loop!)\n";
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests FILENAME_INC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int sim;
  string filename;
  string string1;
  string string2;
  int time_step;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  FILENAME_INC increments a string.\n";
  cout << "\n";
  cout << "  In this example, the string is a file name.\n";
  cout << "  The user is going to carry out several simulations.\n";
  cout << "  For each simulation, a consecutive series of file\n";
  cout << "  names is desired.\n";

  string1 = "file_00";
  string2 = "_000.txt";

  for ( sim = 1; sim < 5; sim++ )
  {
    cout << "\n";
    cout << "  Simulation " << sim << " begins.\n";
    cout << "\n";

    filename_inc ( &string1 );

    filename = string1 + string2;

    for ( time_step = 0; time_step < 4; time_step++ )
    {
      filename_inc ( &filename );

      cout << "  " << setw(11) << filename << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test165 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST165 tests FILENAME_INC_NOWRAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  string filename;

  cout << "\n";
  cout << "TEST165\n";
  cout << "  FILENAME_INC_NOWRAP increments a string\n";
  cout << "\n";
  cout << "     Input             Output\n";

  for ( i = 0; i < 4; i++ )
  {
    if ( i == 0 )
    {
      filename = "file???.dat";
    }
    else if ( i == 1 ) 
    {
      filename = "file072.dat";
    }
    else if ( i == 2 ) 
    {
      filename = "2cat9.dat  ";
    }
    else if ( i == 3 ) 
    {
      filename = "fred98.txt ";
    }
    cout << "\n";
    for ( j = 1; j <= 4; j++ )
    {
      cout << "  " << setw(11) << filename << "  ";

      filename_inc_nowrap ( &filename );

      cout << "  " << setw(11) << filename << "\n";

      if ( s_len_trim ( filename ) <= 0 )
      {
        cout << "  (File name not incrementable.  Quit loop!)\n";
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests FILE_PARA_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "story.txt";
  int para_num;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  FILE_PARA_COUNT counts the paragraphs in a file.\n";
  cout << "\n";
  cout << "  Examining file \"" << filename << "\".\n";
  cout << "\n";

  para_num = file_para_count ( filename );
  cout << "  Number of paragraphs: " << para_num << "\n";

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests FILE_ROW_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "filum_prb_test.txt";

  cout << "\n";
  cout << "TEST22\n";
  cout << "  FILE_ROW_COUNT counts the lines in a file.\n";
  cout << "\n";
  cout << "  Examining file \"" << filename << "\".\n";
  cout << "\n";
  cout << "  Number of lines: " 
       << file_row_count ( filename ) << ".\n";

  return;
}
//****************************************************************************80

void test225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST225 tests FILE_SEQUENCE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int file_dim;
  string filename = "data_100.txt";
  int file_num;

  cout << "\n";
  cout << "TEST225\n";
  cout << "  FILE_SEQUENCE_SIZE \"sizes\" the files in a file sequence.\n";

  cout << "\n";
  cout << "  Examining files in sequence starting with:\n";
  cout << "    \"" << filename << "\".\n";

  file_sequence_size ( filename, &file_dim, &file_num );
  
  cout << "\n";
  cout << "  Number of files:      " << file_num << "\n";
  cout << "  Number of data items: " << file_dim << "\n";

  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests FILE_WORD_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "filum_prb_test.txt";

  cout << "\n";
  cout << "TEST24\n";
  cout << "  FILE_WORD_COUNT counts the words in a file.\n";
  cout << "\n";
  cout << "  Examining file \"" << filename << "\".\n";
  cout << "\n";
  cout << "  Number of words: " 
       << file_word_count ( filename ) << ".\n";

  return;
}
