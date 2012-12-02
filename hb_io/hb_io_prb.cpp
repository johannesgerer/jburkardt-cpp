# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hb_io.hpp"

int main ( );
void test01 ( string input_file );
void test02 ( );
void test03 ( string input_file );
void test04 ( );
void test05 ( string input_file );
void test06 ( );
void test07 ( string input_file );
void test08 ( );
void test09 ( string input_file );
void test10 ( );
void test11 ( );
void test12 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    HB_IO_PRB runs the HB_IO tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HB_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the HB_IO library.\n";

  test01 ( "rua_32.txt" );
  test01 ( "rse_5.txt" );
  test02 ( );
  test03 ( "rua_32.txt" );
  test03 ( "rse_5.txt" );
  test04 ( );
  test05 ( "rua_32.txt" );
  test05 ( "rse_5.txt" );
  test06 ( );
  test07 ( "rua_32.txt" );
  test08 ( );
  test09 ( "rua_32.txt" );
  test10 ( );
  test11 ( );
  test12 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HB_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( string input_file )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests HB_HEADER_READ;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int indcrd;
  char *indfmt = NULL;
  ifstream input;
  char *key = NULL;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  char *rhstyp = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;

  cout << " \n";
  cout << "TEST01\n";
  cout << "  HB_HEADER_READ reads the header of an HB file.\n";
  cout << "\n";
  cout << "  Reading the file '" << input_file << "'.\n";

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_header_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd, 
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt, 
    &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix );

  input.close ( );
//
//  Print out the  header information.
//
  hb_header_print ( title, key, totcrd, ptrcrd, indcrd, valcrd, 
    rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt, 
    rhsfmt, rhstyp, nrhs, nrhsix );

  return;
}
//****************************************************************************80

void test02 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HB_HEADER_WRITE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int indcrd = 8;
  char indfmt[17] = "(16I5)";
  char key[9] = "RUA_32";
  char mxtype[4] = "PUA";
  int ncol = 32;
  int neltvl = 0;
  int nnzero = 126;
  int nrhs = 0;
  int nrhsix = 0;
  int nrow = 32;
  ofstream output;
  string output_file = "rua_32_header.txt";
  int ptrcrd = 3;
  char ptrfmt[17] = "(16I5)";
  int rhscrd = 0;
  char rhsfmt[21] = " ";
  char rhstyp[4] = "   ";
  char title[73] = "1Real unsymmetric assembled matrix based on IBM32";
  int totcrd = 11;
  int valcrd = 0;
  char valfmt[21] = " ";

  cout << "\n";
  cout << "TEST02\n";
  cout << "  HB_HEADER_WRITE writes the header of an HB file.\n";
  cout << "\n";
  cout << "  Writing the file '" << output_file << "'.\n";

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "TEST02 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_header_write ( output, title, key, totcrd, ptrcrd, indcrd, 
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, 
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix );

  output.close ( );

  return;
}
//****************************************************************************80

void test03 ( string input_file )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests HB_STRUCTURE_READ;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int *colptr = NULL;
  int indcrd;
  char *indfmt = NULL;
  ifstream input;
  char *key = NULL;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  char *rhstyp = NULL;
  int *rowind;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  HB_STRUCTURE_READ reads the structure of an HB file.\n";
  cout << "\n";
  cout << "  Reading the file '" << input_file << "'.\n";

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "TEST03 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  cout << "  Reading the header.\n";

  hb_header_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd, 
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt, 
    &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix );

  colptr = new int[ncol+1];

  if ( mxtype[2] == 'A' )
  {
    rowind = new int[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    rowind = new int[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "TEST03 - Fatal error!\n";
    cout << "  Illegal value of MXTYPE character 3 = " << mxtype[2] << "\n";
    exit ( 1 );
  }

  cout << "  Reading the structure.\n";

  hb_structure_read ( input, ncol, mxtype, nnzero, neltvl, 
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

  input.close ( );

  cout << "\n";
  cout << title << "\n";
  cout << "  KEY =    '" << key << "'.\n";
  cout << "\n";
  cout << "  NROW =   " << nrow   << "\n";
  cout << "  NCOL =   " << ncol   << "\n";
  cout << "  NNZERO = " << nnzero << "\n";
  cout << "  NELTVL = " << neltvl << "\n";

  hb_structure_print ( ncol, mxtype, nnzero, neltvl, colptr, rowind );

  delete [] colptr;
  delete [] rowind;

  return;
}
//****************************************************************************80

void test04 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests HB_STRUCTURE_WRITE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NCOL 32
# define NELTVL 0
# define NNZERO 126

  int colptr[NCOL+1] = {
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, 
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, 
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, 
    121, 124, 127 };
  char indfmt[17] = "(16I5)";
  char mxtype[4] = "RUA";
  ofstream output;
  string output_file = "rua_32_structure.txt";
  char ptrfmt[17] = "(16I5)";
  int rowind[NNZERO] = {
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, 
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, 
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, 
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, 
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, 
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, 
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, 
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, 
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, 
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, 
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, 
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, 
   13,   14,   31,   24,   28,   32 };

  cout << "\n";
  cout << "TEST04\n";
  cout << "  HB_STRUCTURE_WRITE writes the structure of an HB file.\n";
  cout << "\n";
  cout << "  Writing the file '" << output_file << "'.\n";

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "TEST04 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_structure_write ( output, NCOL, mxtype, NNZERO, NELTVL, 
    ptrfmt, indfmt, colptr, rowind );

  output.close ( );

  return;
# undef NNZERO
# undef NELTVL
# undef NCOL
}
//****************************************************************************80

void test05 ( string input_file )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests HB_VALUES_READ;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int *colptr = NULL;
  int indcrd;
  char *indfmt = NULL;
  ifstream input;
  char *key = NULL;
  int khi;
  int klo;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  char *rhstyp = NULL;
  int *rowind = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  float *values = NULL;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  HB_VALUES_READ reads the values of an HB file.\n";
  cout << "\n";
  cout << "  Reading the file '" << input_file << "'.\n";

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "TEST05 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  cout << "  Reading the header.\n";

  hb_header_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd, 
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt, 
    &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix );

  colptr = new int[ncol+1];

  if ( mxtype[2] == 'A' )
  {
    rowind = new int[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    rowind = new int[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "TEST05 - Fatal error!\n";
    cout << "  Illegal value of MXTYPE character 3.\n";
    exit ( 1 );
  }

  cout << "  Reading the structure.\n";

  hb_structure_read ( input, ncol, mxtype, nnzero, neltvl, 
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

  if ( mxtype[2] == 'A' )
  {
    values = new float[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    values =  new float[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "TEST05 - Fatal error!\n";
    cout << "  Illegal value of MXTYPE character 3 = " << mxtype[2] << "\n";
    exit ( 1 );
  }

  cout << "  Reading the values.\n";

  hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );

  input.close ( );

  cout << "\n";
  cout << title << "\n";
  cout << "  KEY =    '" << key << "'.\n";
  cout << "\n";
  cout << "  NROW =   " << nrow   << "\n";
  cout << "  NCOL =   " << ncol   << "\n";
  cout << "  NNZERO = " << nnzero << "\n";
  cout << "  NELTVL = " << neltvl << "\n";

  hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values );

  delete [] colptr;
  delete [] rowind;
  delete [] values;

  return;
}
//****************************************************************************80

void test06 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests HB_VALUES_WRITE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NELTVL 0
# define NNZERO 126

  char mxtype[4] = "RUA";
  ofstream output;
  string output_file = "rua_32_values.txt";
  int valcrd = 13;
  char valfmt[21] = "(10F7.1)";
  float values[NNZERO] = {
  101.0,  102.0,  103.0,  104.0,  107.0, 
  126.0,  201.0,  202.0,  209.0,  221.0, 
  228.0,  302.0,  303.0,  306.0,  308.0, 
  309.0,  329.0,  403.0,  404.0,  405.0, 
  412.0,  503.0,  505.0,  523.0,  527.0, 
  601.0,  606.0,  616.0,  703.0,  707.0, 
  714.0,  721.0,  731.0,  801.0,  808.0, 
  812.0,  817.0,  827.0,  907.0,  909.0, 
  910.0,  913.0,  919.0,  923.0,  927.0, 
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, 
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, 
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, 
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, 
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, 
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, 
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, 
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, 
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, 
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, 
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, 
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, 
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, 
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, 
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, 
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, 
 3232.0 };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  HB_VALUES_WRITE writes the values of an HB file.\n";
  cout << "\n";
  cout << "  Writing the file '" << output_file << "'.\n";

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "TEST06 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_values_write ( output, valcrd, mxtype, NNZERO, NELTVL, valfmt, values );

  output.close ( );

  return;
# undef NELTVL
# undef NNZERO
}
//****************************************************************************80

void test07 ( string input_file )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests HB_RHS_READ, HB_GUESS_READ, HB_EXACT_READ;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int *colptr = NULL;
  float *exact = NULL;
  float *guess = NULL;
  int i;
  int indcrd;
  char *indfmt = NULL;
  ifstream input;
  int j;
  char *key = NULL;
  int khi;
  int klo;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  int *rhsind = NULL;
  int *rhsptr = NULL;
  char *rhstyp = NULL;
  float *rhsval = NULL;
  float *rhsvec = NULL;
  int *rowind = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  float *values = NULL;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  HB_RHS_READ reads right hand sides from an HB file.\n";
  cout << "  HB_GUESS_READ reads starting guesses from an HB file.\n";
  cout << "  HB_EXACT_READ reads exact solutions from an HB file.\n";

  cout << "\n";
  cout << "  Reading the file '" << input_file << "'.\n";

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "TEST07 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  cout << "  Reading the header.\n";

  hb_header_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd, 
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt, 
    &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix );

  colptr = new int[ncol+1];

  if ( mxtype[2] == 'A' )
  {
    rowind = new int[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    rowind = new int[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "TEST07 - Fatal error!\n";
    cout << "  Illegal value of MXTYPE character 3 = " << mxtype[2] << "\n";
    exit ( 1 );
  }

  cout << "  Reading the structure.\n";

  hb_structure_read ( input, ncol, mxtype, nnzero, neltvl, 
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

  if ( mxtype[2] == 'A' )
  {
    values = new float[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    values =  new float[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "TEST07 - Fatal error!\n";
    cout << "  Illegal value of MXTYPE character 3 = " << mxtype[2] << "\n";
    exit ( 1 );
  }

  cout << "  Reading the values.\n";

  hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );
//
//  Read the right hand sides.
//
  if ( 0 < rhscrd )
  {
    cout << "  Reading the right hand side.\n";

    if ( rhstyp[0] == 'F' )
    {
      rhsval = new float[nrow*nrhs];
    }
    else if ( rhstyp[0] == 'M' )
    {
      if ( mxtype[2] == 'A' )
      {
        rhsptr = new int[nrhs+1];
        rhsind = new int [nrhsix];
        rhsvec = new float[nrhsix];
      }
      else if ( mxtype[2] == 'E' )
      {
        rhsval = new float [nnzero*nrhs];
      }
    }

    hb_rhs_read ( input, nrow, nnzero, nrhs, nrhsix, 
      rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, 
      rhsind, rhsptr, rhsvec );

    cout << "  Done reading the right hand side.\n";

    if ( rhstyp[0] == 'F' )
    {
      r4mat_print_some ( nrow, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS" );
    }
    else if ( rhstyp[0] == 'M' && mxtype[2] == 'A' )
    {
      i4vec_print_part ( nrhs+1, rhsptr, 10, "  Part of RHSPTR" );
      i4vec_print_part ( nrhsix, rhsind, 10, "  Part of RHSIND" );
      r4vec_print_part ( nrhsix, rhsvec, 10, "  Part of RHSVEC" );
    }
    else if ( rhstyp[0] == 'M' && mxtype[2] == 'E' )
    {
      r4mat_print_some ( nnzero, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS" );
    }
//
//  Read the starting guesses.
//
    if ( rhstyp[1] == 'G' )
    {
      cout << "  Reading the starting guesses.\n";

      guess = new float[nrow*nrhs];

      hb_guess_read ( input, nrow, nrhs, rhscrd, rhsfmt, rhstyp, guess );

      r4mat_print_some ( nrow, nrhs, guess, 1, 1, 5, 5, "  Part of GUESS" );
    }
//
//  Read the exact solutions.
//
    if ( rhstyp[2] == 'X' )
    {
      cout << "  Reading the exact solutions.\n";

      exact = new float[nrow*nrhs];

      hb_exact_read ( input, nrow, nrhs, rhscrd, rhsfmt, rhstyp, exact );

      r4mat_print_some ( nrow, nrhs, exact, 1, 1, 5, 5, "  Part of EXACT" );
    }
  }

  input.close ( );

  if ( colptr )
  {
    delete [] colptr;
  }
  if ( exact )
  {
    delete [] exact;
  }
  if ( guess )
  {
    delete [] guess;
  }
  if ( rhsind )
  {
    delete [] rhsind;
  }
  if ( rhsptr )
  {
    delete [] rhsptr;
  }
  if ( rhsval )
  {
    delete [] rhsval;
  }
  if ( rhsvec )
  {
    delete [] rhsvec;
  }
  if ( rowind )
  {
    delete [] rowind;
  }
  if ( values )
  {
    delete [] values;
  }

  return;
}
//****************************************************************************80

void test08 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests HB_RHS_WRITE, HB_GUESS_WRITE, HB_EXACT_WRITE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NELTVL 0
# define NNZERO 126
# define NRHS 1
# define NRHSIX 0
# define NROW 32

  float exact[NROW*NRHS] = {
    1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0, 
   11.0,  12.0,  13.0,  14.0,  15.0,  16.0,  17.0,  18.0,  19.0,  20.0, 
   21.0,  22.0,  23.0,  24.0,  25.0,  26.0,  27.0,  28.0,  29.0,  30.0, 
   31.0,  32.0 };
  float guess[NROW*NRHS] = {
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0 };
  char indfmt[17] = "(16I5)";
  char mxtype[4] = "RUA";
  ofstream output;
  string output_file = "rua_32_rhs.txt";
  char ptrfmt[17] = "(16I5)";
  int rhscrd = 12;
  char rhsfmt[21] = "(10F7.1)";
  int *rhsind = NULL;
  int *rhsptr = NULL;
  float rhsval[NROW*NRHS] = {
    101.0, 102.0, 103.0, 104.0, 107.0, 126.0, 201.0, 202.0, 209.0, 221.0, 
    228.0, 302.0, 303.0, 306.0, 308.0, 309.0, 329.0, 403.0, 404.0, 405.0, 
    412.0, 503.0, 505.0, 523.0, 527.0, 601.0, 606.0, 616.0, 703.0, 707.0, 
    714.0, 721.0 };
  float *rhsvec = NULL;
  char rhstyp[4] = "FGX";

  cout << "\n";
  cout << "TEST08\n";
  cout << "  HB_RHS_WRITE writes the right hand sides to an HB file.\n";
  cout << "  HB_GUESS_WRITE writes starting guesses to an HB file.\n";
  cout << "  HB_EXACT_WRITE writes exact solutions to an HB file.\n";
  cout << "\n";
  cout << "  Writing the file '" << output_file << "'.\n";

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "TEST08 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }
//
//  Write the right hand sides.
//
  hb_rhs_write ( output, NROW, NNZERO, NRHS, NRHSIX, 
    rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, 
    rhsind, rhsptr, rhsvec );
//
//  Write the right hand sides.
//
  hb_guess_write ( output, NROW, NRHS, rhscrd, rhsfmt, rhstyp, guess );
//
//  Write the right hand sides.
//
  hb_exact_write ( output, NROW, NRHS, rhscrd, rhsfmt, rhstyp, exact );

  output.close ( );

  return;
# undef NELTVL
# undef NNZERO
# undef NRHS
# undef NRHSIX
# undef NROW
}
//****************************************************************************80

void test09 ( string input_file )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests HB_FILE_READ;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int *colptr = NULL;
  float *exact = NULL;
  float *guess = NULL;
  int i;
  int indcrd;
  char *indfmt = NULL;
  ifstream input;
  int j;
  char *key = NULL;
  int khi;
  int klo;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  int *rhsind = NULL;
  int *rhsptr = NULL;
  char *rhstyp = NULL;
  float *rhsval = NULL;
  float *rhsvec = NULL;
  int *rowind = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  float *values = NULL;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  HB_FILE_READ reads all the data in an HB file.\n";
  cout << "  HB_FILE_MODULE is the module that stores the data.\n";

  cout << "\n";
  cout << "  Reading the file '" << input_file << "'.\n";

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "TEST09 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_file_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd, 
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, 
    &ptrfmt, &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix, 
    &colptr, &rowind, &values, &rhsval, &rhsptr, &rhsind, &rhsvec, 
    &guess, &exact );
//
//  Print out the header information.
//
  hb_header_print ( title, key, totcrd, ptrcrd, indcrd, valcrd, 
    rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt, 
    rhsfmt, rhstyp, nrhs, nrhsix );
//
//  Print the structure information.
//
  hb_structure_print ( ncol, mxtype, nnzero, neltvl, colptr, rowind );
//
//  Print the values.
//
  hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values );

  if ( 0 < rhscrd )
  {
//
//  Print a bit of the right hand sides.
//
    if ( rhstyp[0] == 'F' )
    {
      r4mat_print_some ( nrow, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS" );
    }
    else if ( rhstyp[0] == 'M' && mxtype[2] == 'A' )
    {
      i4vec_print_part ( nrhs+1, rhsptr, 10, "  Part of RHSPTR" );
      i4vec_print_part ( nrhsix, rhsind, 10, "  Part of RHSIND" );
      r4vec_print_part ( nrhsix, rhsvec, 10, "  Part of RHSVEC" );
    }
    else if ( rhstyp[0] == 'M' && mxtype[2] == 'E' )
    {
      r4mat_print_some ( nnzero, nrhs, rhsval, 1, 1, 5, 5, "  Part of RHS" );
    }
//
//  Print a bit of the starting guesses.
//
    if ( rhstyp[1] == 'G' )
    {
      r4mat_print_some ( nrow, nrhs, guess, 1, 1, 5, 5, "  Part of GUESS" );
    }
//
//  Print a bit of the exact solutions.
//
    if ( rhstyp[2] == 'X' )
    {
      r4mat_print_some ( nrow, nrhs, exact, 1, 1, 5, 5, "  Part of EXACT" );
    }

  }

  input.close ( );

  if ( colptr )
  {
    delete [] colptr;
  }
  if ( exact )
  {
    delete [] exact;
  }
  if ( guess )
  {
    delete [] guess;
  }
  if ( rhsind )
  {
    delete [] rhsind;
  }
  if ( rhsptr )
  {
    delete [] rhsptr;
  }
  if ( rhsval )
  {
    delete [] rhsval;
  }
  if ( rhsvec )
  {
    delete [] rhsvec;
  }
  if ( rowind )
  {
    delete [] rowind;
  }
  if ( values )
  {
    delete [] values;
  }

  return;
}
//****************************************************************************80

void test10 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests HB_FILE_WRITE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NCOL 32
# define NELTVL 0
# define NNZERO 126
# define NRHS 1
# define NRHSIX 0
# define NROW 32

  int colptr[NCOL+1] = {
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, 
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, 
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, 
    121, 124, 127 };
  float exact[NROW*NRHS] = {
    1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0, 
   11.0,  12.0,  13.0,  14.0,  15.0,  16.0,  17.0,  18.0,  19.0,  20.0, 
   21.0,  22.0,  23.0,  24.0,  25.0,  26.0,  27.0,  28.0,  29.0,  30.0, 
   31.0,  32.0 };
  float guess[NROW*NRHS] = {
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0 };
  int indcrd = 8;
  char indfmt[17] = "(16I5)";
  char key[9] = "RUA_32";
  char mxtype[4] = "RUA";
  string output_file = "rua_32_file.txt";
  ofstream output;
  int ptrcrd = 3;
  char ptrfmt[17] = "(16I5)";
  int rhscrd = 12;
  char rhsfmt[21] = "(10F7.1)";
  int *rhsind = NULL;
  int *rhsptr = NULL;
  float rhsval[NROW*NRHS] = {
    101.0, 102.0, 103.0, 104.0, 107.0, 126.0, 201.0, 202.0, 209.0, 221.0, 
    228.0, 302.0, 303.0, 306.0, 308.0, 309.0, 329.0, 403.0, 404.0, 405.0, 
    412.0, 503.0, 505.0, 523.0, 527.0, 601.0, 606.0, 616.0, 703.0, 707.0, 
    714.0, 721.0 };
  char rhstyp[4] = "FGX";
  float *rhsvec = NULL;
  int rowind[NNZERO] = {
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, 
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, 
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, 
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, 
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, 
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, 
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, 
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, 
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, 
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, 
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, 
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, 
   13,   14,   31,   24,   28,   32 };
  char title[73] = "1Real unsymmetric assembled matrix based on IBM32";
  int totcrd = 36;
  int valcrd = 13;
  char valfmt[21] = "(10F7.1)";
  float values[NNZERO] = { 
  101.0,  102.0,  103.0,  104.0,  107.0, 
  126.0,  201.0,  202.0,  209.0,  221.0, 
  228.0,  302.0,  303.0,  306.0,  308.0, 
  309.0,  329.0,  403.0,  404.0,  405.0, 
  412.0,  503.0,  505.0,  523.0,  527.0, 
  601.0,  606.0,  616.0,  703.0,  707.0, 
  714.0,  721.0,  731.0,  801.0,  808.0, 
  812.0,  817.0,  827.0,  907.0,  909.0,
  910.0,  913.0,  919.0,  923.0,  927.0, 
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, 
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, 
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, 
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, 
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, 
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, 
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, 
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, 
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, 
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, 
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, 
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, 
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, 
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, 
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, 
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, 
 3232.0 };

  cout << "\n";
  cout << "TEST10\n";
  cout << "  HB_FILE_WRITE writes an HB file.\n";
  cout << "\n";
  cout << "  Writing the file '" << output_file << "'.\n";

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "TEST10 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_file_write ( output, title, key, totcrd, ptrcrd, indcrd, 
    valcrd, rhscrd, mxtype, NROW, NCOL, NNZERO, NELTVL, ptrfmt, indfmt, 
    valfmt, rhsfmt, rhstyp, NRHS, NRHSIX, colptr, rowind, values, 
    rhsval, rhsptr, rhsind, rhsvec, guess, exact );

  output.close ( );

  return;
# undef NCOL
# undef NELTVL
# undef NNZERO
# undef NRHS
# undef NRHSIX
# undef NROW
}
//****************************************************************************80

void test11 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests HB_FILE_WRITE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NCOL 32
# define NELTVL 0
# define NNZERO 126
# define NRHS 2
# define NRHSIX 0
# define NROW 32

  int colptr[NCOL+1] = {
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, 
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, 
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, 
    121, 124, 127 };
  float exact[NROW*NRHS] = {
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0, 
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 
    0.0,   0.0,
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0 };
  float guess[NROW*NRHS] = {
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0 };
  int indcrd = 8;
  char indfmt[17] = "(16I5)";
  char key[9] = "RUA_32";
  char mxtype[4] = "RUA";
  string output_file = "rua_32_ax.txt";
  ofstream output;
  int ptrcrd = 3;
  char ptrfmt[17] = "(16I5)";
  int rhscrd = 12;
  char rhsfmt[21] = "(10F7.1)";
  int *rhsind = NULL;
  int *rhsptr = NULL;
  float *rhsval;
  char rhstyp[4] = "FGX";
  float *rhsvec = NULL;
  int rowind[NNZERO] = {
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, 
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, 
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, 
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, 
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, 
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, 
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, 
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, 
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, 
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, 
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, 
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, 
   13,   14,   31,   24,   28,   32 };
  char title[73] = "1Real unsymmetric assembled matrix based on IBM32";
  int totcrd = 36;
  int valcrd = 13;
  char valfmt[21] = "(10F7.1)";
  float values[NNZERO] = { 
  101.0,  102.0,  103.0,  104.0,  107.0, 
  126.0,  201.0,  202.0,  209.0,  221.0, 
  228.0,  302.0,  303.0,  306.0,  308.0, 
  309.0,  329.0,  403.0,  404.0,  405.0, 
  412.0,  503.0,  505.0,  523.0,  527.0, 
  601.0,  606.0,  616.0,  703.0,  707.0, 
  714.0,  721.0,  731.0,  801.0,  808.0, 
  812.0,  817.0,  827.0,  907.0,  909.0,
  910.0,  913.0,  919.0,  923.0,  927.0, 
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, 
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, 
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, 
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, 
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, 
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, 
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, 
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, 
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, 
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, 
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, 
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, 
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, 
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, 
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, 
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, 
 3232.0 };

  cout << "\n";
  cout << "TEST11\n";
  cout << "  HB_MATVEC_A_MEM multiplies a matrix times a vector.\n";
  cout << "\n";
  cout << "  This particular version assumes:\n";
  cout << "  * the matrix is in ""A"" format (assembled),\n";
  cout << "  * the matrix and vectors can fit in memory,\n";
  cout << "  * the matrix and multiplicand have been read into\n";
  cout << "    memory before the routine is called.\n";
  cout << "\n";
  cout << "  For this example, the first vector X is zero except\n";
  cout << "  for a 1 in row 10.  This means A*X should return\n";
  cout << "  column 10 of A.\n";
  cout << "\n";
  cout << "  The second vector X is all 1's.  A*X should be\n";
  cout << "  the sum of the entries of each row.\n";

  rhsval = hb_matvec_a_mem ( NROW, NCOL, NNZERO, NRHS, colptr, rowind, values,
    exact );

  r4mat_print ( NROW, NRHS, rhsval,  "  The product vectors A*X" );

  cout << "\n";
  cout << "  Writing the file '" << output_file << "'.\n";

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "TEST11 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_file_write ( output, title, key, totcrd, ptrcrd, indcrd, 
    valcrd, rhscrd, mxtype, NROW, NCOL, NNZERO, NELTVL, ptrfmt, indfmt, 
    valfmt, rhsfmt, rhstyp, NRHS, NRHSIX, colptr, rowind, values, 
    rhsval, rhsptr, rhsind, rhsvec, guess, exact );

  output.close ( );

  delete [] rhsval;

  return;
# undef NCOL
# undef NELTVL
# undef NNZERO
# undef NRHS
# undef NRHSIX
# undef NROW
}
//****************************************************************************80

void test12 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests HB_FILE_WRITE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NCOL 32
# define NELTVL 0
# define NNZERO 126
# define NRHS 2
# define NRHSIX 0
# define NROW 32

  int colptr[NCOL+1] = {
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, 
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, 
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, 
    121, 124, 127 };
  float exact[NROW*NRHS] = {
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0, 
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 
    0.0,   0.0,
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0 };
  float guess[NROW*NRHS] = {
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, 
    1.0,   1.0 };
  int indcrd = 8;
  char indfmt[17] = "(16I5)";
  char key[9] = "RUA_32";
  char mxtype[4] = "RUA";
  string output_file = "rua_32_xa.txt";
  ofstream output;
  int ptrcrd = 3;
  char ptrfmt[17] = "(16I5)";
  int rhscrd = 12;
  char rhsfmt[21] = "(10F7.1)";
  int *rhsind = NULL;
  int *rhsptr = NULL;
  float *rhsval;
  char rhstyp[4] = "FGX";
  float *rhsvec = NULL;
  int rowind[NNZERO] = {
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, 
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, 
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, 
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, 
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, 
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, 
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, 
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, 
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, 
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, 
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, 
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, 
   13,   14,   31,   24,   28,   32 };
  char title[73] = "1Real unsymmetric assembled matrix based on IBM32";
  int totcrd = 36;
  int valcrd = 13;
  char valfmt[21] = "(10F7.1)";
  float values[NNZERO] = { 
  101.0,  102.0,  103.0,  104.0,  107.0, 
  126.0,  201.0,  202.0,  209.0,  221.0, 
  228.0,  302.0,  303.0,  306.0,  308.0, 
  309.0,  329.0,  403.0,  404.0,  405.0, 
  412.0,  503.0,  505.0,  523.0,  527.0, 
  601.0,  606.0,  616.0,  703.0,  707.0, 
  714.0,  721.0,  731.0,  801.0,  808.0, 
  812.0,  817.0,  827.0,  907.0,  909.0,
  910.0,  913.0,  919.0,  923.0,  927.0, 
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, 
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, 
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, 
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, 
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, 
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, 
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, 
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, 
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, 
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, 
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, 
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, 
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, 
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, 
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, 
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, 
 3232.0 };

  cout << "\n";
  cout << "TEST12\n";
  cout << "  HB_VECMAT_A_MEM multiplies a vector times a matrix.\n";
  cout << "\n";
  cout << "  This particular version assumes:\n";
  cout << "  * the matrix is in ""A"" format (assembled),\n";
  cout << "  * the matrix and vectors can fit in memory,\n";
  cout << "  * the matrix and multiplicand have been read into\n";
  cout << "    memory before the routine is called.\n";
  cout << "\n";
  cout << "  For this example, the first vector X is zero except\n";
  cout << "  for a 1 in row 10.  This means A'*X should return\n";
  cout << "  row 10 of A.\n";
  cout << "\n";
  cout << "  The second vector X is all 1's.  A'*X should be\n";
  cout << "  the sum of the entries of each column.\n";

  rhsval = hb_vecmat_a_mem ( NROW, NCOL, NNZERO, NRHS, colptr, rowind, values,
    exact );

  r4mat_print ( NCOL, NRHS, rhsval,  "  The product vectors A'*X" );

  cout << "\n";
  cout << "  Writing the file '" << output_file << "'.\n";

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "TEST12 - Fatal error!\n";
    cout << "  Error opening the file.\n";
    return;
  }

  hb_file_write ( output, title, key, totcrd, ptrcrd, indcrd, 
    valcrd, rhscrd, mxtype, NROW, NCOL, NNZERO, NELTVL, ptrfmt, indfmt, 
    valfmt, rhsfmt, rhstyp, NRHS, NRHSIX, colptr, rowind, values, 
    rhsval, rhsptr, rhsind, rhsvec, guess, exact );

  output.close ( );

  delete [] rhsval;

  return;
# undef NCOL
# undef NELTVL
# undef NNZERO
# undef NRHS
# undef NRHSIX
# undef NROW
}
