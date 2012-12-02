# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "cnf_io.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CNF_IO_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CNF_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CNF_IO library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CNF_IO_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 calls CNF_WRITE to write a small CNF example to a CNF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define C_NUM 2
# define L_NUM 5

  int c_num = C_NUM;
  int l_num = L_NUM;

  string cnf_file_name = "cnf_io_v3_c2.cnf";
  bool error;
  int l_c_num[C_NUM] = { 2, 3 };
  int l_val[L_NUM] = { 1, -3, 2, 3, -1 };
  int v_num = 3;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  CNF_WRITE can write CNF data to a CNF file.\n";
  cout << "\n";
  cout << "  Here is the data:\n";

  cnf_print ( v_num, c_num, l_num, l_c_num, l_val );

  cout << "\n";
  cout << "  Now we call CNF_WRITE to store this information\n";
  cout << "  in the file \"" << cnf_file_name << "\".\n";

  error = cnf_write ( v_num, c_num, l_num, l_c_num, l_val, cnf_file_name );

  return;
# undef C_NUM
# undef L_NUM
}
//****************************************************************************80

void test02 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 calls CNF_HEADER_READ to read the header of a small CNF example file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2008
//
//  Author:
//
//    John Burkardt
//
{
  int c_num;
  string cnf_file_name = "cnf_io_v3_c2.cnf";
  bool error;
  int l_num;
  int v_num;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  CNF_HEADER_READ reads the header of a CNF file.\n";

  cout << "\n";
  cout << "  Read the header of \"" << cnf_file_name << "\".\n";

  error = cnf_header_read ( cnf_file_name, &v_num, &c_num, &l_num );

  if ( error )
  {
    cout << "\n";
    cout << "  The header information could not be read.\n";
    return;
  }

  cout << "\n";
  cout << "  The number of variables       V_NUM  = " << v_num << "\n";
  cout << "  The number of clauses         C_NUM  = " << c_num << "\n";
  cout << "  The number of signed literals L_NUM  = " << l_num << "\n";

  return;
}
//****************************************************************************80

void test03 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 calls CNF_DATA_READ to read the data of a small CNF example file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2008
//
//  Author:
//
//    John Burkardt
//
{
  int c_num;
  string cnf_file_name = "cnf_io_v3_c2.cnf";
  bool error;
  int *l_c_num;
  int l_num;
  int *l_val;
  int v_num;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  CNF_DATA_READ reads the data of a CNF file.\n";

  cout << "\n";
  cout << "  Read the header of \"" << cnf_file_name << "\".\n";

  error = cnf_header_read ( cnf_file_name, &v_num, &c_num, &l_num );

  if ( error )
  {
    cout << "\n";
    cout << "  The header information could not be read.\n";
    return;
  }

  cout << "\n";
  cout << "  The number of variables       V_NUM  = " << v_num << "\n";
  cout << "  The number of clauses         C_NUM  = " << c_num << "\n";
  cout << "  The number of signed literals L_NUM  = " << l_num << "\n";

  l_c_num = new int[c_num];
  l_val = new int[l_num];

  cnf_data_read ( cnf_file_name, v_num, c_num, l_num, l_c_num, l_val );

  cout << "\n";
  cout << "  Here is the data as read from the file:\n";

  cnf_print ( v_num, c_num, l_num, l_c_num, l_val );

  delete [] l_c_num;
  delete [] l_val;

  return;
}
//****************************************************************************80

void test04 ( void )

//*****************************************************************************80
//
//  Purpose:
//
//    TEST04 calls CNF_WRITE to write a CNF example to a CNF file.
//
//  Discussion:
//
//    This formula is used as an example in the Quinn reference.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define C_NUM 18
# define L_NUM 36

  int c_num = C_NUM;
  int l_num = L_NUM;

  string cnf_file_name ="cnf_io_v16_c18.cnf";
  bool error;
  int l_c_num[C_NUM] = {
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2 };
  int l_val[L_NUM] = {
    1,    2, 
   -2,   -4, 
    3,    4, 
   -4,   -5, 
    5,   -6, 
    6,   -7, 
    6,    7, 
    7,  -16, 
    8,   -9, 
   -8,  -14, 
    9,   10, 
    9,  -10, 
  -10,  -11, 
   10,   12, 
   11,   12, 
   13,   14, 
   14,  -15, 
   15,   16 };
  int v_num = 16;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  CNF_WRITE can write CNF data to a CNF file.\n";
  cout << "\n";
  cout << "  Here is the data to be written to the file:\n";

  cnf_print ( v_num, c_num, l_num, l_c_num, l_val );

  cout << "\n";
  cout << "  Now we call CNF_WRITE to store this information\n";
  cout << "  in the file \"" << cnf_file_name << "\".\n";

  error = cnf_write ( v_num, c_num, l_num, l_c_num, l_val, cnf_file_name );

  return;
# undef C_NUM
# undef L_NUM
}
//****************************************************************************80

void test05 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 calls CNF_HEADER_READ to read the header of a small CNF example file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2008
//
//  Author:
//
//    John Burkardt
//
{
  int c_num;
  string cnf_file_name = "cnf_io_v16_c18.cnf";
  bool error;
  int l_num;
  int v_num;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  CNF_HEADER_READ reads the header of a CNF file.\n";

  cout << "\n";
  cout << "  Read the header of \"" << cnf_file_name << "\".\n";

  error = cnf_header_read ( cnf_file_name, &v_num, &c_num, &l_num );

  if ( error )
  {
    cout << "\n";
    cout << "  The header information could not be read.\n";
  }
  else
  {
    cout << "\n";
    cout << "  The number of variables       V_NUM  = " << v_num << "\n";
    cout << "  The number of clauses         C_NUM  = " << c_num << "\n";
    cout << "  The number of signed literals L_NUM  = " << l_num << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 calls CNF_DATA_READ to read the data of a small CNF example file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2008
//
//  Author:
//
//    John Burkardt
//
{
  int c_num;
  string cnf_file_name = "cnf_io_v16_c18.cnf";
  bool error;
  int *l_c_num;
  int l_num;
  int *l_val;
  int v_num;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  CNF_DATA_READ reads the data of a CNF file.\n";

  cout << "\n";
  cout << "  Read the header of \"" << cnf_file_name << "\".\n";

  error = cnf_header_read ( cnf_file_name, &v_num, &c_num, &l_num );

  if ( error )
  {
    cout << "\n";
    cout << "  The header information could not be read.\n";
    return;
  }

  cout << "\n";
  cout << "  The number of variables       V_NUM  = " << v_num << "\n";
  cout << "  The number of clauses         C_NUM  = " << c_num << "\n";
  cout << "  The number of signed literals L_NUM  = " << l_num << "\n";

  l_c_num = new int[c_num];
  l_val = new int[l_num];

  cnf_data_read ( cnf_file_name, v_num, c_num, l_num, l_c_num, l_val );

  cout << "\n";
  cout << "  Here is the data as read from the file:\n";

  cnf_print ( v_num, c_num, l_num, l_c_num, l_val );

  delete [] l_c_num;
  delete [] l_val;

  return;
}
//****************************************************************************80

void test07 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 calls CNF_EVALUATE to evaluate a formula.
//
//  Discussion:
//
//    This formula is used as an example in the Quinn reference.
//    Here, we seek the logical inputs that make the formula true.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define C_NUM 18
# define L_NUM 36
# define V_NUM 16

  int c_num = C_NUM;
  int l_num = L_NUM;
  int v_num = V_NUM;

  bool f_val;
  int i;
  int ihi;
  int j;
  int l_c_num[C_NUM] = {
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2 };
  int l_val[L_NUM] = {
    1,    2, 
   -2,   -4, 
    3,    4, 
   -4,   -5, 
    5,   -6, 
    6,   -7, 
    6,    7, 
    7,  -16, 
    8,   -9, 
   -8,  -14, 
    9,   10, 
    9,  -10, 
  -10,  -11, 
   10,   12, 
   11,   12, 
   13,   14, 
   14,  -15, 
   15,   16 };
  int solution_num;
  bool v_val[V_NUM];

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Seek inputs to a circuit that produce a 1 (TRUE) output.\n";
  cout << "\n";
  cout << "  Here is the CNF data defining the formula:\n";

  cnf_print ( v_num, c_num, l_num, l_c_num, l_val );
//
//  Initialize the logical vector.
//
  for ( j = 0; j < v_num; j++ )
  {
    v_val[j] = false;
  }
//
//  Compute the number of binary vectors to check.
//
  ihi = i4_power ( 2, v_num );

  cout << "\n";
  cout << "  Number of input vectors to check is " << ihi << "\n";
  cout << "\n";
  cout << "   #       Index    ---------Input Values----------\n";
  cout << "\n";
//
//  Check every possible input vector.
//
  solution_num = 0;

  for ( i = 0; i < ihi; i++ )
  {
    f_val = cnf_evaluate ( v_num, c_num, l_num, l_c_num, l_val, v_val );

    if ( f_val )
    {
      solution_num = solution_num + 1;
      cout << "  " << setw(2) << solution_num
           << "  " << setw(10) << i;
      for ( j = 0; j < v_num; j++ )
      {
        cout << v_val[j];
      }
      cout << "\n";
    }
    lvec_next ( v_num, v_val );
  }
//
//  Report.
//
  cout << "\n";
  cout << "  Number of solutions found was " << solution_num << "\n";

  return;
# undef C_NUM
# undef L_NUM
# undef V_NUM
}
