# include <stdlib.h>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "xyz_io.hpp"

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for XYZ_IO_PRB.
//
//  Discussion:
//
//    XYZ_IO_PRB tests the XYZ_IO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "XYZ_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the XYZ_IO library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "XYZ_IO_PRB:\n";
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
//   TEST01 tests XYZ_EXAMPLE, XYZ_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  string file_name = "helix.xyz";
  int point_num;
  double *xyz;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  XYZ_EXAMPLE sets up sample XYZ data.\n";
  cout << "  XYZ_WRITE writes an XYZ file.\n";

  point_num = xyz_example_size ( );

  cout << "  Example dataset size is " << point_num << "\n";

  xyz = new double[3*point_num];

  xyz_example ( point_num, xyz );

  xyz_write ( file_name, point_num, xyz );

  delete [] xyz;

  cout << "\n";
  cout << "  XYZ_WRITE wrote the header and data for \"" << file_name << "\".\n";
  cout << "  Number of points = " << point_num << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests XYZ_HEADER_READ, XYZ_DATA_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  string file_name = "xyz_io_prb_02.xyz";
  int i;
  int k;
  int point_num;
  double *xyz;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  XYZ_HEADER_READ reads the header of an XYZ file.\n";
  cout << "  XYZ_DATA_READ reads the data of an XYZ file.\n";

  xyz_write_test ( file_name );

  cout << "\n";
  cout << "  XYZ_WRITE_TEST created some data.\n";

  xyz_header_read ( file_name, &point_num );

  cout << "\n";
  cout << "  XYZ_HEADER_READ has read the header.\n";
  cout << "\n";
  cout << "  Number of points = " << point_num << "\n";

  xyz = new double[3*point_num];

  xyz_data_read ( file_name, point_num, xyz );

  cout << "\n";
  cout << "  XYZ_DATA_READ has read the data.\n";

  cout << "\n";
  cout << "  Sample data:\n";
  cout << "\n";

  for ( k =  1; k <= 11; k++ )
  {
    i = ( ( 11 - k     ) * 1 
        + (      k - 1 ) * point_num ) 
        / ( 11     - 1 ) - 1;
    cout << "  " << setw(4) << i
         << "  " << setw(12) << xyz[0+i*3]
         << "  " << setw(12) << xyz[1+i*3]
         << "  " << setw(12) << xyz[2+i*3] << "\n";
  }

  delete [] xyz;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests XYZL_EXAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int *line_data;
  int *line_pointer;
  int line_data_num;
  int line_num;
  int point_num;
  double *xyz;
  string xyz_filename = "cube.xyz";
  string xyzl_filename = "cube.xyzl";

  cout << "\n";
  cout << "TEST03\n";
  cout << "  XYZL_EXAMPLE sets up XYZ and XYZL data.\n";

  xyzl_example_size ( &point_num, &line_num, &line_data_num );

  cout << "\n";
  cout << "  Example has:\n";
  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";
  cout << "  Number of lines      = " << line_num << "\n";
  cout << "  Number of line items = " << line_data_num << "\n";

  line_data = new int[line_data_num];
  line_pointer = new int[line_num+1];
  xyz = new double[3*point_num];

  xyzl_example ( point_num, line_num, line_data_num, xyz, line_pointer, 
    line_data );

  xyz_write ( xyz_filename, point_num, xyz );

  xyzl_write ( xyzl_filename, point_num, line_num, line_data_num, 
    line_pointer, line_data );

  cout << "\n";
  cout << "  Wrote the XYZ file \"" << xyz_filename << "\".\n";
  cout << "  and the XYZL file \"" << xyzl_filename << "\".\n";

  delete [] line_data;
  delete [] line_pointer;
  delete [] xyz;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests XYZL_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int line;
  int *line_data;
  int *line_pointer;
  int line_data_num;
  int line_num;
  int point_num;
  double *xyz;
  string xyz_filename = "cube.xyz";
  string xyzl_filename = "cube.xyzl";

  cout << "\n";
  cout << "TEST04\n";
  cout << "  XYZ_HEADER_READ  reads the header of an XYZ  file.\n";
  cout << "  XYZ_DATA_READ    reads the data   of an XYZ  file.\n";
  cout << "  XYZL_HEADER_READ reads the header of an XYZL file.\n";
  cout << "  XYZL_DATA_READ   reads the data   of an XYZL file.\n";

  cout << "\n";
  cout << "  Examine XYZ file \"" << xyz_filename << "\".\n";

  xyz_header_read ( xyz_filename, &point_num );

  cout << "\n";
  cout << "  Number of points = " << point_num << "\n";

  xyz = new double[3*point_num];

  xyz_data_read ( xyz_filename, point_num, xyz );

  cout << "\n";
  cout << "  Point data:\n";
  cout << "\n";

  for ( i = 0; i < point_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << xyz[0+i*3]
         << "  " << setw(10) << xyz[1+i*3]
         << "  " << setw(10) << xyz[2+i*3] << "\n";
  }

  cout << "\n";
  cout << "  Examine XYZL file \"" << xyzl_filename << "\".\n";

  xyzl_header_read ( xyzl_filename, &line_num, &line_data_num );

  cout << "\n";
  cout << "  Number of lines      = " << line_num << "\n";
  cout << "  Number of line items = " << line_data_num << "\n";

  line_data = new int[line_data_num];
  line_pointer = new int[line_num+1];

  xyzl_data_read ( xyzl_filename, line_num, line_data_num, line_pointer, 
    line_data );

  cout << "\n";
  cout << "  Line pointers:\n";
  cout << "\n";

  for ( line = 0; line < line_num; line++ )
  {
    cout << "  " << setw(4) << line
         << "  " << setw(8) << line_pointer[line]
         << "  " << setw(8) << line_pointer[line+1] - 1 << "\n";
  }

  cout << "\n";
  cout << "  Line data:\n";
  cout << "\n";

  for ( line = 0; line < line_num; line++ )
  {
    cout << "  " << setw(4) << line << "    ";
    for ( j = line_pointer[line]; j <= line_pointer[line+1] - 1; j++ )
    {
      cout << "  " << setw(8) << line_data[j];
    }
    cout << "\n";
  }

  delete [] line_data;
  delete [] line_pointer;
  delete [] xyz;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests XYZF_EXAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *face_data;
  int *face_pointer;
  int face_data_num;
  int face_num;
  int i;
  int point_num;
  double *xyz;
  string xyz_filename = "cube.xyz";
  string xyzf_filename = "cube.xyzf";

  cout << "\n";
  cout << "TEST05\n";
  cout << "  XYZF_EXAMPLE sets up XYZ and XYZF data.\n";

  xyzf_example_size ( &point_num, &face_num, &face_data_num );

  cout << "\n";
  cout << "  Example has:\n";
  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";
  cout << "  Number of faces      = " << face_num << "\n";
  cout << "  Number of face items = " << face_data_num << "\n";

  face_data = new int[face_data_num];
  face_pointer = new int[face_num+1];
  xyz = new double[3*point_num];

  xyzf_example ( point_num, face_num, face_data_num, xyz, face_pointer, 
    face_data );

  xyz_write ( xyz_filename, point_num, xyz );

  xyzf_write ( xyzf_filename, point_num, face_num, face_data_num, 
    face_pointer, face_data );

  cout << "\n";
  cout << "  Wrote the XYZ file \"" << xyz_filename << "\".\n";
  cout << "  and the XYZF file \"" << xyzf_filename << "\".\n";

  delete [] face_data;
  delete [] face_pointer;
  delete [] xyz;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests XYZF_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int face;
  int *face_data;
  int *face_pointer;
  int face_data_num;
  int face_num;
  int i;
  int j;
  int point_num;
  double *xyz;
  string xyz_filename = "cube.xyz";
  string xyzf_filename = "cube.xyzf";

  cout << "\n";
  cout << "TEST06\n";
  cout << "  XYZ_HEADER_READ  reads the header of an XYZ  file.\n";
  cout << "  XYZ_DATA_READ    reads the data   of an XYZ  file.\n";
  cout << "  XYZF_HEADER_READ reads the header of an XYZF file.\n";
  cout << "  XYZF_DATA_READ   reads the data   of an XYZF file.\n";

  cout << "\n";
  cout << "  Examine XYZ file \"" << xyz_filename << "\".\n";

  xyz_header_read ( xyz_filename, &point_num );

  cout << "\n";
  cout << "  Number of points = " << point_num << "\n";

  xyz = new double[3*point_num];

  xyz_data_read ( xyz_filename, point_num, xyz );

  cout << "\n";
  cout << "  Point data:\n";
  cout << "\n";

  for ( i = 0; i < point_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << xyz[0+i*3]
         << "  " << setw(10) << xyz[1+i*3]
         << "  " << setw(10) << xyz[2+i*3] << "\n";
  }

  cout << "\n";
  cout << "  Examine XYZF file \"" << xyzf_filename << "\".\n";

  xyzf_header_read ( xyzf_filename, &face_num, &face_data_num );

  cout << "\n";
  cout << "  Number of faces      = " << face_num << "\n";
  cout << "  Number of face items = " << face_data_num << "\n";

  face_data = new int[face_data_num];
  face_pointer = new int[face_num+1];

  xyzf_data_read ( xyzf_filename, face_num, face_data_num, face_pointer, 
    face_data );

  cout << "\n";
  cout << "  Face pointers:\n";
  cout << "\n";

  for ( face = 0; face < face_num; face++ )
  {
    cout << "  " << setw(4) << face
         << "  " << setw(8) << face_pointer[face]
         << "  " << setw(8) << face_pointer[face+1] - 1 << "\n";
  }

  cout << "\n";
  cout << "  Face data:\n";
  cout << "\n";

  for ( face = 0; face < face_num; face++ )
  {
    cout << "  " << setw(4) << face << "    ";
    for ( j = face_pointer[face]; j <= face_pointer[face+1] - 1; j++ )
    {
      cout << "  " << setw(8) << face_data[j];
    }
    cout << "\n";
  }

  delete [] face_data;
  delete [] face_pointer;
  delete [] xyz;

  return;
}

