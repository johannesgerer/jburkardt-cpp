# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "xy_io.hpp"

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
//    XY_IO_PRB calls the XY_IO test routines.
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
  cout << "XY_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the XY_IO library.\n";

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
  cout << "XY_IO_PRB:\n";
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
//    TEST01 tests XY_EXAMPLE, XY_WRITE.
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
# define POINT_NUM 300

  string file_name = "xy_io_prb_01.xy";
  int point_num = POINT_NUM;
  double xy[2*POINT_NUM];

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  XY_EXAMPLE sets up sample XY data.\n";
  cout << "  XY_WRITE writes an XY file.\n";

  xy_example ( point_num, xy );

  cout << "\n";
  cout << "  XY_EXAMPLE has created the data.\n";

  xy_write ( file_name, point_num, xy );

  cout << "\n";
  cout << "  XY_WRITE wrote the header and data for \"" << file_name << "\"\n";
  cout << "  Number of points = " << point_num << "\n";
  cout << flush;

  return;
# undef POINT_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests XY_READ.
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
  string file_name = "xy_io_prb_02.xy";
  int i;
  int k;
  int point_num;
  double *xy;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  XY_READ reads the header and data of an XY file.\n";

  xy_write_test ( file_name );

  cout << "\n";
  cout << "  XY_WRITE_TEST created data and wrote it to \"" 
       << file_name << "\"\n";
//
//  Now have XY_READ try to read it.
//
  xy_read ( file_name, &point_num, &xy );

  cout << "\n";
  cout << "  XY_READ read the test data successfully.\n";
  cout << "  Number of points = " << point_num << "\n";
  cout << "\n";
  cout << "  Sample data:\n";
  cout << "\n";
  for ( k = 0; k <= 9; k++ )
  {
    i = ( ( 9 - k ) * 0 + k * ( point_num - 1 ) ) / 9;
    cout << setw(4)  << i    << "  "
         << setw(10) << xy[0+i*2] << "  "
         << setw(10) << xy[1+i*2] << "\n";
  }

  delete [] xy;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests XYL_EXAMPLE.
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
  double *xy;
  string xy_filename = "house.xy";
  string xyl_filename = "house.xyl";

  cout << "\n";
  cout << "TEST03\n";
  cout << "  XYL_EXAMPLE sets up XY and XYL data.\n";

  xyl_example_size ( &point_num, &line_num, &line_data_num );

  cout << "\n";
  cout << "  Example has:\n";
  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";
  cout << "  Number of lines      = " << line_num << "\n";
  cout << "  Number of line items = " << line_data_num << "\n";

  line_data = new int[line_data_num];
  line_pointer= new int[line_num+1];
  xy = new double[2*point_num];

  xyl_example ( point_num, line_num, line_data_num, xy, line_pointer, line_data );

  xy_write ( xy_filename, point_num, xy );

  xyl_write ( xyl_filename, point_num, line_num, line_data_num, 
    line_pointer, line_data );

  cout << "\n";
  cout << "  Wrote the XY file \"" << xy_filename << "\",\n";
  cout << "  and the XYL file \"" << xyl_filename << "\".\n";

  delete [] line_data;
  delete [] line_pointer;
  delete [] xy;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests XYL_READ.
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
  double *xy;
  string xy_filename = "house.xy";
  string xyl_filename = "house.xyl";

  cout << "\n";
  cout << "TEST04\n";
  cout << "  XY_HEADER_READ  reads the header of an XY  file.\n";
  cout << "  XY_DATA_READ    reads the data   of an XY  file.\n";
  cout << "  XYL_HEADER_READ reads the header of an XYL file.\n";
  cout << "  XYL_DATA_READ   reads the data   of an XYL file.\n";

  cout << "\n";
  cout << "  Examine XY file \"" << xy_filename << "\".\n";

  xy_header_read ( xy_filename, &point_num );

  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";

  xy = new double[2*point_num];

  xy_data_read ( xy_filename, point_num, xy );

  cout << "\n";
  cout << "  Point data:\n";
  cout << "\n";

  for ( i = 0; i < point_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << xy[0+2*i]
         << "  " << setw(10) << xy[1+2*i] << "\n";
  }

  cout << "\n";
  cout << "  Examine XYL file \"" << xyl_filename << "\".\n";

  xyl_header_read ( xyl_filename, &line_num, &line_data_num );

  cout << "\n";
  cout << "  Number of lines      = " << line_num << "\n";
  cout << "  Number of line items = " << line_data_num << "\n";

  line_data = new int[line_data_num];
  line_pointer = new int[line_num+1];

  xyl_data_read ( xyl_filename, line_num, line_data_num, line_pointer, line_data );

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
    cout << "  " << setw(4) << line << "  ";
    for ( j = line_pointer[line]; j < line_pointer[line+1]; j++ )
    {
      cout << "  " << setw(8) << line_data[j];
    }
    cout << "\n";
  }

  delete [] line_data;
  delete [] line_pointer;
  delete [] xy;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests XYF_EXAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2009
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
  double *xy;
  string xy_filename = "annulus.xy";
  string xyf_filename = "annulus.xyf";

  cout << "\n";
  cout << "TEST05\n";
  cout << "  XYF_EXAMPLE sets up XY and XYF data.\n";

  xyf_example_size ( &point_num, &face_num, &face_data_num );

  cout << "\n";
  cout << "  Example has:\n";
  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";
  cout << "  Number of faces      = " << face_num << "\n";
  cout << "  Number of face items = " << face_data_num << "\n";

  face_data = new int[face_data_num];
  face_pointer= new int[face_num+1];
  xy = new double[2*point_num];

  xyf_example ( point_num, face_num, face_data_num, xy, face_pointer, face_data );

  xy_write ( xy_filename, point_num, xy );

  xyf_write ( xyf_filename, point_num, face_num, face_data_num, 
    face_pointer, face_data );

  cout << "\n";
  cout << "  Wrote the XY file \"" << xy_filename << "\",\n";
  cout << "  and the XYF file \"" << xyf_filename << "\".\n";

  delete [] face_data;
  delete [] face_pointer;
  delete [] xy;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests XYF_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2009
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
  double *xy;
  string xy_filename = "annulus.xy";
  string xyf_filename = "annulus.xyf";

  cout << "\n";
  cout << "TEST06\n";
  cout << "  XY_HEADER_READ  reads the header of an XY  file.\n";
  cout << "  XY_DATA_READ    reads the data   of an XY  file.\n";
  cout << "  XYF_HEADER_READ reads the header of an XYF file.\n";
  cout << "  XYF_DATA_READ   reads the data   of an XYF file.\n";

  cout << "\n";
  cout << "  Examine XY file \"" << xy_filename << "\".\n";

  xy_header_read ( xy_filename, &point_num );

  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";

  xy = new double[2*point_num];

  xy_data_read ( xy_filename, point_num, xy );

  cout << "\n";
  cout << "  Point data:\n";
  cout << "\n";

  for ( i = 0; i < point_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << xy[0+2*i]
         << "  " << setw(10) << xy[1+2*i] << "\n";
  }

  cout << "\n";
  cout << "  Examine XYF file \"" << xyf_filename << "\".\n";

  xyf_header_read ( xyf_filename, &face_num, &face_data_num );

  cout << "\n";
  cout << "  Number of faces      = " << face_num << "\n";
  cout << "  Number of face items = " << face_data_num << "\n";

  face_data = new int[face_data_num];
  face_pointer = new int[face_num+1];

  xyf_data_read ( xyf_filename, face_num, face_data_num, face_pointer, face_data );

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
    cout << "  " << setw(4) << face << "  ";
    for ( j = face_pointer[face]; j < face_pointer[face+1]; j++ )
    {
      cout << "  " << setw(8) << face_data[j];
    }
    cout << "\n";
  }

  delete [] face_data;
  delete [] face_pointer;
  delete [] xy;

  return;
}
