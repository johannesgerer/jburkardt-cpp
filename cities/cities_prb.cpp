# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "cities.hpp"

int main ( );
void test01 ( string prefix );
void test02 ( string prefix );
void test03 ( string prefix );
void test04 ( string prefix );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CITIES_PRB.
//
//  Discussion:
//
//    CITIES_PRB tests routines from the CITIES library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  string filename;

  timestamp ( );

  cout << "\n";
  cout << "CITIES_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CITIES library.\n";

  filename = "wg22";
  test01 ( filename );

  filename = "usca312";
  test02 ( filename );

  test03 ( filename );

  filename = "spaeth2_09";
  test04 ( filename );
//
//  Terminate.
//
  cout << "\n";
  cout << "CITIES_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests POINT_TO_DIST_TABLE.
//
//  Discussion:
//
//    Get the XY coordinates of a set of cities, and compute
//    the city-to-city distance table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, the common file prefix.
//
{
  int dim_num;
  double *dist_table;
  string dist_table_filename;
  string main_filename;
  double *point;
  string point_filename;
  int point_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  POINT_TO_DIST_TABLE computes a distance table from a\n";
  cout << "  list of point locations.\n";

  main_filename = prefix + "_main.txt";
  point_filename = prefix + "_xy.txt";
  dist_table_filename = prefix + "_dist_table.txt";

  cout << "\n";
  cout << "  The main filename is \"" << main_filename << "\".\n";
  cout << "  The point filename is \"" << point_filename << "\".\n";
  cout << "  The distance table filename is \"" << dist_table_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  The spatial dimension is " << dim_num << "\n";
  cout << "  The number of points is  " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );

  r8mat_transpose_print ( dim_num, point_num, point, "  The points:" );

  dist_table = point_to_dist_table ( dim_num, point_num, point );

  r8mat_nint ( point_num, point_num, dist_table );

  r8mat_print_some ( point_num, point_num, dist_table, 1, 1, 5, 5, 
    "  Initial 5x5 distance subtable:" );

  r8mat_write ( dist_table_filename, point_num, point_num, dist_table );

  delete [] dist_table;
  delete [] point;

  return;
}
//****************************************************************************80

void test02 ( string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests MAIN_READ_SIZE, MAIN_READ_DMS, MAIN_READ_NAME.
//
//  Discussion:
//
//    Get the DMS coordinates of a set of cities, and compute
//    the city-to-city distance table, using distances on a sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, the common file prefix.
//
{
  double *dist;
  string dist_filename;
  string dms_filename;
  int *lat_dms;
  int *long_dms;
  string main_filename;
  int n;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Get the DMS coordinates of a set of cities.\n";
  cout << "  Compute the city-to-city distance table,\n";
  cout << "  assuming the cities lie on a sphere (the earth).\n";

  main_filename = prefix + "_main.txt";
  dms_filename = prefix + "_dms.txt";
  dist_filename = prefix + "_dist.txt";

  cout << "\n";
  cout << "  The main filename is \"" << main_filename << "\".\n";
  cout << "  The dms filename is \"" << dms_filename << "\".\n";
  cout << "  The distance filename will be \"" << dist_filename << "\".\n";

  n = main_read_size ( main_filename );

  cout << "\n";
  cout << "  The number of data items is " << n << "\n";

  dist = new double[n*n];
  lat_dms = new int[4*n];
  long_dms = new int[4*n];

  dms_read ( dms_filename, n, lat_dms, long_dms );

  dms_print ( n, lat_dms, long_dms, "  The longitude/latitude data:" );

  dist = dms_to_dist ( n, lat_dms, long_dms );

  r8mat_nint ( n, n, dist );

  cout << "\n";
  cout << "  Distance from Atlanta to Boston  = " << dist[11+33*n] << "\n";
  cout << "  Road distance is 1037\n";
  cout << "  Distance from Boston to Chicago  = " << dist[33+57*n] << "\n";
  cout << "  Road distance is  963\n";
  cout << "  Distance from Chicago to Atlanta = " << dist[57+11*n] << "\n";
  cout << "  Road distance is  674\n";

  r8mat_write ( dist_filename, n, n, dist );

  delete [] dist;
  delete [] lat_dms;
  delete [] long_dms;

  return;
}
//****************************************************************************80

void test03 ( string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests DMS_TO_XY.
//
//  Discussion:
//
//    Get the DMS coordinates of a set of cities, and compute
//    the XY coordinates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2010;
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, the common file prefix.
//
{
  string dms_filename;
  int *lat_dms;
  int *long_dms;
  string main_filename;
  int n;
  double *point_xy;
  string point_xy_filename;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  DMS_TO_XY takes latitude and longitude\n";
  cout << "  information, and assigns pseudo XY coordinates.\n";

  main_filename = prefix + "_main.txt";
  dms_filename = prefix + "_dms.txt";
  point_xy_filename = prefix + "_xy.txt";

  cout << "\n";
  cout << "  The main filename is \"" << main_filename << "\".\n";
  cout << "  The dms filename is \"" << dms_filename << "\".\n";
  cout << "  The point XY filename will be \"" << point_xy_filename << "\".\n";

  n = main_read_size ( main_filename );

  cout << "\n";
  cout << "  The number of data items is " << n << "\n";

  lat_dms = new int[4*n];
  long_dms = new int[4*n];

  dms_read ( dms_filename, n, lat_dms, long_dms );

  dms_print ( n, lat_dms, long_dms, "  The longitude/latidude data:" );

  point_xy = dms_to_xy ( n, lat_dms, long_dms );

  r8mat_transpose_print ( 2, n, point_xy, "  The computed point values:" );

  r8mat_write ( point_xy_filename, 2, n, point_xy );

  delete [] lat_dms;
  delete [] long_dms;
  delete [] point_xy;

  return;
}


//****************************************************************************80

void test04 ( string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests DIST_TABLE_CHECK.
//
//  Discussion:
//
//    Read a distance matrix and check it for consistency.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, the common file prefix.
//
{
  int check;
  double *dist_table;
  string dist_table_filename;
  int dist_num;
  int n;
  int n1;
  int n2;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  DIST_TABLE_CHECK checks a distance table.\n";

  dist_table_filename = prefix + "_dist_table.txt";

  cout << "\n";
  cout << "  The distance table filename is \"" 
       << dist_table_filename << "\"\n";

  r8mat_header_read ( dist_table_filename, &n1, &n2 );

  if ( n1 != n2 )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  The distance table is not square.\n";
    return;
  }

  n = n1;
  cout << "\n";
  cout << "  The number of data items is " << n << "\n";

  dist_table = r8mat_data_read ( dist_table_filename, n, n );

  check = dist_table_check ( n, dist_table );

  cout << "\n";
  if ( check == 0 )
  {
    cout << "  0: The distance table passed all checks.\n";
  }
  else if ( check == 1 )
  {
    cout << "  1: The table failed the nonnegativity check.\n";
  }
  else if ( check == 2 )
  {
    cout << "  2: The table failed the zero self-distance check.\n";
  }
  else if ( check == 3 )
  {
    cout << "  3: The table failed the symmetry check.\n";
  }
  else if ( check == 4 )
  {
    cout << "  4: The table failed the triangle check.\n";
  }

  delete [] dist_table;

  return;
}
