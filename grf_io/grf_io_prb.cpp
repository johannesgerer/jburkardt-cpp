# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "grf_io.hpp"

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_IO_PRB calls the GRF_IO test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "GRF_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the GRF_IO library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "GRF_IO_PRB:\n";
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
//    TEST01 tests GRF_HEADER_WRITE, GRF_DATA_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int edge_num;
  int *edge_data;
  int *edge_pointer;
  int node_num;
  string output_filename = "coxeter.grf";
  double *xy;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  GRF_HEADER_WRITE writes the header of a GRF file.\n";
  cout << "  GRF_DATA_WRITE writes the data of a GRF file.\n";

  grf_example_size ( &node_num, &edge_num );

  grf_header_print ( node_num, edge_num );

  edge_data = new int[edge_num];
  edge_pointer= new int[node_num+1];
  xy = new double[2*node_num];

  grf_example ( node_num, edge_num, edge_pointer, edge_data, xy );

  grf_write ( output_filename, node_num, edge_num, edge_pointer,
    edge_data, xy );

  cout << "\n";
  cout << "  Wrote the GRF file \"" << output_filename << "\",\n";

  delete [] edge_data;
  delete [] edge_pointer;
  delete [] xy;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests GRF_HEADER_READ and GRF_DATA_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int edge_num;
  string input_filename = "coxeter.grf";
  int *edge_data;
  int *edge_pointer;
  int node_num;
  double *xy;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  GRF_HEADER_READ reads the header of a GRF file.\n";
  cout << "  GRF_DATA_READ reads the data of a GRF file.\n";

  cout << "\n";
  cout << "  Reading the GRF file \"" << input_filename << "\"\n";

  grf_header_read ( input_filename, &node_num, &edge_num );

  grf_header_print ( node_num, edge_num );

  edge_pointer = new int[node_num+1];
  edge_data = new int[edge_num];
  xy = new double[2*node_num];

  grf_data_read ( input_filename, node_num, edge_num, edge_pointer,
    edge_data, xy );

  grf_data_print ( node_num, edge_num, edge_pointer, edge_data, xy );

  delete [] edge_data;
  delete [] edge_pointer;
  delete [] xy;

  return;
}
