# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "stla_io.hpp"

int main ( );

void test01 ( string file_name );
void test02 ( string file_name );
void test03 ( string file_name );
void test04 ( string file_name );
void test05 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for STLA_IO_PRB.
//
//  Discussion:
//
//    STLA_IO_PRB runs the tests of the STLA_IO routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "STLA_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the STLA_IO library.\n";

  test01 ( "cube.stla" );
  test02 ( "cube.stla" );
  test03 ( "cube.stla" );
  test04 ( "cube_new.stla" );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "STLA_IO_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( string input_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests STLA_CHECK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{ 
  cout << "\n";
  cout << "TEST01\n";
  cout << "  STLA_CHECK makes some simple checks on a file.\n";

  cout << "\n";
  if ( stla_check ( input_file_name ) )
  {
    cout << "  The file \"" << input_file_name
      << "\" seems to be a legal ASCII STL file.\n";
  }
  else
  {
    cout << "  The file \"" << input_file_name
      << "\" does NOT seem to be a legal ASCII STL file.\n";
  }

  return;
}
//****************************************************************************80

void test02 ( string input_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests STLA_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int face_num;
  int node_num;
  int solid_num;
  int text_num;
  
  cout << "\n";
  cout << "TEST02\n";
  cout << "  STLA_SIZE determines the size of various objects\n";
  cout << "  in an ASCII STL file.\n";

  stla_size ( input_file_name, &solid_num, &node_num, &face_num, &text_num );

  stla_size_print ( input_file_name, solid_num, node_num, face_num, text_num );

  return;
}
//****************************************************************************80

void test03 ( string input_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests STLA_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  bool error;
  int *face_node;
  double *face_normal;
  int face_num;
  int ierror;
  int node_num;
  double *node_xyz;
  int solid_num;
  int text_num;
  
  cout << "\n";
  cout << "TEST03\n";
  cout << "  STLA_READ reads an object in an ASCII STL file.\n";

  stla_size ( input_file_name, &solid_num, &node_num, &face_num, &text_num );

  face_node = new int[3*face_num];
  face_normal = new double[3*face_num];
  node_xyz = new double[3*node_num];

  error = stla_read ( input_file_name, node_num, face_num, node_xyz, 
    face_node, face_normal );

  if ( error )
  {
    cout << "\n";
    cout << "  STLA_READ returned error flag.\n";
    return;
  }

  stla_size_print ( input_file_name, solid_num, node_num, face_num, text_num );

  stla_face_node_print ( face_num, face_node );
  stla_face_normal_print ( face_num, face_normal );
  stla_node_xyz_print ( node_num, node_xyz );

  delete [] face_node;
  delete [] face_normal;
  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test04 ( string output_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests STLA_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define FACE_NUM 12
# define NODE_NUM 8

  int face_node[3*FACE_NUM] = {
    1, 3, 2, 
    2, 3, 4, 
    1, 6, 5,
    1, 2, 6, 
    3, 7, 4, 
    4, 7, 8, 
    5, 6, 8, 
    5, 8, 7, 
    1, 5, 7, 
    1, 7, 3, 
    2, 4, 6, 
    6, 4, 8 };
  double face_normal[3*FACE_NUM] = {
    0.0,  0.0, -1.0, 
    0.0,  0.0, -1.0, 
    0.0, -1.0,  0.0, 
    0.0, -1.0,  0.0, 
    0.0, +1.0,  0.0, 
    0.0, +1.0,  0.0, 
    0.0,  0.0, +1.0, 
    0.0,  0.0, +1.0, 
   -1.0,  0.0,  0.0, 
   -1.0,  0.0,  0.0, 
   +1.0,  0.0,  0.0, 
   +1.0,  0.0,  0.0 }; 
  double node_xyz[3*NODE_NUM] = {
    0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    1.0, 1.0, 0.0, 
    0.0, 0.0, 1.0, 
    1.0, 0.0, 1.0, 
    0.0, 1.0, 1.0, 
    1.0, 1.0, 1.0 };
  int offset;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  STLA_WRITE writes an ASCII STL file.\n";

  offset = 1;
  stla_offset_set ( offset );

  stla_write ( output_file_name, NODE_NUM, FACE_NUM, node_xyz, 
    face_node, face_normal );

  cout << "\n";
  cout << "  Graphics data was written to the STLA file \""
    << output_file_name << "\".\n";

  return;
# undef FACE_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test05 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests STLA_FACE_NORMAL_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define FACE_NUM 12
# define NODE_NUM 8

  double dot_max;
  int face;
  int face_node[3*FACE_NUM] = {
    1, 3, 2, 
    2, 3, 4, 
    1, 6, 5, 
    1, 2, 6, 
    3, 7, 4, 
    4, 7, 8, 
    5, 6, 8, 
    5, 8, 7, 
    1, 5, 7, 
    1, 7, 3, 
    2, 4, 6, 
    6, 4, 8 };
  double face_normal[3*FACE_NUM] = {
    0.0,  0.0, -1.0, 
    0.0,  0.0, -1.0, 
    0.0, -1.0,  0.0, 
    0.0, -1.0,  0.0, 
    0.0, +1.0,  0.0, 
    0.0, +1.0,  0.0, 
    0.0,  0.0, +1.0, 
    0.0,  0.0, +1.0, 
   -1.0,  0.0,  0.0, 
   -1.0,  0.0,  0.0, 
   +1.0,  0.0,  0.0, 
   +1.0,  0.0,  0.0 };
  double *face_normal2;
  double node_xyz[3*NODE_NUM] = {
    0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    1.0, 1.0, 0.0, 
    0.0, 0.0, 1.0, 
    1.0, 0.0, 1.0, 
    0.0, 1.0, 1.0, 
    1.0, 1.0, 1.0 };
  int offset;

  offset = 1;
  stla_offset_set ( offset );

  cout << "\n";
  cout << "TEST05\n";
  cout << "  STLA_FACE_NORMAL_COMPUTE computes the face normal\n";
  cout << "  vectors for an STLA file.\n";
  cout << "\n";
  cout << "  We have an STLA solid, and its exact normals.\n";
  cout << "  We now call STLA_FACE_NORMAL_COMPUTE to\n";
  cout << "  recompute the normals.\n";

  face_normal2 = stla_face_normal_compute ( NODE_NUM, FACE_NUM, node_xyz, 
    face_node );

  cout << "\n";
  cout << "  We print out the maximum error, defined as\n";
  cout << "    |1 - dot ( n1, n2 )|\n";
  cout << "  where n1 and n2 are the exact and computed normals.\n";

  dot_max = 0.0;

  for ( face = 0; face < FACE_NUM; face++ )
  {
    dot_max = r8_max ( dot_max, fabs ( 1.0 - 
      r8vec_dot ( 3, face_normal+face*3, face_normal2+face*3 ) ) );
  }

  cout << "\n";
  cout << "  Maximum error = " << dot_max << "\n";

  delete [] face_normal2;

  return;
# undef FACE_NUM
# undef NODE_NUM
}
