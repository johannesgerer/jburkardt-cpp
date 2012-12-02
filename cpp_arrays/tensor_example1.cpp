# include <iostream>
# include <vector>

using namespace std;

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    Demonstrate one way to create and use 3D arrays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2012
//
//  Author:
//
//    Mauro Perego
//
{
  int m(5);
  int n(2);
  int p(3);

  cout << "\n";
  cout << "TENSOR_EXAMPLE1:\n";
  cout << "  C++ version.\n";
  cout << "  Demonstrate the creation and use of 3D arrays.\n";
//
//  Compact declaration, definition and initialization to 0 for TENSOR1:
//
  vector<vector<vector<double> > > tensor2(m, vector<vector<double> > (n, vector<double>(p, 0.0)));
//
//  Long form declaration for TENSOR1:
//
  vector<vector<vector<double> > > tensor1;
//
//  Definition for TENSOR1:
//
  tensor1.resize ( m );
  for ( int i = 0; i < m; i++ )
  {
    tensor1[i].resize ( n );
    for ( int j = 0; j < n; j++ )
    {
      tensor1[i][j].resize ( p );
      for ( int k = 0; k < p; k++ )
      {
        tensor1[i][j][k] = 0.0;
      }
    }
  }
//
//  Insert values into TESNSOR2.
//
  for ( int i = 0; i < m; i++ )
  {
    for ( int j = 0; j < n; j++ )
    {
      for ( int k = 0; k < p; k++ )
      {
        tensor2[i][j][k] = i+2*j-k;
      }
    }
  }
//
//  Copy Tensor2 in Tensor1.
//
  tensor1 = tensor2;
//
//  Print out some entries.
//
  cout << "\n";
  cout << "T1[3][1][0]: " << tensor1[3][1][0] 
       << ",  T1[4][1][2]: " << tensor1[4][1][2] << endl;
//
//  Terminate.
//
  cout << "\n";
  cout << "TENSOR_EXAMPLE1:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
