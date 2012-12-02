# include <iostream>
# include <vector>

using namespace std;

typedef typename std::vector<double>::size_type   size_tipe;
class Tensor : public std::vector<double>
{
public:
    //constructor
    Tensor(size_tipe _m, size_tipe _n, size_tipe _p, double initial_value=0.0):
        std::vector<double>(_m*_n*_p, initial_value), p(_p), np(_n*_p) {};

    double& operator()(size_tipe i, size_tipe j, size_tipe k) {return (*this)[i*np + p*j + k];}
    const double& operator()(size_tipe i, size_tipe j, size_tipe k) const {return (*this)[i*np + p*j + k];}

private:
    size_tipe p, np;
};

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
  Tensor tensor1(m,n,p);
  Tensor tensor2(m,n,p,0.0);

  cout << "\n";
  cout << "TENSOR_EXAMPLE2:\n";
  cout << "  C++ version.\n";
  cout << "  Demonstrate the creation and use of 3D arrays.\n";

  for ( int i = 0; i < m; i++ )
  {
    for ( int j = 0; j < n; j++ )
    {
      for ( int k = 0; k < p; k++ )
      {
        tensor2(i,j,k) = i+2*j-k;
      }
    }
  }
//
//  Copy one tensor into another.
//
    tensor1 = tensor2;
//
// Print some entries.
//
  cout << "\n";
  cout << "T[3][1][0]: " << tensor1(3,1,0) 
       << ",  T[4][1][2]: " << tensor1(4,1,2) << endl;
//
//  Terminate.
//
  cout << "\n";
  cout << "TENSOR_EXAMPLE2:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
