# include <boost/lambda/lambda.hpp>
# include <iostream>
# include <iterator>
# include <algorithm>

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    BOOST_EXAMPLE1 reads integers, multiplies by 3, and writes them out.
//
//  Discussion:
//
//    The program reads integers from the standard input, and writes out
//    3 times the value of each quantity input.  The program doesn't terminate,
//    and continues to read and triple until killed.
//
//  Modified:
//
//    20 April 2011
//
{
  using namespace boost::lambda;
  typedef std::istream_iterator<int> in;

  std::for_each
  (
    in ( std::cin ), 
    in ( ), 
    std::cout << (_1 * 3) << " "
  );

  return 0;
}
