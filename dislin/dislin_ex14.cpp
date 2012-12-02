# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "dislin.hpp"

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN demonstrates the use of TeX for mathematical formulas.
//
//  Modified:
// 
//    09 April 2011
//
//  Reference:
//
//    Helmut Michels,
//    The Data Plotting Software DISLIN - version 10.4,
//    Shaker Media GmbH, January 2010,
//    ISBN13: 978-3-86858-517-9.
//
{
  char cstr[80];
  int nl;

  cout << "\n";
  cout << "DISLIN_EX14:\n";
  cout << "  C++ version:\n";
  cout << "  Using TeX for mathematical formulas.\n";
//
//  Specify the format of the output file.
//
  metafl ( "png" );
//
//  Indicate that new data overwrites old data.
//
  filmod ( "delete" );
//
//  Specify the name of the output graphics file.
//
  setfil ( "dislin_ex14.png" );
//
//  Choose the page size and orientation.
//
  setpag ( "usal" );
//
//  For PNG output, reverse the default black background to white.
//
  scrmod ( "reverse" );
//
//  Open DISLIN.
//
  disini ( );
//
//  Plot a border around the page.
//
  pagera ( );
//
//  Use the COMPLEX font.
//
  complx ( );
  height ( 40 );

  strcpy ( cstr, "TeX Instructions for Mathematical Formulas" );
  nl = nlmess ( cstr );
  messag ( cstr, ( 2100 - nl ) / 2, 100 );
  
  texmod ( "on" );
  messag ( "$\\frac{1}{x+y}$", 150, 400 );
  messag ( "$\\frac{a^2 - b^2}{a+b} = a - b$", 1200, 400 );
  
  messag ( "$r = \\red{\\sqrt{x^2 + y^2}}", 150, 700);
  messag ( "$\\cos \\phi = \\frac{x}{\\sqrt{x^2 + y^2}}$", 1200, 700 );

  messag ( "$\\Gamma(x) = \\int_0^\\infty e^{-t}t^{x-1}dt$", 150, 1000 );
  messag ( "$\\lim_{x \\to \\infty} (1 + \\frac{1}{x})^x = e$", 1200, 1000 );

  messag ( "$\\mu = \\sum_{i=1}^n x_i p_i$", 150, 1300 );
  messag ( "$\\mu = \\int_{-\\infty}^ \\infty x f(x) dx$", 1200, 1300 );

  messag ( "$\\overline{x} = \\frac{1}{n} \\sum_{i=1}^n x_i$", 150, 1600 );
  messag ( "$s^2 = \\frac{1}{n-1} \\sum_{i=1}^n (x_i - \\overline{x})^2$",
          1200, 1600 );

  messag ( "$\\sqrt[n]{\\frac{x^n - y^n}{1 + u^{2n}}}$", 150, 1900 );  
  messag ( "$\\sqrt[3]{-q + \\sqrt{q^2 + p^3}}$", 1200, 1900 );

  messag ( "$\\int \\frac{dx}{1+x^2} = \\arctan x + C$", 150, 2200 );
  messag ( "$\\int \\frac{dx}{\\sqrt{1+x^2}} = {\\rm arcsinh} x + C$",
          1200, 2200 );

  messag ( "$\\overline{P_1P_2} = \\sqrt{(x_2-x_1)^2 + (y_2-y_1)^2}$",
           150, 2500 );
  messag ( "$x = \\frac{x_1 + \\lambda x_2}{1 + \\lambda}$", 1200, 2500 );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX14:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}  
