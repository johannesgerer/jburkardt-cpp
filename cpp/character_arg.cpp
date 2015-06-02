# include <cstdlib>
# include <iostream>
# include <cstring>

using namespace std;

int main ( );
void fred ( string *name );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CHARACTER_ARG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  string name;

  cout << "\n";
  cout << "CHARACTER_ARG:\n";
  cout << "  Demonstrate how a C function can return character data\n";
  cout << "  through the argument list.\n";
  cout << "\n";
  cout << "  Our main program declares a string:\n";
  cout << "    string name;\n";
  cout << "  then calls function fred() with the ADDRESS of the string:\n";
  cout << "    fred ( &name );\n";
  cout << "  Function fred receives its argument as\n";
  cout << "    void fred ( string *name )\n";
  cout << "  It sets a value to the string:\n";
  cout << "    *name = \"ob_data.txt\" );\n";
  cout << "  The main program now has a string stored in name.\n";

  fred ( &name );

  cout << "\n";
  cout << "  The value of name is now = \"" << name << "\".\n";

  cout << "\n";
  cout << "CHARACTER_ARG:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void fred ( string *name )

//****************************************************************************80
//
//  Purpose:
//
//    FRED returns character data through one of its arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameter:
//
//    Input, string *NAME, the address of the address of a character.
//    The value *NAME is the address of a character, and can be set
//    to the address of a string of interest to the user.
//
{
  *name = "ob_data.txt";

  return;
}

