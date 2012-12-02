# ifndef __c8_complex_h__
# define __c8_complex_h__


class c8_complex

//
//  Modified:
//
//    20 October 2005
//
//  Reference:
//
//    Steve Oualline,
//    Practical C++ Programming,
//    O'Reilly & Associates, 1997,
//    ISBN: 1-56592-139-9
//
{
  private:

    double real_part;
    double imaginary_part;

  public:

    c8_complex ( void )
    {
      real_part = 0.0;
      imaginary_part = 0.0;
    }

    c8_complex ( const c8_complex &other_c8_complex )
    {
      real_part = other_c8_complex.real_part;
      imaginary_part = other_c8_complex.imaginary_part;
    }

    c8_complex ( double init_real, double init_imaginary = 0.0 )
    {
      real_part = init_real;
      imaginary_part = init_imaginary;
    }
    ~c8_complex ( )
    {
    }
    double real ( void ) const
    {
      return ( real_part );
    }
    double imaginary ( void ) const
    {
      return ( imaginary_part );
    }
    void set ( double real, double imaginary )
    {
      real_part = real;
      imaginary_part = imaginary;
    }
    void set_real ( double real )
    {
      real_part = real;
    }
    void set_imaginary ( double imaginary )
    {
      imaginary_part = imaginary;
    }
    c8_complex operator = ( const c8_complex &oper2 )
    {
      set ( oper2.real_part, oper2.imaginary_part );
      return ( *this );
    }
    c8_complex &operator += ( const c8_complex &oper2 )
    {
      real_part = real_part + oper2.real ( );
      imaginary_part = imaginary_part + oper2.imaginary ( );
      return ( *this );
    }
    c8_complex &operator += ( double oper2 )
    {
      real_part = real_part + oper2;
      return ( *this );
    }
    c8_complex &operator -= ( const c8_complex &oper2 )
    {
      real_part = real_part - oper2.real ( );
      imaginary_part = imaginary_part - oper2.imaginary ( );
      return ( *this );
    }
    c8_complex &operator -= ( double &oper2 )
    {
      real_part = real_part - oper2;
      return ( *this );
    }
    c8_complex &operator *= ( const c8_complex &oper2 )
    {
      double real_result = real_part * oper2.real ( )
        - imaginary_part * oper2.imaginary ( );
      imaginary_part = real_part * oper2.imaginary ( )
        + imaginary_part * oper2.real ( );
      real_part = real_result;
      return ( *this );
    }
    c8_complex &operator *= ( double oper2 )
    {
      real_part = real_part * oper2;
      imaginary_part = imaginary_part * oper2;
      return ( *this );
    }
    c8_complex &operator /= ( const c8_complex &oper2 );

    c8_complex &operator /= ( double oper2 )
    {
      real_part /= oper2;
      imaginary_part /= oper2;
      return ( *this );
    }

    c8_complex operator ++ ( int )
    {
      c8_complex result ( *this );
      real_part = real_part + 1.0;
      return ( result );
    }
    c8_complex &operator ++ ( void )
    {
      real_part = real_part + 1.0;
      return ( *this );
    }
    c8_complex operator -- ( int )
    {
      c8_complex result ( *this );
      real_part = real_part - 1.0;
      return ( result );
    }
    c8_complex operator -- ( void )
    {
      real_part = real_part - 1.0;
      return ( *this );
    }
    c8_complex &operator ~ ( void )
    {
      imaginary_part = -imaginary_part;
      return ( *this );
    }
};

inline c8_complex operator + ( const c8_complex &oper1, const c8_complex &oper2 )
{
  return c8_complex ( oper1.real ( ) + oper2.real ( ),
    oper1.imaginary ( ) + oper2.imaginary ( ) );
}

inline c8_complex operator + ( const c8_complex &oper1, double oper2 )
{
  return c8_complex ( oper1.real ( ) + oper2,
    oper1.imaginary ( ) );
}

inline c8_complex operator + ( double oper1, const c8_complex &oper2 )
{
  return c8_complex ( oper1 + oper2.real ( ),
    oper2.imaginary ( ) );
}

inline c8_complex operator - ( const c8_complex &oper1, const c8_complex &oper2 )
{
  return c8_complex ( oper1.real ( ) - oper2.real ( ),
    oper1.imaginary ( ) - oper2.imaginary( ) );
}
inline c8_complex operator - ( const c8_complex &oper1, double oper2 )
{
  return c8_complex ( oper1.real ( ) - oper2,
    oper1.imaginary ( ) );
}

inline c8_complex operator - ( double oper1, const c8_complex &oper2 )
{
  return c8_complex ( oper1 - oper2.real ( ),
    oper2.imaginary ( ) );
}


inline c8_complex operator * ( const c8_complex &oper1, const c8_complex &oper2 )
{
  return c8_complex ( oper1.real ( ) * oper2.real ( ) 
                - oper1.imaginary ( ) * oper2.imaginary ( ),
                 oper1.real ( ) * oper2.imaginary ( ) 
                + oper1.imaginary ( ) * oper2.real ( ) );
}
inline c8_complex operator * ( const c8_complex &oper1, double oper2 )
{
  return c8_complex ( oper1.real ( ) * oper2,
    oper1.imaginary ( ) * oper2 );
}

inline c8_complex operator * ( double oper1, const c8_complex &oper2 )
{
  return c8_complex ( oper1 * oper2.real ( ),
    oper1 * oper2.imaginary ( ) );
}

extern c8_complex operator / ( const c8_complex &oper1, const c8_complex &oper2 );

inline c8_complex operator / ( const double &oper1, const c8_complex &oper2 )
{
  return ( c8_complex ( oper1, 0.0 ) / oper2 );
}
inline c8_complex operator / ( const c8_complex &oper1, const double &oper2 )
{
  return ( oper1 / c8_complex ( oper2, 0.0 ) );
}

inline int operator == ( const c8_complex &oper1, const c8_complex &oper2 )
{
  return ( ( oper1.real ( ) == oper2.real ( ) ) &&
           ( oper1.imaginary ( ) == oper2.imaginary ( ) ) );
}

inline int operator != ( const c8_complex &oper1, const c8_complex oper2 )
{
  return ( ! ( oper1 == oper2 ) );
}

inline c8_complex operator - ( const c8_complex &oper1 )
{
  return ( c8_complex ( - oper1.real ( ), - oper1.imaginary ( ) ) );
}

inline c8_complex operator + ( const c8_complex &oper1 )
{
  return ( c8_complex ( +oper1.real ( ), +oper1.imaginary ( ) ) );
}

inline ostream &operator << ( ostream &out_file, const c8_complex &number )
{
  out_file << '(' << number.real ( ) << ',' << number.imaginary ( ) << ')';
}

extern istream &operator >> ( istream &in_file, c8_complex &number );

# endif

