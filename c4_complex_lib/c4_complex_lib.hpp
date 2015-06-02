# ifndef __c4_complex_h__
# define __c4_complex_h__


class c4_complex

//
//  Modified:
//
//    27 September 2013
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

    float real_part;
    float imaginary_part;

  public:

    c4_complex ( void )
    {
      real_part = 0.0;
      imaginary_part = 0.0;
    }

    c4_complex ( const c4_complex &other_c4_complex )
    {
      real_part = other_c4_complex.real_part;
      imaginary_part = other_c4_complex.imaginary_part;
    }

    c4_complex ( float init_real, float init_imaginary = 0.0 )
    {
      real_part = init_real;
      imaginary_part = init_imaginary;
    }
    ~c4_complex ( )
    {
    }
    float real ( void ) const
    {
      return ( real_part );
    }
    float imaginary ( void ) const
    {
      return ( imaginary_part );
    }
    void set ( float real, float imaginary )
    {
      real_part = real;
      imaginary_part = imaginary;
    }
    void set_real ( float real )
    {
      real_part = real;
    }
    void set_imaginary ( float imaginary )
    {
      imaginary_part = imaginary;
    }
    c4_complex operator = ( const c4_complex &oper2 )
    {
      set ( oper2.real_part, oper2.imaginary_part );
      return ( *this );
    }
    c4_complex &operator += ( const c4_complex &oper2 )
    {
      real_part = real_part + oper2.real ( );
      imaginary_part = imaginary_part + oper2.imaginary ( );
      return ( *this );
    }
    c4_complex &operator += ( float oper2 )
    {
      real_part = real_part + oper2;
      return ( *this );
    }
    c4_complex &operator -= ( const c4_complex &oper2 )
    {
      real_part = real_part - oper2.real ( );
      imaginary_part = imaginary_part - oper2.imaginary ( );
      return ( *this );
    }
    c4_complex &operator -= ( float &oper2 )
    {
      real_part = real_part - oper2;
      return ( *this );
    }
    c4_complex &operator *= ( const c4_complex &oper2 )
    {
      float real_result = real_part * oper2.real ( )
        - imaginary_part * oper2.imaginary ( );
      imaginary_part = real_part * oper2.imaginary ( )
        + imaginary_part * oper2.real ( );
      real_part = real_result;
      return ( *this );
    }
    c4_complex &operator *= ( float oper2 )
    {
      real_part = real_part * oper2;
      imaginary_part = imaginary_part * oper2;
      return ( *this );
    }
    c4_complex &operator /= ( const c4_complex &oper2 );

    c4_complex &operator /= ( float oper2 )
    {
      real_part /= oper2;
      imaginary_part /= oper2;
      return ( *this );
    }

    c4_complex operator ++ ( int )
    {
      c4_complex result ( *this );
      real_part = real_part + 1.0;
      return ( result );
    }
    c4_complex &operator ++ ( void )
    {
      real_part = real_part + 1.0;
      return ( *this );
    }
    c4_complex operator -- ( int )
    {
      c4_complex result ( *this );
      real_part = real_part - 1.0;
      return ( result );
    }
    c4_complex operator -- ( void )
    {
      real_part = real_part - 1.0;
      return ( *this );
    }
    c4_complex &operator ~ ( void )
    {
      imaginary_part = -imaginary_part;
      return ( *this );
    }
};

inline c4_complex operator + ( const c4_complex &oper1, const c4_complex &oper2 )
{
  return c4_complex ( oper1.real ( ) + oper2.real ( ),
    oper1.imaginary ( ) + oper2.imaginary ( ) );
}

inline c4_complex operator + ( const c4_complex &oper1, float oper2 )
{
  return c4_complex ( oper1.real ( ) + oper2,
    oper1.imaginary ( ) );
}

inline c4_complex operator + ( float oper1, const c4_complex &oper2 )
{
  return c4_complex ( oper1 + oper2.real ( ),
    oper2.imaginary ( ) );
}

inline c4_complex operator - ( const c4_complex &oper1, const c4_complex &oper2 )
{
  return c4_complex ( oper1.real ( ) - oper2.real ( ),
    oper1.imaginary ( ) - oper2.imaginary( ) );
}
inline c4_complex operator - ( const c4_complex &oper1, float oper2 )
{
  return c4_complex ( oper1.real ( ) - oper2,
    oper1.imaginary ( ) );
}

inline c4_complex operator - ( float oper1, const c4_complex &oper2 )
{
  return c4_complex ( oper1 - oper2.real ( ),
    oper2.imaginary ( ) );
}


inline c4_complex operator * ( const c4_complex &oper1, const c4_complex &oper2 )
{
  return c4_complex ( oper1.real ( ) * oper2.real ( ) 
                - oper1.imaginary ( ) * oper2.imaginary ( ),
                 oper1.real ( ) * oper2.imaginary ( ) 
                + oper1.imaginary ( ) * oper2.real ( ) );
}
inline c4_complex operator * ( const c4_complex &oper1, float oper2 )
{
  return c4_complex ( oper1.real ( ) * oper2,
    oper1.imaginary ( ) * oper2 );
}

inline c4_complex operator * ( float oper1, const c4_complex &oper2 )
{
  return c4_complex ( oper1 * oper2.real ( ),
    oper1 * oper2.imaginary ( ) );
}

extern c4_complex operator / ( const c4_complex &oper1, const c4_complex &oper2 );

inline c4_complex operator / ( const float &oper1, const c4_complex &oper2 )
{
  return ( c4_complex ( oper1, 0.0 ) / oper2 );
}
inline c4_complex operator / ( const c4_complex &oper1, const float &oper2 )
{
  return ( oper1 / c4_complex ( oper2, 0.0 ) );
}

inline int operator == ( const c4_complex &oper1, const c4_complex &oper2 )
{
  return ( ( oper1.real ( ) == oper2.real ( ) ) &&
           ( oper1.imaginary ( ) == oper2.imaginary ( ) ) );
}

inline int operator != ( const c4_complex &oper1, const c4_complex oper2 )
{
  return ( ! ( oper1 == oper2 ) );
}

inline c4_complex operator - ( const c4_complex &oper1 )
{
  return ( c4_complex ( - oper1.real ( ), - oper1.imaginary ( ) ) );
}

inline c4_complex operator + ( const c4_complex &oper1 )
{
  return ( c4_complex ( +oper1.real ( ), +oper1.imaginary ( ) ) );
}

inline ostream &operator << ( ostream &out_file, const c4_complex &number )
{
  out_file << '(' << number.real ( ) << ',' << number.imaginary ( ) << ')';
}

extern istream &operator >> ( istream &in_file, c4_complex &number );

# endif

