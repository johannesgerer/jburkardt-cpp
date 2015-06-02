#include <boost/serialization/vector.hpp>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef isnan
# define isnan(x) \
	(sizeof (x) == sizeof (long double) ? isnan_ld (x) \
	 : sizeof (x) == sizeof (double) ? isnan_d (x) \
	 : isnan_f (x))
static inline int isnan_f  (float       x) { return x != x; }
static inline int isnan_d  (double      x) { return x != x; }
static inline int isnan_ld (long double x) { return x != x; }
#endif


class pnt {/*{{{*/
	private:
		friend class boost::serialization::access;
		template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
				ar & x;
				ar & y;
				ar & z;
				ar & projTo;
				ar & isBdry;
				ar & idx;
			}
	public:
		double x, y, z;
		int projTo;
		int isBdry;
		int idx;


		pnt(double x_, double y_, double z_, int isBdry_, int idx_, int projTo_)
			:  x(x_), y(y_), z(z_), isBdry(isBdry_), idx(idx_), projTo(projTo_) {	}

		pnt(double x_, double y_, double z_, int isBdry_, int idx_)
			:  x(x_), y(y_), z(z_), isBdry(isBdry_), idx(idx_), projTo(-1) {	}

		pnt(double x_, double y_, double z_)
			: x(x_), y(y_), z(z_), isBdry(0), idx(0), projTo(-1) { }

		pnt(double x_, double y_, double z_, int isBdry_)
			: x(x_), y(y_), z(z_), isBdry(isBdry_), idx(0), projTo(-1) { }

		pnt()
			: x(0.0), y(0.0), z(0.0), isBdry(0), idx(0), projTo(-1) { }

		friend pnt operator*(const double d, const pnt &p);
		friend std::ostream & operator<<(std::ostream &os, const pnt &p);
		friend std::istream & operator>>(std::istream &is, pnt &p);

		pnt& operator=(const pnt &p){/*{{{*/
			x = p.x;
			y = p.y;
			z = p.z;
			isBdry = p.isBdry;
			idx = p.idx;
			projTo = -1;
			return *this;
		}/*}}}*/
		bool operator==(const pnt &p) const {/*{{{*/
			return (x == p.x) & (y == p.y) & (z == p.z);
		}/*}}}*/
		pnt operator-(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x-p.x;
			y_ = y-p.y;
			z_ = z-p.z;

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		pnt operator+(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x+p.x;
			y_ = y+p.y;
			z_ = z+p.z;

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		pnt operator*(double d) const {/*{{{*/
			double x_, y_, z_;

			x_ = x*d;
			y_ = y*d;
			z_ = z*d;

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		pnt operator/(double d) const {/*{{{*/
			double x_, y_, z_;

			if(d == 0.0){
				std::cout << "pnt: operator/" << std::endl;
				std::cout << (*this) << std::endl;
			}

			assert(d != 0.0);
			x_ = x/d;
			y_ = y/d;
			z_ = z/d;
			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		pnt& operator/=(double d){/*{{{*/
			if(d == 0.0){
				std::cout << "pnt: operator /=" << std::endl << (*this) << std::endl;
			}
			assert(d != 0.0);
			x = x/d;
			y = y/d;
			z = z/d;
			return *this;
		}/*}}}*/
		pnt& operator+=(const pnt &p){/*{{{*/
			x += p.x;
			y += p.y;
			z += p.z;
			return *this;
		}/*}}}*/
		double operator[](int i) const {/*{{{*/
			if(i == 0){
				return x;
			} else if(i == 1){
				return y;
			} else {
				return z;
			}
		}/*}}}*/
		void normalize(){/*{{{*/
			double norm;

			norm = x*x + y*y + z*z;
			if(norm == 0){
				std::cout << "pnt: normalize" << std::endl;
				std::cout << x << " " << y << " " << z << " " << isBdry << " " << idx << std::endl;

				assert(norm != 0);
			}	
			norm = sqrt(norm);

			x = x/norm;
			y = y/norm;
			z = z/norm;
		}/*}}}*/
		double dot(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;

			return junk;
		}/*}}}*/
		double dotForDistance(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;

			junk = junk - 1.0;

			return fabs(junk);
		}/*}}}*/
		double dotForAngle(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;
			if(junk > 1.0){
				junk = 1.0;
			}

			if(junk < -1.0){
				junk = -1.0;
			}
			return acos(junk);
		}/*}}}*/
		pnt cross(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = y*p.z - p.y*z;
			y_ = z*p.x - p.z*x;
			z_ = x*p.y - p.x*y;

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		double magnitude() const{/*{{{*/
			return sqrt(x*x + y*y + z*z);
		}/*}}}*/
		double magnitude2() const {/*{{{*/
			return x*x + y*y + z*z;
		}/*}}}*/
		double getLat() const {/*{{{*/
			return asin(z);			
		}/*}}}*/
		double getLon() const {/*{{{*/
			double lon;

			lon = atan2(y,x);

			if(lon < 0){
				return 2.0 * M_PI + lon;
			} else {
				return lon;
			}
		}/*}}}*/
	struct hasher {/*{{{*/
		size_t operator()(const pnt &p) const {
			uint32_t hash; 
			size_t i, key[3] = { (size_t)p.x, (size_t)p.y, (size_t)p.z };
			for(hash = i = 0; i < sizeof(key); ++i) {
				hash += ((uint8_t *)key)[i];
				hash += (hash << 10);
				hash ^= (hash >> 6);
			}
			hash += (hash << 3);
			hash ^= (hash >> 11);
			hash += (hash << 15);
			return hash;
		}
	};/*}}}*/
	struct idx_hasher {/*{{{*/
		size_t operator()(const pnt &p) const {
			return (size_t)p.idx;
		}
	};/*}}}*/
};/*}}}*/
class tri {/*{{{*/
	private:
		friend class boost::serialization::access;
		template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
				ar & vi1;
				ar & vi2;
				ar & vi3;
				ar & idx;
			}
	public:
	int vi1, vi2, vi3;
	int idx;

	tri() : vi1(0), vi2(0), vi3(0), idx(0) { }

	tri(int vi1_, int vi2_, int vi3_)
		: vi1(vi1_), vi2(vi2_), vi3(vi3_), idx(0) { }

	tri(int vi1_, int vi2_, int vi3_, int idx_)
		: vi1(vi1_), vi2(vi2_), vi3(vi3_), idx(idx_) { }


	friend std::ostream & operator<<(std::ostream &os, const tri &t);
	friend std::istream & operator>>(std::istream &is, tri &t);

	tri sortedTri(){/*{{{*/
		int v1, v2, v3, swp_v;
		v1 = vi1;
		v2 = vi2;
		v3 = vi3;

		//Bubble sort on 3 integers
		for(int i = 0; i < 2; i++){
			if(v1 > v2){
				swp_v = v1;
				v1 = v2;
				v2 = swp_v;
			}
			if(v2 > v3){
				swp_v = v2;
				v2 = v3;
				v3 = swp_v;
			}
		}

		return tri(v1,v2,v3);
	}/*}}}*/
	bool operator==(const tri &t) const {/*{{{*/
		return (vi1 == t.vi1) & (vi2 == t.vi2) & (vi3 == t.vi3);
	}/*}}}*/
	struct hasher {/*{{{*/
		size_t operator()(const tri &t) const {
			uint32_t hash; 
			size_t i, key[3] = { t.vi1, t.vi2, t.vi3 };
			for(hash = i = 0; i < sizeof(key); ++i) {
				hash += ((uint8_t *)key)[i];
				hash += (hash << 10);
				hash ^= (hash >> 6);
			}
			hash += (hash << 3);
			hash ^= (hash >> 11);
			hash += (hash << 15);
			return hash;
		}
	};/*}}}*/
};/*}}}*/

inline pnt operator*(const double d, const pnt &p){/*{{{*/
	return pnt(d*p.x, d*p.y, d*p.z, 0, 0);
}/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const pnt &p){/*{{{*/
	//os << '(' << p.x << ", " << p.y << ", " << p.z << ") bdry=" << p.isBdry << " idx=" << p.idx << ' ';
	os  << std::setiosflags(std::ios::fixed) << std::setprecision(16) << p.x << " " << p.y << " " << p.z;
	return os;
}/*}}}*/
inline std::istream & operator>>(std::istream &is, pnt &p){/*{{{*/
	//is >> p.x >> p.y >> p.z >> p.isBdry >> p.idx;
	is >> p.x >> p.y >> p.z;
	return is;
}/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const tri &t){/*{{{*/
	//return os << t.vi1 << " " << t.vi2 << " " << t.vi3;
	return os << t.vi1 << " " << t.vi2 << " " << t.vi3;
}/*}}}*/
inline std::istream & operator>>(std::istream &is, tri &t){/*{{{*/
	//return is >> t.vi1 >> t.vi2 >> t.vi3;
	return is >> t.vi1 >> t.vi2 >> t.vi3;
}/*}}}*/

void circumcenter(const pnt &A,const pnt &B,const pnt &C, pnt &cent){/*{{{*/
	double a, b, c;
	double pbc, apc, abp;

	a = (B-C).magnitude2();
	b = (C-A).magnitude2();
	c = (A-B).magnitude2();

	pbc = a*(-a + b + c);
	apc = b*( a - b + c);
	abp = c*( a + b - c);

	cent = (pbc*A + apc*B + abp*C)/(pbc + apc + abp);

}/*}}}*/
double circumradius(const pnt &A, const pnt &B, const pnt &C){/*{{{*/

	pnt ccenter;  
	pnt ac;

	circumcenter(A,B,C,ccenter);
	ac = A - ccenter;

	return ac.magnitude();

/*
	pnt a(0.0,0.0,0.0,0);
	pnt b(0.0,0.0,0.0,0);
	pnt sub(0.0,0.0,0.0,0);
	pnt cross(0.0,0.0,0.0,0);
	double top, bot;

	a = A-C;
	b = B-C;

	sub = a-b;
	cross = a.cross(b);

	top = a.magnitude() * b.magnitude() * sub.magnitude();
	bot = 2.0 * cross.magnitude();

	return top/bot; // */

/*	dx = a.x - b.x;
	dy = a.y - b.y;
	dz = a.z - b.z;

	sa = sqrt(dx*dx + dy*dy + dz*dz);

	dx = b.x - c.x;
	dy = b.y - c.y;
	dz = b.z - c.z;

	sb = sqrt(dx*dx + dy*dy + dz*dz);

	dx = c.x - a.x;
	dy = c.y - a.y;
	dz = c.z - a.z;

	sc = sqrt(dx*dx + dy*dy + dz*dz);

	dx = (sa+sb+sc)/2.0; // Semiperimeter
	dy = 4.0*sqrt(dx*(sa+sb-dx)*(sa+sc-dx)*(sb+sc-dx)); // Bottom of circumradius computation

	return (sa*sb*sc)/dy;*/
}/*}}}*/
double triArea(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	/**************************************************************************
	 * - This function calculates the area of the triangle A,B,C on the
	 *   surface of a sphere.
	 *
	 *   Input: A, B, C
	 *        A: vertex 1 of triangle
	 *        B: vertex 2 of triangle
	 *        C: vertex 3 of triangle
	 *   Output: (returned value area)
	 *   	area: surface area of triangle on sphere.
	 **************************************************************************/
	pnt u12, u23, u31;
	double a, b, c, s, tanqe, area;	
	double sign;

	area = 0.0;

	//Compute Surface normal for triangle to get "sign" of area
	u12 = B - A;
	u23 = C - A;

	u31 = u12.cross(u23);

	u31.normalize();
	sign = u31.magnitude();
	assert(sign != 0.0);

	//dot the surface norm with one of the points to get sign.
	sign = u31.dot(A);
	sign = sign/fabs(sign);

	a = A.dotForAngle(B);
	b = B.dotForAngle(C);
	c = C.dotForAngle(A);

	s = 0.5*(a+b+c);

	tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)));

	area = sign*4.0*atan(tanqe);
//	area = 4.0*atan(tanqe);

	if(isnan(area))
		area = 0.0;

	return area;
}/*}}}*/
int isCcw(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	double sign;
	pnt ab, ac, cross;

	ab = B-A;
	ac = C-A;
	cross = ab.cross(ac);

	sign = cross.dot(A);
	if(sign > 0){
		return 1;
	} else {
		return 0;
	}
}/*}}}*/
pnt pntFromLatLon(const double &lat, const double &lon){/*{{{*/
	pnt temp;
	temp.x = cos(lon) * cos(lat);
	temp.y = sin(lon) * cos(lat);
	temp.z = sin(lat);
	return temp;
}/*}}}*/
