#ifndef VEC2_H_
#define VEC2_H_

#include <iostream>

template<typename T>
class Vec2
{
  public:
    T x,y;

    Vec2(T xi = 0, T yi = 0) : x(xi), y(yi) { };
    Vec2(const Vec2 & v) : x(v.x), y(v.y) { };

    T dot(const Vec2<T> rhs)
    {
        return x * rhs.x + y * rhs.y;
    }

    void normalize()
    {
        const T L = length();
        if (L) {
            (*this) * (1.0/L);
        }
    }

    T length() const
    {
        return sqrt(x * x + y * y);
    }

    void operator=(const Vec2<T> rhs)
    {
        x = rhs.x;
        y = rhs.y;
    }

    Vec2<T> operator+(const Vec2<T> & rhs)
    {
        return Vec2<T>(x + rhs.x, y + rhs.y);
    }

    void operator+=(const Vec2<T> & rhs)
    {
        x += rhs.x;
        y += rhs.y;
    }

    Vec2<T> operator-(const Vec2<T> & rhs) const
    {
        return Vec2<T>(x - rhs.x, y - rhs.y);
    }

    Vec2<T> operator*(T rhs) const
    {
        return Vec2<T>(x * rhs, y * rhs);
    }

    void operator*=(T rhs)
    {
        x*= rhs;
        y*= rhs;
    }
    
  protected:
};

typedef Vec2<float> Vec2f;
typedef Vec2<double> Vec2d;
typedef Vec2<unsigned int> Vec2i;

template<typename T>
inline Vec2<T> operator*(T lhs, const Vec2<T> & rhs)
{
    return Vec2<T>(lhs * rhs.x, lhs * rhs.y);
}

template<typename T>
std::ostream & operator << (std::ostream & os, const Vec2<T> & v)
{
    return os << "(" << v.x << "," << v.y << ")";
}

#endif
