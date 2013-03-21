#ifndef UTIL_H_
#define UTIL_H_

#include <stdlib.h>

inline float random(float min, float max)
{
    return static_cast<float>(rand()) / RAND_MAX * (max - min) + min;
}

template<class T>
inline T sqr(const T &x){
    return x*x;
}

template<class T>
inline T min(const T &a1, const T &a2)
{
    if (a1 < a2) {
        return a1;
    } else {
        return a2;
    }
}

template<class T>
inline T min(const T &a1, const T &a2, const T &a3)
{
    if(a1<a2) {
        return min(a1,a3);
    } else {
        return min(a2,a3);
    }
}

template<class T>
inline T max(const T &a1, const T &a2)
{
    if (a1 > a2) {
        return a1;
    } else {
        return a2;
    }
}

template<class T>
inline T max(const T &a1, const T &a2, const T &a3)
{
    if(a1 > a2) {
        return max(a1,a3);
    } else {
        return max(a2,a3);
    }
}

template<class T>
inline T clamp(T a, T lower, T upper)
{
    if(a < lower) {
        return lower;
    } else if(a > upper) {
        return upper;
    } else {
        return a;
    }
}

#endif
