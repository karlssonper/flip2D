#ifndef ARRAY_H_
#define ARRAY_H_

#include "vec2.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cmath>

template<typename T>
class Array2
{
  public:
    Array2() : _nx(0), _ny(0), _dx(1.0) {}
    
    Array2(size_t nx, size_t ny, float dx) :
            _nx(nx), _ny(ny), _dx(dx)
    {
        _data.resize(_nx*_ny);
    }

    const T & operator()(size_t i, size_t j) const
    {
        assert(i < _nx);
        assert(j < _ny);
        return _data[_idx(i,j)];
    }

    T & operator()(size_t i, size_t j)
    {
        assert(i < _nx);
        assert(j < _ny);
        return _data[_idx(i,j)];
    }

    Vec2f pos(size_t i, size_t j) const
    {
        return _dx * Vec2<T>(i + 0.5, j + 0.5);
    }

    void bary(float x,
              float y,
              size_t & i,
              size_t & j,
              float & tx,
              float & ty) const
    {
        // Scale from world coordinates to index coordinates
        x/= _dx;
        y/= _dx;

        if (x < 0.5) {
            i = 0;
            tx = 0;
        } else if (x > _nx - 0.5) {
            i = _nx - 2;
            tx = 1;
        } else {
            i = floorf(x - 0.5);
            tx = x - 0.5 - i;
        }

        if (y < 0.5) {
            j = 0;
            ty = 0;
        } else if (y > _ny - 0.5) {
            j = _ny - 2;
            ty = 1;
        } else {
            j = floorf(y - 0.5);
            ty = y - 0.5 - j;
        }
    }

    T bilerp(const Vec2f & pos) const
    {
        return bilerp(pos.x, pos.y);
    }
    
    T bilerp(float x, float y) const
    {
        size_t i, j;
        float tx, ty;
        bary(x, y, i, j, tx, ty);
        return bilerp(i,j,tx,ty);
    }

    T bilerp(size_t i, size_t j, float tx, float ty) const
    {
        return (1-tx) * ((1-ty) * (*this)(i,j) + 
               ty * (*this)(i,j+1)) + 
               tx * ((1-ty) * (*this)(i+1,j) +
               ty * (*this)(i+1,j+1)); 
    }

    void reset()
    {
        for (size_t i = 0; i < _data.size(); ++i) {
            _data[i] = 0;
        }
    }

    void resize(size_t nx, size_t ny, float dx = 1.0)
    {
        _nx = nx;
        _ny = ny;
        _dx = dx;
        _data.resize(_nx*_ny);
    }

    void copy(const Array2<T> & src)
    {
        assert(src._data.size() == _data.size());
        std::copy(src._data.begin(),
                  src._data.end(),
                  _data.begin());
    }

    void swap(Array2<T> & src)
    {
        _data.swap(src._data);
    }

    T dot(const Array2<T> & rhs)
    {
        assert(rhs._data.size() == _data.size());
        T sum = 0;
        for (size_t i = 0; i < _data.size(); ++i) {
            sum += _data[i] * rhs._data[i];
        }
        return sum;
    }

    void set(T value)
    {
        for (size_t i = 0; i < _data.size(); ++i) {
            _data[i] = value;
        }
    }

    void add(T t)
    {
        for (size_t i = 0; i < _data.size(); ++i) {
            _data[i] += t;
        }
    }
    
    void add(const Array2<T> & rhs, T scale = 1.0)
    {
        for (size_t i = 0; i < _data.size(); ++i) {
            _data[i] += rhs._data[i] * scale;
        }
    }

    void scaleAndAdd(T scale, const Array2<T> & addArray)
    {
        for (size_t i = 0; i < _data.size(); ++i) {
            _data[i] = _data[i] * scale + addArray._data[i];
        }
    }

    void multiply(const Array2<T> & rhs)
    {
        for (size_t i = 0; i < _data.size(); ++i) {
                _data[i] *= rhs._data[i];
        }
    }

    void multiply(T  rhs)
    {
        for (size_t i = 0; i < _data.size(); ++i) {
            _data[i] *= rhs;
        }
    }
    
    void divide(const Array2<T> & rhs)
    {
        for (size_t i = 0; i < _data.size(); ++i) {
            if (rhs._data[i]) {
                _data[i] /= rhs._data[i];
            }
        }
    }

    /**
        Return the magnitude of the largest or smallest value
    */
    T infNorm() const
    {
        T norm = 0;
        for (size_t i = 0; i < _data.size(); ++i) {
            if (std::fabs(_data[i]) > norm) {
                norm = std::fabs(_data[i]);
            }
        }
        return norm;
    }

    void writeMatlab(const char * filename, int frame) const
    {
        std::stringstream ss;
        ss << filename << frame << ".txt";
        std::ofstream file(ss.str().c_str());
        if (file.is_open()) {
            for(int j = _ny - 1; j >= 0; --j) {
                for(int i = 0; i < _nx - 1; ++i) {
                    file << (double)(*this)(i,j) << ", ";
                }
                file << (double)(*this)(_nx - 1,j) << "\n";
            }
        } else {
            std::cerr << "Can't write matlab file " << filename << std::endl;
        }
    }
    
    size_t nx() const { return _nx; }
    size_t ny() const { return _ny; }
    float dx() const { return _dx; }
    
  protected:
    size_t _nx, _ny;
    float _dx;
    std::vector<T> _data;

    size_t _idx(size_t i, size_t j) const
    {
        assert(i < _nx);
        assert(j < _ny);
        return j + _ny * i;
    }
};

typedef Array2<float> Array2f;
typedef Array2<double> Array2d;
typedef Array2<char> Array2c;
typedef Array2<int> Array2i;

enum Faces
{
    CENTER = 0,
    TOP = 1,
    BOTTOM = 2,
    LEFT = 3,
    RIGHT = 4,
};

template <typename T, size_t T_DIMX, size_t T_DIMY>
class FaceArray2 : public Array2<T>
{
  public:
    FaceArray2() {}
    FaceArray2(size_t nx, size_t ny, float dx) :
            Array2<T>(nx + T_DIMX, ny + T_DIMY, dx)
    {
        assert(T_DIMX || T_DIMY);
    }

    Vec2f pos(size_t i, size_t j) const
    {
        return _dx * Vec2<T>(i + 0.5*T_DIMY, j + 0.5*T_DIMX);
    }
    
    template<Faces T_FACE>
    Vec2<float> pos(size_t i, size_t j) const
    {
        if (T_DIMX) {
            assert(T_FACE == LEFT || T_FACE == RIGHT);
        } else if (T_DIMY) {
            assert(T_FACE == BOTTOM || T_FACE == TOP);
        }
        
        if (T_FACE == LEFT) {
            return _dx * Vec2<T>(i + 0.5 * T_DIMY, j + 0.5 * T_DIMX);
        } else if (T_FACE == RIGHT) {
            return _dx * Vec2<T>(i + 1 + 0.5 * T_DIMY, j + 0.5 * T_DIMX);
        } else if (T_FACE == TOP) {
            return _dx * Vec2<T>(i + 0.5 * T_DIMY, j + 1 + 0.5 * T_DIMX);
        } else if (T_FACE == BOTTOM) {
            return _dx * Vec2<T>(i + 0.5 * T_DIMY, j + 0.5 * T_DIMX);
        }
    }

    template<Faces T_FACE>
    T & face(size_t i, size_t j)
    {
        if (T_FACE == LEFT ||
            T_FACE == BOTTOM) {
            return _data[_idx(i,j)];
        } else if (T_FACE == RIGHT ||
                   T_FACE == TOP) {
            return _data[_idx(i + T_DIMX,j + T_DIMY)];
        }   
    }
    
    template<Faces T_FACE>
    T face(size_t i, size_t j) const
    {        
        if (T_FACE == LEFT ||
            T_FACE == BOTTOM) {
            return _data[_idx(i,j)];
        } else if (T_FACE == RIGHT ||
                   T_FACE == TOP) {
            return _data[_idx(i + T_DIMX,j + T_DIMY)];
        }  
    }

    T center(size_t i, size_t j) const
    {
        return 0.5*(_data[_idx(i,j)] +
                    _data[_idx(i+T_DIMX,j+T_DIMY)]);
        
    }

    void bary(float x,
              float y,
              size_t & i,
              size_t & j,
              float & tx,
              float & ty) const
    {
        // Scale from world coordinates to index coordinates
        x/= _dx;
        y/= _dx;

        if (x < 0.5 * T_DIMY) {
            i = 0;
            tx = 0;
        } else if (x >= _nx - T_DIMX - 0.5 * T_DIMY) {
            i = _nx - 2;
            tx = 1;
        } else {
            i = floorf(x - 0.5 * T_DIMY);
            tx = x - 0.5 * T_DIMY - i;
        }

        if (y < 0.5 * T_DIMX) {
            j = 0;
            ty = 0;
        } else if (y >= _ny - T_DIMY - 0.5 * T_DIMX) {
            j = _ny - 2;
            ty = 1;
        } else {
            j = floorf(y - 0.5 * T_DIMX);
            ty = y - 0.5 * T_DIMX - j;
        }
    }

    T bilerp(Vec2f pos) const
    {
        return bilerp(pos.x, pos.y);
    }
    
    T bilerp(float x, float y) const
    {
        size_t i,j;
        float tx,ty;
        bary(x, y, i, j, tx, ty);
        return Array2<T>::bilerp(i, j, tx, ty);
    }
    
    void resize(size_t nx, size_t ny, float dx = 1.0)
    {
        Array2<T>::resize(nx + T_DIMX, ny + T_DIMY, dx);
    }

    size_t nx() const { return _nx - 1 * T_DIMX; }
    size_t ny() const { return _ny - 1 * T_DIMY; }
    
  protected:
    // To make it compile in earlier gcc versions.
    // This bug is fixed in 4.7 and forward it seems
    using Array2<T>::_data;
    using Array2<T>::_idx;
    using Array2<T>::_nx;
    using Array2<T>::_ny;
    using Array2<T>::_dx;
};

typedef FaceArray2<float,1,0> FaceArray2Xf;
typedef FaceArray2<float,0,1> FaceArray2Yf;
typedef FaceArray2<double,1,0> FaceArray2Xd;
typedef FaceArray2<double,0,1> FaceArray2Yd;
typedef FaceArray2<char,1,0> FaceArray2Xc;
typedef FaceArray2<char,0,1> FaceArray2Yc;
typedef FaceArray2<int,1,0> FaceArray2Xi;
typedef FaceArray2<int,0,1> FaceArray2Yi;

enum Corner
{
    BOTTOM_LEFT, BOTTOM_RIGHT, TOP_LEFT, TOP_RIGHT
};

template<typename T>
class CornerArray2 : public Array2<T>
{
  public:
    CornerArray2() {} 
    CornerArray2(size_t nx, size_t ny, float dx)
            : Array2<T>(nx+1,nx+1,dx) {}

    Vec2f pos(size_t i, size_t j) const
    {
        return _dx * Vec2<T>(i, j);
    }
    
    template<Corner T_CORNER>
    Vec2<float> pos(size_t i, size_t j) const
    {
        if (T_CORNER == BOTTOM_LEFT) {
            return _dx * Vec2<T>(i, j);
        } else if (T_CORNER == BOTTOM_RIGHT) {
            return _dx * Vec2<T>(i+1, j);
        } else if (T_CORNER == TOP_LEFT) {
            return _dx * Vec2<T>(i, j+1);
        } else if (T_CORNER == TOP_RIGHT) {
            return _dx * Vec2<T>(i+1, j+1);
        }
    }

    T center(size_t i, size_t j) const
    {
        return  0.25 * (_data[_idx(i,j)] +
                        _data[_idx(i+1,j)] +
                        _data[_idx(i,j+1)] +
                        _data[_idx(i+1,j+1)]);
    }
    
    template<Corner T_CORNER>
    T & corner(size_t i, size_t j)
    {
        assert(i < _nx - 1);
        assert(j < _ny - 1);
        if (T_CORNER == BOTTOM_LEFT) {
            return _data[_idx(i,j)];
        } else if (T_CORNER == BOTTOM_RIGHT) {
            return _data[_idx(i+1,j)];
        } else if (T_CORNER == TOP_LEFT) {
            return _data[_idx(i,j+1)];
        } else if (T_CORNER == TOP_RIGHT) {
            return _data[_idx(i+1,j+1)];
        }
    }

    template<Corner T_CORNER>
    T corner(size_t i, size_t j) const
    {
        assert(i < _nx - 1);
        assert(j < _ny - 1);
        if (T_CORNER == BOTTOM_LEFT) {
            return _data[_idx(i,j)];
        } else if (T_CORNER == BOTTOM_RIGHT) {
            return _data[_idx(i+1,j)];
        } else if (T_CORNER == TOP_LEFT) {
            return _data[_idx(i,j+1)];
        } else if (T_CORNER == TOP_RIGHT) {
            return _data[_idx(i+1,j+1)];
        }
    }

    void bary(float x,
              float y,
              size_t & i,
              size_t & j,
              float & tx,
              float & ty) const
    {
        // Scale from world coordinates to index coordinates
        x/= _dx;
        y/= _dx;
            
        if (x < 0) {
            i = 0;
            tx = 0;
        } else if (x >= _nx - 1) {
            i = _nx - 2;
            tx = 1;
        } else {
            i = floorf(x);
            tx = x - i;
        }

        if (y < 0) {
            j = 0;
            ty = 0;
        } else if (y >= _ny - 1) {
            j = _ny - 2;
            ty = 1;
        } else {
            j = floorf(y);
            ty = y - j;
        }
    }

    T bilerp(const Vec2f & pos) const
    {
        return bilerp(pos.x, pos.y);
    }
    
    T bilerp(float x, float y) const
    {
        size_t i,j;
        float tx,ty;
        bary(x, y, i, j, tx, ty);
        return Array2<T>::bilerp(i,j,tx,ty);
    }
    
    void resize(size_t nx, size_t ny, float dx)
    {
        Array2<T>::resize(nx + 1, ny + 1, dx);
    }

    size_t nx() const { return _nx - 1; }
    size_t ny() const { return _ny - 1; }
        
  protected:
    // To make it compile in earlier gcc versions.
    // This bug is fixed in 4.7 and forward it seems
    using Array2<T>::_data;
    using Array2<T>::_idx;
    using Array2<T>::_nx;
    using Array2<T>::_ny;
    using Array2<T>::_dx;
};

typedef CornerArray2<float> CornerArray2f;
typedef CornerArray2<double> CornerArray2d;
typedef CornerArray2<char> CornerArray2c;
typedef CornerArray2<int> CornerArray2i;


#endif
