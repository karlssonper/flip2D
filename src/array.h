#ifndef ARRAY_H_
#define ARRAY_H_

#include "vec2.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>

template<typename T>
class Array2
{
  public:
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

    Vec2<T> pos(size_t i, size_t j) const
    {
        return _dx * Vec2<T>(i + 0.5, j + 0.5);
    }

    
    T operator()(float x, float y) const
    {
        // Check if on or outside the border
        if (x < 0.5 * _dx ||
            y < 0.5 * _dx ||
            x > _dx * _nx ||
            y > _dx * _ny) {

            int i,j;
            if (x < 0.5 * _dx) {
                i = 0;
            }


            
        } else {
            T ii = x * _dx;
            T jj = y * _dx;
            int i = floorf(ii);
            int j = floorf(jj);
            return bilerp(i, j, ii - i, jj - j);
        }
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
        for (size_t i = 0; i < _size; ++i) {
            _data[i] = 0;
        }
    }

    void resize(size_t nx, size_t ny)
    {
        _nx = nx;
        _ny = ny;
        _data.resize(_nx*_ny);
    }

    void copy(const Array2<T> & src)
    {
        assert(src._data.size() == _data.size());
        std::copy(src._data.begin(), src._data.end(), _data.begin());
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

    /**
        Return the magnitude of the largest or smallest value
    */
    T infNorm() const
    {
        T norm = 0;
        for (size_t i = 0; i < _size; ++i) {
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
    
  protected:
    size_t _nx, _ny;
    float _dx;
    std::vector<T> _data;

    size_t _idx(size_t i, size_t j) const
    {
        return j + _ny * i;
    }
};

typedef Array2<float> Array2f;
typedef Array2<double> Array2d;
typedef Array2<char> Array2c;
typedef Array2<int> Array2i;

#endif
