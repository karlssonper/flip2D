#ifndef SPARSE_H_
#define SPARSE_H_

#include "array.h"

template <typename T>
class SparseLaplacianMatrix
{
  public:
    SparseLaplacianMatrix();

    SparseLaplacianMatrix(size_t nx, size_t ny)
    {
        resize(nx,ny);
    }

    void resize(size_t nx, size_t ny)
    {
        _Adiag.resize(nx,ny);
        _Aplusi.resize(nx,ny);
        _Aplusj.resize(nx,ny);
    }

    void reset()
    {
        _Adiag.reset();
        _Aplusi.reset();
        _Aplusj.reset();
    }
    
    template<int T_ELEMENT>
    T& value(size_t i, size_t j)
    {
        if (T_ELEMENT == CENTER) {
            return _Adiag(i,j);
        } else if (T_ELEMENT == RIGHT) {
            return _Aplusi(i,j);
        } else if (T_ELEMENT == LEFT) {
            return _Aplusi(i-1,j);
        } else if (T_ELEMENT == BOTTOM) {
            return _Aplusj(i,j-1);
        } else if (T_ELEMENT == TOP) {
            return _Aplusj(i,j);
        }
    }

    template<int T_ELEMENT>
    T value(size_t i, size_t j) const
    {
        if (T_ELEMENT == CENTER) {
            return _Adiag(i,j);
        } else if (T_ELEMENT == RIGHT) {
            return _Aplusi(i,j);
        } else if (T_ELEMENT == LEFT) {
            return _Aplusi(i-1,j);
        } else if (T_ELEMENT == TOP) {
            return _Aplusj(i,j);
        } else if (T_ELEMENT == BOTTOM) {
            return _Aplusj(i,j-1);
        }
    }

    T mult(const Array2<T> & input, size_t i, size_t j) const
    {
        return value<CENTER>(i,j) * input(i,j) +
                (i > 0 ? value<LEFT>(i,j) * input(i-1,j) : 0) +
                (i < input.nx()-1 ? value<RIGHT>(i,j) * input(i+1,j) : 0) +
                (j > 0 ? value<BOTTOM>(i,j) * input(i,j-1) : 0) +
                (j < input.ny()-1 ? value<TOP>(i,j)* input(i,j+1) : 0);
    }
  protected:
    Array2<T> _Adiag;
    Array2<T> _Aplusi;
    Array2<T> _Aplusj;
};

#endif
