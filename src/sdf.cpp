#include "sdf.h"
#include <cassert>
#include <cmath>

SolidSDF::SolidSDF()
{
    
}

SolidSDF::SolidSDF(int nx, int ny, float dx)
{
    _phi.resize(nx,ny,dx);
}

void SolidSDF::createWeights(FaceArray2Xf & uw, FaceArray2Yf & vw)
{
    int iEnd = uw.nx() - 1;
    for (int j = 0; j < uw.ny(); ++j) {
        for (int i = 0; i < uw.nx(); ++i) {
            uw.face<LEFT>(i,j) = _weight(_phi.corner<TOP_LEFT>(i,j),
                                         _phi.corner<BOTTOM_LEFT>(i,j));
        }
        uw.face<RIGHT>(iEnd,j) = _weight(_phi.corner<TOP_RIGHT>(iEnd,j),
                                         _phi.corner<BOTTOM_RIGHT>(iEnd,j));
    }
        
    int jEnd = vw.ny() - 1;
    for (int i = 0; i < vw.nx(); ++i) {
        for (int j = 0; j < vw.ny(); ++j) {
            vw.face<BOTTOM>(i,j) = _weight(_phi.corner<BOTTOM_RIGHT>(i,j),
                                           _phi.corner<BOTTOM_LEFT>(i,j));
        }
        vw.face<TOP>(i,jEnd) = _weight(_phi.corner<TOP_RIGHT>(i,jEnd),
                                       _phi.corner<TOP_LEFT>(i,jEnd));
    }
}

void SolidSDF::initBoxBoundary(int width)
{
    _phi.set(-1.5f + width + 3);
    float x = -1.5f;
    int m = 0;
    while (m < width + 3) {
        for (int i = m; i < _phi.nx() - m; ++i) {
            _phi.corner<BOTTOM_LEFT>(i,m) = x * _phi.dx();
            _phi.corner<TOP_RIGHT>(i,_phi.ny()-1-m) = x * _phi.dx();
        }
    
        for (int j = m; j < _phi.ny() - m; ++j) {
            _phi.corner<TOP_LEFT>(m,j) = x * _phi.dx();;
            _phi.corner<BOTTOM_RIGHT>(_phi.nx()-1-m,j) = x * _phi.dx();;
        }
        ++m;
        x+= 1.0f;
    }
}

FluidSDF::FluidSDF()
{
}

FluidSDF::FluidSDF(int nx, int ny, float dx)
{
    _phi.resize(nx,ny,dx);
    _sum.resize(nx,ny,dx);
    _pAvg.resize(nx,ny,dx);
}

void FluidSDF::reconstructSurface(Particles::Ptr particles, float R, float r)
{
    assert(R && r);
    _sum.reset();
    _pAvg.reset();
    
    int nx = _phi.nx();
    int ny = _phi.ny();
    
    for (int p = 0; p < particles->numParticles(); ++p) {
        const int ci = particles->pos(p).x / _phi.dx();
        const int cj = particles->pos(p).y / _phi.dx();
        assert(ci < _phi.nx() && cj < _phi.ny());

        for (int i = std::max(0,ci-2); i < std::min(nx,ci+2); ++i) {
            for (int j = std::max(0,cj-2); j < std::min(ny,cj+2); ++j) {
                const Vec2f xd = particles->pos(p) - _phi.pos(i,j);
                const float k = _kernel(xd.length() / R);
                _sum(i,j) +=  k;
                _pAvg(i,j) += k * particles->pos(p);
            }
        }
    }
    
    for (int i = 0; i < _phi.nx(); ++i) {
        for (int j = 0; j < _phi.ny(); ++j) {
            if (_sum(i,j) > 0) {
                _pAvg(i,j) *= (1.0f / _sum(i,j));
                const Vec2f xd = _phi.pos(i,j) - _pAvg(i,j);
                _phi(i,j) = xd.length() - r;
            } else {
                _phi(i,j) = _phi.dx() * 1e15;
            }
        }
    }
}

void FluidSDF::reinitialize(int numSwepIterations)
{
    for (int i = 0; i < numSwepIterations; ++i) {
        _sweep(1,_phi.nx(), 1,_phi.ny());
        _sweep(_phi.nx() - 2, -1,_phi.ny() - 2, -1);
        _sweep(1, _phi.nx(), _phi.ny() - 2, -1);
        _sweep(_phi.nx()-2, -1, 1, _phi.ny());
    }
}

void FluidSDF::_sweep(int i0, int i1, int j0, int j1)
{
    const int di = i0 < i1 ? 1 : -1;
    const int dj = j0 < j1 ? 1 : -1;
    for (int i = i0; i != i1; i+= di) {
        for (int j = j0; j != j1; j+= dj) {
            if (!isFluid(i,j)) {
                _solveEikonal(_phi(i - di,j), _phi(i,j - dj), _phi(i,j));
            }
        }
    }
}

void FluidSDF::extrapolateIntoSolid(SolidSDF::Ptr solid)
{
    for (int i = 0; i < _phi.nx(); ++i) {
        for (int j = 0; j < _phi.ny(); ++j) {
            if (_phi(i,j) < 0.5 * _phi.dx() && solid->center(i,j) < 0) {
                _phi(i,j) = -0.5f * _phi.dx();
            }
        }
    }
}
