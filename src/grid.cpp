#include "grid.h"
#include "util.h"

Grid::Grid(Settings::Ptr s)
{
    _u.resize(s->nx,s->ny,s->dx);
    _v.resize(s->nx,s->ny,s->dx);
    _uSum.resize(s->nx,s->ny,s->dx);
    _vSum.resize(s->nx,s->ny,s->dx);
    _uWeights.resize(s->nx,s->ny,s->dx);
    _vWeights.resize(s->nx,s->ny,s->dx);
}

void Grid::sampleVelocities(Particles::Ptr p)
{
    _reset();
    size_t i,j;
    float tx,ty;
    for (int idx = 0; idx < p->numParticles(); ++idx) {
        _u.bary(p->pos(idx).x, p->pos(idx).y, i, j, tx, ty);
        _accumulate(_u, _uSum, _uWeights, p->vel(idx).x, i, j, tx, ty);
        _v.bary(p->pos(idx).x, p->pos(idx).y, i, j, tx, ty);
        _accumulate(_v, _vSum, _vWeights, p->vel(idx).y, i, j, tx, ty);
    }
    _u.divide(_uSum);
    _v.divide(_vSum);
}

void Grid::applyGravity(const Vec2f & g, float dt)
{
    _u.add(g.x * dt);
    _v.add(g.y * dt);
}

float Grid::CFL() const
{
    float x = sqr(_u.infNorm()) + sqr(_v.infNorm());
    if (x < 1e-16){
        x = 1e-16;
    }
    return _u.dx() / sqrt(x);
}

void Grid::extrapolateVelocities(FluidSDF::Ptr f, int numSweepIterations)
{
    for (int i = 0; i < numSweepIterations; ++i) {
        _sweep<RIGHT, true>(f, true, true);
        _sweep<RIGHT, true>(f, true, false);
        _sweep<LEFT, true>(f, false, true);
        _sweep<LEFT, true>(f, false, false);

        _sweep<TOP, false>(f, true, true);
        _sweep<TOP, false>(f, true, false);
        _sweep<BOTTOM, false>(f, false, true);
        _sweep<BOTTOM, false>(f, false, false);
     }
}

void Grid::pressureProjection(const Array2f & p,
                              FluidSDF::Ptr f,
                              float dt)
{
    float scale = dt / p.dx();
    float theta;
    _u.reset();
    for (int i = 1; i < p.nx(); ++i) {
        for (int j = 0; j < p.ny(); ++j) {
            const float uw = _uWeights.face<LEFT>(i,j);
            if (_theta(uw, f->phi(i,j), f->phi(i-1,j), theta)) {
                _u.face<LEFT>(i,j) -= scale * (p(i,j) - p(i-1,j)) / theta;
            }
        }
    }

    _v.reset();
    for (int i = 0; i < p.nx(); ++i) {
        for (int j = 1; j < p.ny(); ++j) {
            const float vw = _vWeights.face<BOTTOM>(i,j);
            if (_theta(vw, f->phi(i,j), f->phi(i,j-1), theta)) {
                _v.face<BOTTOM>(i,j) -= scale * (p(i,j) - p(i,j-1)) / theta;
            }
        }
    }
}

void Grid::enforceBoundaryConditions(SolidSDF::Ptr s)
{
    _uSum.reset();
    _vSum.reset();
        
    for (int i = 1; i < _u.nx(); ++i) {
        for (int j = 0; j < _u.ny(); ++j) {
            if (_uWeights.face<LEFT>(i,j) < 0) {
                const Vec2f pos = _uWeights.pos<LEFT>(i,j);
                Vec2f gradient = s->gradient(pos);
                gradient.normalize();
                const Vec2f vel(_u.face<LEFT>(i,j), _v.bilerp(pos));
                _uSum.face<LEFT>(i,j) = vel.x - gradient.dot(vel);
            } else {
                _uSum.face<LEFT>(i,j) = _u.face<LEFT>(i,j);
            }
        }
    }
    for (int i = 0; i < _v.nx(); ++i) {
        for (int j = 1; j < _v.ny(); ++j) {
            if (_vWeights.face<BOTTOM>(i,j) < 0) {
                const Vec2f pos = _vWeights.pos<BOTTOM>(i,j);
                Vec2f gradient = s->gradient(pos);
                gradient.normalize();
                const Vec2f vel(_u.bilerp(pos), _v.face<BOTTOM>(i,j));
                _vSum.face<BOTTOM>(i,j) = vel.y - gradient.dot(vel);
            } else {
                _vSum.face<BOTTOM>(i,j) = _v.face<BOTTOM>(i,j);
            }
        }
    }

    _u.swap(_uSum);
    _v.swap(_vSum);
}

template<Faces T_FACE, bool T_U>
void Grid::_sweep(FluidSDF::Ptr f, bool upsweepX, bool upsweepY)
{
    const int i0 = upsweepX ? 0 : _u.nx() - 1;
    const int i1 = upsweepX ? _u.nx() : -1;
    const int j0 = upsweepY ? 0 : _u.ny() - 1;
    const int j1 = upsweepY ? _u.ny() : -1;
    const int di = i0 < i1 ? 1 : -1;
    const int dj = j0 < j1 ? 1 : -1;
    for (int j = j0; j != j1; j += dj) {
        for (int i = i0; i != i1; i += di) {
            if (T_U ? _uWeights.face<T_FACE>(i,j) :
                      _vWeights.face<T_FACE>(i,j)) {
                const Vec2f &p =T_U ? _u.pos<T_FACE>(i,j) : _v.pos<T_FACE>(i,j);
                const Vec2f grad = f->gradient(p);
                // Only interested in upwinding. If any of the derivates is
                // negative it means its propegating in the wrong direction
                if (grad.x < 0.0f || grad.y < 0.0f) {
                    continue;
                }

                // Interpolate between the derivates. 
                // Special case if the denominator is zero.
                const float sum = grad.x + grad.y;
                const float a = sum ? grad.x / sum : 0.5f;

                if (T_U) {
                     _u.face<T_FACE>(i,j) = a * _u.face<T_FACE>(i-di,j) +
                                            (1.0f-a)  * _u.face<T_FACE>(i,j-dj);
                } else {
                    _v.face<T_FACE>(i,j) = a * _v.face<T_FACE>(i-di,j) +
                                            (1.0f-a)  * _v.face<T_FACE>(i,j-dj);
                }
            }
        }
    }
}

template <typename T_ARRAY>
void Grid::_accumulate(T_ARRAY & array,
                       T_ARRAY & sum,
                       T_ARRAY & marker,
                       float q,
                       int i,
                       int j,
                       float tx,
                       float ty)
{
    float weight = (1.0f - tx) * (1.0f - ty);
    array(i,j) += weight * q;
    sum(i,j) += weight;
    marker(i,j) = 1;

    weight = tx * (1.0f - ty);
    array(i + 1,j) += weight * q;
    sum(i + 1,j) += weight;
    marker(i + 1,j) = 1;

    weight = (1.0f - tx) * ty;
    array(i,j + 1) += weight * q;
    sum(i,j + 1) += weight;
    marker(i,j + 1) = 1;

    weight = tx * ty;
    array(i + 1,j + 1) += weight * q;
    sum(i + 1,j + 1) += weight;
    marker(i + 1,j + 1) = 1;
}

void Grid::_reset()
{
    _u.reset();
    _uSum.reset();
    _uWeights.reset();
    _v.reset();
    _vSum.reset();
    _vWeights.reset();
}
