#include "grid.h"
#include "util.h"

Grid::Grid(int nx, int ny, float dx)
{
    _u.resize(nx,ny,dx);
    _v.resize(nx,ny,dx);
    _uSum.resize(nx,ny,dx);
    _vSum.resize(nx,ny,dx);
    _uWeights.resize(nx,ny,dx);
    _vWeights.resize(nx,ny,dx);
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
    
}

template<int T_FACEDIR>
void Grid::_sweep(FluidSDF::Ptr f)
{
    const int di = i0 < i1 ? 1 : -1;
    const int dj = j0 < j1 ? 1 : -1;
    for (int j = j0; j != j1; j += dj) {
        for (int i = i0; i != i1; i += di) {
            if (!marker<T_FACEDIR>(i,j)) {



                
            }
        }
    }
}

template<bool T_SWEEP_X>
void Grid::_sweepVel(Array2f & vel, int i0, int i1, int j0, int j1)
{
    // Sweep direction variables
    const int si = T_SWEEP_X ? -1 : 0, sj = T_SWEEP_X ? 0: -1;
    const int di = i0 < i1 ? 1 : -1, dj = j0 < j1 ? 1 : -1;
    for (int j = j0; j != j1; j += dj) {
        for (int i = i0; i != i1; i += di) {
            // Only sweep in cells not covered by particles
            if (!marker(i,j) && !marker(i + si, j + sj)) {
                // Spatial partial derivatives of the signed distance
                const float dphidx = T_SWEEP_X ? 
                        di*(_phi(i,j)-_phi(i-1,j)) :
                        0.5f*(_phi(i,j)+_phi(i,j-1)-_phi(i-di,j)-_phi(i-di,j-1));
                const float dphidy = !T_SWEEP_X ? 
                        dj*(_phi(i,j)-_phi(i,j-1)) :
                        0.5f*(_phi(i,j)+_phi(i-1,j)-_phi(i,j-dj)-_phi(i-1,j-dj));

                // Only interested in upwinding. If any of the derivates is
                // negative it means its propegating in the wrong direction
                if (dphidx < 0.0f || dphidy < 0.0f) {
                    continue;
                }

                // Interpolate between the derivates. 
                // Special case if the denominator is zero.
                const float sum = dphidx + dphidy;
                const float alpha = sum ? dphidx / sum : 0.5f;
                vel(i,j) = alpha * vel(i-di,j) + (1.0f-alpha) * vel(i,j-dj);
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
