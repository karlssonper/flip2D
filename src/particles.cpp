#include "particles.h"

Particles::Particles()
{
    
}

void Particles::addParticle(const Vec2f & pos, Vec2f vel)
{
    _pos.push_back(pos);
    _vel.push_back(vel);
}

void Particles::addParticles(const std::vector<Vec2f> & pos,
                             const std::vector<Vec2f> & vel)
{
    _pos = pos;
    _vel = vel;
}

void Particles::updateVelocities(const Array2f & u, const Array2f & v)
{
    for(int i = 0; i < _pos.size(); ++i) {
        _vel[i] = Vec2f(u.bilerp(_pos[i]), v.bilerp(_pos[i]));
    }
}

void Particles::advect(float dt)
{
    for(int i = 0; i < _pos.size(); ++i) {
        _pos[i] += dt * _vel[i];
    }
}
