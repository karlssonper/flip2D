#include "particles.h"
#include "util.h"
#include "log.h"

Particles::Particles()
{
    
}

void Particles::initSphere(const Array2f & solidPhi,
                           const Vec2f & center,
                           float radius,
                           int particlesPerCell,
                           Vec2f vel)
{
    LOG_OUTPUT("Initating the fluid as a sphere at " << vel);
    const float r2 = sqr(radius);
    for (int i = 0; i < solidPhi.nx(); ++i) {
        for (int j = 0; j < solidPhi.ny(); ++j) {
            for (int n = 0; n < particlesPerCell; ++n) {
                const Vec2f offset(random(-0.495, 0.495),random(-0.495, 0.495));
                const Vec2f pos = solidPhi.pos(i,j) + solidPhi.dx() * offset;
                const Vec2f d = pos - center;

                // Inside test
                if (sqr(d.x) + sqr(d.y) - r2 < 0 && solidPhi.bilerp(pos) >= 0) {
                    addParticle(pos, vel);
                }
            }
        }
    }
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

void Particles::updateVelocities(const FaceArray2Xf & u, const FaceArray2Yf & v)
{
    LOG_OUTPUT("Updating particle velocities from grid.");
    for(int i = 0; i < _pos.size(); ++i) {
        _vel[i] = Vec2f(u.bilerp(_pos[i]), v.bilerp(_pos[i]));
    }
}

void Particles::advect(const FaceArray2Xf & u, const FaceArray2Yf & v, float dt)
{
    LOG_OUTPUT("Advecting particles position in the grid velocity field");
    Vec2f mid;
    for(int i = 0; i < _pos.size(); ++i) {
        mid = _pos[i] + 0.5f * dt * Vec2f(u.bilerp(_pos[i]),v.bilerp(_pos[i]));
        _pos[i] = mid + 0.5f * dt * Vec2f(u.bilerp(mid), v.bilerp(mid));
    }
}

void Particles::write(std::ofstream & out) const
{
    int N = numParticles();
    int size = sizeof(Vec2f) * N;
    out.write(reinterpret_cast<const char*>(&N), sizeof(int));
    out.write(reinterpret_cast<const char*>(&_pos[0]), size);
    out.write(reinterpret_cast<const char*>(&_vel[0]), size);
}

void Particles::read(std::ifstream & in)
{
    int N;
    in.read(reinterpret_cast<char *>(&N), sizeof(int));
    _pos.resize(N);
    _vel.resize(N);
    int size = sizeof(Vec2f) * N;
    in.read(reinterpret_cast<char*>(&_pos[0]), size);
    in.read(reinterpret_cast<char*>(&_vel[0]), size);
}

