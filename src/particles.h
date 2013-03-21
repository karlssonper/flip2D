#ifndef PARTICLES_H_
#define PARTICLES_H_

#include "ptr.h"
#include "vec2.h"
#include "array.h"

class Particles : public SmartPtrInterface<Particles>
{
  public:
    typedef SmartPtr<Particles> Ptr;

    static Ptr create() { return new Particles(); }

    void addParticle(const Vec2f & pos, Vec2f vel = Vec2f());
    void addParticles(const std::vector<Vec2f> & pos,
                      const std::vector<Vec2f> & vel);
    
    int numParticles() const { return _pos.size(); }
    const Vec2f & pos(int particleIdx) { return _pos[particleIdx]; }
    const Vec2f & vel(int particleIdx) { return _vel[particleIdx]; }

    void updateVelocities(const Array2f & u, const Array2f & v);
    
    void advect(float dt);

  protected:
    std::vector<Vec2f> _pos;
    std::vector<Vec2f> _vel;
    
    Particles();
    Particles(const Particles &);
    void operator=(const Particles &);
};

#endif
