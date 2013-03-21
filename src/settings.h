#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "ptr.h"
#include "vec2.h"

class Settings : public SmartPtrInterface<Settings>
{
  public:
    typedef SmartPtr<Settings> Ptr;

    static Ptr create() { return new Settings(); }
    
    // Resolution and size
    int nx;
    int ny;
    float dx;

    // SDF
    float R;
    float r;
    int numPhiSweepIterations;

    // Grid
    Vec2f gravity;
    int numVelSweepIterations;
    
    // PCG
    float tolerance;
    int maxIterations;

  protected:
    Settings() {}
    Settings(const Settings &);
    void operator=(const Settings &);
};

#endif
