#ifndef FLIP2D_H_
#define FLIP2D_H_

#include "ptr.h"
#include "settings.h"
#include "grid.h"
#include "pressure.h"

class FLIP2D : public SmartPtrInterface<FLIP2D>
{
  public:
    typedef SmartPtr<FLIP2D> Ptr;
    
    static FLIP2D::Ptr create(Settings::Ptr s)
    {
        return new FLIP2D(s);
    }

    void step(float dt);

    bool write(const char * filename) const;

    static bool read(const char * filename, Settings::Ptr s, Particles::Ptr p);
    
  protected:
    Settings::Ptr _settings;
    Grid::Ptr _grid;
    Particles::Ptr _particles;
    FluidSDF::Ptr _fluid;
    SolidSDF::Ptr _solid;
    PressureSolver::Ptr _pressureSolver;
    
    FLIP2D(Settings::Ptr s);
    FLIP2D();
    FLIP2D(const FLIP2D &);
    void operator=(const FLIP2D &);
};

#endif
