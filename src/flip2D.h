#ifndef FLIP2D_H_
#define FLIP2D_H_

#include "ptr.h"
#include "grid.h"

class FLIP2D : public SmartPtrInterface<FLIP2D>
{
  public:
    typedef SmartPtr<FLIP2D> Ptr;
    
    static FLIP2D::Ptr create(int nx, int ny, float dx);

    void step(float dt);

    bool write(const char * filename);
    
  protected:
    Grid::Ptr _grid;
    
    FLIP2D(int nx, int ny, float dx);
    FLIP2D();
    FLIP2D(const FLIP2D &);
    void operator=(const FLIP2D &);
};

#endif
