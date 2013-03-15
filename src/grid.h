#ifndef GRID_H_
#define GRID_H_

#include "ptr.h"

class Grid : public SmartPtrInterface<Grid>
{
  public:
    typedef SmartPtr<Grid> Ptr;

    static Ptr create(int nx, int ny, float dx);

    int nx() const { return _nx; }
    int ny() const { return _ny; }
    float dx() const { return _dx; }
    float invDx() const { return _invDx; }
    
  protected:
    int _nx;
    int _ny;
    float _dx;
    float _invDx;
    
    Grid(int nx, int ny, float dx);
    Grid();
    Grid(const Grid &);
    void operator=(const Grid&);
};

#endif
