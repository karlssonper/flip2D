#include "flip2D.h"

FLIP2D::FLIP2D(int nx, int ny, float dx)
{
    _grid = Grid::create(nx,ny,dx);
}

FLIP2D::Ptr FLIP2D::create(int nx, int ny, float dx)
{
    return new FLIP2D(nx,ny,dx);
}

void FLIP2D::step(float dt)
{
    
}

bool FLIP2D::write(const char * filename)
{

}

