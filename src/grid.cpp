#include "grid.h"

Grid::Grid(int nx, int ny, float dx) : _nx(nx), _ny(ny), _dx(dx)
{
    _invDx = 1.0f / _dx;
}

Grid::Ptr Grid::create(int nx, int ny, float dx)
{
    return new Grid(nx,ny,dx);
}
