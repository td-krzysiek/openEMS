// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPAQUE_SURFACE_UNION_H
#define OPAQUE_SURFACE_UNION_H
#include <cstddef>
#include <stdexcept>
#include <vector>

// Occupancy of open primary cells. Shared primitive faces are not boundaries
// of the physical conductor union. Coordinates are supplied as mesh indices.
class OpaqueSurfaceUnion
{
public:
    OpaqueSurfaceUnion(unsigned int nx, unsigned int ny, unsigned int nz)
    {
        size[0]=nx; size[1]=ny; size[2]=nz;
        cells.assign(static_cast<size_t>(nx)*ny*nz,-1);
    }
    int Get(int x,int y,int z) const
    {
        if (x<0 || y<0 || z<0 || static_cast<unsigned int>(x)>=size[0] ||
            static_cast<unsigned int>(y)>=size[1] || static_cast<unsigned int>(z)>=size[2]) return -1;
        return cells[(static_cast<size_t>(x)*size[1]+y)*size[2]+z];
    }
    void Add(const unsigned int* lo,const unsigned int* hi,unsigned int fit)
    {
        for (int n=0;n<3;++n)
            if (lo[n]>=hi[n] || hi[n]>size[n]) throw std::runtime_error("Invalid opaque union box bounds");
        for (unsigned int x=lo[0];x<hi[0];++x)
        for (unsigned int y=lo[1];y<hi[1];++y)
        for (unsigned int z=lo[2];z<hi[2];++z) {
            int& value=cells[(static_cast<size_t>(x)*size[1]+y)*size[2]+z];
            if (value>=0 && value!=static_cast<int>(fit)) throw std::runtime_error("Overlapping opaque fits differ");
            value=fit;
        }
    }
    bool InteriorH(unsigned int n,const unsigned int* pos) const
    {
        int low[3]={static_cast<int>(pos[0]),static_cast<int>(pos[1]),static_cast<int>(pos[2])};
        int high[3]={low[0],low[1],low[2]}; --low[n];
        return Get(low[0],low[1],low[2])>=0 && Get(high[0],high[1],high[2])>=0;
    }
    void Quadrants(unsigned int n,const unsigned int* pos,int fit[2][2]) const
    {
        const unsigned int p=(n+1)%3,q=(n+2)%3;
        for (int a=0;a<2;++a) for (int b=0;b<2;++b) {
            int cell[3]={static_cast<int>(pos[0]),static_cast<int>(pos[1]),static_cast<int>(pos[2])};
            cell[p]+=a-1; cell[q]+=b-1;
            fit[a][b]=Get(cell[0],cell[1],cell[2]);
        }
    }
private:
    unsigned int size[3];
    std::vector<int> cells;
};
#endif
