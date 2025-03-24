#ifndef ICOSAHEDRON_H
#define ICOSAHEDRON_H

#include "particle.h"
#include "surface.h"
#include <array>

using namespace std;

class IcoVertice : public array<Atom,12>
{
public:
    IcoVertice()
    {
        //
        // Vertices of icosahedron centered on (0, 0, 0), with edge length 1.0.
        //
        (*this)[0] = Atom(0.0,                                   0.0,                               -5/( sqrt(50-10*sqrt(5)) ));
        (*this)[1] = Atom(0.0,                                   0.0,                                5/( sqrt(50-10*sqrt(5)) ));
        (*this)[2] = Atom(-sqrt(2/(5-sqrt(5))),                  0.0,                               -1/( sqrt(10-2*sqrt(5)) ));

        (*this)[3] = Atom(sqrt(2/(5-sqrt(5))),                   0.0,                                1/( sqrt(10-2*sqrt(5)) ));
        (*this)[4] = Atom((1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),   -0.5,                               -1/( sqrt(10-2*sqrt(5)) ));
        (*this)[5] = Atom((1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),    0.5,                               -1/( sqrt(10-2*sqrt(5)) ));

        (*this)[6] = Atom(-(1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),  -0.5,                                1/( sqrt(10-2*sqrt(5)) ));
        (*this)[7] = Atom(-(1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),   0.5,                                1/( sqrt(10-2*sqrt(5)) ) );
        (*this)[8] = Atom(-(-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))), -0.5*sqrt((5+sqrt(5))/(5-sqrt(5))), -1/( sqrt(10-2*sqrt(5)) ) );

        (*this)[9] = Atom(-(-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),  0.5*sqrt((5+sqrt(5))/(5-sqrt(5))), -1/( sqrt(10-2*sqrt(5)) ) );
        (*this)[10] = Atom((-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))), -0.5*sqrt((5+sqrt(5))/(5-sqrt(5))),  1/( sqrt(10-2*sqrt(5)) ));
        (*this)[11] = Atom((-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),  0.5*sqrt((5+sqrt(5))/(5-sqrt(5))),  1/( sqrt(10-2*sqrt(5)) ));
    }
};

template<class Surface>
class Icosahedron : public Particle {
    Surface surface;
public:
    IcoVertice edge;
    vector<int> types;
    int molTag_offset;
    int ligandModulo;
    int offset;

    Icosahedron(string name) : Particle(name), surface(edge, beads) {}

    virtual void generate( Data& data )
    {
        //bparam.push_back();
        molTag_offset = data.getMaxMolTag();
        offset = data.all_beads.size();
        this->types = data.in.types;
        ligandModulo = data.in.num_lig;

        this->icoFrame( data.in.num_of_beads, data.in.scale, data.in.com_pos); // generate icosahedron
        setUnitSize();
    }

protected:
    void setUnitSize()
    {
        double len = edge[0].size();
        for(int i=0; i<beads.size(); ++i) {
            beads[i] *= 1.0/len;
        }
    }

    bool isEdge(int i, int j, int size)
    {
        if(j == 0) return true;
        if(i == size-1) return true;
        if(i == j) return true;
        return false;
    }

    void icoFrame(int size, double scale=1.0, Atom com_pos = Atom(1.0,1.0,1.0))
    {
        Atom push;
        int count=0;
        vector<int> face;

        //
        // Loop over faces -> 20 Faces
        //
        for(int a=0; a<12; a++) { // 12 Vertices
            for(int b=0; b<12; b++) {
                if(edge[a].isNeighbor(edge[b]) ) {
                    for(int c=0; c<12; c++) {
                        if(edge[b].isNeighbor(edge[c]) && edge[a].isNeighbor(edge[c]) ) // for face => triangular area:
                        {
                            ++count;

                            for(int i=0; i<size; i++) { // separate into smaller triangles
                                //
                                // i==0 and i==size-1 on j==0 and j==size-1 -> 3 vertices, each is duplicate 4 times
                                //

                                for(int j=0; j<=i; j++) {
                                    // This is executed as:
                                    //                          .       i sets depth as well as number of repetitions
                                    //                          ..
                                    //                          ...
                                    //                          ....
                                    // i - number of lines, j = number of columns

                                    push = getSubPoint(edge[a], edge[b], edge[c], size, i, j);

                                    if( isVertice(i,j,size) )
                                    {
                                        // Push the point inside into the nanoparticle
                                        // Particle centered at (0.0, 0.0, 0.0) -> simply scale appropriately by given radiuses of vertice and next sphere in line
                                        push *= ( getSubPoint(edge[a], edge[b], edge[c], size, 1, 0).size() ) / push.size();
                                    }

                                    push *= scale;
                                    push += com_pos;

                                    if( isSame(beads, push) )
                                    {
                                        if( !isEdge(i, j, size) )
                                        {
                                            bool exist = false;
                                            for(auto item : face)
                                            {
                                                if(count == item)
                                                {
                                                    exist = true;
                                                }
                                            }
                                            if( !exist )
                                            {
                                                face.push_back(count);
                                            }

                                        }
                                        setType( size,i,j,a,b,c, push, face.size() );
                                        push.N = offset + beads.size() +1;
                                        beads.push_back(push);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void setType(int size, int i, int j, int a, int b, int c, Atom& push, int face)
    {
        if(i%ligandModulo == 0 && j%ligandModulo ==0)
        {
            push.type = types[1]; // typeLig
            if( isEdge(i, j, size))
            {
                push.mol_tag = molTag_offset+2;
            }
            else
            {
                push.mol_tag = molTag_offset+2+face;
            }
        }
        else
        {
            push.type = types[0]; // type nano
            push.mol_tag = molTag_offset+1;
        }
    }

    bool isVertice(int i, int j, int size) {
        if( ( j == 0 && j == i ) || ( j == i && (i == size-1) ) || ( (i == size-1) && ( j == 0 ) ) )
            return true;
        return false;
    }

    Atom getSubPoint(Atom& a, Atom& b, Atom& c, int size, int i, int j) {
        Atom vecAB( b.x - a.x, b.y - a.y, b.z - a.z); // vector from a to b
        Atom vecBC( c.x - b.x, c.y - b.y, c.z - b.z); // vector from b to c
        size-=1;

        vecAB *= (1.0/size);
        vecBC *= (1.0/size);

        return Atom(a.x + i*vecAB.x + j*vecBC.x, a.y + i*vecAB.y + j*vecBC.y, a.z + i*vecAB.z + j*vecBC.z);
    }
};


#endif // ICOSAHEDRON_H
