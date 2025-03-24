#ifndef PENTAMER_H
#define PENTAMER_H

#include "icosahedron.h"
#include "surface.h"

template<class Surface>
class Pentamer : public Icosahedron<Surface> {
public:
    using Icosahedron<Surface>::edge;
    using Icosahedron<Surface>::beads;
    using Icosahedron<Surface>::bonds;

    Pentamer(string name) : Icosahedron<Surface>(name)  {}

    virtual void generate( Data& data ) override
    {
        double len = edge[0].size();
        double limit = 1.0 / (data.in.num_of_beads-1);
        this->icoFrame(data.in.num_of_beads);
        for(int i=0; i<beads.size(); ++i) {
            beads[i] *= 1.0/len;
        }
        int start = beads.size();
        this->icoFrame(data.in.num_of_beads);
        for(int i=start; i<beads.size(); ++i) {
            beads[i] *= 1.143/len;
        }

        for(int i=0; i<start; ++i) {
            for(int a=0; a<12; ++a) {
                if( beads[i].type == a ) {
                    for(int j=i+1; j<start; ++j) {
                        if(beads[j].type == a && (beads[j]-beads[i]).size() < limit+0.01) {
                            bonds.push_back( Bond(bonds.size(), bonds.size(), i+1, j+1, (beads[j]-beads[i]).size() ) );
                        }
                    }
                }
            }
        }
        for(int i=start; i<beads.size(); ++i) {
            for(int a=0; a<12; ++a) {
                if( beads[i].type == a ) {
                    for(int j=i+1; j<beads.size(); ++j) {
                        if(beads[j].type == a && (beads[j]-beads[i]).size() < limit*1.143+0.01) {
                            bonds.push_back( Bond(bonds.size(), bonds.size(), i+1, j+1, (beads[j]-beads[i]).size() ) );
                        }
                    }
                }
            }
        }
        for(int i=0; i<start; ++i) {
            for(int a=0; a<12; ++a) {
                if( beads[i].type == a ) {
                    for(int j=start; j<beads.size(); ++j) {
                        if(beads[j].type == a && (beads[j]-beads[i]).size() < 0.125) {
                            bonds.push_back( Bond(bonds.size(), bonds.size(), i+1, j+1, (beads[j]-beads[i]).size() ) );
                        }
                    }
                }
            }
        }

        cerr << "Bonds: " << bonds.size() << endl;

        setTypes(0, start, 1.0, 12, 13);
        setTypes(start, beads.size(), 1.0/len, 12, 13);

        for(int i=0; i<beads.size(); ++i) {
            if( beads[i].type < 12 )
                beads[i].type = 1;
            if( beads[i].type == 12 )
                beads[i].type = 2;
            if( beads[i].type == 13 )
                beads[i].type = 3;
        }
    }

private:
    void setTypes(int start, int stop, double scale, int type1=2, int type2=3)
    {
        // Generate Interaction A, B types
        for(int i=start; i<stop; ++i) {
            // identify edge
            for(int a=0; a<12; ++a) {
                if(beads[i].type == a && (beads[i] - edge[a]).size() > 0.45*scale ) { // closest to edge[0]
                    //
                    // Identified pentamer edge
                    //
                    for(int b=0; b<12; ++b) {
                        if( edge[a].isNeighbor(edge[b]) ) {
                            for(int c=0; c<12; ++c) {
                                if( edge[a].isNeighbor(edge[c]) && edge[b].isNeighbor(edge[c]) ) {
                                    //
                                    // a,b,c vertices of triangle
                                    //
                                    Atom center(edge[a].x + edge[b].x + edge[c].x, edge[a].y + edge[b].y + edge[c].y, edge[a].z + edge[b].z + edge[c].z);
                                    Atom ab = edge[a] + edge[b];
                                    Atom ac = edge[a] + edge[c];
                                    Atom bc = edge[b] + edge[c];
                                    center *= scale*1.0/3.0;
                                    ab *= scale*0.5;
                                    ac *= scale*0.5;
                                    bc *= scale*0.5;

                                    if(( (  (edge[a] - center).cross( (edge[b] - center) )  ) + center ).size() < 1.0) {

                                        double lenCenAB = (center-ab).size();

                                        if( (beads[i] - center).size() < lenCenAB && (beads[i] - ab).size() < lenCenAB )
                                            beads[i].type = type1;
                                        if( (beads[i] - center).size() < lenCenAB && (beads[i] - ac).size() < lenCenAB )
                                            beads[i].type = type2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif // PENTAMER_H
