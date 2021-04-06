#ifndef SLAB_H
#define SLAB_H

//#include "sphere.h"
#include "particle.h"
#include "cow.h"

class Slab1D : public Particle
{
public:
    Slab1D() : Particle("Slab") {}
    Slab1D(string name) : Particle(name) {}

    void generate( Data& data )
    {
        int row = sqrt(data.in.num_of_beads);
        int bead_size = data.in.scale;
        double factor = 1.0; // length were we generate beads


        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                for(int k=0; k<1; ++k) {
                     beads.push_back(Atom( factor*i, factor*k, factor*j, 0));
                }
            }
        }

        sigma_size = 1;
        sigma[0][0] = bead_size*1.0;

        add(data);
    }
};

class Slab : public Particle
{
public:
    Slab() : Particle("Slab") {}
    Slab(string name) : Particle(name) {}

    void generate( Data& data )
    {
        double factor = data.in.scale;

        sigma[0][0] = 5*factor;
        sigma[0][1] = 3*factor;
        sigma[1][1] = factor;
        sigma_size = 2;
        factor /= 4.0;

        vector<Atom> bVec;
        int size = data.in.num_of_beads;
        for(int a=0; a<size; ++a) {
            for(int b=0; b<size; ++b) {
                for(int c=0; c<size; ++c) {
                    if((a-size/2.0)*(a-size/2.0) + (b-size/2.0)*(b-size/2.0) + (c-size/2.0)*(c-size/2.0) < size*size*0.25)
                        bVec.push_back(Atom(a*factor*2.0 - size/2.0*factor, b*factor*2.0 - size/2.0*factor, c*factor*2.0 - size/2.0*factor, 1, 1));
                }
            }
        }

        Atom move = Atom( 14, 14, 14, 1, 1);
        for(Atom& item : bVec)
            item += move;

        beads.insert(this->beads.end(), bVec.begin(), bVec.end());

        add(data);
    }
private:
    void gen_bonds(int row) {
        int actual;
        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                actual = i*row + j;

                if(j > 0)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i)*row + j-1, actual, beads[(i)*row + j-1].dist(beads[actual]) ) );
                if(i > 0)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i-1)*row + j, actual, beads[(i-1)*row + j].dist(beads[actual]) ) );
                if(j < row-1)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i)*row + j+1, actual, beads[(i)*row + j+1].dist(beads[actual]) ) );
                if(i < row-1)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i+1)*row + j, actual, beads[(i+1)*row + j].dist(beads[actual]) ) );
            }
        }

        // Correct for offset -> lammps starts at 1 not 0
        for(int i=0; i<bonds.size(); ++i) {
            ++bonds[i].at1;
            ++bonds[i].at2;
        }
    }
};


class NettedSlab : public Particle
{
public:
    NettedSlab() : Particle("Slab") {}
    NettedSlab(string name) : Particle(name) {}

    void generate( Data& data )
    {
        int row = sqrt(data.in.num_of_beads);
        int bead_size = data.in.scale;
        double factor = 1.0; // length were we generate beads


        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                for(int k=0; k<1; ++k) {
                    if( j==0 || i == 0 || i == row-1 || j == row-1 )
                        beads.push_back(Atom( factor*i, factor*k, factor*j, 1));
                    else
                        beads.push_back(Atom( factor*i, factor*k, factor*j, 0));
                }
            }
        }

        gen_bonds(row);

        cerr << "Bonds: " << bonds.size() << endl;
        Sphere sp;
        vector<Atom> bVec;
        factor *= 2.5/2.0;



        for(int i=0; i<data.in.num_lig; ++i) {

            Atom next;
            bool clash=true;
            while(clash) {
                double x = ((0.20+ran()/10*6)*data.in.boxp.x - data.in.com_pos.x) / data.in.scale;
                double y = ( (0.20+ran()/4) *data.in.boxp.y - data.in.com_pos.y) / data.in.scale;
                double z = ((0.20+ran()/10*6)*data.in.boxp.z - data.in.com_pos.z) / data.in.scale;

                next = Atom(x,y,z,2);

                /*sp.beads.clear();
                sp.fibonacci_sphere( 300, 2, i+1);
                sp.rescale(data.c);
                sp.move(next);*/

                bVec.clear();
                for(int a=0; a<11; ++a) {
                    for(int b=0; b<11; ++b) {
                        for(int c=0; c<11; ++c) {
                            if((a-5.5)*(a-5.5) + (b-5.5)*(b-5.5) + (c-5.5)*(c-5.5) < 11.0*11.0*0.25)
                                bVec.push_back(Atom(x+a*factor*2.0 - 11.0*factor, y+b*factor*2.0 - 11.0*factor, z+c*factor*2.0 - 11.0*factor, 2, i+1));
                        }
                    }
                }

                clash=false;

                for(int j=0; j<beads.size(); ++j) {
                    if(beads[j].dist(next) < data.in.c ) {
                        clash = true;
                    }
                }

                if(!clash) {
                    //beads.insert(this->beads.end(), sp.beads.begin(), sp.beads.end());
                    beads.insert(this->beads.end(), bVec.begin(), bVec.end());
                }
            }
        }

        sigma[0][0] = bead_size*1.0;
        sigma[1][1] = bead_size*1.0;
        sigma[2][2] = bead_size*2.5;
        sigma_size = 3;

        mixing_rules();
    }

private:
    void mixing_rules() {
        for(int i = 0; i< sigma_size; ++i) {
            for(int j = 0; j< sigma_size; ++j) {
                if(i != j)
                    sigma[i][j] = 0.5*(sigma[i][i] + sigma[j][j]);
            }
        }
    }

    void gen_bonds(int row) {
        int actual;
        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                actual = i*row + j;
                if(j > 0)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i)*row + j-1, actual, beads[(i)*row + j-1].dist(beads[actual]) ) );
                if(i > 0)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i-1)*row + j, actual, beads[(i-1)*row + j].dist(beads[actual]) ) );
                if(j < row-1)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i)*row + j+1, actual, beads[(i)*row + j+1].dist(beads[actual]) ) );
                if(i < row-1)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i+1)*row + j, actual, beads[(i+1)*row + j].dist(beads[actual]) ) );
            }
        }

        // Correct for offset -> lammps starts at 1 not 0
        for(int i=0; i<bonds.size(); ++i) {
            ++bonds[i].at1;
            ++bonds[i].at2;
        }
    }
};

#endif // SLAB_H
