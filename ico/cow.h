#ifndef COW_H
#define COW_H

#include <fstream>
#include <sstream>
#include "atom.h"
#include "particle.h"

class Cow : public Particle
{
public:
    Cow() : Particle("Cow")  {}
    Cow(string name) : Particle(name) {}

    void generate( Data& data  )
    {
        int type=2;
        int mol_tag=2;
        std::ifstream fs( "cow3.xyz" );
        int size=0;
        char line[256];
        string temp;
        double x,y,z;

        fs >> size;
        fs.getline(line, 256);
        for(int i=0; i<size; ++i) {
            fs >> temp >> x >> y >> z;
            beads.push_back( Atom(x, y, z, type, mol_tag) );
        }

        fs.close();
    }
};

#endif // COW_H
