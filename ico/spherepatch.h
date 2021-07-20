#ifndef SPHEREPATCH_H
#define SPHEREPATCH_H

#include "sphere.h"


class SpherePatch : public Sphere
{
public:
    SpherePatch() : Sphere("SpherePatch") {}
    SpherePatch(string name) : Sphere(name) {}

    void generate( Data& data )
    {
        vector<Atom> ligand;

        // generate nano
        int nano_start = beads.size();
        fibonacci_sphere(beads, data.in.num_of_beads, typeNano);
        int nano_end = beads.size();

        // generate temporary ideal ligand placement
        fibonacci_sphere(ligand, data.in.num_lig, typeTemp); // second fib. sphere

        // Erase part of ligand placement
        int num_lig2 = createPatch(nano_end, data.in.c, typeTemp);

        // change type nano to type lig based on placement of type temp
        Atom patch = Atom(1,1,1,typeLig);
        gen_ligands(data, ligand, patch, typeNano); // find closest spheres on first fib. sphere and change their type

        beads.erase(beads.begin()+nano_end, beads.end()); // erase second fib sphere

        for(int i=nano_start; i<nano_end; ++i) {
            if(beads[i].type == typeLig)
                gen_polymer(6,beads[i], data.in.scale);
        }
    }


    int createPatch(int nano_end, double c, int type) {
        for(int i = beads.size()-1; i >= nano_end; --i ) {
            if( beads[i].type == type && !isPatchZ(i,c) ) {
                beads.erase(beads.begin()+i);
            }
        }

        return beads.size() - nano_end;
    }

    bool isPatchX(int i, double c) {
        if( fabs(beads[i].x) > c*0.5 )
            return false;
        return true;
    }

    bool isPatchY(int i, double c) {
        if( fabs(beads[i].y) > c*0.5 )
            return false;
        return true;
    }

    bool isPatchZ(int i, double c) {
        if( fabs(beads[i].z) > c*0.5 )
            return false;
        return true;
    }

    void gen_polymer(int size, Atom start, double scale) {
        Atom unit = start;
        unit.normalise();
        unit *= 1.0/scale;

        for(int i=0; i<size; ++i) {
            start += unit;
            beads.push_back(start);
        }
    }
};

#endif // SPHEREPATCH_H
