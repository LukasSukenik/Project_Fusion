#ifndef FLATMEMBRANE_H
#define FLATMEMBRANE_H

#include "membrane.h"


class FlatMembrane : public Membrane {
public:
    FlatMembrane();
    FlatMembrane(string in_file) : Membrane(in_file) {}

    void trim(double trim, int multiple) {

        int b_size = 0;
        int size;
        int trim_size = lipids.size();
        Bond bond;
        Particle part;


        cout << "Trimming existing membrane by " << trim << endl;
        // TRIM the membrane edges -> avoid lipids over boundary conditions
        for(int i=lipids.size()-1; i>=0; --i) {
            if(lipids[i].part[0].x < box[0]+trim || lipids[i].part[0].x > box[1]-trim || lipids[i].part[0].y < box[2]+trim || lipids[i].part[0].y > box[3]-trim ||
               lipids[i].part[1].x < box[0]+trim || lipids[i].part[1].x > box[1]-trim || lipids[i].part[1].y < box[2]+trim || lipids[i].part[1].y > box[3]-trim ||
               lipids[i].part[2].x < box[0]+trim || lipids[i].part[2].x > box[1]-trim || lipids[i].part[2].y < box[2]+trim || lipids[i].part[2].y > box[3]-trim) {
                lipids.erase(lipids.begin()+i);
            }
        }

        cout << "Trimmed " << trim_size - lipids.size() << endl;

        // change box
        double boxsize[2];
        box[0] += (trim-0.6);
        box[1] -= (trim-0.6);
        box[2] += (trim-0.6);
        box[3] -= (trim-0.6);
        boxsize[0] = box[1] - box[0];
        boxsize[1] = box[3] - box[2];

        // generate new vec and bonds
        convertLipids();

        // COPY membrane + generate new bonds
        size = vec.size();
        for(int x=0; x < multiple; x++) {
            for(int y=0; y < multiple; y++) {
                if(!(x == 0 && y == 0)) {
                    b_size = bonds.size();
                    for(int i=0; i<size; i++) {
                        part = vec[i];
                        part.x += x*boxsize[0];
                        part.y += y*boxsize[1];

                        part.N += b_size;

                        bond = bonds[i];
                        bond.N +=   b_size;
                        bond.at1 += b_size;
                        bond.at2 += b_size;

                        bonds.push_back(bond);
                        vec.push_back(part);
                    }
                }
            }
        }

        // change box size
        for(int i=1; i<multiple; i++) {
            box[1] += boxsize[0];
            box[3] += boxsize[1];
        }
    }
};

#endif // FLATMEMBRANE_H
