#ifndef CHAIN_H
#define CHAIN_H

#include "particle.h"



class Chain : public Particle
{
public:
    Chain() : Particle("Dodecahedron chain") {}
    Chain(string name) : Particle(name) {}

    const bool angle_pot = false;
    const int chain_size= 3705;
    const int moltype = 12;
    const double bead_size = 12.0;

    void loadFF( Data& data )
    {
        sigma_size = data.all_sigma_size;
        for(int i=0; i<data.all_sigma_size; ++i)
        {
            for(int j=0; j<data.all_sigma_size; ++j)
            {
                sigma[i][j] = data.all_sigma[i][j];
                sigma_cosatt[i][j] = data.all_sigma_cosatt[i][j];
            }
        }
    }

    void generate( Data& data )
    {
        // Chain parameters
        int res_num = data.in.num_lig;
        int type = data.all_sigma_size;

        // Force field stuff - set sigma
        loadFF(data);
        ++sigma_size;
        sigma[type][type] = bead_size;

        if(res_num > 0) {
            ++sigma_size;
            sigma[type+1][type+1] = bead_size;
            mixing_rules(type+1);
        }
        mixing_rules(type);
        sigma_cosatt[11][11] = true;

        // Bond parameters
        double bond_size = 11.0;
        int offset = data.all_beads.size()+1;
        int bondOffset = data.all_bonds.size()+1;
        double convert = bond_size; // / data.scale;

        cerr << "Chain type:" << type << endl;
        printSigma(type);

        //
        // Starting point
        //
        Atom next;
        for(auto& a : data.all_beads) {
            if(a.type == 5) {
                next = a;
                break;
            }
        }
        // First chain BEAD
        next.N = data.all_beads.size()+1;
        next.type = type+1;
        next.z += 100.0;
        next.mol_tag = moltype;
        beads.push_back(next);

        //
        // Generating Chain
        //
        sigma[type][type] = 6.5; // needs to be below 11
        bool clashExist = true;
        int tick=0;
        for(int i=1; i<chain_size; ++i) {
            clashExist=true;
            while( clashExist ) {
                if(tick % 10000 == 0)
                {
                    cerr << "Try:" << tick << ", chain bead:" << i << endl;
                }
                if(tick > 100000)
                {
                    cerr << "Error, chain not generated" << endl;
                    exit(0);
                }
                clashExist=false;

                // Generate new chain bead
                next.randomUnitSphere();
                next = next*convert + beads.back(); // move atom to last generated atom (+ convert * random unit dist)

                clashExist = ( clash(next, beads) || clash(next, data.all_beads) );

                // Add bead to structure
                if( !clashExist ) {
                   next.mol_tag = moltype;
                   if(i >= chain_size-res_num) {
                       next.type=type+2;
                   } else {
                       next.type=type+1;
                   }
                   next.N = beads.back().N + 1;
                   beads.push_back(next);
                   bonds.push_back( Bond(bonds.size() + bondOffset, bonds.size() + bondOffset, i+offset, i+offset-1, bond_size ) );
                   if(i >= chain_size-res_num && angle_pot) {
                       //                                       current      -1        -2
                       angles.push_back( Angle(angles.size(), 1, i+offset, i+offset-1, i+offset-2 ) );
                   }
                } else {
                    ++tick;
                }
            }
        }

        sigma[type][type] = bead_size;
        cout << beads.size() << endl;
    }

    bool clash(Atom& a, vector<Atom>& others)
    {
        bool clashExist = false;
        for( auto& o : others)
        {
            //cerr << a.dist(o) << " " << sigma[a.type-1][o.type-1] << " " << beads.size() << endl;
            if( a.dist(o) < sigma[a.type-1][o.type-1] ) // type inde in memory from 0, but lammps need types to start at 1
            {
                clashExist = true;
            }
        }
        return clashExist;
    }

    void mixing_rules(int i) {
        for(int j = 0; j< sigma_size; ++j) {
            if(i != j) {
                sigma[i][j] = 0.5*(sigma[i][i] + sigma[j][j]);
                sigma[j][i] = 0.5*(sigma[i][i] + sigma[j][j]);
            }
        }

    }
};

#endif // CHAIN_H
