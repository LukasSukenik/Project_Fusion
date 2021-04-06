#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <cstdlib>
#include <algorithm>
#include <random>

#include "rng.h"
#include "icosahedron.h"
#include "oblatespheroid.h"
#include "sphere.h"
#include "tennisball.h"
#include "spherepatch.h"
#include "pentamer.h"
#include "dodecahedron.h"
#include "chain.h"
#include "slab.h"

void helpMessage();



int main(int argc, char* argv[]) // // $num of beads per edge, box dimensions X(same as Y) $beg $end, $position in Z, $offset
{
    if( argc == 1 ) {
        cout << "No input file specified, number of parameters: " << argc << endl;
        helpMessage();
        exit(1);
    }

    Data data;
    vector<Particle*> nano;
    nano.push_back( new Icosahedron<Surface>("Icosahedron") );  // 1
    nano.push_back( new Sphere() );                                                     // 2
    nano.push_back( new TennisBall() );                                                 // 3
    nano.push_back( new OblateSpheroid() );                                             // 4
    nano.push_back( new SpherePatch() );                                                // 5
    nano.push_back( new Pentamer<PentaSurface>("Icosahedron + PentaSurface" ) );      // 6
    nano.push_back( new Dodecahedron("Dodecahedron") );  // 7
    nano.push_back( new Chain() );                                                      // 8
    nano.push_back( new Slab() );                                                       // 9
    nano.push_back( new NettedSlab());                                                   // 10
    nano.push_back( new SphereJanus());                                                   // 11
    nano.push_back( new Cow());                                                   // 12

    for(int i=1; i<argc; ++i)
    {
        data.loadInput(argv[i]);

        cerr << data.in.toString() << endl;

        if( !data.in.infile.empty() )
        {
            data.load(data.in.infile);
            data.rescale(data.in.scale);
            data.move(data.in.com_pos);
            if( data.in.mol_tag != -1)
            {
                data.mol_tag(data.in.mol_tag);
            }
            if(data.in.nanoup)
            {
                data.rotate();
            }
            if(data.in.z_stack)
            {
                data.z_stack();
            }
            data.offset(data.all_beads.size());
            data.add();


            if(data.in.center)
            {
                data.center();
            }
        }
        else
        {
            if(data.in.nano_type < 1 || data.in.nano_type > nano.size()) {
                cerr << "Bad nanoparticle type selected" << endl;
                exit(2);
            }

            cerr << "Generating: " << nano[data.in.nano_type-1]->name << endl;
            nano[data.in.nano_type-1]->generate( data );
            nano[data.in.nano_type-1]->rescale( data.in.scale );
            nano[data.in.nano_type-1]->move( data.in.com_pos );
            nano[data.in.nano_type-1]->add(data);
        }
    }

    vector<string> coeff;
    vector<double> dist = data.createBondGroups(coeff);
    data.roundBondDist(dist);
    //data.removeDuplicateBond();

    if( data.in.out_type == 1)
    {
        data.printLammps();
        cerr << "Lammps out" << endl;
    }
    else
    {
        data.printXYZ();
        cerr << "XYZ out" << endl;
    }

    for(int i=0; i<nano.size(); ++i)
        delete nano[i];

    cerr << data.toString() << endl;

    return 0;
}

void helpMessage() {
    cout << "Run binary with filenames that contain the following parameters:" << endl;

    cout << "Particle_type: 1 = Icosahedron" << endl;
    cout << "               2 = Sphere" << endl;
    cout << "               3 = Tennis ball" << endl;
    cout << "               4 = Oblate spheroid" << endl;
    cout << "               5 = Sphere patch" << endl;
    cout << "               6 = Icosahedron pentamer" << endl;
    cout << "               7 = Dodecahedron pentamer, param c -> scale of second dodecahedron" << endl;
    cout << "               8 = chain of particles, num_lig = type 2 chain beads" << endl;
    cout << "               9 = slab of particles" << endl;
    cout << "               10 = netted slab of particles" << endl;
    cout << "               11 = Janus sphere" << endl;
    cout << "               12 = xyz file - name cow2.xyz curently" << endl;

    cout << "Output_type: 1 = Lammps" << endl;
    cout << "             other = XYZ" << endl;

    cout << "Num_of_beads: XXX - for 1,6,7 nano type - number of beads per edge" << endl;
    cout << "                 - others - number of beads for entire nanoparticle/structure" << endl;

    cout << "Scale: XXX - nanoparticle radius" << endl;
    cout << "Lammps_offset: 0 - if previous lammps structure exist" << endl;
    cout << "c: XXX - (oblate spheroid param || tennis ball param || width of patch (0.0 to 2.0)) " << endl;
    cout << "Number_of_ligands: XXX" << endl;
    cout << "box: XXX XXX XXX" << endl;
    cout << "position_shift: XXX XXX XXX" << endl;

    cout << "Load_file:" << endl;
    cout << "Center:" << endl;
    cout << "Seed:" << endl;
    cout << "Beads_lj/cut:" << endl;

    cout << "Seed: XXX - for random generator" << endl;
}



