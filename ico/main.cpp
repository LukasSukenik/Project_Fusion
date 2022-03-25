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
#include "xtcanalysis.h"
#include "thermo_analsis.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <cstdlib>
#include <algorithm>

#include <vector>

// Sofia first commit 2

using namespace std;

void helpMessage();

int main(int argc, char* argv[]) // // $num of beads per edge, box dimensions X(same as Y) $beg $end, $position in Z, $offset
{
    Thermo_analysis t;
    /*t.timestep_analyze();
    return 0;*/

    if( argc == 1 || strcmp(argv[1], "-h") == 0 ) {
        cout << "No input file specified" << endl;
        helpMessage();
        exit(1);
    }

    Data data;
    vector<Particle*> nano;
    nano.push_back( new Icosahedron<Surface>("Icosahedron") );                      // 1
    nano.push_back( new Sphere() );                                                 // 2
    nano.push_back( new TennisBall() );                                             // 3
    nano.push_back( new OblateSpheroid() );                                         // 4
    nano.push_back( new SpherePatch() );                                            // 5
    nano.push_back( new Pentamer<PentaSurface>("Icosahedron + PentaSurface" ) );    // 6
    nano.push_back( new Dodecahedron("Dodecahedron") );                             // 7
    nano.push_back( new Chain() );                                                  // 8
    nano.push_back( new Slab() );                                                   // 9
    nano.push_back( new NettedSlab());                                              // 10
    nano.push_back( new SphereJanus());                                             // 11
    nano.push_back( new Cow());                                                     // 12

    XTCAnalysis xtc_analyzator;

    for(int i=1; i<argc; ++i)
    {
        data.loadInput(argv[i]);

        cerr << data.in.toString() << endl;

        if( data.isDefined() )
        {
            data.load(data.in.infile);      // Load Data from file "data.in.infile"
            if(data.in.isScale())
                data.rescale(data.in.scale);    // Rescale atom positions by data.in.scale - usually 1
            if(data.in.isCOM_pos())
                data.move(data.in.com_pos);     // Move by vector defined in input file
            if(data.in.is_mol_tag())
                data.mol_tag(data.in.mol_tag);  // Change mol_tag of all particles to one set by input
            if(data.in.is_mtag_12())
                data.align(data.in.mtag_1, data.in.mtag_2); // align mol_tag particles in z axis and XY plane
            if(data.in.is_mtag_12z())
                data.align_z(data.in.mtag_1z, data.in.mtag_2z); // align mol_tag particles in z axis and XY plane

            data.offset(data.all_beads.size());

            if( data.in.is_fit() )
                data.fit(data.in.fit_x, data.in.fit_y,data.in.fit_z);
            if( data.in.fit_lipo )
                data.fit_lipo();
            if( data.in.center )
                data.center();

            if( !data.in.analize_infile.empty() )
                xtc_analyzator.analyze_histogram( data.in.analize_infile, data );
            if( data.in.timestep )
                t.timestep_analyze();

            data.add();
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
            nano[data.in.nano_type-1]->rotate( data.in.rotation, data );
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
        double COM=0.0; double min=9999; double max=-9999;
        for(Atom& item : data.all_beads)
        {
            if(item.z > max)
                max = item.z;
            if(item.z < min)
                min = item.z;
            COM += item.z;
        }

        ofstream myfile;
        myfile.open ("wall_properties");
        myfile << COM/data.all_beads.size() << " " << min << " " << max << endl;
        myfile.close();

        data.printLammps();
        cerr << "Lammps out" << endl;
    }
    if( data.in.out_type == 2)
    {
        data.printXYZ();
        cerr << "XYZ out" << endl;
    }

    for(int i=0; i<nano.size(); ++i)
        delete nano[i];

    //cerr << data.toString() << endl;

    return 0;
}

void helpMessage() {
    cout << "./ico filename_1 filename_2" << endl;

    cout << "Keywords:" << endl;
    cout << "Ouput/Input category:" << endl;
    cout << "Output_type: 1 = Lammps" << endl;
    cout << "             0 = XYZ" << endl;
    cout << "Lammps_offset: integer" << endl;
    cout << " - offset the generated structure for manual insertion into another lammps structure file" << endl;
    cout << "Load_file: filename" << endl;

    cout << "\nGenerated structure properties category:" << endl;
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
    cout << "Num_of_beads: integer" << endl;
    cout << " - Particle_type: 1,6,7 = number of beads per edge" << endl;
    cout << " - other Particle_type = number of beads for entire nanoparticle/structure" << endl;
    cout << "Scale: floating_point_number - nanoparticle radius" << endl;
    cout << "c: float " << endl;
    cout << " - Particle_type: 4 = oblate spheroid < 1, prolate spheroid > 1, 1.0 - ERROR, not defined" << endl;
    cout << " - Particle_type: 3 = width of patch (0.0 to 2.0)" << endl;
    cout << "Number_of_ligands: integer" << endl;
    cout << "Mol_tag: integer" << endl;
    cout << " - change mol_tag of generated/loaded structure " << endl;
    cout << "Atom_type: integer" << endl;
    cout << " - Atom_type of generated structure, if structure has more atom_types they are incremented from provided value " << endl;
    cout << "Janus: float float float" << endl;

    cout << "\nPosition/Box properties category:" << endl;
    cout << "Box: float float float float float float" << endl;
    cout << "position_shift: float float float" << endl;
    cout << "Center = centers particles at 0.0" << endl;
    cout << "Align: mol_tag_1 mol_tag_2 integer" << endl;
    cout << " - align mol_tag_1 with x-axis" << endl;
    cout << " - center mol_tag_2 around z axis" << endl;
    cout << "Impact_vector: float float float" << endl;
    cout << "Fit" << endl;
    cout << " - positions the loaded/generated structure next to previosly generated/loadedd structure" << endl;
    cout << " - used for ideal collision position of two liposomes and a nanoparticle" << endl;

    cout << "\nForce-Field category:" << endl;
    cout << "Beads_lj/cut:" << endl;

    cout << "\nSeed: intege = radom generator" << endl;
}



