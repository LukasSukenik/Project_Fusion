#ifndef MEMBRANE_H
#define MEMBRANE_H

#include <sstream>
#include <vector>
#include "xdrfile-1.1.4/include/xdrfile.h"
#include "xdrfile-1.1.4/include/xdrfile_xtc.h"
#include <stdio.h>
#include <cmath>
#include <array>

#include "input.h"

#define HEAD 1
#define TAIL 2
#define RECEPTOR_HEAD 3
#define RECEPTOR_TAIL 4

#define RECEPTOR_BOND_TYPE 4
#define RECEPTOR_ANGLE_TYPE 1

#define TAIL_HEAD_BOND 1
#define TAIL_TAIL_BOND 2
#define HARMONIC_BOND 3

using namespace std;

class Particle {
public:
    Particle() : N(-1) {}
    Particle(int N, int atom_type, double x, double y, double z) : N(N), mol_tag( 0 ), atom_type(atom_type), q(0), x(x), y(y), z(z), nx(0), ny(0), nz(0) {}
    Particle(int N, int mol_tag, int atom_type, double x, double y, double z) : N(N), mol_tag( mol_tag ), atom_type(atom_type), q(0), x(x), y(y), z(z), nx(0), ny(0), nz(0) {}

    int N; // 1 - infi
    int mol_tag; // (N+1)%1000
    int atom_type; // 1 head, 2 tail
    double q,x,y,z; // 0, x, y, z
    int nx,ny,nz; // 0, 0, 0

    bool lipidLoaded = false; // variable for fast loading, not lammps format related

    double dist(Particle& o) {
        return sqrt( (this->x-o.x)*(this->x-o.x) + (this->y-o.y)*(this->y-o.y) + (this->z-o.z)*(this->z-o.z) );
    }
};

//N bond-type atom-1 atom-2
class Bond {
public:
    Bond() : N(-1), type(-1), at1(-1), at2(-1) {}
    int N, type, at1, at2;
};

// N angle-type atom-1 atom-2 atom-3, atom 2 is center
class Angle {
public:
    Angle(){}
    int N, type, at1, at2, at3;
};

class My_string{
public:
    My_string(char str[256] ) {
        strcpy(this->str,str);
    }
    char str[256];
};

class Lipid {
public:
    Lipid(){}

    Lipid(Particle p1, Particle p2, Particle p3) {
        part[0] = p1;
        part[1] = p2;
        part[2] = p3;

        changeN(p1.N, 1);
        bond[0].type = TAIL_HEAD_BOND;
        bond[1].type = TAIL_TAIL_BOND;
        bond[2].type = HARMONIC_BOND;
    }

    Particle part[3];  // consecutive:  0:HEAD, 1:TAIL, 2:TAIL (the further one)
    Bond bond[3];      // 0: between 1-2 TAIL_HEAD_BOND, 1: 2-3 TAIL_TAIL_BOND, 2: 1-3 HARMONIC_BOND

    void changeN(int N, int mol_tag) {
        part[0].N = N;
        part[1].N = N+1;
        part[2].N = N+2;

        part[0].mol_tag = mol_tag;
        part[1].mol_tag = mol_tag;
        part[2].mol_tag = mol_tag;

        bond[0].N = N;
        bond[1].N = N+1;
        bond[2].N = N+2;

        bond[0].at1 = N;
        bond[0].at2 = N+1;

        bond[1].at1 = N+1;
        bond[1].at2 = N+2;

        bond[2].at1 = N;
        bond[2].at2 = N+2;
    }
};


bool sortN(const Particle& i, const Particle& j) {
    return i.N < j.N;
}


bool isTop(int i, vector<Particle>& particles) {
    //
    // search particles, if another particle of type HEAD found under i, then i is top layer
    //
    for(unsigned int j=0; j<particles.size(); ++j) {
        if(particles[j].atom_type == HEAD) {
            if( ( (particles[j].x - particles[i].x)*(particles[j].x - particles[i].x)
                  + (particles[j].y - particles[i].y)*(particles[j].y - particles[i].y)
                  + (particles[j].z - particles[i].z)*(particles[j].z - particles[i].z) ) < 36.0 ) {
                if( (particles[j].z - 4.0 > particles[i].z) ) {
                    return true;
                }
            }
        }
    }
    return false;
}

class ParticleSort{
public:
    bool operator() (Particle& i, Particle& j) {return i.N < j.N;}
};

class BondsSort{
public:
    bool operator() (Bond& i, Bond& j) {return i.at1 < j.at1;}
};

class Membrane {
public:
    Membrane() {}
    Membrane(string in_file) {
        load(in_file);
    }

    bool loaded = false;
    vector<My_string > file_head;
    vector<Bond> bonds;
    vector<Particle> vec;
    vector<Angle> angles;
    vector<Lipid> lipids;

    double last_cm;

    double box[6];


    //
    // Not pretty inheritance
    //
    virtual void trim(double param1, int param2) {}
    virtual void generate(Input& in) {}


    void center(string in_file) {

        cout << "Centering membrane" << endl;

        double move = box[0] + 0.5*(box[1]-box[0]);

        box[0] -= move;
        box[1] -= move;
        box[2] -= move;
        box[3] -= move;

        for(unsigned int i=0; i<vec.size(); ++i) {
            vec[i].x -= move;
            vec[i].y -= move;
        }
    }


    void analyze(string in_file) {

        cout << "Analyzing xtc" << endl;

        const char* name = in_file.c_str();
        XDRFILE* xd;
        char mode;

        int natoms;
        int step;
        float time;
        float prec;
        int status;
        matrix box;

        char fileName[] = "file.xtc";

        status = read_xtc_natoms(fileName, &natoms);

        cout << "Succesful read of natoms: " << natoms << " " << status << " " << exdrOK << endl;

        if (status == exdrOK) {
            XDRFILE* xfp = xdrfile_open(fileName, "r");

            cout << "xdrfile_open" << endl;


            if (xfp != NULL) {
                rvec k[natoms];
                status = read_xtc(xfp, natoms, &step, &time, box, k, &prec);
                last_cm = 10.0;

                //
                // Analyze frame, track types 4 and 2
                //
                for(int j=0; j<1400; ++j) { // step is 25000, have done 37M
                    read_xtc(xfp, natoms, &step, &time, box, k, &prec);
                    cout << step << " " << cm(2, k, box) << " " << cm(4, k, box) << endl;
                    last_cm = cm(4, k, box);
                }


                xdrfile_close(xfp);
            } else {
                cout << "File not opened:" << fileName << endl;
            }

        } else {
            fprintf(stderr, "read_xtc_natoms failure; return code %d", status);
        }
    }

    double cm(int type, rvec* k, matrix box) {
        double cm=0.0;
        int count = 0.0;
        for(unsigned int i=0; i<vec.size(); ++i) {
            if(vec[i].atom_type == type) {
                cm += pbc_z(k[i][2], box[2][2]); // z axis
                ++count;
            }
        }
        return cm/count;
    }

    double pbc_z(double value, double z) {
        if( fabs(value - last_cm) > 50.0 ) {
            return value - z;
        }
        return value;
    }

    void load(string in_file) {
        if(in_file.empty()) {
            in_file = string("data.mem100k");
        }

        cout << "Loading file " << in_file << endl;
        loadFileHeadAndPart(in_file);
        cout << "Loading box parameters" << endl;
        loadBox();

        // THROW AWAY VELOCITIES

        cout << "Parsing bonds parameters" << endl;
        loadBonds(in_file);

        // TEST FOR CONSECCUTIVE INDEXES
        cout << "Sorting lipids by index" << endl;
        std::sort(vec.begin(), vec.end(), sortN);
        cout << "Testing index duplicity" << endl;
        for(int i=0; i<vec.size(); i++) {
            if(i+1 != vec[i].N) {
                cout << "ERROR, missing index, " << i << " != " << vec[i].N << endl;
            }
        }

        cout << "Making lipids whole" << endl;

        loadLipids();
        loaded = true;
    }



    void copyZ(double dist) {
        Lipid temp;
        int size=lipids.size();
        lipids.reserve(2*lipids.size());

        // COPY
        for(int i=0; i<size; ++i) {

            temp = lipids[i];
            temp.part[0].z += dist;
            temp.part[1].z += dist;
            temp.part[2].z += dist; // displace in z dir

            lipids.push_back(temp);
        }

        cerr << lipids[10].part[0].mol_tag << endl;

        convertLipids();

        for(int i=size*3; i<vec.size(); ++i) {
            ++vec[i].mol_tag;
        }
    }


    void printOutSC() {
        std::fstream out ("config.init", std::fstream::out);

        if(out.is_open()) {
            out << box[1] - box[0] << "  " << box[3] - box[2] << "  " << box[5] - box[4] << endl;

            double minX=200,maxX=-200;
            double minY=200,maxY=-200;


            for(unsigned int i=0; i<vec.size(); i++) {
                if(minX > vec[i].x)
                    minX = vec[i].x;

                if(maxX < vec[i].x)
                    maxX = vec[i].x;

                if(minY > vec[i].y)
                    minY = vec[i].y;

                if(maxY < vec[i].y)
                    maxY = vec[i].y;
            }

            for(unsigned int i=0; i<vec.size(); i++) {
                out << vec[i].x - minX << " " << vec[i].y - minY << " " << vec[i].z << "  0 0 0  0 0 0  0 0\n";
            }

            cout << minX << " " << maxX << endl;
        }

        out.close();
    }

    void printOutLammps(Input& in) {

        std::fstream out (in.out_file, std::fstream::out);

        // PRINT OUT
        if(out.is_open())
        {
            out << "LAMMPS data file via gen_membrane\n" << endl;

            out << vec.size() << " atoms\n";
            out << num_types() << " atom types\n";
            out << bonds.size() << " bonds\n";
            out << 3 << " bond types\n";

            if(!angles.empty()) {
                out << angles.size() << " angles\n";
                out << "1 angle types\n";
            }
            out << "\n";

            if(loaded)
            {
                out << box[0] << " " << box[1] << " xlo xhi\n";
                out << box[2] << " " << box[3] << " ylo yhi\n";
                out << box[4] << " " << box[5] << " zlo zhi\n";
            } else
            {
                out << in.box_xm << " " << in.box_xp << " xlo xhi\n";
                out << in.box_ym << " " << in.box_yp << " ylo yhi\n";
                out << in.box_zm << " " << in.box_zp << " zlo zhi\n";
            }

            out << "\nMasses\n\n";
            for(int i=1; i<=num_types(); ++i )
                out << i << " 1.0\n";

            out << "\nAtoms # full\n\n";

            for(unsigned int i=0; i<vec.size(); i++) {
                out << vec[i].N << " " << vec[i].mol_tag << " " << vec[i].atom_type << " ";
                out << vec[i].q << " " << vec[i].x << " " << vec[i].y << " " << vec[i].z << "\n";
            }

            out << "\nBonds\n\n";

            for(unsigned int i=0; i<bonds.size(); i++) {
                out << bonds[i].N << " " << bonds[i].type << " " << bonds[i].at1 << " " << bonds[i].at2 << "\n";
            }

            if(!angles.empty()) {
                out << "\nAngles\n\n";

                for(unsigned int i=0; i<angles.size(); i++) {
                    out << angles[i].N << " " << angles[i].type << " " << angles[i].at1 << " " << angles[i].at2 << " " << angles[i].at3 << "\n";
                }
            }
        }

        out.close();
    }



    /**
     * @brief populate_receptors Populate membrane with receptors -> change head group type
     * @param frac
     */
    void populate_receptors(double frac, int lenght)
    {
        cout << "Populating membrane receptors, fraction: " << frac << ", receptor length: " << lenght << endl;
        bool both = false;
        srand(999);
        int count_tail = 0;
        int count_receptor_head = 0;
        int count_receptor_tail = 0;
        int count_head = 0;
        if(lenght <= 0)
        {
            for(unsigned int i=0; i<vec.size(); ++i)
            {
                if(vec[i].atom_type == HEAD)
                {
                    if( (rand()) < frac*RAND_MAX )
                    {
                        vec[i].atom_type = RECEPTOR_TAIL;
                    }
                }
            }

        }
        else
        {
            cout << "GENERATING BEADS, BONDS, ANGLES OF RECEPTORS" << endl;
            Angle angle;
            Particle part;
            Bond bond;
            const int head_size = 1;
            double receptor_size = 1.0;
            int num_rec = 0;
            int receptor[lenght + 2];
            int receptor_index[lenght + 2];
            double x,y,z;
            bool cont = false;
            for(unsigned int i=0; i<vec.size(); i++)
            {
                if(vec[i].atom_type == HEAD && (isTop(i, vec) || both) && (rand() < frac*RAND_MAX) ) // head bead on correct side
                {
                    cont = false;

                    for(unsigned int j=0; j<bonds.size(); j++)
                    {
                        if(vec[i].N == bonds[j].at1 && bonds[j].type == 1)
                        {
                            receptor[0] = bonds[j].at2;
                            break;
                        }
                        if(vec[i].N == bonds[j].at2 && bonds[j].type == 1)
                        {
                            receptor[0] = bonds[j].at1;
                            break;
                        }

                    }
                    receptor[1] = vec[i].N;
                    receptor_index[1] = i;

                    double upDown = 0.0;
                    if(isTop(i, vec))
                        upDown = 1.0;
                    else upDown = -1.0;

                    // check overlaps
                    for(int q=0; q<lenght; q++)
                    {
                        for(unsigned int j=0; j<vec.size(); j++)
                        {
                            x = vec[i].x - vec[j].x;
                            x*=x;
                            y = vec[i].y - vec[j].y;
                            y*=y;
                            z = (vec[i].z - (1.0+q)*upDown*receptor_size ) - vec[j].z;
                            z *= z;
                            if((x+y+z) < receptor_size*receptor_size)
                            {
                                cont = true;
                                break;
                            }
                        }
                    }
                    if(cont)
                        continue;

                    num_rec++;

                    // Generate receptors
                    for(int j=0; j<lenght; j++)
                    {
                        part.N = vec.size()+1;
                        receptor[j+2] = part.N;
                        receptor_index[j+2] = vec.size();
                        part.mol_tag = 1;
                        if(j >= lenght-head_size) {
                            part.atom_type = RECEPTOR_HEAD;
                        } else {
                            part.atom_type = RECEPTOR_TAIL;
                        }
                        part.q = 0.0;
                        part.x = vec[i].x;
                        part.y = vec[i].y;
                        part.z = vec[i].z - (j+1.0)*upDown*receptor_size;
                        part.nx = 0;
                        part.ny = 0;
                        part.nz = 0;

                        vec.push_back(part);
                    }
                    for(int j=0; j<lenght; j++)
                    {
                        bond.N = bonds.size()+1;
                        bond.type = RECEPTOR_BOND_TYPE;
                        bond.at1 = receptor[j+1];
                        bond.at2 = receptor[j+2];
                        bonds.push_back(bond);
                    }
                    // GENERATE ANGLES For receptors
                    for(int j=0; j<lenght; j++)
                    {
                        angle.N = angles.size()+1;
                        angle.type = RECEPTOR_ANGLE_TYPE;
                        angle.at1 = receptor[j];
                        angle.at2 = receptor[j+1];
                        angle.at3 = receptor[j+2];
                        angles.push_back(angle);
                    }
                }

            }
        }

        for(unsigned int i=0; i<vec.size(); ++i) {
            if(vec[i].atom_type == HEAD) {
                ++count_head;
            }
            if(vec[i].atom_type == RECEPTOR_HEAD) {
                ++count_receptor_head;
            }
            if(vec[i].atom_type == TAIL) {
                ++count_tail;
            }
            if(vec[i].atom_type == RECEPTOR_TAIL) {
                ++count_receptor_tail;
            }
        }
        cout << "Lipid head: " << count_head << ", lipid tail: " << count_tail << ", receptor head: " << count_receptor_head << ", receptor tail: " << count_receptor_tail << endl;
    }

protected:
    void convertLipids() {
        vec.clear();
        bonds.clear();

        cout << "converting lipids -> particles + bonds\nParticles: " << vec.size() << ", Bonds: " << bonds.size() << ", Lipids: " << lipids.size() << endl;

        for(int i=0; i<lipids.size(); i++) {
            if(i%1000 == 0)
                cout << "converting lipid " << i << endl;

            lipids[i].changeN( (3*i) + 1, 1);

            for(int j=0; j<3; j++) {
                vec.push_back(lipids[i].part[j]);
                bonds.push_back(lipids[i].bond[j]);
            }
        }

        cout << "converting done...\nParticles: " << vec.size() << ", Bonds: " << bonds.size() << ", Lipids: " << lipids.size() << endl;
    }


    void loadBox() {

        stringstream str;
        stringstream str2;
        stringstream str3;

        for(unsigned int i=0; i<file_head.size(); i++) {
            //cout << file_head[i].str << endl;
            if(strstr(file_head[i].str, "xlo xhi") != NULL) {
                str << file_head[i].str;
                str >> box[0] >> box[1];
                continue;
            }

            if(strstr(file_head[i].str, "ylo yhi") != NULL) {
                str2 << file_head[i].str;
                str2 >> box[2] >> box[3];
                continue;
            }

            if(strstr(file_head[i].str, "zlo zhi") != NULL) {
                str3 << file_head[i].str;
                str3 >> box[4] >> box[5];
                continue;
            }
        }

        //cout << "BOX: " << box[0] << " " << box[1] << " " << box[2] << " " << box[3] << " " << box[4] << " " << box[5] << endl;
    }

    void search(Lipid& lipid, bool& next, array<int, 3>& at, array<int, 2>& hORt2, array<int, 2>& hORt, array<int, 2>& tORt2) {

        for(int j=0; j<bonds.size(); j++) {
            if( bonds[j].at1 == at[0] && next ) { // this is the bond
                if(bonds[j].type == TAIL_HEAD_BOND) {
                    lipid.bond[0] = bonds[j];
                    hORt[0] = at[0];
                    hORt[1] = bonds[j].at2;
                    at[2] = bonds[j].at2;
                }

                if(bonds[j].type == TAIL_TAIL_BOND) {
                    lipid.bond[1] = bonds[j];
                    tORt2[0] = at[0];
                    tORt2[1] = bonds[j].at2;
                    at[2] = bonds[j].at2;
                }

                if(bonds[j].type == HARMONIC_BOND) {
                    lipid.bond[2] = bonds[j];
                    hORt2[0] = at[0];
                    hORt2[1] = bonds[j].at2;
                    at[2] = bonds[j].at2;
                }
                break;
            }
            if( bonds[j].at1 == at[0] && !next ) { // this is the bond
                next=true;
                if(bonds[j].type == TAIL_HEAD_BOND) {
                    lipid.bond[0] = bonds[j];
                    hORt[0] = at[0];
                    hORt[1] = bonds[j].at2;
                    at[1] = bonds[j].at2;
                }

                if(bonds[j].type == TAIL_TAIL_BOND) {
                    lipid.bond[1] = bonds[j];
                    tORt2[0] = at[0];
                    tORt2[1] = bonds[j].at2;
                    at[1] = bonds[j].at2;
                }

                if(bonds[j].type == HARMONIC_BOND) {
                    lipid.bond[2] = bonds[j];
                    hORt2[0] = at[0];
                    hORt2[1] = bonds[j].at2;
                    at[1] = bonds[j].at2;
                }
            }
        } // NOW We know 2 bonds and 3 atoms, but not their order
    }


    void loadLipids() {
        Lipid lipid;
        ParticleSort partSort;
        BondsSort bondSort;

        std::sort(vec.begin(), vec.end(), partSort); // Sort vec on particle N
        std::sort(bonds.begin(), bonds.end(), bondSort); // Sort bonds on at1

        array<int, 3> at;
        int otherAt;
        bool next;

        array<int, 2> hORt2;
        array<int, 2> hORt;
        array<int, 2> tORt2; // for determining atoms from bonds

        int head,tail,t2;

        for(int i=0; i< vec.size(); i++) {
            if(i % 3000 == 0) // show every 1000 lipids
                cout << "Making lipid number " << i/3 << endl;

            if(vec[i].lipidLoaded) // is particle loded in some lipid?
                continue;

            // generate empty bonds
            lipid.bond[0] = Bond();
            lipid.bond[1] = Bond();
            lipid.bond[2] = Bond();

            at[0] = i+1; // we are indexing from 1 in lammps
            head= -1; tail=-1; t2=-1;

            // Find other particles and bonds, each atom in lipid has 2 bonds
            hORt2[0] = -1; hORt[0] = -1; tORt2[0] = -1;
            hORt2[1] = -1; hORt[1] = -1; tORt2[1] = -1;
            next = false;

            // NOW We know 2 bonds and 3 atoms, but not their order
            search(lipid, next, at, hORt2, hORt, tORt2);

            otherAt = lipid.bond[0].at2;
            if(otherAt == -1)
                otherAt = lipid.bond[1].at2;

            if(hORt[0] == hORt2[0])  head = hORt[0];
            if(tORt2[0] == hORt[0])  tail = hORt[0];
            if(tORt2[0] == hORt2[0]) t2 = hORt2[0];

            if(hORt[1] == hORt2[1])  head = hORt[1];
            if(tORt2[1] == hORt[1])  tail = hORt[1];
            if(tORt2[1] == hORt2[1]) t2 = hORt2[1];



            next = false;
            for(int j=0; j<bonds.size(); j++) {
                if( bonds[j].at1 == otherAt && next ) { // this is the bond
                    if(bonds[j].type == TAIL_HEAD_BOND) {
                        lipid.bond[0] = bonds[j];
                        hORt[0] = otherAt;
                        hORt[1] = bonds[j].at2;
                    }

                    if(bonds[j].type == TAIL_TAIL_BOND) {
                        lipid.bond[1] = bonds[j];
                        tORt2[0] = otherAt;
                        tORt2[1] = bonds[j].at2;
                    }

                    if(bonds[j].type == HARMONIC_BOND) {
                        lipid.bond[2] = bonds[j];
                        hORt2[0] = otherAt;
                        hORt2[1] = bonds[j].at2;
                    }
                    break;
                }
                if( bonds[j].at1 == otherAt && !next ) { // this is the bond
                    next=true;
                    if(bonds[j].type == TAIL_HEAD_BOND) {
                        lipid.bond[0] = bonds[j];
                        hORt[0] = otherAt;
                        hORt[1] = bonds[j].at2;
                    }

                    if(bonds[j].type == TAIL_TAIL_BOND) {
                        lipid.bond[1] = bonds[j];
                        tORt2[0] = otherAt;
                        tORt2[1] = bonds[j].at2;
                    }

                    if(bonds[j].type == HARMONIC_BOND) {
                        lipid.bond[2] = bonds[j];
                        hORt2[0] = otherAt;
                        hORt2[1] = bonds[j].at2;
                    }
                }
            } // NOW We know 3 bonds and 3 atoms

            if(hORt[0] == hORt2[0])  head = hORt[0];
            if(tORt2[0] == hORt[0])  tail = hORt[0];
            if(tORt2[0] == hORt2[0]) t2 = hORt2[0];

            if(hORt[1] == hORt2[1])  head = hORt[1];
            if(tORt2[1] == hORt[1])  tail = hORt[1];
            if(tORt2[1] == hORt2[1]) t2 = hORt2[1];

            for(int x=0; x<3; ++x) {
                if( !(at[x] == head || at[x] == tail || at[x] == t2) ) { // at1 is not used
                    if(head == -1)  head = at[x];
                    if(tail == -1)  tail = at[x];
                    if(t2 == -1) t2 = at[x];
                }
            }

            lipid.part[0] = vec[head-1];
            lipid.part[1] = vec[tail-1];
            lipid.part[2] = vec[t2-1];

            vec[at[0]-1].lipidLoaded = true;
            vec[at[1]-1].lipidLoaded = true;
            vec[at[2]-1].lipidLoaded = true;

            if(tail==-1 || t2==-1 || head==-1) {
                cerr << "ERROR, failure to load lipid, head: " << head << ", tail: " << tail << ", tail2: " << t2 << endl;
                cerr << at[0] << " " << at[1] << " " << at[2] << endl;
                exit(1);
            }

            lipids.push_back(lipid);
        }
    }


    void loadFileHeadAndPart(string filename) {

        char str[256];
        Particle part;
        std::fstream in;
        std::stringstream ss;
        bool first=true;

        in.open(filename, std::fstream::in);

        if(in.is_open()) {
            while(in.good() && strstr(str, "Atoms") == NULL) {
                in.getline(str, 256);
                file_head.push_back(My_string(str));
            }
            // N molecule-tag atom-type q x y z nx ny nz
            while(in.good()) {
                part.N = -1;
                part.nx = 0;
                part.ny = 0;
                part.nz = 0;
                in.getline(str, 256);
                ss.str( str );
                ss >> part.N >> part.mol_tag >> part.atom_type >> part.q >> part.x >> part.y >> part.z >> part.nx >> part.ny >> part.nz;
                ss.flush();
                ss.clear();

                if( part.N == -1 ) {
                    if(vec.empty())
                        continue;
                    else
                        break;
                }
                vec.push_back(part);
            }
            cout << "Loaded " << vec.size() << " particles" << endl;
        }
        in.close();
    }

    void loadBonds(string filename) {

        char str[256];
        Bond bond;
        std::fstream in;

        in.open(filename, std::fstream::in);
        if(in.good()) {
            while(in.good()) {
                in.getline(str, 256);
                if(strstr(str, "Bonds") != NULL)
                    break;
            }

            while(in.good()) {
                in >> bond.N >> bond.type >> bond.at1 >> bond.at2;
                if(bonds.size() == 0 || bond.N != bonds[bonds.size()-1].N) {
                    bonds.push_back(bond);
                }
            }
        }
        in.close();
    }

    /**
     * @brief num_types - return number of types used in sim, NOTE: we assing consecutive integers from 1 as types
     * @return
     */
    int num_types() {
        int num=0;
        for(unsigned int i=0; i<vec.size(); ++i) {
            if(num < vec[i].atom_type)
                num = vec[i].atom_type;
        }
        return num;
    }

    void min_dist(double scale, int start, int stop) {
        double min_head = 1000;
        double min_tail1 = 1000;
        double min_tail2 = 1000;
        for(int i=start; i < stop; i++) {
            for(int j=start; j < stop; j++) {
                if(i != j) {
                    if( lipids[i].part[0].dist(lipids[j].part[0]) < min_head )
                        min_head = lipids[i].part[0].dist(lipids[j].part[0]);

                    if( lipids[i].part[1].dist(lipids[j].part[1]) < min_tail1 )
                        min_tail1 = lipids[i].part[1].dist(lipids[j].part[1]);

                    if( lipids[i].part[2].dist(lipids[j].part[2]) < min_tail2 )
                        min_tail2 = lipids[i].part[2].dist(lipids[j].part[1]);
                }
            }
        }

        // Verify mindist on heads -> 1.06633894589390435
        // TAILS -> 1.122462048309373
        cerr << "Minimal head dist   :  " << scale*min_head << " > " << 1.06633894589390435 << endl;
        cerr << "Minimal tail 1 dist :  " << scale*min_tail1 << " > " << 1.122462048309373 << endl;
        cerr << "Minimal tail 2 dist :  " << scale*min_tail2 << " > " << 1.122462048309373 << endl;
    }

};

#endif // MEMBRANE_H
