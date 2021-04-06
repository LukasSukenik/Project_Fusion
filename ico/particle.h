#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "xdrfile-1.1.4/include/xdrfile.h"
#include "xdrfile-1.1.4/include/xdrfile_xtc.h"

#include<stdio.h>

#include "data.h"
#include "atom.h"


using namespace std;

/**
 * @brief The Cuboid class - Represents standart simulation box - orthogonal with PBC
 */
class Cuboid {
public:
    Atom box; // starts at 0

    Cuboid() {}
    Cuboid(Atom box) : box(box) {}

    inline double volume() {
        return box.x*box.y*box.z;
    }

    Atom usePBC(Atom& pos_orig, double scale) const {
        Atom pos = pos_orig;

        while (pos.x < 0.0) {
            pos.x += box.x/scale;
        }
        while (pos.x > box.x/scale) {
            pos.x -= box.x/scale;
        }
        while (pos.y < 0.0) {
            pos.y += box.y/scale;
        }
        while (pos.y > box.y/scale) {
            pos.y -= box.y/scale;
        }
        while (pos.z < 0.0) {
            pos.z += box.z/scale;
        }
        while (pos.z > box.z/scale) {
            pos.z -= box.z/scale;
        }
        return pos;
    }
};


/**
 * @brief The Particle class - Abstract class of Particle - a assembly made from several beads
 */
class Particle{
public:
    Cuboid pbc;
    string name;
    const double degToRad = 0.0174532925;

    /// Local data for structure, need to save to Data class for output
    vector<Atom> beads;
    vector< Bond > bonds;
    vector<Angle> angles;
    vector<BeadParam> bparam;

    int typeNano = 1;
    int typeLig = 2;
    const int typeTemp = 99;

    const int orientX = 101;
    const int orientY = 102;
    const int orientZ = 103;

    /// Force-field Data
    int sigma_size = 0;
    array<array<double, 100>, 100> sigma;   
    array< array<bool, 100>, 100> sigma_cosatt;

    //
    // Class stuff
    //
    Particle(string name) : name(name) {}
    virtual ~Particle() {}


    /**
     * @brief generate - Generate particle
     * - Particle need to be generated in scale and in position (data.in.scale and data.in.com_pos)
     * - Position overlap check to already existing particles in data.all_beads
     *
     * @param data
     */
    virtual void generate( Data& data )=0;

    void printSigma()
    {
        cerr << "printSigma:" << endl;
        for(unsigned int i=0; i<sigma_size; ++i)
        {
            cerr << "[" << i+1 << "][" << i+1 << "] = "  << sigma[i][i] << endl;
        }
    }

    void printSigma(int j)
    {
        cerr << "printSigma type " << j << ":" << endl;
        for(unsigned int i=0; i<sigma_size; ++i)
        {
            cerr << "[" << j+1 << "][" << i+1 << "] = "  << sigma[j][i] << " == " << sigma[i][j] << endl;
        }
    }

    void move(Atom move)
    {
        for(Atom& item : this->beads)
            item += move;
    }

    void rescale(double rescale)
    {
        for(Atom& item : this->beads)
            item *= rescale;
    }

    bool isSame (vector<Atom>& container, Atom& push)
    {
        bool same = false;

        for(unsigned int q=0; q< container.size(); q++) {
            if(container[q].isAproxSame(push)) {
                same = true;
            }
        }

        return !same;
    }

    void add(Data& data)
    {
        bool overlap = false;
        for(auto& item : this->beads)
        {
            if( data.isOverlap(item) )
            {
                overlap = true;
            }
        }

        data.all_beads.insert(data.all_beads.end(), this->beads.begin(), this->beads.end());
        data.all_bonds.insert(data.all_bonds.end(), this->bonds.begin(), this->bonds.end());
        data.all_angles.insert(data.all_angles.end(), this->angles.begin(), this->angles.end());

        data.all_sigma = this->sigma;
        data.all_sigma_size = this->sigma_size;
        data.all_sigma_cosatt = this->sigma_cosatt;

        if(overlap)
        {
            cerr << "WARNING: overlap detected" << endl;
        }
    }


    void clusterRotate_random(vector<Atom>& cluster, double max_angle) {
        double vc,vs;
        Atom newaxis;

        Atom cluscm = clusterCM(cluster);

        // create rotation quaternion
        newaxis.randomUnitSphere(); // random axes for rotation
        vc = cos(max_angle * ( 0.0001 * (rand()%10000) ) );
        if (( 0.0001 * (rand()%10000) ) <0.5) vs = sqrt(1.0 - vc*vc);
        else vs = -sqrt(1.0 - vc*vc); // randomly choose orientation of direction of rotation clockwise or counterclockwise

        Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

        //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

        //shift position to geometrical center
        for(unsigned int i=0; i<cluster.size(); ++i) {
            //shift position to geometrical center
            cluster[i].x -= cluscm.x;
            cluster[i].y -= cluscm.y;
            cluster[i].z -= cluscm.z;
            //do rotation
            cluster[i].rotate(newquat);
            //shift positions back
            cluster[i].x += cluscm.x;
            cluster[i].y += cluscm.y;
            cluster[i].z += cluscm.z;
        }
    }

    void clusterRotate(vector<Atom>& cluster, double angle, Atom axis) {
        double vc,vs;

        Atom cluscm = clusterCM(cluster);

        axis.normalise();

        vc = cos( angle );
        vs = sqrt(1.0 - vc*vc);

        Quat newquat(vc, axis.x*vs, axis.y*vs, axis.z*vs);

        //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

        //shift position to geometrical center
        for(unsigned int i=0; i<cluster.size(); ++i) {
            //shift position to geometrical center
            cluster[i].x -= cluscm.x;
            cluster[i].y -= cluscm.y;
            cluster[i].z -= cluscm.z;
            //do rotation
            cluster[i].rotate(newquat);
            //shift positions back
            cluster[i].x += cluscm.x;
            cluster[i].y += cluscm.y;
            cluster[i].z += cluscm.z;
        }
    }


    Atom clusterCM(vector<Atom>& cluster) {

        Atom cluscm(0.0, 0.0, 0.0);

        for(unsigned int i=0; i<cluster.size(); i++) {
            cluscm.x += cluster[i].x;
            cluscm.y += cluster[i].y;
            cluscm.z += cluster[i].z;
        }

        cluscm.x /= cluster.size();
        cluscm.y /= cluster.size();
        cluscm.z /= cluster.size();

        return cluscm;
    }
};


class LoadParticle : public Particle
{
public:
    LoadParticle() : Particle("LoadParticle") {}
    LoadParticle(string name) : Particle(name) {}
    void generate( Data& data ) {}
};

#endif // PARTICLE_H
