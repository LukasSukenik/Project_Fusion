#ifndef SPHERE_H
#define SPHERE_H

#include "particle.h"

using namespace std;

class Sphere: public Particle {
public:
    Sphere() : Particle("Sphere") {}
    Sphere(string name) : Particle(name) {}

    void generate( Data& data )
    {
        vector<Atom> ligand;
        int nano_start = beads.size();

        fibonacci_sphere( beads, data.in.num_of_beads, typeNano);
        fibonacci_sphere_z_distrib_linear( ligand, data.in.num_lig, data.in.c, typeTemp); // second fib. sphere
        Atom patch = Atom(1,1,1,typeLig);
        gen_ligands( data, ligand, patch, typeNano); // find closest spheres on first fib. sphere and change their type


        beads.erase(beads.begin()+data.in.num_of_beads, beads.end()); // erase second fib sphere

        int nano_end = beads.size();
        for(int i=nano_start; i<nano_end; ++i) {
            beads[i] = ((beads[i]*data.in.scale) + data.in.com_pos);
            beads[i].mol_tag = 2;
        }
    }


    void fibonacci_sphere(vector<Atom>& container, int samples, int type, int mol_tag=0) {
        const double PI = 3.141592653589793;
        double offset = 2.0/samples;
        double increment = PI * (3.0 - sqrt(5.0));
        double x,y,z,r,phi;

        for(int i=0; i<samples; ++i) {
            z = ((i * offset) - 1) + (offset / 2);
            r = sqrt(1 - z*z);

            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;

            container.push_back(Atom(x,y,z,type, mol_tag));
        }
    }

    void move(Atom move) {
        for(Atom& item : this->beads)
            item += move;
    }

    void rescale(double rescale) {
        for(Atom& item : this->beads)
            item *= rescale;
    }

protected:    
    void gen_ligands( Data& data, vector<Atom>& ligand, Atom patch, int type_from, int type_temp=-1)
    {
        for(auto& lig : ligand)
        {
            if(lig.x < patch.x && lig.x > patch.vx &&
               lig.y < patch.y*data.in.c && lig.y > patch.vy*data.in.c &&
               lig.z < patch.z && lig.z > patch.vz)
            {
                Atom* select = &beads[0]; // Stupid C++, for some reason reference dont work
                for(auto& item : beads)
                {
                    if(select->dist(lig) > item.dist(lig) && item.type == type_from )
                    {
                        select = &item;
                    }
                }
                if( select->type == type_from )
                    select->type = patch.type;
            }
        }
    }

    /**
     * @brief fibonacci_sphere_z_distrib_linear - We are transforming z axis coordinates where we generate particles. In sphere its linear - slope 1
     *                                          - We want x from -1 to 0
     * @param samples
     * @param type
     */
    void fibonacci_sphere_z_distrib_linear(vector<Atom>& container, int samples, double S, int type) {
        const double PI = 3.141592653589793;
        double offset = 2.0/samples;
        double increment = PI * (3.0 - sqrt(5.0));
        double x,y,z,r,phi;

        for(int i=0; i<samples; ++i) {

            z = ((i * offset) - 1) + (offset / 2); // z linearly distributed from (-1, -1) to (-1, -1)
            z = transform(z, S);
            r = sqrt(1 - z*z);

            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;

            container.push_back(Atom(x,y,z,type));
        }
    }


    void testDistribution() {
        int size = 20;
        int hist[size] = {0};
        for(int i=0; i<beads.size(); ++i) {
            if(beads[i].type == typeLig) {
                cout << beads[i].z << endl;
                ++hist[ (int)(beads[i].z*20) ];
            }
        }
        for(int i=0; i< size; ++i) {
            //cout << hist[i] << endl;
        }
        exit(1);
    }

private:
    double transform(double z, double param) {
        if(z < 0)
            return -( param * pow(fabs(z), 0.5)) - (1-param);
        else
            return ( param * pow(fabs(z), 0.5)) + (1-param);
    }
};

class SphereJanus : public Sphere
{
public:
    SphereJanus() : Sphere("SphereJanus") {}
    SphereJanus(string name) : Sphere(name) {}

    void generate( Data& data ) {

        // generate nano
        int nano_start = beads.size();
        fibonacci_sphere(beads, data.in.num_of_beads, typeNano);
        int nano_end = beads.size();

        // change type nano to type lig based on placement
        for(int i=nano_start; i<nano_end; ++i) {
            if(beads[i].x > 0.0) {
                beads[i].type = typeLig;
            }
        }

        for(int i=nano_start; i<nano_end; ++i) {
            beads[i] = ((beads[i]*data.in.scale) + data.in.com_pos);
            beads[i].mol_tag = 2;
        }

        /*nano_start = beads.size();
        fibonacci_sphere(data.num_of_beads, typeNano);
        nano_end = beads.size();

        // change type nano to type lig based on placement
        for(int i=nano_start; i<nano_end; ++i) {
            if(beads[i].x < 0.0) {
                beads[i].type = typeLig;
            }
        }

        data.com_pos.x = 50.0-data.com_pos.x;

        for(int i=nano_start; i<nano_end; ++i) {
            beads[i] = ((beads[i]*data.scale) + data.com_pos);
            beads[i].mol_tag = 3;
        }*/
    }
};

#endif // SPHERE_H
