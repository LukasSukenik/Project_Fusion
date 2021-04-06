#ifndef VESICLE_H
#define VESICLE_H

#include "membrane.h"
#include "input.h"

class Vesicle : public Membrane {
public:
    Vesicle() {}
    Vesicle(string in_file) : Membrane(in_file) {}

    /**
     * @brief generate - Generate lipid vesicle using Fibonacci sphere algorithm
     * @param in
     */
    void generate(Input& in)
    {
        vec.clear();
        bonds.clear();

        cout << "Vesicle radius: " << in.radius << ", number of lipids: " << in.num_lipids << endl;
        int samples = in.num_lipids;

        const double frac = 0.1;
        int old_samples;
        const double PI = 3.141592653589793;
        double offset = 2.0/samples;
        double increment = PI * (3.0 - sqrt(5.0));
        double x,y,z,r,phi;
        int num=1;

        for(int i=0; i<samples; ++i) {
            z = ((i * offset) - 1) + (offset / 2);
            r = sqrt(1 - z*z);

            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;

            lipids.push_back( Lipid(Particle(num + in.offset,   in.mol_tag, HEAD, x,            y,            z),
                                    Particle(num+1 + in.offset, in.mol_tag, TAIL, x + x*frac,   y + y*frac,   z + z*frac),
                                    Particle(num+2 + in.offset, in.mol_tag, TAIL, x + 2*x*frac, y + 2*y*frac, z + 2*z*frac) ) );
            num += 3;
        }

        int sum=lipids.size();
        int count=0;
        int random=0;
        while(count < in.num_rec)
        {
            random = (int)(rand()/(RAND_MAX + 1.) * samples);
            if(random < lipids.size())
            {
                if(lipids[random].part[0].atom_type == HEAD)
                {
                    lipids[random].part[0].atom_type = RECEPTOR_HEAD;
                    count++;
                }
            }
        }

        //min_dist(scale, 0, samples);
        old_samples = samples;

        samples = samples*(1.0+3*frac)*(1.0+3*frac)*(1.06633894589390435/1.122462048309373)*(1.06633894589390435/1.122462048309373);
        offset = 2.0/samples;

        // inner leaflet
        for(int i=0; i<samples; ++i) {
            z = ((i * offset) - 1) + (offset / 2);
            r = sqrt(1 - z*z);

            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;

            lipids.push_back( Lipid(Particle(num + in.offset,   in.mol_tag, HEAD, x + 5*x*frac, y + 5*y*frac, z + 5*z*frac),
                                    Particle(num+1 + in.offset, in.mol_tag, TAIL, x + 4*x*frac, y + 4*y*frac, z + 4*z*frac),
                                    Particle(num+2 + in.offset, in.mol_tag, TAIL, x + 3*x*frac, y + 3*y*frac, z + 3*z*frac) ) );
            num += 3;
        }

        count=0;
        random=0;
        while(count < in.num_rec)
        {
            random = (int)(rand()/(RAND_MAX + 1.) * samples) + sum;
            if(random < lipids.size())
            {
                if(lipids[random].part[0].atom_type == HEAD)
                {
                    lipids[random].part[0].atom_type = RECEPTOR_HEAD;
                    count++;
                }
            }
        }


        // generate bonds
        for(int i=0; i<lipids.size(); i++) {
            lipids[i].changeN(3*i+1, 1); // Head, tail, tail

            for(int j=0; j<3; j++) {
                vec.push_back(lipids[i].part[j]);
                bonds.push_back(lipids[i].bond[j]);
            }
        }



        for(unsigned int i=0; i<vec.size(); i++)
        {
            vec[i].x = in.radius*vec[i].x + in.com_x;
            vec[i].y = in.radius*vec[i].y + in.com_y;
            vec[i].z = in.radius*vec[i].z + in.com_z;
        }
    }
};

#endif // VESICLE_H
