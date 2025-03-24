#ifndef THERMO_ANALSIS_H
#define THERMO_ANALSIS_H

#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "welford.h"
#include "atom.h"
#include "data.h"

class Thermo_analysis
{
public:
    Thermo_analysis (){}

//
// constants
//
double kb=1.38e-23; // m^2 . kg / (s^2 K)
double T=310; // K
double pi=3.141592;
double vis=0.00089; // J.s/m^3
double r=5e-9; // m


    int timestep_analyze()
    {
        std::fstream fs( "thermo", std::fstream::in );

        double temp, pe, ke, gyr, cm_x, cm_y, cm_z;
        int step;
        int jump_limit = 150;

        // Einstein-Stokes
        ofstream file;
        file.open ("D_constant");
        double D=kb*T/(6*pi*vis*r); // m^2 / s
        D=D*10e12; //um2/s
        file<<D<<endl;
        file.close();

        Welford_Algo aver[jump_limit];
        Atom origin(0,0,0);
        Atom current;

        double dist_sq;
        vector<int> steps;
        vector<double> cx, cy, cz;
        int persistence_length=5000; // sim setup is thermo every 5000

        while( !fs.eof() ) // Lines in input, For (steps), eof = end of file
        {
            // step temp pe ke c_gyr1 c_cm1[1] c_cm1[2] c_cm1[3]
            fs >> step >> temp >> pe >> ke >> gyr >> cm_x >> cm_y >> cm_z;
            if( (step % persistence_length ) == 0 )
            {
                steps.push_back(step);
                cx.push_back(cm_x);
                cy.push_back(cm_y);
                cz.push_back(cm_z);
            }
        }
        fs.close();

        for(int jj = 1; jj < jump_limit; ++jj) // For ( Jump length in persistence length )
        {
            origin = Atom(0,0,0);
            for (int i=0; i<steps.size(); ++i) // For ( Steps )
            {
                if( (steps[i] % ( persistence_length*(jj) ) ) == 0)
                {
                    current.x = cx[i];
                    current.y = cy[i];
                    current.z = cz[i];

                    dist_sq = origin.distSQ(current);
                    //cout << dist_sq<<endl;
                    aver[jj].push(dist_sq);

                    origin = current;
                }
            }
            cout << jj << " " << aver[jj].getNum() << " "<< aver[jj].getAverage() << endl;
        }
        return 0;
    }
};

#endif // THERMO_ANALSIS_H
