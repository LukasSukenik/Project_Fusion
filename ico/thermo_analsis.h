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

        // Einstein-Stokes
        double D=kb*T/(6*pi*vis*r); // m^2 . s

        Welford_Algo aver[4];
        Atom origin(0,0,0);
        Atom current;

        double dist_sq;
        vector<int> steps;
        vector<double> cx, cy, cz;
        int jump=5000;

        while( !fs.eof() ) // Lines in input, For (steps), eof = end of file
        {
            // step temp pe ke c_gyr1 c_cm1[1] c_cm1[2] c_cm1[3]
            fs >> step >> temp >> pe >> ke >> gyr >> cm_x >> cm_y >> cm_z;
            if( (step % jump ) == 0 )
            {
                steps.push_back(step);
                cx.push_back(cm_x);
                cy.push_back(cm_y);
                cz.push_back(cm_z);
            }
        }
        fs.close();


        for(int jj = 1; jj <= step; ++jj) // For ( Jump )
        {
            origin = Atom(0,0,0);
            for (int i=0; i<steps.size(); ++i) // For ( Steps )
            {
                if( (steps[i] % ( jump*(jj) ) ) == 0)
                {
                    current.x = cx[i];
                    current.y = cy[i];
                    current.z = cz[i];

                    dist_sq = origin.distSQ(current);
                    aver[jj].push(dist_sq);

                    origin = current;

                }
            }
            cout << aver[jj].getAverage() << endl;
        }


        return 0;
    }
};

#endif // THERMO_ANALSIS_H


/*
 * BEGIN{
        # constants
        kb=1.38e-23 # m^2 . kg / (s^2 K)
        T=310 # K
        pi=3.141592
        vis=0.00089 # J.s/m^3
        r=5e-9 # m

        # simulation settings
        every=5000

        # Einstein-Stokes
        D=kb*T/(6*pi*vis*r) # m^2 . s

        cumdis=0
        x=-1
        y=-1
        z=-1
        bool_beg_set=0 # 0=false, 1=true
}
{

        #expected=6*45*0.01*$1

        if($1==0 && bool_beg_set == 0)
        {
          x=$6
          y=$7
          z=$8
        }

        if( ($1 % (10*every)) == 0 )
       {
          cumdis=(($6-x)**2+($7-y)**2+($8-z)**2) # cumulative distance from beginning
          #expected=6*D*$1

          print $1, cumdis
        }
}

 * */
