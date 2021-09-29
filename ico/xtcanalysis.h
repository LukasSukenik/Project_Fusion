#ifndef XTCANALYSIS_H
#define XTCANALYSIS_H

#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>

#include "welford.h"

#include "xdrfile-1.1.4/include/xdrfile.h"
#include "xdrfile-1.1.4/include/xdrfile_xtc.h"

using namespace std;



class XTCAnalysis
{
private:
    char fileName[64];
    int natoms;

    int step;
    int release_step = numeric_limits<int>::max();
    float time;
    matrix box;
    float prec;

public:
    Welford_Algo average;

    XTCAnalysis()
    {
        // -1 for all fold indexes
    }


    /**
     * @brief potential - cos^2 potential
     * @param dist
     * @param epsilon
     * @param sigma
     * @param rm
     * @return
     */
    double potential(double dist, double epsilon, double sigma, double rm)
    {
        if(dist > rm+sigma)
            return 0.0;
        if(dist > rm)
        {
            const double pih = std::atan(1.0)*2;
            return -epsilon * cos( pih * (dist - rm) / sigma ) * cos( pih * (dist - rm) / sigma );
        }
        else
        {
            // https://en.wikipedia.org/wiki/Lennard-Jones_potential
            // -epsilon from cos^2 potential, eps for LJ is 1.0
            double ssigma = rm / pow(2.0, 1.0/6.0);
            return 4*epsilon*( pow(ssigma/dist, 12) - pow(ssigma/dist,6) );
        }
    }

    int analyse(string inName, int stop=-1)
    {
        int status=exdrOK;
        double dist,max=0;

        int first_step;

        // convert string to char
        int i;
        for (i = 0; inName[i] != '\0'; ++i) {
          fileName[i] = inName[i];
        }
        fileName[i] = '\0';

        status = read_xtc_natoms(fileName, &natoms);
        if (status == exdrOK)
        {
            XDRFILE* xfp = xdrfile_open(fileName, "r");
            if (xfp != NULL)
            {
                rvec k[natoms];
                status = read_xtc(xfp, natoms, &first_step, &time, box, k, &prec);

                while(status == exdrOK && ( stop == -1 || step < first_step+stop) ) // Analyze frame
                {
                    status = read_xtc(xfp, natoms, &step, &time, box, k, &prec);

                    i = 0; // atom index
                    double x = k[i][0];
                    double y = k[i][1];
                    double z = k[i][2];
                }
                xdrfile_close(xfp);
            }
            else
            {
                cout << "File not opened:" << fileName << endl;
            }
        }
        else
        {
            fprintf(stderr, "read_xtc_natoms failure; return code %d", status);
            cerr << inName << endl;
        }
        return -1;
    }
};

#endif // XTCANALYSIS_H
