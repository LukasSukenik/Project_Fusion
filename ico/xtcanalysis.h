#ifndef XTCANALYSIS_H
#define XTCANALYSIS_H

#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <stack>
#include <tuple>
#include <iterator>

#include "welford.h"
#include "atom.h"
#include "data.h"

#include "xdrfile-1.1.4/include/xdrfile.h"
#include "xdrfile-1.1.4/include/xdrfile_xtc.h"

using namespace std;
using namespace chrono;

class Histogram
{
private:
    int bin_size = 150;
public:
    Histogram(){}

    void init_frame()
    {
        bin.push_back( vector<int>(bin_size, 0) );
    }

    void increment_last(int index)
    {
        ++bin.back()[index]; // bin.back() refers to frame, the [index] refers to particle == bin[bin.size()-1][index]
    }

    void print()
    {
        return;
    }

    vector< vector<int> > bin;
};

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

    const double degToRad = 1.0/57.2958;
    const double radToDeg = 57.2958; // 180deg = 3.1415 radians, 1rad = 57.2958deg

    enum{MOLTAG_NANO=2};

    int grid_size=200;

public:
    Welford_Algo average;

    XTCAnalysis()
    {

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

    int symmetry( Atom test, Atom cm, rvec* frame, Data& data)
    {
        // do the radial symmetry test
        //align structure to z axis

        Atom axis = test-cm;                          // Axis between com and particle at the end
        axis.normalise();                                 // normalise axis vector for correct rotation
        Atom z_axis = Atom(0,0,1);
        Atom rot_axis = axis.cross(z_axis);
        rot_axis.normalise();                                  // normalise for rotation algo
        double angle = acos( axis.dot(z_axis) );

        Atom x2=axis;

        x2.rotate(rot_axis, angle);
        cm.rotate(rot_axis, angle);

        int posx=0;
        int negx=0;
        int posy=0;
        int negy=0;
        int diff1=0;
        int diff2=0;
        int diff3 = 0;
        int poso=0;
        int nego=0;

        for(int j=0; j<data.temp_beads.size(); ++j)
        {
            test.x = frame[j][0];
            test.y = frame[j][1];
            test.z = frame[j][2];
            test.mol_tag = data.temp_beads[j].mol_tag;

            // determine whether you want +angle or - angle
            if( radToDeg*acos( x2.dot(z_axis) ) < 0.000001)
            {
                 test.rotate(rot_axis, angle);
            }
            else{
                test.rotate(rot_axis, -angle);
            }

            //data.temp_beads[j] = test; // just for testing xyz in vmd

            if(test.mol_tag != MOLTAG_NANO && test.x-cm.x > 0)
            {
                ++posx;
            }
            if(test.mol_tag != MOLTAG_NANO && test.x-cm.x < 0)
            {
                ++negx;
            }
            if(test.mol_tag != MOLTAG_NANO && test.y-cm.y > 0)
            {
                ++posy;
            }
            if(test.mol_tag != MOLTAG_NANO && test.y-cm.y < 0)
            {
                ++negy;
            }
            if(test.mol_tag != MOLTAG_NANO && test.y-cm.y > test.x-cm.x)
            {
                ++poso;
            }
            if(test.mol_tag != MOLTAG_NANO && test.y-cm.y < test.x-cm.x)
            {
                ++nego;
            }
        }
        diff1 = abs(posx - negx);
        diff2 = abs(posy - negy);
        diff3 = abs(poso - nego);

        // void printXYZ(Data& data, Atom axis, Atom com)
        //printXYZ(data, axis, cm);

        return diff1+diff2+diff3;

    }


    double median(vector<double>& arr)
    {
        int n = arr.size();
        int count=0;
        //sorting - ASCENDING ORDER
        for(int i=0;i<n;i++)
        {
            ++count;
            for(int j=i+1;j<n;j++)
            {
                if(arr[i]>arr[j])
                {
                    int temp  =arr[i];
                    arr[i]=arr[j];
                    arr[j]=temp;
                }
            }
        }
        int median=count/2;
        return arr[median];
    }


    Atom com(rvec* frame, Data& data)
    {
        Atom cm;
        int count=0;
        for(int i=0; i<data.temp_beads.size(); ++i)
        {
            if(data.temp_beads[i].mol_tag != 2)
            {
                cm.x += frame[i][0];
                cm.y += frame[i][1];
                cm.z += frame[i][2];
                ++count;
            }
        }
        cm *= 1.0/count;
        return cm;
    }
    Atom com_ves1(rvec* frame, Data& data)
    {
        Atom cm;
        int count=0;
        for(int i=0; i<data.temp_beads.size(); ++i)
        {
            if(data.temp_beads[i].mol_tag == 1)
            {
                cm.x += frame[i][0];
                cm.y += frame[i][1];
                cm.z += frame[i][2];
                ++count;
            }
        }
        cm *= 1.0/count;
        return cm;
    }
    Atom com_ves2(rvec* frame, Data& data)
    {
        Atom cm;
        int count=0;
        for(int i=0; i<data.temp_beads.size(); ++i)
        {
            if(data.temp_beads[i].mol_tag == 3)
            {
                cm.x += frame[i][0];
                cm.y += frame[i][1];
                cm.z += frame[i][2];
                ++count;
            }
        }
        cm *= 1.0/count;
        return cm;
    }

    int distance(int i, rvec* frame, Atom cm)
    {
        Atom part;
        int COM;
        int centre;

        part.x = frame[i][0];
        part.y = frame[i][1];
        part.z = frame[i][2];

        COM=sqrt(part.x*part.x + part.y*part.y + part.z*part.z);
        centre=sqrt(cm.x*cm.x + cm.y*cm.y + cm.z*cm.z);
        int rest = COM -centre;
        if(rest<0)
        {
            rest=rest*-1;
        }
        return rest;
    }

    void printXYZ(Data& data, Atom axis, Atom com)
    {
        std::fstream fout( "text.xyz", std::fstream::out );
        fout << data.temp_beads.size()+100+100 << "\nparticle\n";
        for(unsigned int i=0; i<data.temp_beads.size(); ++i) {
            fout << "C" << data.temp_beads[i].type <<  " " << data.temp_beads[i] << endl;
        }
        axis *= 0.1;
        for(int i=-50; i<50; ++i) {
            fout << "H" << 1 <<  " " << com+(axis*i) << endl;
        }
        Atom z_axis(0,0,2);
        for(int i=-50; i<50; ++i) {
            fout << "P" << 1 <<  " " << com+(z_axis*i) << endl;
        }
        fout.close();
    }

    /**
     * @brief get_cutoff - determine the distance cutoff - median distance, sorting algorithm, vector<double> test_array; // double test_array[500];
     * @param data
     * @return
     */
    double get_cutoff(Data& data)
    {
        vector<double> cmarray(data.temp_beads.size() , 0); //calculate median for cutoff
        Atom cm;
        //
        for(int a=0; a<data.temp_beads.size(); ++a) // loop over particles to build cm array
        {
            if(data.temp_beads[a].mol_tag != 2)
            {
                cmarray[a] = cm.distSQ( data.temp_beads[a] );
            }
        }
        return 1.5*sqrt( median(cmarray) );
    }

    Atom get_axis(Data& data, rvec* frame, double cutoff)
    {
        int temp=10000;
        Atom id, cm, test;

        for(int i=0; i<data.temp_beads.size(); ++i) // loop over particles to construct axis through
        {
            test.x = frame[i][0];
            test.y = frame[i][1];
            test.z = frame[i][2];
            test.mol_tag = data.temp_beads[i].mol_tag;

            if(data.temp_beads[i].mol_tag != MOLTAG_NANO && (test-cm).size() > cutoff) // look at vesicle caps only
            {
                int diff = symmetry(test, cm, frame, data);
                if(diff<temp)
                {
                    temp=diff;
                    id=test;
                }
            }
        }
        return id;
    }

    Atom right_angle(Atom x3, Atom z_axis, Atom part, Atom rot_axis, double angle){
        if( radToDeg*acos( x3.dot(z_axis) ) < 0.000001)
        {
            part.rotate(rot_axis, angle);
        }
        else{
            part.rotate(rot_axis, -angle);
        }
        return part;
    }


    int analyze_histogram(string inName, Data& data, int stop=-1)
    {
        int status=exdrOK;
        double dist,max=0;

        Atom cm;
        Atom id;
        double cutoff;

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
                rvec frame[natoms];
                rvec **pointer;
                status = read_xtc(xfp, natoms, &first_step, &time, box, frame, &prec);

                Histogram hist;
                //auto t1 = high_resolution_clock::now();
                cout << "build start " << endl;

                while(status == exdrOK && ( stop == -1 || step < first_step+stop) ) // Analyze frame
                {
                    //cout << "#step: " << step << endl;
                    status = read_xtc(xfp, natoms, &step, &time, box, frame, &prec); // reads 1 entire frame from xtc
                    cm = com(frame, data);

                    cutoff = get_cutoff(data);
                    id = get_axis(data, frame, cutoff);

                    //rotate structure by cm-id axis
                    Atom z_axis = Atom(0,0,1);
                    Atom axis = id-cm;                          // Axis between com and selected particle (id)
                    Atom start=com_ves1(frame, data);
                    Atom start2=com_ves2(frame, data);

                    //printXYZ(data, axis, cm);

                    axis.normalise();                                 // normalise axis vector for correct rotation

                    Atom rot_axis = axis.cross(z_axis);
                    rot_axis.normalise();                                  // normalise for rotation algo
                    cm.normalise();
                    Atom particle;
                    double angle = acos( axis.dot(z_axis) );

                    Atom x3=axis;
                    //cout << axis << " " << radToDeg*acos( axis.dot(z_axis) ) << endl;
                    x3.rotate(rot_axis, angle);

                    //histogram, bins
                    int histo_size=150;
                    vector<int> bin(histo_size, 0);

                    //vector<vector<int> > grid(grid_size+1, vector<int> (grid_size+1, 0));
                    //
                    // int > tuple<int,bool>, get<0>(name) -> int get<1>(name) -> bool
                    //
                    vector<vector<vector< tuple<int,bool> >>> grid(grid_size, vector<vector<tuple<int,bool>>>(grid_size, vector<tuple<int,bool>>(grid_size)));
                    //
                    // build the grid
                    //
                    for (int c=0; c<grid_size;++c)
                    {
                        for (int b=0; b<grid_size;++b)
                        {
                            for (int a=0; a<grid_size;++a)
                            {
                                get<0>(grid[c][b][a])=0;
                                get<1>(grid[c][b][a])=false;
                            }
                        }
                    }


                    hist.init_frame();


                    // set grid origin
                    int x_offset = 20;
                    int y_offset = 20;
                    int z_offset = 20;
                    double bin_size=0.25;
                    Atom grid_point;
                    right_angle(x3, z_axis, start, rot_axis, angle);

                    start.x=round( (start.x +x_offset)/bin_size );
                    start.y=round( (start.y +y_offset)/bin_size );
                    start.z=round( (start.z +z_offset)/bin_size );

                    right_angle(x3, z_axis, start2, rot_axis, angle);

                    start2.x=round( (start2.x +x_offset)/bin_size );
                    start2.y=round( (start2.y +y_offset)/bin_size );
                    start2.z=round( (start2.z +z_offset)/bin_size );


                    for(int j=0; j<data.temp_beads.size(); ++j) // ~10 000
                    {
                        particle.x = frame[j][0];
                        particle.y = frame[j][1];
                        particle.z = frame[j][2];
                        particle.type = data.temp_beads[j].type;


                        // determine whether you want +angle or - angle and rotate

                        right_angle(x3, z_axis, particle, rot_axis, angle);
                        right_angle(x3, z_axis, cm, rot_axis, angle);


                        build_grid(particle,grid);

                        //fill histo
                        for (int b=0; b<histo_size;b++)
                        {
                            if (particle.type == 2 && (particle.z - cm.z + 40) >= 0.6*b && (particle.z - cm.z + 40) < 0.6*(b+1))
                            {
                                bin[b]++;
                                hist.increment_last(b);
                            }
                        }
                    }


                    //Build neighbor list. Nneigh-->list with number of neigh for each particle; neigh_index --> list of neigh of each particle.
                    Atom at;
                    vector <int> Nneigh;
                    vector <int> neigh_index;
                    vector <bool> used(data.temp_beads.size(), false);
                    vector <bool> added(data.temp_beads.size(), false);

                    for(int i=0; i<data.temp_beads.size(); ++i) // ~10 000
                    {
                        particle.x = frame[i][0];
                        particle.y = frame[i][1];
                        particle.z = frame[i][2];
                        particle.type = data.temp_beads[i].type;

                        int neigh=0;

                        for(int j=0; j<data.temp_beads.size(); ++j) // ~10 000
                        {
                            at.x = frame[j][0];
                            at.y = frame[j][1];
                            at.z = frame[j][2];
                            at.type = data.temp_beads[j].type;

                            if(at.distSQ(particle) < pow(2.0, 1.0/3.0) ){
                                neigh++;
                                neigh_index.push_back(j);
                            }
                        }
                        Nneigh.push_back(neigh);
                    }


                    //add particles to cluster
                    int end;
                    vector <int> cluster;

                    cluster.push_back(1);

                    for (int j=0;j<cluster.size();++j) {

                        //seguir solo si j no ha sido comprobado antes(if not visited)
                        if(data.temp_beads[cluster[j]].type == 2 && used[cluster[j]]==false){
                            used[cluster[j]]=true;
                            int begin=0;
                            for (int a=0;a<cluster[j];++a) {
                                begin=begin+Nneigh[a];
                            }

                            end=begin+Nneigh[cluster[j]]-1;

                            for (int n=begin;n<=end;++n) {
                                //cout<<Nneigh[cluster[j]]<<" "<<neigh_index[n]<<endl;

                                //if neigh is not already in cluster, then add
                                if(added[neigh_index[n]]==false){
                                    added[neigh_index[n]]=true;
                                    cluster.push_back(neigh_index[n]);
                                }
                            }
                        }
                    }
                    if(cluster.size() < 0.8 * data.temp_beads.size()){
                        cout<<"No stalk"<<endl;

                    }else{
                        cout<<"Stalk"<<endl;
                    }

                    //
                    // analyze grid
                    //
                    if( is_pore(start, grid,start2) )
                        cout << step << " Pore" << endl;

                    else {
                        cout << step << " NO PORE" << endl;
                    }

                    //print grid
                    /*for (int n=0;n<grid_size;++n){
                        for (int m=0;m<grid_size;++m){
                            for (int l=0;l<grid_size;++l){
                                if (l!=grid_size-1)
                                {
                                    if (get<1>(grid[l][m][n])==true){
                                        cout << "\033[1;31m"<<get<0>(grid[l][m][n])<<"\033[0m"<<" ";
                                    }else{
                                        cout <<get<0>(grid[l][m][n])<<" ";
                                    }
                                }
                                else{
                                    if (get<1>(grid[l][m][n])==true){
                                        cout << "\033[1;31m"<<get<0>(grid[l][m][n])<<"\033[0m"<<endl;
                                    }else{
                                        cout <<get<0>(grid[l][m][n])<<endl;
                                    }
                                }
                            }
                        }
                    }*/

                    //int frame = 0;
                    /*for (int n=0;n<histo_size;++n){   //print histo
                        cout << n << " " << bin[n] << /*" " << hist.bin[frame][n] << endl;
                    }*/

                }
                xdrfile_close(xfp);
                /*auto t2 = high_resolution_clock::now();
                auto ms_int = duration_cast<milliseconds>(t2 - t1);
                std::cout << "duration" << ms_int.count() << "ms\n";*/
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

    int build_grid(Atom particle, vector<vector<vector< tuple<int,bool> >>>& grid){
        Atom grid_point;
        int x_offset = 20;
        int y_offset = 20;
        int z_offset = 20;
        double bin_size=0.25;
        int cut, lowx, highx, lowy, highy, lowz, highz;
        grid_point.x = (particle.x +x_offset)/bin_size ;
        grid_point.y = (particle.y +y_offset)/bin_size ;
        grid_point.z = (particle.z +z_offset)/bin_size ;

        // Find grip point surrounding
        cut=(pow(2.0, 1.0/6.0))/bin_size +1;
        lowx=round(grid_point.x -cut);
        if(lowx < 0)
            lowx = 0; //lowx = (lowx < 0) ? 0 : lowx;
        highx=round(grid_point.x +cut);
        if(highx >= grid_size)
            highx = grid_size-1;
        lowy=round(grid_point.y -cut);
        if(lowy < 0)
            lowy = 0;
        highy=round(grid_point.y +cut);
        if(highy >= grid_size)
            highy = grid_size-1;
        lowz=round(grid_point.z -cut);
        if(lowz < 0)
            lowz = 0;
        highz=round(grid_point.z +cut);
        if(highz >= grid_size)
            highz = grid_size-1;

        // assing gridpoing to particle
       // if (particle.type == 2){
            for (int pointx=lowx;pointx<=highx;pointx++){
                for (int pointy=lowy;pointy<=highy;pointy++){
                    for (int pointz=lowz;pointz<=highz;pointz++)
                    {
                        grid_point = Atom( pointx*bin_size - x_offset , pointy*bin_size - y_offset, pointz*bin_size - z_offset);
                        //cout << grid_point.x << "~=" << particle.x << "   " << grid_point.y << "~=" << particle.y << "   " << grid_point.z << "~=" << particle.z << "   " << endl;
                        if ( grid_point.distSQ(particle) < pow(2.0, 1.0/3.0)+pow(bin_size, 2.0) )
                        {
                            ++get<0>(grid[pointx][pointy][pointz]);
                            //cout << pointx << " " << pointy << " " << pointz << endl;
                        }
                    }

                }
            }
       // }
        return 0;
    }

    bool is_pore(Atom start, vector<vector<vector< tuple<int,bool> >>>& grid, Atom start2)
    {
        //
        // analyze grid
        //
        stack<int> instack;
        int pore_out=0;
        // push i * x_size*y_size + j*y_size + k
        // pop(), only removes last element, use top() then pop()
        int indx=start.x+start.y*grid_size+start.z*grid_size*grid_size;
        int x,y,z;
        instack.push(indx);

        int dx[] = {0, 1,  0, -1,  0,  0};
        int dy[] = {0, 0,  1,  0, -1,  0};
        int dz[] = {1, 0,  0,  0,  0, -1};
        while(!instack.empty()){

            z = instack.top() / (grid_size * grid_size);
            instack.top() -= (z * grid_size * grid_size);
            y = instack.top() / grid_size;
            x=instack.top() % grid_size;

            instack.pop(); // remove that index from stack
            //check neighbors
            for (int d=0; d < 6; ++d)
            {
                int neigh_x=x+dx[d];
                int neigh_y=y+dy[d];
                int neigh_z=z+dz[d];
                //
                // Check boundaries
                //
                if(neigh_x>=0 && neigh_y>=0 && neigh_z>=0 && neigh_x<grid_size && neigh_y<grid_size && neigh_z<grid_size) {
                    if(get<0>(grid[neigh_x][neigh_y][neigh_z])==0 && get<1>(grid[neigh_x][neigh_y][neigh_z])==false)
                    {
                        //replace neighbors with new neigh whose value is 0
                        instack.push(neigh_x+neigh_y*grid_size+neigh_z*grid_size*grid_size); //add new neighbors to stack
                        get<1>(grid[neigh_x][neigh_y][neigh_z])=true;
                    }
                }
                else
                {
                    return true;
                }
            }
        }

        if(get<1>(grid[start2.x][start2.y][start2.z])==true){
            return true;
        }
        return false;
    }
};

#endif // XTCANALYSIS_H
