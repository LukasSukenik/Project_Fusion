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

const double degToRad = 1.0/57.2958;
const double radToDeg = 57.2958; // 180deg = 3.1415 radians, 1rad = 57.2958deg
enum{MOLTAG_NANO=2};

Atom right_angle(Atom x3, Atom z_axis, Atom part, Atom rot_axis, double angle)
{
    if( radToDeg*acos( x3.dot(z_axis) ) < 0.000001)
    {
        part.rotate(rot_axis, angle);
    }
    else{
        part.rotate(rot_axis, -angle);
    }
    return part;
}

Atom com(rvec* frame, Data& data, vector<int> mol_tag)
{
    Atom cm;
    int count=0;
    for(int i=0; i<data.temp_beads.size(); ++i)
    {
        for(int& mtag : mol_tag)
        {
            if(data.temp_beads[i].mol_tag == mtag)
            {
                cm.x += frame[i][0];
                cm.y += frame[i][1];
                cm.z += frame[i][2];
                ++count;
            }
        }
    }
    cm *= 1.0/count;
    return cm;
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
                double temp  =arr[i];
                arr[i]=arr[j];
                arr[j]=temp;
            }
        }
    }
    int median=count/2;
    return arr[median];
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




int symmetry( Atom test, Atom cm, rvec* frame, Data& data);



Atom get_axis(Data& data, rvec* frame, double cutoff)
{
    int temp=10000;
    int diff=-1;
    Atom id, test;
    Atom cm=com(frame, data, vector<int>{1,3});

    for(int i=0; i<data.temp_beads.size(); ++i) // loop over particles to construct axis through
    {
        test.x = frame[i][0];
        test.y = frame[i][1];
        test.z = frame[i][2];
        test.mol_tag = data.temp_beads[i].mol_tag;

        if(test.mol_tag != MOLTAG_NANO && (test-cm).size()  > cutoff) // look at vesicle caps only
        {
            diff = symmetry(test, cm, frame, data);
            //cout<<diff<<endl;

            if(diff<temp)
            {
                temp=diff;
                id=test;
            }
        }
    }
      return id;
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
    //cm.rotate(rot_axis, angle);

    int  posx1=0, posx3=0;
    int  negx1=0, negx3=0;
    int  posy1=0, posy3=0;
    int  negy1=0, negy3=0;
    int df1=0, df3=0, df5=0, df7=0;
    int df2=0, df4=0, df6=0, df8=0;
    Atom particle;

    for(int j=0; j<data.temp_beads.size(); ++j)
    {
        particle.x = frame[j][0];
        particle.y = frame[j][1];
        particle.z = frame[j][2];
        particle.mol_tag = data.temp_beads[j].mol_tag;

        // determine whether you want +angle or - angle
        Atom part=particle-cm;
        part=right_angle(x2, z_axis, part, rot_axis, angle);
        Atom cm_rot=cm-cm;
        cm_rot=right_angle(x2, z_axis, cm_rot, rot_axis, angle);

        //data.temp_beads[j] = test; // just for testing xyz in vmd

        if(particle.mol_tag == 1 && part.x-cm_rot.x > 0)
        {
                ++posx1;
        }
        if(particle.mol_tag == 3 && part.x-cm_rot.x > 0){
                ++posx3;
        }

        if(particle.mol_tag ==1 && part.x-cm_rot.x < 0)
        {
            ++negx1;
        }
        if(particle.mol_tag == 3 && part.x-cm_rot.x < 0){
            ++negx3;
        }


        if(particle.mol_tag ==1 && part.y-cm_rot.y > 0)
        {
            ++posy1;
        }
        if(particle.mol_tag == 3 && part.y-cm_rot.y > 0){
            ++posy3;
        }

        if(particle.mol_tag ==1 && part.y-cm_rot.y < 0)
        {
            ++negy1;
        }
        if(particle.mol_tag == 3 && part.y-cm_rot.y < 0){
            ++negy3;
        }
    }

    df1=abs(posx1-posx3);
    df2=abs(posx1-negx1);
    df7=abs(posx3-negx3);
    df3=abs(posy1-posy3);
    df4=abs(posy1-negy1);
    df8=abs(posy3-negy3);
    df5=abs(negx1-negx3);
    df6=abs(negy1-negy3);
     //void printXYZ(Data& data, Atom axis, Atom com);
    //printXYZ(data, axis, cm);

    return df1+df2+df3+df4+df5+df6+df7+df8;

}



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
        ofstream myfile;
        myfile.open ("histogram");
        for (int i=0; i<bin.size(); ++i)
        {
            for (int j=0; j<bin[i].size(); ++j)
            {
                myfile << i << " " << j << " " << bin[i][j] << endl;
            }
        }
        myfile.close();
    }

    void frame_histo(Data& data, rvec* frame)
    {
        Atom cm = com(frame, data, vector<int>{1,3});
        double cutoff = get_cutoff(data);
        Atom id = get_axis(data, frame, cutoff);
        //rotate structure by cm-id axis
        Atom z_axis = Atom(0,0,1);
        Atom axis = id-cm;                          // Axis between com and selected particle (id)
        axis.normalise();  // normalise axis vector for correct rotation

        Atom rot_axis = axis.cross(z_axis);
        rot_axis.normalise();                                  // normalise for rotation algo
        Atom particle;
        double angle = acos( axis.dot(z_axis) );

        Atom x3=axis;
        //cout << axis << " " << radToDeg*acos( axis.dot(z_axis) ) << endl;
        x3.rotate(rot_axis, angle);

        //histogram, bins
        int histo_size=150;
        init_frame();

        for(int j=0; j<data.temp_beads.size(); ++j) // ~10 000
        {
            particle.x = frame[j][0];
            particle.y = frame[j][1];
            particle.z = frame[j][2];
            particle.type = data.temp_beads[j].type;

            // determine whether you want +angle or - angle and rotate

            Atom part=right_angle(x3, z_axis, particle, rot_axis, angle);
            Atom cm_rot=right_angle(x3, z_axis, cm, rot_axis, angle);

            //fill histo
            for (int b=0; b<histo_size;b++)
            {
                if (particle.type == 2 && (part.z - cm_rot.z + 40) >= 0.6*b && (part.z - cm_rot.z + 40) < 0.6*(b+1))
                {
                    increment_last(b);
                }
            }
        }

        /*for (int j=0; j<bin.back().size(); ++j)
        {
            cout << j << " " << bin.back()[j] << endl;
        }*/
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

    int grid_size=200;

public:
    Welford_Algo average;

    XTCAnalysis() {}


    void grid_analysis(Data& data, rvec* frame)
    {
        Atom cm = com(frame, data, vector<int>{1,3});
        double cutoff = get_cutoff(data);
        Atom id = get_axis(data, frame, cutoff);
        Atom z_axis = Atom(0,0,1);
        Atom axis = id-cm;                          // Axis between com and selected particle (id)
        axis.normalise();  // normalise axis vector for correct rotation
        Atom rot_axis = axis.cross(z_axis);
        rot_axis.normalise();                                  // normalise for rotation algo
        Atom particle;
        double angle = acos( axis.dot(z_axis) );
        Atom x3=axis;
        x3.rotate(rot_axis, angle);

        Atom start=com(frame, data, vector<int>{1} );
        Atom start2=com(frame, data, vector<int>{3});

        //vector<vector<int> > grid(grid_size+1, vector<int> (grid_size+1, 0));
        //
        // int > tuple<int,bool>, get<0>(name) -> int get<1>(name) -> bool
        //
        vector<vector<vector< tuple<int,bool> >>> grid(grid_size, vector<vector<tuple<int,bool>>>(grid_size, vector<tuple<int,bool>>(grid_size)));
        init_grid(grid);

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
        }

        //
        // analyze grid
        //
        if( is_pore(start, grid,start2) )
            cout <<step<<" 1 "; // 1 means pore, 0 means no pore

        else {
            cout <<step<<" 0 ";
        }
        //print_grid(grid);
    }

    int get_init_stalk(Data& data, vector<int>& Nneigh)
    {
        int init_stalk=-1;
        for(int p=0; p<data.temp_beads.size(); ++p) // ~10 000
        {
            if(Nneigh[p]>3 && data.temp_beads[p].type == 2) { // choose a particle with 3< neighbors
                init_stalk=p;
                break;
            }
        }
        return init_stalk;
    }

    int analyze_histogram(string inName, Data& data, int stop=-1)
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
                rvec frame[natoms];
                rvec **pointer;
                status = read_xtc(xfp, natoms, &first_step, &time, box, frame, &prec);

                Histogram hist;

                //auto t1 = high_resolution_clock::now();
                //cout << "build start " << endl;
                cout << "#step  pore (0-no, 1-yes) stalk (0-no, 1-yes)"<<endl;

                while(status == exdrOK && ( stop == -1 || step < first_step+stop) ) // Analyze frame
                {
                    status = read_xtc(xfp, natoms, &step, &time, box, frame, &prec); // reads 1 entire frame from xtc

                    //
                    // Histogram - Fusion detection
                    //
                    if( data.in.histo )
                    {
                        hist.frame_histo( data, frame );
                        hist.print();
                    }

                    // Grid for pore search
                    //
                    // Cluster analysis
                    // - Identify Stalk
                    // - Two vesicles comparison sims
                    // -- Average distance of vesicles
                    grid_analysis( data, frame );

                    //
                    //Build neighbor list. Nneigh-->list with number of neigh for each particle; neigh_index --> list of neigh of each particle.
                    //
                    vector <int> Nneigh, neigh_index;
                    generate_pairlist(data, frame, Nneigh, neigh_index);

                    //
                    //  Cluster analysis - stalk
                    //
                    vector <int> cluster, ves1, ves2;
                    int init_stalk=get_init_stalk(data, Nneigh);
                    clusterize(data, cluster, Nneigh, neigh_index, init_stalk, vector<int>{2});

                    if(cluster.size() < 0.8 * data.temp_beads.size())
                    {
                        cout<<" 0 "<<endl; // 0 means no stalk, 1 means stalk
                    }else{
                        cout<<" 1 "<<endl;;
                    }

                    //
                    //Cluster analysis - distances
                    //
                    vector<double> distances;
                    if( data.in.dist ){

                        int init_head1=select_init(data, Nneigh, 1);
                        clusterize_mtag(data, ves1, Nneigh, neigh_index, init_head1, 1);

                        int init_head2=select_init2(data, Nneigh, 3, ves1);
                        clusterize_mtag(data, ves2, Nneigh, neigh_index, init_head2, 3);
                        get_distances(data, frame, ves1, ves2, distances, 50);
                        order(distances);
                        int sum=0;
                        int avg_dist=-1;
                        int num=100;

                        //
                        // Distance average
                        //
                        if(distances.size()!=0)
                        {
                            for (int n=0;n<num;++n) {
                                sum+=distances[n];
                            }
                            avg_dist=sum/num;
                        }else {
                            avg_dist=0;
                        }

                        ofstream file;
                        file.open ("distances");

                        file << step<< " " << avg_dist <<endl;

                        file.close();
                    }

                } // end for freme while cycle
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

    /**
     * @brief init_grid
     * @param grid
     */
    void init_grid(vector<vector<vector< tuple<int,bool> >>>& grid)
    {
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
    }

    void get_distances(Data& data,rvec* frame,vector <int>& ves1, vector <int>& ves2, vector <double>& distances, int coff){

        Atom cent = com(frame, data, vector<int>{1,3}); //center of mass between ves1 y ves2
        double min, current;
        for(int i=0; i<ves2.size(); ++i) // ~10 000
        {

            Atom ves_2;
            ves_2.x = frame[ves2[i]][0];
            ves_2.y = frame[ves2[i]][1];
            ves_2.z = frame[ves2[i]][2];
            ves_2.type=data.temp_beads[i].type;

            if ((ves_2.type == 1 || ves_2.type==3) && ves_2.distSQ(cent)< coff){

                min = 99999.0;
                for (int n=0;n<ves1.size();n++)     //calculate distance of heads
                {
                    Atom ves;
                    ves.x = frame[ves1[n]][0];
                    ves.y = frame[ves1[n]][1];
                    ves.z = frame[ves1[n]][2];
                    ves.type=data.temp_beads[i].type;


                    if((ves.type == 1 || ves.type==3) && ves.distSQ(cent)< coff){
                        current = ves_2.distSQ(ves);
                        if(current < min)
                            min = current;
                    }
                }
                distances.push_back(min);
            }
        }
    }

    void clusterize(Data& data, vector <int>& cluster, vector <int>& Nneigh, vector <int>& neigh_index, int init, vector<int> tp)
    {
        vector <bool> used(data.temp_beads.size(), false);
        vector <bool> added(data.temp_beads.size(), false);
        int end;

        cluster.push_back(init);

        for (int j=0;j<cluster.size();++j)
        {
            //seguir solo si j no ha sido comprobado antes(if not visited)
            for(int& type : tp)
            {
                if(data.temp_beads[cluster[j]].type == type && used[cluster[j]]==false){
                    used[cluster[j]]=true;
                    int begin=0;
                    for (int a=0;a<cluster[j];++a) {
                        begin=begin+Nneigh[a];
                    }

                    end=begin+Nneigh[cluster[j]]-1;

                    for (int n=begin;n<=end;++n)
                    {
                        //cout<<Nneigh[cluster[j]]<<" "<<neigh_index[n]<<endl;

                        //if neigh is not already in cluster, then add
                        if(added[neigh_index[n]]==false){
                            added[neigh_index[n]]=true;
                            cluster.push_back(neigh_index[n]);
                        }
                    }
                }
            }
        }
    }
    void clusterize_mtag(Data& data, vector <int>& cluster, vector <int>& Nneigh, vector <int>& neigh_index, int init, int moltag)
    {
        vector <bool> used(data.temp_beads.size(), false);
        vector <bool> added(data.temp_beads.size(), false);
        int end;

        cluster.push_back(init);

        for (int j=0;j<cluster.size();++j)
        {
            //seguir solo si j no ha sido comprobado antes(if not visited)
                if(data.temp_beads[cluster[j]].mol_tag == moltag && used[cluster[j]]==false){
                    used[cluster[j]]=true;
                    int begin=0;
                    for (int a=0;a<cluster[j];++a) {
                        begin=begin+Nneigh[a];
                    }

                    end=begin+Nneigh[cluster[j]]-1;

                    for (int n=begin;n<=end;++n)
                    {
                        //cout<<Nneigh[cluster[j]]<<" "<<neigh_index[n]<<endl;

                        //if neigh is not already in cluster, then add
                        if(added[neigh_index[n]]==false){
                            added[neigh_index[n]]=true;
                            cluster.push_back(neigh_index[n]);
                        }
                    }
                }

        }
    }

    void generate_pairlist(Data& data, rvec* frame, vector <int>& Nneigh, vector <int>& neigh_index)
    {
        Atom at, particle;
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

                if(at.distSQ(particle) < pow(2.0, 1.0/3.0) +1.6 ){ // cos^2 attractive potential of Deserno lipid
                    neigh++;
                    neigh_index.push_back(j);
                }
            }
            Nneigh.push_back(neigh);
        }
    }

    void print_grid(vector<vector<vector< tuple<int,bool> >>>& grid)
    {
        for (int n=0;n<grid_size;++n){
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
        }
    }

    int build_grid(Atom particle, vector<vector<vector< tuple<int,bool> >>>& grid)
    {
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

                        if ( grid_point.distSQ(particle) < pow(2.0, 1.0/3.0)+pow(bin_size, 2.0) )
                        {
                            ++get<0>(grid[pointx][pointy][pointz]);
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
                    return false;
                }
            }
        }

        if(get<1>(grid[start2.x][start2.y][start2.z])==true){
            return true;
        }
        return false;
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

    void order(vector<double>& arr)
    {
        int n = arr.size();
        //sorting - ASCENDING ORDER
        for(int i=0;i<n;i++)
        {
            for(int j=i+1;j<n;j++)
            {
                if(arr[i]>arr[j])
                {
                    double temp  =arr[i];
                    arr[i]=arr[j];
                    arr[j]=temp;
                }
            }
        }
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




    int select_init(Data& data, vector <int>& Nneigh, int moltag){
        int init;
        for(int p=0; p<data.temp_beads.size(); ++p) // ~10 000
        {
            Atom part;
            part.mol_tag = data.temp_beads[p].mol_tag;
            if(Nneigh[p]>3 && part.mol_tag == moltag){
                init=p;
                break;
            }
        }
        return init;
    }

    int select_init2(Data& data, vector <int>& Nneigh, int moltag, vector <int>& ves1){
        int init;
        for(int p=0; p<data.temp_beads.size(); ++p) // ~10 000
        {
            bool isves1=false;

            for(int i=0; i<ves1.size(); ++i)    // check if i belongs to ves1, we want it to not be true
            {
                if(p==ves1[i]){
                    isves1=true;
                }
            }
            if(isves1==false){
                Atom part;
                part.mol_tag = data.temp_beads[p].mol_tag;
                if(Nneigh[p]>3 && part.mol_tag == moltag)
                {
                    init=p;
                    break;
                }
            }
        }
        return init;
    }

};

#endif // XTCANALYSIS_H
