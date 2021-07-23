#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <cstdlib>
#include <algorithm>
#include <random>
#include <sstream>
#include <string.h>

#include "atom.h"


using namespace std;

class BeadParam {
public:
    BeadParam(){}
    BeadParam(int type, double epsilon, double sigma, double cutoff) : type(type), epsilon(epsilon), sigma(sigma), cutoff(cutoff) {}

    int type=-1;
    double epsilon=1.0;
    double sigma=1.0;
    double cutoff=25.0;

    string toString()
    {
        stringstream ss;
        ss << type << " " << epsilon << " " << sigma << " " << cutoff << endl;
        return ss.str();
    }
};

class CosParam {
public:
    CosParam(){}
    CosParam(int type1, int type2, double epsilon, double start_dis, double range) : type1(type1), type2(type2), epsilon(epsilon), start_dis(start_dis), range(range){}

    int type1;
    int type2;
    double epsilon;
    double start_dis;
    double range;

    string toString()
    {
        stringstream ss;
        ss << type1 << " " << type2 << " " << epsilon << " " << start_dis << " " << range << endl;
        return ss.str();
    }
};

/**
 * @brief The Input class - Parameters for particle generation
 */
class Input{
public:
    Input() {}

    int nano_type;
    int out_type;
    int num_of_beads;
    myFloat scale;
    int offset = 1;
    int seed=0;
    myFloat c;
    int num_lig;
    int mol_tag=-1;
    int atom_type=1;
    bool center=false;
    int mtag_1=-1;
    int mtag_2=-1;
    Atom ivx = Atom(0.0, 0.0, 0.0);
    bool fit = false;
    Atom boxm = Atom(-1,-1,-1);
    Atom boxp = Atom(-1,-1,-1);
    Atom com_pos = Atom(0.0, 0.0, 0.0);
    Atom patch_1 = Atom(1,1,1,0);
    Atom patch_2 = Atom(1,1,1,0);
    string infile;

    vector<BeadParam> bparam;
    vector<CosParam> cparam;
    vector<int> types;

    string toString()
    {
        stringstream ss;

        ss << "Generating type: " << nano_type << endl;
        ss << "Out type: " << out_type << endl;
        ss << "Number of beads: " << num_of_beads << endl;
        ss << "Scale: " << scale << endl;
        ss << "Offset: " << offset << endl;
        ss << "c: " << c << endl;
        ss << "Number of ligands: " << num_lig << endl;
        ss << "Box: (" << boxm.x << ", " << boxp.x << ", " << boxm.y << ", " << boxp.y << ", " << boxm.z << ", " << boxp.z << ")" << endl;
        ss << "Position: (" << com_pos.x << ", " << com_pos.y << ", " << com_pos.z << ")" << endl;
        ss << "Patch_1: (" << patch_1.x << "-" << patch_1.vx << ", " << patch_1.y << "-" << patch_1.vy << ", " << patch_1.z << "-" << patch_1.vz << ", " << patch_1.type << ")" << endl;
        ss << "Patch_1: (" << patch_2.x << ", " << patch_2.y << ", " << patch_2.z << ", " << patch_2.type << ")" << endl;

        ss << "Beads_lj/cut:" << endl;
        for(auto i : bparam)
        {
            ss << i.toString();
        }

        ss << "Cosatt:" << endl;
        for(auto i : cparam)
        {
            ss << i.toString();
        }

        return ss.str();
    }

    bool loadInput(string input)
    {
        std::fstream fs( input, std::fstream::in );
        string line, what;
        stringstream ss;
        int len=0;

        while( !fs.eof() ) // Lines in input
        {
            ss.flush();
            ss.clear();
            getline(fs, line);
            ss.str(line);

            ss >> what;
            if( what.compare("Particle_type:") == 0 )   { ss >> nano_type; }
            if( what.compare("Output_type:") == 0 )     { ss >> out_type; }
            if( what.compare("Num_of_beads:") == 0 )    { ss >> num_of_beads; }
            if( what.compare("Scale:") == 0 )           { ss >> scale; }
            if( what.compare("Lammps_offset:") == 0 )   { ss >> offset; }
            if( what.compare("c:") == 0 )               { ss >> c; }
            if( what.compare("Number_of_ligands:") == 0 ) { ss >> num_lig; }
            if( what.compare("Box:") == 0 )             { ss >> boxm.x >> boxp.x >> boxm.y >> boxp.y >> boxm.z >> boxp.z; }
            if( what.compare("Position_shift:") == 0 )  { ss >> com_pos.x >> com_pos.y >> com_pos.z; }
            if( what.compare("Load_file:") == 0 )  { ss >> infile; }
            if( what.compare("Center") == 0 )  { center=true; }
            if( what.compare("Fit") == 0 )  { fit=true; }
            if( what.compare("Align:") == 0 )  { ss >> mtag_1 >> mtag_2;  }
            if( what.compare("Impact_vector:") == 0 )  { ss >> ivx.x >> ivx.y >> ivx.z; }
            if( what.compare("Patch_1:") == 0 )  { ss >> patch_1.vx >> patch_1.x >> patch_1.vy >> patch_1.y >> patch_1.vz >> patch_1.z >> patch_1.type; }
            if( what.compare("Patch_2:") == 0 )  { ss >> patch_2.vx >> patch_2.x >> patch_2.vy >> patch_2.y >> patch_2.vz >> patch_2.z >> patch_2.type; }
            if( what.compare("Seed:") == 0 )  { ss >> seed; rng.seed(seed); }

            if( what.compare("Mol_tag:") == 0 ) { ss >> mol_tag; }
            if( what.compare("Atom_type:") == 0 ) { ss >> atom_type; }

            if( what.compare("Beads_lj/cut:") == 0 )
            {
                ss >> len;
                for(int i=0; i < len; ++i)
                {
                    ss.flush();
                    ss.clear();
                    getline(fs, line);
                    ss.str(line);
                    bparam.push_back(BeadParam());
                    ss >> bparam.back().type >> bparam.back().epsilon >> bparam.back().sigma >> bparam.back().cutoff;
                }
                len=0;
            }

            if( what.compare("Cosatt:") == 0 )
            {
                ss >> len;
                for(int i=0; i < len; ++i)
                {
                    ss.flush();
                    ss.clear();
                    getline(fs, line);
                    ss.str(line);
                    cparam.push_back(CosParam());
                    ss >> cparam.back().type1 >> cparam.back().type2 >> cparam.back().epsilon >> cparam.back().start_dis >> cparam.back().range;
                }
            }

            what.clear();
        }
        fs.close();

        for(auto i : bparam)
        {
            types.push_back(i.type);
        }
        return true;
    }

    void clear()
    {
        nano_type=-1;
        out_type=-1;
        num_of_beads=-1;
        scale=0.0;
        offset=1;
        c=0;
        num_lig=0;
        mol_tag=-1;
        atom_type=1;

        center=false;
        fit = false;
        mtag_1=-1;
        mtag_2=-1;

        boxm=Atom(-1,-1,-1);
        boxp=Atom(-1,-1,-1);
        com_pos=Atom(0.0, 0.0, 0.0);
        patch_1=Atom(1,1,1,0);
        patch_2=Atom(1,1,1,0);
        ivx=Atom(0.0, 0.0, 0.0);

        infile.clear();
        bparam.clear();
        cparam.clear();
        types.clear();
    }
};

class My_string{
public:
    My_string(char str[256] ) {
        strcpy(this->str,str);
    }
    char str[256];
};

bool sortN(const Atom& i, const Atom& j) {
    return i.N < j.N;
}

/**
 * @brief The Data class - Generated data from class particle
 */
class Data
{
public:
    /// lammps data stuff
    vector< Atom > all_beads;
    vector< Bond > all_bonds;
    vector<Angle> all_angles;
    vector<BeadParam> all_bparam;
    vector<CosParam> all_cparam;

    vector< Atom > temp_beads;
    vector< Bond > temp_bonds;
    vector<Angle> temp_angles;

    /// force field stuff
    array<array<double, 100>, 100> all_sigma;
    int all_sigma_size = 0;
    array< array<bool, 100>, 100> all_sigma_cosatt;

    /// For lammps file IO
    vector<My_string > file_head;
    double box[6] = {0};

    bool first_file=true;

    myFloat x_min = 0.0;
    myFloat x_max = 0.0;
    myFloat y_min = 0.0;
    myFloat y_max = 0.0;
    myFloat z_min = 0.0;
    myFloat z_max = 0.0;

    Input in;

    Data()
    {
        for(int i=0; i<99; ++i)
            for(int j=0; j<99; ++j)
                all_sigma_cosatt[i][j] = false;
        for(int i=0; i<99; ++i)
            for(int j=0; j<99; ++j)
                all_sigma[i][j] = 0.0;
    }

    bool isDefined()
    {
        if(!in.infile.empty())
            return true;
        return false;
    }

    int getMaxMolTag()
    {
        int max=0;
        for(auto& a : all_beads)
        {
            if(a.mol_tag > max)
            {
                max = a.mol_tag;
            }
        }
        return max;
    }

    void move(Atom move)
    {
        for(Atom& item : temp_beads)
            item += move;
        cerr << "move " << move.x << " " << move.y << " " << move.z << " done" << endl;
    }

    /**
     * @brief center_of_mass - function computes Center-Of-Mass (COM) of particles with a givem mol_tag
     * @param mtag - mol_tag of particles for COM calculation, -1 = all particles regarles of mol_tag
     */
    Atom center_of_mass(int mtag=-1, int start=-1, int stop=-1)
    {
        int count=0;
        int total=0;
        Atom cm;
        for(Atom& item : temp_beads)
        {
            if( (item.mol_tag == mtag || mtag == -1) && total >= start && (total < stop || stop == -1) )
            {
                cm += item;
                ++count;
            }

            if(item.mol_tag == mtag || mtag == -1)
                ++total;
        }
        cm *= 1.0/count;
        return cm;
    }

    /**
     * @brief impact - Deprecated, don't use
     * @param ivx
     */
    void impact(Atom ivx)
    {
        if(ivx.size() > 0.1)
        {
            Atom impact = Atom(0.0, 0.0, 0.0);

            // move the liposome to COM
            Atom com = center_of_mass();
            com *= -1;
            move(com);

            // Define nanoparticle and liposomes dimensions
            double z_min_lipo3 = 9999; // liposome we are adding
            double y_min_lipo3 = 9999; // liposome we are adding
            double y_max_lipo3 = -9999; // liposome we are adding
            double z_max_nano = -9999;  // nanoparticle
            double z_max_lipo1 = -9999; // bound liposome
            double y_max_lipo1 = -9999; // bound liposome
            double y_min_lipo1 = 9999; // bound liposome


            for(Atom& item : temp_beads)
            {
                if(z_min_lipo3 > item.z)
                    z_min_lipo3 = item.z;
                if(y_min_lipo3 > item.y)
                    y_min_lipo3 = item.y;
                if(y_max_lipo3 < item.y)
                    y_max_lipo3 = item.y;
            }

            for(Atom& item : all_beads)
            {
                if(y_min_lipo1 > item.y && item.mol_tag == 1)
                    y_min_lipo1 = item.y;
                if(y_max_lipo1 < item.y && item.mol_tag == 1)
                    y_max_lipo1 = item.y;
                if(z_max_lipo1 < item.z && item.mol_tag == 1)
                    z_max_lipo1 = item.z;
                if(z_max_nano < item.z && item.mol_tag == 2)
                    z_max_nano = item.z;
            }

            // z impact
            //z_impact.z = z_max_nano - z_min_lipo3 + 3;

            // y impact
            //y_impact.y = - y_max_lipo3 - z_max_nano -1;
            //y_impact.z = z_max_lipo1 - z_min_lipo3 + 0.5;

            ivx.normalise();
            impact.y = (- y_max_lipo3 - z_max_nano -1) * ivx.y;
            impact.z = (z_max_lipo1 - z_min_lipo3 + 0.5)*(1-ivx.z) + ivx.z*(z_max_nano - z_min_lipo3 + 3);
            move(impact);

            cerr << "Liposome moved by " << impact.x << " " << impact.y << " " << impact.z << endl;
        }
    }

    /**
     * @brief fit - positions the loaded/generated structure next to previosly generated/loadedd structure
     * - used for ideal collision position of two liposomes and a nanoparticle
     * see Fit_function.blend, need blender 2.9
     */
    void fit()
    {
        //
        // TODO: Construct an algorithm for positioning the second liposome in an ideal collision position.
        // - I think something as shown in Fit_function.blend will work fine, but feel free to innovate
        //
        // Moving the structure is already implemented in move(Atom vec) function, example below
        // Calculating overlap is simple as well, example below
        //

        Atom displace = Atom(0, 0, 23); // class Atom works as a vector as well.
        move(displace); // displace liposome2 by vector displace
        return;

        //for(Atom& item : temp_beads) // loop over liposome2
        //for(Atom& item : all_beads) // loop over liposome+nanoparticle structure
        double distance_squared; // we are using distance squared because it takes less resources -> we are not calculating the square root
        double too_small = 1.4; // maybe bigger, smaller? Best to eyeball it for vmd once you make it semi-functional
        for(Atom& lip2 : temp_beads)
        {
            for(Atom& nano_lip : all_beads)
            {
                if(nano_lip.mol_tag == 2) // 2 is nanoparticle in our use case, 1 is the liposome
                {
                    distance_squared = lip2.distSQ(nano_lip); // assume overlap if distance squared too small
                    if(distance_squared < too_small)
                    {
                        cerr << "Overlap!" << endl;
                    }
                }
            }
        }

        //
        // Rotating a structure is shown in align function
        //

        //
        // Finally Calculate impact vector and print it out into a file, then load in in prep.sh script
        // tutorial for input/output in c++ https://www.cplusplus.com/doc/tutorial/files/
        //
    }

    /**
     * @brief align - align a liposome and nanoparticle
     * @param mtag - nanoparticle mol_tag
     * @param mtag2 - liposome mol_tag
     */
    void align(int mtag, int mtag2)
    {
        // Test for empty mol_tags
        if( mtag != -1 && mtag2 != -1 )
        {
            Atom x_axis = Atom(1,0,0);
            Atom x_axis_negative = Atom(-1,0,0);
            Atom z_axis = Atom(0,0,-1);

            //
            // Move nanoparticle to center (0,0,0)
            //
            center(mtag);

            //
            // Rotate structure so that mtag beads (nanoparticle) align with x_axis
            // - nanoparticle generated from poles = tips in prolate form, same as in oblate form
            // -- 1/4 beads from each end identify the poles (tips)
            //
            int count = countMoltag(mtag, temp_beads);             // number of mtag (nanoparticle beads)
            Atom nano1 = center_of_mass(mtag, 0, count/4);         // first 1/4 COM of mtag beads
            Atom nano2 = center_of_mass(mtag, 1+3*count/4, count); // last 1/4 COM of mtag beads
            Atom nano_axis = nano1-nano2;                          // Axis of mtag beads
            nano_axis.normalise();                                 // normalise axis vector for correct rotation
            //
            // Look at vector_magic.blend for visual example, need blender 2.9
            // - rot_axis defined plane is between the nano_axis and x_axis vector
            // - by rotating the nano_axis vector 180deg in this plane we align in at x_axis_negative exactly
            //
            Atom rot_axis = nano_axis-x_axis;
            rot_axis.normalise();                                  // normalise for rotation algo
            for(Atom& item : temp_beads)
            {
                item.rotate(rot_axis, 3.14159265359);       // rotate 180deg = 3.1415 radians, rad to deg = 57.2958
            }

            //
            // Rotate mtag so that COM of mtag2 is (*,*,0) = centered around z axis
            // - to keep it aligned with x axis, we rotate only around x axis
            //
            Atom com_mtag2 = center_of_mass(mtag2);
            com_mtag2.x = 0.0;
            com_mtag2.normalise();
            double angle = acos( com_mtag2.dot(z_axis) );
            double clockwise = (com_mtag2.cross(z_axis)).x ;

            if(clockwise > 0.0)
            {
                for(Atom& item : temp_beads)
                {
                    item.rotate(x_axis, angle);       //
                }
            } else
            {
                for(Atom& item : temp_beads)
                {
                    item.rotate(x_axis_negative, angle);       //
                }
            }

            //
            // TODO: Construct nanoparticle patch, rotate structure so patch points to in +y axis
            //
            /*Atom patch_vec;
            Atom rotate_axis;
            double rotate_angle;*/

            //
            // Class Atom has variable x,y,z. You can access them via . (dot)
            // example: patch_vec.y
            //
            // there are a number of function within class Atom that you can use. -, +, *, /, dot, cross.
            // If you are unsure what they do look at how they are programmed in class Atom
            //

            //
            // function center_of_mass(mol_tag) returns position of center of mass. This is stored in class Atom
            //
            //Atom mtag2_COM = center_of_mass(mtag2);

            //
            // Rotates structure around rotate_axis by angle rotate_angle
            //
            /*rotate_axis.normalise();
            for(Atom& item : temp_beads)
            {
                item.rotate(rotate_axis, rotate_angle);       //
            }*/

            cerr << "Aligned to x axis and z axis" << endl;
        }
    }

    void add()
    {
        all_beads.insert(all_beads.end(), temp_beads.begin(), temp_beads.end());
        all_bonds.insert(all_bonds.end(), temp_bonds.begin(), temp_bonds.end());
        all_angles.insert(all_angles.end(), temp_angles.begin(), temp_angles.end());

        temp_beads.clear();
        temp_bonds.clear();
        temp_angles.clear();
    }

    void offset(int offs)
    {
        for(Atom& item : temp_beads)
        {
            item.N += offs;
        }

        for(Bond& item : temp_bonds)
        {
            item.N += offs;
            item.at1 += offs;
            item.at2 += offs;
        }

        for(Angle& item : temp_angles)
        {
            item.N += offs;
            item.at1 += offs;
            item.at2 += offs;
            item.at3 += offs;
        }
    }

    void mol_tag(int mtag)
    {
        if(in.mol_tag != -1)
        {
            for(Atom& item : temp_beads)
            {
                item.mol_tag = mtag;
            }
            cerr << "All beads changed to mol_tag = " << mtag << endl;
        }
    }

    void printAllSigma()
    {
        cerr << "printAllSigma:" << endl;
        for(unsigned int i=0; i<all_sigma_size; ++i)
        {
            cerr << "[" << i+1 << "][" << i+1 << "] = "  << all_sigma[i][i] << endl;
        }
    }


    /**
     * @brief center - moves structure, so that center_of_mass is (0,0,0)
     * @param mtag
     */
    void center(int mtag=-1)
    {
        if(! temp_beads.empty())
        {
            Atom cm = center_of_mass(mtag);
            cm*=-1.0;
            move( cm );
            cerr << "center of " << mtag << " done" << endl;
        }
        else
        {
            cerr << "Error: No atoms loaded" << endl;
        }
    }

    void rescale(double rescale)
    {
        for(Atom& item : temp_beads)
            item *= rescale;
        cerr << "rescale " << rescale << " done" << endl;
    }

    bool loadInput(string input)
    {
        in.clear();
        in.loadInput(input);
        return true;
    }

    bool isOverlap(Atom& a)
    {
        for(Atom& item : all_beads)
        {
            if( item.dist(a) > getBeadParam(item.type).cutoff * 2)
            {
                return true;
            }
        }
        return false;
    }

    BeadParam getBeadParam(int type)
    {
        for(auto item : all_bparam)
        {
            if( type == item.type)
            {
                return item;
            }
        }
        return BeadParam();
    }

    string toString()
    {
        stringstream ss;
        vector<int> moltags = getMolTypes();
        vector<int> types = getAtomTypes();

        ss << "Beads: " << all_beads.size() << endl;
        ss << types.size() << " atom types:" << endl;
        for(int atp=0; atp< types.size(); ++atp)
        {
            ss << "Atom type " << types[atp] << " " << countAtomType(types[atp]) << endl;
        }

        ss << moltags.size() << " molTypes:" << endl;
        for(int moltg=0; moltg< moltags.size(); ++moltg)
        {
            ss << "Molecule " << moltags[moltg] << " " << countMoltag(moltags[moltg], all_beads) << endl;
        }

        return ss.str();
    }

    int countAtomType( int atype)
    {
        int count = 0;
        for(Atom& item : all_beads)
        {
            if( item.type == atype)
            {
                ++count;
            }
        }
        return count;
    }

    int countMoltag(int mTag, vector< Atom >& container)
    {
        int count = 0;
        for(Atom& item : container)
        {
            if( item.mol_tag == mTag)
            {
                ++count;
            }
        }
        return count;
    }

    vector<int> getMolTypes()
    {
        vector<int> moltags;
        bool exist = false;
        for(Atom& item : all_beads)
        {
            exist = false;
            for(auto mtype : moltags)
            {
                if( item.mol_tag == mtype )
                {
                    exist = true;
                }
            }
            if(!exist)
            {
                moltags.push_back(item.mol_tag);
            }
        }
        return moltags;
    }


    /// Loading Lammps file - taylored to membrane ///
    void load(string in_file);
    void loadFileHeadAndPart(string filename);
    void loadBox();
    void loadBonds(string filename);
    void loadAngles(string filename)
    {
        char str[256];
        Angle angle;
        std::fstream in;

        in.open(filename, std::fstream::in);
        if(in.good())
        {
            while(in.good())
            {
                in.getline(str, 256);
                if(strstr(str, "Angles") != NULL)
                    break;
            }

            while(in.good())
            {
                in >> angle.N >> angle.type >> angle.at1 >> angle.at2 >> angle.at3;
                if( temp_angles.empty() || angle.N != temp_angles.back().N)
                {
                    temp_angles.push_back(angle);
                }
            }
        }
        in.close();
    }
    //////////////////////////////////////////////////

    vector<int> getAtomTypes() const
    {
        vector<int> atom_types;
        bool exist = false;

        atom_types.push_back(all_beads[0].type);
        for(int i=0; i<all_beads.size(); ++i) {
            exist = false;
            for(int j=0; j<atom_types.size(); ++j) {
                if(atom_types[j] == all_beads[i].type)
                    exist = true;
            }
            if(!exist)
                atom_types.push_back( all_beads[i].type );
        }
        return atom_types;
    }

    /**
     * @brief calcBondTypes - count all bond types
     * @return
     */
    int calcBondTypes() const
    {
        vector<int> bond_types;
        bool exist = false;

        if(!all_bonds.empty())
        {
            bond_types.push_back(all_bonds[0].type);
            for(auto& bond : all_bonds) {

                exist = false;
                for(auto btype : bond_types) {
                    if(btype == bond.type) {
                        exist = true;
                    }
                }
                if(!exist) {
                    bond_types.push_back( bond.type );
                }
            }
        }

        return bond_types.size();
    }



    void printForceField(vector<double>& dist, vector<string> &dist_coeff, double scale ) const
    {
        fstream force_field("force_field", fstream::out);

        force_field << "variable coeff_1 string 1.0" << endl;
        force_field << "variable coeff_2 string 1.0" << endl;
        force_field << "variable coeff_3 string 1.0" << endl;
        force_field << "variable coeff_bond string 50" << endl;
        force_field << endl;

        for(int i=0; i<all_sigma_size; ++i) {
            for(int j=i; j<all_sigma_size; ++j) {
                force_field << "pair_coeff " << i+1 << " " << j+1 << " lj/cut 1.0 " << 0.5 * (all_sigma[i][j]+all_sigma[i][j]) / 1.122462048 << " " << 0.5 * (all_sigma[i][j]+all_sigma[i][j]) << endl;
            }
            force_field << endl;
        }

        force_field << endl;
        int count = 1;
        for(int i=0; i<all_sigma_size; ++i) {
            for(int j=i; j<all_sigma_size; ++j) {
                if(all_sigma_cosatt[i][j]) {
                    force_field << "pair_coeff " << i+1 << " " << j+1 << " cosatt ${coeff_" << count << "} " << 0.5 * (all_sigma[i][j]+all_sigma[i][j]) << " 1.0"<< endl; // 2 3
                    ++count;
                }
            }
        }

        force_field << "\n" << endl;

        for(int j=0; j<dist.size(); ++j) {
            force_field << "bond_coeff " << j+1 << " harmonic ${" << dist_coeff[j] << "} " << dist[j] << endl;
        }

        force_field.close();
    }

    void roundBondDist(vector<double>& dist, double precision=1000.0)
    {
        for(int j=0; j<dist.size(); ++j) // each j+1 is bond type
        {
            for(auto& bond : all_bonds)
            {
                if( !bond.typelock && isAproxSame( bond.r0, dist[j], 1.0/precision) )
                {
                    bond.type = j+1;
                }
            }
        }
    }

    void removeDuplicateBond()
    {
        int count = 0;
        for(auto& a : all_bonds)
        {
            for(auto& b : all_bonds)
            {
                if(a.at1 == b.at1 && a.at2 == b.at2 && a.N != b.N)
                {
                    cerr << "Duplicate bond " << a.N << "==" << b.N << endl;
                    ++count;
                }
            }
        }
        cerr << count << endl;
    }

    vector<double> createBondGroups(vector<string> &bond_coeff, double precision=1000.0) const
    {
        vector<double> dist;
        if(!all_bonds.empty())
        {
            bool exist;
            for(auto& bond : all_bonds)
            {
                if( !bond.typelock )
                {
                    exist = false;
                    for(int j=0; j<dist.size(); ++j)
                    {
                        if( isAproxSame( bond.r0, dist[j], 1.0/precision) )
                        {
                            exist=true;
                        }
                    }
                    if(!exist)
                    {
                        dist.push_back( round(precision * bond.r0) / precision );
                        bond_coeff.push_back( bond.coeff_name );
                    }
                }
            }
        }
        return dist;
    }

    void printLammps() const
    {
        if(all_beads.empty()) {
            cerr << "No beads generated" << endl;
            return;
        }

        vector<double> dist;
        vector<string> dist_coeff;

        dist = createBondGroups(dist_coeff);
        printForceField(dist, dist_coeff, in.scale);

        int bond_types = calcBondTypes();
        int num_a_types = getAtomTypes().size();

        //
        // Print Head
        //
        cout << "LAMMPS data file via gen_membrane\n" << endl;
        cout << all_beads.size() << " atoms\n";
        cout << num_a_types << " atom types\n";

        //
        // Print Bond types
        //
        if(!all_bonds.empty()) {
            cout << all_bonds.size() << " bonds\n";
            cout << bond_types << " bond types\n";
        }
        //
        // Print Angle types
        //
        if(!all_angles.empty()) {
            cout << all_angles.size() << " angles\n";
            cout << "1 angle types\n";
        }

        //
        // Print Box
        //
        cout << "\n";
        if(in.boxm.x != -1) {
            cout << in.boxm.x << " " << in.boxp.x << " xlo xhi\n";
            cout << in.boxm.y << " " << in.boxp.y << " ylo yhi\n";
            cout << in.boxm.z << " " << in.boxp.z << " zlo zhi\n";
        } else {
            cout << ( box[0]-((x_min+x_max)/2.0) ) << " " << ( box[1]-((x_min+x_max)/2.0) ) << " xlo xhi\n";
            cout << ( box[2]-((x_min+x_max)/2.0) ) << " " << ( box[3]-((x_min+x_max)/2.0) ) << " ylo yhi\n";
            cout << box[4] << " " << box[5] << " zlo zhi\n";
        }

        //
        // Print Masses
        //
        cout << "\nMasses\n\n";
        for(int i=0; i<num_a_types; ++i )
            cout << i+1 << " mass_" << i+1 << "\n";

        //
        // Print Atoms
        //
        cout <<"\nAtoms # full\n" << endl;
        for(auto& a : all_beads) {
            cout << a.N << " " << a.mol_tag << " " << a.type << " " << 0 << " " << a << " 0 0 0" << endl;
        }

        //
        // Print Bonds
        //
        if(!all_bonds.empty())
        {
            cout << "\nBonds\n\n";

            for(auto& bond : all_bonds)
            {
                cout << bond.N << " " << bond.type << " " << bond.at1 << " " << bond.at2 << "\n";
            }
        }

        //
        // Print Angles
        //
        if(!all_angles.empty()) {
            cout << "\nAngles\n\n";

            for(auto & a : all_angles) {
                cout << a.N << " " << a.type << " " << a.at1 << " " << a.at2 << " " << a.at3 << "\n";
            }
        }
    }

    void printXYZ() const
    {
        cout << all_beads.size() << "\nparticle\n";
        for(unsigned int i=0; i<all_beads.size(); i++) {
            cout << "C" << all_beads[i].type <<  " " << (all_beads[i]*in.scale) + in.com_pos << endl;
        }
    }
};



void Data::load(string in_file)
{
    if(in_file.empty()) {
        return;
    }

    cerr << "Loading file " << in_file << endl;
    loadFileHeadAndPart(in_file);
    cerr << "Loading box parameters" << endl;
    loadBox();

    // THROW AWAY VELOCITIES

    cerr << "Parsing bonds parameters" << endl;
    loadBonds(in_file);

    // TEST FOR CONSECCUTIVE INDEXES
    cerr << "Sorting lipids by index" << endl;
    std::sort(temp_beads.begin(), temp_beads.end(), sortN);
    cerr << "Testing index duplicity" << endl;
    for(int i=0; i<temp_beads.size(); i++) {
        if(i%1000 == 0)
        {
            cerr << i << " : " << temp_beads.size() << endl;
        }
        if(i+1 != temp_beads[i].N) {
            cerr << "ERROR, missing index, " << i << " != " << temp_beads[i].N << endl;
        }
    }
    cerr << "Load done" << endl;

    loadAngles(in_file);
}

void Data::loadBox()
{
    stringstream str;
    stringstream str2;
    stringstream str3;

    for(unsigned int i=0; i<file_head.size(); i++) {
        //cout << file_head[i].str << endl;
        if(strstr(file_head[i].str, "xlo xhi") != nullptr) {
            str << file_head[i].str;
            str >> box[0] >> box[1];
            continue;
        }

        if(strstr(file_head[i].str, "ylo yhi") != nullptr) {
            str2 << file_head[i].str;
            str2 >> box[2] >> box[3];
            continue;
        }

        if(strstr(file_head[i].str, "zlo zhi") != nullptr) {
            str3 << file_head[i].str;
            str3 >> box[4] >> box[5];
            continue;
        }
    }
}

void Data::loadBonds(string filename)
{
    char str[256];
    Bond bond;
    bond.typelock = true;
    std::fstream in;

    in.open(filename, std::fstream::in);
    if(in.good())
    {
        while(in.good())
        {
            in.getline(str, 256);
            if(strstr(str, "Bonds") != NULL)
                break;
        }

        while(in.good())
        {
            in >> bond.N >> bond.type >> bond.at1 >> bond.at2;
            temp_bonds.push_back(bond);
        }
    }
    in.close();
}

void Data::loadFileHeadAndPart(string filename)
{
    char str[256];
    Atom part;
    std::fstream in;
    std::stringstream ss;
    int size=temp_beads.size();

    in.open(filename, std::fstream::in);

    if(in.is_open())
    {
        // Load until "Atoms" occurs
        while(in.good() && strstr(str, "Atoms") == NULL)
        {
            in.getline(str, 256);
            if(first_file)
            {
                file_head.push_back(My_string(str));
            }
        }
        in.getline(str, 256); // 1 empty line after atoms

        // Load: N molecule-tag atom-type q x y z nx ny nz
        while(in.good())
        {
            part.N = -1;
            part.nx = 0;
            part.ny = 0;
            part.nz = 0;
            in.getline(str, 256);
            ss.str( str );
            ss >> part.N >> part.mol_tag >> part.type >> part.q >> part.x >> part.y >> part.z >> part.nx >> part.ny >> part.nz;
            ss.flush();
            ss.clear();

            if( part.N == -1 ) {
                if(temp_beads.empty())
                    continue;
                else
                    break;
            }
            temp_beads.push_back(part);
        }
        cerr << "  Added " << temp_beads.size() - size << " beads" << endl;
    } else {
        cerr << "File " << filename << " file not opened" << endl;
    }

    first_file = false;
    in.close();
}

#endif // DATA_H
