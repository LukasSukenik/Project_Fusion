#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <cstdlib>
#include <algorithm>
#include <random>

using namespace std;

class Input
{
public:
    Input(){}

    int op=-1;
    double radius;
    double trim;
    int out_type;
    int num_lipids;
    int multiple;
    string in_file;
    string out_file;
    int mol_tag;
    int offset;
    int num_rec;

    double com_x;
    double com_y;
    double com_z;

    double box_xm;
    double box_ym;
    double box_zm;
    double box_xp;
    double box_yp;
    double box_zp;

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
            if( what.compare("Operation_type:") == 0 ) { ss >> op; }
            if( what.compare("Radius:") == 0 )         { ss >> radius; }
            if( what.compare("Trim:") == 0 )           { ss >> trim; }
            if( what.compare("Out_type:") == 0 )       { ss >> out_type; }
            if( what.compare("Num_lipids:") == 0 )     { ss >> num_lipids; }
            if( what.compare("Multiple:") == 0 )       { ss >> multiple; }
            if( what.compare("In_file:") == 0 )        { ss >> in_file; }
            if( what.compare("Out_file:") == 0 )       { ss >> out_file; }
            if( what.compare("Mol_tag:") == 0 )        { ss >> mol_tag; }

            if( what.compare("Lammps_offset:") == 0 )   { ss >> offset; }
            if( what.compare("Number_of_receptors:") == 0 ) { ss >> num_rec; }
            if( what.compare("Box:") == 0 )             { ss >> box_xm >> box_xp >> box_ym >> box_yp >> box_zm >> box_zp; }
            if( what.compare("Position_shift:") == 0 )  { ss >> com_x >> com_y >> com_z; }

            what.clear();
        }
        fs.close();

        return true;
    }

    void helpMessage()
    {
        cout << "Run binary with filenames that contain the following parameters:" << endl;

        cout << "Operation_type: 1 - trim and multiple,  (param1 = trim size, param2 = multiply)" << endl;
        cout << "                2 - generate vesicle,   (param1 = vesicle radius)" << endl;
        cout << "                3 - copy in z dir,      (param1 = z dir distance)" << endl;
        cout << "                4 - populate receptors, (param1 = fraction of receptors, param2=receptor_lenght)" << endl;
        cout << "                5 - Analyze trajectory, (param1 = dummy, param2=dummy)" << endl;
        cout << "                5 - Analyze trajectory, (input_file = lammps_input_data, output_file=xtc_file)" << endl;
        cout << "                6 - center membrane around 0 in xy plane, rest dummy params" << endl;

        cout << "Radius: xxx - floating point, radius of generated vesicle" << endl;
        cout << "Trim: xxx - floating point, trim edges of input membrane" << endl;

        cout << "Out_type: 1 - Lammps" << endl;
        cout << "          2 - SC" << endl;
        cout << "          other - do nothing" << endl;

        cout << "Num_lipids: XXX - number of lipids in vesicle outer-leaflet" << endl;
        cout << "Multiple: XXX - number of copies of input membrane" << endl;

        cout << "In_file: xxx" << endl;
        cout << "Out_file: xxx" << endl;

        cout << "Mol_tag: XXX - numeric label for generated molecule" << endl;
        cout << "Lammps_offset: XXX " << endl;
        cout << "Number_of_receptors: XXX " << endl;
        cout << "Box: XX XX XX XX XX XX" << endl;
        cout << "Position_shift: XXX XXX XXX" << endl;

    }
};

#endif // INPUT_H
