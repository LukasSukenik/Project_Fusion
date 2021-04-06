#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <random>

using namespace std;

class Vector{
public:
    int type;
    double x,y,z;
    Vector(){}
    Vector(double x, double y, double z, int type=0): x(x), y(y), z(z), type(type) {}

    bool operator==(const Vector& o) {
        if(x != o.x) return false;
        if(y != o.y) return false;
        if(z != o.z) return false;
        /*const double aprox = 0.0000000000001;
        if(x < o.x+aprox && x > o.x - aprox)
            return false;
        if(y < o.y+aprox && y > o.y - aprox)
            return false;
        if(z < o.z+aprox && z > o.z - aprox)
            return false;*/

        return true;
    }

    bool isAproxSame(const Vector& o) {
        const double aprox = 0.000001;
        if(!(x < o.x+aprox && x > o.x - aprox))
            return false;
        if(!(y < o.y+aprox && y > o.y - aprox))
            return false;
        if(!(z < o.z+aprox && z > o.z - aprox))
            return false;

        return true;
    }

    double size() {
        return sqrt(x*x + y*y + z*z);
    }

    Vector operator*(double a) {
        return Vector(x*a,y*a,z*a);
    }

    void operator*=(double a) {
        x*=a;
        y*=a;
        z*=a;
    }

    double dot(Vector& o) {
        return this->x*o.x + this->y*o.y + this->z*o.z;
    }

    Vector operator+(Vector& o) {
        return Vector(x+o.x, y+o.y, z+o.z);
    }

    Vector operator-(Vector& o) {
        return Vector(x-o.x, y-o.y, z-o.z);
    }

    double dist(Vector& o) {
        return sqrt( (this->x-o.x)*(this->x-o.x) + (this->y-o.y)*(this->y-o.y) + (this->z-o.z)*(this->z-o.z) );
    }

    inline Vector cross(const Vector& B) const {
        return Vector(this->y*B.z - this->z*B.y, -this->x*B.z + this->z*B.x, this->x*B.y - this->y*B.x);
    }

    bool isNeighbor(const Vector& o) {
        Vector vec(o.x - x, o.y - y, o.z - z);
        double dist = vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
        return (dist < 1.00000000000001 && dist > 0.9999999999999 );
    }
};

std::ostream& operator<<(std::ostream& os, const Vector& vec) {
  os << vec.x << " " << vec.y << " " << vec.z;
  return os;
}

Vector cut(Vector& a, Vector& b, Vector& c, int size, int i, int j) {
    Vector vecAB( b.x - a.x, b.y - a.y, b.z - a.z); // vector from a to b
    Vector vecBC( c.x - b.x, c.y - b.y, c.z - b.z); // vector from b to c
    size-=1;

    vecAB *= (1.0/size);
    vecBC *= (1.0/size);

    return Vector(a.x + i*vecAB.x + j*vecBC.x, a.y + i*vecAB.y + j*vecBC.y, a.z + i*vecAB.z + j*vecBC.z);
}

bool myfunction (Vector i,Vector j) { return (i.size()<j.size()); }

class Particle {
public:
    const int typeNano = 5;
    const int typeLig = 6;
    vector< Vector> beads;

    virtual void generate(const int size)=0;

    void printXYZ(double scale, Vector& com_pos)
    {
        cout << beads.size() << "\nparticle\n";
        for(unsigned int i=0; i<beads.size(); i++) {
            if(beads[i].type == typeNano) // other
                cout << "C2 " << (beads[i]*scale) + com_pos << endl;
            if(beads[i].type == typeLig) // ligands
                cout << "C1 " << (beads[i]*scale) + com_pos << endl;
        }
    }

    void printLammps(double scale, Vector& com_pos, int offset) {

    }
};

class Icosahedron : public Particle {
public:
    Vector edge[12];

    void printXYZ(double scale, Vector& com_pos);
    void printLammps(double scale, Vector& com_pos, int offset);

    Icosahedron() {
        // save edges
        edge[0] = Vector(0.0,                                   0.0,                               -5/( sqrt(50-10*sqrt(5)) ));
        edge[1] = Vector(0.0,                                   0.0,                                5/( sqrt(50-10*sqrt(5)) ));
        edge[2] = Vector(-sqrt(2/(5-sqrt(5))),                  0.0,                               -1/( sqrt(10-2*sqrt(5)) ));

        edge[3] = Vector(sqrt(2/(5-sqrt(5))),                   0.0,                                1/( sqrt(10-2*sqrt(5)) ));
        edge[4] = Vector((1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),   -0.5,                               -1/( sqrt(10-2*sqrt(5)) ));
        edge[5] = Vector((1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),    0.5,                               -1/( sqrt(10-2*sqrt(5)) ));

        edge[6] = Vector(-(1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),  -0.5,                                1/( sqrt(10-2*sqrt(5)) ));
        edge[7] = Vector(-(1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),   0.5,                                1/( sqrt(10-2*sqrt(5)) ) );
        edge[8] = Vector(-(-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))), -0.5*sqrt((5+sqrt(5))/(5-sqrt(5))), -1/( sqrt(10-2*sqrt(5)) ) );

        edge[9] = Vector(-(-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),  0.5*sqrt((5+sqrt(5))/(5-sqrt(5))), -1/( sqrt(10-2*sqrt(5)) ) );
        edge[10] = Vector((-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))), -0.5*sqrt((5+sqrt(5))/(5-sqrt(5))),  1/( sqrt(10-2*sqrt(5)) ));
        edge[11] = Vector((-1+sqrt(5))/(2*sqrt(10-2*sqrt(5))),  0.5*sqrt((5+sqrt(5))/(5-sqrt(5))),  1/( sqrt(10-2*sqrt(5)) ));
    }

    void generate(const int size) {
        // Identify neighbors -> 20 times 3 neighbors
        bool same = false;
        Vector push;
        int aa=0;
        int bb=0;
        for(int a=0; a<12; a++) {
            for(int b=a; b<12; b++) {
                if(edge[a].isNeighbor(edge[b]) ) {
                    for(int c=b; c<12; c++) {
                        if(edge[b].isNeighbor(edge[c]) && edge[a].isNeighbor(edge[c]) ) {
                            for(int i=0; i<size; i++) {
                                for(int j=0; j<=i; j++) {
                                    // This is executed as:
                                    //                          .       i sets depth as well as number of repetitions
                                    //                          ..
                                    //                          ...
                                    //                          ....
                                    // and so on
                                    same = false;
                                    push = cut(edge[a], edge[b], edge[c], size, i, j);
                                    for(unsigned int q=0; q< beads.size(); q++) {
                                        if(beads[q].isAproxSame(push)) {
                                            same = true;
                                        }
                                    }
                                    if(!same) {

                                        /*if(i%4 == 0 && j%4 ==0) { //  5,9,13,17,21,25, 29,....
                                            push.type = 6;
                                            bb++;
                                        }
                                        else {
                                            push.type = 5;
                                            aa++;
                                        }*/

                                        if(i%3 == 0 && j%3 ==0) { // 11% coverage, 4,7,10,13,16,....
                                            push.type = 6;
                                            bb++;
                                            if( ( j == 0 && j == i ) || ( j == i && (i == size-1) ) || ( (i == size-1) && ( j == 0 ) ) )
                                                push.type = 8;
                                            else if(j == 0 || j == i || (i == size-1))
                                                push.type = 7;
                                        }
                                        else {
                                            push.type = 5;
                                            if(j == 0 || j == i || (i == size-1))
                                                push.type = 9;
                                            aa++;
                                        }

                                        /*if(i%2 == 0 && j%2 ==0) { // 11% coverage, 3,5,7,...
                                            push.type = 6;
                                            bb++;
                                        }
                                        else {
                                            push.type = 5;
                                            aa++;
                                        }*/

                                        beads.push_back(push);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

};

class OblateSpheroid: public Particle {
public:
    /**
     * @brief generate (x^2+y^2)/(a^2) + (z^2/c^2) = 1, c > 1 prolate, a>1 oblate
     * @param size
     */
    void generate(const int size) {
        beads.clear();
        Vector sample(0.0, 0.0, 0.0, typeNano);
        double c = 0.2;
        double tolerance = 0.00001;
        double step = 1.0/(double(size));
        /*double spacing = 0.1;
        double angle = 0.0;
        double var = 0.0; // spacing*0.5/tan_angle
        double curr_radius = 3.0;
        double temp_x = 0.0;*/

        //
        // Starting point x = a, y=0.0, z=0.0
        //
        /*sample.x = -1;
        sample.y = -1;
        sample.z = -c;
        for(int i=0; i <= size*2; ++i) { // -1 to 1 X
            sample.y = -1;
            for(int i=0; i <= size*2; ++i) { // -1 to 1 Y

                if( ( ((sample.x*sample.x) + (sample.y*sample.y)) ) < (1.0+tolerance) ) {
                    // calc z coordinates

                    sample.z = sqrt(c*c*(1-(sample.x*sample.x + sample.y*sample.y)));
                    if(!isnan(sample.z)) {
                        beads.push_back(sample);
                        sample.z *= -1.0;
                        beads.push_back(sample);
                    }
                }
                sample.y += step;
            }
            sample.x += step;
        }*/
        fibonacci_spheroid(size);
        cerr << beads.size() << endl;
        /*for(int i=0; i< 10; ++i) {
            //if(beads[i].type == typeLig)
                cout << beads[i].x << " " << beads[i].y << " " << beads[i].z << endl;
        }

        exit(0);*/
    }

    void fibonacci_spheroid(int samples) {
        //
        // r(z) = a * sqrt(1-(z*z)/(c*c))
        //
        double c = 0.8;
        const double PI = 3.141592653589793;
        double tolerance = 0.00001;
        //beads.clear();

        double offset = c*2.0/samples;
        double sum = 0.0;
        double increment = PI * (3.0 - sqrt(5.0));

        double x,y,z,r,phi;

        // NOTE sinh^-1 x = ln(x + sqrt( x^2 + 1 ) )
        double transform = ( - 1*c) + (offset / 2);
        z = transform;


        transform = transform / ( 1.0/( c*c * 1.0/log(( sqrt(1-c*c) * z/(c*c) ) + sqrt( ( sqrt(1-c*c) * z/(c*c) )*( sqrt(1-c*c) * z/(c*c) ) + 1 ) ) / sqrt(1-c*c) + z*sqrt( z*z/(c*c*c*c) - z*z/(c*c) + 1 ) ) );

        for(int i=0; i<samples; ++i) {

            z = ((i * offset) - 1*c) + (offset / 2);
            cout << z;
            z = ( c*c * 1.0/log(( sqrt(1-c*c) * z/(c*c) ) + sqrt( ( sqrt(1-c*c) * z/(c*c) )*( sqrt(1-c*c) * z/(c*c) ) + 1 ) ) / sqrt(1-c*c) + z*sqrt( z*z/(c*c*c*c) - z*z/(c*c) + 1 ) );
            cout << " " << z << endl;

            r = sqrt( 1 - ((z*z)/(c*c))  ); // r = sqrt(1 - z*z); -> c=1.0 -> sqrt(1-z*z) - z*z

            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;


            if( ( ( (x*x)+(y*y) ) + ((z*z) / (c*c)) ) > (1 - tolerance) && ( ( (x*x)+(y*y) ) + ((z*z) / (c*c)) ) < (1 + tolerance))
                beads.push_back(Vector(x,y,z,typeNano));
        }
        exit(1);
    }
};

class Sphere: public Particle {
public:
    void generate(const int size) {
        fibonacci_sphere(size);
    }

    void fibonacci_sphere(int samples) {
        const double PI = 3.141592653589793;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        //std::normal_distribution<double> distribution(0.0,0.5);

        double serge=0.1; // separate or merge the distribution at zero

        double rand;
        beads.clear();

        double offset = 2.0/samples;

        double increment = PI * (3.0 - sqrt(5.0));

        double x,y,z,r,phi;
        for(int i=0; i<samples; ++i) {
            y = ((i * offset) - 1) + (offset / 2);
            r = sqrt(1 - pow(y,2));

            phi = i * increment;

            x = cos(phi) * r;
            z = sin(phi) * r;

            // min = -1, max = 1
            rand = distribution(generator);
            if((fabs(z) - serge) > fabs(rand))
                beads.push_back(Vector(x,y,z,typeLig));
            else
                beads.push_back(Vector(x,y,z,typeNano));
        }

        int count = 0;
        for(int i=0; i<beads.size(); ++i) {
            if(beads[i].type == typeLig)
                count++;
        }
        cerr << "Number of Ligands " << count << endl;


        //testDistribution();
    }

    void testDistribution() {

        for(int i=0; i<beads.size(); ++i) {
            if(beads[i].type == typeLig)
                cout << beads[i].z << endl;
        }

        exit(1);
    }
};

//
// Icosahedron
//
int main(int argc, char* argv[]) // // $num of beads per edge, box dimensions X(same as Y) $beg $end, $position in Z, $offset
{
    //
    // Process input parameters
    //
    if(argc == 1) {
        cout << "Parameters:\n" << "num_beads_edge xy_pos_beg xy_pos_end z_pos particle_offset_lammps" << endl;
    }
    const int size = atoi(argv[1]); // num of beads per edge
    string str1(argv[2]);
    string str2(argv[3]);
    string str3(argv[4]);
    string str4(argv[5]);

    double beg = std::stod(str1);  // Box size -> in Lammps defined as begining pos and ending pos
    double end = std::stod(str2);  // Box size -> in Lammps defined as begining pos and ending pos
    double pos = std::stod(str3);  // Center of mass particle position in Z axis
    int offset = std::stoi(str4);  // Lammps format particle numbering offset - adding at the end of file - GIVE PARTICLE COUNT
    const double scale = /*0.5*size*/ 10.0; // Generated beads scale relative to unit length

    Vector com_pos( (-beg + end)*0.5 + beg, (-beg + end)*0.5 + beg, pos); // Center of mass of particle Vector
    offset++;

    //
    // Generate Particle
    //
    OblateSpheroid part;
    part.generate(size);

    //exit(1);

    /*cout << beads.size() << endl;
    cout << "Icosahedron " << size << endl;
    cout << std::setprecision(10);
    for(unsigned int i=0; i<beads.size(); i++) {
        if(rand()%100 > 15)
            cout << beads[i].type << " " << beads[i]*scale << endl;
        else
            cout << beads[i].type << " " << beads[i]*scale << endl;
    }*/

    //cout << " %:" << 100.0*bb/(aa+bb) << endl;

    // Distance of neighbors
    /*double min = 50.0;
    for(int i=0; i < beads.size(); ++i) {
        for(int j=0; j < beads.size(); ++j) {
            if(i != j && min > beads[i].dist(beads[j])) {
                min = beads[i].dist(beads[j]);
            }
        }
    }
    cout << min*scale << endl;
    exit(1);*/


    //part.printLammps( scale, com_pos, offset);
    part.printXYZ( scale, com_pos);

    //cout << "Total: " << beads.size() << ", aa: " << aa << ", bb: " << bb << endl;
    //cout << "Diameter of circumscribed sphere: " << (10/( sqrt(50-10*sqrt(5)) )) * scale  << endl;


    return 0;
}

void Icosahedron::printXYZ(double scale, Vector& com_pos) {
    // XYZ FORMAT
    //shift = beads[0]*scale;
    //shift = shift * -1.0;
    int other=-1;
    std::vector< Vector> neighbors;
    Vector normal;
    Vector dist;
    cout << beads.size() << "\nicosahedron\n";
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 5) // other
            cout << "C2 " << (beads[i]*scale) + com_pos << endl;
        if(beads[i].type == 6) // ligands
            cout << "C1 " << (beads[i]*scale) + com_pos << endl;
        if(beads[i].type == 7) { // edge ligands
            cout << "Type 7" << endl;
            neighbors.clear();
            for(unsigned int j=0; j < beads.size(); ++j) { // search for 4 neighbors
                if(beads[j].type == 5 && i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    neighbors.push_back( beads[j] );
                    cout << "neighbors push " << beads[j] << endl;
                }
            }
            normal = (neighbors[1]-neighbors[0]).cross(neighbors[2]-neighbors[0]); // calculate normal
            dist = beads[i] - neighbors[0];
            normal = normal * (-scale*(dist.dot(normal) / (normal.size() * normal.size() ) ) );
            cout << "C1 " << (beads[i]*scale) +com_pos +normal << endl;
        }
        if(beads[i].type == 8) { // Tops
            for(unsigned int j=0; j < beads.size(); ++j) {
                if(i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    other = j;
                }
            }
            cout << "C1 " << (beads[i]*(scale/( beads[i].size() * beads[i].size() ) * beads[i].dot(beads[other]) ) ) + com_pos << endl;
        }
    }
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 8) { // Tops
            for(unsigned int j=0; j < beads.size(); ++j) {
                if(i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    beads[j].type=10;
                }
            }
        }
    }
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 9) {// edge others
            neighbors.clear();
            for(unsigned int j=0; j < beads.size(); ++j) { // search for 4 neighbors
                if(beads[j].type == 5 && i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    neighbors.push_back( beads[j] );
                }
            }
            //normal = Vector(0.0, 0.0, 0.0);
            //normal = (neighbors[1]-neighbors[0]).cross(neighbors[2]-neighbors[0]); // calculate normal
            //dist = beads[i] - neighbors[0];
            //normal = normal * (-scale*(dist.dot(normal) / (normal.size() * normal.size() ) ) );
            cout << "C2 " << (beads[i]*scale) +com_pos /*+normal*/ << endl;
        }
    }
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 10) {// edge others
            cout << "C2 " << (beads[i]*scale) + com_pos << endl;
        }

    }
}

void Icosahedron::printLammps(double scale, Vector& com_pos, int offset) {
    Vector a,b;

    int other=-1;
    std::vector< Vector> neighbors;
    Vector normal;
    Vector dist;
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 8) { // Tops
            for(unsigned int j=0; j < beads.size(); ++j) {
                if(i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    beads[j].type=10;
                }
            }
        }
    }
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 9) {// edge others
            neighbors.clear();
            for(unsigned int j=0; j < beads.size(); ++j) { // search for 4 neighbors
                if(beads[j].type == 5 && i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    neighbors.push_back( beads[j] );
                }
            }
            normal = Vector(0.0, 0.0, 0.0);
            normal = (neighbors[1]-neighbors[0]).cross(neighbors[2]-neighbors[0]); // calculate normal
            dist = beads[i] - neighbors[0];
            normal = normal * (-scale*(dist.dot(normal) / (normal.size() * normal.size() ) ) );

            beads[i] = (beads[i]*scale) +com_pos +normal;
            beads[i].type = 9;
            //cout << "C2 " << (beads[i]*scale) +shift +normal << endl;
        }
    }
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 5) { // other
            beads[i] = (beads[i]*scale) +com_pos;
            beads[i].type = 5;
            //cout << "C2 " << (beads[i]*scale) + shift << endl;
        }
        if(beads[i].type == 6) { // ligands
            beads[i] = (beads[i]*scale) +com_pos;
            beads[i].type = 6;
            //cout << "C1 " << (beads[i]*scale) + shift << endl;
        }
        if(beads[i].type == 7) { // edge ligands
            neighbors.clear();
            for(unsigned int j=0; j < beads.size(); ++j) { // search for 4 neighbors
                if(beads[j].type == 5 && i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    neighbors.push_back( beads[j] );
                }
            }
            normal = (neighbors[1]-neighbors[0]).cross(neighbors[2]-neighbors[0]); // calculate normal
            dist = beads[i] - neighbors[0];
            normal = normal * (-scale*(dist.dot(normal) / (normal.size() * normal.size() ) ) );

            beads[i] = (beads[i]*scale) +com_pos +normal;
            beads[i].type = 7;
            //cout << "C1 " << (beads[i]*scale) +shift +normal << endl;
        }
        if(beads[i].type == 8) { // Tops
            for(unsigned int j=0; j < beads.size(); ++j) {
                if(i != j && ( 0.53/scale > beads[i].dist(beads[j]) ) ) {
                    other = j;
                }
            }
            beads[i] = (beads[i]*(scale/( beads[i].size() * beads[i].size() ) * beads[i].dot(beads[other]) ) ) + com_pos;
            //cout << "C1 " << (beads[i]*(scale/( beads[i].size() * beads[i].size() ) * beads[i].dot(beads[other]) ) ) + shift << endl;
            beads[i].type = 8;
        }
    }
    for(unsigned int i=0; i<beads.size(); i++) {
        if(beads[i].type == 10) {// edge others
            beads[i] = (beads[i]*scale) + com_pos;
            //cout << "C2 " << (beads[i]*scale) + shift << endl;
            beads[i].type = 10;
        }

    }
    for(unsigned int i=0; i < beads.size(); i++) {
        if(beads[i].type == 10) beads[i].type = 5;
        if(beads[i].type == 9) beads[i].type = 5;
        if(beads[i].type == 8) beads[i].type = 6;
        if(beads[i].type == 7) beads[i].type = 6;
    }


    for(unsigned int i=0; i < beads.size(); i++) {
        cout << i+offset << " " << ((i+offset)%1000) +1 << " " << beads[i].type << " " << 0 << " " << beads[i] << " 0 0 0 #ico" << endl;
    }
}

