#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <complex>

#include "sphere.h"

using namespace std;

class Ellipsoid_Surf
{
public:
    Ellipsoid_Surf(double a, double b, double c, int samples)
    {
        this->a=a;
        this->b=b;
        this->c=c;
        this->samples = samples;
        surf_distrib_Z.reserve(samples);

        //
        // construct surface distribution in Z direstion
        //
        elipsoid_surface_distrib();

        //
        // Construct points in Z axis, with constant surface area increment between points
        //
        generate_const_surf_distrib();

        //
        // Construct azimuth values, with constant perimeter distance increment between points
        //
        generate_const_perimeter_distrib();
    }

    const double PI = 3.141592653589793;
    const double deg_to_rad = 0.0174532925199;
    const double rad_to_deg = 57.2957795131;
    const int count = 36000;
    double a=1.0;
    double b=1.0;
    double c=1.0;
    int samples=0;

    vector<double> surf_distrib_Z;
    vector<double> const_surf_Z_distrib; // z coordinate in <-c; c>, samples, each representing constant increment in surf area
    vector<double> const_perimeter_Phi_distrib;




    void generate_const_perimeter_distrib()
    {
        Atom a1,a2;
        double deg_increment = 0.00001;
        double sum=0.0;
        double p_ellipse = perimeter_ellipse(a, b, deg_increment);
        double perimeter_increment = p_ellipse / count;

        const_perimeter_Phi_distrib.push_back(0.0);
        for ( double i = 0.0; i <= 360.0-deg_increment; i += deg_increment )
        {
            a1 = get_atom(i, 0.0);
            a2 = get_atom(i+deg_increment, 0.0);

            sum += (a1-a2).size();

            if(sum > perimeter_increment)
            {
                const_perimeter_Phi_distrib.push_back( i );
                sum = 0.0;
            }
        }
    }




    /**
     * @brief azimuth_conversion
     * @param phi - distribution from sphere
     * @return - distribution from ellipse
     */
    double azimuth_conversion(double phi)
    {
        phi = fmod( phi * rad_to_deg, 360.0);
        int index = phi / (360.0/count);
        return deg_to_rad * const_perimeter_Phi_distrib[index];
    }




    double surface_Area(Atom from = Atom(0,0,-999), Atom to = Atom(0,0,999))
    {
        double sum=0.0;
        double z_increment = c/samples;
        //double tolerance = 0.1*z_increment;
        double z;

        for ( int i = -samples; i < samples; i+=2 )
        {
            z = i*z_increment;
            if( (z >= from.z && z < to.z) )
                sum += surf_distrib_Z[(i+samples)/2];
        }
        return sum;
    }




    double surface_Area_approx()
    {
        if(a==1.0 && b==1.0 && c==1.0)
            return 4*PI;

        double value;
        double aa = pow(a, 1.6075);
        double bb = pow(b, 1.6075);
        double cc = pow(c, 1.6075);
        value = aa*bb + aa*cc + bb*cc;
        value /= 3.0;
        return 4*PI*pow(value, 1.0/1.6075); // aproximate, error 1%
    }




private:

    void generate_const_surf_distrib()
    {
        int real_samples = samples*10000;

        long double cumul_surf = 0.0;
        long double cumul_surf_prev = 0.0;
        long double total_surf = 0.0;

        long double z_increment = c/real_samples;
        long double j_z;

        for ( int j = -real_samples+2; j <= real_samples; j+=2 )
        {
            j_z = j*z_increment;
            total_surf += lateral_area_elliptic_frustrum_approx_2(j_z, (j-2)*z_increment, 2.0*z_increment);
        }

        long double surf_increment = total_surf / (samples-1);

        const_surf_Z_distrib.push_back(-c);
        for ( int j = -real_samples+2; j <= real_samples; j+=2 )
        {
            j_z = j*z_increment;
            cumul_surf += lateral_area_elliptic_frustrum_approx_2(j_z, (j-2)*z_increment, 2.0*z_increment);
            if(cumul_surf > cumul_surf_prev + surf_increment && const_surf_Z_distrib.size() < samples)
            {
                cumul_surf_prev = cumul_surf;
                const_surf_Z_distrib.push_back(j_z);
            }
        }
        const_surf_Z_distrib.push_back(c);

        if(const_surf_Z_distrib.size() != samples)
        {
            cerr << "Generated " << const_surf_Z_distrib.size() << "instead of " << samples << endl;
        }

        /*for(auto item : const_surf_Z_distrib)
        {
            cout << item << endl;
        }*/
        //exit(1);
    }




    /**
     * @brief construct a general elipsoid in a naive manner, save surface distribution
     * @param container
     */
    void elipsoid_surface_distrib()
    {
        double z_increment = c/samples;

        double j_z;
        for ( int j = -samples+2; j <= samples; j+=2 )
        {
            j_z = j*z_increment;
            surf_distrib_Z.push_back( lateral_area_elliptic_frustrum_approx_2(j_z, (j-2)*z_increment, 2.0*z_increment) );
        }
    }




    /**
     * @brief lateral_area_elliptic_frustrum_approx
     * Deprecated - use lateral_area_elliptic_frustrum_approx_2
     *
     * Calculated the lateral area of an eliptical frustrum
     *
     * @return
     */
    double lateral_area_elliptic_frustrum_approx(double z, double z2, double h)
    {
        Atom a1,a2,a3,a4;
        double val=0.0;
        double deg_increment = 0.01;
        double sum=0.0;

        for ( double i = 0.0; i <= 360.0-deg_increment; i += deg_increment )
        {
            a1 = get_atom(i, z2);
            a2 = get_atom(i, z);

            a3 = get_atom(i+deg_increment, z);
            a4 = get_atom(i+deg_increment, z2);

            val=0.0;
            val += 0.5*(a1-a2).cross(a1-a3).size();
            val += 0.5*(a1-a4).cross(a1-a3).size();
            sum += val;
        }

        return sum;
    }




    double lateral_area_elliptic_frustrum_approx_2(double z, double z2, double h)
    {
        if(z > c)
            z = c;
        if(z < -c)
            z = -c;
        if(z2 > c)
            z2 = c;
        if(z2 < -c)
            z2 = -c;

        double x1 = a * sqrt( 1 - ((z*z) / (c*c)) );
        double y1 = b * sqrt( 1 - ((z*z) / (c*c)) );

        double x2 = a * sqrt( 1 - ((z2*z2) / (c*c)) );
        double y2 = b * sqrt( 1 - ((z2*z2) / (c*c)) );

        if(x1 != x1 && y1 != y1 && x2 != x2 && y2 != y2)
        {
            cerr << "Ellipsoid::lateral_area_elliptic_frustrum_approx_2 :: x1 is nan" << endl;
            exit(1);
        }

        double h_total = 0.0;

        //
        // double has 15 decimal points
        // a very small |x2-x1| produces bad results. We approximate those cases by eliptical cylinder
        //
        if( abs(x2-x1) < 1e-10 )
        {
            return perimeter_ellipse(x1,y1)*h;
        }

        if( x2 > x1)
        {
            h_total = h * x2 / (x2-x1);
            return lateral_area_elliptic_cone_approx(x2, y2, h_total) - lateral_area_elliptic_cone_approx(x1, y1, h_total-h);
        }
        else
        {
            h_total = h * x1 / (x1-x2);

            return lateral_area_elliptic_cone_approx(x1, y1, h_total) - lateral_area_elliptic_cone_approx(x2, y2, h_total-h);
        }

        return 0.0;
    }




    double lateral_area_elliptic_cone_approx(double x, double y, double h)
    {
        return PI/2.0 * (x*sqrt(y*y + h*h) + y*sqrt(x*x + h*h) );
    }




    double perimeter_ellipse(double x, double y, double deg_increment = 0.1)
    {
        Atom a1,a2;
        double sum=0.0;

        for ( double i = 0.0; i <= 360.0-deg_increment; i += deg_increment )
        {
            a1 = get_atom(i, 0.0);
            a2 = get_atom(i+deg_increment, 0.0);

            sum += (a1-a2).size();
        }

        return sum;
    }




    /**
     * @brief get_atom
     *
     * lambda - azimuth
     *
     * x = a cos(phi) cos(lambda)
     * y = b cos(phi) sin(lambda)
     * z = c sin(phi)
     *
     * cos(phi) = sqrt( 1 - ((z*z) / (c*c)) )
     * phi - reduced latitude
     *
     * @return
     */
    Atom get_atom(double deg, double z)
    {
        Atom add;

        add.x = a * sqrt( 1 - ((z*z) / (c*c)) ) * cos (deg*deg_to_rad);
        add.y = b * sqrt( 1 - ((z*z) / (c*c)) ) * sin (deg*deg_to_rad);
        add.z = z;

        return add;
    }
};




/**
 * @brief The OblateSpheroid class - generate oblate spheroid with surface area of 4*PI (equivalent to sphere of unit radius)
 */
class Ellipsoid: public Sphere {
public:   
    Ellipsoid() : Sphere("Ellipsoid") {}

    /**
     * @brief generate (x^2+y^2)/(a^2) + (z^2/c^2) = 1, c > 1 prolate, a>1 oblate
     * @param size
     */
    void generate( Data& data )
    {
        vector<Atom> ligand;

        typeNano = data.in.atom_type;
        typeLig = data.in.atom_type + 1;

        int orientations = orientZ;

        //
        // construct base shape
        //
        fibonacci_tri_axial_spheroid(beads, data.in.num_of_beads, data.in.a, data.in.b, data.in.c, typeNano, orientations);

        //
        // construct ligand shape for selection purposes
        //
        fibonacci_tri_axial_spheroid(ligand, data.in.num_lig, data.in.a, data.in.b, data.in.c, typeTemp, orientations);

        //
        // Validation
        //
        /*if( validate(beads, data.in.a, data.in.b, data.in.c) )
            cout << "Shape validation OK" << endl;
        else
            cout << "Shape validation NOK" << endl;*/

        //
        // Test distribution vs surface area
        //
        /*distribution(beads, data.in.a, data.in.b, data.in.c);
        exit(1);*/

        //
        // Generate ligands
        //
        gen_ligands( data, ligand, data.in.patch_1, typeNano, typeTemp);
        gen_ligands( data, ligand, data.in.patch_2, typeNano, typeTemp);

        //
        // Set moltag, and offset N for lammps output style
        //
        int i=0;
        for(auto& item : beads)
        {
            item.mol_tag = data.in.mol_tag;
            item.N = i+1+data.in.offset+data.all_beads.size();
            ++i;
        }
    }

protected:

    inline double linear_Z(int i, int samples, double c)
    {
        double z_increment = c*2.0/samples;
        return - c + (i * z_increment) + (z_increment * 0.5); // z linearly distributed from -c to c
    }

    void fibonacci_tri_axial_spheroid(vector<Atom>& container, int samples, double a, double b, double c, const int type, int orientation)
    {
    	const double PI = 3.141592653589793;
        Ellipsoid_Surf ellipsoid(a,b,c,samples);
        double angle_increment = PI * ( 3.0 - sqrt(5.0));

        //vector<double> cumul_fce = elipsoid_cumulative_function(1.0, b, c);

        double x,y,z,phi;

        for(int i=0; i<samples; ++i)
        {
            //z = linear_Z(i, samples, c);
            z = ellipsoid.const_surf_Z_distrib[i];
            phi = i * angle_increment;

            x = a * sqrt( 1 - ((z*z) / (c*c)) ) * cos ( ellipsoid.azimuth_conversion(phi) );
            y = b * sqrt( 1 - ((z*z) / (c*c)) ) * sin ( ellipsoid.azimuth_conversion(phi) );

            //if( ( ((x*x)+(y*y)) + ((z*z)/(c*c)) ) > (1 - tolerance) && ( ((x*x)+(y*y)) + ((z*z)/(c*c)) ) < (1 + tolerance)) {
                if(orientX == orientation)
                    container.push_back(Atom(z,x,y,type));
                if(orientY == orientation)
                    container.push_back(Atom(x,z,y,type));
                if(orientZ == orientation)
                    container.push_back(Atom(x,y,z,type));
            //}
        }
    }

    /**
     * @brief validate general elipsoid
     * x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
     * @return
     */
    bool validate(vector<Atom>& container, double a, double b, double c)
    {
        double limit=0.0001;
        double value;
        for (auto& item : container)
        {
            value = (item.x*item.x) / (a*a);
            value += (item.y*item.y) / (b*b);
            value += (item.z*item.z) / (c*c);

            //cout << item.x << " " << item.y << " " << item.z << " " << value << endl;

            if( value > 1+limit || value < 1-limit)
                return false;
        }
        return true;
    }


    void distribution(vector<Atom>& container, double a, double b, double c)
    {
        cout << "Testing distribution, unit size" << endl;
        const int size = 10;
        int x[size] = {0};
        int y[size] = {0};
        int z[size] = {0};
        double x_r, y_r, z_r;
        double x_i, y_i, z_i;
        double total_surf;
        Ellipsoid_Surf ellipsoid_1(c,b,a, 1000);
        Ellipsoid_Surf ellipsoid_2(a,c,b, 1000);
        Ellipsoid_Surf ellipsoid_3(a,b,c, 1000);

        for(auto& item : container)
        {
            ++x[ (int) (  (item.x+a)*(size/(2*a))  ) ];
            ++y[ (int) (  (item.y+b)*(size/(2*b))  ) ];
            ++z[ (int) (  (item.z+c)*(size/(2*c))  ) ];

            //cout << (int) ((beads[i].z+1)*(size/2)) << endl;
        }

        /*cout << "X " << get_min_x() << " to " << get_max_x() << endl;
        cout << "Y " << get_min_y() << " to " << get_max_y() << endl;
        cout << "Z " << get_mim_z() << " to " << get_max_z() << endl;*/

        total_surf = ellipsoid_1.surface_Area(Atom(0,0, -c), Atom(0,0, c) );
        for(int i=0; i<size; ++i)
        {
            x_r = 100.0 * x[i]/beads.size();
            y_r = 100.0 * y[i]/beads.size();
            z_r = 100.0 * z[i]/beads.size();

            x_i = ellipsoid_1.surface_Area(Atom(0,0, -a + (i)*(2.0*a/size) ), Atom(0,0, -a + (i+1)*(2.0*a/size) ) );
            x_i = 100.0 * x_i / total_surf;

            y_i = ellipsoid_2.surface_Area(Atom(0,0, -b + (i)*(2.0*b/size) ), Atom(0,0, -b + (i+1)*(2.0*b/size) ) );
            y_i = 100.0 * y_i / total_surf;

            z_i = ellipsoid_3.surface_Area(Atom(0,0, -c + (i)*(2.0*c/size) ), Atom(0,0, -c + (i+1)*(2.0*c/size) ) );
            z_i = 100.0 * z_i / total_surf;
            std::cout << std::fixed << std::setprecision(2);
            cout << "real: "   << std::setw(5) << x_r << " " << std::setw(5) << y_r << " " << std::setw(5) << z_r;
            cout << " ideal: " << std::setw(5) << x_i << " " << std::setw(5) << y_i << " " << std::setw(5) << z_i << endl;
        }
    }

private:

    double get_max_x()
    {
        double max = beads[0].x;
        for(int i=0; i<beads.size(); ++i)
        {
            if(beads[i].x > max)
                max = beads[i].x;
        }
        return max;
    }

    double get_max_y()
    {
        double max = beads[0].y;
        for(int i=0; i<beads.size(); ++i)
        {
            if(beads[i].y > max)
                max = beads[i].y;
        }
        return max;
    }

    double get_max_z()
    {
        double max = beads[0].z;
        for(int i=0; i<beads.size(); ++i)
        {
            if(beads[i].z > max)
                max = beads[i].z;
        }
        return max;
    }

    double get_min_x()
    {
        double min = beads[0].x;
        for(int i=0; i<beads.size(); ++i)
        {
            if(beads[i].x < min)
                min = beads[i].x;
        }
        return min;
    }

    double get_min_y()
    {
        double min = beads[0].y;
        for(int i=0; i<beads.size(); ++i)
        {
            if(beads[i].y < min)
                min = beads[i].y;
        }
        return min;
    }

    double get_mim_z()
    {
        double min = beads[0].z;
        for(int i=0; i<beads.size(); ++i)
        {
            if(beads[i].z < min)
                min = beads[i].z;
        }
        return min;
    }
};

#endif // ELLIPSOID_H
