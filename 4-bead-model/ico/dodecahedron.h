#ifndef DODECAHEDRON_H
#define DODECAHEDRON_H

#include "pentamer.h"
#include "icosahedron.h"



class Dodecahedron : public Particle
{
public:
    IcoVertice edge;

    myFloat shiftTop;
    myFloat shiftMid;
    myFloat shiftBot;
    myFloat fscale;
    myFloat scale = 1.0;
    Atom com_pos = Atom(0.0, 0.0, 0.0);
    myFloat c;
    int size;
    int offset;

    enum interface: unsigned short { INNER=0, OUTER=1, UNIFORM=2 };
    enum type: unsigned short { BASIC, INTERFACE_SHIFT, INTERFACE_SHIFT_ONLY };

    interface interface_type = UNIFORM;

    Dodecahedron(string name) : Particle(name) {}

    virtual void generate( Data& data ) override
    {
        // LOAD DATA
        this->size = data.in.num_of_beads;
        this->c = data.in.c;
        this->pbc.box.x = data.in.boxp.x - data.in.boxm.x;
        this->pbc.box.y = data.in.boxp.y - data.in.boxm.y;
        this->pbc.box.z = data.in.boxp.z - data.in.boxm.z;
        this->offset = data.all_beads.size();

        init(size);// set shift, fscale
        info();

        this->scale = data.in.scale;
        this->com_pos = data.in.com_pos;

        //
        // Define sigma values for forcefield
        //
        set_sigma();
        cerr << "c: " << c << ", scale: " << scale << ", fscale: " << fscale << endl;
        printSigma();

        // 12 pentamers in dodecahedron
        //vector<int> penta = {0,1,2,3,4,5,6,7,8,9,10,11}; // Entire capside
        //vector<int> penta = {0,1,2,3,4, 5,6,7,8,9, 10,11,0,1,2, 3,4,5,6,7, 8,9,10,11,0, 1,2,3,4,5 }; // Entire capside
        //vector<int> penta = {0,1}; // patological case
        //vector<int> penta = {0,1,2}; // patological case 2
        vector<int> penta = {0}; // monomer case
        //vector<int> penta = {0,2}; // dimer
        //vector<int> penta = {0,4,5}; // trimer

        vector<Atom> penta_particles;
        //
        // Generate pentamers + bonds
        //
        for(int i=0; i<penta.size(); ++i)
        {
            this->offset = data.all_beads.size() + beads.size();
            generatePenta(penta_particles, penta[i]);
            //random_Translate_Rotate(penta_particles);
            /*if(i==1) { // patological case 1 rotate
                for(int j=0; j<penta_particles.size(); ++j)
                    penta_particles[j].z -= 1.4*i;

                clusterRotate(penta_particles, 10*degToRad, Atom(0,0,1) );
            }*/

           beads.insert(beads.end(), penta_particles.begin(), penta_particles.end());
           penta_particles.clear();
        }

        // set atoms N
        for(unsigned int i=0; i< beads.size(); ++i)
        {
            beads[i].N = i+1;
        }
        cerr << "Bonds: " << bonds.size() << endl;
    }

private:
    myFloat separation_same, separation_other, separation_other_shifted, separation_same_shifted;

    void mixing_rules() {
        for(int i = 0; i< sigma_size; ++i) {
            for(int j = 0; j< sigma_size; ++j) {
                if(i != j) {
                    sigma[i][j] = 0.5*(sigma[i][i] + sigma[j][j]);
                    sigma[j][i] = 0.5*(sigma[i][i] + sigma[j][j]);
                }
            }
        }
    }

    /**
     * @brief set_sigma - Exclude Intra-pentamer interactions - only bonds make pentamer shape
     */
    void set_sigma() {
        sigma_size=11;

        sigma[0][0] = separation_other* c*scale/fscale; // outer body
        sigma[1][1] = separation_other* c*scale/fscale; // interface outer
        sigma[2][2] = separation_other* c*scale/fscale; // interface outer
        sigma[3][3] = separation_other* c*scale/fscale; // interface outer

        sigma[4][4] = separation_other*scale/fscale; // inner body
        sigma[5][5] = separation_other*scale/fscale; // interface inner
        sigma[6][6] = separation_other*scale/fscale; // interface inner
        sigma[7][7] = separation_other*scale/fscale; // interface inner

        sigma[8][8] = separation_other*(1.0 + 0.5*(c-1.0))*scale/fscale; // interface middle
        sigma[9][9] = separation_other*(1.0 + 0.5*(c-1.0))*scale/fscale; // interface middle
        sigma[10][10] = separation_other*(1.0 + 0.5*(c-1.0))*scale/fscale; // interface middle

        mixing_rules();

        sigma_cosatt[1][2] = sigma_cosatt[2][1] = true;
        sigma_cosatt[5][6] = sigma_cosatt[6][5] = true;
        sigma_cosatt[8][9] = sigma_cosatt[9][8] = true;
    }




    /**
     * @brief gen_bonds
     * @param dodeca
     * @return particles in pentamer
     */
    void generatePenta(vector<Atom>& particles, int a)
    {
        if(!particles.empty()) {
            cerr << "particles not empty" << endl;
            exit(1);
        }

        //
        // Generate dodecahedron just for bonds
        //
        int inner_start, inner_end, outer_start, outer_end, middle_start, middle_end;

        inner_start = particles.size();
        pentaFrame(particles, a, 1.0/fscale, type::INTERFACE_SHIFT, 5);
        inner_end = particles.size();

        outer_start = particles.size();
        pentaFrame(particles, a, c/fscale, type::INTERFACE_SHIFT, 1);
        outer_end = particles.size();

        middle_start = particles.size();
        pentaFrame(particles, a, (1.0 + 0.5*(c-1.0))/fscale, type::INTERFACE_SHIFT_ONLY, 9);
        middle_end = particles.size();

        //
        // Generate bonds
        //
        /*generateSame(particles, 1.0, inner_start, inner_end);
        generateSame(particles, c, outer_start, outer_end);
        generateSame(particles, (1.0 + 0.5*(c-1.0)), middle_start, middle_end);*/

        generateCross2(particles, (1.0 + 0.53*(c-1.0)), middle_start, middle_end, outer_start, outer_end); // adds 1 bond to middle and outer
        generateCross2(particles, (1.0 + 0.53*(c-1.0)), inner_start, inner_end, middle_start, middle_end); // adds 1 bond to middle and inner
        generateCross3(particles, inner_start, inner_end, outer_start, outer_end); // adds 1 bond to inner and outer
        generateCross4(particles, inner_start, inner_end, outer_start, outer_end); // adds bonds to inner and outer interface

        /*generateCrossSpecial(particles);
        generateSameSpecial(particles);*/

        //
        // update 2019 18.6. -> adding more bonds to interface to prevent excessive interface beads movement
        //

        for(int i=0; i<particles.size(); ++i) {
            particles[i].mol_tag=a;
            if(particles[i].to_type != -1)
                particles[i].type = particles[i].to_type;
        }
    }

    /**
     * @brief pentaFrame - generate pentamer
     * @param beads - stl container for generated beads
     * @param a - identifier of pentamer
     * @param scale
     * @param shift
     * @param tri_type
     * @return
     */
    myFloat pentaFrame(vector<Atom>& particles, int a, myFloat factor, unsigned short tri_type, int types)
    {
        myFloat ss=0.0;

        for(int b=0; b<12; b++) {
            if(edge[a].isNeighbor(edge[b]) ) {
                for(int c=0; c<12; c++) {
                    if(edge[b].isNeighbor(edge[c]) && edge[a].isNeighbor(edge[c]) ) { // for face => triangular area:
                        for(int d=0; d<12; d++) {
                            if(edge[d].isNeighbor(edge[a]) && edge[d].isNeighbor(edge[b]) && d != c ) {
                                //
                                // a,b,c,d vertices of triangle
                                //
                                Atom center = (edge[a] + edge[b] + edge[c]) * (1.0/3.0);
                                Atom center2 = (edge[a] + edge[b] + edge[d]) * (1.0/3.0);
                                Atom cc = (center + center2) * 0.5;

                                myFloat dist = cc.dot(edge[a]) / edge[a].size();
                                Atom edge_A_scaled = edge[a] * ( dist / cc.size() );

                                ss = edge_A_scaled.size();

                                trianle(particles, edge_A_scaled, center, center2, a,b,c, tri_type, types, factor);
                            }
                        }
                    }
                }
            }
        }

        return ss;
    }

    void trianle(vector<Atom>& particles, Atom va, Atom vb, Atom vc, int a, int b, int c, int tri_type, int types, myFloat factor) { // edge_A_scaled, center, center2, size, a,b,c

        Atom cc = (vb + vc) * 0.5;
        Atom vec = cc - va;
        vec.normalise();

        if(( (  (edge[a] - vb).cross( (edge[b] - vb) )  ) + vb ).size() < 1.0) {

            for(int i=0; i<size-1; i++) { // separate into smaller triangle
                for(int j=0; j<=i; j++) {

                    switch( tri_type ) {
                        case type::BASIC : basic(particles, va, vb, vc, i, j, types, factor);
                        break;
                        case type::INTERFACE_SHIFT : interface(particles, va, vb, vc, i, j, types, factor);
                        break;
                        case type::INTERFACE_SHIFT_ONLY : interface_only(particles, va, vb, vc, i, j, types, factor);
                        break;
                    }
                }
            }
        }
    }

    Atom getPosition( Atom va, Atom vb, Atom vc, int size, int i, int j, myFloat factor)
    {
        Atom push;
        push = getSubPoint(va, vb, vc, size, i, j);
        push *= factor;
        push *= scale;
        push += com_pos;

        return push;
    }

    void basic(vector<Atom>& particles, Atom va, Atom vb, Atom vc, int i, int j, int types, myFloat factor)
    {
        Atom push = getPosition(va, vb, vc, size, i, j, factor);

        if( isSame(particles, push) )
        {
            push.type = types; // Upper and bottom pentamer body, type 1 and 6
            particles.push_back(push);
        }
    }

    bool isInterface(int i)
    {
        return (i == size-2);
    }

    bool isFirstPatch(int j)
    {
        return (j >= 1 && j < size/2);
    }

    bool isInterfaceEdge(int j) // Interface, shift type from non-interacting to interacting, j=0 start, j=size-1 pentamer edge
    {
        return (j==0 || j == size-2); // Pentamer general
    }

    void interface(vector<Atom>& particles, Atom va, Atom vb, Atom vc, int i, int j, int types, myFloat factor) // type 1 and  6
    {
        Atom push = getPosition(va, vb, vc, size, i, j, factor);

        if( isInterface(i) )
        {
            if( isInterfaceEdge(j) )
            {
                if( isSame(particles, push) )
                {
                    push.type = types+3;  // bead on pentamer interface edge
                    particles.push_back(push);
                }
                return;
            }

            if( isFirstPatch(j) )
            {
                push = getSubPoint(va, vb, vc, size, i, j) + getSubPoint(va, vb, vc, size, i, j-1);
                push *= 0.5;

                Atom cc = (vb + vc) * 0.5;
                Atom vec = cc - va;
                vec.normalise();

                push = push + vec*shiftTop*factor;

                push *= factor;
                push *= scale;
                push += com_pos;

                if( isSame(particles, push) )
                {
                    push.type = types+2; // interface 2
                    switch( interface_type )
                    {
                        case INNER: if(j == 1 || j == 2) {push.to_type=types+3;} break;
                        case OUTER: if(j == 3) {push.to_type=types+3;} break;
                        case UNIFORM: if(j == 1 ) {push.to_type=types+3;} break;
                        default: push.to_type=types+3; break;
                    }
                    particles.push_back(push);
                }
                return;
            }

            if( !isFirstPatch(j) )
            {
                if( isSame(particles, push) )
                {
                    push.type = types+1; // interface 1
                    switch( interface_type )
                    {
                        case INNER: if(j == 6) {push.to_type=types+3;} break;
                        case OUTER: if(j == 4 || j == 5) {push.to_type=types+3;} break;
                        case UNIFORM: if( false ) {push.to_type=types+3;} break;
                        default: push.to_type=types+3; break;
                    }
                    particles.push_back(push);
                }
                return;
            }

        }

        if( !isInterface(i) )
        {
            if( isSame(particles, push) )
            {
                push.type = types; // general body
                particles.push_back(push);
            }
        }
    }

    void interface_only(vector<Atom>& particles, Atom va, Atom vb, Atom vc, int i, int j, int types, myFloat factor) // types = 11
    {
        Atom push = getPosition(va, vb, vc, size, i, j, factor);

        if(i == size -2)
        {
            if(j==0 || j == size-2)
            {
                if( isSame(particles, push) )
                {
                    push.type = types+2;
                    particles.push_back(push);
                }
                return;
            }

            if(j >= 1 && j < size/2)
            {
                push = getSubPoint(va, vb, vc, size, i, j) + getSubPoint(va, vb, vc, size, i, j-1);
                push *= 0.5;

                Atom cc = (vb + vc) * 0.5;
                Atom vec = cc - va;
                vec.normalise();

                push = push + vec*shiftTop*factor;

                push *= factor;
                push *= scale;
                push += com_pos;

                if( isSame(particles, push) )
                {
                    push.type = types; // type 9 - interacting bead
                    switch( interface_type )
                    {
                        case INNER: if(j == 1 || j == 2) {push.to_type = types+2;} break;
                        case OUTER: if(j == 3) {push.to_type = types+2;} break;
                        case UNIFORM: if(j == 1 ) {push.to_type = types+2;} break;
                        default: push.to_type = types+2; break;
                    }
                    particles.push_back(push);
                }
            }
            else
            {
                if( isSame(particles, push) )
                {
                    push.type = types+1;
                    switch( interface_type )
                    {
                        case INNER: if(j == 6) {push.to_type = types+2;} break;
                        case OUTER: if(j == 4 || j == 5) {push.to_type = types+2;} break;
                        case UNIFORM: if( false ) {push.to_type = types+2;} break;
                        default: push.to_type = types+2; break;
                    }
                    particles.push_back(push);
                }
            }
        }
    }

    //
    // BONDS
    //

    /**
     * @brief generateSame - generate bonds between particles in container particles
     * @param particles
     * @param ccc
     * @param fscale
     * @param start
     * @param end
     */
    void generateSame(vector<Atom> particles, myFloat ccc, int start, int end)
    {
        myFloat limit = (ccc*separation_same/fscale*scale)+0.1; // ~13.33
        cerr << "Limit: " << limit << endl;
        for(int i=start; i<end; ++i)
        {
            for(int j=i+1; j<end; ++j)
            {
                if( particles[i].dist(particles[j]) < limit )
                {
                    //                    N             type          at1       at2       dist
                    bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_same" ) );
                }
            }
        }
    }

    void generateCross2(vector<Atom> particles, myFloat c, int inner_start, int inner_end, int outer_start, int outer_end)
    {
        for(int i=inner_start; i<inner_end; ++i)
        {
            for(int j=outer_start; j<outer_end; ++j)
            {
                if(particles[i].dist(particles[j]) < (c-1.01)*scale)
                {
                    bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_cross2" ) );
                }
            }
        }
    }

    void generateCross3(vector<Atom> particles, int inner_start, int inner_end, int outer_start, int outer_end)
    {
        for(int i=inner_start; i<inner_end; ++i)
        {
            for(int j=outer_start; j<outer_end; ++j)
            {
                if(i!=j && particles[i].type == 5 && particles[j].type == 1 && particles[i].dist(particles[j]) < 20.5  ) // (c-0.995)*scale
                {
                    bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_cross" ) );
                }
            }
        }
    }

    void generateCross4(vector<Atom> particles, int inner_start, int inner_end, int outer_start, int outer_end)
    {
        for(int i=inner_start; i<inner_end; ++i)
        {
            for(int j=outer_start; j<outer_end; ++j)
            {
                if(i!=j && (particles[i].type == 6 || particles[i].type == 7 || particles[i].type == 8) && particles[j].type == 1
                        && particles[i].dist(particles[j]) < 22.0 ) // (c-0.987)*scale
                {
                    bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_cross" ) );
                }
                if(i!=j && particles[j].type == 2 && particles[i].type == 5 && particles[i].dist(particles[j]) < (c-0.9975)*scale )
                {
                    bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_cross" ) );
                }
                if(i!=j && particles[j].type == 3 && particles[i].type == 5 && particles[i].dist(particles[j]) < 21.0 ) //(c-0.995)*scale
                {
                    bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_cross" ) );
                }
                if(i!=j && particles[j].type == 4 && particles[i].type == 5 && particles[i].dist(particles[j]) < (c-0.980)*scale )
                {
                    bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_cross" ) );
                }
            }
        }
    }

    /**
     * @brief generateCross - Deprecated
     * @param particles
     * @param inner_start
     * @param inner_end
     * @param outer_start
     * @param outer_end
     * @param offset
     */
    void generateCross(vector<Atom> particles,  int inner_start, int inner_end, int outer_start, int outer_end)
    {
        int  j = outer_start;
        for(int i=inner_start; i<inner_end; ++i) {
            bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+1+offset, j+1+offset, particles[i].dist(particles[j]), "coeff_bond_cross" ) );
            ++j;
        }
    }


    void generateCrossSpecial(vector<Atom> particles)
    {
        myFloat limit = scale * separation_other / fscale;

        for(int i=0; i<particles.size(); ++i)
        {
            if( particles[i].type == 10 || particles[i].type == 11 )
            {
                for(int j=0; j<particles.size(); ++j)
                {
                    if( ( particles[j].type == 1 || particles[j].type == 5 ) && particles[i].dist(particles[j]) < limit*0.9 ) {
                        bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_special" ) );
                    }
                }
            }
            if( particles[i].type == 9 )
            {
                for(int j=0; j<particles.size(); ++j)
                {
                    if( ( particles[j].type == 1 || particles[j].type == 5 )
                            && particles[i].dist(particles[j]) < limit*1.1
                            && particles[i].dist(particles[j]) > limit*0.8 )
                    {
                        bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_special" ) );
                    }
                }
            }
        }
    }

    void generateSameSpecial(vector<Atom> particles)
    {
        myFloat limit = scale * separation_same / fscale;

        for(int i=0; i<particles.size(); ++i)
        {
            if( particles[i].type == 3 )
            {
                for(int j=0; j<particles.size(); ++j)
                {
                    if( ( particles[j].type == 1 ) &&
                            particles[i].dist(particles[j]) > 1.0*limit  &&
                            particles[i].dist(particles[j]) < 1.8*limit )
                    {
                        bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_same_special" ) );
                    }
                }
            }
            if( particles[i].type == 7 )
            {
                for(int j=0; j<particles.size(); ++j)
                {
                    if( ( particles[j].type == 5 ) &&
                            particles[i].dist(particles[j]) > 1.0*limit  &&
                            particles[i].dist(particles[j]) < 1.6*limit )
                    {
                        bonds.push_back( Bond(bonds.size()+1, bonds.size()+1, i+offset+1, j+offset+1, particles[i].dist(particles[j]), "coeff_bond_same_special" ) );
                    }
                }
            }
            /*if( particles[i].type == 2 ) {
                for(int j=0; j<particles.size(); ++j) {
                    if( ( particles[j].type == 7 ) &&
                            particles[i].dist(particles[j]) > 1.0*limit  &&
                            particles[i].dist(particles[j]) < 2.5*limit ) {
                        bonds.push_back( Bond(bonds.size(), bonds.size(), i+1+offset, j+1+offset, particles[i].dist(particles[j]) ) );
                    }
                }
            }
            if( particles[i].type == 3 ) {
                for(int j=0; j<particles.size(); ++j) {
                    if( ( particles[j].type == 6 ) &&
                            particles[i].dist(particles[j]) > 1.0*limit  &&
                            particles[i].dist(particles[j]) < 2.5*limit ) {
                        bonds.push_back( Bond(bonds.size(), bonds.size(), i+1+offset, j+1+offset, particles[i].dist(particles[j]) ) );
                    }
                }
            }
            if( particles[i].type == 3 ) {
                for(int j=0; j<particles.size(); ++j) {
                    if( ( particles[j].type == 2 ) &&
                            particles[i].dist(particles[j]) > 1.0*limit  &&
                            particles[i].dist(particles[j]) < 2.0*limit ) {
                        bonds.push_back( Bond(bonds.size(), bonds.size(), i+1+offset, j+1+offset, particles[i].dist(particles[j]) ) );
                    }
                }
            }
            if( particles[i].type == 7 ) {
                for(int j=0; j<particles.size(); ++j) {
                    if( ( particles[j].type == 6 ) &&
                            particles[i].dist(particles[j]) > 1.0*limit  &&
                            particles[i].dist(particles[j]) < 2.0*limit ) {
                        bonds.push_back( Bond(bonds.size(), bonds.size(), i+1+offset, j+1+offset, particles[i].dist(particles[j]) ) );
                    }
                }
            }

            if( particles[i].type == 9 ) {
                for(int j=0; j<particles.size(); ++j) {
                    if( ( particles[j].type == 10 ) &&
                            particles[i].dist(particles[j]) > 1.0*limit  &&
                            particles[i].dist(particles[j]) < 2.0*limit ) {
                        bonds.push_back( Bond(bonds.size(), bonds.size(), i+1+offset, j+1+offset, particles[i].dist(particles[j]) ) );
                    }
                }
            }*/
        }
    }

    void init(int size)
    {
        vector<Atom> temp;

        shiftTop = 0.0;
        shiftMid = 0.0;
        shiftBot = 0.0;
        getInfo(size, 1.0, separation_same, separation_other);
        shiftTop = separation_other - sqrt(separation_other*separation_other - 0.25*separation_same*separation_same);
        getInfo(size, 1.0, separation_same_shifted, separation_other_shifted);

        this->fscale = pentaFrame(temp, 0, 1.0, type::BASIC, 0);
    }

    void info()
    {
        cerr << "Params:" << endl;
        cerr << "Beads separation on interface same type inner = " << separation_same * scale/fscale << " : " << separation_same / fscale << endl;
        cerr << "Beads separation on interface same type outer = " << separation_same * c*scale/fscale << " : " << separation_same * c / fscale << endl;

        cerr << "Shift interface same type Top = " << shiftTop * c*scale/fscale << endl;

        cerr << "Beads separation on interface different type inner = " << separation_other * scale/fscale << endl;
        cerr << "Beads separation on interface different type outer = " << separation_other * c*scale/fscale << endl;

        cerr << "Beads separation on interface different type shifted outer = " << separation_other_shifted << " : " << separation_other_shifted * c*scale/fscale << endl;

        cerr << "Beads radius should be inner(6,7,8,9,10) be = " << separation_other*0.5 << " : " << separation_other*0.5 * scale/fscale << endl;
        cerr << "Beads radius should be outer(1,2,3,4,5) be = " << separation_other*0.5 << " : " << separation_other*0.5 * c*scale/fscale << endl;
    }

    void random_Translate_Rotate(vector<Atom>& particles)
    {
        vector<Atom> temp;
        temp.resize( particles.size() );
        // Random translate
        int q=0;
        myFloat box_len = pbc.box.x / scale;

        bool clash = true;
        while(clash) {
            q++;
            Atom pos( 0.0001 * (rand()%10000) * box_len,
                        0.0001 * (rand()%10000) * box_len,
                        0.0001 * (rand()%10000) * box_len );

            Atom translate = particles[0] + ( pos - particles[0] );

            for(int i=0; i< particles.size(); ++i) {
                temp[i] = ( particles[i] + translate );
                temp[i].mol_tag = particles[i].mol_tag;
            }

            clusterRotate_random(temp, 180);

            for(int i=0; i< temp.size(); ++i) {
                temp[i] = pbc.usePBC( temp[i], scale );
            }

            clash = false;
            for(int i=0; i<temp.size(); ++i) {
                for(int j=0; j<beads.size(); ++j) {
                    if( ( temp[i].dist( beads[j] ) < 1.2 * sigma[temp[i].type][beads[j].type] / scale ) ) {
                        clash = true;
                    }
                }
            }
            cerr << q << " clash!" << endl;
        }
        swap(temp, particles);
    }

    void getInfo(int size, myFloat iscale, myFloat& separation_same, myFloat& separation_other)
    {
        Atom push;

        int a=0, b=2, c=8, d=9;

        vector<Atom> veca;
        vector<Atom> vecb;

        //
        // a,b,c,d vertices of triangle
        //
        Atom center = (edge[a] + edge[b] + edge[c]) * (iscale/3.0);
        Atom center2 = (edge[a] + edge[b] + edge[d]) * (iscale/3.0);
        Atom cc = (center + center2) * 0.5;

        myFloat dist = cc.dot(edge[a]) / edge[a].size();
        Atom edge_A_scaled = edge[a] * ( dist / cc.size() );

        Atom vec = cc - edge_A_scaled;
        vec *= 1.0/(vec.size());

        push = getSubPoint(edge_A_scaled, center, center2, size, 1, 1) - getSubPoint(edge_A_scaled, center, center2, size, 1, 0);
        separation_same = push.size();

        pentaFrame(veca, 0, iscale, type::INTERFACE_SHIFT, 0 );
        pentaFrame(vecb, 2, iscale, type::INTERFACE_SHIFT, 0 );

        myFloat min = 99999;

        for( auto& a : veca)
        {
            for( auto& b : vecb)
            {
                if(min > a.dist(b) )
                {
                    min = a.dist(b);
                }
            }
        }

        separation_other = min;
    }

    myFloat dist_from_plane(Atom point, Atom plane_normal, Atom point_of_plane)
    {
        myFloat d = -1.0 * ( plane_normal.x*point_of_plane.x + plane_normal.y*point_of_plane.y + plane_normal.z*point_of_plane.z );
        return  (plane_normal.x*point.x + plane_normal.y*point.y + plane_normal.z*point.z + d) / plane_normal.size();
    }

    Atom getSubPoint(Atom& a, Atom& b, Atom& c, int size, int i, int j)
    {
        Atom vecAB( b.x - a.x, b.y - a.y, b.z - a.z); // vector from a to b
        Atom vecBC( c.x - b.x, c.y - b.y, c.z - b.z); // vector from b to c
        size-=1;

        vecAB *= (1.0/size);
        vecBC *= (1.0/size);

        return Atom(a.x + i*vecAB.x + j*vecBC.x, a.y + i*vecAB.y + j*vecBC.y, a.z + i*vecAB.z + j*vecBC.z);
    }




    //
    // DEPRECATED STUFF
    //



    /**
     * @brief set_sigma_explicit - Deprecated
     */
    void set_sigma_explicit() {
        //
        // MARK 2 MODEL, Worse than MARK 1 with fluctuation
        //

        sigma_size=11;

        sigma[0][0] = 12.46; // outer body
        sigma[1][1] = separation_other* c*scale/fscale; // interface outer
        sigma[2][2] = separation_other* c*scale/fscale; // interface outer
        sigma[3][3] = separation_other* c*scale/fscale; // interface outer

        sigma[4][4] = 10.9; // inner body
        sigma[5][5] = separation_other*scale/fscale; // interface inner
        sigma[6][6] = separation_other*scale/fscale; // interface inner
        sigma[7][7] = separation_other*scale/fscale; // interface inner

        sigma[8][8] = separation_other*(1.0 + 0.5*(c-1.0))*scale/fscale; // interface middle
        sigma[9][9] = separation_other*(1.0 + 0.5*(c-1.0))*scale/fscale; // interface middle
        sigma[10][10] = separation_other*(1.0 + 0.5*(c-1.0))*scale/fscale; // interface middle

        mixing_rules();

        sigma[0][1] = sigma[1][0] = 12.46;
        sigma[0][2] = sigma[2][0] = 11.6;
        sigma[0][3] = sigma[3][0] = 12.46;
        sigma[0][4] = sigma[4][0] = 17.56;
        sigma[0][8] = sigma[8][0] = 14.73;
        sigma[0][9] = sigma[9][0] = 14.64;
        sigma[0][10] = sigma[10][0] = 14.4;

        sigma[1][3] = sigma[3][1] = 13.34;
        sigma[1][5] = sigma[5][1] = 17.5;
        sigma[1][7] = sigma[7][1] = 19.0;
        sigma[1][8] = sigma[8][1] = 18.33;
        sigma[1][9] = sigma[9][1] = 8.75;
        sigma[1][10] = sigma[10][1] = 13.6;

        sigma[2][3] = sigma[3][2] = 6.75;
        sigma[2][6] = sigma[6][2] = 17.54;
        sigma[2][7] = sigma[7][2] = 17.48;
        sigma[2][8] = sigma[8][2] = 8.77;
        sigma[2][10] = sigma[10][2] = 9.53;

        sigma[3][7] = sigma[7][3] = 18.42;
        sigma[3][8] = sigma[8][3] = 12.73;
        sigma[3][9] = sigma[9][3] = 17.72;
        sigma[3][10] = sigma[10][3] = 9.21;

        sigma[4][5] = sigma[5][4] = 10.9;
        sigma[4][6] = sigma[6][4] = 10.16;
        sigma[4][7] = sigma[7][4] = 10.9;
        sigma[4][8] = sigma[8][4] = 13.25;
        sigma[4][9] = sigma[9][4] = 13.02;
        sigma[4][10] = sigma[10][4] = 15.24;

        //sigma[5][6] = sigma[6][5] = 15.21; // need to be set to separation_other
        sigma[5][7] = sigma[7][5] = 11.67;
        sigma[5][8] = sigma[8][5] = 18.89;
        sigma[5][9] = sigma[9][5] = 8.74;
        sigma[5][10] = sigma[10][5] = 17.0;

        sigma[6][7] = sigma[7][6] = 5.91;
        sigma[6][8] = sigma[8][6] = 8.77;
        sigma[6][9] = sigma[9][6] = 17.39;
        sigma[6][10] = sigma[10][6] = 12.42;

        sigma[7][8] = sigma[8][7] = 9.37;
        sigma[7][9] = sigma[9][7] = 12.96;
        sigma[7][10] = sigma[10][7] = 9.21;

        //sigma[8][9] = sigma[9][8] = 16.3;  // need to be set to separation_other
        sigma[8][10] = sigma[10][8] = 6.33;

        sigma[9][10] = sigma[10][9] = 12.5;

        sigma_cosatt[1][2] = sigma_cosatt[2][1] = true;
        sigma_cosatt[5][6] = sigma_cosatt[6][5] = true;
        sigma_cosatt[8][9] = sigma_cosatt[9][8] = true;
    }
};

#endif // DODECAHEDRON_H
