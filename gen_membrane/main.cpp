#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <string>

#include "membrane.h"
#include "flatmembrane.h"
#include "vesicle.h"
#include "input.h"


using namespace std;



int main(int argc, char* argv[])
{
    Input in;
    if(argc < 2 ) {
        in.helpMessage();
        return 1;
    }

    for(int i=1; i<argc; ++i)
    {
        in.loadInput(argv[i]);
    }

    Membrane* mem;

    if(in.op == 1) { // trim and multiply
        mem = new FlatMembrane(in.in_file);
        mem->trim(in.trim, in.multiple); // size of edges to cut off
    }

    if(in.op == 2) {
        mem = new Vesicle();
        mem->generate(in);
    }

    if(in.op == 3) {  // copy membrane in z dir
        mem = new Membrane(in.in_file);
        mem->copyZ(in.trim);
    }

    if(in.op == 4) { // populate receptors
        mem = new Membrane(in.in_file);
        mem->populate_receptors(in.trim, in.multiple);
    }

    if(in.op == 5) { // analyze
        mem = new Membrane(in.in_file);
        mem->analyze(in.out_file);
    }

    if(in.op == 6) { // analyze
        mem = new Membrane(in.in_file);
        mem->center(in.out_file);
    }



    if(in.out_type == 1) {
        mem->printOutLammps(in);
    }
    if(in.out_type == 2) {
        mem->printOutSC();
    }

    delete mem;

    return 0;
}





