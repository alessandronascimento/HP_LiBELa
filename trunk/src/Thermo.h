#include "Lattice.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>

#ifndef THERMO_H
#define THERMO_H

using namespace std;


class Thermo
{
public:
    Thermo(Lattice* _Binding_Lattice, double _ti, double _tf, double _dt);
    Lattice* Binding_Lattice;
    double ti, tf, dt;
    vector<double> vkt;
    vector<double> vt;
    vector<double> vQ;
    vector<double> vlnQ;
    vector<double> dlnQdT;
    vector<double> U;
    vector<double> dUdT;
    vector<double> vprob_energies;
    vector<double> F;
    vector<double> S;
    vector<double> minusTS;
    vector<double> Cv;
};

#endif // THERMO_H
