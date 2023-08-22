#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "Lattice.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

using namespace std;


class MonteCarlo
{
public:
    gsl_rng * r;
    MonteCarlo(Lattice* _lattice, double T, int nsteps);
    Lattice* lattice;
    Lattice::Pose* pose;

    int run_MC(double kT, int nsteps);
    void sample(void);
    double score_pose(Lattice::Pose* pose);
};

#endif // MONTECARLO_H
