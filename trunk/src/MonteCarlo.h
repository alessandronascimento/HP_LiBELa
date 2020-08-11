#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "Lattice.h"
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

using namespace std;


class MonteCarlo
{
public:
    gsl_rng * r;
    MonteCarlo(Lattice* _lattice);
    Lattice* lattice;

    int run_MC(double kT);
    Lattice::Pose* sample(void);
    double score_pose(Lattice::Pose* pose);
};

#endif // MONTECARLO_H
