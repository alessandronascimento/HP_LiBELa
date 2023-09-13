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
    vector<int> steps;
    vector<double> scores;
    vector<double> frac_accept;
    vector<double> frac_bound;
    vector<double> distances;

    int run_MC(double kT, int nsteps);
    void sample(void);
    double score_pose(Lattice::Pose* pose);
    double distance(Lattice::Pose* Pose1, Lattice::Pose Pose2);
};

#endif // MONTECARLO_H
