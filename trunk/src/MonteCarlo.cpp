#include "MonteCarlo.h"

MonteCarlo::MonteCarlo(Lattice* _lattice)
{
    this->lattice = _lattice;
    srand(rand());
    r = gsl_rng_alloc (gsl_rng_ranlxs2);
    double seed = 30.0 * (rand()/(RAND_MAX + 1.0));
    gsl_rng_set(r, seed);

     this->run_MC(2.5);

}

int MonteCarlo::run_MC(double kT){
    Lattice::Pose* pose;
    double score, new_score, p, rnumber;
    int nsteps=1000000;
    int accepted=0, step=0;
    score = lattice->ligand_energies[0];
    while(accepted < nsteps){
        pose = this->sample();
        new_score = this->score_pose(pose);
        if (new_score < score){
            score = new_score;
            accepted++;
        }
        else {
            p = exp(-(new_score-score)/kT);
            rnumber = gsl_rng_uniform(r);
            if (p > rnumber){
                score = new_score;
                accepted++;
            }
        }
        step++;
        if (step % 1000 == 0){
            printf("%10d %10.2f %10.6f\n", step, score, (100.*accepted/step));
        }
    }
    return 0;

}

Lattice::Pose* MonteCarlo::sample(){
    Lattice::Pose* pose = new Lattice::Pose;
    pose = &(lattice->ligand_slots[gsl_rng_uniform_int(r, this->lattice->ligand_slots.size())]);
    return (pose);
}

double MonteCarlo::score_pose(Lattice::Pose* pose){
    return (lattice->score_pair(pose));
}

