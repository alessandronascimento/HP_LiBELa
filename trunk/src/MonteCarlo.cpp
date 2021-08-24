#include "MonteCarlo.h"

MonteCarlo::MonteCarlo(Lattice* _lattice, double T, int nsteps)
{
    this->lattice = _lattice;
    srand(rand());
    r = gsl_rng_alloc (gsl_rng_ranlxs2);
    double seed = 30.0 * (rand()/(RAND_MAX + 1.0));
    gsl_rng_set(r, seed);

     this->run_MC(0.001985875*T, nsteps);

    pose = new Lattice::Pose;

}

int MonteCarlo::run_MC(double kT, int nsteps){
    double score, new_score, p, rnumber;
    int accepted=0, step=0;
    int is_bound=0;
    long double sum_ene = 0.0L ;
    lattice->search_lattice_mc(kT);
    score = lattice->ligand_energies[0];
    printf("MC [ %10s ] %10s %10s %10s %24s\n", "step", "<score>", "%accept", "%bound" , "current configuration");
    lattice->print_line();
    while(accepted < nsteps){
        step++;
        this->sample();
        new_score = this->score_pose(pose);
        if (new_score < score){
            score = new_score;
            accepted++;
            sum_ene += score;
            if (score == lattice->lowest_energy){
                is_bound++;
            }
        }
        else {
            p = exp(-(new_score-score)/kT);
            rnumber = gsl_rng_uniform(r);
            if (p > rnumber){
                score = new_score;
                accepted++;
                sum_ene += score;
                if (score == lattice->lowest_energy){
                    is_bound++;
                }
            }
            else {
                sum_ene += score;
            }
        }
        if (step % 1000000 == 0){
            printf("MC [ %10d ] %10.2f %10.6f %10.6f %3d %3d %3d %3d %3d %3d\n", accepted, double(sum_ene/step), (100.*accepted/step), (100.*is_bound/accepted) , this->pose->ijk[0][0], this->pose->ijk[0][1], this->pose->ijk[1][0],
                    this->pose->ijk[1][1], this->pose->ijk[2][0], this->pose->ijk[2][1]);

        }
    }
    return 0;
}

void MonteCarlo::sample(void){
    int rnumber = gsl_rng_uniform_int(r, lattice->poses_found);
    this->pose = &(lattice->ligand_slots[rnumber]);
}

double MonteCarlo::score_pose(Lattice::Pose* pose){
    return (lattice->score_pair(pose));
}

