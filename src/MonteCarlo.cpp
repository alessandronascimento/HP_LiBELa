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
//    printf("MC [ %10s ] %10s %10s %10s %24s\n", "step", "<score>", "%accept", "%bound" , "current configuration");
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
        if (step % 10000 == 0){
            this->steps.push_back(accepted);
            this->scores.push_back(double(sum_ene/step));
            this->frac_accept.push_back(100.*accepted/step);
            this->frac_bound.push_back(100.*is_bound/accepted);
            this->distances.push_back(this->distance(this->pose, lattice->best_pose));
//            printf("MC [ %10d ] %10.2f %10.6f %10.6f %3d %3d %3d %3d %3d %3d\n", accepted, double(sum_ene/step), (100.*accepted/step), (100.*is_bound/accepted) , this->pose->ijk[0][0], this->pose->ijk[0][1], this->pose->ijk[1][0],
//                    this->pose->ijk[1][1], this->pose->ijk[2][0], this->pose->ijk[2][1]);
        }
    }
    return 0;
}

void MonteCarlo::sample(void){
    int rnumber = gsl_rng_uniform_int(r, lattice->poses_found);
    this->pose = &(lattice->ligand_slots[rnumber]);
}

double MonteCarlo::score_pose(Lattice::Pose* pose){
    return (lattice->score_pair_pointer(pose));
}

double MonteCarlo::distance(Lattice::Pose* Pose1, Lattice::Pose Pose2){
    double distance=0.0;
    if (Pose1->ijk.size() == Pose2.ijk.size()){
        for (unsigned i=0; i< Pose1->ijk.size(); i++){
            distance += pow(Pose1->ijk[i][0]-Pose2.ijk[i][0],2) + pow(Pose1->ijk[i][1]-Pose2.ijk[i][1],2);
        }
        distance = sqrt(distance);
    }
    else{
        distance=10000.0;
    }
        return(distance);
}


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


BOOST_PYTHON_MODULE(pyMonteCarlo)
{

/*    
    class_< vector<double> >("vectorDouble")
    .def(vector_indexing_suite<vector<double> >())
    ;
*/
    class_< vector<double> >("vectorInt")
    .def(vector_indexing_suite<vector<int> >())
    ;

    class_<MonteCarlo>("MonteCarlo", init<Lattice*, double, int >())
    .def("run_MC", & MonteCarlo::run_MC)
    .def("sample", & MonteCarlo::sample)
    .def("score_pose", & MonteCarlo::score_pose)
    .def("distance", & MonteCarlo::distance)
    .def_readwrite("r", & MonteCarlo::r)
    .def_readwrite("lattice", & MonteCarlo::lattice)
    .def_readwrite("pose", & MonteCarlo::pose)
    .def_readwrite("steps", & MonteCarlo::steps)
    .def_readwrite("scores", & MonteCarlo::scores)
    .def_readwrite("distances", & MonteCarlo::distances)        
    .def_readwrite("frac_accept", & MonteCarlo::frac_accept)
    .def_readwrite("frac_bound", & MonteCarlo::frac_bound)
    ;
}
