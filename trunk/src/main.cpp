#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "Lattice.h"
#include "MonteCarlo.h"

using namespace std;

int main(int argc, char **argv)
{

    printf("**********************************************************************\n");
    printf("***                                                                ***\n");
    printf("***                          HP - LiBELA                           ***\n");
    printf("***                                                                ***\n");
    printf("***  A minimalist lattice model for ligand binding thermodynamics  ***\n");
    printf("***                                                                ***\n");
    printf("***         Written by Alessandro S. Nascimento - IFSC/USP         ***\n");
    printf("***                       asnascimento@ifsc.usp.br                 ***\n");
    printf("***                                                                ***\n");
    printf("***                            August/2020                         ***\n");
    printf("***                                                                ***\n");
    printf("**********************************************************************\n");

    vector<double> vkt;
    vector<double> vt;
    vector<double> vQ;
    vector<double> vlnQ;
    double Q;
    bool verbose = false;
    double k=0.001985875;
    int c;
    double epsilon=1.0, polar_epsilon=1.5;
    double ti=25., tf=3000., dt=25. ;

    if (argc < 5){
        printf("Usage: %s -t <ti> -f <tf> -d <dt> -p <polar_epsilon> -e <epsilon>\n", argv[0]);
        exit(1);
    }

    while ((c = getopt (argc, argv, "e:p:i:f:d:v")) != -1)
        switch (c)
        {
        case 'i':
            ti = (atof(optarg));
            break;
        case 'f':
            tf = (atof(optarg));
            break;
        case 'd':
            dt = (atof(optarg));
            break;
        case 'p':
            polar_epsilon = (atof(optarg));
            break;
        case 'e':
            epsilon = (atof(optarg));
            break;
        case 'v':
            verbose = true;
            break;
        case '?':
            if (optopt == 'c') {
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                printf("Usage: %s -t <ti> -f <tf> -d <dt> -p <polar_epsilon> -e <epsilon>\n", argv[0]);
            }
            else if (isprint (optopt)) {
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                printf("Usage: %s -t <ti> -f <tf> -d <dt> -p <polar_epsilon> -e <epsilon>\n", argv[0]);
            }
            else {
                fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                printf("Usage: %s -t <ti> -f <tf> -d <dt> -p <polar_epsilon> -e <epsilon>\n", argv[0]);
            }
            return 1;
        default:
            abort ();
        }

    printf("***                                                                ***\n");
    printf("*** Parameters:                                                    ***\n");
    printf("*** Ti = %5.3f                                                    ***\n", ti);
    printf("*** Tf = %5.3f                                                  ***\n", tf);
    printf("*** dT = %5.3f                                                    ***\n", dt);
    printf("*** <epsilon> = %5.3f                                             ***\n", epsilon);
    printf("*** <epsilon_polar> = %5.3f                                       ***\n", polar_epsilon);
    printf("***                                                                ***\n");

    Lattice* Binding_Lattice = new Lattice(epsilon, polar_epsilon, verbose);

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    Binding_Lattice->create_binding_lattice();

    Binding_Lattice->print_lattice();

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    for (double kt=(0.001985875*ti); kt<=(0.001985875*tf); kt+=(0.001985875*dt)){
        vkt.push_back(kt);
        vt.push_back(kt/0.001985875);
        Q = Binding_Lattice->search_lattice2(kt);
        vQ.push_back(Q);
        vlnQ.push_back(log(Q));
    }

    vector<double> dlnQdT;
    vector<double> U;
    for (unsigned i=0; i<vlnQ.size()-1; i++){
        dlnQdT.push_back((vlnQ[i+1]-vlnQ[i])/(vt[i+1]-vt[i]));
        U.push_back(vkt[i]*vt[i]*dlnQdT[i]);
    }

    vector<double> dUdT;
    for (unsigned i=0; i<U.size()-1; i++){
        dUdT.push_back((U[i+1]-U[i])/(vt[i+1]-vt[i]));
    }

    printf ("#Thermo: %10.5s %10.5s %10s %10s %10s %10s %10s\n", "T", "ln(Q)", "dlnQ/dT", "U", "S", "F", "CV");
    for (unsigned i=0; i<vlnQ.size()-1; i++){
        printf("Thermo: %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", vt[i], vlnQ[i], dlnQdT[i] , (k*vt[i]*vt[i]*dlnQdT[i]), ((k*vlnQ[i])+(vkt[i]*dlnQdT[i])), (-vkt[i]*vlnQ[i]), dUdT[i]);
    }

    vector<double> vprob_energies;
    vprob_energies.push_back(-Binding_Lattice->lowest_energy);
    double dene=abs(Binding_Lattice->lowest_energy/10);
    for (int i=1; i<10; i++){
        vprob_energies.push_back(-(Binding_Lattice->lowest_energy+(i*dene)));
    }

    printf("#Probs: %4s ", "T");
    for (unsigned i=0; i<10; i++){
        printf("%4s ", string("p("+to_string(vprob_energies[i])+")").c_str());
    }
    printf("%4s \n", "p(0.0)");

    for (unsigned i=0; i<vQ.size()-1; i++){
        printf("Probs: %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n", vt[i], (exp(vprob_energies[0]/vkt[i])/vQ[i]), (exp(vprob_energies[1]/vkt[i])/vQ[i]), (exp(vprob_energies[2]/vkt[i])/vQ[i]),
                (exp(vprob_energies[3]/vkt[i])/vQ[i]), (exp(vprob_energies[4]/vkt[i])/vQ[i]), (exp(vprob_energies[5]/vkt[i])/vQ[i]), (exp(vprob_energies[6]/vkt[i])/vQ[i]), (exp(vprob_energies[7]/vkt[i])/vQ[i]),
                (exp(vprob_energies[8]/vkt[i])/vQ[i]), (exp(vprob_energies[9]/vkt[i])/vQ[i]), (exp(0/vkt[i])/vQ[i]));
    }

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

//    MonteCarlo* MC = new MonteCarlo(Binding_Lattice);
//    delete MC;

    delete Binding_Lattice;

    printf("***                                                                ***\n");
    printf("***                                                                ***\n");
    printf("**********************************************************************\n");

    return 0;
}

