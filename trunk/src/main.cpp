//#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "Lattice.h"
#include "MonteCarlo.h"
#include "Thermo.h"

//using namespace std;

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

    bool verbose = false;
    int c;
    bool mc = false;
    int mc_steps = 1000000;
    bool single_temp = false;
    double temperature = 300.;

    double epsilon=1.0, polar_epsilon=1.5;

    double ti=25., tf=3000., dt=25. ;

    if (argc < 5){
        printf("Usage: %s -i <initial temperature> -f <final temperature> -d <dT> -p <polar_epsilon> -e <epsilon> -t <temperature> [ -m <mc_steps> ] [-v ]\n ", argv[0]);
        exit(1);
    }

    while ((c = getopt (argc, argv, "e:p:i:f:d:m:t:v")) != -1)
        switch (c)
        {
        case 'i':
            ti = (atof(optarg));
            single_temp = false;
            break;
        case 'f':
            tf = (atof(optarg));
            single_temp = false;
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
        case 't':
            temperature = double(atof(optarg));
            single_temp = true;
            break;
        case 'm':
            mc = true;
            mc_steps = (atoi(optarg));
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
    printf("*** T  = %5.3f                                                    ***\n", temperature);
    printf("*** <epsilon> = %5.3f                                             ***\n", epsilon);
    printf("*** <epsilon_polar> = %5.3f                                       ***\n", polar_epsilon);
    printf("***                                                                ***\n");

    Lattice* Binding_Lattice = new Lattice(epsilon, polar_epsilon, verbose);

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    Binding_Lattice->create_binding_lattice();
    //    Binding_Lattice->create_empty_binding_lattice();

    Binding_Lattice->print_lattice();

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    Thermo* ThermoData;

    if (single_temp){
        ThermoData = new Thermo(Binding_Lattice, (0.001985875*temperature));
    }
    else {
        ThermoData = new Thermo(Binding_Lattice, ti, tf, dt);
    }

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    Lattice* Empty_Lattice = new Lattice(epsilon, polar_epsilon, verbose);

    Empty_Lattice->print_line();
    Empty_Lattice->print_line();

    Empty_Lattice->create_empty_binding_lattice();

    Empty_Lattice->print_lattice();

    Empty_Lattice->print_line();
    Empty_Lattice->print_line();

    Thermo* ThermoData2;

    if (single_temp){
        ThermoData2 = new Thermo(Empty_Lattice, (0.001985875*temperature));
    }
    else {
        ThermoData2 = new Thermo(Empty_Lattice, ti, tf, dt);
    }

    Empty_Lattice->print_line();
    Empty_Lattice->print_line();

    double dF, dU, dS;

    if (! single_temp){
        printf("#Binding_Data: %10s %10s %10s %10s %10s\n", "Temp(K)", "dF", "dU", "dS", "-TdS");
        for (unsigned i=0; i<ThermoData->vlnQ.size()-1; i++){
            dF = ThermoData->F[i] - ThermoData2->F[i];
            dU = ThermoData->U[i] - ThermoData2->U[i];
            dS = ThermoData->S[i] - ThermoData2->S[i];
            printf("Binding_Data: %10.5f %10.5f %10.5f %10.5f %10.5f\n", ThermoData->vt[i], dF, dU, dS, -ThermoData->vt[i]*dS);
        }
        Binding_Lattice->print_line();
        Binding_Lattice->print_line();
    }
    else {
        printf("#Binding_Data: %10s %10s %10s %10s %10s\n", "Temp(K)", "dF", "dU", "dS", "-TdS");
        dF = ThermoData->single_F - ThermoData2->single_F;
        dU = ThermoData->single_U - ThermoData2->single_U;
        dS = ThermoData->single_S - ThermoData2->single_S;
        printf("Binding_Data: %10.5f %10.5f %10.5f %10.5f %10.5f\n", temperature, dF, dU, dS, -temperature*dS);
        Binding_Lattice->print_line();
        Binding_Lattice->print_line();
    }

    if (mc){
        MonteCarlo* MC = new MonteCarlo(Binding_Lattice, 300.0, mc_steps);
        delete MC;
    }

    delete ThermoData;
    delete ThermoData2;
    delete Binding_Lattice;
    delete Empty_Lattice;

    printf("***                                                                ***\n");
    printf("***                                                                ***\n");
    printf("**********************************************************************\n");

    return 0;
}

