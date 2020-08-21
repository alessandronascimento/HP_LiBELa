#include <iostream>
#include "Lattice.h"
#include "MonteCarlo.h"

using namespace std;

int main()
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


    Lattice* Binding_Lattice = new Lattice;
    vector<double> vkt;
    vector<double> vt;
    vector<double> vQ;
    vector<double> vlnQ;
    double Q;
    double k=0.001985875;

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    Binding_Lattice->create_binding_lattice();

    Binding_Lattice->print_lattice();

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    for (double kt=0.05; kt<=3.0; kt+=0.01){
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

    printf("#Probs: %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "T", "p(-8.0)", "p(-5.0)", "p(-4.5)", "p(-3.5)", "p(-3.0)", "p(-2.0)", "p(-1.5)", "p(0.0)");
    for (unsigned i=0; i<vQ.size()-1; i++){
        printf("Probs: %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", vt[i], (exp(8/vkt[i])/vQ[i]), (exp(5/vkt[i])/vQ[i]), (exp(4.5/vkt[i])/vQ[i]), (exp(3.5/vkt[i])/vQ[i]), (exp(3/vkt[i])/vQ[i]), (exp(2/vkt[i])/vQ[i]),
               (exp(1.5/vkt[i])/vQ[i]), (exp(0/vkt[i])/vQ[i]));
    }

    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

//    MonteCarlo* MC = new MonteCarlo(Binding_Lattice);

//    delete MC;
    delete Binding_Lattice;
    return 0;
}

