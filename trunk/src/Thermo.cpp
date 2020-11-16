#include "Thermo.h"

Thermo::Thermo(Lattice* _Binding_Lattice, double _ti, double _tf, double _dt)
{
    this->Binding_Lattice = _Binding_Lattice;
    this->ti = _ti;
    this->tf = _tf;
    this->dt = _dt;

    double k=0.001985875;

    double Q;

    for (double kt=(k*ti); kt<=(k*tf); kt+=(k*dt)){
        vkt.push_back(kt);
        vt.push_back(kt/0.001985875);
        Q = Binding_Lattice->search_lattice3(kt);
        vQ.push_back(Q);
        vlnQ.push_back(log(Q));
    }

    for (unsigned i=0; i<vlnQ.size()-1; i++){
        dlnQdT.push_back((vlnQ[i+1]-vlnQ[i])/(vt[i+1]-vt[i]));
        U.push_back(vkt[i]*vt[i]*dlnQdT[i]);
    }

    for (unsigned i=0; i<U.size()-1; i++){
        dUdT.push_back((U[i+1]-U[i])/(vt[i+1]-vt[i]));
    }

    printf("#Thermo: %10.5s %10.5s %10s %10s %10s %10s %10s\n", "T", "ln(Q)", "dlnQ/dT", "U", "S", "F", "CV");
    for (unsigned i=0; i<vlnQ.size()-1; i++){
        S.push_back((k*vlnQ[i])+(U[i]/vt[i]));
        minusTS.push_back(((k*vlnQ[i])+(U[i]/vt[i]))*(-vt[i]));
        F.push_back(-vkt[i]*vlnQ[i]);
        Cv.push_back(dUdT[i]);
        printf("Thermo: %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", vt[i], vlnQ[i], dlnQdT[i] , (U[i]), (S[i]), (F[i]), Cv[i]);
    }

    Binding_Lattice->print_line();
    printf("*** Probabilities:                                                 ***\n");
    Binding_Lattice->print_line();

    vprob_energies.push_back(Binding_Lattice->lowest_energy);
    double dene=abs(Binding_Lattice->lowest_energy/10);
    for (int i=1; i<10; i++){
        vprob_energies.push_back((Binding_Lattice->lowest_energy+(i*dene)));
    }

    printf("#Probs: %4s ", "T");
    for (unsigned i=0; i<10; i++){
        printf("%4s ", string("p("+to_string(vprob_energies[i])+")").c_str());
    }
    printf("%4s \n", "p(0.0)");

    for (unsigned i=0; i<vQ.size()-1; i++){
        printf("Probs: %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n", vt[i], (exp(-vprob_energies[0]/vkt[i])/vQ[i]), (exp(-vprob_energies[1]/vkt[i])/vQ[i]), (exp(-vprob_energies[2]/vkt[i])/vQ[i]),
                (exp(-vprob_energies[3]/vkt[i])/vQ[i]), (exp(-vprob_energies[4]/vkt[i])/vQ[i]), (exp(-vprob_energies[5]/vkt[i])/vQ[i]), (exp(-vprob_energies[6]/vkt[i])/vQ[i]), (exp(-vprob_energies[7]/vkt[i])/vQ[i]),
                (exp(-vprob_energies[8]/vkt[i])/vQ[i]), (exp(-vprob_energies[9]/vkt[i])/vQ[i]), (exp(0/vkt[i])/vQ[i]));
    }
}
