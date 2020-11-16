#include "Thermo.h"

Thermo::Thermo(Lattice* _Binding_Lattice, double _ti, double _tf, double _dt)
{
    this->Binding_Lattice = _Binding_Lattice;
    this->ti = _ti;
    this->tf = _tf;
    this->dt = _dt;
    this->Run_Temp_Scan();
}

Thermo::Thermo(Lattice* _Binding_Lattice, double kt){
    this->Binding_Lattice = _Binding_Lattice;
    this->Single_Temp(kt);
}


int Thermo::Run_Temp_Scan(void){
    double k=0.001985875;
    double Q;

/* Scanning the partition function Q as a function of the temperature
 * Also storing ln (Q) in C++ vectors.
 */

    for (double kt=(k*this->ti); kt<=(k*this->tf); kt+=(k*this->dt)){
        vkt.push_back(kt);
        vt.push_back(kt/0.001985875);
        Q = Binding_Lattice->search_lattice(kt);
        vQ.push_back(Q);
        vlnQ.push_back(log(Q));
    }

/*
 * Computing d ln(Q) / dT and
 * U = kt^2* [ d ln(Q)/dT ]
 */

    for (unsigned i=0; i<vlnQ.size()-1; i++){
        dlnQdT.push_back((vlnQ[i+1]-vlnQ[i])/(vt[i+1]-vt[i]));
        U.push_back(vkt[i]*vt[i]*dlnQdT[i]);
    }

/*
 * Computing dU/dT
 */

    for (unsigned i=0; i<U.size()-1; i++){
        dUdT.push_back((U[i+1]-U[i])/(vt[i+1]-vt[i]));
    }

/*
 * Printing the results:
 * S = (k*lnQ) + (U/T)
 * F = -kT*ln(Q)
 * Cv = dU/dT
 */

    printf("#Thermo: %10.5s %10.5s %10s %10s %10s %10s %10s\n", "T", "ln(Q)", "dlnQ/dT", "U", "S", "F", "CV");
    for (unsigned i=0; i<vlnQ.size()-1; i++){
        S.push_back((k*vlnQ[i])+(U[i]/vt[i]));
        minusTS.push_back(((k*vlnQ[i])+(U[i]/vt[i]))*(-vt[i]));
        F.push_back(-vkt[i]*vlnQ[i]);
        Cv.push_back(dUdT[i]);
        printf("Thermo: %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", vt[i], vlnQ[i], dlnQdT[i] , (U[i]), (S[i]), (F[i]), Cv[i]);
    }

/*
 * Printing the probabilities...
 */

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

    return 0;
}

int Thermo::Single_Temp(double kt){
    double k=0.001985875;
    double Q = Binding_Lattice->search_lattice_mc(kt);
    double lnQ = log(Q);

    double U=0.0;
    for (unsigned i=0; i<Binding_Lattice->ligand_energies.size(); i++){
        U += Binding_Lattice->ligand_energies[i]*exp(-Binding_Lattice->ligand_energies[i]/kt) / Q;
    }
    double t = kt/k;
    this->single_S = k*lnQ + (U/t);
    this->single_F = -kt*lnQ;
    this->single_U = U;
    this->single_minusTS = -t*this->single_S;

    printf("#Thermo: %10.5s %10s %10s %10s %10s %10s\n", "T", "ln(Q)", "U", "S", "-TS", "F");
    printf("Thermo: %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", t, lnQ, this->single_U, this->single_S, -t*this->single_S ,this->single_F);

    return 0;

}

