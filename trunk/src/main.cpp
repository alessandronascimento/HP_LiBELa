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
    Binding_Lattice->print_line();
    Binding_Lattice->print_line();
    Binding_Lattice->create_binding_lattice();
    Binding_Lattice->print_lattice();
    Binding_Lattice->print_line();
    Binding_Lattice->print_line();
    Binding_Lattice->search_lattice();
    Binding_Lattice->print_line();
    Binding_Lattice->print_line();

    MonteCarlo* MC = new MonteCarlo(Binding_Lattice);

    delete MC;
    delete Binding_Lattice;
    return 0;
}
