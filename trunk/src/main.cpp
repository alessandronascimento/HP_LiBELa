#include <iostream>
#include "Lattice.h"

using namespace std;

int main()
{
    Lattice* Binding_Lattice = new Lattice;
    Binding_Lattice->create_binding_lattice();
    Binding_Lattice->print_lattice();
    Binding_Lattice->search_lattice();

    delete Binding_Lattice;
    return 0;
}
