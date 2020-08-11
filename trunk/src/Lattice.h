#ifndef LATTICE_H
#define LATTICE_H

#include<vector>

using namespace std;


class Lattice
{
public:
    vector<vector<int> > lattice;
    Lattice();
    void create_lattice(int m, int n);
    void create_binding_lattice(void);
    void print_lattice(void);
    void search_lattice(void);
    bool is_empty(unsigned i, unsigned j);
};

#endif // LATTICE_H
