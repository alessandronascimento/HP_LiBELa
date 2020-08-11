#ifndef LATTICE_H
#define LATTICE_H

#include<vector>

using namespace std;


class Lattice
{
public:
    struct  coord{
        int x;
        int y;
    };

    struct Pose {
        int n;
        vector<vector<int> > ijk;
    };
    vector<int> ligand_types;
    vector<vector<int> > lattice;
    double epsilon;
    double epsilon_polar;
    Lattice();
    void create_lattice(int m, int n);
    void create_binding_lattice(void);
    void print_lattice(void);
    void search_lattice(void);
    bool is_empty(int i, int j);
    vector<coord> find_contacts(int i, int j);
    bool is_occupied(int i, int j);
    double score_pair(Pose binding_pose);
    void print_line(void);
};

#endif // LATTICE_H
