#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

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
    double kt;
    bool verbose = false;
    vector<Pose> ligand_slots;
    vector<double> ligand_energies;
    vector<int> ligand_types;
    vector<vector<int> > lattice;
    double epsilon;
    double epsilon_polar;
    double lowest_energy;
    int poses_found;
    Lattice();
    Lattice(double epsilon, double polar, bool _verbose);
    void create_lattice(int m, int n);
    void create_binding_lattice(void);
    void create_ligand(void);
    void create_empty_binding_lattice(void);
    void read_lattice_from_file(string filename);
    void print_lattice(void);
    void print_ligand(void);
    double search_lattice(double this_kt);
    vector<double> find_lowest(vector<double> v);
    double search_lattice_mc(double kt);
    bool is_empty(int i, int j);
    vector<coord> find_contacts(int i, int j);
    bool is_occupied(int i, int j);
    double score_pair(Pose binding_pose);
    double score_pair(Pose* binding_pose);
    void print_line(void);
    bool is_triangle(int i, int j, int k, int l, int m, int n);
    bool exists(int i, int j);

    double search_triangle_a(int i, int j);
    double search_triangle_b(int i, int j);
    double search_triangle_c(int i, int j);
    double search_triangle_d(int i, int j);
    double search_triangle_e(int i, int j);
    double search_triangle_f(int i, int j);
    double search_triangle_g(int i, int j);
    double search_triangle_h(int i, int j);
    double search_triangle_i(int i, int j);
    double search_triangle_j(int i, int j);
    double search_triangle_k(int i, int j);
    double search_triangle_l(int i, int j);
};

#endif // LATTICE_H
