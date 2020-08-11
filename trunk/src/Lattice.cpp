#include "Lattice.h"

Lattice::Lattice()
{
    this->epsilon = -1.0;
    this->epsilon_polar = -1.5;
}


void Lattice::create_lattice(int m, int n){
    vector<int> v;
    this->lattice.clear();
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            v.push_back(0);
        }
        this->lattice.push_back(v);
        v.clear();
    }
}

void Lattice::create_binding_lattice(void){

    // create a 8 x 8 lattice

    this->create_lattice(8,8);

/* fill it with some predefined pattern. Here +2 describes a positively charge
 * and -2 describes a negative charge, 1 means occupied and hydrophobic site
 * 0 means empyt site.
 */

    for (unsigned i=0; i<4; i++){
        for (unsigned j=0; j<8; j++){
            this->lattice[j][i] = 1;
        }
    }
    this->lattice[3][3] = 0;
    this->lattice[1][4] = -2;
    this->lattice[4][4] = 2;
    this->ligand_types.push_back(0);
    this->ligand_types.push_back(2);
    this->ligand_types.push_back(-2);

}

void Lattice::print_lattice(void){
    int lines = int(this->lattice.size());
    int cols = int(this->lattice[0].size());

    printf("%2s ", "#");

    for (int i=0; i<cols; i++){
        printf("%2d ", i);
    }
    printf("\n");

    for (int i=0; i<lines; i++){
        printf("%2d ", i);
        for (int j=0; j<cols; j++){
            switch (this->lattice[i][j]) {
            case 0:
                printf("%2s ", "-");
                break;
            case 1:
                printf("%2s ", "*");
                break;
            case -2:
                printf("\033[0;31m");
                printf("%2s ", "*");
                printf("\033[0m");
                break;
            case 2:
                printf("\033[1;34m");
                printf("%2s ", "*");
                printf("\033[0m");
                break;
            }
        }
        printf("\n");
    }
}

void Lattice::search_lattice(void){
    double score;
    vector<Pose> poses_found;
    vector<int> tmp(2);

    for (int i=1; i<this->lattice.size()-1; i++){
        for (int j=1; j< this->lattice[i].size()-1; j++){
            if (this->is_empty(i, j)){
                if (this->is_empty(i, j+1)){
                    if (this->is_empty(i+1, j+1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j+1);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j+1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j+1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j+1);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        score = this->score_pair(inv_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);
                    }
                    if (this->is_empty(i-1, j+1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j+1);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j+1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        poses_found.push_back(this_pose);
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j+1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j+1);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        poses_found.push_back(inv_pose);
                        score = this->score_pair(inv_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);
                    }
                }
                if (this->is_empty(i+1, j)){
                    if (this->is_empty(i+1, j+1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j+1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        poses_found.push_back(this_pose);
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j+1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        poses_found.push_back(inv_pose);
                        score = this->score_pair(inv_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);
                    }
                    if (this->is_empty(i+1, j-1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j-1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        poses_found.push_back(this_pose);
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j-1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i+1);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        poses_found.push_back(inv_pose);
                        score = this->score_pair(inv_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);
                    }
                }
                if (this->is_empty(i, j-1)){
                    if (this->is_empty(i-1, j-1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j-1);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j-1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        poses_found.push_back(this_pose);
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j-1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j-1);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        poses_found.push_back(inv_pose);
                        score = this->score_pair(inv_pose);printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);
                    }
                    if (this->is_empty(i-1, j+1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j-1);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j+1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        poses_found.push_back(this_pose);
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j+1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i);
                        tmp[1] = (j-1);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        poses_found.push_back(inv_pose);
                        score = this->score_pair(inv_pose);
                        score = this->score_pair(inv_pose);printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);
                    }
                }
                if (this->is_empty(i-1, j)){
                    if (this->is_empty(i-1, j-1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j-1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        poses_found.push_back(this_pose);
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j-1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        poses_found.push_back(inv_pose);
                        score = this->score_pair(inv_pose);
                        score = this->score_pair(inv_pose);printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);

                    }
                    if (this->is_empty(i-1, j+1)){
                        Pose this_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j);
                        this_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j+1);
                        this_pose.ijk.push_back(tmp);
                        poses_found.push_back(this_pose);
                        this_pose.n=3;
                        poses_found.push_back(this_pose);
                        score = this->score_pair(this_pose);
                        printf("%2d %2d %2d %2d %2d %2d %8.2f\n", this_pose.ijk[0][0], this_pose.ijk[0][1], this_pose.ijk[1][0], this_pose.ijk[1][1],
                                this_pose.ijk[2][0], this_pose.ijk[2][1], score);
                        Pose inv_pose;
                        tmp[0] = (i);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j+1);
                        inv_pose.ijk.push_back(tmp);
                        tmp[0] = (i-1);
                        tmp[1] = (j);
                        inv_pose.ijk.push_back(tmp);
                        poses_found.push_back(inv_pose);
                        inv_pose.n=3;
                        poses_found.push_back(inv_pose);
                        score = this->score_pair(inv_pose);
                        score = this->score_pair(inv_pose);printf("%2d %2d %2d %2d %2d %2d %8.2f\n", inv_pose.ijk[0][0], inv_pose.ijk[0][1], inv_pose.ijk[1][0], inv_pose.ijk[1][1],
                                inv_pose.ijk[2][0], inv_pose.ijk[2][1], score);
                    }
                }
            }
        }
    }
    printf("Number of poses found: %5d.\n", int(poses_found.size()));
}

bool Lattice::is_empty(int i, int j){
    bool ret;
    if ((this->lattice[i][j] == 0) and (i < this->lattice.size()) and (j < this->lattice[i].size())){
        ret = true;
    }
    else {
        ret = false;
    }
    return (ret);
}

bool Lattice::is_occupied(int i, int j){
    bool ret = false;
    if ( (i < this->lattice.size()) and (i >= 0) and (j < this->lattice[i].size()) and (j>=0)){
        if (this->lattice[i][j] != 0){
            ret = true;
        }
    }
    return (ret);
}

vector<Lattice::coord> Lattice::find_contacts(int i, int j){
    vector<coord> coords;
    for(int l=i-1; l <= i+1; l++){
        for (int c=j-1; c<=j+1; c++){
            if (this->is_occupied(l, c)){
                coord tmp;
                tmp.x=l;
                tmp.y=c;
                coords.push_back(tmp);
            }
        }
    }
    return (coords);
}

double Lattice::score_pair(Pose binding_pose){
    double score=0.0;
    int type_R, type_L;

    for (unsigned atom=0; atom < binding_pose.ijk.size(); atom++){
        vector<coord> coords = this->find_contacts(binding_pose.ijk[atom][0], binding_pose.ijk[atom][1]);
        for (unsigned i=0; i<coords.size(); i++){
            type_R = this->lattice[coords[i].x][coords[i].y];
            type_L = this->ligand_types[atom];

            switch (type_R) {
            case 1:
                score += epsilon;
                break;
            case 2:
                if (type_L == -2){
                    score += epsilon_polar;
                }
                break;
            case -2:
                if (type_L == 2){
                    score += epsilon_polar;
                }
                break;
            }
        }
    }
    return (score);
}


void Lattice::print_line(void){
printf("*************************\n");
}
