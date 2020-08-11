#include "Lattice.h"

Lattice::Lattice()
{

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
    int found_poses = 0;

    for (unsigned i=1; i<this->lattice.size()-1; i++){
        for (unsigned j=1; j< this->lattice[i].size()-1; j++){
            if (this->is_empty(i, j)){
                if (this->is_empty(i, j+1)){
                    if (this->is_empty(i+1, j+1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i, j+1, i+1, j+1);
                        found_poses++;
                    }
                    if (this->is_empty(i-1, j+1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i, j+1, i-1, j+1);
                        found_poses++;
                    }
                }
                if (this->is_empty(i+1, j)){
                    if (this->is_empty(i+1, j+1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i+1, j, i+1, j+1);
                        found_poses++;
                    }
                    if (this->is_empty(i+1, j-1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i+1, j, i+1, j-1);
                        found_poses++;
                    }
                }
                if (this->is_empty(i, j-1)){
                    if (this->is_empty(i-1, j-1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i, j-1, i-1, j-1);
                        found_poses++;
                    }
                    if (this->is_empty(i-1, j+1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i, j-1, i-1, j+1);
                        found_poses++;
                    }
                }
                if (this->is_empty(i-1, j)){
                    if (this->is_empty(i-1, j-1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i-1, j, i-1, j-1);
                        found_poses++;
                    }
                    if (this->is_empty(i-1, j+1)){
                        printf("Solution found at [%2d %2d] [%2d %2d] [%2d %2d]\n", i, j, i-1, j, i-1, j+1);
                        found_poses++;
                    }
                }
            }
        }
    }
    printf("Number of solutions found = %d\n", found_poses);
}

bool Lattice::is_empty(unsigned i, unsigned j){
    bool ret;
    if ((this->lattice[i][j] == 0) and (i < this->lattice.size()) and (j < this->lattice[i].size())){
        ret = true;
    }
    else {
        ret = false;
    }
    return (ret);
}
