#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
# include <cstdlib>
#include <ctime>
# include "mpi.h"


using namespace std;

class Bound
{
    int endRow;
    int endColunm;
    int endDepth;
    int startRow;
    int startColunm;
    int startDepth;
public:
    int getStartRow (int idBlock, int nbProcessor, int n) {
        startRow = idBlock*(n/nbProcessor);
        return startRow;
    }

    int getStartColunm (int idBlock, int nbProcessor, int n) {
        startColunm = idBlock*(n/nbProcessor);
        return startColunm;
    }

    int getStartDepth (int idBlock, int nbProcessor, int n) {
        startDepth = idBlock*(n/nbProcessor);
        return startDepth;
    }

    int getEndRow (int nbProcessor, int n) {
        endRow = startRow + (n/nbProcessor)-1;
        return endRow;
    }

    int getEndColunm (int nbProcessor, int n) {
        endColunm = startColunm + (n/nbProcessor)-1;
        return endColunm;
    }

    int getEndDepth (int nbProcessor, int n) {
        endDepth = startDepth + (n/nbProcessor)-1;
        return endDepth;
    }

};


class Block
{
public:
    int id;
    int idPlan;
    
    void buildBound (int idBlock, int nbProcessor, int str_size_1, int str_size_2, int nb_str_1) {
        Bound currentBlockBound;
        currentBlockBound.getStartRow(idBlock, nbProcessor, str_size_1);
        currentBlockBound.getStartColunm(idBlock, nbProcessor, nb_str_1);
        currentBlockBound.getStartDepth(idBlock, nbProcessor, str_size_2);

        currentBlockBound.getEndRow(nbProcessor, str_size_1);
        currentBlockBound.getEndColunm(nbProcessor, nb_str_1);
        currentBlockBound.getEndDepth(nbProcessor, str_size_2);
    }
};


class ParallelPlan
{
public:
    int id;
    int row;
    int colunm;
    int rank;
    int nextRank;

    void buildBock(int row, int colunm, int idPlan, int nbProcessor, int str_size_1, int str_size_2, int nb_str_1) {
        Block block;
        for(int i=0; i<row; i++) {
            for(int j=0; j<colunm; j++) {
                block.id = idPlan;
                block.buildBound(block.id, nbProcessor, str_size_1, str_size_2, nb_str_1);
            }
        }
    }
};


class Cube
{
public:
    void buildCube (int nbProcessor, int str_size_1, int str_size_2, int nb_str_1) {

        for(int i=0; i<nbProcessor; i++) {
            ParallelPlan plan;
            plan.id = i;
            plan.row = 1;
            plan.colunm = 1;
            plan.rank = i;

            if(nbProcessor != i+1) {
                plan.nextRank = i+1;
            } else {
                plan.nextRank = -1;
            }

            for(int a=0; a<nbProcessor; a++) {
                for (int b=0; b<nbProcessor; b++) {
                    plan.buildBock(a, b, i, nbProcessor, str_size_1, str_size_2, nb_str_1);
                 }
            }
        }
    }

};


class Dag {
public:
    Cube* buildDag (int nbProcessor, int str_size_1, int str_size_2, int nb_str_1, int nb_str_2) {
        int nbCube = nb_str_2;
        Cube cubeTab[nbCube];
        for (int i=0; i<= nbCube; i++) {
            cubeTab[i].buildCube(nbProcessor, str_size_1,str_size_2, nb_str_1);
        }
        return cubeTab;
    }
};


// Classique sequences edit distance method
int classic_seq_edit_distance(const string &a, const string &b)
{
    int m = a.size();
    int n = b.size();

    int **scores = new int *[m + 1];
    scores[0] = new int[(m + 1) * (n + 1)];
    for (int i = 1; i < m + 1; i++)
    {
        scores[i] = scores[0] + i * (n + 1);
    }
    

    for (int i = 0; i <= m; ++i)
        scores[i][0] = i;

    for (int j = 0; j <= n; ++j)
        scores[0][j] = j;

    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            if (a[i - 1] == b[j - 1])
            {
                scores[i][j] = scores[i - 1][j - 1];
                continue;
            }

            int x = scores[i - 1][j] + 1;     // Deletion
            int y = scores[i][j - 1] + 1;     // Insertion
            int z = scores[i - 1][j - 1] + 1; // Substitution

            int t = min(x, y);

            scores[i][j] = min(t, z);
        }
    }
    int result =  scores[m][n];

    for (int i= 0; i < m+1; i++){
        delete scores[i];
    }
    delete[] scores;
    return result;
}

// Building patterns called "motif"
 string* first_pattern_construction (int str_size_1, int nb_str_1) {
      ifstream datasets("input/datasets.txt");
      string* motif_1;
      string str1;
      motif_1 = new string[nb_str_1];
      if (datasets) {
         string dataset;
         datasets >> dataset;
         ifstream input(dataset.c_str());
         if (input) {
             for (int i = 0; i < nb_str_1; i++) {
                 input >> str1;
                 if (str_size_1 > str1.size()) {
                     if (str_size_1 > str1.size()) {
                         cout << "[ERROR]the max size of the first string is" << str1.size() << "in the dataset" << dataset << endl;
                     }

                     exit(0);
                 }
                 str1 = str1.substr(0, str_size_1);
                 motif_1[i] = str1;
             }

             input.close();
         }
         datasets.close();
     }

     return motif_1;

}

 string* second_pattern_construction (int str_size_2, int nb_str_2) {
      ifstream datasets("input/datasets.txt");
      string* motif_2;
      string str2;
      motif_2 = new string[nb_str_2];
      if (datasets) {
         string dataset;
         datasets >> dataset;
         ifstream input(dataset.c_str());
         if (input) {
             for (int i = 0; i < nb_str_2; i++) {
                 input >> str2;
                 if (str_size_2 > str2.size()) {
                     if (str_size_2 > str2.size()) {
                         cout << "[ERROR]the max size of the first string is" << str2.size() << "in the dataset" << dataset << endl;
                    }
                     exit(0);
                }

                 str2 = str2.substr(0, str_size_2);
                 motif_2[i] = str2;
             }

             input.close();
         }
         datasets.close();
     }

     return motif_2;
}
// Filling of Dr table
 int** dr_matrice_func(int str_size_1, int nb_str_1){
     int** dr_matrice;
     dr_matrice = new int* [str_size_1];
     for (int i= 0; i <str_size_1; i++){
         dr_matrice[i] = new int[nb_str_1];
     }

     for (int i = 0; i < str_size_1; i++) {
         for (int j = 0; j < nb_str_1; j++) {
             dr_matrice[i][j] = (j + 1);
         }
     }
   return dr_matrice;

}

// Filling of Dc table
int** dc_matrice_func(int str_size_1, int nb_str_1){
    int** dc_matrice;
    dc_matrice = new int* [str_size_1];
    for (int i= 0; i <str_size_1; i++){
        dc_matrice[i] = new int[nb_str_1];
    }

    for (int i = 0; i < str_size_1; i++) {
         for (int j = 0; j <  nb_str_1; j++) {
             dc_matrice[i][j] = (i + 1);
         }
     }
    return dc_matrice;
}

// Filling of Ir table
 int** ir_matrice_func(int str_size_2, int nb_str_2){
     int** ir_matrice;
     ir_matrice = new int* [str_size_2];
     for (int i= 0; i <str_size_2; i++){
         ir_matrice[i] = new int[nb_str_2];
     }
     for (int i = 0; i < str_size_2; i++) {
         for (int j = 0; j < nb_str_2; j++) {
             ir_matrice[i][j] = (j + 1);
         }
     }
    return ir_matrice;
}

// Filling of Ic table
int** ic_matrice_func(int str_size_2, int nb_str_2){
    int** ic_matrice;
    ic_matrice = new int* [str_size_2];
    for (int i= 0; i < str_size_2; i++){
        ic_matrice[i] = new int[nb_str_2];
    }
     for (int i = 0; i < str_size_2; i++) {
         for (int j = 0; j < nb_str_2; j++) {
             ic_matrice[i][j] = (i + 1);
        }
     }
     return ic_matrice;
}

// Filling of R table
int**** r_matrice_func(int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string* motif_1, string* motif_2){
     int**** r_matrice;
     r_matrice = new int*** [nb_str_1];
     for (int i= 0; i < nb_str_1; i++){
         r_matrice[i] = new int**[str_size_1];
         for (int j= 0; j < str_size_1; j++){
             r_matrice[i][j] = new int*[nb_str_2];
             for (int k= 0; k < nb_str_2; k++){
                 r_matrice[i][j][k] = new int[str_size_2];
             }
         }
     }


     for (int i = 0; i < nb_str_1; i++) {
         for (int j = 0; j < str_size_1; j++) {
             for (int k = 0; k < nb_str_2; k++) {
               for (int l = 0; l < str_size_2; l++) {
                    r_matrice[i][j][k][l] = classic_seq_edit_distance(motif_1[i].substr(0, j), motif_2[k].substr(0, l));
                 }

             }

         }

     }

     return r_matrice;
}


// Filling of C table
int**** c_matrice_func(int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string* motif_1, string* motif_2) {

    int**** c_matrice;
    c_matrice = new int*** [nb_str_1];
    for (int i= 0; i < nb_str_1; i++){
        c_matrice[i] = new int**[str_size_1];
        for (int j= 0; j < str_size_1; j++){
            c_matrice[i][j] = new int*[nb_str_2];
            for (int k= 0; k < nb_str_2; k++){
                c_matrice[i][j][k] = new int[str_size_2];
            }
        }
    }

    string* motif_trans_1 = new string[str_size_1];
    for (int j = 0; j < str_size_1; j++) {
        string str_trans_1(nb_str_1, '0');
        for (int i = 0; i< nb_str_1 ; i++) {
                 str_trans_1[i]= motif_1[i][j];
             }
        motif_trans_1[j] = str_trans_1;
     }

    string* motif_trans_2 = new string[str_size_2];
    for (int j = 0; j < str_size_2; j++) {
        string str_trans_2(nb_str_2, '0');
        for (int i = 0; i< nb_str_2 ; i++) {
                 str_trans_2[i]= motif_2[i][j];
             }
        motif_trans_2[j] = str_trans_2;
     }

     for (int i = 0; i < nb_str_1; i++) {
         for (int j = 0; j < str_size_1; j++) {
             for (int k = 0; k < nb_str_2; k++) {
                 for (int l = 0; l < str_size_2; l++) {
                     c_matrice[i][j][k][l] = classic_seq_edit_distance(motif_trans_1[i].substr(0, j), motif_trans_2[k].substr(0, l));
                 }

             }

         }

    }

     return c_matrice;
}



//First dimensional transformtion

int f(int i, int j, int k, int nb_str_1, int str_size_1)
{
    //initialisation of e
    int e = i * nb_str_1 + j;
    return k * (str_size_1 * nb_str_1) + e;
}

//Filling of T table
void t_matrice_func(int argc, char *argv[], int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string* motif_1, string* motif_2){

    int rank;
    int nproc;
    // First call MPI_Init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   // Get my rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int taille_motif = motif_1->size();
    cout << rank << " La talle du motif  est de : " << taille_motif << endl;

    int nombre_de_ligne_par_processeur = nb_str_1/ nproc;
    cout << rank << " Mon nombre de ligne est : " << nombre_de_ligne_par_processeur << endl;

    string* new_motif = new string[nombre_de_ligne_par_processeur];
    string** all_new_motifs = new string* [nproc];
    for (int i= 0; i < nproc; i++){
                 all_new_motifs[i] = new string[nombre_de_ligne_par_processeur];
             }
    int cpt = 0;

    while (cpt != taille_motif - 1) {
        int i = 0;
        if (new_motif->size() != nombre_de_ligne_par_processeur){
            (*new_motif).append(motif_1[cpt]);
            cpt += 1;
        }
        else
        {
            all_new_motifs[i] = new_motif;
            string * new_motif = new string[nombre_de_ligne_par_processeur];
            cpt += 1;
            i = i+ 1;
        }
        
    }

    string * my_motif = new string [nombre_de_ligne_par_processeur];
    for (int rk = 0; rk < nproc; rk++){
        if (rank == rk){
            my_motif  = all_new_motifs[rk];
        }
    }

     int** dr_mat = dr_matrice_func(str_size_1, nb_str_1);
     int** dc_mat = dc_matrice_func(str_size_1, nb_str_1);
     int** ir_mat = ir_matrice_func(str_size_2, nb_str_2);
     int** ic_mat = ic_matrice_func(str_size_2, nb_str_2);
     int**** r_mat = r_matrice_func(my_motif[0].size(), str_size_2, my_motif->size(), nb_str_2, my_motif, motif_2);
     int**** c_mat = c_matrice_func(my_motif[0].size(), str_size_2, my_motif->size(), nb_str_2, my_motif, motif_2);
    
    //prepare data partitioning data structure.
    Dag my_dag;
    Cube* my_cube = my_dag.buildDag(nproc, str_size_1, str_size_2, nb_str_1, nb_str_2);

    cout <<'my cube'<< my_cube << endl;

    for (int l = 0; l < nb_str_2; l++)
    {

    }


    int**** t_matrice;
        t_matrice = new int*** [nb_str_1];
            for (int i= 0; i < nb_str_1; i++){
            t_matrice[i] = new int**[str_size_1];
                for (int j= 0; j < str_size_1; j++){
                t_matrice[i][j] = new int*[nb_str_2];
                    for (int k= 0; k < nb_str_2; k++){
                    t_matrice[i][j][k] = new int[str_size_2];
             }
         }
     }






///Initalization maginales of T table
         for (int i = 0; i < my_motif->size(); i++) {
             for (int j = 0; j < my_motif[0].size(); j++) {
                 for (int k = 0; k < nb_str_2; k++) {
                     for (int l = 0; l < str_size_2; l++) {
                         t_matrice[i][j][k][l] = 0;
                         t_matrice[0][j][k][l] = (k+1)*(l+1);
                         t_matrice[i][0][k][l] = (k+1)*(l+1);
                         t_matrice[i][j][0][l] = (i+1)*(j+1);
                         t_matrice[i][j][k][0] = (i+1)*(j+1);
                     }

                 }

             }
         }

         cout <<"nb_str_1->"<<nb_str_1 <<"str_size_1->"<<str_size_1 <<"nb_str_2->"<<nb_str_2 <<"str_size_2->"<<str_size_2 <<endl;

//Filling of T table
         for (int i = 1; i < my_motif->size(); i++) {
             for (int j = 1; j < my_motif[0].size(); j++) {
                 for (int k = 1; k < nb_str_2; k++) {
                     for (int l = 1; l < str_size_2; l++) {
                         int val1 = max(t_matrice[i-1][j][k][l] + dr_mat[i][nb_str_1-1], t_matrice[i][j-1][k][l] + dc_mat[i][nb_str_1-1]);
                         int val2 = max(val1, t_matrice[i][j][k-1][l] + ir_mat[i][nb_str_2-1]);
                         int val3 = max(val2, t_matrice[i][j][k-1][l-1] + ic_mat[i][nb_str_2-1]);
                         int val4 = max(val3, t_matrice[i-1][j][k-1][l] + r_mat[i][j][k][l]);
                         int val5 =  max(val4, t_matrice[i][j-1][k][l-1] + c_mat[i][j][k][l]);
                         int val6 = max(val5, t_matrice[i-1][j-1][k-1][l-1] + c_mat[i-1][j][k-1][l] + r_mat[i][j][k][l]);
                         t_matrice[i][j][k-1][l] = max (val6, t_matrice[i-1][j-1][k-1][l-1] + c_mat[i][j][k][l] + r_mat[i][j-1][k][l-1]);
                     }

                }

             }

         }


//Display of T table
    // cout <<"\n Display T table \n" << endl;
    // for (int i = 0; i < nb_str_1; i++) {
    //      for (int j = 0; j < str_size_1; j++) {
    //          for (int k = 0; k < nb_str_2; k++) {
    //              for (int l = 0; l < str_size_2; l++) {
    //                  cout << i << " " << j << " " << k << " " << l << " -> " << t_matrice[i][j][j][l] << endl;
    //              }
    //          }
    //     }
    //  }

     MPI_Finalize();
   
}



int main(int argc, char *argv[]) {
	time_t first_timer;
    time_t second_timer;
    int start = time( &first_timer);


    unsigned int str_size_1 = 10 , str_size_2 = 10;
    string str1 = "tessa";
    string str2 = "Grace";

    int nb_str_1 = 10, nb_str_2 = 10;
    ifstream datasets("input/datasets.txt");
    string* motif_1= first_pattern_construction(str_size_1, nb_str_1);
    string* motif_2 = second_pattern_construction(str_size_2, nb_str_2);
    t_matrice_func(argc, argv, str_size_1, str_size_2, nb_str_1, nb_str_2, motif_1, motif_2);

 	int end = time(& second_timer);
    	cout << "Treatment start period is about "<< start << endl;
    	cout << "Treatment end period is about "<< end << endl;
     	cout << "Treatment end period is about "<< end-start<< endl;


}
