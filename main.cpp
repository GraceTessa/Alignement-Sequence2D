#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "mpi.h"



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

    int** scores;
    scores = new int* [m + 1];
    for (int i= 0; i < m+1; i++){
        scores[i] = new int[n+1];
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

            scores[i][j] = min({x, y, z});
        }
    }
    int result =  scores[m][n];

     // cout << "TESSA"<< endl;

    for (int i= 0; i < m+1; i++){
        delete scores[i];
    }
    delete[] scores;

    return result;
}

// Building patterns called "motif"
 vector<string> first_pattern_construction (int str_size_1, int nb_str_1) {
      ifstream datasets("input/datasets.txt");
      vector<string> motif_1;
      string str1;
      // motif_1 = new string[nb_str_1];
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
                 motif_1.push_back(str1);
             }

             input.close();
         }
         datasets.close();
     }

     return motif_1;

}

 vector<string> second_pattern_construction (int str_size_2, int nb_str_2) {
      ifstream datasets("input/datasets.txt");
      vector<string> motif_2;
      string str2;
      // motif_2 = new string[nb_str_2];
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
                 motif_2.push_back(str2);
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
int**** r_matrice_func(int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, vector<string> motif_1, vector<string> motif_2){
     int**** r_matrice;
     r_matrice = new int*** [nb_str_1];
     for (int i= 0; i < nb_str_1; i++){
         r_matrice[i] = new int**[str_size_1];
         for (int j= 0; j < str_size_1; j++){
             r_matrice[i][j] = new int*[nb_str_2];
             for (int k= 0; k < nb_str_2; k++){
                 r_matrice[i][j][k] = new int[str_size_2];
                 // cout << "TESSA"<< endl;
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
int**** c_matrice_func(int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, vector<string> motif_1, vector<string> motif_2) {

    int**** c_matrice;
    c_matrice = new int*** [nb_str_1];
    for (int i= 0; i < nb_str_1; i++){
        c_matrice[i] = new int**[str_size_1];
        for (int j= 0; j < str_size_1; j++){
            c_matrice[i][j] = new int*[nb_str_2];
            for (int k= 0; k < nb_str_2; k++){
                c_matrice[i][j][k] = new int[str_size_2];
                // cout << "TESSA"<< endl;
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
                     // cout << "TESSA"<< endl;
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
void t_matrice_func(int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, vector<string> motif_1, vector<string> motif_2){
    MPI::Init();
    time_t first_timer;
    time_t second_timer;
    double start = time( &first_timer);

    int nproc = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();
    int etiquette = 10;
    

    time_t first_pre_timer;
    time_t second_pre_timer;
    time_t f_timer;
    time_t s_timer;

    double pre_start = time(&first_pre_timer);
    int taille_motif = motif_1.size();
    cout << rank << " motif length = "<< taille_motif << endl;

    int nombre_de_ligne_par_processeur = taille_motif / nproc;

    // Partition list of motifs for all processors
    vector<vector<string>> all_new_motifs;
    
    /*Begin motif partitionnig */
    int step = 0;
    for (int i = 0; i < nproc; i++)
    {
        vector<string> new_motif;
        while (new_motif.size() != nombre_de_ligne_par_processeur)
        {
            new_motif.push_back(motif_1[step]);
            step++;
        }
        
        // cout << rank << " new motif length = "<< new_motif.size() << endl;
        all_new_motifs.push_back(new_motif);
        
    }
    cout << "i completed the partitionning function with success !!!" << endl;
    cout << "all motifs length = "<< all_new_motifs.size() << endl;

    /* Each processor takes it own motif*/
    vector<string>  my_motif;
    for (int rk = 0; rk < nproc; rk++){
        if (rank == rk){
            my_motif  = all_new_motifs[rk];
        }
    }
    cout << rank << " i took my own motif" << endl;
    /* End motif partitionning*/

    MPI::COMM_WORLD.Barrier();
    // return exit(0);

    int** dr_mat = dr_matrice_func(str_size_1, nb_str_1);
    int** dc_mat = dc_matrice_func(str_size_1, nb_str_1);
    int** ir_mat = ir_matrice_func(str_size_2, nb_str_2);
    int** ic_mat = ic_matrice_func(str_size_2, nb_str_2);
    int**** r_mat = r_matrice_func(my_motif[0].size(), str_size_2, my_motif.size(), nb_str_2, my_motif, motif_2);
    cout << rank << " finish r matrix computations " << endl;

    int**** c_mat = c_matrice_func(my_motif[0].size(), str_size_2, my_motif.size(), nb_str_2, my_motif, motif_2);
    cout << rank << " finish c matrix computations " << endl;
    
    double pre_end = time(&second_pre_timer);

    cout << rank << " finish pre-computations " << endl;
    

    // MPI::Datatype mpi_vector;
    // int dest, source = 0;

    int my_nb_of_lines = my_motif.size();
    int my_nb_of_columns = my_motif[0].size();
    int matrix_of_matrix_length = my_nb_of_lines*my_nb_of_columns*nb_str_2*str_size_2;

    
    int r_mat_prime[my_nb_of_lines][my_nb_of_columns][nb_str_2][str_size_2];
    int c_mat_prime[my_nb_of_lines][my_nb_of_columns][nb_str_2][str_size_2];

    cout << rank << " i am ready to copy and line = " << my_nb_of_lines << " with column = " << my_nb_of_columns  << endl;

    for (int i = 0; i < my_nb_of_lines; i++)
    {
        for (int j = 0; j < my_nb_of_columns; j++)
        {
            for (int k = 0; k < nb_str_2; k++)
            {
                for (int l = 0; l < str_size_2; l++)
                {
                    // cout << r_mat[i][j][k][l] << endl;
                    r_mat_prime[i][j][k][l] = r_mat[i][j][k][l];
                    c_mat_prime[i][j][k][l] = c_mat[i][j][k][l];
                }
                
            }
            
        }
        
    }
    
    cout << rank << " i have finish to copy and line = " << my_nb_of_lines << " with column = " << my_nb_of_columns  << endl;
    

    MPI::COMM_WORLD.Barrier();
    double end = time( &first_timer);
    cout << rank << " pre-computation time = " << end-start << endl;

    time_t first_comm_timer;
    time_t second_comm_timer;
    
    /* Start Communication*/
    if (rank != 0){
        /*count = ;
        blocklengths = bd.core_data.d.column;
        stride = 1;
        MPI_Type_vector(count, blocklengths, stride, MPI_INT, &mpi_vector);
        MPI_Type_commit(&mpi_vector);
        */
       
        double begin_comm = time( &first_comm_timer);
        cout << rank << " Entering to sending matrix" << endl;
        //MPI_Send(&r_mat, matrix_of_matrix_length, MPI_INT, dest, etiquette, MPI_COMM_WORLD);
        MPI::COMM_WORLD.Send(&r_mat_prime, matrix_of_matrix_length,  MPI_INT, 0 , etiquette);
        cout << rank << " sending r matrix to rank = 0" << endl;
        
        MPI::COMM_WORLD.Send(&c_mat_prime, matrix_of_matrix_length,  MPI_INT, 0 , etiquette);
        cout << rank << " i am going to send data to rank = 0" << endl;
        // return exit(0);
        double end_comm = time( &second_comm_timer);
        cout << rank << " communication time = " << end_comm-begin_comm << endl;
    }
    else
    {
        int offset = 0;
        int all_r_matrice[nb_str_1][str_size_1][nb_str_2][str_size_2];
        int all_c_matrice[nb_str_1][str_size_1][nb_str_2][str_size_2];
        
        for (int i = 0; i < my_nb_of_lines; i++)
        {
            for (int j = 0; j < my_nb_of_columns; j++)
            {
                for (int k = 0; k < nb_str_2; k++)
                {
                    for (int l = 0; l < str_size_2; l++)
                    {
                        // cout << r_matrice[i][j][k][l] << endl;
                        all_r_matrice[offset][j][k][l] = r_mat_prime[i][j][k][l];
                        all_c_matrice[offset][j][k][l] = c_mat_prime[i][j][k][l];
                    }   
                }
            }
            offset += 1;
        }
        cout << rank << " offset = " << offset << endl;

        for (int rk = 1; rk < nproc; rk++)
        {
            int r_matrice[my_nb_of_lines][my_nb_of_columns][nb_str_2][str_size_2];
            int c_matrice[my_nb_of_lines][my_nb_of_columns][nb_str_2][str_size_2];

            int source = rk;
            cout << rank << " Entering to receiving matrix" << endl;
            // MPI_Recv(&r_matrice, matrix_of_matrix_length, MPI_INT, source, etiquette, MPI_COMM_WORLD, &status);
            MPI::COMM_WORLD.Recv(&r_matrice, matrix_of_matrix_length,  MPI_INT, source, etiquette);
            MPI::COMM_WORLD.Recv(&c_matrice, matrix_of_matrix_length,  MPI_INT, source, etiquette);
            
            for (int i = 0; i < my_nb_of_lines; i++)
            {
                for (int j = 0; j < my_nb_of_columns; j++)
                {
                    for (int k = 0; k < nb_str_2; k++)
                    {
                        for (int l = 0; l < str_size_2; l++)
                        {
                            // cout << r_matrice[i][j][k][l] << endl;
                            all_r_matrice[offset][j][k][l] = r_matrice[i][j][k][l];
                            all_c_matrice[offset][j][k][l] = c_matrice[i][j][k][l];
                        }   
                    }
                }
                offset += 1;
            }
            cout << rank << " offset = " << offset << endl;
            /*for (int i = 0; i < nb_str_1; i++)
            {
                for (int j = 0; j < str_size_1; j++)
                {
                    for (int k = 0; k < nb_str_2; k++)
                    {
                        for (int l = 0; l < str_size_2; l++)
                        {
                            cout << all_c_matrice[i][j][k][l] << endl;
                        }
                        
                    }
                    
                }
        
            }*/
            cout << rank << " i am going to receive data from all the other ranks" << endl;
            
            
            //return exit(0);
        }

    
    /* End Communication*/
       
    //    MPI_Bcast(&r_matrice, pre_computation_final_size, MPI_INT, 0, MPI_COMM_WORLD);
   

        //prepare data partitioning data structure.
        Dag my_dag;
        Cube* my_cube = my_dag.buildDag(nproc, str_size_1, str_size_2, nb_str_1, nb_str_2);

        cout <<"my cube"<< my_cube << endl;

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
        for (int i = 0; i < nb_str_1; i++) {
            for (int j = 0; j < str_size_1; j++) {
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

        

        //Filling of T table
        for (int i = 1; i < nb_str_1; i++) {
            for (int j = 1; j < str_size_1; j++) {
                for (int k = 1; k < nb_str_2; k++) {
                    for (int l = 1; l < str_size_2; l++) {
                        int val1 = max(t_matrice[i-1][j][k][l] + dr_mat[i][nb_str_1-1], t_matrice[i][j-1][k][l] + dc_mat[i][nb_str_1-1]);
                        int val2 = max(val1, t_matrice[i][j][k-1][l] + ir_mat[i][nb_str_2-1]);
                        int val3 = max(val2, t_matrice[i][j][k-1][l-1] + ic_mat[i][nb_str_2-1]);
                        int val4 = max(val3, t_matrice[i-1][j][k-1][l] + all_r_matrice[i][j][k][l]);
                        int val5 =  max(val4, t_matrice[i][j-1][k][l-1] + all_c_matrice[i][j][k][l]);
                        int val6 = max(val5, t_matrice[i-1][j-1][k-1][l-1] + all_c_matrice[i-1][j][k-1][l] + all_r_matrice[i][j][k][l]);
                        t_matrice[i][j][k-1][l] = max (val6, t_matrice[i-1][j-1][k-1][l-1] + all_c_matrice[i][j][k][l] + all_r_matrice[i][j-1][k][l-1]);
                    }
                }
            }
        }


        //Display of T table
        /*cout <<"\n Display T table \n" << endl;
        for (int i = 0; i < nb_str_1; i++) {
            for (int j = 0; j < str_size_1; j++) {
                for (int k = 0; k < nb_str_2; k++) {
                    for (int l = 0; l < str_size_2; l++) {
                        cout << i << " " << j << " " << k << " " << l << " -> " << t_matrice[i][j][j][l] << endl;
                    }
                }
            }
        }*/
    }     
    
    
    
    MPI::Finalize();
   
}



int main(int argc, char *argv[]) {
	time_t first_timer;
    time_t second_timer;
    double start = time( &first_timer);
    cout << " i entered the main function " << endl;

    unsigned int str_size_1 = 50 , str_size_2 = 50;
    
    int nb_str_1 = 50, nb_str_2 = 50;
    ifstream datasets("input/datasets.txt");
    vector<string> motif_1= first_pattern_construction(str_size_1, nb_str_1);
    // cout << " i got the first motif " << endl;

    vector<string> motif_2 = second_pattern_construction(str_size_2, nb_str_2);
    // cout << " i got the second motif " << endl;

    t_matrice_func(str_size_1, str_size_2, nb_str_1, nb_str_2, motif_1, motif_2);

 	double end = time(& second_timer);
    cout << "Treatment start period is about "<< start << endl;
    cout << "Treatment end period is about "<< end << endl;
    cout << "Treatment end period is about "<< end-start<< endl;


}
