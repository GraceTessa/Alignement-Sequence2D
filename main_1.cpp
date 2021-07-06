#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>

using namespace std;

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
    return scores[m][n];
}

// Building patterns called "motif"
 string* first_pattern_construction (int str_size_1, int nb_str_1) {
      ifstream datasets("input/datasets.txt");
      string* motif_1;
      motif_1 = new string[nb_str_1];
      if (datasets) {
         string dataset;
         datasets >> dataset;
         ifstream input(dataset);
         if (input) {
             for (int i = 0; i < nb_str_1; i++) {
                 string str1;
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
      motif_2 = new string[nb_str_2];
      if (datasets) {
         string dataset;
         datasets >> dataset;
         ifstream input(dataset);
         if (input) {
             for (int i = 0; i < nb_str_2; i++) {
                string str2;
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



//Filling of T table
void t_matrice_func( int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string* motif_1, string* motif_2){

     time_t first_timer;
     time_t second_timer;

     int** dr_mat = dr_matrice_func(str_size_1, nb_str_1);
     int** dc_mat = dc_matrice_func(str_size_1, nb_str_1);
     int** ir_mat = ir_matrice_func(str_size_2, nb_str_2);
     int** ic_mat = ic_matrice_func(str_size_2, nb_str_2);
     int**** r_mat = r_matrice_func(str_size_1, str_size_2, nb_str_1, nb_str_2, motif_1, motif_2);
     int**** c_mat = c_matrice_func(str_size_1, str_size_2, nb_str_1, nb_str_2, motif_1, motif_2);

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


     int start = time(nullptr)*1000;


//Initalization maginales of T table
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

         cout <<"nb_str_1->"<<nb_str_1 <<"str_size_1->"<<str_size_1 <<"nb_str_2->"<<nb_str_2 <<"str_size_2->"<<str_size_2 <<endl;

//Filling of T table
         for (int i = 1; i < nb_str_1; i++) {
             for (int j = 1; j < str_size_1; j++) {
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


        int end = time(nullptr)*1000;

//Display of T table
    cout <<"\n Display T table \n" << endl;
    for (int i = 0; i < nb_str_1; i++) {
         for (int j = 0; j < str_size_1; j++) {
             for (int k = 0; k < nb_str_2; k++) {
                 for (int l = 0; l < str_size_2; l++) {
                     cout << i << " " << j << " " << k << " " << l << " -> " << t_matrice[i][j][j][l] << endl;
                 }
             }
        }
     }
    cout << "Treatment start period is about "<< start << endl;
    cout << "Treatment end period is about "<< end << endl;
}



int main(int argc, char *argv[]) {

    unsigned int str_size_1 = 256, str_size_2 = 256;
    string str1 = "tessa";
    string str2 = "Grace";

    int nb_str_1 = 256, nb_str_2 = 256;
    ifstream datasets("input/datasets.txt");
    string* motif_1= first_pattern_construction(str_size_1, nb_str_1);
    string* motif_2 = second_pattern_construction(str_size_2, nb_str_2);
    t_matrice_func(str_size_1, str_size_2, nb_str_1, nb_str_2, motif_1, motif_2);




}
