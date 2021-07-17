
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;

// Classique sequences edit distance method
int classic_seq_edit_distance(const string &a, const string &b)
{
    int m = a.size();
    int n = b.size();

    vector<vector<int>> scores(m + 1, vector<int>(n + 1));

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
vector<vector<string>> pattern_construction (unsigned int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string str1, string str2, vector<string> motif_1, vector<string> motif_2) {
     ifstream datasets("input/datasets.txt");
     if (datasets) {
        string dataset;
        datasets >> dataset;
        ifstream input(dataset);
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


    vector<vector<string>> motifs;
    motifs.push_back(motif_1);
    motifs.push_back(motif_2);

    return motifs;
   
}

// Filling of Dr table
void dr_matrice_func(vector<vector<int>> & dr_matrice){
    for (int i = 0; i < dr_matrice.size(); i++) {
        for (int j = 0; j < dr_matrice.size(); j++) {
            dr_matrice[i][j] = (j + 1);
        }
    }
}


// Filling of Dc table
void dc_matrice_func(vector<vector<int>> & dc_matrice){
    for (int i = 0; i < dc_matrice.size(); i++) {
        for (int j = 0; j < dc_matrice.size(); j++) {
            dc_matrice[i][j] = (i + 1);
        }
    }
}

// Filling of Ir table
void ir_matrice_func(vector<vector<int>> & ir_matrice){
    for (int i = 0; i < ir_matrice.size(); i++) {
        for (int j = 0; j < ir_matrice.size(); j++) {
            ir_matrice[i][j] = (j + 1);
        }
    }
}

// Filling of Ic table 
void ic_matrice_func(vector<vector<int>> & ic_matrice){
    for (int i = 0; i < ic_matrice.size(); i++) {
        for (int j = 0; j < ic_matrice.size(); j++) {
            ic_matrice[i][j] = (i + 1);
        }
    }
}

// Filling of R table
void r_matrice_func(vector<vector<vector<vector<int>>>> &r_matrice,unsigned int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string str1, string str2, vector<string> motif_1, vector<string> motif_2){

    vector<vector<string>> motifs = pattern_construction (str_size_1, str_size_2, nb_str_1, nb_str_2, str1, str2, motif_1, motif_2);

    for (int i = 0; i < motifs.size(); i++) { 
        for (int j = 0; j < motifs[i].size(); j++) {
            cout << motifs[i][j] << endl; 
            }
    } 

    motif_1 = motifs[0];
    motif_2 = motifs[1];

    for (int i = 0; i < nb_str_1; i++) {
        for (int j = 0; j < str_size_1; j++) {
            for (int k = 0; k < nb_str_2;k++) {
                for (int l = 0; l < str_size_2; l++) {
                   r_matrice[i][j][k][l] = classic_seq_edit_distance(motif_1.at(i).substr(0, j), motif_2.at(k).substr(0, l));
                }
              
            }
            
        }

    }
}


// Filling of C table
void c_matrice_func(vector<vector<vector<vector<int>>>> &c_matrice, unsigned int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string str1, string str2, vector<string> motif_1, vector<string> motif_2){

    vector<vector<string>> motifs = pattern_construction ( str_size_1, str_size_2, nb_str_1, nb_str_2, str1, str2, motif_1, motif_2);

    for (int i = 0; i < motifs.size(); i++) { 
        for (int j = 0; j < motifs[i].size(); j++) {
            cout << motifs[i][j] << endl; 
            }
    } 

    motif_1 = motifs[0];
    motif_2 = motifs[1];


    int  m = motif_1.size();
    int  n = m >0 ? motif_1[0].size() : 0;
    vector<string> motif_trans_1(n, string(m, '\0'));

    for (int j = 0; j< n ; j++){
            for (int i = 0; i< m; i++){
                motif_trans_1[j][i] = motif_1[i][j];
            }  
    }


    int  m2 = motif_2.size();
    int  n2 = m2 >0 ? motif_2[0].size() : 0;
    vector<string> motif_trans_2(n2, string(m2, '\0'));

    for (int j = 0; j< n2 ; j++){
            for (int i = 0; i< m2; i++){
                motif_trans_2[j][i] = motif_2[i][j];
            }  
    }
        
    for (int i = 0; i < nb_str_1; i++) {
        for (int j = 0; j < str_size_1; j++) {
            for (int k = 0; k < nb_str_2; k++) {
                for (int l = 0; l < str_size_2; l++) {
                    c_matrice[i][j][k][l] = classic_seq_edit_distance(motif_trans_1.at(i).substr(0, j), motif_trans_2.at(k).substr(0, l));
                }
              
            }
            
        }

    }
}

// Filling of T table
void t_matrice_func(vector<vector<vector<vector<int>>>> & t_matrice, ifstream & datasets, unsigned int str_size_1, int str_size_2, int nb_str_1, int nb_str_2, string str1, string str2, vector<string> motif_1, vector<string> motif_2){
    time_t first_timer;
    time_t second_timer;
    vector<vector<int>> dr_mat(nb_str_1, vector<int>(str_size_1, 0));
    vector<vector<int>> dc_mat(nb_str_1, vector<int>(str_size_1, 0));
    vector<vector<int>> ir_mat(nb_str_2, vector<int>(str_size_2, 0));
    vector<vector<int>> ic_mat(nb_str_2, vector<int>(str_size_2, 0));
    vector<vector<vector<vector<int>>>> r_mat(nb_str_1, vector<vector<vector<int>>>(str_size_1, vector<vector<int>>(nb_str_2, vector<int>(str_size_2, 0))));
    vector<vector<vector<vector<int>>>> c_mat(nb_str_1, vector<vector<vector<int>>>(str_size_1, vector<vector<int>>(nb_str_2, vector<int>(str_size_2, 0))));
    
    int start = time( &first_timer);
    dr_matrice_func(dr_mat);
    dc_matrice_func(dc_mat);
    ir_matrice_func(ir_mat);
    ic_matrice_func(ic_mat);
    r_matrice_func(r_mat, str_size_1, str_size_2, nb_str_1, nb_str_2, str1, str2, motif_1, motif_1);
    c_matrice_func(c_mat, str_size_1, str_size_2, nb_str_1, nb_str_2, str1, str2, motif_1, motif_1);
    

    // Initalization maginales of T table
        for (int i = 0; i < t_matrice.size(); i++) {
            for (int j = 0; j < t_matrice.size(); j++) {
                for (int k = 0; k < t_matrice.size(); k++) {
                    for (int l = 0; l < t_matrice.size(); l++) {
                        t_matrice[0][j][k][l] = (k+1)*(l+1);
                        t_matrice[i][0][k][l] = (k+1)*(l+1);
                        t_matrice[i][j][0][l] = (i+1)*(j+1);
                        t_matrice[i][j][k][0] = (i+1)*(j+1);
                    }
                
                }
                
            }

        }

        // Filling of T table
        for (int i = 1; i < t_matrice.size(); i++) {
            for (int j = 1; j < t_matrice.size(); j++) {
                for (int k = 1; k < t_matrice.size(); k++) {
                    for (int l = 1; l < t_matrice.size(); l++) {
                        int val1 = max(t_matrice[i-1][j][k][l] + dr_mat[i][nb_str_1-1], t_matrice[i][j-1][k][l] + dc_mat[i][nb_str_1-1]);
                        int val2 = max(val1, t_matrice[i][j][k-1][l] + ir_mat[i][nb_str_2-1]);
                        int val3 = max(val2, t_matrice[i][j][k-1][l-1] + ic_mat[i][nb_str_2-1]);
                        int val4 = max (val3, t_matrice[i-1][j][k-1][l] + r_mat[i][j][k][l]);
                        int val5 =  max (val4, t_matrice[i][j-1][k][l-1] + c_mat[i][j][k][l]);
                        int val6 = max (val5, t_matrice[i-1][j-1][k-1][l-1] + c_mat[i-1][j][k-1][l] + r_mat[i][j][k][l]);
                        t_matrice[i][j][k-1][l] = max (val5, t_matrice[i-1][j-1][k-1][l-1] + c_mat[i][j][k][l] + r_mat[i][j-1][k][l-1]);
                        
                    }
                
                }
                
            }

        }


        int end = time(& second_timer);
        
        // Display of T table
        cout <<"\n Display T table \n" << endl;
        for (int i = 0; i < t_matrice.size(); i++) {
            for (int j = 0; j < t_matrice.size(); j++) {
                for (int k = 0; k < t_matrice.size(); k++) {
                    for (int l = 0; l < t_matrice.size(); l++) {
                        cout << i << " " << j << " " << k << " " << l << " -> " << t_matrice[i][j][j][l] << endl;
                    }
                }
            }
        }
        cout << "Treatment start period is about "<< start << endl;
        cout << "Treatment end period is about "<< end << endl;
}



int main(int argc, char *argv[]) {

    unsigned int str_size_1 = 30, str_size_2 = 30;
    string str1 = "tessa";
    string str2 = "Grace";
    
    vector<string> motif_1;
    vector<string> motif_2;

    int nb_str_1 = 30, nb_str_2 = 30;
    ifstream datasets("input/datasets.txt");
    //pattern_construction(datasets, str_size_1, str_size_2, nb_str_1, nb_str_2, str1, str2, motif_1, motif_2);
    vector<vector<vector<vector<int>>>> t_matrice(nb_str_1, vector<vector<vector<int>>>(str_size_1, vector<vector<int>>(nb_str_2, vector<int>(str_size_2, 0))));
    t_matrice_func(t_matrice, datasets, str_size_1, str_size_2, nb_str_1, nb_str_2, str1, str2, motif_1, motif_2);
   
}
