#include <iostream>
#include "Matrix_diag.cpp"
using namespace std;


int main(){
    //example 1:  Complex Matrix
    int n;
    n=3;
    vector <vector <comp>> a(n,vector <comp>(n,0.0));
    a={{-1.0+3.0i, -8.0+2.0i,  0.0+1.0i},
    {-1.0-2.0i,  1.0-3.0i, -5.0-5.0i},
    {3.0-9.0i,  0.0-7.0i,  2.0-6.0i}};
    diagonalization di(true,a);
    di.diag();
    vector <comp> v=di.data_complex.eigenvalues;
    cout << "--- Eigenvalues---" << endl;
    int i=1;
    for (comp t:v){
        cout<<"Eigenvalue\t"<<i<<":\t"<<t<<endl;
        i++;
    }
    vector <vector <comp>> u=di.data_complex.eigenvectors;
    cout << "--- Eigenvectors---" << endl;
    for (int i=0;i<n;i++){
        cout<<"eigenvector"<<'\t'<<i+1<<':'<<'\t';
        for (int j=0;j<n;j++){
            cout<<u[j][i]<<'\t';
        }
        cout<<endl;
    }

    //example 2:  Real Symmetric Matrix
    int n2;
    n2=2;
    vector <vector <double>> a2;
    a2={{1.0,2.0},{2.0,4.0}};
    diagonalization di_re(true,a2,1);
    di_re.diag();
    vec eigv_r=di_re.data_real_symm.eigenvalues;
    vector <vec> eigvec_r=di_re.data_real_symm.eigenvectors;

    cout << "--- Eigenvalues---" << endl;
    i=1;
    for (double ev:eigv_r){
        cout<<"Eigenvalue\t"<<i<<":\t"<<ev<<endl;
        i++;
    }
    cout << "--- Eigenvectors---" << endl;
    for (int i=0;i<n2;i++){
        cout<<"eigenvector"<<'\t'<<i+1<<':'<<'\t';
        for (int j=0;j<n2;j++){
            cout<<eigvec_r[j][i]<<'\t';
        }
        cout<<endl;
    }

    //example 3:  Real NonSymmetric Matrix
    int n3;
    n3=5;
    vector <vector <double>> a3;
    a3={{-1.01, 3.98, 3.30, 4.43, 7.31},
    { 0.86, 0.53, 8.26, 4.96, -6.43},
    {-4.60, -7.04, -3.89, -7.66, -6.16},
    { 3.31, 5.29, 8.20, -7.33, 2.47},
    {-4.81, 3.55, -1.51, 6.18, 5.58}};

    diagonalization di_nons(true,a3,2);
    di_nons.diag();
    vector<comp> eigv_nonr=di_nons.data_real_nonsym.eigenvalues;
    vector<vector<comp>> eigvec_nonsr=di_nons.data_real_nonsym.eigenvectors;


    cout << "--- Eigenvalues---" << endl;
    i=1;
    for (comp ev:eigv_nonr){
        cout<<"Eigenvalue\t"<<i<<":\t"<<ev<<endl;
        i++;
    }
    cout << "--- Eigenvectors---" << endl;
    for (int i=0;i<n3;i++){
        cout<<"eigenvector"<<'\t'<<i+1<<':'<<'\t';
        for (int j=0;j<n3;j++){
            cout<<eigvec_nonsr[j][i]<<'\t';
        }
        cout<<endl;
    }


    return 0;


}
