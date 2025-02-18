#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream> // For file handling
#include <iomanip> // For formatting
#include <complex>
#include <typeinfo>

using namespace std;
using vec=vector <double>;
using comp=complex <double>;

/*
This first part declares lapack functions for matrix diagonalization
*/

// zgeev_ is a symbol in the LAPACK library files
//Diagonalizes complex matrices either symmetric or nonsymmetric
extern "C"{
    extern int zgeev_(char*,char*,int*,comp*,int*,comp*,comp*,int*,
        comp*,int*,comp*,int*,double*,int*);
    }

// dsyev_ is a symbol in the LAPACK library files
//Diagonalizes real symmetric matrices 
extern "C" {  
    extern int dsyev_(char*,char*,int*, double*,
        int*,double*,double*,int*,int* 
        );
    }	
// dgeev_ is a symbol in the LAPACK library files
//Diagonalizes real nonsymmetric matrices  
  extern "C" {
    extern int dgeev_(char*,char*,int*,double*,int*,double*,
         double*, double*, int*, double*, int*, double*, int*, int*);
    }



class diagonalization{
    vector<vec> m_b;
    vector<vector <comp>> m_bc;
    int ind; //flag for determining the nature of the matrix to be diagonalized 1 for real and  2 for complex
    struct comple{
        vector<vector<comp>> eigenvectors;
        vector<comp> eigenvalues;
    }; 
    struct reall{
        vector<vec> eigenvectors;
        vec eigenvalues;
    };

    public:
    bool eigv;
    int inds; //flag for determing if the real matrix is symmetric (1) or nonsymmetric (2)
    //variables where is going to be storage the eigenvalues and eigenvectors
    comple data_complex;
    reall data_real_symm;
    comple data_real_nonsym;

    diagonalization(
        bool eigvec, //this define if eigenvector matrix is return or no
        vector<vec> a, //Matrix to be diagonalized
        int iindi=0  // flag for matrix symmetry 

    )
    : m_b(a),ind(1)
    {
        eigv=eigvec;
        if (iindi==0){
            cout<<"ERROR, you must indicate the symmetry of the marix. 1 for symmetric and 2 for nonsymmetric"<<endl;
            exit(-1);
        }
        if ((iindi==1) &&(a[0][1]!=a[1][0])){
            cout<<"ERROR, The matrix provided is not symmetric."<<endl;
            cout<<"Use the flag 2 for nonsymmetric matrices."<<endl;
            exit(-1);
        }
        if((iindi==2)&&(a[0][1]==a[1][0])){
            cout<<"ERROR, The matrix provided is symmetric."<<endl;
            cout<<"Use the flag 1 for symmetric matrices."<<endl;
            exit(-1);
        }
        else{inds=iindi;}        
    }

    diagonalization(
        bool eigvec, //this define if eigenvector matrix is return or no
        vector<vector <comp>> a //Matrix to be diagonalized
    ): m_bc(a),ind(2)
    {eigv=eigvec;}



    void diag(){
        if(ind==1){
            if(inds==1){
                diag_real_sym(m_b);
            }
            else{
                diag_real_nonsym(m_b);
            }
        }
        else{ 
            diag_compl(m_bc);            
    }
    }

    private:

    //functions to prepare the matrix in the lapack's format 
    //complex matrices
    comp* transf_c(vector <vector <comp>> a,int n){
        int l;
        l=0;
        comp* b=new comp [n*n];
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                b[l]=a[j][i];
                l++;
            }
        }
        return b;
    }
    //real matrices
    double* transf_r(vector <vec> a,int n){
        int l;
        l=0;
        double* b=new double [n*n];
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                b[l]=a[j][i];
                l++;
            }
        }
        return b;
    }

    //functions to transform eigenvector lapack's format into a readable matrix
    //complex matrices
    vector<vector<comp>> eigevect_c(comp* a,int n){
        vector<vector<comp>> b(n,vector<comp> (n,0.0));
        int i,j,l;
        l=0;
        for(i=0;i<n;i++){
            for (j=0;j<n;j++){
                b[j][i]=a[l];
                l++;
            }
        }
        return b;   
    }

    //symmetric real matrices
    vector <vec> eigevect_sr(double* a,int n){
        int i,j,l;
        vector<vec> b(n, vec (n, 0.0));
        l=0;
        for (i=0;i<n;i++){
          for (j=0;j<n;j++){
            b[j][i]=a[l];
            l++;
          }
        }   
        return b;
    }
    //nonsymmetric real matrices
    vector <vector <comp>> eigevect_ar(double* eig,double* wr,double* wi,int n){
        int i,j,l,l2;
        double temp;
        vector<vector<double>> b(2 * n, vector<double>(n, 0.0));
        vector<vector <comp>>c(n, vector<comp>(n, 0.0));
        //determining the number of Real and complex eigenvectors and storage them in the matrix of eigenvalues
        l=0;
        l2=0;
        for (i=0;i<n;i++){
            for (j=0;j<2*n;j++){
                if(wi[i]!=0){//complex eigenvectors
                    if(fabs(wi[i])!=fabs(wi[i-1])){
                        b[j][i]=eig[l];
                        l2=0;
                        l2=l-2*(n)+1;
                        l++;                    
                    }
                    else{
                        if(j<n){
                            b[j][i]=eig[l2];
                            l2++;
                        }
                        else{
                            b[j][i]=eig[l2]*(wi[i]/fabs(wi[i]));
                            l2++;
                        }
                    }
                }
                else{//real eigenvectors
                    if(j<n){
                        b[j][i]=eig[l];
                        l++;
                    }
                    else{
                        b[j][i]=0.0;
                    }
                }

            }
        }
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){ 
                c[j][i]={b[j][i],b[j+n][i]};
            }
        }
        return c;
    }


    void diag_compl(vector <vector <comp>> a){
        int n;
        char JOBVR1;
        n=a.size();
        comp* b;
        b=transf_c(a,n);
        //Definining variables for diagonalization
        if (eigv){
            JOBVR1='V';
        }
        else{
            JOBVR1='N';
        }
        char 	JOBVL='N';
        char JOBVR=JOBVR1;
        comp* W=new comp [n];
        int LDVL=1;
        comp* VL;
        int LDVR=n;
        comp* VR=new comp [n*n];
        int LWORK=6*n;
        comp* 	WORK=new comp [LWORK];
        double* RWORK=new double [2*n];
        int INFO;
        zgeev_(&JOBVL,&JOBVR,&n,b,&n,W,VL,&LDVL,
            VR,&LDVR,WORK,&LWORK,RWORK,&INFO);
        if(INFO!=0){
            cout<<"Error: zgeev returned error code " << INFO << endl;
            exit(-1);
            }
        vector<comp> val;
        vector<vector<comp>>u=eigevect_c(VR,n);
        for (int i=0;i<n;i++){
            val.push_back(W[i]);
        }
        data_complex.eigenvectors=u;
        data_complex.eigenvalues=val;

    }

    void diag_real_sym(vector<vec> a){
        int n;
        char jobzz;
        n=a.size();
        double* b;
        b=transf_r(a,n);
        //Definining variables for diagonalization
        if(eigv){
            jobzz='V';
        }
        else{
            jobzz='N';
        }
        char	JOBZ=jobzz;
        char 	UPLO='U';
        double* W=new double [n];
        int	LWORK=6*n;
        double *WORK=new double [LWORK*6];
        int INFO;
        dsyev_(&JOBZ,&UPLO,&n,b,&n,W,WORK,&LWORK,&INFO);
        if(INFO!=0){
            cout<<"Error: dsyev returned error code " << INFO << endl;
            exit(-1);
            } 
        vector<vec> u(n, vec(n,0.0));
        vec v;

        if(eigv){
            u=eigevect_sr(b,n);
        } 
        for (int i=0;i<n;i++){
            v.push_back(W[i]);
        }
        data_real_symm.eigenvalues=v;
        data_real_symm.eigenvectors=u;

    }

    void diag_real_nonsym(vector <vec> a){
    int n;
    char jobzz;
    n=a.size();
    double* b=transf_r(a,n);
    //Definining variables for diagonalization
    if(eigv){
        jobzz='V';
    }
    else{
        jobzz='N';
    }
    char	JOBVL='N';
    char 	JOBVR=jobzz;
    double* WR=new double [n];
    double*	WI=new double [n];
    int LDVL=1;
    double* VL;
    int LDVR=n;
    double*	VR=new double [n*n];
    int LWORK=6*n;
    double*	WORK=new double[LWORK];
    int INFO;
    dgeev_(&JOBVL,&JOBVR,&n,b,&n,WR,WI,VL,&LDVL,VR,&LDVR,WORK,&LWORK,&INFO);
    if(INFO!=0){
        cout<<"Error: dgeev returned error code " << INFO << endl;
        exit(-1);
        }
    vector<comp> v; 
    for(int i=0;i<n;i++){
        v.push_back({WR[i],WI[i]});

    }
    vector <vector <comp>> u=eigevect_ar(VR,WR,WI,n);

    data_real_nonsym.eigenvalues=v;
    data_real_nonsym.eigenvectors=u;
    
    }
};
