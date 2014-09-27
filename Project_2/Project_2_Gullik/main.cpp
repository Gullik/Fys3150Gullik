#include <iostream>
#include <fstream>
#include "armadillo"
#include <math.h>
#include <time.h>
#include <Proj2Lib.h>

using namespace std;
using namespace arma;

// ------------------------------------------------
// On the last assignment I sent in the indents dissappeared
// send me a mail if it happens again so I can fix it so you don't need to read the code unindented
// gullik.killie@gmail.com
// ------------------------------------------------


int main()
{
    int N = 200;                      //Number of steps
    int Iterations = 100000;
    double epsilon = 0.00001;
    cout << "Enter wanted number of steps? "<< endl;
//    cin >> n;
    cout << "Enter wanted number of iterations" << endl;
//    cin >> Iterations;

    double rho_min = 0;
    double rho_max = 10;

// Creating the vectors to populate the matrix A
    double h = (rho_max - rho_min)/N;

    vec d = zeros(N);
    vec e = zeros(N);

    for(int i = 0; i < N - 2; i++)
    {
//        d[i] = 2; //2/pow(h,2) + i*h;
//        e[i] = 1; //-1/pow(h,2);
        d[i] = 2/pow(h,2) + pow((i+1)*h,2);
        e[i] = -1/pow(h,2);
    }

cout << "The stepsize is " << h << endl;

    mat A = zeros<mat>(N-2 ,N-2);

// Populate the matrix that should be solved with the Jacobi method
    for(int i=0; i < N - 2 ; i++)
    {
       A(i,i) = d(i);
       if(i > 0)
           A(i, i -1) = e(i-1);
       if(i < N - 3)
           A(i, i +1) = e(i);
    }

// Constructing the eigenvector

    mat eigenvector = eye(N-2,N-2);

//    cout << "The operator matrix, A is:" << endl << A
//         << endl << "The eigenvector is:" << endl << eigenvector << endl;

// ///////////////////////////////////////////////////////////////////////
// Test of the algorithm compared with the armadillo algorithms
    mat B = A;
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, B);

    cout <<  "Test version gave eigenvalues: ";
    for(int i = 0 ; i < 4 ; i++)
    {
        cout << eigval(i) << ", ";
    }
    cout << endl;

//         << "and eigenvectors " << endl << eigvec ;

// ///////////////////////////////////////////////////////////////////////
    JacobiRot(&A, &eigenvector, Iterations, epsilon);

    vec Eigenvalues = sort(A.diag());       // Sorts the eigenvalues extracted from the diagonal of A

//    cout << "Finished." << endl << "The operator matrix, A is:" << endl << A
//         << endl << "The eigenvector is:" << endl << eigenvector << endl;

    cout <<  "My eigenvalues : ";
    for(int i = 0 ; i < 4 ; i++)
    {
        cout << Eigenvalues(i) << ", ";
    }
    cout << endl;

/*
    ofstream myfile;
    myfile.open ("../Project_2_Gullik/Matrix_A.csv");
        for(int i = 0; i<N - 1 ; i++)
        {
            for(int j = 0; j <N - 1; j++)
            {
                myfile << A(i,j) << ",";
            }
            myfile << endl;
        }
    myfile.close();


    myfile.open ("../Project_2_Gullik/Eigenvectors.csv");
        for(int i = 0; i<eigenvector.n_rows ; i++)
        {
            for(int j = 0; j <eigenvector.n_rows; j++)
            {
                myfile << eigenvector(i,j) << ",";
            }
            myfile << endl;
        }
    myfile.close();
*/

    return 0;
}


