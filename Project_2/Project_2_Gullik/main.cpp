#include <iostream>
#include <fstream>
#include "armadillo"
#include <math.h>
#include <time.h>
#include <Proj2Lib.h>
#include <time.h>
#include <string>

using namespace std;
using namespace arma;


// ------------------------------------------------
// On the last assignment I sent in the indents dissappeared
// send me a mail if it happens again so I can fix it so you don't need to read the code unindented
// gullik.killie@gmail.com
// ------------------------------------------------



int main()
{
    int N = 750;                      //Number of steps
    int Iterations = 100000000;
    double epsilon = 0.001;
    double w_r = 0.5;                   // Used in task 2 to calculate the new potential
    int Task = 2;                     // Chooses which task to perform, 1, corresponds to task a and b, with a single electron,
                                      // 2 is Task 2 with two interacting electrons in the potential well
    cout << "Enter wanted number of steps? "<< endl;
//    cin >> n;
    cout << "Enter wanted number of iterations" << endl;
//    cin >> Iterations;

    double rho_min = 0.0;
    double rho_max = 20.0;

    double Start, Finish, MyTime, ArmaTime; // Used for calculating how long algorithms use

// Creating the vectors to populate the matrix A
    double h = (rho_max - rho_min)/N;
    vec rho = zeros(N);

    for(int i=0; i <N ; i++)
        rho(i) = rho_min + i*h;

    vec d = zeros(N);
    vec e = zeros(N);

    for(int i = 0; i < N - 2; i++)
    {
        if(Task == 1)
        {
            d[i] = 2/pow(h,2) + pow(rho(i + 1),2);
            e[i] = -1/pow(h,2);
        }
        else
        {
           d[i] = 2/pow(h,2) + pow(w_r*rho(i + 1),2) + 1.0/(rho(i + 1));
           e[i] = -1/pow(h,2);
        }

    }

cout << "The stepsize is " << h << endl;

    mat A = zeros<mat>(N-2 ,N-2); // Constructing the eigenvector and matrix
    mat eigenvector = eye(N-2,N-2);

// Populate the matrix that should be solved with the Jacobi method
    for(int i=0; i < N - 2 ; i++)
    {
       A(i,i) = d(i);
       if(i > 0)
           A(i, i -1) = e(i-1);
       if(i < N - 3)
           A(i, i +1) = e(i);
    }

//    cout << "The operator matrix, A is:" << endl << A
//         << endl << "The eigenvector is:" << endl << eigenvector << endl;

// ///////////////////////////////////////////////////////////////////////
// Test of the algorithm compared with the armadillo algorithms
    mat B = A;
    vec eigval;
    mat eigvec;

    cout << "Running Armadillo function; eig_sym" << endl;

    Start = clock();

    eig_sym(eigval, eigvec, B);               // Armadillo's symmetric eigenvalue problem solver, used to check mine against

    Finish = clock();

    ArmaTime = (Finish - Start)/ CLOCKS_PER_SEC;

    cout << "Running my implementation of the Jacobi method" << endl;

    Start = clock();

    JacobiRot(&A, &eigenvector, Iterations, epsilon);       // My symmetric eigenvalue problem solver

    Finish = clock();

    MyTime = (Finish - Start)/ CLOCKS_PER_SEC;

    vec Eigenvalues = sort(A.diag());       // Sorts the eigenvalues extracted from the diagonal of A


//    cout << "Finished." << endl << "The operator matrix, A is:" << endl << A
//         << endl << "The eigenvector is:" << endl << eigenvector << endl;



    cout <<  "Test version gave eigenvalues: ";         //Printing the eigenvalues from the armadillo function and from my function
    for(int i = 0 ; i < 4 ; i++)
    {
        cout << eigval(i) << ", ";
    }
    cout << endl;

    cout <<  "My eigenvalues : ";
    for(int i = 0 ; i < 4 ; i++)
    {
        cout << Eigenvalues(i) << ", ";
    }
    cout << endl;

    cout << "The Armadillo eig_sym used " << ArmaTime << "s, while my implementation of the Jacobi method used " <<
            MyTime << "s, to solve a  " << N << "x" << N << "eigenvalue problem" << endl;

// Checking my eigenvectors compared to the armadillo function's
//    cout << "Armadillo got eigvec:" << endl <<eigvec << endl << "Jacobi_Rot got eigenvector " << endl <<
//                  eigenvector << endl; //<< "The difference is" << endl << eigvec - eigenvector << endl;


    // ///////////////////////////////////////////////////////////////////////////
    // To plot the eigenvector corresponding to the lowest eigenvalue it first finds
    // the position of the lowest and then stores it in a csv file
    // which is plotted by a python script
    // //////////////////////////////////////////////////////////////////////////

    vec temp_eigenvalues = diagvec(A);
    double min, temp_eigenvalue;;
    vec temp_row;

    for(int i = 0; i < temp_eigenvalues.n_elem ;i++)        //Algorithm to sort the eigenvectors corresponding to ascending eigenvalues
    {
        min = 10000;

        for(int j = i; j < temp_eigenvalues.n_elem ;j++)
        {


            if(temp_eigenvalues(j) < min)
            {
                min = temp_eigenvalues(j);

                temp_row = eigenvector.col(i);
                eigenvector.col(i) = eigenvector.col(j);
                eigenvector.col(j) = temp_row;

//                temp_eigenvalue =
            }
        }

    }

    ofstream myfile;        //Storing the eigenvectors as a csv file to be plotted in a python script

    myfile.open ("../Project_2_Gullik/Eigenvectors_w_r_0.5.csv");
        for(int i = 0; i<eigenvector.n_rows ; i++)
        {
            for(int j = 0; j <eigenvector.n_rows; j++)
            {
                myfile << eigenvector(i,j) << ",";
            }
            myfile << endl;
        }
    myfile.close();

    return 0;
}


