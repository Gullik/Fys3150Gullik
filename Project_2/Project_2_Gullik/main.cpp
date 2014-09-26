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
    int N = 4;
    int Iterations = 10;
    double epsilon = 0.0001;
    cout << "Enter wanted number of steps? "<< endl;
//    cin >> N;
    cout << "Enter wanted number of iterations" << endl;
//    cin >> Iterations;


    double rho_min = 0;
    double rho_max = 100;

// Creating the vectors to populate the matrix A
    double h = (rho_max - rho_min)/N;

    vec d = zeros(N);
    vec e = zeros(N);

    for(int i = 0; i < N; i++)
    {
        d[i] = 2; //2/pow(h,2) + i*h;
        e[i] = 1; //-1/pow(h,2);
    }

cout << "The stepsize is " << h << endl;

    mat A = zeros<mat>(N,N);

// Populate the matrix that should be solved with the Jacobi method
    for(int i=0; i < N ; i++)
    {
       A(i,i) = d(i);
       if(i > 0)
           A(i, i -1) = e(i-1);
       if(i < N-1)
           A(i, i +1) = e(i);
     }

// Constructing the eigenvector

    mat eigenvector = eye(N,N);



    JacobiRot(&A, &eigenvector, Iterations, epsilon);


    return 0;
}


