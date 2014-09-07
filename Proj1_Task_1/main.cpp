#include <iostream>
#include <fstream>
#include "armadillo"
#include <math.h>
#include <time.h>





using namespace std;
using namespace arma;


vec Tridim_Solver(int N)
{
    int a = -1; // The three constants for the three diagonals in the
    int b = 2;  // matrix
    int c = -1;

    vec divisor_temp(N);    // Used in the forward subsitution step,
                            // predefined here to avoid redefining it several times in the for loop

    vec u(N);
    int i;

    //Creating some points of the source term f(x) = 100e^{-10x}*h²
    vec tilde_b = 100*exp(-10 * linspace(0,1, N))/pow((1+N), 2);



//    cout << f << endl;

    //The forward subsitution step, which gives u(i) + cst*u(i+1) = cst*f(i)

    double btemp = b;
    u(0) = tilde_b(0)/b; //f(0) not being 0 is a problem




    for(i=1 ; i < N ; i++)
    {
        divisor_temp(i) = c/btemp;
        btemp = (b - (a*divisor_temp(i)));
        u(i) = (tilde_b(i) - (a*u(i-1)))/btemp;
    }

    // The backwards step;
    for ( i = N-2 ; i>=0 ; i--)
    {
        u(i) -= divisor_temp(i)*u(i+1);
    }


// Checking that the values that get out makes any sense
//    for(int i=0 ; i<N ; i++)
//    {
//        cout << "u(" << i << ") = " << u(i)   <<endl;
//    }

    return u;
}

mat RelError(int N)
{

    //The matrix should only use int for the first row and double for the second row, to be fixed later
    mat ErrorTable = zeros(N/10,2);

    for(int i = 10; i <= N; i+=10)
    {
        vec x = linspace(0,1,i);

        vec u = Tridim_Solver(i);

        vec Exact = 1 - (1-exp(-10))*x - exp(-10*x);

        // The first value and last value of the error is cut out since it results in an infinity because the exact solution
        // is 0 at u(0) and u(N)

       vec epsilon = log10(abs((u-Exact)/Exact));
       double MaxRelError = max(epsilon.subvec(1, epsilon.n_elem -2));

        ErrorTable(i/10 - 1, 1) = MaxRelError;
        ErrorTable(i/10 - 1, 0) = i;


    }

    return ErrorTable;
}

int Time_Comp(int N)
{
    //Since the specialised tridimensional solver doesn't need the matrix to be
    // created the matrix creation is included for the general armadillo solver.

clock_t Start, Finish;

Start = clock();

cout << "Starttiden er " <<Start << endl;

    //Creating the A matrix,
    mat A = zeros<mat>(N,N);

    for(int i = 0 ; i<N ; i++)
    {
        for(int j = 0 ; j < N ; j++)
        {
            if(j==i)
                A(i,j) = 2;
            if(j==i+1)
                A(i,j) = -1;
            if(j==i-1)
                A(i,j)= -1;
        }
    }


    //Creating source term
    vec tilde_b = 100*exp(-10 * linspace(0,1, N))/pow((1+N), 2);

    vec u_solver = zeros<vec>(N);

    u_solver = solve(A, tilde_b);

Finish = clock();

cout << "Slutttid " << Finish << endl;

double GenSolv = (Finish - Start)/CLOCKS_PER_SEC;


Start = clock();
    vec u_mysolver = Tridim_Solver(N);
Finish = clock();

double MySolv = (Finish - Start)/CLOCKS_PER_SEC;

cout << GenSolv << " vs  " << MySolv << endl;


    return 0;
}

int main()
{


   int N;

   cout << "Number of steps:";

   cin >> N;

   vec x = linspace(0,1,N);

   vec u = Tridim_Solver(N);

   vec Exact = 1 - (1-exp(-10))*x - exp(-10*x);

//   cout << u << endl;
//   cout << Exact << endl;

   ofstream myfile;
   myfile.open ("../Proj1_Task_1/Results.csv");
       for(int i = 0; i<N ; i++)
       {
            myfile <<x(i) << "," << u(i) << "," <<x(i) << "," << Exact(i) << endl ;
       }
   myfile.close();

// Error calculation
cout << "Want to do relative error?" << endl;

    int Ans;
    cin >> Ans;
    if (Ans == 1)
     {   mat ErrorTable = RelError(N);



//   ofstream myfile;
   myfile.open ("../Proj1_Task_1/ErrorTable.csv");
       for(uint i = 0; i< ErrorTable.n_rows ; i++)
       {
            myfile <<ErrorTable(i,0) << "," << ErrorTable(i,1) << endl ;
       }
   myfile.close();
    }
//Done with error calculation part

   Time_Comp(N);


  
    return 0;
}
