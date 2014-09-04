#include <iostream>
#include "armadillo"


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

    //Creating some points of the source term f(x) = 100e^{-10x}
    vec f = exp(-10 * linspace(0,1, N));

//    cout << f << endl;

    //The forward subsitution step, which gives u(i) + cst*u(i+1) = cst*f(i)

    double btemp = b;
    u(0) = f(0)/b;

    for(i=1 ; i < N ; i++)
    {
        divisor_temp(i) = c/btemp;
        btemp = (b - (a*divisor_temp(i)));
        u(i) = (f(i) - (a*u(i-1)))/btemp;
    }

    // The backwards step;
    for ( i = N-2 ; i>=0 ; i--)
    {
        u(i) -= divisor_temp(i)*u(i+1);
    }



//    for(int i=0 ; i<N ; i++)
//    {
//        cout << "u(" << i << ") = " << u(i)   <<endl;
//    }

    return u;
}

int main()
{
   int N;

   cout << "Number of steps:";

   cin >> N;

   vec x = linspace(0,1,N);

   vec u = Tridim_Solver(N);

   vec Exact = 1 - (1-exp(-10))*x - exp(-10*x);







//   cout << "Vectoren er " << endl << a  << endl << "Lengden er " << b << endl;

//   return 0;




    return 0;
}

