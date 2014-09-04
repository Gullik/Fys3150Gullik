#include <iostream>
#include "armadillo"

using namespace std;
using namespace arma;

void Tridim_Solver(int N)
{
    ivec a = -1*ones<ivec>(N);
    ivec b = 2*ones<ivec>(N);
    ivec c = -1*ones<ivec>(N);

    vec u(N);

    cout << u.n_elem << endl;




//    cout << "Vectoren a er " << endl << a  << endl;
//    cout << "Vectoren b er " << endl << b  << endl;
//    cout << "Vectoren c er " << endl << c  << endl;

    int i;

    for(i=0 ; i<N ; i++)
    {

        cout << "u(" << i << ") = " << u(i) << endl;

    }

    return;
}

int main()
{
   int N;

   cout << "Number of steps:";

   cin >> N;

   Tridim_Solver(N);





//   cout << "Vectoren er " << endl << a  << endl << "Lengden er " << b << endl;

//   return 0;




    return 0;
}

