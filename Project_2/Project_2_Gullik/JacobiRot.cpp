#include <Proj2Lib.h>

uvec MaxOffDiag(mat * A)
{

   uvec Loc = zeros<uvec>(2);
   uint MatLength = A->n_rows;
   double Max = 0.0;

   mat B = mat(A->memptr(),A->n_rows,A->n_cols,false);  //This creates a matrix B that shares the same memory as the matrix A,
                                                        //so changes to that is also done to A and it is easier to get matrix elements out of it



    for(uint i = 0 ; i < MatLength; i++)
    {
        for(uint j = i + 1 ; j < MatLength; j++)
        {
            if(fabs(B(i,j)) >= Max)
            {
                Max = fabs(B(i,j));
//                cout << "i er: "<< i << "       j er: "<<j << endl;
                Loc(0) = i;
                Loc(1) = j;
            }
        }

    }

    return Loc;
}

int Rotation(mat * A, mat * eigenvector, uvec Location)      // Do this by hand for a simple matrix and compare results
{
    int x = Location(0);    //Row-number
    int y = Location(1);    //Column-number

    mat B = mat(A->memptr(),A->n_rows,A->n_cols,false);
    mat R = mat(eigenvector->memptr(),eigenvector->n_rows,eigenvector->n_cols,false);
                                                        //This creates a matrix B that shares the same memory as the matrix A,
                                                        //so changes to that is also done to A and it is easier to get matrix elements out of it
                                                        //and all operations done to B also applies to A.
                                                        //This also avoids any copying of data to save time

    // Calculating the angle
    double c, s,t, tau;
    if(B(x,y) != 0)
    {
        tau = (B(y,y) - B(x,x))/(2.0*B(x,y));
        if(tau > 0)
            t = -tau + sqrt(1.0 + pow(tau,2));
        else
            t = -tau - sqrt(1.0 + pow(tau,2));
    }

    c = 1.0/sqrt(1+pow(t,2));
    s = t*c;

    cout << "tau er: " << tau << ", t er: " << t << ", c er: " << c << ", og s er: " << s << endl;

//Making the rotation matrix S, a matrix that should rotate the elements xx,xy,yx, and yy by tau
//Just for testing purposes to see that the quicker algorithm is correct, since this is quite slow
//It also produced some numerical small errors, probably because of approximations in the inversion

    mat Rotation = eye(B.n_rows, B.n_rows);

    Rotation(x,x) = c;
    Rotation(x,y) = -s;
    Rotation(y,x) = s;
    Rotation(y,y) = c;

//    cout << "S is:" << endl << Rotation << endl;

    mat Test = Rotation*B*inv(Rotation);

    if(trace(inv(Rotation) ) > trace(Rotation.t()) + 0.01  ||       // Checks to se if the angles are correct
            trace(inv(Rotation) ) < trace(Rotation.t()) - 0.01)     // by checking that the transposed and inverse have the same trace
    {
        cout << "Error: the similarity matrix is not a similarity matrix" << endl;
        exit(0);
        return 0;
    }

    cout << "A is : " << endl << *A << endl << "After Rotation:"  <<  endl <<Test << endl;

// //////////////////////////// Testing done

    double bxx = B(x,x);
    double byy = B(y,y);
    double bxy = B(x,y);

    B(x,x) = byy*pow(c,2) + bxx*pow(s,2) - 2.0*bxy*c*s;    // Assuming that it is a symmetric matrix so A(x,y)=A(y,x)
    B(y,y) = byy*pow(c,2) + bxx*pow(s,2) + 2.0*bxy*c*s;    // x is the row position and y is the column position,
    B(x,y) = (bxx-byy)*s*c + bxy*(pow(c,2)-pow(s,2));      // These can be
    B(y,x) = (bxx- byy)*s*c + bxy*(pow(c,2)-pow(s,2));     // set to 0

    // Changing the remaining elements by first calculating the uppr triangular part
    // then since we started with a symmetric A, and the similarity transformation keeps
    // the symmetri copying the lower triangular part

    double bix,biy;

    for(int i=0; i < abs(B.n_rows); i++)
    {
        if(i!=x && i!=y)
        {
 //           cout << "i is:" << i <<"; x is: " << x << "; y is: " << y << endl;
            bix = B(i,x);
            biy = B(i,y);
            B(i,x) = c*bix - s*biy;
            B(i,y) = c*biy + s*bix;
            B(x,i) = B(i,x);
            B(y,i) = B(i,y);
        }
    }

    cout << B << endl;

// Time to calculate the new eigenvector that has been operated on by S^(-1)

    return 0;
}






int JacobiRot(mat * A, mat * eigenvector, int Iterations, double epsilon)
{

    uvec MaxLoc = zeros<uvec>(2);

    for(int i = 0; i < Iterations; i++)
    {
        cout << "Doing Jacobi rotation number " << i << endl;

        MaxLoc = MaxOffDiag (A);

        cout << "Biggest matrix element is at position: " << MaxLoc[0] << "," << MaxLoc[1] << " and has the value " << A[0,0](MaxLoc(0),MaxLoc(1)) <<endl;

        if(fabs(A[0,0](MaxLoc(0),MaxLoc(1))) < epsilon)         //The A[0,0](i,j) gets goes to the matrix stored in the pointer and gets element ij
        {
            cout << "The max offdiagonal element is below the threshold, " << epsilon << ", at rotation " << i << endl;
            return 0;
        }
        Rotation(A, eigenvector, MaxLoc);
    }



    return 0;

}
