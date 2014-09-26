#include <iostream>
#include <fstream>
#include "armadillo"
#include <math.h>
#include <time.h>

using namespace arma;
using namespace std;

#ifndef PROJ2_H
#define PROJ2_H

int JacobiRot(mat *, mat *, int , double epsilon);

uvec MaxOffDiag(mat *);

int Rotation(mat *,mat *, uvec);

#endif
