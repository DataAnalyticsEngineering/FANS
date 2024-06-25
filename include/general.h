

#ifndef GENERAL_H_
#define GENERAL_H_

#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstring>
#include <vector>
#include <string>
#include <map>

using namespace std;

// Packages
#include "fftw3.h"  // this includes the serial fftw as well as mpi header files! See https://fftw.org/doc/MPI-Files-and-Data-Types.html
#include "fftw3-mpi.h"

#include "H5Cpp.h"

#include "reader.h"


#include <Eigen/Dense>
using namespace Eigen;

#include "sys/stat.h"
#include "mpi.h"

#endif


#ifndef FANS_MALLOC_H
#define FANS_MALLOC_H

/* Usage: V *data = FANS_malloc<V>(n); */

template <class V>
inline V* FANS_malloc(size_t n){
    //V* out = new V[n];
    V* out = (V*) aligned_alloc(4* sizeof(V), n*sizeof(V));

    return out;
}
#endif //FANS_MALLOC_H

#define VERBOSITY 0


//#define EIGEN_RUNTIME_NO_MALLOC



