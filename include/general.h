

#ifndef GENERAL_H_
#define GENERAL_H_

#include <cmath>
#include <complex>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// Packages
#include "fftw3-mpi.h"
#include "fftw3.h" // this includes the serial fftw as well as mpi header files! See https://fftw.org/doc/MPI-Files-and-Data-Types.html

#include "H5Cpp.h"

#include "reader.h"

#include <Eigen/Dense>
using namespace Eigen;

#include "mpi.h"
#include "sys/stat.h"

#endif

#ifndef FANS_MALLOC_H
#define FANS_MALLOC_H

/* Usage: V *data = FANS_malloc<V>(n); */

template <class V>
inline V *FANS_malloc(size_t n)
{
    // V* out = new V[n];
    V *out = (V *) aligned_alloc(4 * sizeof(V), n * sizeof(V));

    return out;
}
#endif // FANS_MALLOC_H

#define VERBOSITY 0

//#define EIGEN_RUNTIME_NO_MALLOC
