//====================================================================================================================
//
//  CSpack_main.cpp
//
//  Created by Abir Sarkar on 21/01/2020 with modifications from Jens Chluba.
//  This code illustrates some of the computations with CSpack
//
//====================================================================================================================
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;
using namespace CSpack_kernel_moments;

//====================================================================================================================
int main(int narg, char *args[])
{
    double omega0=0.01, p0=0.1;

    cout << compute_exact_moments_analytical(omega0, p0, 0) << endl;
    cout << compute_exact_moments_analytical(omega0, p0, 1) << endl;
    cout << compute_exact_moments_analytical(omega0, p0, 2) << endl;

    return 0;
}

//====================================================================================================================
//====================================================================================================================
