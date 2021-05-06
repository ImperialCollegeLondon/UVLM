#pragma once
#include "EigenInclude.h"
#include "types.h"
#include "constants.h"
#include "geometry.h"
#include "steady.h"
#include "unsteady.h"
#include "struct_utils.h"
#include <iostream>
//#ifdef _OPENMP
    //#include "omp.h"
//#endif

#define DLLEXPORT extern "C"
#pragma warning disable 1017


namespace UVLMlin{

  extern "C" void call_der_biot_panel(double p_DerP[9],
                    double p_DerVertices[36],
                    double p_zetaP[3],
                    double p_ZetaPanel[12],
                    const double& gamma,
                    double& vortex_radius);


  extern "C" void call_biot_panel(double p_vel[3],
                  double p_zetaP[3],
                  double p_ZetaPanel[12],
                  const double& gamma,
                  double& vortex_radius);
}
