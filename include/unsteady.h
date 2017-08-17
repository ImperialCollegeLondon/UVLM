#pragma once

#include "types.h"
#include "triads.h"
#include "constants.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"
#include "matrix.h"
#include "wake.h"
#include "postproc.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Unsteady
    {
        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_rbm_velocity,
                  typename t_forces>
        void solver
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_rbm_velocity& rbm_velocity,
            t_forces& forces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
    }
}
