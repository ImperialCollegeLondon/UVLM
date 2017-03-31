#pragma once

#include "EigenInclude.h"
#include "types.h"

namespace UVLM
{
    namespace BiotSavart
    {
        // DECLARATIONS
        template <typename t_zeta,
                  typename t_gamma,
                  typename t_ttriad,
                  typename t_uout>
        void surface
        (
            const t_zeta&       zeta,
            const t_gamma&      gamma,
            const t_ttriad&     target_triad,
            const t_uout&       uout,
            unsigned int        Mstart = 0,
            unsigned int        Nstart = 0,
            unsigned int        Mend   = -1,
            unsigned int        Nend   = -1,
            const bool&         image_method = false
        );
    }
}



// SOURCE CODE
template <typename t_zeta,
          typename t_gamma,
          typename t_ttriad,
          typename t_uout>
void UVLM::BiotSavart::surface
(
    const t_zeta&       zeta,
    const t_gamma&      gamma,
    const t_ttriad&     target_triad,
    const t_uout&       uout,
    unsigned int        Mstart = 0,
    unsigned int        Nstart = 0,
    unsigned int        Mend   = -1,
    unsigned int        Nend   = -1,
    const bool&         image_method = false
)
{
    // If Mend or Nend are == -1, their values are taken as the surface M and N
    if (Mend == -1) {Mend = zeta[0].rows();}
    if (Nend == -1) {Nend = zeta[0].cols();}

    






}
