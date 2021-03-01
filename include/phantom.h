#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "debugutils.h"

#include <fstream>

namespace UVLM
{
    namespace Phantom
    {
        template<typename t_zeta,
                 typename t_zeta_phantom,
                 typename t_flag_zeta_phantom>
        void create_phantom_zeta
        (
            t_zeta& zeta,
            t_zeta_phantom& zeta_phantom,
            t_flag_zeta_phantom& flag_zeta_phantom
        );
        template<typename t_flag_zeta_phantom>
        bool check_for_true_in_bool_vec_mat
        (
            t_flag_zeta_phantom& flag_zeta_phantom
        );
    }
}

template<typename t_zeta,
         typename t_zeta_phantom,
            typename t_flag_zeta_phantom>
void UVLM::Phantom::create_phantom_zeta
(
    t_zeta& zeta,
    t_zeta_phantom& zeta_phantom,
    t_flag_zeta_phantom& flag_zeta_phantom
)
{
    const uint n_surf = zeta.size();
    uint N_col;
    float phantom_dy=0.0;
    zeta_phantom.resize(n_surf);
    unsigned int N_phantom_cell, N_row_phantom, index_phantom;
    for(uint i_surf=0; i_surf<n_surf; i_surf++)
    {
        N_phantom_cell = (flag_zeta_phantom[i_surf] > 0).count();
        N_row_phantom = zeta[i_surf][0].rows();
 
        zeta_phantom[i_surf].resize(3);
        if(N_phantom_cell != 0)
        {
            std::vector<Eigen::Index> idxs;
            for(Eigen::Index i=0; i<flag_zeta_phantom[i_surf].size(); ++i)
            {
                if(flag_zeta_phantom[i_surf](i))
                {
                    // To-DO: Here only one junction is assumeds
                    index_phantom = i;
                    break;
                }

            }
            
            // Required phantom panels to extend phantom panels to symmetry plane
           N_col = abs(round(zeta[i_surf][1](0, index_phantom)/(zeta[i_surf][1](0, index_phantom+1)-zeta[i_surf][1](0, index_phantom))));
            if(N_col==1)
            {
                phantom_dy=zeta[i_surf][1](0, index_phantom+1)-zeta[i_surf][1](0, index_phantom);
            }
            else
            {
                phantom_dy=zeta[i_surf][1](0, index_phantom)/N_col;
            }
            for(uint i_dim=0; i_dim<3; ++i_dim)
            {
                zeta_phantom[i_surf][i_dim].resize(N_row_phantom, N_col+1);
                
                zeta_phantom[i_surf][i_dim].block(0, 0,N_row_phantom, 1) = zeta[i_surf][i_dim].block(0, index_phantom,N_row_phantom, 1);
                
                for(uint i_col = 1;  i_col<N_col+1;++i_col)
                {
                    
                    if(i_dim==1)
                    {
                       //TO-DO:
                       // - Symmetry plane is used for end of panel
                       // - Here symmetry plane is assumed to be at y = 0 --> adjust e.g. if sideslip angle != 0
                       
                       zeta_phantom[i_surf][i_dim].col(i_col) = zeta_phantom[i_surf][i_dim].col(i_col-1);
                       for(uint i_row=0;i_row<N_row_phantom;++i_row)
                       {
                           zeta_phantom[i_surf][i_dim](i_row, i_col)-=phantom_dy;
                        }
                    }
                    else
                    {      
                        // TO-DO:
                        // - here only unswept wings are covered
                        zeta_phantom[i_surf][i_dim].col(i_col) = zeta_phantom[i_surf][i_dim].col(index_phantom);
                    }
                }
            }
        
        }
        else
        {
            for(uint i_dim=0; i_dim<3; ++i_dim)
            {
            zeta_phantom[i_surf][i_dim].resize(0,0);
            }
        }

    }
}


template<typename t_flag_zeta_phantom>
bool UVLM::Phantom::check_for_true_in_bool_vec_mat
(
    t_flag_zeta_phantom& flag_zeta_phantom
)
{
    for(uint i_surf=0; i_surf<flag_zeta_phantom.size(); ++i_surf)
    {
        for(uint i_row=0; i_row < flag_zeta_phantom[i_surf].size();++i_row)
        {
            if (flag_zeta_phantom[i_surf][i_row])
            {
                return true;
            }
        }
    }
    return false;
}