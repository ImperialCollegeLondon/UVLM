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
        template<typename t_zeta_phantom,
                 typename t_zeta_star,
                 typename t_zeta_phantom_star>
        void create_phantom_zeta_star
        (
            t_zeta_phantom& zeta_phantom,
            t_zeta_star& zeta_star,
            t_zeta_phantom_star& zeta_phantom_star
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
    // Parameter initialisation
    const uint n_surf = zeta.size();
    int N_col, i_col;
    double phantom_dy=0.0;  
    double ratio_vector_length = 0.0;
    bool phantom_surface = false;
    bool phantom_surface_already_set = false;
    zeta_phantom.resize(n_surf);
    uint N_phantom_cell, N_row_phantom, index_phantom_0, index_phantom_1;
    uint i_surf_partner_junction = 0;
    UVLM::Types::Vector3 junction_point_0, junction_point_1, junction_connection_vector;
    UVLM::Types::Vector3 interpolated_point_0, interpolated_point_1;

    // Resize surfaces beforehand to avoid problems with empty phantom surfaces
    for(uint i_surf=0; i_surf<n_surf; i_surf++)
    {
        zeta_phantom[i_surf].resize(3);
        for(uint i_dim=0; i_dim<3; ++i_dim)
        {
            zeta_phantom[i_surf][i_dim].resize(0, 0);
        }
    }
    for(uint i_surf=0; i_surf<n_surf; i_surf++)
    {
        phantom_surface = false;        
        N_row_phantom = zeta[i_surf][0].rows();
        //To-Do: Check why -1 necessary (bug probably in datastructres.py or aerogrid.py)
        for(uint index_phantom=0;index_phantom<zeta[i_surf][0].cols()-1;++index_phantom)
        {            
            // if partner junnction surface in flag lower than current surface
            // than the surface has already been created
            if (flag_zeta_phantom[i_surf](index_phantom, 0) > i_surf)
            {
                i_surf_partner_junction = flag_zeta_phantom[i_surf](index_phantom, 0);
                phantom_dy= zeta[i_surf][1](0, index_phantom+1)-zeta[i_surf][1](0, index_phantom);
                N_col = abs(round((zeta[i_surf][1](0, index_phantom)-zeta[i_surf_partner_junction][1](0, index_phantom))/(2.0*phantom_dy)));
                
                for(uint i_dim=0; i_dim<3; ++i_dim)
                {
                    zeta_phantom[i_surf][i_dim].resize(N_row_phantom, N_col+1);
                    zeta_phantom[i_surf_partner_junction][i_dim].resize(N_row_phantom, N_col+1);
                    zeta_phantom[i_surf][i_dim].block(0, 0,N_row_phantom, 1) = zeta[i_surf][i_dim].block(0, index_phantom,N_row_phantom, 1);
                    zeta_phantom[i_surf_partner_junction][i_dim].block(0, 0,N_row_phantom, 1) = zeta[i_surf_partner_junction][i_dim].block(0, index_phantom,N_row_phantom, 1);
                }
                phantom_surface = true;
                break;
            }
        }
        if (phantom_surface)
        {
            // Dynamic phantom surface creation depending on vectors between junction surfaces
            for(uint i_row=0;i_row<N_row_phantom;++i_row)
            {
                // get vector between boundaries
                junction_point_0 << zeta_phantom[i_surf][0](i_row, 0), zeta_phantom[i_surf][1](i_row, 0), zeta_phantom[i_surf][2](i_row, 0);
                junction_point_1 << zeta_phantom[i_surf_partner_junction][0](i_row, 0), zeta_phantom[i_surf_partner_junction][1](i_row, 0), zeta_phantom[i_surf_partner_junction][2](i_row, 0);
                junction_connection_vector = junction_point_0 - junction_point_1;
                for(uint i_col=0;i_col<N_col+1;++i_col)
                {
                    ratio_vector_length =double(N_col - i_col) / double(N_col);
                    interpolated_point_0 = junction_point_0 - (0.5 * ratio_vector_length ) * junction_connection_vector;
                    interpolated_point_1 = junction_point_1 + (0.5 * ratio_vector_length ) * junction_connection_vector;
                    
                    for(uint i_dim=0; i_dim<3; ++i_dim)
                    {
                        zeta_phantom[i_surf][i_dim](i_row, i_col) = interpolated_point_0(i_dim);                    
                        zeta_phantom[i_surf_partner_junction][i_dim](i_row, i_col) = interpolated_point_1(i_dim);
                    }
                }
            }
        }
    }
}

template<typename t_zeta_phantom,
            typename t_zeta_star,
            typename t_zeta_phantom_star>
void UVLM::Phantom::create_phantom_zeta_star
(
    t_zeta_phantom& zeta_phantom,
    t_zeta_star& zeta_star,
    t_zeta_phantom_star& zeta_phantom_star
)
{
    // allocate zeta phantom star
    uint N_rows, N_cols, N_rows_zeta_phantom;
    const uint n_surf_phantom = zeta_phantom.size();
    zeta_phantom_star.resize(n_surf_phantom);
    for (uint i_surf=0; i_surf <n_surf_phantom; ++i_surf)
    {
        zeta_phantom_star[i_surf].resize(UVLM::Constants::NDIM);
        N_cols = zeta_phantom[i_surf][0].cols();
        N_rows = zeta_phantom[i_surf][0].rows();
        for(uint i_dim=0; i_dim<UVLM::Constants::NDIM;++i_dim)
        {
            zeta_phantom_star[i_surf][i_dim].resize(N_rows, N_cols);
        }
        if (zeta_phantom_star[i_surf][0].size() > 0)
        {
            // set coordinates of zeta phantom star
            for (uint i_row=0; i_row<N_rows;++i_row)
            {
                for (uint i_col=0; i_col<N_cols;++i_col)
                {
                    // TO-DO: Generalize case, or outsource whole phantom panels generation to SHARPy
                    zeta_phantom_star[i_surf][1](i_row, i_col) = zeta_phantom[i_surf][1](N_rows-1, i_col);
                    zeta_phantom_star[i_surf][0](i_row, i_col) = zeta_star[i_surf][0](i_row, 0);
                    zeta_phantom_star[i_surf][2](i_row, i_col) = zeta_star[i_surf][2](i_row, 0);
                }
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
        for(uint i_row=0; i_row < flag_zeta_phantom[i_surf].rows();++i_row)
        {
            if (flag_zeta_phantom[i_surf](i_row, 0))
            {
                return true;
            }
        }
    }
    return false;
}