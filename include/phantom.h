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
        template<typename t_flag_phantom,
                 typename t_zeta_phantom,
                 typename t_zeta_star,
                 typename t_zeta_phantom_star>
        void create_phantom_zeta_star
        (
            const t_flag_phantom& flag_zeta_phantom,
            t_zeta_phantom& zeta_phantom,
            t_zeta_star& zeta_star,
            t_zeta_phantom_star& zeta_phantom_star
        );
        template<typename t_flag_phantom,
                 typename t_flag_zeta_phantom>
        bool check_for_true_in_bool_vec_mat
        (
            t_flag_zeta_phantom& flag_zeta_phantom
        );
        template<typename t_flag_phantom>
        void get_parameter_phantom_setup
        (
            const t_flag_phantom& flag_zeta_phantom,
            const uint& i_surf,
            uint& i_surf_partner_junction,
            uint& idx_junction,
            bool& phantom_surface
        );
        template<typename t_zeta_out,
                 typename t_zeta_in>
        void interpolate_geometry_coordinates
        (
            t_zeta_out& zeta_out,
            t_zeta_in& zeta_in,
            const uint& N_row,
            const uint& N_col,
            const uint& idx_junction,
            const uint& i_surf,
            const uint& i_surf_partner_junction
        );
        template<typename t_gamma_out,
                typename t_gamma_in,
                typename t_zeta_col_out,
                typename t_zeta_col_in>
        void interpolate_circulation_strength
        (
            t_gamma_out& gamma_out,
            const t_gamma_in& gamma_in,
            const t_zeta_col_out& zeta_col_out,
            const t_zeta_col_in& zeta_col_in,
            const uint& idx_junction,
            const uint& i_surf,
            const uint& i_surf_partner_junction
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
    uint idx_junction= 0;
    double phantom_dy=0.0;  
    bool phantom_surface = false;
    uint  N_row_phantom;
    uint i_surf_partner_junction = 0;

    // Allocate surfaces    
    zeta_phantom.resize(n_surf);

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
        UVLM::Phantom::get_parameter_phantom_setup(flag_zeta_phantom, i_surf,i_surf_partner_junction,idx_junction,phantom_surface);
        if (phantom_surface)
        {
            N_row_phantom = zeta[i_surf][0].rows();
            // ToDo: Fix since following line only works if junction is not in last column        
            phantom_dy= zeta[i_surf][1](0, idx_junction+1)-zeta[i_surf][1](0, idx_junction);
            N_col = abs(round((zeta[i_surf][1](0, idx_junction)-zeta[i_surf_partner_junction][1](0, idx_junction))/(2.0*phantom_dy)));
            interpolate_geometry_coordinates(zeta_phantom, zeta, N_row_phantom, N_col+1, idx_junction, i_surf, i_surf_partner_junction);
        }
        
    }
}

template<typename t_flag_phantom,
        typename t_zeta_phantom,
            typename t_zeta_star,
            typename t_zeta_phantom_star>
void UVLM::Phantom::create_phantom_zeta_star
(
    const t_flag_phantom& flag_zeta_phantom,
    t_zeta_phantom& zeta_phantom,
    t_zeta_star& zeta_star,
    t_zeta_phantom_star& zeta_phantom_star
)
{
    // allocate zeta phantom star
    uint N_rows, N_cols, N_rows_zeta_phantom, i_surf_partner_junction, idx_junction = 0;
    const uint n_surf_phantom = zeta_phantom.size();
    zeta_phantom_star.resize(n_surf_phantom);
    
    bool uninitiliased_phantom_surface = false;
    // Surface allocation has to be done before the coordinate interpolation
    // because of the simultaneous generation of the partner surface
    for (uint i_surf=0; i_surf <n_surf_phantom; ++i_surf)
    {
        zeta_phantom_star[i_surf].resize(UVLM::Constants::NDIM);
    }
    for (uint i_surf=0; i_surf <n_surf_phantom; ++i_surf)
    {
        N_cols = zeta_phantom[i_surf][0].cols();
        N_rows_zeta_phantom = zeta_phantom[i_surf][0].rows();
        N_rows = zeta_star[i_surf][0].rows();
        if (zeta_phantom[i_surf][0].size() > 0)
        {   
            UVLM::Phantom::get_parameter_phantom_setup(flag_zeta_phantom, i_surf,i_surf_partner_junction,idx_junction, uninitiliased_phantom_surface);
            // Note: phantom surface can be false, if phantom surface is already initiliased
            if (uninitiliased_phantom_surface)
            {              
                // set coordinates of zeta phantom star
                interpolate_geometry_coordinates(zeta_phantom_star, zeta_star, N_rows, N_cols, idx_junction, i_surf, i_surf_partner_junction);
            
                // coordinates of phantom surface and wake are coincident at trailing edge
                for (uint i_dim=0; i_dim<3; i_dim++)
                {
                    zeta_phantom_star[i_surf][i_dim].topRows(1) = zeta_phantom[i_surf][i_dim].bottomRows(1);
                }
            }
        }
        else
        {
            for (uint i_dim=0; i_dim<3; i_dim++)
            {
                zeta_phantom_star[i_surf][i_dim].resize(0,0);
    }
}
    }

}

template<typename t_flag_phantom>
void UVLM::Phantom::get_parameter_phantom_setup
(
    const t_flag_phantom& flag_zeta_phantom,
    const uint& i_surf,
    uint& i_surf_partner_junction,
    uint& idx_junction,
    bool& phantom_surface
)
{

    // if partner junnction surface in flag lower than current surface
    // then the surface has already been created
    phantom_surface = false;
    if (flag_zeta_phantom(0, i_surf) >= 0)
    {
        if (flag_zeta_phantom(0, i_surf) > i_surf)
        {
            i_surf_partner_junction = flag_zeta_phantom(0, i_surf);
            idx_junction = 0;
            phantom_surface = true;
        }
    }

}


template<typename t_zeta_out,
         typename t_zeta_in>
void UVLM::Phantom::interpolate_geometry_coordinates
(
    t_zeta_out& zeta_out,
    t_zeta_in& zeta_in,
    const uint& N_row,
    const uint& N_col,
    const uint& idx_junction,
    const uint& i_surf,
    const uint& i_surf_partner_junction

)
{
    UVLM::Types::Vector3 junction_point_0, junction_point_1, junction_connection_vector,
                         interpolated_point_0, interpolated_point_1;
    double ratio_vector_length = 0.0;
    // Allocate
    for(uint i_dim=0; i_dim<3; ++i_dim)
    {
        zeta_out[i_surf][i_dim].resize(N_row, N_col);
        zeta_out[i_surf_partner_junction][i_dim].resize(N_row, N_col);
    }
    // Dynamic phantom surface creation depending on vectors between junction surfaces
    for(uint i_row=0;i_row<N_row;++i_row)
    {
        // get vector between boundaries
        // ToDO: Create function to read Vector from i_surf, i_row, i_col
        junction_point_0 << zeta_in[i_surf][0](i_row, idx_junction), zeta_in[i_surf][1](i_row, idx_junction), zeta_in[i_surf][2](i_row, idx_junction);
        junction_point_1 << zeta_in[i_surf_partner_junction][0](i_row, idx_junction), zeta_in[i_surf_partner_junction][1](i_row, idx_junction), zeta_in[i_surf_partner_junction][2](i_row, idx_junction);
        junction_connection_vector = junction_point_0 - junction_point_1;
        for(uint i_col=0;i_col<N_col;++i_col)
        {
            ratio_vector_length =double(N_col - 1 - i_col) / double(N_col - 1);
            interpolated_point_0 = junction_point_0 - (0.5 * ratio_vector_length ) * junction_connection_vector;
            interpolated_point_1 = junction_point_1 + (0.5 * ratio_vector_length ) * junction_connection_vector;
            
            for(uint i_dim=0; i_dim<3; ++i_dim)
            {
                zeta_out[i_surf][i_dim](i_row, i_col) = interpolated_point_0(i_dim);                    
                // ToDo: Check if at symmetry plane coordinates of both surfaces are the same!                                 
                zeta_out[i_surf_partner_junction][i_dim](i_row, i_col) = interpolated_point_1(i_dim);
            }
        }
    }
}
template<typename t_gamma_out,
         typename t_gamma_in,
         typename t_zeta_col_out,
         typename t_zeta_col_in>
void UVLM::Phantom::interpolate_circulation_strength
(
    t_gamma_out& gamma_out,
    const t_gamma_in& gamma_in,
    const t_zeta_col_out& zeta_col_out,
    const t_zeta_col_in& zeta_col_in,
    const uint& idx_junction,
    const uint& i_surf,
    const uint& i_surf_partner_junction
)
{
    
    // To-Do: - final interpolation in function to avoid double code for i_surf and partner surface
    //        - Create function to read Vector from i_surf, i_row, i_col
    //        - avoid abs(gamma) and then multpzing by -1 (reargannge surface vectors?)
    UVLM::Types::Vector3 junction_point_0, junction_point_1,
                         collocation_point_phantom, collocation_point_phantom_partner;
       uint N_col, N_row, updated_surf, partner_surf;
    double distance_0, distance_1 = 0.0;

    // Allocate
    N_row = zeta_col_out[i_surf][0].rows();
    N_col =  zeta_col_out[i_surf][0].cols();
    gamma_out[i_surf].resize(N_row, N_col);
    gamma_out[i_surf_partner_junction].resize(N_row, N_col);
    
    // Dynamic phantom surface creation depending on vectors between junction surfaces
    for(uint i_row=0;i_row<N_row;++i_row)
    {
        // get vector between boundaries
        junction_point_0 << zeta_col_in[i_surf][0](i_row, idx_junction), zeta_col_in[i_surf][1](i_row, idx_junction), zeta_col_in[i_surf][2](i_row, idx_junction);
        junction_point_1 << zeta_col_in[i_surf_partner_junction][0](i_row, idx_junction), zeta_col_in[i_surf_partner_junction][1](i_row, idx_junction), zeta_col_in[i_surf_partner_junction][2](i_row, idx_junction);
        
        for(uint i_col=0;i_col<N_col;++i_col)
        {
            for (uint idx_surf=0; idx_surf<2; idx_surf++)
            {
                // extra loop to avoid duplicate code for i_surf and i_surf_partner junction interpolation
                if (idx_surf==0)
                {
                    updated_surf = i_surf;
                    partner_surf = i_surf_partner_junction;
                }
                else
                {
                    updated_surf = i_surf_partner_junction;
                    partner_surf = i_surf;
                }
                collocation_point_phantom << zeta_col_out[updated_surf][0](i_row, i_col), zeta_col_out[updated_surf][1](i_row, i_col), zeta_col_out[updated_surf][2](i_row, i_col);                  
                distance_0 = (collocation_point_phantom - junction_point_0).norm();
                distance_1 = (collocation_point_phantom - junction_point_1).norm();
                
                gamma_out[updated_surf](i_row, i_col) = (abs(gamma_in[i_surf](i_row,idx_junction)) * distance_0 +abs(gamma_in[i_surf_partner_junction](i_row,idx_junction)) * distance_1)
                                                /(distance_0 + distance_1);
                                                
                if (gamma_in[updated_surf](i_row, idx_junction) < 0.0)
                {
                    gamma_out[updated_surf](i_row, i_col) *= -1;
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
    for(uint i_surf=0; i_surf<flag_zeta_phantom.cols(); ++i_surf)
    {
        for(uint i_row=0; i_row < flag_zeta_phantom.rows();++i_row)
        {
            if (flag_zeta_phantom(i_row, i_surf))
            {
                return true;
            }
        }
    }
    return false;
}