#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "biotsavart.h"
#include "unitsourcedensity.h"

#include <fstream>

namespace UVLM
{
    namespace Matrix
    {
        // DECLARATIONS
        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_zeta_star,
                //   typename t_zeta_star_col,
                  typename t_uext_col,
                  typename t_normals,
                  typename t_aic>
        void AIC
        (
            const uint& Ktotal,
            const t_zeta& zeta,
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            // const t_zeta_star_col& zeta_star_col,
            const t_uext_col& uext_col,
            const t_normals& normals,
            const UVLM::Types::VMopts& options,
            const bool horseshoe,
            t_aic& aic
        );

        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_surf_vec_panel,
                  typename t_surf_vec_col,
                  typename t_aic>
        void AIC_sources
        (
            const t_zeta& zeta,
            const t_zeta_col& zeta_col,
            const t_surf_vec_panel& longitudinals_panel,
            const t_surf_vec_panel& perpendiculars_panel,
            const t_surf_vec_panel& normals_panel,
            const t_surf_vec_col& longitudinals_col,
            const t_surf_vec_col& perpendiculars_col,
            const t_surf_vec_col& normals_col,
            const UVLM::Types::VMopts& options,
            t_aic& aic_sources_x,
            t_aic& aic_sources_y,
            t_aic& aic_sources_z,
            const bool& same_body = true
        );
        template <typename t_zeta_lifting,
                typename t_zeta_col_lifting,
                typename t_zeta_star_lifting,
                typename t_uext_col_lifting,
                typename t_surface_vec_lifting,
                typename t_zeta_nonlifting,
                typename t_zeta_col_nonlifting,
                typename t_uext_col_nonlifting,
                typename t_surface_vec_nonlifting,
                typename t_zeta_phantom,
                typename t_zeta_star_phantom,
                typename t_surface_vec_phantom,
                typename t_flag_zeta_phantom>
        void aic_combined
            (
                const uint& Ktotal,
                const uint& Ktotal_lifting,
                const uint& Ktotal_nonlifting, 
                const uint& Ktotal_phantom, 
                const UVLM::Types::VMopts& options,
                const UVLM::Types::VMopts& options_nonlifting,
                const t_zeta_lifting& zeta,
                const t_zeta_col_lifting& zeta_col,
                const t_zeta_star_lifting& zeta_star,
                const t_uext_col_lifting& uext_col,
                const t_surface_vec_lifting& normals,
                const t_surface_vec_lifting& longitudinals,
                const t_surface_vec_lifting& perpendiculars,
                UVLM::Types::MatrixX& aic_lifting,
                const t_zeta_nonlifting& zeta_nonlifting,
                const t_zeta_col_nonlifting& zeta_col_nonlifting,
                const t_uext_col_nonlifting& uext_col_nonlifting,
                const t_surface_vec_nonlifting& normals_nonlifting,
                const t_surface_vec_nonlifting& longitudinals_nonlifting,
                const t_surface_vec_nonlifting& perpendiculars_nonlifting,
                const t_zeta_phantom& zeta_phantom,
                const t_zeta_star_phantom& zeta_star_phantom,
                const t_surface_vec_phantom& normals_phantom,
                const t_surface_vec_phantom& longitudinals_phantom,
                const t_surface_vec_phantom& perpendiculars_phantom,
                const t_flag_zeta_phantom& flag_zeta_phantom,
                UVLM::Types::MatrixX& aic_nonlifting_x,
                UVLM::Types::MatrixX& aic_nonlifting_y,      
                UVLM::Types::MatrixX& aic_nonlifting_z,
                UVLM::Types::MatrixX& aic
            );
            template<typename t_zeta_col, 
                    typename t_zeta_phantom_col,
                    typename t_aic_phantom_out>
            void aic_phantom_interp_condition
            (
                const uint& Ktotal_lifting,
                const uint& Ktotal_nonlifting, 
                const uint& Ktotal_phantom, 
                t_zeta_col& zeta_col,
                t_zeta_phantom_col& zeta_phantom_col,
                t_aic_phantom_out& aic_phantom_out
            );
        template<typename t_aic_phantom,
                typename t_aic,
                typename t_zeta,
                typename t_flag_zeta_phantom>
        void add_phantom_AIC_to_AIC
        (
            t_aic_phantom& aic_phantom,
            t_aic& aic,
            t_zeta& zeta,
            t_flag_zeta_phantom& flag_zeta_phantom
        );
        template <typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_uext_col,
                //   typename t_zeta_dot_col,
                  typename t_gamma_star,
                  typename t_normal>
        void RHS
        (
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_uext_col& uext_col,
            // const t_zeta_dot_col& zeta_dot_col,
            const t_gamma_star& gamma_star,
            const t_normal& normal,
            const UVLM::Types::VMopts& options,
            UVLM::Types::VectorX& rhs,
            const uint& Ktotal
        );


        template <typename t_uext_col,
                  typename t_normal>
        void RHS_nonlifting_body
        (
            const t_uext_col& uinc_col,
            const t_normal& normal,
            UVLM::Types::VectorX& rhs,
            const uint& Ktotal,
            const uint n_surf
        );


        void generate_assembly_offset
        (
            const UVLM::Types::VecDimensions& dimensions,
            const UVLM::Types::VecDimensions& dimensions_star,
            UVLM::Types::MatrixX& offset,
            const bool steady
        );


        template <typename t_gamma,
                  typename t_zeta_col>
        void reconstruct_gamma
        (
            const UVLM::Types::VectorX& gamma_flat,
            t_gamma& gamma,
            const t_zeta_col& zeta_col
        );


        template <typename t_gamma,
                  typename t_zeta_col>
        void deconstruct_gamma
        (
            const t_gamma& gamma,
            UVLM::Types::VectorX& gamma_flat,
            const t_zeta_col& zeta_col
        );
		
		template <typename t_vec_mat,
		          typename t_zeta_col>
		void reconstruct_MatrixX
		(
			const UVLM::Types::MatrixX& mat_flat,
			t_vec_mat& vec_mat,
			const t_zeta_col& zeta_col
		);
		void build_offsets
		(
			const uint& n_surf,
			const UVLM::Types::VecDimensions& dimensions,
			std::vector<uint>& offset
		);
        
        uint get_total_VecVecMat_size(UVLM::Types::VecVecMatrixX mat_in);
    }
}
// SOURCE CODE
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/

template <typename t_zeta,
          typename t_zeta_col,
          typename t_zeta_star,
        //   typename t_zeta_star_col,
          typename t_uext_col,
          typename t_normals,
          typename t_aic>
void UVLM::Matrix::AIC
(
    const t_zeta& zeta,
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    // const t_zeta_star_col& zeta_star_col,
    const t_uext_col& uext_col,
    const t_normals& normals,
    const UVLM::Types::VMopts& options,
    const bool horseshoe,
    t_aic& aic
)
{
    const uint n_surf_panel = zeta.size();
    const uint n_surf_col = zeta_col.size();
    UVLM::Types::VecDimensions dimensions_panel, dimensions_col, dimensions_star;
    UVLM::Types::generate_dimensions(zeta, dimensions_panel, - 1);
    UVLM::Types::generate_dimensions(zeta_star, dimensions_star, -1);
    UVLM::Types::generate_dimensions(zeta_col, dimensions_col, 0);

    // build the offsets beforehand
    // (parallel variation)
    std::vector<uint> offset_panel;
    std::vector<uint> offset_col;
	UVLM::Matrix::build_offsets(n_surf_panel,
								dimensions_panel,
								offset_panel);
	UVLM::Matrix::build_offsets(n_surf_col,
								dimensions_col,
								offset_col);

    // fill up AIC
    for (uint icol_surf=0; icol_surf<n_surf_col; ++icol_surf)
    {
        uint k_surf_col = dimensions_col[icol_surf].first*
                          dimensions_col[icol_surf].second;

        // uint ii_offset = 0;
        for (uint jpanel_surf=0; jpanel_surf<n_surf_panel; ++jpanel_surf)
        {
            uint k_surf_panel = dimensions_panel[jpanel_surf].first*
								dimensions_panel[jpanel_surf].second;
            UVLM::Types::MatrixX dummy_gamma;
            UVLM::Types::MatrixX dummy_gamma_star;
            UVLM::Types::Block block = aic.block(offset_col[icol_surf], offset_panel[jpanel_surf], k_surf_col, k_surf_panel);
            // steady wake coefficients
            dummy_gamma.setOnes(dimensions_panel[jpanel_surf].first,
                                dimensions_panel[jpanel_surf].second);
            if (options.Steady)
            {
                dummy_gamma_star.setOnes(dimensions_star[jpanel_surf].first,
                                        dimensions_star[jpanel_surf].second);

                UVLM::BiotSavart::multisurface_steady_wake
                (
                    zeta[jpanel_surf],
                    zeta_star[jpanel_surf],
                    dummy_gamma,
                    dummy_gamma_star,
                    zeta_col[icol_surf],
                    horseshoe,
                    block,
                    options.ImageMethod,
                    normals[icol_surf],
                    options.vortex_radius
                );
            } else // unsteady case
            {
                dummy_gamma_star.setOnes(1,
                                         dimensions_star[jpanel_surf].second);
                UVLM::BiotSavart::multisurface_unsteady_wake
                (
                    zeta[jpanel_surf],
                    zeta_star[jpanel_surf],
                    dummy_gamma,
                    dummy_gamma_star,
                    zeta_col[icol_surf],
                    block,
                    options.ImageMethod,
                    normals[icol_surf],
                    0,
                    options.vortex_radius
                );
            }
        }
    }
}
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/

void UVLM::Matrix::build_offsets
(
    const uint& n_surf,
    const UVLM::Types::VecDimensions& dimensions,
    std::vector<uint>& offset
)
{
    uint i_offset = 0;
    uint k_surf;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        offset.push_back(i_offset);
        k_surf = dimensions[i_surf].first*
                 dimensions[i_surf].second;
       i_offset += k_surf;
    }
}
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/

template <typename t_zeta,
          typename t_zeta_col,
          typename t_surf_vec_panel,
          typename t_surf_vec_col,
          typename t_aic>

void UVLM::Matrix::AIC_sources
(
    const t_zeta& zeta,
    const t_zeta_col& zeta_col,
    const t_surf_vec_panel& longitudinals_panel,
    const t_surf_vec_panel& perpendiculars_panel,
    const t_surf_vec_panel& normals_panel,
    const t_surf_vec_col& longitudinals_col,
    const t_surf_vec_col& perpendiculars_col,
    const t_surf_vec_col& normals_col,
    const UVLM::Types::VMopts& options,
    t_aic& aic_sources_x,
    t_aic& aic_sources_y,
    t_aic& aic_sources_z,
    const bool& same_body
)
{
    const uint n_surf_panel = zeta.size();
    const uint n_surf_col = zeta_col.size();
    bool same_surface= false;
    UVLM::Types::VecDimensions dimensions_panel, dimensions_col;
    UVLM::Types::generate_dimensions(zeta, dimensions_panel, - 1);
    UVLM::Types::generate_dimensions(zeta_col, dimensions_col, 0);


    // build the offsets beforehand
    // (parallel variation)
    std::vector<uint> offset_panel;
    std::vector<uint> offset_col;
	UVLM::Matrix::build_offsets(n_surf_panel,
								dimensions_panel,
								offset_panel);
	UVLM::Matrix::build_offsets(n_surf_col,
								dimensions_col,
								offset_col);

    // fill up AIC
    for (uint icol_surf=0; icol_surf<n_surf_col; ++icol_surf)
    {
        uint k_surf_col_i = dimensions_col[icol_surf].first*
                            dimensions_col[icol_surf].second;

        // uint ii_offset = 0;
        for (uint jpanel_surf=0; jpanel_surf<n_surf_panel; ++jpanel_surf)
        {
            uint k_surf_panel_j = dimensions_panel[jpanel_surf].first*
                                  dimensions_panel[jpanel_surf].second;
            UVLM::Types::Block block_aic_x = aic_sources_x.block(offset_col[icol_surf], offset_panel[jpanel_surf], k_surf_col_i, k_surf_panel_j);
            UVLM::Types::Block block_aic_y = aic_sources_y.block(offset_col[icol_surf], offset_panel[jpanel_surf], k_surf_col_i, k_surf_panel_j);
            UVLM::Types::Block block_aic_z = aic_sources_z.block(offset_col[icol_surf], offset_panel[jpanel_surf], k_surf_col_i, k_surf_panel_j);

            if (options.Steady)
            {
                same_surface=((icol_surf==jpanel_surf)&&(same_body));
                UVLM::UnitSourceDensity::get_influence_coefficient(zeta[jpanel_surf],
                                                                zeta_col[icol_surf],
                                                                block_aic_x,
                                                                block_aic_y,
                                                                block_aic_z,
                                                                longitudinals_panel[jpanel_surf],
                                                                perpendiculars_panel[jpanel_surf],
                                                                normals_panel[jpanel_surf],
                                                                longitudinals_col[icol_surf],
                                                                perpendiculars_col[icol_surf],
                                                                normals_col[icol_surf],
                                                                same_surface
                                                                );
            }
        }
    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_col,
          typename t_gamma_star,
          typename t_normal>
void UVLM::Matrix::RHS
(
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_col& uinc_col,
    const t_gamma_star& gamma_star,
    const t_normal& normal,
    const UVLM::Types::VMopts& options,
    UVLM::Types::VectorX& rhs,
    const uint& Ktotal
)
{
    const uint n_surf = options.NumSurfaces;

    rhs.setZero(Ktotal);

    // filling up RHS
    int ii = -1;
    int istart = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = uinc_col[i_surf][0].rows();
        uint N = uinc_col[i_surf][0].cols();

        if (!options.Steady)
        {
            #pragma omp parallel for collapse(2)
            for (uint i=0; i<M; ++i)
            {
                for (uint j=0; j<N; ++j)
                {
                    UVLM::Types::Vector3 v_ind;
                    UVLM::Types::Vector3 collocation_coords;
                    UVLM::Types::Vector3 u_col;

                    u_col << uinc_col[i_surf][0](i,j),
                             uinc_col[i_surf][1](i,j),
                             uinc_col[i_surf][2](i,j);
                    // we have to add the wake effect on the induced velocity.
                    collocation_coords << zeta_col[i_surf][0](i,j),
                                          zeta_col[i_surf][1](i,j),
                                          zeta_col[i_surf][2](i,j);

                    v_ind.setZero();
                    for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                    {
                        v_ind += UVLM::BiotSavart::whole_surface(zeta_star[ii_surf],
                                                                      gamma_star[ii_surf],
                                                                      collocation_coords,
                                                                      0,
                                                                      0,
                                                                      -1,
                                                                      -1,
                                                                      options.ImageMethod,
                                                                      options.vortex_radius);
                    }
                    u_col += v_ind;

                    // dot product of uinc and panel normal
                    uint counter = istart + j + i*N;
                    rhs(counter) =
                    -(
                        u_col(0)*normal[i_surf][0](i,j) +
                        u_col(1)*normal[i_surf][1](i,j) +
                        u_col(2)*normal[i_surf][2](i,j)
                    );
                }
            }
            istart += M*N;
        } else {
            for (uint i=0; i<M; ++i)
            {
                for (uint j=0; j<N; ++j)
                {
                    rhs(++ii) =
                    -(
                        uinc_col[i_surf][0](i,j)*normal[i_surf][0](i,j) +
                        uinc_col[i_surf][1](i,j)*normal[i_surf][1](i,j) +
                        uinc_col[i_surf][2](i,j)*normal[i_surf][2](i,j)
                    );
                }
            }
        }
    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_uext_col,
          typename t_normal>
void UVLM::Matrix:: RHS_nonlifting_body
(
    const t_uext_col& uinc_col,
    const t_normal& normal,
    UVLM::Types::VectorX& rhs,
    const uint& Ktotal,
    const uint n_surf
)
{
    rhs.setZero(Ktotal);

    // filling up RHS
    int ii = -1;
    int istart = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = uinc_col[i_surf][0].rows();
        uint N = uinc_col[i_surf][0].cols();
        for (uint i=0; i<M; ++i)
        {
            for (uint j=0; j<N; ++j)
            {
                rhs(++ii) =
                -(
                    uinc_col[i_surf][0](i,j)*normal[i_surf][0](i,j) +
                    uinc_col[i_surf][1](i,j)*normal[i_surf][1](i,j) +
                    uinc_col[i_surf][2](i,j)*normal[i_surf][2](i,j)
                );
            }
        }
    }
}


/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/

void UVLM::Matrix::generate_assembly_offset
(
    const UVLM::Types::VecDimensions& dimensions,
    const UVLM::Types::VecDimensions& dimensions_star,
    UVLM::Types::MatrixX& offset,
    const bool steady
)
{
    std::cerr << "Not sure this is correct, dont use (or debug)!!!" << std::endl;
    uint n_surf = dimensions.size();
    offset.setZero(n_surf, 2);

    uint counter = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        offset(i_surf, 0) = counter;
        counter += dimensions[i_surf].first*
                   dimensions[i_surf].second;

        offset(i_surf, 1) = counter;
        if (!steady)
        {
            counter += dimensions_star[i_surf].first*
                       dimensions_star[i_surf].second;
        }
    }
}

template <typename t_gamma,
          typename t_zeta_col>
void UVLM::Matrix::reconstruct_gamma
(
    const UVLM::Types::VectorX& gamma_flat,
    t_gamma& gamma,
    const t_zeta_col& zeta_col
)
{
    const uint n_surf = gamma.size();
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);

    uint i_flat = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        for (uint i=0; i<dimensions[i_surf].first; ++i)
        {
            for (uint j=0; j<dimensions[i_surf].second; ++j)
            {
                gamma[i_surf](i, j) = gamma_flat(i_flat++);
            }
        }
    }

}
template <typename t_vec_mat,
		  typename t_zeta_col>
void UVLM::Matrix::reconstruct_MatrixX
(
    const UVLM::Types::MatrixX& mat_flat,
    t_vec_mat& vec_mat,
    const t_zeta_col& zeta_col
)
{
    const uint n_surf = zeta_col.size();
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);

    uint i_flat = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        for (uint i=0; i<dimensions[i_surf].first; ++i)
        {
            for (uint j=0; j<dimensions[i_surf].second; ++j)
            {
				for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; i_dim++)
				{
					vec_mat[i_surf][i_dim](i, j) = mat_flat(i_dim,i_flat);
				}
				i_flat++;
            }
        }
    }

}
template <typename t_gamma,
          typename t_zeta_col>
void UVLM::Matrix::deconstruct_gamma
(
    const t_gamma& gamma,
    UVLM::Types::VectorX& gamma_flat,
    const t_zeta_col& zeta_col
)
{
    const uint n_surf = gamma.size();
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);

    uint n_total = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        n_total += dimensions[i_surf].first*dimensions[i_surf].second;
    }
    gamma_flat.resize(n_total);

    uint i_flat = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        for (uint i=0; i<dimensions[i_surf].first; ++i)
        {
            for (uint j=0; j<dimensions[i_surf].second; ++j)
            {
                gamma_flat(i_flat++) = gamma[i_surf](i, j);
            }
        }
    }

}

template <typename t_zeta_lifting,
        typename t_zeta_col_lifting,
        typename t_zeta_star_lifting,
        typename t_uext_col_lifting,
        typename t_surface_vec_lifting,
        typename t_zeta_nonlifting,
        typename t_zeta_col_nonlifting,
        typename t_uext_col_nonlifting,
        typename t_surface_vec_nonlifting,
        typename t_zeta_phantom,
        typename t_zeta_star_phantom,
        typename t_surface_vec_phantom,
        typename t_flag_zeta_phantom>
void UVLM::Matrix::aic_combined
    (
        const uint& Ktotal,
        const uint& Ktotal_lifting,
        const uint& Ktotal_nonlifting, 
        const uint& Ktotal_phantom, 
        const UVLM::Types::VMopts& options,
        const UVLM::Types::VMopts& options_nonlifting,
        const t_zeta_lifting& zeta,
        const t_zeta_col_lifting& zeta_col,
        const t_zeta_star_lifting& zeta_star,
        const t_uext_col_lifting& uext_col,
        const t_surface_vec_lifting& normals,
        const t_surface_vec_lifting& longitudinals,
        const t_surface_vec_lifting& perpendiculars,
        UVLM::Types::MatrixX& aic_lifting,
        const t_zeta_nonlifting& zeta_nonlifting,
        const t_zeta_col_nonlifting& zeta_col_nonlifting,
        const t_uext_col_nonlifting& uext_col_nonlifting,
        const t_surface_vec_nonlifting& normals_nonlifting,
        const t_surface_vec_nonlifting& longitudinals_nonlifting,
        const t_surface_vec_nonlifting& perpendiculars_nonlifting,
        const t_zeta_phantom& zeta_phantom,
        const t_zeta_star_phantom& zeta_star_phantom,
        const t_surface_vec_phantom& normals_phantom,
        const t_surface_vec_phantom& longitudinals_phantom,
        const t_surface_vec_phantom& perpendiculars_phantom,
        const t_flag_zeta_phantom& flag_zeta_phantom,
        UVLM::Types::MatrixX& aic_nonlifting_x,
        UVLM::Types::MatrixX& aic_nonlifting_y,       
        UVLM::Types::MatrixX& aic_nonlifting_z,
        UVLM::Types::MatrixX& aic
    )
    {
        
	//TO-DO: delete unnessary inputs like aic_lifting
    //UVLM::Types::MatrixX aic_lifting = UVLM::Types::MatrixX::Zero(Ktotal_lifting, Ktotal_lifting);
    UVLM::Matrix::AIC(zeta,
                      zeta_col,
                      zeta_star,
                      uext_col,
                      normals,
                      options,
                      false,
                      aic_lifting);
        
    // LNonfting on nonlifting surfaces
    UVLM::Matrix::AIC_sources(zeta_nonlifting,
                            zeta_col_nonlifting,
                            longitudinals_nonlifting,
                              perpendiculars_nonlifting,
                              normals_nonlifting,
                              longitudinals_nonlifting, //collocation
                              perpendiculars_nonlifting, //collocation
                              normals_nonlifting, //collocation
                              options_nonlifting,
							  aic_nonlifting_x,
							  aic_nonlifting_y,
                              aic_nonlifting_z,
                              true);

    // Lifting on nonlifting surfaces
	UVLM::Types::MatrixX aic_lifting_on_nonlifting = UVLM::Types::MatrixX::Zero(Ktotal_nonlifting, Ktotal_lifting);
    UVLM::Matrix::AIC(zeta,
                      zeta_col_nonlifting,
                      zeta_star,
                      uext_col_nonlifting,
                      normals_nonlifting,
                      options,
                      false,
                      aic_lifting_on_nonlifting);

    // Nonlifting on lifting surfaces
	UVLM::Types::MatrixX aic_nonlifting_on_lifting_z = UVLM::Types::MatrixX::Zero(Ktotal_lifting, Ktotal_nonlifting);
    UVLM::Types::MatrixX aic_nonlifting_on_lifting_x = UVLM::Types::MatrixX::Zero(Ktotal_lifting, Ktotal_nonlifting);
    UVLM::Types::MatrixX aic_nonlifting_on_lifting_y = UVLM::Types::MatrixX::Zero(Ktotal_lifting, Ktotal_nonlifting);
    UVLM::Matrix::AIC_sources(zeta_nonlifting,
                              zeta_col,
                              longitudinals_nonlifting,
                              perpendiculars_nonlifting,
                              normals_nonlifting,
                              longitudinals, //collocation
                              perpendiculars, //collocation
                              normals, //collocation
                              options_nonlifting,
							  aic_nonlifting_on_lifting_x,
							  aic_nonlifting_on_lifting_y,
                              aic_nonlifting_on_lifting_z,
                              false);    
                 
    //Phantom panels        
    UVLM::Types::MatrixX aic_phantom_on_lifting = UVLM::Types::MatrixX::Zero(Ktotal_lifting, Ktotal_phantom);    
    UVLM::Matrix::AIC(zeta_phantom,
                      zeta_col,
                      zeta_star_phantom,
                      uext_col,
                      normals,
                      options,
                      false,
                      aic_phantom_on_lifting);

    UVLM::Types::MatrixX aic_phantom_on_nonlifting = UVLM::Types::MatrixX::Zero(Ktotal_nonlifting, Ktotal_phantom);   
    UVLM::Matrix::AIC(zeta_phantom,
                      zeta_col_nonlifting,
                      zeta_star_phantom,
                      uext_col_nonlifting,
                      normals_nonlifting,
                      options,
                      false,
                      aic_phantom_on_nonlifting);
    
    // Get matrix to enforce linear interpolated circulation on phantom panels
    UVLM::Types::VecVecMatrixX zeta_col_phantom;          
    UVLM::Geometry::generate_colocationMesh(zeta_phantom, zeta_col_phantom);     
    UVLM::Types::MatrixX circulation_bc_phantom = UVLM::Types::MatrixX::Zero(Ktotal_phantom, Ktotal+Ktotal_phantom);
    UVLM::Matrix::aic_phantom_interp_condition(Ktotal_lifting,
                                               Ktotal_nonlifting, 
                                               Ktotal_phantom, 
                                               zeta_col,
                                               zeta_col_phantom,
                                               circulation_bc_phantom);
    UVLM::Types::copy_Mat_to_block(aic_lifting, aic, 0, 0);   

    UVLM::Types::copy_Mat_to_block(aic_nonlifting_z, aic, Ktotal_lifting, Ktotal_lifting);
    UVLM::Types::copy_Mat_to_block(aic_lifting_on_nonlifting, aic, Ktotal_lifting,0);
    UVLM::Types::copy_Mat_to_block(aic_nonlifting_on_lifting_z, aic, 0, Ktotal_lifting);

    UVLM::Types::copy_Mat_to_block(aic_phantom_on_nonlifting, aic, Ktotal_lifting, Ktotal);

    UVLM::Types::copy_Mat_to_block(aic_phantom_on_lifting, aic, 0, Ktotal);
    UVLM::Types::copy_Mat_to_block(circulation_bc_phantom, aic, Ktotal, 0);
}

uint UVLM::Matrix::get_total_VecVecMat_size(UVLM::Types::VecVecMatrixX mat_in)
{
    uint n_surf = mat_in.size();
    uint ii = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = mat_in[i_surf][0].rows();
        uint N = mat_in[i_surf][0].cols();

        ii += M*N;
    }
    return ii;
}



template<typename t_aic_phantom,
            typename t_aic,
            typename t_zeta,
            typename t_flag_zeta_phantom>
void UVLM::Matrix::add_phantom_AIC_to_AIC
(
    t_aic_phantom& aic_phantom,
    t_aic& aic,
    t_zeta& zeta,
    t_flag_zeta_phantom& flag_zeta_phantom
)
{
    /* Function adds the AIC of the phantom panel on the collocation points to the adjacent wing panel
       --> Keeps the combined AIC matrix squared
       --> enforces the circulation BC at the junction (Gamma_phantom == Gamma_wing_junction)*/
    const uint n_surf_panel = zeta.size();
    const uint N_rows_aic = aic.rows();
    uint N_phantom_cell, N_row_phantom, index_phantom, offset_aic_phantom = 0;
    // build the offsets beforehand
    std::vector<uint> offset_panel;
    UVLM::Types::VecDimensions dimensions_panel;
    UVLM::Types::generate_dimensions(zeta, dimensions_panel, - 1);
	UVLM::Matrix::build_offsets(n_surf_panel,
								dimensions_panel,
								offset_panel);

    for (uint jpanel_surf=0; jpanel_surf<n_surf_panel; ++jpanel_surf)
    {
        // Get index
        for(Eigen::Index i=0; i<flag_zeta_phantom.rows(); ++i)
        {
            if(flag_zeta_phantom(i, jpanel_surf))
            {
                index_phantom = i;
                break;
            }
        }
        N_row_phantom = dimensions_panel[jpanel_surf].first;
        
        for(Eigen::Index j_panel=0; j_panel<flag_zeta_phantom.rows(); ++j_panel)
        {
            if(flag_zeta_phantom(j_panel, jpanel_surf))
            {
                for(uint i_phantom_row=0; i_phantom_row<N_row_phantom; ++i_phantom_row)
                {
                aic.col(offset_panel[jpanel_surf]+j_panel*N_row_phantom+i_phantom_row)
                        +=aic_phantom.col(offset_aic_phantom+i_phantom_row);
                }
                offset_aic_phantom+=N_row_phantom;
                }
            }
    }
}

template<typename t_zeta_col, 
         typename t_zeta_phantom_col,
         typename t_aic_phantom_out>
void UVLM::Matrix::aic_phantom_interp_condition
(
    const uint& Ktotal_lifting,
    const uint& Ktotal_nonlifting, 
    const uint& Ktotal_phantom, 
    t_zeta_col& zeta_col,
    t_zeta_phantom_col& zeta_phantom_col,
    t_aic_phantom_out& aic_phantom_out
)
{
    uint N_phantom_panels, M_phantom_panels, counter_row, N_spanwise_panels;
    double yL0_minus_yR0, interpolated_value;
    const uint n_surf = zeta_col.size();
    // set influence of phantom panel on own collocation point
    aic_phantom_out.topRightCorner(Ktotal_phantom, Ktotal_phantom).setIdentity();
    aic_phantom_out *= -1;

    // get offsets beforehand    
    std::vector<uint> offset_panel, offset_phantom_panel;
    UVLM::Types::VecDimensions dimensions_panel, dimensions_phantom_panel;
    UVLM::Types::generate_dimensions(zeta_col, dimensions_panel, 0);
    UVLM::Types::generate_dimensions(zeta_phantom_col, dimensions_phantom_panel, 0);
	UVLM::Matrix::build_offsets(n_surf, dimensions_panel,offset_panel);
	UVLM::Matrix::build_offsets(n_surf, dimensions_phantom_panel,offset_phantom_panel);

    for(uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        counter_row = 0;
        N_spanwise_panels = zeta_col[i_surf][0].cols();
        N_phantom_panels= zeta_phantom_col[i_surf][0].rows();
        M_phantom_panels= zeta_phantom_col[i_surf][0].cols();
        
        for (uint i_row=0; i_row<N_phantom_panels; ++i_row)
        {
            yL0_minus_yR0 = zeta_col[1][1](i_row, 0) - zeta_col[0][1](i_row, 0);
            for (uint i_col=0; i_col<M_phantom_panels; ++i_col)
            {
                interpolated_value = (zeta_phantom_col[i_surf][1](i_row, i_col)-zeta_col[0][1](i_row, 0))/yL0_minus_yR0;
                aic_phantom_out(offset_phantom_panel[i_surf]+counter_row,offset_panel[0]+i_row*N_spanwise_panels) = 1-interpolated_value;
                aic_phantom_out(offset_phantom_panel[i_surf]+counter_row,offset_panel[1]+i_row*N_spanwise_panels) = interpolated_value;
                if (i_surf==0)
                {
                    aic_phantom_out(offset_phantom_panel[i_surf]+counter_row,offset_panel[1]+i_row*N_spanwise_panels)*=-1;
                }
                else
                {
                    aic_phantom_out(offset_phantom_panel[i_surf]+counter_row,offset_panel[0]+i_row*N_spanwise_panels)*=-1;
                }
                counter_row++;
            }
        }
    }
}
