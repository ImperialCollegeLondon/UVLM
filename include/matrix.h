#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "biotsavart.h"
#include "sources.h"

#include <fstream>

/**
 * @brief Namespace for the UVLM library.
 */
namespace UVLM
{
    /**
     * @brief Namespace for functions involved to generate the linear set of equation to be solved-state.
     */
    namespace Matrix
    {
    /**
     * \brief Computes the Aerodynamic Influence Coefficients (AIC) matrix for vortex ring singularities.
     * 
     * This function calculates the AIC matrix that represents the influence of vortex filaments
     * on panel collocation points. The matrix relates the induced velocities at collocation points
     * to the strength of vortices associated with vortex filaments. The influence is computed using
     * the Biot-Savart law, considering both steady and unsteady flow conditions.
     * 
     * \tparam t_zeta Type representing the panel vortex corner point coordinates.
     * \tparam t_zeta_col Type representing the collocation points on panels.
     * \tparam t_zeta_star Type representing the vortex sheet vertices.
     * \tparam t_uext_total_col Type representing the total external velocities at collocation points.
     * \tparam t_normals Type representing the panel normals.
     * \tparam t_aic Type representing the computed AIC matrix.
     * \param zeta The panel vortex corner point coordinates.
     * \param zeta_col The collocation points on panels.
     * \param zeta_star The vortex sheet vertices.
     * \param uext_total_col The total external velocities at collocation points.
     * \param normals The panel normals.
     * \param options Struct containing simulation options.
     * \param horseshoe Flag indicating whether horseshoe vortices are used.
     * \param aic Computed AIC matrix (output).
    **/
        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_uext_total_col,
                  typename t_normals,
                  typename t_aic>
        void AIC
        (
            const t_zeta& zeta,
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_uext_total_col& uext_total_col,
            const t_normals& normals,
            const UVLM::Types::VMopts& options,
            const bool horseshoe,
            t_aic& aic
        );

        /**
         * @brief Calculate the AIC matrices for source panel singularities.
         *
         * This function computes the AIC matrices induced by a linear source panel on a set of collocation points
         * of an arbitary surface (lifting, nonlifting, or panel).
         *
         * @tparam t_zeta Type representing the source panel corner point coordinates.
         * @tparam t_zeta_col Type of the input matrix for zeta_col (collocation points).
         * @tparam t_surf_vec_panel Type of the input matrix for surface vectors on panels.
         * @tparam t_surf_vec_col Type of the input matrix for surface vectors at collocation points.
         * @tparam t_aic Type representing the computed AIC matrix.
         * @param zeta Matrix of coordinates of the non-lifting body.
         * @param zeta_col Matrix of collocation points on the non-lifting body.
         * @param longitudinals_panel Surface vectors in the longitudinal direction on panels.
         * @param perpendiculars_panel Surface vectors perpendicular to the panels.
         * @param normals_panel Surface normals on panels.
         * @param longitudinals_col Surface vectors in the longitudinal direction at collocation points.
         * @param perpendiculars_col Surface vectors perpendicular to the collocation points.
         * @param normals_col Surface normals at collocation points.
         * @param aic_sources_x Output AIC matrix for x-direction.
         * @param aic_sources_y Output AIC matrix for y-direction.
         * @param aic_sources_z Output AIC matrix for z-direction.
         * @param same_body Flag indicating whether the sources are on the same body.
         */
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
            t_aic& aic_sources_x,
            t_aic& aic_sources_y,
            t_aic& aic_sources_z,
            const bool& same_body = true
        );

        /**
         * @brief Computes the AIC matrix in case multiple type of singularities (i.e. vortex rings and source panels) exist.
         * 
         * The AIC matrices are computed separately for each surface and singularity type, and are merged finally in a combined AIC matrix.
         *
         * @tparam t_struct_lifting_surfaces Type of the struct representing lifting surfaces.
         * @tparam t_struct_nl_body Type of the struct representing non-lifting bodies.
         * @tparam t_struct_phantom_surf Type of the struct representing phantom surfaces.
         *
         * @param lifting_surfaces Struct representing lifting surfaces.
         * @param nl_body Struct representing non-lifting bodies.
         * @param phantom_surfaces Struct representing phantom surfaces.
         * @param aic Output combined AIC matrix.
         * @param options Options for the UVLM solver.
         */
        template <typename t_struct_lifting_surfaces,
                typename t_struct_nl_body,
                typename t_struct_phantom_surf>
        void aic_combined
            (
                t_struct_lifting_surfaces& lifting_surfaces,
                t_struct_nl_body& nl_body,
                t_struct_phantom_surf& phantom_surfaces,
                UVLM::Types::MatrixX& aic,
                const UVLM::Types::VMopts& options
            );
        /**
         * @brief Computes the values within the AIC matrix to set the boundary conditions for the circulation strength of the phantom vortex panels. 
         *
         * @tparam t_zeta_col Type of the zeta_col (collocation points) data structure.
         * @tparam t_zeta_phantom_col Type of the zeta_phantom_col (collocation points of phantom surfaces) data structure.
         * @tparam t_aic_phantom_out Type of the output phantom AIC data structure.
         * @tparam t_flag_zeta_phantom Type of the flag_zeta_phantom data structure.
         *
         * @param Ktotal_lifting Total number of collocation points for lifting surfaces.
         * @param Ktotal_phantom Total number of collocation points for phantom surfaces.
         * @param zeta_col Collocation points located on the lifting surfaces.
         * @param zeta_phantom_col Collocation points located on the phantom surfaces.
         * @param aic_phantom_out Output phantom AIC matrix.
         * @param flag_zeta_phantom Flag indicating phantom AIC conditions.
         * @param only_for_update Flag indicating whether to update phantom AIC.
         */
        template<typename t_zeta_col, 
                typename t_zeta_phantom_col,
                typename t_aic_phantom_out,
                typename t_flag_zeta_phantom>
        void aic_phantom_interp_condition
        (
            const uint& Ktotal_lifting,
            const uint& Ktotal_phantom, 
            t_zeta_col& zeta_col,
            t_zeta_phantom_col& zeta_phantom_col,
            t_aic_phantom_out& aic_phantom_out,
            t_flag_zeta_phantom& flag_zeta_phantom,
        const bool only_for_update = false
    );
        /**
         * @brief Calculates the right-hand side (RHS)/ non-penetrating boundary conditions for the given set of collocation points.
         *
         * @tparam t_zeta_col Type of the zeta_col (surface panel collocation points) data structure.
         * @tparam t_zeta_star Type of the zeta_star (wake surface corner points) data structure.
         * @tparam t_uext_total_col Type of the uext_total_col (total external velocities at collocation points) data structure.
         * @tparam t_gamma_star Type of the gamma_star (circulation strength on wake panels) data structure.
         * @tparam t_normal Type of the normal vectors data structure.
         *
         * @param zeta_col Collocation points of the lifting surface.
         * @param zeta_star Vortex ring panel vorner points for wake surfaces.
         * @param uext_total_col Total external velocities at collocation points.
         * @param gamma_star Circulation strength on the vortex ring wake panels of the lifting surfaces.
         * @param normal Normal vectors of the panel associated with the given collocation points.
         * @param options Options for the UVLM solver.
         * @param rhs Output boundary condition vector.
         * @param Ktotal Total number of collocation points.
         */
        template <typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_uext_total_col,
                  typename t_gamma_star,
                  typename t_normal>
        void RHS
        (
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_uext_total_col& uext_total_col,
            const t_gamma_star& gamma_star,
            const t_normal& normal,
            const UVLM::Types::VMopts& options,
            UVLM::Types::VectorX& rhs,
            const uint& Ktotal
        );
        /**
         * @brief Calculates the non-penetrating boundary conditions for the given set of collocation points for an unsteady case including phantom panels.
         *
         * @tparam t_zeta_col Type of the zeta_col (surface panel collocation points) data structure.
         * @tparam t_zeta_star Type of the zeta_star (wake surface corner points) data structure.
         * @tparam t_uext_total_col Type of the uext_total_col (total external velocities at collocation points) data structure.
         * @tparam t_gamma_star Type of the gamma_star (circulation strength on wake panels) data structure.
         * @tparam t_normal Type of the normal vectors data structure.
         *
         * @param zeta_col Collocation points of the lifting surface.
         * @param zeta_star Vortex ring panel corner points for lifting wake surfaces.
         * @param uext_total_col Total external velocities (flow, structure, and rigid body motions) at collocation points.
         * @param gamma_star Circulation strength on the vortex ring wake panels of the lifting surfaces.
         * @param normal Normal vectors of the panel associated with the given collocation points.
         * @param options Options for the UVLM solver.
         * @param rhs Output boundary condition vector.
         * @param Ktotal Total number of collocation points.
         * @param phantom_gamma_star Circulation strength on the vortex ring wake panels of the phantom surfaces.
         * @param phantom_zeta_star Vortex ring panel corner points for phantom wake surfaces.
         * */
        template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_col,
          typename t_gamma_star,
          typename t_normal,
          typename t_phantom_gamma_star,
          typename t_phantom_zeta_star>
        void RHS_unsteady_phantom_unsteady
        (
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_uext_col& uinc_col,
            const t_gamma_star& gamma_star,
            const t_normal& normal,
            const UVLM::Types::VMopts& options,
            UVLM::Types::VectorX& rhs,
            const uint& Ktotal,
            const t_phantom_gamma_star& phantom_gamma_star,
            const t_phantom_zeta_star& phantom_zeta_star
        );

        /**
         * @brief Calculates the non-penetrating boundary conditions for the given set of collocation points in the presence of source panels.
         *
         * @tparam t_uext_col Type of the uext_col (external velocities at collocation points) data structure.
         * @tparam t_normal Type of the normal vectors data structure.
         *
         * @param uinc_col External velocities at collocation points.
         * @param normal Normal vectors on the panels associated with the collocation points.
         * @param rhs Output boundary condition vector.
         * @param Ktotal Total number of collocation points.
         * @param n_surf Number of considered surfaces.
         */
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

        /**
         * @brief Copies values from vector to a VecMatrix with the same size but different shape.
         *
         * @tparam t_gamma Type of the gamma (circulation) data structure.
         * @tparam t_zeta_col Type of the zeta_col (collocation points) data structure.
         *
         * @param gamma_flat Gamma vector.
         * @param gamma Gamma matrix.
         * @param zeta_col Collocation points of the surfaces.
         */
        template <typename t_gamma,
                  typename t_zeta_col>
        void reconstruct_VecMatrixX_values_from_vector
        (
            const UVLM::Types::VectorX& gamma_flat,
            t_gamma& gamma,
            const t_zeta_col& zeta_col
        );

        /**
         * @brief Copies values from VecMatrixX to a vector with the same size but different shape.
         *
         * @tparam t_gamma Type of the gamma (circulation) data structure.
         * @tparam t_zeta_col Type of the zeta_col (collocation points) data structure.
         *
         * @param gamma Gamma matrix.
         * @param gamma_flat Gamma vector.
         * @param zeta_col Collocation points of the surfaces.
         */
        template <typename t_gamma,
                  typename t_zeta_col>
        void reconstruct_vector_values_from_VecMatrixX
        (
            const t_gamma& gamma,
            UVLM::Types::VectorX& gamma_flat,
            const t_zeta_col& zeta_col
        );
		
        /**
         * @brief Copies values from a MatriX to a VecVecMatrix with the same size but different shape.
         *
         * @tparam t_vec_mat ype for VecVecMatrixX (sometimes mapped matrix).
         * @tparam t_zeta_col Type of the zeta_col (collocation points) data structure.
         *
         * @param vec_mat Output matrix with values to be filled in this function.
         * @param zeta_col Collocation points on the surfaces (for dimensions).
         */
        template <typename t_vec_mat,
		          typename t_zeta_col>
		void reconstruct_VecVecMatrixX_values_from_MatrixX
		(
			const UVLM::Types::MatrixX& mat_flat,
			t_vec_mat& vec_mat,
			const t_zeta_col& zeta_col
		);
        /**
         * @brief Builds assembly offsets for UVLM matrices.
         *
         * @param n_surf Number of surfaces.
         * @param dimensions Dimensions of the surfaces.
         * @param offset Output assembly offset vector.
         */
		void build_offsets
		(
			const uint& n_surf,
			const UVLM::Types::VecDimensions& dimensions,
			std::vector<uint>& offset
		);
        
        /**
         * @brief Calculates the total size of a VecVecMat data structure.
         *
         * @tparam t_mat_in Type of the input matrix data structure.
         *
         * @param mat_in Input matrix data structure.
         *
         * @return Total size of the VecVecMat data structure.
         */
        template<typename t_mat_in>
        uint get_total_VecVecMat_size(t_mat_in mat_in);
    }
}
// SOURCE CODE
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/

template <typename t_zeta,
          typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_total_col,
          typename t_normals,
          typename t_aic>
void UVLM::Matrix::AIC
(
    const t_zeta& zeta,
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_total_col& uext_total_col,
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

            same_surface=((icol_surf==jpanel_surf)&&(same_body));
            UVLM::Sources::get_influence_coefficient(zeta[jpanel_surf],
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

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_total_col,
          typename t_gamma_star,
          typename t_normal>
void UVLM::Matrix::RHS
(
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_total_col& uext_total_col,
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
        uint M = uext_total_col[i_surf][0].rows();
        uint N = uext_total_col[i_surf][0].cols();
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

                    u_col << uext_total_col[i_surf][0](i,j),
                             uext_total_col[i_surf][1](i,j),
                             uext_total_col[i_surf][2](i,j);
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
                        uext_total_col[i_surf][0](i,j)*normal[i_surf][0](i,j) +
                        uext_total_col[i_surf][1](i,j)*normal[i_surf][1](i,j) +
                        uext_total_col[i_surf][2](i,j)*normal[i_surf][2](i,j)
                    );
                }
            }
        }
    }
}
template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_col,
          typename t_gamma_star,
          typename t_normal,
          typename t_phantom_gamma_star,
          typename t_phantom_zeta_star>
void UVLM::Matrix::RHS_unsteady_phantom_unsteady
(
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_col& uinc_col,
    const t_gamma_star& gamma_star,
    const t_normal& normal,
    const UVLM::Types::VMopts& options,
    UVLM::Types::VectorX& rhs,
    const uint& Ktotal,
    const t_phantom_gamma_star& phantom_gamma_star,
    const t_phantom_zeta_star& phantom_zeta_star
)
{
    // TODO: Merge with RHS Function above
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
                                                                    options.ImageMethod,
                                                                    options.vortex_radius);
                                                                    
                        v_ind += UVLM::BiotSavart::whole_surface(phantom_zeta_star[ii_surf],
                                                                    phantom_gamma_star[ii_surf],
                                                                    collocation_coords,
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
            istart += M*N; //???
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
template <typename t_gamma,
          typename t_zeta_col>
void UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector
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
void UVLM::Matrix::reconstruct_VecVecMatrixX_values_from_MatrixX
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
void UVLM::Matrix::reconstruct_vector_values_from_VecMatrixX
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

template <typename t_struct_lifting_surfaces,
        typename t_struct_nl_body,
        typename t_struct_phantom_surf>
void UVLM::Matrix::aic_combined
    (
        t_struct_lifting_surfaces& lifting_surfaces,
        t_struct_nl_body& nl_body,
        t_struct_phantom_surf& phantom_surfaces,
        UVLM::Types::MatrixX& aic,
        const UVLM::Types::VMopts& options
    )
    {
    const uint Ktotal_lifting_and_nonlifting = lifting_surfaces.Ktotal + nl_body.Ktotal;
    if (!options.phantom_wing_test)
    {
    // Lifting on nonlifting surfaces
	UVLM::Types::MatrixX aic_lifting_on_nonlifting = UVLM::Types::MatrixX::Zero(nl_body.Ktotal, lifting_surfaces.Ktotal);
    UVLM::Matrix::AIC(lifting_surfaces.zeta,
                      nl_body.zeta_col,
                      lifting_surfaces.zeta_star,
                      nl_body.uext_col,
                      nl_body.normals,
                      options,
                      options.horseshoe,
                      aic_lifting_on_nonlifting);
    // Nonlifting on lifting surfaces
	UVLM::Types::MatrixX aic_nonlifting_on_lifting_z = UVLM::Types::MatrixX::Zero(lifting_surfaces.Ktotal, nl_body.Ktotal);
    UVLM::Types::MatrixX aic_nonlifting_on_lifting_x = UVLM::Types::MatrixX::Zero(lifting_surfaces.Ktotal, nl_body.Ktotal);
    UVLM::Types::MatrixX aic_nonlifting_on_lifting_y = UVLM::Types::MatrixX::Zero(lifting_surfaces.Ktotal, nl_body.Ktotal);
    UVLM::Matrix::AIC_sources(nl_body.zeta,
                              lifting_surfaces.zeta_col,
                              nl_body.longitudinals,
                              nl_body.perpendiculars,
                              nl_body.normals,
                              lifting_surfaces.longitudinals, //collocation
                              lifting_surfaces.perpendiculars, //collocation
                              lifting_surfaces.normals, //collocation
							  aic_nonlifting_on_lifting_x,
							  aic_nonlifting_on_lifting_y,
                              aic_nonlifting_on_lifting_z,
                              false);   
                              
    UVLM::Types::MatrixX aic_phantom_on_nonlifting = UVLM::Types::MatrixX::Zero(nl_body.Ktotal, phantom_surfaces.Ktotal);   
    UVLM::Matrix::AIC(phantom_surfaces.zeta,
                      nl_body.zeta_col,
                      phantom_surfaces.zeta_star,
                      nl_body.uext_col,
                      nl_body.normals,
                      options,
                      options.horseshoe,
                      aic_phantom_on_nonlifting);
                      
    UVLM::Types::copy_Mat_to_block(nl_body.aic_sources_z, aic, lifting_surfaces.Ktotal, lifting_surfaces.Ktotal);
    UVLM::Types::copy_Mat_to_block(aic_lifting_on_nonlifting, aic, lifting_surfaces.Ktotal,0);
    UVLM::Types::copy_Mat_to_block(aic_nonlifting_on_lifting_z, aic, 0, lifting_surfaces.Ktotal);
    UVLM::Types::copy_Mat_to_block(aic_phantom_on_nonlifting, aic, lifting_surfaces.Ktotal, Ktotal_lifting_and_nonlifting);
    
    }
    
    
    //Phantom panels     
    UVLM::Types::MatrixX aic_phantom_on_lifting = UVLM::Types::MatrixX::Zero(lifting_surfaces.Ktotal, phantom_surfaces.Ktotal);    
    phantom_surfaces.update_wake(lifting_surfaces.zeta_star);
    UVLM::Matrix::AIC(phantom_surfaces.zeta,
                      lifting_surfaces.zeta_col,
                      phantom_surfaces.zeta_star,
                      lifting_surfaces.uext_col,
                      lifting_surfaces.normals,
                      options,
                      options.horseshoe,
                      aic_phantom_on_lifting);
    // Get matrix to enforce linear interpolated circulation on phantom panels   
    UVLM::Types::MatrixX circulation_bc_phantom = UVLM::Types::MatrixX::Zero(phantom_surfaces.Ktotal, aic.cols());
    UVLM::Matrix::aic_phantom_interp_condition(lifting_surfaces.Ktotal,
                                               phantom_surfaces.Ktotal, 
                                               lifting_surfaces.zeta_col,
                                               phantom_surfaces.zeta_col,
                                               circulation_bc_phantom,
                                               phantom_surfaces.flag_zeta_phantom);
    // Merge all matrices into combined aic matrix
    UVLM::Types::copy_Mat_to_block(lifting_surfaces.aic, aic, 0, 0); 
    UVLM::Types::copy_Mat_to_block(aic_phantom_on_lifting, aic, 0, Ktotal_lifting_and_nonlifting);
    UVLM::Types::copy_Mat_to_block(circulation_bc_phantom, aic, Ktotal_lifting_and_nonlifting, 0);
}

template<typename t_mat_in>
uint UVLM::Matrix::get_total_VecVecMat_size(t_mat_in mat_in)
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

template<typename t_zeta_col, 
         typename t_zeta_phantom_col,
         typename t_aic_phantom_out,
         typename t_flag_zeta_phantom>
void UVLM::Matrix::aic_phantom_interp_condition
(
    const uint& Ktotal_lifting,
    const uint& Ktotal_phantom, 
    t_zeta_col& zeta_col,
    t_zeta_phantom_col& zeta_phantom_col,
    t_aic_phantom_out& aic_phantom_out,
    t_flag_zeta_phantom& flag_zeta_phantom,
    const bool only_for_update
)
{
     uint N_phantom_panels, M_phantom_panels, counter_row, N_spanwise_panels;
    double yL0_minus_yR0, interpolated_value;
    const uint n_surf = zeta_col.size();
    // set influence of phantom panel on own collocation point
    if (!only_for_update)
    {
        // if matrix is calculated due to update gamma (e.g. for unsteady wake convection)
        // we don`t need the identy matrix
        aic_phantom_out.topRightCorner(Ktotal_phantom, Ktotal_phantom).setIdentity();
        aic_phantom_out *= -1;
    }
    

    // get offsets beforehand 
    // Check if better dimensions and offsets are better to be stored in structs   
    std::vector<uint> offset_panel, offset_phantom_panel;
    UVLM::Types::VecDimensions dimensions_panel, dimensions_phantom_panel;
    UVLM::Types::generate_dimensions(zeta_col, dimensions_panel, 0);
    UVLM::Types::generate_dimensions(zeta_phantom_col, dimensions_phantom_panel, 0);
    UVLM::Matrix::build_offsets(n_surf, dimensions_panel,offset_panel);
	UVLM::Matrix::build_offsets(n_surf, dimensions_phantom_panel,offset_phantom_panel);
    UVLM::Types::Vector3 col_point_lifting, col_point_phantom;
    uint partner_surface = 0;
     for(uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        if (flag_zeta_phantom(0, i_surf) > i_surf)
        {
            partner_surface = flag_zeta_phantom(0, i_surf);
            counter_row = 0;
            N_spanwise_panels = zeta_col[i_surf][0].cols();
            N_phantom_panels= zeta_phantom_col[i_surf][0].rows();
            M_phantom_panels= zeta_phantom_col[i_surf][0].cols();
            
            for (uint i_row=0; i_row<N_phantom_panels; ++i_row)
            {
                yL0_minus_yR0 = zeta_col[1][1](i_row, 0) - zeta_col[0][1](i_row, 0);
                for (uint i_col=0; i_col<M_phantom_panels; ++i_col)
                {
                    interpolated_value = (zeta_phantom_col[i_surf][1](i_row, i_col)-zeta_col[i_surf][1](i_row, 0))/yL0_minus_yR0;
                    aic_phantom_out(offset_phantom_panel[i_surf]+counter_row,offset_panel[i_surf]+i_row*N_spanwise_panels) = 1-interpolated_value;
                    aic_phantom_out(offset_phantom_panel[i_surf]+counter_row,offset_panel[partner_surface]+i_row*N_spanwise_panels) = -interpolated_value;
                    interpolated_value = (zeta_phantom_col[partner_surface][1](i_row, i_col)-zeta_col[i_surf][1](i_row, 0))/yL0_minus_yR0;
                    // ToDO; find a better solution than multiplying values by minues 1 because the circulation is negativ there
                    aic_phantom_out(offset_phantom_panel[partner_surface]+counter_row,offset_panel[i_surf]+i_row*N_spanwise_panels) = (1-interpolated_value)*(-1.);
                    aic_phantom_out(offset_phantom_panel[partner_surface]+counter_row,offset_panel[partner_surface]+i_row*N_spanwise_panels) = interpolated_value;
                    
                    counter_row++;
                }                
            }
        }
    }
}
