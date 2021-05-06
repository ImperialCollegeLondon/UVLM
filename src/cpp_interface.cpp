#include "cpp_interface.h"
#include <fenv.h>


DLLEXPORT void run_VLM
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    double** p_zeta,
    double** p_zeta_star,
    double** p_zeta_dot,
    double** p_u_ext,
    double** p_gamma,
    double** p_gamma_star,
    double** p_forces,
    double* p_rbm_vel_g,
    double* p_centre_rot_g
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif
    struct UVLM::StructUtils::lifting_surface Lifting_surfaces = UVLM::StructUtils::lifting_surface
            (options.NumSurfaces,
            p_dimensions,
            p_zeta,
            p_u_ext,
            p_forces,
            p_zeta_star,
            p_zeta_dot,
            p_gamma,
            p_gamma_star,
            p_dimensions_star);

    UVLM::Steady::solver(Lifting_surfaces,
                         options,
                         flightconditions);
}

DLLEXPORT void run_VLM_nonlifting_body
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    double** p_zeta,
    double** p_u_ext,
    double** p_sigma,
    double** p_forces
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif

    struct UVLM::StructUtils::nonlifting_body nl_body = UVLM::StructUtils::nonlifting_body(options.NumSurfaces,
                                                    p_dimensions,
                                                    p_zeta,
                                                    p_u_ext,
                                                    p_forces,
                                                    p_sigma);
    UVLM::Steady::solver_nonlifting_body(nl_body,
                                         options,
                                         flightconditions);
}

DLLEXPORT void run_VLM_lifting_and_nonlifting_bodies
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    double** p_zeta,
    double** p_zeta_star,
    double** p_zeta_dot,
    double** p_u_ext,
    double** p_gamma,
    double** p_gamma_star,
    double** p_forces,
    int* p_flag_zeta_phantom,
    double* p_rbm_vel_g,
    unsigned int** p_dimensions_nonlifting,
    double** p_zeta_nonlifting,
    double** p_u_ext_nonlifting,
    double** p_sigma_nonlifting,
    double** p_forces_nonlifting
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif
    // Setup Lifting Surfaces
    struct UVLM::StructUtils::lifting_surface Lifting_surfaces = UVLM::StructUtils::lifting_surface
            (options.NumSurfaces,
            p_dimensions,
            p_zeta,
            p_u_ext,
            p_forces,
            p_zeta_star,
            p_zeta_dot,
            p_gamma,
            p_gamma_star,
            p_dimensions_star);
    
    // Setup Nonlifting Body  
    struct UVLM::StructUtils::nonlifting_body nl_body = UVLM::StructUtils::nonlifting_body(options.NumSurfacesNonlifting,
                                                    p_dimensions_nonlifting,
                                                    p_zeta_nonlifting,
                                                    p_u_ext_nonlifting,
                                                    p_forces_nonlifting,
                                                    p_sigma_nonlifting);
    // Setup Phantom Surfaces
    struct UVLM::StructUtils::phantom_surface phantom_surfaces = UVLM::StructUtils::phantom_surface(p_flag_zeta_phantom,
                                                                                                    Lifting_surfaces.n_surf,
                                                                                                    Lifting_surfaces.zeta,
                                                                                                    Lifting_surfaces.zeta_star,
                                                                                                    Lifting_surfaces.dimensions);
    
    //Start solver
	UVLM::Steady::solver_lifting_and_nonlifting_bodies
	(
		Lifting_surfaces,
		phantom_surfaces,
		options,
		flightconditions,
		nl_body
	);
}

DLLEXPORT void run_UVLM
(
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    unsigned int i_iter,
    double** p_u_ext,
    double** p_uext_star,
    double** p_zeta,
    double** p_zeta_star,
    double** p_zeta_dot,
    double*  p_rbm_vel_g,
    double** p_gamma,
    double** p_gamma_star,
    double** p_dist_to_orig,
    // double** p_previous_gamma,
    double** p_normals,
    double** p_forces,
    double** p_dynamic_forces
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif
    // Setup Lifting Surfaces
    struct UVLM::StructUtils::lifting_surface_unsteady Lifting_surfaces_unsteady = UVLM::StructUtils::lifting_surface_unsteady
            (options.NumSurfaces,
            p_dimensions,
            p_zeta,
            p_u_ext,
            p_forces,
            p_zeta_star,
            p_zeta_dot,
            p_gamma,
            p_gamma_star,
            p_dimensions_star,
            p_dist_to_orig,
            p_dynamic_forces,
            p_uext_star);

    UVLM::Unsteady::solver
    (
        i_iter,
        Lifting_surfaces_unsteady,
        options,
        flightconditions
    );
}

DLLEXPORT void run_UVLM_lifting_and_nonlifting
(
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    unsigned int i_iter,
    double** p_uext,
    double** p_uext_star,
    double** p_zeta,
    double** p_zeta_star,
    double** p_zeta_dot,
    double*  p_rbm_vel,
    double*  p_centre_rot,
    double** p_gamma,
    double** p_gamma_star,
    double** p_dist_to_orig,
    // double** p_previous_gamma,
    double** p_normals,
    double** p_forces,
    double** p_dynamic_forces,
    int* p_flag_zeta_phantom,
    unsigned int** p_dimensions_nonlifting,
    double** p_zeta_nonlifting,
    double** p_u_ext_nonlifting,
    double** p_sigma,
    double** p_forces_nonlifting
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif
    struct UVLM::StructUtils::lifting_surface_unsteady Lifting_surfaces_unsteady = UVLM::StructUtils::lifting_surface_unsteady
        (options.NumSurfaces,
        p_dimensions,
        p_zeta,
        p_uext,
        p_forces,
        p_zeta_star,
        p_zeta_dot,
        p_gamma,
        p_gamma_star,
        p_dimensions_star,
        p_dist_to_orig,
        p_dynamic_forces,
        p_uext_star);

    struct UVLM::StructUtils::nonlifting_body nl_body = UVLM::StructUtils::nonlifting_body
    (
        options.NumSurfacesNonlifting,
        p_dimensions_nonlifting,
        p_zeta_nonlifting,
        p_u_ext_nonlifting,
        p_forces_nonlifting,
        p_sigma
    );
                                                    
    struct UVLM::StructUtils::phantom_surface phantom_surfaces = UVLM::StructUtils::phantom_surface(p_flag_zeta_phantom,
                                                                                                    Lifting_surfaces_unsteady.n_surf,
                                                                                                    Lifting_surfaces_unsteady.zeta,
                                                                                                    Lifting_surfaces_unsteady.zeta_star,
                                                                                                    Lifting_surfaces_unsteady.dimensions);
    UVLM::Unsteady::solver_lifting_and_nonlifting
    (
        i_iter,
        Lifting_surfaces_unsteady,
        nl_body,
        phantom_surfaces,
        options,
        flightconditions
    );
}

DLLEXPORT void calculate_unsteady_forces
(
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    double** p_zeta,
    double** p_zeta_star,
    double*  p_rbm_vel,
    double** p_gamma,
    double** p_gamma_star,
    double** p_gamma_dot,
    double** p_normals,
    double** p_dynamic_forces
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif
    uint n_surf = options.NumSurfaces;
    UVLM::Types::VecDimensions dimensions;
    UVLM::Mapping::transform_dimensions(n_surf,
                                             p_dimensions,
                                             dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::Mapping::transform_dimensions(n_surf,
                                             p_dimensions_star,
                                             dimensions_star);

    UVLM::Types::VecVecMapX zeta;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_zeta,
                                      zeta,
                                      1);

    UVLM::Types::VecVecMapX zeta_star;
    UVLM::Mapping::map_VecVecMat(dimensions_star,
                                      p_zeta_star,
                                      zeta_star,
                                      1);

    UVLM::Types::MapVectorX rbm_velocity (p_rbm_vel, 2*UVLM::Constants::NDIM);

    UVLM::Types::VecMapX gamma;
    UVLM::Mapping::map_VecMat(dimensions,
                                   p_gamma,
                                   gamma,
                                   0);

    UVLM::Types::VecMapX gamma_star;
    UVLM::Mapping::map_VecMat(dimensions_star,
                                   p_gamma_star,
                                   gamma_star,
                                   0);

    UVLM::Types::VecMapX gamma_dot;
    UVLM::Mapping::map_VecMat(dimensions,
                                   p_gamma_dot,
                                   gamma_dot,
                                   0);

    UVLM::Types::VecVecMapX normals;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_normals,
                                      normals,
                                      0);

    UVLM::Types::VecVecMapX dynamic_forces;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_dynamic_forces,
                                      dynamic_forces,
                                      1,
                                      2*UVLM::Constants::NDIM);


    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);

    //std::cout << "Dynamic forces being calculated, new routine" << std::endl;
    UVLM::PostProc::calculate_dynamic_forces
    (
        zeta,
        zeta_star,
        zeta_col,
        gamma,
        gamma_star,
        gamma_dot,
        normals,
        dynamic_forces,
        options,
        flightconditions
    );
}

DLLEXPORT void UVLM_check_incidence_angle
(
    uint& n_surf,
    unsigned int** p_dimensions,
    double** p_uext,
    double** p_zeta,
    double** p_zeta_dot,
    double** p_normals,
    double*  p_rbm_vel,
    double** p_incidence_angle,
    const UVLM::Types::UVMopts& options
)
{
    UVLM::Types::VecDimensions dimensions;
    UVLM::Mapping::transform_dimensions(n_surf,
                                             p_dimensions,
                                             dimensions);

    UVLM::Types::VecVecMapX u_ext;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_uext,
                                      u_ext,
                                      1);

    UVLM::Types::VecVecMapX zeta;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_zeta,
                                      zeta,
                                      1);

    UVLM::Types::VecVecMapX zeta_dot;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_zeta_dot,
                                      zeta_dot,
                                      1);

    UVLM::Types::VecVecMapX normals;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_normals,
                                      normals,
                                      0);

    UVLM::Types::MapVectorX rbm_velocity (p_rbm_vel, 2*UVLM::Constants::NDIM);

    UVLM::Types::VecMapX incidence_angle;
    UVLM::Mapping::map_VecMat(dimensions,
                                   p_incidence_angle,
                                   incidence_angle,
                                   0);

    UVLM::PostProc::calculate_incidence_angle
    (
        u_ext,
        zeta,
        zeta_dot,
        normals,
        options,
        incidence_angle
    );
}

DLLEXPORT void total_induced_velocity_at_points
(
    const UVLM::Types::UVMopts& options,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    double** p_zeta,
    double** p_zeta_star,
    double** p_gamma,
    double** p_gamma_star,
    double* p_target_triads,
    double* p_uout,
    unsigned int npoints
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif
    uint n_surf = options.NumSurfaces;
    UVLM::Types::VecDimensions dimensions;
    UVLM::Mapping::transform_dimensions(n_surf,
                                             p_dimensions,
                                             dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::Mapping::transform_dimensions(n_surf,
                                             p_dimensions_star,
                                             dimensions_star);

    UVLM::Types::VecVecMapX zeta;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_zeta,
                                      zeta,
                                      1);

    UVLM::Types::VecVecMapX zeta_star;
    UVLM::Mapping::map_VecVecMat(dimensions_star,
                                      p_zeta_star,
                                      zeta_star,
                                      1);

    UVLM::Types::MapMatrixX uout(p_uout,
                                 npoints,
                                 UVLM::Constants::NDIM);

    UVLM::Types::MapMatrixX target_triads(p_target_triads,
                                          npoints,
                                          UVLM::Constants::NDIM);

    UVLM::Types::VecMapX gamma;
    UVLM::Mapping::map_VecMat(dimensions,
                                   p_gamma,
                                   gamma,
                                   0);

    UVLM::Types::VecMapX gamma_star;
    UVLM::Mapping::map_VecMat(dimensions_star,
                                   p_gamma_star,
                                   gamma_star,
                                   0);

    #pragma omp parallel for
    for (uint ipoint=0; ipoint<npoints; ipoint++)
    {
        UVLM::Types::Vector3 target_triad;
        UVLM::Types::Vector3 aux_uout;
        target_triad << target_triads(ipoint, 0),
                        target_triads(ipoint, 1),
                        target_triads(ipoint, 2);
        aux_uout = UVLM::BiotSavart::total_induced_velocity_on_point(target_triad,
                        zeta,
                        zeta_star,
                        gamma,
                        gamma_star,
                        options.ImageMethod,
                        options.vortex_radius);
        uout(ipoint, 0) = aux_uout(0);
        uout(ipoint, 1) = aux_uout(1);
        uout(ipoint, 2) = aux_uout(2);
    }

}

DLLEXPORT void multisurface
(
    const UVLM::Types::UVMopts& options,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_target,
    unsigned int** p_dimensions_uout,
    double** p_zeta,
    double** p_gamma,
    double** p_target_surface,
    double** p_uout
    // double  vortex_radius
)
{
#if defined(_OPENMP)
    omp_set_num_threads(options.NumCores);
#endif
    uint n_surf = options.NumSurfaces;
    uint n_surf_target = 1;

    UVLM::Types::VecDimensions dimensions;
    UVLM::Mapping::transform_dimensions(n_surf,
                                             p_dimensions,
                                             dimensions);

    // UVLM::Types::VecDimensions dimensions_star;
    // UVLM::Mapping::transform_dimensions(n_surf,
    //                                          p_dimensions_star,
    //                                          dimensions_star);

    UVLM::Types::VecDimensions dimensions_target;
    UVLM::Mapping::transform_dimensions(n_surf_target,
                                             p_dimensions_target,
                                             dimensions_target);

    UVLM::Types::VecDimensions dimensions_uout;
    UVLM::Mapping::transform_dimensions(n_surf_target,
                                             p_dimensions_uout,
                                             dimensions_uout);

    UVLM::Types::VecVecMapX zeta;
    UVLM::Mapping::map_VecVecMat(dimensions,
                                      p_zeta,
                                      zeta,
                                      1);

    // UVLM::Types::VecMapX gamma;
    // UVLM::Mapping::map_VecMat(dimensions,
    //                               p_gamma,
    //                               gamma,
    //                               0);
    //
    // UVLM::Types::VecMatrixX target_surface;
    // UVLM::Mapping::map_VecMatrixX(dimensions_target,
    //                                   p_target_surface,
    //                                   target_surface,
    //                                   1);

    UVLM::Types::VecMapX gamma;
    UVLM::Mapping::map_VecMat(dimensions,
                                  p_gamma,
                                  gamma,
                                  0);

    // std::cout << dimensions[0].first << std::endl;
    // std::cout << dimensions[0].second << std::endl;
    // std::cout << dimensions[1].first << std::endl;

    UVLM::Types::VecMapX uout;
    UVLM::Mapping::map_VecMat(dimensions_uout,
                                  p_uout,
                                  uout,
                                  0);

    UVLM::Types::VecVecMapX target_surface;
    UVLM::Mapping::map_VecVecMat(dimensions_target,
                                      p_target_surface,
                                      target_surface,
                                      1);
    // UVLM::Types::VecMatrixX target_surface_col;
    // UVLM::Geometry::generate_colocationMesh(target_surface, target_surface_col);
    //
    // UVLM::Types::VecMatrixX normal;
    // UVLM::Geometry::generate_surfaceNormal(zeta, normal);

    UVLM::Types::VecVecMatrixX target_surface_col;
    // UVLM::Types::allocate_VecVecMat
    //     (
    //         UVLM::Types::VecVecMatrixX& mat,
    //         const unsigned int& n_surf,
    //         const unsigned int& n_dim,
    //         const unsigned int& M,
    //         const unsigned int& N
    //     )
    UVLM::Geometry::generate_colocationMesh(target_surface, target_surface_col);

    UVLM::Types::VecVecMatrixX normal;
    UVLM::Types::allocate_VecVecMat(normal, target_surface_col);
    UVLM::Geometry::generate_surfaceNormal(target_surface, normal);

    UVLM::BiotSavart::multisurface(
        zeta[0],
        gamma[0],
        target_surface_col[0],
        uout[0],
        options.ImageMethod,
        normal[0],
        options.vortex_radius
    );
}


// linear UVLM interface

DLLEXPORT void call_der_biot_panel(double p_DerP[9],
                    double p_DerVertices[36],// 4x9
                    double p_zetaP[3],
                    double p_ZetaPanel[12],
                    const double& gamma,
                    double& vortex_radius)
  { /*
    To interface with python, matrices need to be mapped into 1d arrays.
    */

    int vv;

    UVLMlin::map_Mat3by3 DerP(p_DerP);
    UVLMlin::Vec_map_Mat3by3 DerVertices; // requires push_back to assigns
    const UVLMlin::map_RowVec3 zetaP(p_zetaP);
    const UVLMlin::map_Mat4by3 ZetaPanel(p_ZetaPanel);

    // initialise DerVertices - all done by reference
    for(vv=0;vv<4;vv++){
      DerVertices.push_back( UVLMlin::map_Mat3by3(p_DerVertices+9*vv) );
    }

    UVLMlin::der_biot_panel_map( DerP, DerVertices, zetaP, ZetaPanel, gamma, vortex_radius );
  }



DLLEXPORT void call_biot_panel(double p_vel[3],
                  double p_zetaP[3],
                  double p_ZetaPanel[12],
                  const double& gamma,
                  double& vortex_radius){
    /*
    To interface Eigen based routines with python, matrices need to be mapped
    into 1d arrays.
    */

    UVLMlin::map_RowVec3 velP(p_vel);
    const UVLMlin::map_RowVec3 zetaP(p_zetaP);
    const UVLMlin::map_Mat4by3 ZetaPanel(p_ZetaPanel);

    UVLMlin::biot_panel_map(velP, zetaP, ZetaPanel, gamma, vortex_radius);
  }



DLLEXPORT void call_dvinddzeta(double p_DerC[9],
                  double p_DerV[],
                  double p_zetaC[3],
                  double p_ZetaIn[],
                  double p_GammaIn[],
                  int& M_in,
                  int& N_in,
                  bool& IsBound,
                  int& M_in_bound, // M of bound surf associated
                  double& vortex_radius)
  {
    int cc;
    int Kzeta_in=(M_in+1)*(N_in+1);
    int Kzeta_in_bound=(M_in_bound+1)*(N_in+1);

    // interface
    UVLMlin::map_Mat3by3 DerC(p_DerC);
    const UVLMlin::map_RowVec3 zetaC(p_zetaC);

    UVLMlin::map_Mat DerV(p_DerV,3,3*Kzeta_in_bound);
    UVLMlin::map_Mat GammaIn(p_GammaIn,M_in,N_in);

    UVLMlin::Vec_map_Mat ZetaIn;
    for(cc=0;cc<3;cc++){
      ZetaIn.push_back( UVLMlin::map_Mat(p_ZetaIn+cc*Kzeta_in, M_in+1, N_in+1) );
    }

    UVLMlin::dvinddzeta( DerC,DerV,
          zetaC,ZetaIn,GammaIn,
          M_in,N_in,Kzeta_in,
          IsBound,M_in_bound,Kzeta_in_bound,
          vortex_radius);
  }



DLLEXPORT void call_aic3(  double p_AIC3[],
                double p_zetaC[3],
                double p_ZetaIn[],
                int& M_in,
                int& N_in,
                double& vortex_radius)
  {
    int cc;
    int K_in=M_in*N_in;

    UVLMlin::map_Mat AIC3(p_AIC3,3,K_in);
    const UVLMlin::map_RowVec3 zetaC(p_zetaC);

    int Kzeta_in=(M_in+1)*(N_in+1);
    UVLMlin::Vec_map_Mat ZetaIn;
    for(cc=0;cc<3;cc++){
      ZetaIn.push_back( UVLMlin::map_Mat(p_ZetaIn+cc*Kzeta_in, M_in+1, N_in+1) );
    }

    UVLMlin::aic3(AIC3, zetaC, ZetaIn, M_in, N_in, vortex_radius);
  }



DLLEXPORT void call_ind_vel(
                double p_vel[3],
                double p_zetaC[3],
                double p_ZetaIn[],
                double p_GammaIn[],
                int& M_in,
                int& N_in,
                double& vortex_radius)
  {
    int cc;

    UVLMlin::map_RowVec3 velC(p_vel);
    const UVLMlin::map_RowVec3 zetaC(p_zetaC);

    UVLMlin::map_Mat GammaIn(p_GammaIn,M_in,N_in);

    int Kzeta_in=(M_in+1)*(N_in+1);
    UVLMlin::Vec_map_Mat ZetaIn;
    for(cc=0;cc<3;cc++){
      ZetaIn.push_back( UVLMlin::map_Mat(p_ZetaIn+cc*Kzeta_in, M_in+1, N_in+1) );
    }

    UVLMlin::ind_vel(velC, zetaC, ZetaIn, GammaIn, M_in, N_in, vortex_radius);
  }
