#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "Eigen/IterativeLinearSolvers"

namespace UVLM
{
    namespace LinearSolver
    {
        template <typename t_a,
                  typename t_b,
                  typename t_x,
                  typename t_options>
        void solve_system
        (
            t_a& a,
            t_b& b,
            t_options& options,
            t_x& x
        )
        {
            if (options.iterative_solver)
            {
                // we use iterative solver
                // if (options.iterative_precond)
                // {
                //     Eigen::BiCGSTAB<t_a, Eigen::IncompleteLUT<UVLM::Types::Real, int>> solver;
                //     // Eigen::ConjugateGradient<t_a, Eigen::IncompleteLUT<UVLM::Types::Real>> solver;
                //     solver.preconditioner().setDroptol(options.iterative_tol);
                //     solver.compute(a);
                //     x = solver.solve(b);
                // } else
                // {
                    // std::cout << "iterative " << std::endl;
                    Eigen::ConjugateGradient<t_a, Eigen::Lower|Eigen::Upper> solver;
                    // Eigen::BiCGSTAB<t_a> solver;
                    solver.compute(a);
                    x = solver.solve(b);
                // }
            } else
            {
                    // std::cout << "direct" << std::endl;
                x = a.partialPivLu().solve(b);
            }
        }
    }
}
