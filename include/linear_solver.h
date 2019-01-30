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
                    Eigen::BiCGSTAB<t_a> solver;
                    solver.compute(a);
                    x = solver.solveWithGuess(b, x);
            } else
            {
                x = a.partialPivLu().solve(b);
            }
        }
    }
}
