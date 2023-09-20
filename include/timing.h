#pragma once

#include <ctime>
#include <chrono>
#include <iostream>
#include <stdio.h>

/**
 * @file timing.h
 * @brief Contains a Timer class for measuring CPU and wall clock times.
 */

namespace UVLM
{
    namespace Timing
    {
        /**
         * @brief A Timer class for measuring CPU and wall clock times.
         */
        class Timer
        {
        private:
            std::clock_t cpu_tic; ///< CPU time at the start of timing.
            std::clock_t cpu_toc; ///< CPU time at the end of timing.
            double cpu_time; ///< CPU time elapsed in seconds.

            std::chrono::high_resolution_clock::time_point wall_tic; ///< Wall clock time at the start of timing.
            std::chrono::high_resolution_clock::time_point wall_toc; ///< Wall clock time at the end of timing.
            std::chrono::duration<float> wall_duration; ///< Wall clock time duration.

            /**
             * @brief Output timing results to the console.
             */
            void output()
            {
                printf("Finished in %6.3f seconds [CPU time]\n", this->cpu_time);
                printf("            %6.3f seconds [wall clock time]\n", this->wall_duration.count());
            };

        public:
            /**
             * @brief Default constructor for the Timer class.
             */
            Timer()
            { };

            /**
             * @brief Default destructor for the Timer class.
             */
            ~Timer()
            { };

            /**
             * @brief Initialize the Timer and start timing.
             */
            void tic()
            {
                this->cpu_tic = std::clock();
                this->wall_tic = std::chrono::high_resolution_clock::now();
            };

            /**
             * @brief Stop timing and display the results.
             */
            void toc()
            {
                this->cpu_toc = std::clock();
                this->cpu_time = (this->cpu_toc - this->cpu_tic) / (double)(CLOCKS_PER_SEC);

                this->wall_toc = std::chrono::high_resolution_clock::now();
                this->wall_duration = this->wall_toc - this->wall_tic;

                this->output();
            };
        };
    }
}
