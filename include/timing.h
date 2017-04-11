#pragma once

#include<ctime>
#include<chrono>
#include<iostream>
#include<stdio.h>

namespace UVLM
{
    namespace Timing
    {
        class Timer
        {
        private:
            std::clock_t cpu_tic;
            std::clock_t cpu_toc;
            double cpu_time;

            std::chrono::high_resolution_clock::time_point wall_tic;
            std::chrono::high_resolution_clock::time_point wall_toc;
            std::chrono::duration<float> wall_duration;

            void output()
            {
                // std::cout << "Finished in " << this->cpu_time << "seconds [CPU time]" << std::endl;
                printf("Finished in %6.3f seconds [CPU time]\n", this->cpu_time);
                // std::cout << "            " << this->wall_duration.count() << "seconds [wall time]" << std::endl;
                printf("            %6.3f seconds [wall clock time]\n", this->wall_duration.count());
            };
        public:
            // constructor
            Timer()
            { };
            // destructor
            ~Timer()
            { };
            // init Timer
            void tic()
            {
                this->cpu_tic = std::clock();
                this->wall_tic = std::chrono::high_resolution_clock::now();
            };
            void toc()
            {
                this->cpu_toc = std::clock();
                this->cpu_time = (this->cpu_toc - this->cpu_tic)/(double)(CLOCKS_PER_SEC);

                this->wall_toc = std::chrono::high_resolution_clock::now();
                this->wall_duration = this->wall_toc - this->wall_tic;

                this->output();
            };
        };
    }
}
