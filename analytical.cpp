#include <iostream>
#include <cstdlib>
#include <chrono>
#include "function.hpp"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " L grid_size" << std::endl;
        return 1;
    }
    double Lx, Ly, Lz;
    Lx = Ly = Lz = atof(argv[1]);
    int grid_size = atoi(argv[2]);
    int time_grid_size = 21;
    double T = static_cast<float>(time_grid_size) / static_cast<float>(grid_size) / 5.0;
    double dx = Lx / static_cast<float>(grid_size-1), dy, dz;
    dy = dz = dx;
    double dt = T / time_grid_size;

    Grid grid(grid_size, dx, dy, dz);
    Function function(grid_size, dx, dy, dz, Lx, Ly, Lz);
    
    std::chrono::high_resolution_clock::time_point start_time, end_time;
    start_time = std::chrono::high_resolution_clock::now();

    for (double t = 0; t < time_grid_size; ++t) {
        #pragma omp parallel for
        for (int i = 0; i < grid_size; i++) {
            #pragma omp parallel for
            for (int j = 0; j < grid_size; j++) {
                #pragma omp parallel for
                for (int k = 0; k < grid_size; k++) {
                    grid.set(i, j, k, function.u_from_idx(i, j, k, t * dt));
                }
            }
        }
        grid.save(t*dt, Lx, "analytical");
    }
    end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    double seconds = duration.count();

    std::cout << "Time taken: " << seconds << " seconds" << std::endl;

    return 0;
}