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

    Grid grid_cur(grid_size, dx, dy, dz), grid_prev(grid_size, dx, dy, dz), grid_next(grid_size, dx, dy, dz);
    Function function(grid_size, dx, dy, dz, Lx, Ly, Lz);
    
    function.fill_grid(grid_prev, 0);
    function.fill_grid(grid_cur, dt);
    
    double max_err = 0;
    double val, err;

    std::chrono::high_resolution_clock::time_point start_time, end_time;

    start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for reduction(max:max_err)
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            for (int k = 0; k < grid_size; k++) {
                double val = grid_prev.get(i, j, k) + (dt * dt) * grid_prev.get_diff(i, j, k) / 2.0;
                grid_cur.set(i, j, k, val);
                err = std::fabs(function.u_from_idx(i, j, k, dt) - val);
                max_err = std::max(err, max_err);
            }
        }
    }

    // grid_prev.save(0, 1.0, "_computational");
    // grid_cur.save(dt, 1.0, "_computational");

    for (double t = 2; t < time_grid_size; ++t) {
        #pragma omp parallel for reduction(max:max_err)
        for (int i = 0; i < grid_size; i++) {
            for (int j = 0; j < grid_size; j++) {
                for (int k = 0; k < grid_size; k++) {
                    val = (dt * dt) * grid_cur.get_diff(i, j, k) + 2 * grid_cur.get(i, j, k) - grid_prev.get(i, j, k);
                    grid_next.set(i, j, k, val);
                    err = std::fabs(function.u_from_idx(i, j, k, t * dt) - val);
                    max_err = std::max(err, max_err);
                }
            }
        }
        grid_prev.swap(grid_cur);
        grid_cur.swap(grid_next);
        // grid_next.save(t*dt, 1.0, "_computational");
        std::cout << "Maximum error: " << max_err << std::endl;
    }

    std::cout << "Maximum error: " << max_err << std::endl;

    end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    double seconds = duration.count();

    std::cout << "Time taken: " << seconds << " seconds" << std::endl;


    return 0;
}