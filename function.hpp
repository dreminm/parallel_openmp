#include <cmath>
#include "grid.hpp"

class Function {
public:
    Function(int grid_size, double dx, double dy, double dz, double Lx, double Ly, double Lz)
    : grid_size(grid_size), dx(dx), dy(dy), dz(dz), Lx(Lx), Ly(Ly), Lz(Lz) {
        a_t = M_PI * std::sqrt(4.0 / std::pow(Lx, 2) + 16.0 / std::pow(Ly, 2) + 36.0 / std::pow(Lz, 2));
    }

    double u(double x, double y, double z, double t) const {
        return std::sin(x * 2.0 * M_PI / Lx) * std::sin(y * 4.0 * M_PI / Ly) * std::sin(z * 6.0 * M_PI / Lz) * std::cos(a_t * t);
    }

    double u_from_idx(int i, int j, int k, double t) const {
        //граничные условия
        if (i == grid_size-1)
            i = 0;
        if (j == grid_size-1)
            j = 0;
        if (k == grid_size-1)
            k = 0;
        return u(i * dx, j * dy, k * dz, t);
    }

    void fill_grid(Grid& grid, double t) const {
        for (int i = 0; i < grid_size; i++) {
            for (int j = 0; j < grid_size; j++) {
                for (int k = 0; k < grid_size; k++) {
                    grid.set(i, j, k, u_from_idx(i, j, k, t));
                }
            }
        }
    }

private:
    int grid_size;
    double dx, dy, dz;
    double Lx, Ly, Lz;
    double a_t;
};
