#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <functional>

class Grid {
public:
    Grid(int grid_size, double dx, double dy, double dz)
        : Nx(grid_size), Ny(grid_size), Nz(grid_size), dx(dx), dy(dy), dz(dz), array(Nx * Ny * Nz, 0) {}

    ~Grid() = default;

    double get(int i, int j, int k) const {
        checkBounds(i, j, k);
        return array[index(i, j, k)];
    }

    void set(int i, int j, int k, double value) {
        checkBounds(i, j, k);
        array[index(i, j, k)] = value;
    }

    //семиточечный разностный аналог оператора Лапласа
    double get_diff(int i, int j, int k) const {
        checkBounds(i, j, k);

        auto diff = [&](int coord, int size, double delta, std::function<double(int)> getNeighbor) {
            double diff = 0.0;
            if (coord == 0 || coord == size - 1) {
                diff = (getNeighbor(size - 2) - 2 * getNeighbor(size - 1) + getNeighbor(1)) / (delta * delta);
            } else {
                diff = (getNeighbor(coord - 1) - 2 * getNeighbor(coord) + getNeighbor(coord + 1)) / (delta * delta);
            }
            return diff;
        };

        double diff_x = diff(i, Nx, dx, [&](int idx) { return this->get(idx, j, k); });
        double diff_y = diff(j, Ny, dy, [&](int idx) { return this->get(i, idx, k); });
        double diff_z = diff(k, Nz, dz, [&](int idx) { return this->get(i, j, idx); });

        return diff_x + diff_y + diff_z;
    }

    void swap(Grid& other) {
        std::swap(this->array, other.array);
    }

    void save(float T, float L, std::string postfix) {
        std::ostringstream formattedStream;
        formattedStream << "results/L_" << L << "_T_" << T << "_" << postfix << ".txt";
        std::string fileName = formattedStream.str();
        std::ofstream outFile(fileName);
        if (outFile.is_open()) {
            for (int i = 0; i < this->array.size(); i++) {
                outFile << this->array[i] << " ";
            }
            outFile.close();
            std::cout << "Vector saved to " << fileName << std::endl;
        } else {
            std::cerr << "Unable to open file: " << fileName << std::endl;
        }
    }

private:
    int Nx, Ny, Nz;
    std::vector<double> array;
    double dx, dy, dz;

    int index(int i, int j, int k) const {
        return i * (Ny * Nz) + j * Nz + k;
    }

    void checkBounds(int i, int j, int k) const {
        if (i < 0 || i >= Nx || j < 0 || j >= Ny || k < 0 || k >= Nz) {
            std::cout << i << " " << j << " " << k << std::endl << std::flush;
            throw std::out_of_range("Index out of bounds.");
        }
    }
};