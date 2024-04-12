#include "mpi.h"
#include <string>
#include <fstream>


int idx(int x, int y, int width) {
    return y * (width + 2) + x;
}

// TODO if number of processes > nof_rows -> do sth about it

std::tuple<int, int> get_dimensions(char * filename) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        throw std::invalid_argument("Error while opening a file (first argument).");
    }

    int width, height = 1;

    std::string line;
    std::getline(f, line);
    width = line.size();

    while (std::getline(f, line)) {
        height++;
        if (line.size() != width) {
            perror("Input file doesn't represent rectangular grid.");
            return std::make_tuple(0, 0);
        }
    }

    return std::make_tuple(width, height);
}

/** Get number of alive neighbours for a cell given by coordinates x and y.
 *
 * @param array representing grid
 * @param x
 * @param y
 * @param width of the grid
 * */
int get_number_of_alive_neighbours(int *array, int x, int y, int width) {
    int count = 0;
    for (int i = x - 1; i < x + 2; i++) {
        for (int j = y - 1; j < y + 2; j++) {
            if (j == y && i == x) {
                continue;
            }
            if (array[idx(i, j, width)] == 1) {
                count++;
            }
        }
    }
    return count;
}

void calculate_point(int *array, int *new_array, int x, int y, int width, int height) {
    int alive_neigh = get_number_of_alive_neighbours(array, x, y, width);
    if (array[idx(x, y, width)] == 1) {
        if (alive_neigh < 2 || alive_neigh > 3) {
            new_array[idx(x, y, width)] = 0;
        } else {
            new_array[idx(x, y, width)] = 1;
        }
    } else {
        if (alive_neigh == 3) {
            new_array[idx(x, y, width)] = 1;
        } else {
            new_array[idx(x, y, width)] = 0;
        }
    }
}

void calculate_step(int *array, int *new_array, int width, int rows_per_process, int height) {
    // TODO what if last process
    for (int h = 1; h < rows_per_process + 1; h++) {
        for (int w = 1; w < width + 1; w++) {
            calculate_point(array, new_array, w, h, width, height);
        }
    }
}

void load_array(char *filename, int *array, int width, int height) {
    std::ifstream grid_definition(filename);
    std::string line;

    int h = 1;
    while (std::getline(grid_definition, line)) {
        array[idx(0, h, width)] = 0;
        // TODO raise exception printf("line size: %zu, width %d\n", line.size(), width);
        for (int w = 0; w < line.size(); w++) {
            int x = line[w] - '0';
            array[idx(w + 1, h, width)] = x;
        }
        h++;
    }
}

void print_grid(int *grid, int width, int height, int rows_per_process) {
    int row = 0;
    int rows = 0;
    for (int i = 1; i < height + 1; i++) {
        if (rows == rows_per_process) {
            row++;
            rows = 0;
        }
        printf("%d: ", row);
        rows++;

        for (int j = 1; j < width + 1; j++) {
            printf("%d", grid[idx(j, i, width)]);
        }
        printf("\n");
    }
}

std::tuple<char *, int> parse_args(int argc, char *argv[]) {
    if (argc != 3) {
        throw std::invalid_argument("Illegal number of arguments.");
    }

    int iterations;
    try {
        iterations = std::stoi(argv[2]);
    } catch (int exception) {
        throw std::invalid_argument("Second argument must be an integer.");
    }
    return std::make_tuple(argv[1], iterations);
}

int run(int argc, char **argv) {
    int rank, nof_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    printf("nof_processes %d\n", nof_processes);
    if (rank == 0) {
        int iterations;
        char *filename;
        std::tie(filename, iterations) = parse_args(argc, argv);

        int width, height;
        std::tie(width, height) = get_dimensions(filename);
        int grid_even[(width + 2) * (height + 2)] = {0};
        int grid_odd[(width + 2) * (height + 2)] = {0};
        load_array(argv[1], grid_even, width, height);

//        printf("Beginning state:\n");
//        for (int i = 0; i < height + 2; i++) {
//            for (int j = 0; j < width + 2; j++) {
//                printf("%d", grid_even[idx(j, i, width)]);
//            }
//            printf("\n");
//        }

        int rows_per_process = height / nof_processes;
        printf("rows per process %d\n", rows_per_process);

        // send data to each process
        for (int i = 1; i < nof_processes; i++) {
            if (i == nof_processes - 1) {
                MPI_Send(
                        grid_even + i * rows_per_process, (height - rows_per_process * i) * width, MPI_INT, i, 0,
                        MPI_COMM_WORLD
                );
            } else {
                MPI_Send(grid_even + i * rows_per_process, rows_per_process * width, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }

        for (int i = 0; i < iterations; i++) {
            if (i % 2 == 0) {
                calculate_step(grid_even, grid_odd, width, rows_per_process, height);
            } else {
                calculate_step(grid_odd, grid_even, width, rows_per_process, height);
            }
            // todo send the boundary row to the next process and receive the boundary row from the next process
        }

        // gather the results
        // MPI_Gather();
        // send it to all processes
//        for (int i = 1; i < height + 1; i++) {
//            for (int j = 1; j < width + 1; j++) {
//                printf("%d", grid_even[i][j]);
//            }
//        }
        if (iterations % 2 == 1) {
            print_grid(grid_odd, width, height, rows_per_process);
        } else {
            print_grid(grid_even, width, height, rows_per_process);
        }
    } else {

        //MPI_Recv();

        // for each step
        // calculate step, receive the boundary row from the previous process, send boundary row to the previous process,
        // send boundary row to the next process, receive boundary row from the next process

        // odd even do it in opposite order
    }
    return 0;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int ret = run(argc, argv);
    MPI_Finalize();
    return ret;
}
