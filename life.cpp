/** PRL 2nd project: Game of Life
 * Author: Milan Tichavský (xticha09)
 *
 * Implementation with a fixed boundary (the cells behind the boundary are considered to have 0 values).
 *
 * Design:
 * The work of processes is split by rows. At first, the 0th process loads the file with the initial grid definition
 * and sends each process the same amount of rows to calculate. In case the number of rows cannot be fully divided by
 * the number of processes, the last process is responsible for calculating all the remaining rows.
 *
 * The internal grid, which has certain `width` and `height`, is padded internally from each side by zeros
 * to simplify the calculations (therefore, you'll see a lot of usages of PADDING macro, which says, that for
 * the internal representation of the grid: width_internal = width + PADDING and height_internal = height + PADDING).
 *
 * Also, this means that the boundaries of the grid are defined by the dimensions of the grid on input.
 * The grid is stored in a 1D array, where the 2D indexing is possible via the `idx` function.
 *
 * During calculation, the neighbouring processes exchange the information about boundary rows between each other.
 * For more information, see `exchange_data method`. After all iterations, the results are sent back to the zeroth
 * process, which prints out the final grid.
 */

#include "mpi.h"
#include <string>
#include <fstream>

#define PADDING 2


#ifdef DEBUG
#define debug(...) \
    do {           \
        if (DEBUG) { \
            fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);  \
            fprintf(stderr, __VA_ARGS__);                    \
            fflush(stderr);                                  \
        } \
    } while (0);
#else
#define debug(...) do {} while(0);
#endif

/* Enable two 2D indexing in 1D array representing the grid */
int idx(int x, int y, int width) {
    return y * (width + 2) + x;
}

#ifdef DEBUG
/** Function used for debugging purposes only */
void debug_print_local_grid(int *local_grid, int width, int rows_per_process, int rank) {
    for (int i = 0; i < (rows_per_process) + 2; i++) {
        printf("%d: ", rank);
        for (int j = 0; j < width + 2; j++) {
            printf("%d", local_grid[idx(j, i, width)]);
        }
        printf("\n");
        fflush(stdout);
    }
}
#endif

/* Get dimensions of the grid stored in file called `filename` */
std::tuple<int, int> get_dimensions(char *filename) {
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

/* Calculate a new value of point located at `x`, `y` in `array` and store it into `new_array` */
void calculate_point(int *array, int *new_array, int x, int y, int width) {
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

/* Get last process local grid height padding included */
int get_last_process_grid_height(int height, int rows_per_process, int nof_processes) {
    return (height - (rows_per_process * (nof_processes - 1))) + PADDING;
}

int is_last_proc(int rank, int nof_processes) {
    return rank == nof_processes - 1;
}

void calculate_step(int *array, int *new_array, int width, int rows_per_process, int height, int rank, int nof_processes) {
    int height_iterate;
    if (is_last_proc(rank, nof_processes)) {
        height_iterate = get_last_process_grid_height(height, rows_per_process, nof_processes) - 1;
    } else {
        height_iterate = rows_per_process + 1;
    }

    for (int h = 1; h < height_iterate; h++) {
        for (int w = 1; w < width + 1; w++) {
            calculate_point(array, new_array, w, h, width);
        }
    }
}

/** Load grid into `array` from file `filename` */
void load_grid(char *filename, int *array, int width) {
    std::ifstream grid_definition(filename);
    std::string line;

    int h = 1;
    while (std::getline(grid_definition, line)) {
        array[idx(0, h, width)] = 0;
        for (int w = 0; w < line.size(); w++) {
            int x = line[w] - '0';
            array[idx(w + 1, h, width)] = x;
        }
        h++;
    }
}

/** Print grid (ignores padding) */
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
        fflush(stdout);
    }
}

/** Returns filename of file containing the grid and number of iterations */
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

/** Exchange data between processes
 *
 * Each process needs to exchange data with the two neighbouring processes. The data being exchanged are boundary rows,
 * that are needed for proper calculation of the rows on the edges of the local grid (=part of the grid that given
 * process calculates).
 *
 * The sending is always done from even process to odd process first, and than followed the other way around.
 * In the first step, even processes exchange data with the following process (rank + 1), and in the second step they
 * exchange data with the previous process (rank - 1).
 *
 * Zeroth and the last process exchange data from only one side (the other doesn't make sense).
 * */
void exchange_data(int *grid, int width, int nof_processes, int rows_per_process, int rank) {
    int payload_size = (width + PADDING);
    if (rank % 2 == 0) {
        if (rank != nof_processes - 1) {
            // send the boundary row to the next process
            MPI_Send(grid + rows_per_process * (width + PADDING), payload_size, MPI_INT, rank + 1, 0,
                     MPI_COMM_WORLD);
            MPI_Recv(grid + (rows_per_process + 1) * (width + PADDING), payload_size, MPI_INT, rank + 1, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // receive the boundary row from the previous process
        if (rank != 0) {
            MPI_Send(grid + (width + PADDING), payload_size, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(grid, payload_size, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        if (rank != 0) {
            // send the boundary row to the next process
            MPI_Recv(grid, payload_size, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(grid + (width + PADDING), payload_size, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
        }
        // receive the boundary row from the previous process
        if (rank != nof_processes - 1) {
            // debug_print(rank, rank + 1, grid + rows_per_process * (width + PADDING), payload_size);
            MPI_Recv(grid + (rows_per_process + 1) * (width + PADDING), payload_size, MPI_INT, rank + 1, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(grid + rows_per_process * (width + PADDING), payload_size, MPI_INT, rank + 1, 0,
                     MPI_COMM_WORLD);
        }
    }
}

/** Calculate the grid for `iterations` iterations
 *
 * @param grid_even, grid_odd this method needs two arrays for doing calculations
 * @param width, height dimensions of the grid
 * @param rows_per_process
 * @param iterations
 * @param rank the rank of the process executing this method
 * @param nof_processes total number of processes
 */
int *calculate(int *grid_even, int *grid_odd, int width, int height, int rows_per_process, int iterations, int rank,
               int nof_processes) {
    for (int i = 0; i < iterations; i++) {
        if (i % 2 == 0) {
            calculate_step(grid_even, grid_odd, width, rows_per_process, height, rank, nof_processes);
            exchange_data(grid_odd, width, nof_processes, rows_per_process, rank);
        } else {
            calculate_step(grid_odd, grid_even, width, rows_per_process, height, rank, nof_processes);
            exchange_data(grid_even, width, nof_processes, rows_per_process, rank);
        }
    }

    if (iterations % 2 == 1) {
        return grid_odd;
    } else {
        return grid_even;
    }
}

/** Send parts of the grid to each process for computation */
void distribute_data_across_processes(int *grid, int width, int height, int rows_per_process, int nof_processes,
                                      int iterations) {
    for (int i = 1; i < nof_processes; i++) {
        MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&height, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&rows_per_process, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&iterations, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

        if (i == (nof_processes - 1)) {
            int *send_start = grid + (i * rows_per_process) * (width + PADDING);
            int send_length = ((height - (rows_per_process * i)) + PADDING) * (width + PADDING);
            MPI_Send(
                    send_start,
                    send_length,
                    MPI_INT, i, 0,
                    MPI_COMM_WORLD
            );
        } else {
            int *send_start = grid + (i * rows_per_process) * (width + PADDING);
            int send_length = (rows_per_process + 2) * (width + PADDING);
            MPI_Send(send_start, send_length, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
}

int get_local_grid_size(int nof_processes, int width,  int height, int rows_per_process, int rank) {
    if (rank == (nof_processes - 1)) {
        return get_last_process_grid_height(height, rows_per_process, nof_processes) *
                          (width + PADDING);
    } else {
        return (width + PADDING) * (rows_per_process + PADDING);
    }
}

/** Gather results into `result` array
 *
 * @param result
 * @param local_grid the local grid of 0th process
 * @param width, height dimensions of the grid
 * @param nof_processes
 * @param rows_per_process
 */
void gather_results(int* result, int *local_grid, int width, int height, int nof_processes, int rows_per_process) {
    int local_grid_size = get_local_grid_size(nof_processes, width, height, rows_per_process, 0);
    for (int i = 0; i < local_grid_size; i++) {
        result[i] = local_grid[i];
    }
    for (int i = 1; i < nof_processes; i++) {
        local_grid_size = get_local_grid_size(nof_processes, width, height, rows_per_process, i);
        MPI_Recv(result + (i * rows_per_process + 1) * (width + PADDING),
                 local_grid_size - (width + PADDING) * PADDING,
                 MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );
    }
}


int run(int argc, char **argv) {
    int rank, nof_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    if (rank == 0) {
        int iterations;
        char *filename;
        std::tie(filename, iterations) = parse_args(argc, argv);

        int width, height;
        std::tie(width, height) = get_dimensions(filename);

        int grid_plus_padding_size = (width + 2) * (height + 2);
        int grid_even[grid_plus_padding_size];
        memset(grid_even, 0, grid_plus_padding_size * sizeof(int));
        int grid_odd[grid_plus_padding_size];
        memset(grid_odd, 0, grid_plus_padding_size * sizeof(int));

        load_grid(filename, grid_even, width);

        int rows_per_process = height / nof_processes;

        distribute_data_across_processes(grid_even, width, height, rows_per_process, nof_processes, iterations);

        // The start of computation of the split for 0th process
        int local_grid_size = (rows_per_process + PADDING) * (width + PADDING);
        int local_grid_even[local_grid_size];
        int local_grid_odd[local_grid_size];
        memset(local_grid_odd, 0, local_grid_size * sizeof(int));

        for (int i = 0; i < local_grid_size; i++) {
            local_grid_even[i] = grid_even[i];
        }

        int *grid = calculate(local_grid_even, local_grid_odd, width, height, rows_per_process, iterations, rank,
                              nof_processes);

        gather_results(grid_even, grid, width, height, nof_processes, rows_per_process);
        print_grid(grid_even, width, height, rows_per_process);
    } else {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int width, height, rows_per_process, iterations;
        MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rows_per_process, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&iterations, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int local_grid_size;
        if (rank == (nof_processes - 1)) {
            local_grid_size = ((height - (rows_per_process * (nof_processes - 1))) + PADDING) * (width + PADDING);
        } else {
            local_grid_size = (width + PADDING) * (rows_per_process + PADDING);
        }

        int local_grid_even[local_grid_size];
        memset(local_grid_even, 0, local_grid_size * sizeof(int));
        int local_grid_odd[local_grid_size];
        memset(local_grid_odd, 0, local_grid_size * sizeof(int));

        MPI_Recv(&local_grid_even, local_grid_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int *grid = calculate(local_grid_even, local_grid_odd, width, height, rows_per_process, iterations, rank,
                              nof_processes);

        MPI_Send(grid + (width + PADDING), local_grid_size - (width + PADDING) * PADDING, MPI_INT, 0, 0,
                 MPI_COMM_WORLD);
    }
    return 0;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int ret = run(argc, argv);
    MPI_Finalize();
    return ret;
}
