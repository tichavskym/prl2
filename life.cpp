#include "mpi.h"
#include <string>
#include <fstream>

#define PADDING 2

#define DEBUG 1
#ifdef DEBUG
#define debug(...) \
    do {           \
        if (DEBUG) { \
            fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);  \
            fprintf(stderr, __VA_ARGS__);                    \
        } \
    } while (0);
#else
#define debug(...) do {} while(0);
#endif

struct Metadata {
    int width;
    int height;
    int nof_processes;
    int rows_per_process;
};

int idx(int x, int y, int width) {
    return y * (width + 2) + x;
}

// TODO if number of processes > nof_rows -> do sth about it

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

void print_local_grid(int* local_grid,int  width,int  rows_per_process, int rank) {
    for (int i = 0; i < (rows_per_process) + 2; i++) {
        printf("%d: ", rank);
        for (int j = 0; j < width + 2; j++) {
            printf("%d", local_grid[idx(j, i, width)]);
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

// Exchange data between processes
void exchange_data(int *grid, int width, int nof_processes, int rows_per_process, int rank) {
    if (rank % 2 == 0) {
        if (rank != nof_processes - 1){
            debug("%d: Exchanging data with process %d\n", rank, rank + 1);
            // send the boundary row to the next process
            MPI_Send(grid + rows_per_process * (width + PADDING), (width + PADDING), MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(grid + (rows_per_process + 1) * (width + PADDING) , (width + PADDING), MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // receive the boundary row from the previous process
        if (rank != 0) {
            debug("%d: Exchanging data with process %d\n", rank, rank - 1);
            MPI_Send(grid + (width + PADDING), (width + PADDING), MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(grid, (width + PADDING), MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        if (rank  != 0){
            debug("%d: Exchanging data with process %d\n", rank, rank - 1);
            // send the boundary row to the next process
            MPI_Recv(grid, (width + PADDING), MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(grid + (width + PADDING), (width + PADDING), MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
        }
        // receive the boundary row from the previous process
        if (rank != nof_processes - 1) {
            debug("%d: Exchanging data with process %d\n", rank, rank + 1);
            MPI_Recv(grid + (rows_per_process + 1) * (width + PADDING), (width + PADDING), MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(grid + rows_per_process * (width + PADDING), (width + PADDING), MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
        }
    }
}

void calculate(int *grid_even, int *grid_odd, int width, int height, int rows_per_process, int iterations, int rank, int nof_processes) {
    for (int i = 0; i < iterations; i++) {
        if (i % 2 == 0) {
            calculate_step(grid_even, grid_odd, width, rows_per_process, height);
            exchange_data(grid_odd, width, nof_processes, rows_per_process, rank);
        } else {
            calculate_step(grid_odd, grid_even, width, rows_per_process, height);
            exchange_data(grid_even, width, nof_processes, rows_per_process, rank);
        }
        debug("%d: Calculate iteration\n", rank);
        // todo send the boundary row to the next process and receive the boundary row from the next process
    }
}

void distribute_data_across_processes(int *grid, int width, int height, int rows_per_process, int nof_processes, int iterations ) {
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
            MPI_Send(
                    grid + (i * rows_per_process - 1) * (width + PADDING),
                    (rows_per_process + 2) * (width + PADDING),
                    MPI_INT,
                    i,
                    0,
                    MPI_COMM_WORLD
            );
        }
    }
}


int run(int argc, char **argv) {
    int rank, nof_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    if (rank == 0) {
        debug("Number of processes: %d\n", nof_processes);
        int iterations;
        char *filename;
        std::tie(filename, iterations) = parse_args(argc, argv);

        int width, height;
        std::tie(width, height) = get_dimensions(filename);
        int grid_even[(width + 2) * (height + 2)] = {0};
        int grid_odd[(width + 2) * (height + 2)] = {0};
        load_array(filename, grid_even, width, height);

        int rows_per_process = height / nof_processes;

        distribute_data_across_processes(grid_even, width, height, rows_per_process, nof_processes, iterations);
        int local_grid_even[(rows_per_process + 2) * (width + PADDING)];
        int local_grid_odd[(rows_per_process + 2) * (width + PADDING)] = {0};
        for (int i = 0; i < (rows_per_process + 2) * (width + PADDING); i++) {
            local_grid_even[i] = grid_even[i];
        }

        calculate(local_grid_even, local_grid_odd, width, height, rows_per_process, iterations, rank, nof_processes);

        int *grid;
        if (iterations % 2 == 1) {
            grid = local_grid_odd;
        } else {
            grid = local_grid_even;
        }

        MPI_Gather(grid + (width + PADDING), rows_per_process * (width + PADDING), MPI_INT, grid_even + (width + PADDING), rows_per_process * (width + PADDING), MPI_INT, 0, MPI_COMM_WORLD);

        debug("Result:\n");
        print_grid(grid_even, width, height, rows_per_process);
    } else {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int width, height, rows_per_process, iterations;
        MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rows_per_process, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&iterations, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int local_grid_size = (width + PADDING) * (rows_per_process + PADDING);
        int local_grid_even[local_grid_size] = {0};
        int local_grid_odd[local_grid_size] = {0};
        MPI_Recv(&local_grid_even, local_grid_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        calculate(local_grid_even, local_grid_odd, width, height, rows_per_process, iterations, rank, nof_processes);

        int *grid;
        if (iterations % 2 == 1) {
            grid = local_grid_odd;
        } else {
            grid = local_grid_even;
        }
        MPI_Gather(grid + (width + PADDING), rows_per_process * (width + PADDING), MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
    return 0;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int ret = run(argc, argv);
    MPI_Finalize();
    return ret;
}
