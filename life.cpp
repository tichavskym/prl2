#include "mpi.h"
#include <string>
#include <fstream>


int idx(int x, int y, int width) {
    return y * (width + 2) + x;
}

// TODO if number of processes > nof_rows -> do sth about it

std::tuple<int, int> get_dimensions(std::ifstream &f) {
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
            new_array[idx(x, y, width)] = array[idx(x, y, width)];
        }
    } else {
        if (alive_neigh == 3) {
            new_array[idx(x, y, width)] = 1;
        } else {
            new_array[idx(x, y, width)] = array[idx(x, y, width)];
        }
    }
}

void calculate_step(int *array, int *new_array, int width, int rows_per_process, int height) {
    // TODO what if last process
    for (int h = 1; h < rows_per_process + 1; h++) {
        printf("Calculating step for row %d\n", h);
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

int run(int argc, char **argv) {
    int rank, nof_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    printf("nof_processes %d\n", nof_processes);
    if (rank == 0) {
        if (argc != 3) {
            perror("Illegal number of arguments.");
            return 1;
        }

        std::ifstream grid_definition(argv[1]);
        if (!grid_definition.is_open()) {
            perror("Error while opening a file (first argument).");
            return 1;
        }

        int iterations;
        try {
            iterations = std::stoi(argv[2]);
        } catch (int exception) {
            perror("Error while converting second argument to integer.\n");
            return 1;
        }

        int width, height;
        std::tie(width, height) = get_dimensions(grid_definition);
        int array[(width + 2) * (height + 2)] = {0};
        int new_array[(width + 2) * (height + 2)] = {0};
        load_array(argv[1], array, width, height);

        printf("Beginning state:\n");
        for (int i = 0; i < height + 2; i++) {
            for (int j = 0; j < width + 2; j++) {
                printf("%d", array[idx(j, i, width)]);
            }
            printf("\n");
        }

        int rows_per_process = height / nof_processes;
        printf("rows per process %d\n", rows_per_process);

        // send data to each process
        for (int i = 1; i < nof_processes; i++) {
            if (i == nof_processes - 1) {
                MPI_Send(
                        array + i * rows_per_process, (height - rows_per_process * i) * width, MPI_INT, i, 0,
                        MPI_COMM_WORLD
                );
            } else {
                MPI_Send(array + i * rows_per_process, rows_per_process * width, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        printf("iterations %d\n", iterations);

        for (int i = 0; i < iterations; i++) {
            calculate_step(array, new_array, width, rows_per_process, height);
//            array = new_array;
            // send the boundary row to the next process and receive the boundary row from the next process
        }
        // calculate step, send the boundary row to the next process and receive the boundary row from the next process

        // gather the results
        // MPI_Gather();
        // send it to all processes
//        for (int i = 1; i < height + 1; i++) {
//            for (int j = 1; j < width + 1; j++) {
//                printf("%d", array[i][j]);
//            }
//        }
        printf("the result:\n");
        for (int i = 1; i < height + 1; i++) {
            for (int j = 1; j < width + 1; j++) {
                printf("%d", new_array[idx(j, i, width)]);
            }
            printf("\n");
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
