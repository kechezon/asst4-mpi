#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>

#include <mpi.h>

#include "wireroute.h"
#define ROOT 0

void print_stats(const std::vector<std::vector<int>>& occupancy) {
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto& row : occupancy) {
    for (const int count : row) {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}

void write_output(const std::vector<Wire>& wires, const int num_wires, const std::vector<std::vector<int>>& occupancy, const int dim_x, const int dim_y, const int nproc, std::string input_filename) {
  if (std::size(input_filename) >= 4 && input_filename.substr(std::size(input_filename) - 4) == ".txt") {
    input_filename.resize(std::size(input_filename) - 4);
  }

  const std::string occupancy_filename = input_filename + "_occupancy_" + std::to_string(nproc) + ".txt";
  const std::string wires_filename = input_filename + "_wires_" + std::to_string(nproc) + ".txt";

  std::ofstream out_occupancy(occupancy_filename, std::fstream::out);
  if (!out_occupancy) {
    std::cerr << "Unable to open file: " << occupancy_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_occupancy << dim_x << ' ' << dim_y << '\n';
  for (const auto& row : occupancy) {
    for (const int count : row) {
      out_occupancy << count << ' ';
    }
    out_occupancy << '\n';
  }

  out_occupancy.close();

  std::ofstream out_wires(wires_filename, std::fstream:: out);
  if (!out_wires) {
    std::cerr << "Unable to open file: " << wires_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n' << num_wires << '\n';

  for (const auto& [start_x, start_y, end_x, end_y, bend1_x, bend1_y] : wires) {
    out_wires << start_x << ' ' << start_y << ' ' << bend1_x << ' ' << bend1_y << ' ';

    if (start_y == bend1_y) {
    // first bend was horizontal

      if (end_x != bend1_x) {
        // two bends

        out_wires << bend1_x << ' ' << end_y << ' ';
      }
    } else if (start_x == bend1_x) {
      // first bend was vertical

      if (end_y != bend1_y) {
        // two bends

        out_wires << end_x << ' ' << bend1_y << ' ';
      }
    }
    out_wires << end_x << ' ' << end_y << '\n';
  }

  out_wires.close();
}

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();
  int pid;
  int nproc;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  std::string input_filename;
  double SA_prob = 0.1;
  int SA_iters = 5;
  char parallel_mode = '\0';
  int batch_size = 1;

  // Read command line arguments
  int opt;
  while ((opt = getopt(argc, argv, "f:p:i:m:b:")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;
      case 'p':
        SA_prob = atof(optarg);
        break;
      case 'i':
        SA_iters = atoi(optarg);
        break;
      case 'm':
        parallel_mode = *optarg;
        break;
      case 'b':
        batch_size = atoi(optarg);
        break;
      default:
        if (pid == ROOT) {
          std::cerr << "Usage: " << argv[0] << " -f input_filename [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
        }

        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || SA_iters <= 0 || (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0) {
    if (pid == ROOT) {
      std::cerr << "Usage: " << argv[0] << " -f input_filename [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  if (pid == ROOT) {
    std::cout << "Number of processes: " << nproc << '\n';
    std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
    std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
    std::cout << "Input file: " << input_filename << '\n';
    std::cout << "Parallel mode: " << parallel_mode << '\n';
    std::cout << "Batch size: " << batch_size << '\n';
  }

  int dim_x, dim_y, num_wires;
  std::vector<Wire> wires;
  std::vector<std::vector<int>> occupancy;

  if (pid == ROOT) {
      std::ifstream fin(input_filename);

      if (!fin) {
        std::cerr << "Unable to open file: " << input_filename << ".\n";
        exit(EXIT_FAILURE);
      }

      /* Read the grid dimension and wire information from file */
      fin >> dim_x >> dim_y >> num_wires;

      wires.resize(num_wires);
      for (auto& wire : wires) {
        fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
        wire.bend1_x = wire.start_x;
        wire.bend1_y = wire.start_y;
      }
  }

  /* Initialize any additional data structures needed in the algorithm */

  if (pid == ROOT) {
    const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
    std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
  }

  const auto compute_start = std::chrono::steady_clock::now();

  /** 
   * (TODO)
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Use MPI to parallelize the algorithm. 
   */
  
  //this will broadcast the grid dimensions and wire data to all of the processors
  MPI_Bcast(&dim_x, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&dim_y, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&num_wires, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  if (pid != ROOT) {
    wires.resize(num_wires);
  }

  //broadcast all the wire data
  for (int i = 0; i < num_wires; i++) {
    MPI_Bcast(&wires[i].start_x, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&wires[i].start_y, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&wires[i].end_x, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&wires[i].end_y, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&wires[i].bend1_x, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&wires[i].bend1_y, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  }

  //initialize the occcupancy matrix for all processors
  occupancy.resize(dim_y, std::vector<int>(dim_x, 0));


  //helper functions must be outside so all processors can use them
   //helper functions from asst3
    // INCLUSIVE
    auto horizontal_line = [&](int x1, int x2, int y, int val) {
      if (x1 < 0 || x2 >= dim_x || x2 < 0 || x2 >= dim_x || y < 0 || y >= dim_y) {
        printf("Tried to draw horizontal line out of bounds! (%i->%i, %i) in (%i, %i)\n", x1, x2, y, dim_x, dim_y);
        abort();
      }

      if (x1 <= x2) {
        for (int i = x1; i <= x2; i++) occupancy[y][i] += val;
      }
      else {
        for (int i = x1; i >= x2; i--) occupancy[y][i] += val;
      }

      return;
    };

    // INCLUSIVE
    auto vertical_line = [&](int y1, int y2, int x, int val) {
      if (y1 < 0 || y2 >= dim_y || y2 < 0 || y2 >= dim_y || x < 0 || x >= dim_x) {
          printf("Tried to draw vertical line out of bounds! (%i, %i->%i) in (%i, %i)\n", x, y1, y2, dim_x, dim_y);
          abort();
      }

      if (y1 <= y2) {
        for (int i = y1; i <= y2; i++) occupancy[i][x] += val;
      }
      else {
        for (int i = y1; i >= y2; i--) occupancy[i][x] += val;
      }

      return;
    };


    // To draw wires that have at least one bend
  auto draw_wire = [&](Wire wire, int val) {
    if (wire.start_y == wire.bend1_y) { // horizontal first
      if (wire.start_y == wire.end_y) {
          printf("Attempted to draw horizontal line... Please only use the draw_wire function for wires with at least one bend!");
          abort();
      }
      horizontal_line(wire.start_x, wire.bend1_x, wire.start_y, val);

      int line_start = wire.bend1_y + (wire.end_y > wire.start_y ? 1 : -1);
      vertical_line(line_start, wire.end_y, wire.bend1_x, val);

      if (wire.bend1_x != wire.end_x) {
        line_start = wire.bend1_x + (wire.end_x > wire.start_x ? 1 : -1);
        horizontal_line(line_start, wire.end_x, wire.end_y, val);
      }
    }

    else { assert(wire.start_x == wire.bend1_x); // vertical first
      if (wire.start_x == wire.end_x) {
          printf("Attempted to draw vertical line... Please only use the draw_wire function for wires with at least one bend!");
          abort();
      }
      vertical_line(wire.start_y, wire.bend1_y, wire.start_x, val);

      int line_start = wire.bend1_x + (wire.end_x > wire.start_x ? 1 : -1);
      horizontal_line(line_start, wire.end_x, wire.bend1_y, val);

      if (wire.bend1_y != wire.end_y) {
        line_start = wire.bend1_y + (wire.end_y > wire.start_y ? 1 : -1);
        vertical_line(line_start, wire.end_y, wire.end_x, val);
      }
    }
  };


  auto route_cost_help = [&](Wire wire, int bend1_x, int bend1_y) {
    int cost = 0;

    // first move is horizontal
    if (wire.start_y == bend1_y) {
      if (wire.start_x < wire.end_x)
        for (int i = wire.start_x; i <= bend1_x; i++) {
          cost += 2*(occupancy[wire.start_y][i]) + 1;
        }
      else
        for (int i = wire.start_x; i >= bend1_x; i--) {
          cost += 2*(occupancy[wire.start_y][i]) + 1;
        }

      // First bend
      if (wire.start_y < wire.end_y)
        for (int i = bend1_y + 1; i <= wire.end_y; i++) {
          cost += 2*(occupancy[i][bend1_x]) + 1;
        }
      else
        for (int i = bend1_y - 1; i >= wire.end_y; i--) {
          cost += 2*(occupancy[i][bend1_x]) + 1;
        }

      if (bend1_x == wire.end_x) return cost;

      // Second bend
      if (wire.start_x < wire.end_x)
        for (int i = bend1_x + 1; i <= wire.end_x; i++) {
          cost += 2*(occupancy[wire.end_y][i]) + 1;
        }
      else
        for (int i = bend1_x - 1; i >= wire.end_x; i--) {
          cost += 2*(occupancy[wire.end_y][i]) + 1;
        }

      return cost;
    }

    // first move is vertical
    else {
      if (wire.start_y < wire.end_y)
        for (int i = wire.start_y; i <= bend1_y; i++) {
          cost += 2*(occupancy[i][wire.start_x]) + 1;
        }
      else
        for (int i = wire.start_y; i >= bend1_y; i--) {
          cost += 2*(occupancy[i][wire.start_x]) + 1;
        }

      // First bend
      if (wire.start_x < wire.end_x)
        for (int i = bend1_x + 1; i <= wire.end_x; i++) {
          cost += 2*(occupancy[bend1_y][i]) + 1;
        }
      else
        for (int i = bend1_x - 1; i >= wire.end_x; i--) {
          cost += 2*(occupancy[bend1_y][i]) + 1;
        }

      if (bend1_y == wire.end_y) return cost;

      // Second bend
      if (wire.start_y < wire.end_y)
        for (int i = bend1_y + 1; i <= wire.end_y; i++) {
          cost += 2*(occupancy[i][wire.end_x]) + 1;
        }
      else
        for (int i = bend1_y - 1; i >= wire.end_y; i--) {
          cost += 2*(occupancy[i][wire.end_x]) + 1;
        }

      return cost;
    }
  };

  auto route_cost_iteration = [&](Wire wire, int route_idx, int dx, int dy, int num_routes) {
    int bend1_x = wire.start_x;
    int bend1_y = wire.start_y;

    // Determine minimum horizontal-first path
    if (route_idx < abs(dx)) { // Determine minimum horizontal-first path
      int shift = route_idx + 1;
      int dir = dx >= 0 ? 1 : -1;
      bend1_x += shift * dir;

      return route_cost_help(wire, bend1_x, wire.start_y);
    }
    else if (route_idx < num_routes) { // 3) Determine minimum vertical-first path
      int shift = route_idx - abs(dx) + 1;
      int dir = dy >= 0 ? 1 : -1;
      bend1_y += shift * dir;

      return route_cost_help(wire, wire.start_x, bend1_y);
    }

    printf("Attempting to find cost of invalid route %i on wire with %i routes\n", route_idx, num_routes);
    abort();
  };

  // Initialize occupancy with all wires, all processors will do this


  //main loop SA anealing
  for (int i = 0; i < SA_iters; i++) {

    //assign each processor to work on its share of wires
    int wires_per_process;
    int start_wire;
    int end_wire;

    for (int j = start_wire; j < end_wire; j += batch_size){
      int end_of_batch;

      //loop to route the wires in this batch
      for (int k = j; k < end_of_batch; k++){

        //find the best route for wire k using routing logic from asst3

        //if the route has changed, add the updates to the batch

      }
      // Propagate batch updates and also probably do .clear() on the batch (maybe fact check this)
    }

    


    //
  }

  //note that only root should do the updates
  if (pid == ROOT) {
    const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
    std::cout << "Computation time (sec): " << std::fixed << std::setprecision(10) << compute_time << '\n';

    /* Write wires and occupancy matrix to files */
    print_stats(occupancy);
    write_output(wires, num_wires, occupancy, dim_x, dim_y, nproc, input_filename);
  }

  // Cleanup
  MPI_Finalize();
}
