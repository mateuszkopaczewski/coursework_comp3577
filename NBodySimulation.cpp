#include "NBodySimulation.h"
#include <vector>

NBodySimulation::NBodySimulation () :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0) {};

NBodySimulation::~NBodySimulation () {
  if (x != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] x[i];
    delete [] x;
  }
  if (v != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] v[i];
    delete [] v;
  }
  if (mass != nullptr) {
    delete [] mass;
  }
}

void NBodySimulation::checkInput(int argc, char** argv) {
    if (argc==1) {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:       interval after how many time units to plot."
                 " Use 0 to switch off plotting" << std::endl
              << "  final-time:      simulated time (greater 0)" << std::endl
              << "  dt:              time step size (greater 0)" << std::endl
              << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
              << std::endl
              << "Examples of arguments:" << std::endl
              << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body spiralling around the other" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup from first lecture" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
              << std::endl;

    throw -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each body is given by seven entries"
                 " (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments"
                 " (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    throw -2;
  }
}

void NBodySimulation::setUp (int argc, char** argv) {

  checkInput(argc, argv);

  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies"
            << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta
              << " time units" << std::endl;
    tPlot = 0.0;
  }
}


double NBodySimulation::force_calculation (int i, int j, int direction){
  // Euclidean distance
  const double distance = sqrt(
                               (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                               (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                               (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                               );
  const double distance3 = distance * distance * distance;
  minDx = std::min( minDx,distance );

  return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
}

void NBodySimulation::updateBody () {

  timeStepCounter++;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  // Allocate arrays for forces acting on each body along x,y,z directions.
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];

  const double cutoff = 2.0;
  const double skin = 0.2;
  const double cutoff2 = (cutoff + skin) * (cutoff + skin);

  std::vector<std::vector<int>> neighbors(NumberOfBodies);

  // Construct the neighbor list using the variable Linked-Cell Algorithm with linked list.
  std::vector<int> head(NumberOfCells, -1);
  std::vector<int> next(NumberOfBodies, -1);
  std::vector<int> cellIndex(NumberOfBodies, -1);

  for (int i = 0; i < NumberOfBodies; i++) {
    // Calculate the index of the cell that contains the body.
    int ix = static_cast<int>(floor((x[i][0] - xmin) / cellSize));
    int iy = static_cast<int>(floor((x[i][1] - ymin) / cellSize));
    int iz = static_cast<int>(floor((x[i][2] - zmin) / cellSize));
    int cell = ix + nx*(iy + ny*iz);

    // Link the body to the head of the corresponding cell.
    next[i] = head[cell];
    head[cell] = i;
    cellIndex[i] = cell;
  }

  // Calculate the forces acting on each body.
  for (int i=0; i<NumberOfBodies; i++) {
    force0[i] = 0.0;
    force1[i] = 0.0;
    force2[i] = 0.0;

    // Calculate the index of the cell that contains the body.
    int ix = static_cast<int>(floor((x[i][0] - xmin) / cellSize));
    int iy = static_cast<int>(floor((x[i][1] - ymin) / cellSize));
    int iz = static_cast<int>(floor((x[i][2] - zmin) / cellSize));
    int cell = ix + nx*(iy + ny*iz);

    // Loop over the neighboring cells.
    for (int cx = ix-1; cx <= ix+1; cx++) {
      for (int cy = iy-1; cy <= iy+1; cy++) {
        for (int cz = iz-1; cz <= iz+1; cz++) {
          if (cx < 0 || cx >= nx || cy < 0 || cy >= ny || cz < 0 || cz >= nz){
            continue; } 
          int cell2 = cx + nx*(cy + ny*cz);
              // Loop over the bodies in the neighboring cell.
              int j = head[cell2];
              while (j != -1) {
                // Compute the distance between body i and body j.
                if (i != j) {
                  double dx = x[j][0] - x[i][0];
                  double dy = x[j][1] - x[i][1];
                  double dz = x[j][2] - x[i][2];
                  double r2 = dx*dx + dy*dy + dz*dz;

                  // Check if j is within the cutoff distance from i.
                  if (r2 < cutoff2) {
                    // Add j to the neighbor list of i.
                    neighbors[i].push_back(j);

                    // Calculate the force acting on i due to j.
                    double r = sqrt(r2);
                    double rinv = 1.0 / r;
                    double rinv2 = rinv * rinv;
                    double rinv3 = rinv2 * rinv;
                    double f = G * m[i] * m[j] * rinv2;
                    force0[i] += f * dx * rinv;
                    force1[i] += f * dy * rinv;
                    force2[i] += f * dz * rinv;
                  }
                }
                j = next[j]; // Move to the next body in the neighboring cell.
              }
            }
          }
        }

        // Update the position and velocity of body i using the computed forces.
        double dt = timeStepSize;
        double invm = 1.0 / m[i];
        x[i][0] += dt * v[i][0] + 0.5 * dt * dt * force0[i] * invm;
        x[i][1] += dt * v[i][1] + 0.5 * dt * dt * force1[i] * invm;
        x[i][2] += dt * v[i][2] + 0.5 * dt * dt * force2[i] * invm;
        v[i][0] += dt * force0[i] * invm;
        v[i][1] += dt * force1[i] * invm;
        v[i][2] += dt * force2[i] * invm;

        // Update the maximum velocity and minimum distance.
        double v2 = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (v2 > maxV) {
          maxV = v2;
        }
        for (int j=0; j<neighbors[i].size(); j++) {
          int k = neighbors[i][j];
          double dx = x[k][0] - x[i][0];
          double dy = x[k][1] - x[i][1];
          double dz = x[k][2] - x[i][2];
          double r2 = dx*dx + dy*dy + dz*dz;
          if (r2 < minDx) {
            minDx = r2;
          }
        }

        // Deallocate the arrays.
        delete[] force0;
        delete[] force1;
        delete[] force2;
        }
      
  

/**
 * Check if simulation has been completed.
 */
bool NBodySimulation::hasReachedEnd () {
  return t > tFinal;
}

void NBodySimulation::takeSnapshot () {
  if (t >= tPlot) {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}


void NBodySimulation::openParaviewVideoFile () {
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
    " version=\"0.1\""
    " byte_order=\"LittleEndian\""
    " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

void NBodySimulation::closeParaviewVideoFile () {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulation::printParaviewSnapshot () {
  static int counter = -1;
  counter++;
  std::stringstream filename, filename_nofolder;
  filename << "paraview-output/result-" << counter <<  ".vtp";
  filename_nofolder << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\""
    " NumberOfComponents=\"3\""
    " format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  out.close();

  videoFile << "<DataSet timestep=\"" << counter
            << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
            << "\"/>" << std::endl;
}

void NBodySimulation::printSnapshotSummary () {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
}

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}
