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

  const double cutoff = 2.0
  const double skin = 0.2
  const double cutoff2 = (cutoff + skin) * (cutoff + skin);
  const double C = (1/100)/NumberOfBodies
  std::vector<std::vector<int>> neighbors(NumberOfBodies);
  
  for (int i = 0; i < NumberOfBodies; i++) {
    for (int j = i + 1; j < NumberOfBodies; j++) {
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 <= C * (mass[i]+mass[j])) {
        x[i][0] = (mass[i]*x[i][0] + mass[j]*x[j][0]) / (mass[i] + mass[j]);
        x[i][1] = (mass[i]*x[i][1] + mass[j]*x[j][1]) / (mass[i] + mass[j]);
        x[i][2] = (mass[i]*x[i][2] + mass[j]*x[j][2]) / (mass[i] + mass[j]);
        v[i][0] = (mass[i]*v[i][0] + mass[j]*v[j][0]) / (mass[i] + mass[j]);
        v[i][1] = (mass[i]*v[i][1] + mass[j]*v[j][1]) / (mass[i] + mass[j]);
        v[i][2] = (mass[i]*v[i][2] + mass[j]*v[j][2]) / (mass[i] + mass[j]);
        mass[i] += mass[j];
        x[j][0] = x[NumberOfBodies-1][0];
        x[j][1] = x[NumberOfBodies-1][1];
        x[j][2] = x[NumberOfBodies-1][2];
        v[j][0] = v[NumberOfBodies-1][0];
        v[j][1] = v[NumberOfBodies-1][1];
        v[j][2] = v[NumberOfBodies-1][2];
        mass[j] = mass[NumberOfBodies-1];
        NumberOfBodies--;
        j--;
        continue
      }
      if (r2 <= cutoff2) {
        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
      }
    }
  }

  // Calculate the forces acting on each body.
  for (int i=0; i<NumberOfBodies; i++) {
    force0[i] = 0.0;
    force1[i] = 0.0;
    force2[i] = 0.0;

    for (int j : neighbors[i]) {
      // Calculate the force acting on body i due to body j.
      force0[i] += force_calculation(i,j,0);
      force1[i] += force_calculation(i,j,1);
      force2[i] += force_calculation(i,j,2);
    }

    // Half-step update of the velocity.
    v[i][0] += 0.5 * timeStepSize * force0[i] / mass[i];
    v[i][1] += 0.5 * timeStepSize * force1[i] / mass[i];
    v[i][2] += 0.5 * timeStepSize * force2[i] / mass[i];

    // Full-step update of the position.
    x[i][0] += timeStepSize * v[i][0];
    x[i][1] += timeStepSize * v[i][1];
    x[i][2] += timeStepSize * v[i][2];

    // Calculate the forces acting on the updated position.
    double force0_new = 0.0;
    double force1_new = 0.0;
    double force2_new = 0.0;

    for (int j : neighbors[i]) {
      if (i == j) continue; // Skip self-interaction.

      // Calculate the force acting on body i due to body j.
      force0_new += force_calculation(i,j,0);
      force1_new += force_calculation(i,j,1);
      force2_new += force_calculation(i,j,2);
    }

    // Half-step update of the velocity.
    v[i][0] += 0.5 * timeStepSize * force0_new / mass[i];
    v[i][1] += 0.5 * timeStepSize * force1_new / mass[i];
    v[i][2] += 0.5 * timeStepSize * force2_new / mass[i];

    for ()
    // Update the max velocity.
    double vel = std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] );
    maxV = std::max(maxV, vel);

  }

  t += timeStepSize;

  // Free the memory allocated for the forces.
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
