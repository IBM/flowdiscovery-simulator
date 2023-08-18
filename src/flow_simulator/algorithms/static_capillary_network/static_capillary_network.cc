/**
 * \file src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.cc
 * \brief Contains the implementation of \c StaticCapillaryNetworkAlgorithm class methods.
 *
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2016
 *
 * This source file contains the implementation of \c StaticCapillaryNetworkAlgorithm class methods.
 * The \c Initialise() method creates and initialise the \c context_, \c model_, and \c experiment_ objects.
 * The \c BuildGeometricalRepresentation() method is responsible for loading the appropriate file
 *  from disk and creating the geometrical representation required by the algorithm.
 * The \c BuildPhysicalEquations() method builds the \f$ A \vec{x} = \vec{b} \f$ equation.
 * The \c SolvePhysicalEquations() method solves the previous equation for \f$ x \f$ (pressure).
 * The \c CalculateDerivedQuantities() method extracts flow fields from pressure fields.
 * The \c SaveResultsToDisk() method writes all calculated results to files for later visualisation.
 */

#include <array>
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network.h"
#include "src/flow_simulator/algorithms/static_capillary_network/static_capillary_network_context.h"
#include "src/flow_simulator/algorithms/static_capillary_network/models/include_models.h"
#include "src/flow_simulator/algorithms/static_capillary_network/experiments/include_experiments.h"
#include "src/exec_manager/simulation_config.h"

namespace simulator {

std::pair<int, std::string> StaticCapillaryNetworkAlgorithm::Initialise(
  SimulationConfig &simulation_cfg) {
/**
 * Initialise all necessary objects for the simulation, and verify the parameters, in the
 * \c simulation_cfg object, for inconsistencies not restricted by the JSON schema.
 *
 * \param[in]  simulation_cfg  Encapsulates all the input parameter information.
 */
  for (auto fluid : simulation_cfg.fluids_json) {
    if (fluid.GetViscosityBehaviourName() == "constant") {
      if (!fluid.HasProperty("dynamic_viscosity")) {
        return std::pair<int, std::string>(-3,
          "All fluids must have `dynamic_viscosity` property.");
      }
    }
  }

  context_ = std::make_shared<StaticCapillaryNetworkContext>(
    simulation_cfg.voxel_size,
    simulation_cfg.experiment_json.GetFlowAxis(),
    simulation_cfg.shape,
    simulation_cfg.experiment_json.GetBoundaryThickness());

  // Set Model
  std::string model_name = simulation_cfg.algorithm_json.GetModel();
  if (model_name == "poiseuille") {
    model_ = std::make_unique<PoiseuilleModel>(simulation_cfg.folder,
                                               context_,
                                               simulation_cfg.fluids_json,
                                               simulation_cfg.experiment_json);
  } else if (model_name == "hydrostatic") {
    model_ = std::make_unique<HydrostaticModel>(simulation_cfg.folder,
                                                context_,
                                                simulation_cfg.fluids_json,
                                                simulation_cfg.experiment_json);
  } else {
    std::fprintf(stderr, "Please select a valid flow simulation physical model.\n");
    std::exit(-1);
  }

  // Set Experiment
  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();
  double absolute_pressure = simulation_cfg.experiment_json.GetAbsolutePressure();

  if (bc.first == "flow_rate_open") {
    experiment_ = std::make_unique<SinglePhaseFlowRateOpen>(absolute_pressure,
                                                            bc.second,
                                                            context_);
  } else if (bc.first == "flow_rate_closed") {
    experiment_ = std::make_unique<SinglePhaseFlowRateClosed>(absolute_pressure,
                                                              bc.second,
                                                              context_);
  } else if (bc.first == "pressure_gradient") {
    experiment_ = std::make_unique<SinglePhasePressureGradient>(absolute_pressure,
                                                                bc.second,
                                                                context_);
  } else {
    std::fprintf(stderr, "Please select a valid flow simulation experiment.\n");
    std::exit(-1);
  }

  // Checks that the temperature and pressure are within the limits accepted by `ViscosityBehaviour`
  if (!model_->GetPhysics()->IsViscosityBehaviourDesignedForInjectedFluid(absolute_pressure)) {
    LOG(FATAL) << "Injected fluid has pressure or temperature beyond `ViscosityBehaviour` limits";
  }

  return std::pair<int, std::string>(0, "OKay");
}  // end of StaticCapillaryNetworkAlgorithm::Initialise() method



void StaticCapillaryNetworkAlgorithm::BuildGeometricalRepresentation(void) {
/**
 * The \c BuildGeometricalRepresentation() method calls \c ReadNetworkFile() to load
 * \c centerlines.json from disk into member variables of the \c geometry and \c context_ objects.
 */
  // Get geometry and physics members
  auto geometry = model_->GetGeometry();
  auto physics = model_->GetPhysics();

  // Load centerlines file from disk into \c geometry and \c context_ members
  geometry->ReadNetworkFile();

  // Build sparse matrix representing connectivity between nodes
  geometry->BuildConnectivityMatrix();

  // Build diagonal matrix representing link geometrical factors
  geometry->BuildGeometryMatrix(physics->GetFluidViscosity(experiment_->GetAbsolutePressure()));

  // Locate nodes at the edges of the sample for applying boundary conditions
  geometry->LocateBoundaryNodes();

  // Calculate midpoint along flow axis
  context_->midpoint_ = geometry->CalculateMidpointAlongFlowAxis();
}  // end of StaticCapillaryNetworkAlgorithm::BuildGeometricalRepresentation() method



void StaticCapillaryNetworkAlgorithm::BuildPhysicalEquations(void) {
/**
 * The \c BuildPhysicalEquations() builds the matrix equation
 * \f$ A \vec{x} = \vec{b} \f$ of the Capillary Network Algorithm by calling
 * \c StaticCapillaryNetworkAlgorithm::CalculateLeftHandSideMatrix() and
 * \c StaticCapillaryNetworkAlgorithm::CalculateRightHandSideVector().
 */
  // Calculate left hand side matrix
  CalculateLeftHandSideMatrix();

  // Calculate right hand side vector
  CalculateRightHandSideVector();
}  // end of StaticCapillaryNetworkAlgorithm::BuildPhysicalEquations() method



void StaticCapillaryNetworkAlgorithm::CalculateLeftHandSideMatrix(void) {
/**
 * The \c CalculateLeftHandSideMatrix() calculates the \f$ A = C G C^{T} \f$ matrix that goes into
 * the Capillary Network Algorithm matrix equation \f$ A \vec{x} = \vec{b} \f$.
 */
  // Get geometry member
  auto geometry = model_->GetGeometry();

  // Get connectivity and geometry matrices from geometry
  arma::SpMat<double> C = geometry->GetConnectivityMatrix();
  arma::SpMat<double> G = geometry->GetGeometryMatrix();

  // Calculate matrix from C and G
  lhs_matrix_ = C * G * C.t();
}  // end of StaticCapillaryNetworkAlgorithm::CalculateLeftHandSideMatrix() method



void StaticCapillaryNetworkAlgorithm::CalculateRightHandSideVector(void) {
/**
 * The \c CalculateRightHandSideVector() calculates the \f$ \vec{b} \f$ vector that goes into the
 * Capillary Network Algorithm matrix equation \f$ A \vec{x} = \vec{b} \f$,
 */
  // Get physics member
  auto physics = model_->GetPhysics();

  // Calculate right hand side vector from pressure boundary conditions along flow direction
  arma::vec rhs_vector_term = physics->CalculateRightHandSideVectorTerm();
  rhs_vector_ = experiment_->ApplyBoundaryConditions(lhs_matrix_, rhs_vector_term);
}  // end of StaticCapillaryNetworkAlgorithm::CalculateRightHandSideVector() method



void StaticCapillaryNetworkAlgorithm::SolvePhysicalEquations(void) {
/**
 * The \c SolvePhysicalEquations() solves the matrix equation
 * \f$ A \vec{x} = \vec{b} \f$ of the Capillary Network Algorithm.
 *
 */
  // Configure SuperLU solver
  arma::superlu_opts solver_config;
  solver_config.equilibrate = true;

  // Solve matrix equation
  pressures_ = arma::spsolve(lhs_matrix_, rhs_vector_, "superlu", solver_config);
}  // end of StaticCapillaryNetworkAlgorithm::SolvePhysicalEquations() method



void StaticCapillaryNetworkAlgorithm::CalculateDerivedQuantities(void) {
/**
 * After computing the pressures, a number of derived quantities can be calculated, such as
 * flow rate, flow speed and permeability.
 *
 * The local flow speeds are calculated by \c StaticCapillaryExperimentBase::CalculateFlowSpeed,
 * the permeability is calculated by \c StaticCapillaryExperimentBase::CalculatePermeability ,
 * while the flow rate is calculated differently according to each model that implements
 * \c StaticCapillaryGeometryBase::CalculateFlowRate.
 */
  // Get physics member
  auto physics = model_->GetPhysics();

  // Calculate flow rate from node pressures
  flow_rate_ = physics->CalculateFlowRate(pressures_);

  // Calculate flow speeds from flow rates
  flow_speed_ = physics->CalculateFlowSpeed(flow_rate_);

  // Calculate permeability from flow rates
  permeability_ = physics->CalculatePermeability(flow_rate_, pressures_);
}  // end of StaticCapillaryNetworkAlgorithm::CalculateDerivedQuantities() method



void StaticCapillaryNetworkAlgorithm::SaveResultsToDisk(void) {
/**
 * The \c SaveResultsToDisk() method saves the resulting pressures, flow rates, flow speeds and
 * permeability to disk.
 */
  // Get geometry member
  auto geometry = model_->GetGeometry();

  // Get node indexes from geometry
  arma::Row<IndexType> node_index = context_->linked_nodes_.row(0);

  // Define auxiliary variables
  arma::uword u = context_->flow_axis_;
  std::string file = geometry->folder_ + "/static_results.h5";
  arma::vec permeability = { permeability_ };
  std::array<const char*, 3> axis = { {"x", "y", "z"} };

  // Save all resulting quantities to disk
  pressures_.save(arma::hdf5_name(file, "pressures_" + std::string(axis[u]),
                                  arma::hdf5_opts::append));
  flow_rate_.save(arma::hdf5_name(file, "flow_rate_" + std::string(axis[u]),
                                        arma::hdf5_opts::append));
  flow_speed_.save(arma::hdf5_name(file, "flow_speed_" + std::string(axis[u]),
                                         arma::hdf5_opts::append));
  permeability.save(arma::hdf5_name(file, "permeability_" + std::string(axis[u]),
                                    arma::hdf5_opts::append));
}  // end of StaticCapillaryNetworkAlgorithm::SaveResultsToDisk() method



void StaticCapillaryNetworkAlgorithm::Finalize(void) {
}  // end of StaticCapillaryNetworkAlgorithm::Finalize() method

}  // namespace simulator
