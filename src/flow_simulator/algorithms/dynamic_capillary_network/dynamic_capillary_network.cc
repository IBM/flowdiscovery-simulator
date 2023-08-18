/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.cc
 * \brief Contains the implementation of \c DynamicCapillaryNetwork class methods.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This source file contains the implementation of \c DynamicCapillaryNetworkAlgorithm class methods.
 * The \c Initialise() method creates and initialise the \c context_, \c model_, and \c experiment_ objects.
 * The \c BuildGeometricalRepresentation() method is responsible for loading the appropriate file
 *  from disk and creating the geometrical representation required by the algorithm.
 * The \c SolvePhysicalEquations() method solves the system of DAEs required for the dynamic simulation through time.
 * The \c CalculateDerivedQuantities() method extracts derived quantities from the simulation results.
 * The \c SaveResultsToDisk() method writes all calculated results to files for later visualisation.
 * The \c SetNumberOfEquations() method calculates the number of equations in the DAE system.
 * The \c SolveDAE() method solves the DAE system using the SUNDIALS IDA solver.
 * The \c StoreStepResults() method stores the results of each time step in of the solver.
 * The \c StoreResumeInfo() method stores the state of variables (resume).
 * The \c LoadResumeInfo() method loads the state of variables (resume).
 * The \c FinishSolveDAE() method clear memory allocation when the solver finishes.
 * The \c InitialiseIDASolver() method initialises the variables for the \c IDASolve.
 * The \c InitialiseIDAState() method initialises the state of the system before the next iteration in time.
 * The \c ComputeInitialValues() method compute the values of the initial state of the simulation based on input.
 * The \c ComputeDifferentialEquationForAllCapillaries() method compute the values of the differential equations of the simulation.
 * The \c ComputeDifferentialEquationForSelectedCapillaries() method compute the values of the differential for selected equations of the simulation.
 * The \c ComputeDifferentialEquation() method compute the values of one differential equation of the simulation.
 * The \c SetIDValues() method set the values to identify which equations are differential or algebraic.
 * The \c UpdateSundialsVariables() method update the values of x and dx before the next time iteration.
 * The \c UpdateInterfacesInXArray() method update the values of x before the next time iteration.
 * The \c StopCriteriaMet() method check if any stop criteria is met.
 * The \c CheckIDAState() method check if there is any state change, so the new state need to be updated and initialised.
 * The \c StoreInteriorPressures() method stores the current pressures.
 * The \c RevertPressures() method Revert the pressures to a previous state.
 * The \c StoreState() method stores the current IDA Sundials state.
 * The \c RevertState() method Revert the pressures to a previous IDA Sundials state.
 * The \c RemoveBubblesFromOrigin() method removes pairs of very close interfaces that are beyond the limits of capillaries with roots.
 * The \c PlugCapillariesInterfacesJumpedThroughInletNodes() method plug capillaries witch interfaces jumped through inlet-nodes
 */

#include <glog/logging.h>
#include <ida/ida.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <ida/ida_direct.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/experiments/include_experiments.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/linear_molecular_kinetics/linear_molecular_kinetics_model.h"
#include "src/exec_manager/simulation_config.h"

namespace simulator {

std::ostream& operator<<(std::ostream& os, const Side& a) {
  switch (a) {
  case Side::none:
    return os << "none";
    break;
  case Side::source:
    return os << "source";
    break;
  case Side::target:
    return os << "target";
    break;
  default:
    return os << "<error:undefined_value>";
    break;
  }
}

std::pair<int, std::string> DynamicCapillaryNetworkAlgorithm::Initialise(
  SimulationConfig &simulation_cfg) {
/**
 * Initialise all necessary objects for the simulation, and verify the parameters, in the
 * \c Settings object, for inconsistencies not restricted by the JSON schema.
 *
 * \param[in]  simulation_cfg  Encapsulate all the input parameter information.
 */
  if (simulation_cfg.fluids_json.size() < 2)
    return std::pair<int, std::string>(-2,
        "Parameter `fluids` must have, at least, two fluids.");

  for (auto fluid : simulation_cfg.fluids_json) {
    if (fluid.GetViscosityBehaviourName() == "constant") {
      if (!fluid.HasProperty("dynamic_viscosity")) {
        return std::pair<int, std::string>(-3,
          "All fluids with constant viscosity behaviour must have `dynamic_viscosity` property.");
      }
    }
  }

  initial_time_ = simulation_cfg.algorithm_json.GetInitialTime();
  final_time_ = simulation_cfg.algorithm_json.GetFinalTime();
  resume_ = simulation_cfg.algorithm_json.GetResume();

  if (final_time_ <= initial_time_)
    return std::pair<int, std::string>(-1,
        "Parameter `final_time` must be greater than `initial_time`.");

  time_step_size_ = simulation_cfg.algorithm_json.GetTimeStepSize();
  relative_tolerance_ = simulation_cfg.algorithm_json.GetRelativeTolerance();
  absolute_link_tolerance_ = simulation_cfg.algorithm_json.GetAbsoluteLinkTolerance();
  absolute_node_tolerance_ = simulation_cfg.algorithm_json.GetAbsoluteNodeTolerance();

  // Set Context
  context_ = std::make_shared<DynamicCapillaryNetworkContext>(
    simulation_cfg.voxel_size,
    simulation_cfg.experiment_json.GetFlowAxis(),
    simulation_cfg.shape,
    simulation_cfg.experiment_json.GetBoundaryThickness());

  // Set Interface
  interface_ = std::make_shared<CapillaryInterface>(context_);

  // Set Model
  if (simulation_cfg.algorithm_json.GetModel() == "linear_molecular_kinetics") {
    model_ = std::make_unique<LinearMolecularKineticsModel>(simulation_cfg.folder,
                                                            context_,
                                                            simulation_cfg.fluids_json,
                                                            simulation_cfg.wettability_json,
                                                            simulation_cfg.fluid_interface_json,
                                                            simulation_cfg.experiment_json);
  }

  // Set Experiment
  std::pair<std::string, double> bc = simulation_cfg.experiment_json.GetBoundaryCondition();
  double absolute_pressure = simulation_cfg.experiment_json.GetAbsolutePressure();

  if (simulation_cfg.experiment_json.GetBoundaryCondition().first == "pressure_gradient_open") {
    experiment_ = std::make_unique<TwoPhasePressureGradientOpen>(absolute_pressure,
                                                                 bc.second,
                                                                 context_);
  } else if (simulation_cfg.experiment_json.GetBoundaryCondition()
      .first == "pressure_gradient_closed") {
    experiment_ = std::make_unique<TwoPhasePressureGradientClosed>(absolute_pressure,
                                                                   bc.second,
                                                                   context_);

  } else {
    std::fprintf(stderr, "Please select a valid flow simulation experiment.\n");
    std::exit(-1);
  }

  // Logging some inputted parameters

  LOG(INFO) << std::endl
      << "Initial time: " << initial_time_ << std::endl
      << "Final time: " << final_time_ << std::endl
      << "Time step size" << time_step_size_ << std::endl
      << "Relative tolerance: " << relative_tolerance_ << std::endl
      << "Absolute Link tolerance: " << absolute_link_tolerance_ << std::endl
      << "Absolute Node tolerance: " << absolute_node_tolerance_;


  // Checks that the temperature and pressure are within the limits accepted by `ViscosityBehaviour`
  DynamicCapillaryPhysicsBase *physics = model_->GetPhysics().get();
  double voxel_size = context_->voxel_size_;
  double number_of_voxel = static_cast<double>(context_->shape_(context_->flow_axis_));
  double boundary_condition_value = bc.second;
  double p_min = absolute_pressure - (voxel_size * number_of_voxel / 2)
               * boundary_condition_value;
  double p_max = absolute_pressure + (voxel_size * number_of_voxel / 2)
               * boundary_condition_value;

  if (!physics->IsViscosityBehaviourDesignedForResidentFluid(p_min, p_max)) {
    LOG(INFO) << "Pressure: (" << p_min << ", " << p_max << "), check Temp and other parameters";
    LOG(FATAL) << "Resident fluid has pressure or temperature beyond `ViscosityBehaviour` limits";
  }

  if (!physics->IsViscosityBehaviourDesignedForInjectedFluid(p_min, p_max)) {
    LOG(INFO) << "Pressure: (" << p_min << ", " << p_max << "), check Temp and other parameters";
    LOG(FATAL) << "Injected fluid has pressure or temperature beyond `ViscosityBehaviour` limits";
  }

  // Logging important informations

  for (auto wet_props : simulation_cfg.wettability_json.GetProperties()) {
    LOG(INFO) << wet_props.first << ": " << wet_props.second;
  }

  for (auto fluid_int_props : simulation_cfg.fluid_interface_json.GetProperties()) {
    LOG(INFO) << fluid_int_props.first << ": " << fluid_int_props.second;
  }

  LOG(INFO) << "Properties from fluid " << simulation_cfg.fluids_json[0].GetName() << ":";
  for (auto res_fluid_props : simulation_cfg.fluids_json[0].GetProperties()) {
    LOG(INFO) << res_fluid_props.first << ": " << res_fluid_props.second;
  }

  LOG(INFO) << "Properties from fluid " << simulation_cfg.fluids_json[1].GetName() << ":";
  for (auto inj_fluid_props : simulation_cfg.fluids_json[1].GetProperties()) {
    LOG(INFO) << inj_fluid_props.first << ": " << inj_fluid_props.second;
  }

  return std::pair<int, std::string>(0, "OKay");
}  // end of DynamicCapillaryNetworkAlgorithm::Initialise() method

int DynamicCapillaryNetworkAlgorithm::DetectJumps(const double t,
                                                  const N_Vector x,
                                                  const N_Vector dx,
                                                  double *gout,
                                                  void *user_data) {
/**
 * The jump of an interface from a capillary to its neighbours is defined by the
 * by the reach of a position 0 or 1 in the capillary. This functions is called by
 * \c IDASolve to find the roots of f(x) as defined below.
 *
 * For each capillary with an interface
 *
 * \f[
 *    f(x) = [x - (eff_i - min_i)] [x - (eff_i + (1 - max_i))] \, .
 * \f]
 *
 * Refer to SUNDIALS documentation for additional info.
 *
 * \param[in]      t           Current solution time of the DAE.
 * \param[in]      x           SUNDIALS N_Vector containing the capillaries effective interface position.
 * \param[in]      dx          SUNDIALS N_Vector containing the capillaries interface movement speed.
 * \param[out]     gout        Double array where the root function current value is stored
 * \param[in]      user_data   DynamicCapillaryNetworkAlgorithm object
 * \retval         error_code  Integer containing zero for success or non zero for error.
 */
  // The line below is required to avoid unused variable warning messages. This function is called
  // by IDASolve from SUNDIALS.
  static_cast<void>(t);  // t is not used for this model

  // Convert SUNDIALS's N_Vectors into Armadillo's arma::vec
  const arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);
  const arma::vec dx_vec(NV_DATA_S(dx), static_cast<arma::uword>(NV_LENGTH_S(dx)), false, true);

  // Vectorise output C array
  const auto interface = static_cast<DynamicCapillaryNetworkAlgorithm*>(user_data)->GetInterface();
  const auto capillaries = static_cast<DynamicCapillaryNetworkAlgorithm*>(user_data)
      ->caps_with_interfaces_not_plugged_;
  const auto caps_residual_limits = static_cast<DynamicCapillaryNetworkAlgorithm*>(user_data)
      ->caps_residual_limits_;

  // Convert SUNDIALS's double array (double*) into Armadillo's arma::vec
  arma::vec gout_vec(gout, caps_residual_limits.n_elem, false, true);

  // TODO(ralves): parallelisation is possible
  for (IndexType i = 0U; i < capillaries.n_elem; ++i) {
    const arma::uword cap_i = capillaries(i);

    // Compute the root values that determines if there will be an interface jump
    const IndexType index = 2U * i;
    gout_vec(index) = x_vec(cap_i) - caps_residual_limits(index);
    gout_vec(index + 1) = x_vec(cap_i) - caps_residual_limits(index + 1);
  }

  return 0;
}  // end of DynamicCapillaryNetworkAlgorithm::DetectJumps() method



int DynamicCapillaryNetworkAlgorithm::Residual(const double t,
                                               const N_Vector x,
                                               const N_Vector dx,
                                               N_Vector resval,
                                               void *user_data) {
/**
 * The residual function is called by \c IDASolver. This is where the physical equations that govern
 * the fluid movement and interaction are applied. Refer to SUNDIALS documentation for more info.
 *
 * The residual of the physical equation being solved is computed as follows:
 *
 * For \f$ 0 \le i \lt M \f$, the residual represents the difference between the calculated and target speeds
 *
 * \f[
 *    r_i = F(x_i) - \dot{x_i} \, .
 * \f]
 *
 * For \f$ M \le i \lt E \f$, the residual represents the mass continuity equation at each interior node \f$ k \f$
 *
 * \f[
 *    r_i = \sum\limits_{j \rightarrow k} \alpha_j \dot{x_j} R_j^2 \, .
 * \f]
 *
 * Where:
 *  - \f$ M \f$ is the number of capillaries.
 *  - \f$ N_{\text{int}} \f$ is the number of interior nodes.
 *  - \f$ E = N_{\text{int}} + M \f$ is the number of equations.
 *  - \f$ r_i \f$ is the residual of the \f$ i \f$-th equation in the DAE system.
 *  - \f$ x_i \f$ is the effective interface position in capillary \f$ i \f$.
 *  - \f$ F(x_i) \f$ is the effective interface speed in capillary \f$ i \f$ as calculated by ComputeDifferentialEquation().
 *  - \f$ \dot{x_i} \f$ is the effective interface speed of the capillary \f$ i \f$ as obtained by \c IDASolver.
 *  - \f$ R_i^2 \f$ is the squared radius of capillary \f$ i \f$.
 *  - \f$ j \rightarrow k \f$ is the set of capillaries \f$ j \f$ connected to interior node \f$ k \f$ where \f$ k = i - M \f$.
 *  - \f$ \alpha_i \f$ defines the direction of the flow in the capillary \f$ i \f$.
 *
 * \param[in]     t           Current solution time of the DAE.
 * \param[in]     x           SUNDIALS N_Vector containing the capillaries effective interface position.
 * \param[in]     dx          SUNDIALS N_Vector containing the capillaries interface movement speed.
 * \param[out]    resval      SUNDIALS N_Vector containing the residual of the physical equations of the DAE
 * \param[in]     user_data   DynamicCapillaryNetworkAlgorithm object
 * \retval        error_code  Integer containing SUNDIALS error code.
 */
  // The line below is required to avoid unused variable warning messages. This function is called
  // by IDASolve from SUNDIALS.
  static_cast<void>(t);  // t is not used for this model

  //  First step is to populate pressures_ with the values from x
  const auto data = static_cast<DynamicCapillaryNetworkAlgorithm*>(user_data);
  const auto interface = data->GetInterface();
  const auto context = data->GetContext();
  const auto geometry = data->model_->GetGeometry();

  // Allocate local variables and Armadillo vectors
  const auto pi = arma::datum::pi;
  const IndexType n_links = context->number_of_links_;
  const IndexType n_equations = data->GetNumberOfEquations();
  const arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);
  const arma::vec dx_vec(NV_DATA_S(dx), static_cast<arma::uword>(NV_LENGTH_S(dx)), false, true);
  arma::vec residual(NV_DATA_S(resval), static_cast<arma::uword>(NV_LENGTH_S(resval)), false, true);

  // Copy the pressure in the interiors to the complete vector of pressures
  data->SetPressures(context->interior_nodes_, x_vec.tail(context->interior_nodes_.n_elem));

  // For each differential equation, check if there is an interface and use the according equation
  data->ComputeDifferentialEquationForAllCapillaries(residual, x_vec);
  residual.head(context->number_of_links_) -= dx_vec.head(context->number_of_links_);
  residual.tail(context->interior_nodes_.n_elem).zeros();

  #pragma omp parallel for
  for (IndexType i = n_links; i < n_equations; ++i) {
    // Find the capillaries connected to \c interior_node, either as a source or as a target
    const IndexType interior_node = context->interior_nodes_(i - n_links);
    const arma::uvec &capillaries_source = geometry->GetCapillariesWhoseSourceIs(interior_node);
    const arma::uvec &capillaries_target = geometry->GetCapillariesWhoseTargetIs(interior_node);

    // Accumulate negative flow contributions from source capillaries connected to interior_node
    for (const auto &cap : capillaries_source) {
      residual(i) = residual(i) - (pi *
                                   context->link_squared_radius_(cap) *
                                   context->link_length_(cap) *
                                   dx_vec(cap));
    }

    // Accumulate positive flow contributions from target capillaries connected to interior_node
    for (const auto &cap : capillaries_target) {
      residual(i) = residual(i) + (pi *
                                   context->link_squared_radius_(cap) *
                                   context->link_length_(cap) *
                                   dx_vec(cap));
    }
  }

  return 0;
}  // end of DynamicCapillaryNetworkAlgorithm::Residual() method



void DynamicCapillaryNetworkAlgorithm::BuildGeometricalRepresentation(void) {
/**
 * The \c BuildGeometricalRepresentation() method calls \c ReadNetworkFile() to load
 * \c centerlines.json from disk into member variables of the \c geometry and \c context_ objects.
 */
  // Get geometry and physics members
  auto geometry = model_->GetGeometry();

  // Load centerlines file from disk into \c geometry and \c context_ members
  geometry->ReadNetworkFile();

  // Locate nodes at the edges of the sample for applying boundary conditions
  geometry->LocateBoundaryNodes();

  // Calculate midpoint along flow axis
  geometry->CalculateMidpointAlongFlowAxis();

  // Identify the list of capillaries connected to each node
  geometry->LocateCapillariesConnectedToNodes();
}  // end of DynamicCapillaryNetworkAlgorithm::BuildGeometricalRepresentation() method



void DynamicCapillaryNetworkAlgorithm::BuildPhysicalEquations(void) {
  // Set sizes of interface-related variables
  interface_->InitialiseCapillaryInterface();

  // Apply boundary conditions to nodes
  pressures_ = experiment_->ApplyBoundaryConditions(interface_);

  // Set number of equations private member
  SetNumberOfEquations();
}  // end of DynamicCapillaryNetworkAlgorithm::BuildPhysicalEquations() method



void DynamicCapillaryNetworkAlgorithm::SolvePhysicalEquations(void) {
/**
 * The \c SolvePhysicalEquations method solves the DAE system as a function of time using SUNDIALS IDA.
 * The method begins by initialising the SUNDIALS IDA solver based on the current simulation.
 * After initialising, it loops into initialising the new DAE state, adjusting the interface
 * information, changing the equations. Then the innermost cycle calls IDASolve increasing the time by
 * \c time_step_size. This cycle is broken whenever there is a change of state in the system, due
 * an interface jump, causing changes in DAE system equations.
 */
  // Create folder to store snapshot files
  CreateSnapshotFolder();

  // Create SUNDIALS variables
  int retval;
  int retvalr;
  N_Vector x, dx;
  SUNMatrix A = NULL;
  SUNLinearSolver LS = NULL;

  // Initialise containers for state variables
  previous_interior_pressures_.zeros(context_->interior_nodes_.n_elem);
  previous_x_.zeros(number_of_equations_);
  previous_dx_.zeros(number_of_equations_);
  previous_A_.zeros(number_of_equations_ * number_of_equations_);

  IndexType output_time_index = 0U;  // Stores the index of the current output time step

  InitialiseIDASolver(&mem_, &x, &dx, &A, &LS, initial_time_);
  const arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);

  // initial data
  current_time_ = initial_time_;
  double next_target_time = current_time_ + time_step_size_;

  if (resume_) {
    LoadResumeInfo(x, dx, A, &output_time_index, &next_target_time);
  }

  InitialiseIDAState(mem_, x, dx, A, LS,
                     current_time_, next_target_time, output_time_index);

  while (1) {
    if (current_time_ >= final_time_) break;

    retval = IDASolve(mem_, next_target_time, &current_time_, x, dx, IDA_NORMAL);
    if (check_retval(&retval, "IDASolve", 1)) exit(1);
    interface_->ComputeDeltas(x_vec);
    interface_->UpdateInterfaceInformation();

    if (retval == IDA_ROOT_RETURN) {  // Process Jumps
      retvalr = IDAGetRootInfo(mem_, roots_found_.memptr());
      check_retval(&retvalr, "IDAGetRootInfo", 1);
      InitialiseIDAState(mem_, x, dx, A, LS,
                        initial_time_, next_target_time, output_time_index);
    } else if (retval == IDA_SUCCESS) {  // Process Time Out
      StoreStepResults(x, dx, static_cast<IndexType>(output_time_index));
      next_target_time += RCONST(time_step_size_);
      output_time_index++;
    } else {
      LOG(FATAL) << "Error return value from SUNDIALS inexpected!";
    }

    // Statistics

    // Calculate number of Capillaries Plugged for logs purposes
    IndexType ncp = 0;
    const arma::Col<IndexType> capillaries = interface_->GetIndexOfCapillariesWithInterfaces();
    for (const auto &capillary : capillaries) {
      if (interface_->IsCapillaryPlugged(capillary)) {
        ncp++;
      }
    }
    LOG(INFO) << "### Begin ###";
    LOG(INFO) << "next_target_time: " << next_target_time;
    LOG(INFO) << "current_time_: " << current_time_;
    LOG(INFO) << "Number of capillaries with interfaces: "
        << interface_->GetNumberOfCapillariesWithInterfaces();
    LOG(INFO) << "Number of interfaces: " << interface_->GetNumberOfInterfaces();
    LOG(INFO) << "Number of capillaries plugged: " << ncp;
    LOG(INFO) << "@@@ End @@@";
  }

  // Compressing the resume file
  std::string resume_file = GetFolder() + "/resume.h5";
  if (std::filesystem::exists(resume_file)) {
    std::string compress_command = "tar -czf " + resume_file + ".tgz -C" + GetFolder() +
                                   " resume.h5";
    system(compress_command.c_str());
    std::filesystem::remove(resume_file);
  }
}  // end of DynamicCapillaryNetworkAlgorithm::SolvePhysicalEquations() method



void DynamicCapillaryNetworkAlgorithm::CalculateDerivedQuantities(void) {
}  // end of DynamicCapillaryNetworkAlgorithm::CalculateDerivedQuantities() method



void DynamicCapillaryNetworkAlgorithm::SaveResultsToDisk(void) {
}  // end of DynamicCapillaryNetworkAlgorithm::SaveResultsToDisk() method



void DynamicCapillaryNetworkAlgorithm::Finalize(void) {
  int retval;
  int64_t nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  retval = IDAGetNumSteps(mem_, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(mem_, &nre);
  check_retval(&retval, "IDAGetNumResEvals", 1);
  retval = IDADlsGetNumJacEvals(mem_, &nje);
  check_retval(&retval, "IDADlsGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(mem_, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(mem_, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem_, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDAGetNumResEvals(mem_, &nreLS);
  check_retval(&retval, "IDADlsGetNumResEvals", 1);
  retval = IDAGetNumGEvals(mem_, &nge);
  check_retval(&retval, "IDAGetNumGEvals", 1);

  LOG(INFO) << "\nFinal Run Statistics: \n\n"
    << "Number of steps                    = " << nst << std::endl
    << "Number of residual evaluations     = " << nre + nreLS << std::endl
    << "Number of Jacobian evaluations     = " << nje << std::endl
    << "Number of nonlinear iterations     = " << nni << std::endl
    << "Number of error test failures      = " << netf << std::endl
    << "Number of nonlinear conv. failures = " << ncfn << std::endl
    << "Number of root fn. evaluations     = " << nge << std::endl;
}  // end of DynamicCapillaryNetworkAlgorithm::Finalize() method



void DynamicCapillaryNetworkAlgorithm::CreateSnapshotFolder(void) {
/**
 * The \c SetupResultsStorage method creates the folder where the results will be stored and
 * stores the folder name to be used in the results storage methods.
 */

  // Defines the folder where results for this simulation will be saved
  std::string results_folder = GetFolder() + snapshots_folder_;

  // Creates the folder if it does not exists
  if (!std::filesystem::is_directory(results_folder)) {
    std::filesystem::create_directory(results_folder);
  } else {
    // Not resuming = starting a new analysis
    if (!resume_) {
      // Delete all files in the snapshots folder
      for (auto const& file : std::filesystem::directory_iterator(results_folder)) {
        std::filesystem::remove(file);
      }
    }
  }
}  // end of DynamicCapillaryNetworkAlgorithm::SetupResultsStorage() method



void DynamicCapillaryNetworkAlgorithm::StoreStepResults(const N_Vector x,
                                                        const N_Vector dx,
                                                        IndexType output_time_index) {
/**
 * The \c StoreStepResults method stores the DAE system state at every output time step
 * defined by \c time_step_size. Each output time step is saved in a different HDF5 file, containing
 * the following data sets:
 * - \c pressures
 * - \c permeability
 * - \c current_time
 * - \c fluid_saturation
 * - \c fluid_at_source
 * - \c interface_offsets
 * - \c interface_positions
 *
 * \param[in] x                   SUNDIALS N_Vector containing effective interface positions.
 * \param[in] dx                  SUNDIALS N_Vector containing effective interface speed.
 * \param[in] output_time_index   Unsigned int that holds the output time index.
 */
  /// Creates the filename based on the output time index
  std::string temp_str = std::to_string(output_time_index);
  temp_str.insert(temp_str.begin(), 6 - temp_str.size(), '0');
  std::string file = GetFolder() + snapshots_folder_ + "/simres_output_time_" + temp_str + ".h5";

  // Writes the current time of simulation
  const arma::vec current_time = {current_time_};

  // Copy the pressure in the interiors to the complete vector of pressures
  const arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);
  const arma::vec dx_vec(NV_DATA_S(dx), static_cast<arma::uword>(NV_LENGTH_S(dx)), false, true);
  SetPressures(context_->interior_nodes_, x_vec.tail(context_->interior_nodes_.n_elem));

  // Calculates all interface positions and index offsets for accessing them
  arma::vec interface_positions(interface_->GetNumberOfInterfaces());
  arma::uvec interface_offsets(context_->number_of_links_ + 1U, arma::fill::zeros);
  interface_->ComputeInterfacePositionsAndOffsets(interface_positions, interface_offsets);

  // Calculate fluid saturations
  const arma::Col<IndexType> fluid_at_source = interface_->GetFluidAtSource();
  const arma::mat capillary_saturation = interface_->ComputeCapillarySaturation();

  // Calculate fluid permeability
  auto physics = model_->GetPhysics();
  const arma::vec flow_speed = physics->CalculateFlowSpeed(dx_vec);
  const arma::vec flow_rate = physics->CalculateFlowRate(flow_speed);
  const arma::vec permeability = physics->CalculatePhasePermeability(flow_rate, pressures_);
  const arma::vec fluid_saturation = physics->CalculateFluidSaturation(capillary_saturation);

  // Write network variables to disk
  pressures_.save(arma::hdf5_name(file, "pressure", arma::hdf5_opts::append));
  flow_rate.save(arma::hdf5_name(file, "flow_rate", arma::hdf5_opts::append));
  flow_speed.save(arma::hdf5_name(file, "flow_speed", arma::hdf5_opts::append));
  permeability.save(arma::hdf5_name(file, "permeability", arma::hdf5_opts::append));
  current_time.save(arma::hdf5_name(file, "current_time", arma::hdf5_opts::append));
  fluid_at_source.save(arma::hdf5_name(file, "fluid_at_source", arma::hdf5_opts::append));
  fluid_saturation.save(arma::hdf5_name(file, "fluid_saturation", arma::hdf5_opts::append));
  interface_offsets.save(arma::hdf5_name(file, "interface_offsets", arma::hdf5_opts::append));
  interface_positions.save(arma::hdf5_name(file, "interface_positions", arma::hdf5_opts::append));

  LOG(INFO) << "Permeability[" << permeability(0) << ", " << permeability(1)
            << "] Fluid Saturation[" << fluid_saturation(0) << ", " << fluid_saturation(1) << "]";
}  // end of DynamicCapillaryNetworkAlgorithm::StoreStepResults() method


void DynamicCapillaryNetworkAlgorithm::StoreResumeInfo(const N_Vector x,
                                                        const N_Vector dx,
                                                        const SUNMatrix A,
                                                        IndexType output_time_index,
                                                        double next_target_time) {
/**
 * The \c StoreResumeInfo method stores the system parameters every time a root is found. These
 * can be later loaded to resume the analysis. The outputs are saved in a HDF5 file
 * (resume.h5), containing the following data sets:
 * - \c current_time
 * - \c next_target_time
 * - \c output_time_index
 * - \c resume_next_target_time
 * - \c last_jump_time
 * - \c x
 * - \c dx
 * - \c A
 * - \c previous_x
 * - \c previous_dx
 * - \c previous_A
 * - \c all_interfaces
 * - \c deltas
 * - \c capillaries_with_interfaces
 * - \c plugs
 * - \c last_plug_info
 * - \c last_plug_inlet_info
 * - \c last_jump_info_size
 * - \c last_jump_info_data
 * - \c caps_with_interfaces_not_plugged
 * - \c caps_with_residual_limits
 * - \c roots_found
 */
  /// Resume file
  std::string resume_file = GetFolder() + "/resume.h5";

  // Gets time info of simulation
  const arma::vec time_info = {current_time_, next_target_time, last_jump_time_};

  // Gets the current output time index of simulation
  const arma::uvec current_time_step = {output_time_index};

  // Get references for the x and dx vectors
  const arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);
  const arma::vec dx_vec(NV_DATA_S(dx), static_cast<arma::uword>(NV_LENGTH_S(dx)), false, true);

  // Get reference for the elements of the current A matrix of the simulation
  const arma::vec A_vec(NV_DATA_S(A), static_cast<arma::uword>(SM_LDATA_D(A)), false, true);

  // Get interfaces information
  const arma::field<arma::Col<double>> all_interfaces = interface_->GetInterfacePositions();
  const arma::Col<double> deltas = interface_->GetDeltas();
  const arma::Col<IndexType> capillaries_with_interfaces =
                                            interface_->GetIndexOfCapillariesWithInterfaces();
  const arma::Col<IndexType> plugs = interface_->GetPlugs();
  const arma::Col<IndexType> last_plug_info = interface_->GetLastPlugInfo();
  const arma::Col<IndexType> last_plug_inlet_info = interface_->GetLastPlugInletInfo();
  InterfaceJumpEventInfo last_jump_info = interface_->GetLastJumpInfo();

  // Calculate fluid saturations
  const arma::Col<IndexType> fluid_at_source = interface_->GetFluidAtSource();

  // Write variables required to resume to disk
  time_info.save(arma::hdf5_name(resume_file, "time_info"));
  current_time_step.save(arma::hdf5_name(resume_file,
                                         "output_time_index",
                                         arma::hdf5_opts::append));

  x_vec.save(arma::hdf5_name(resume_file, "x_vec", arma::hdf5_opts::append));
  dx_vec.save(arma::hdf5_name(resume_file, "dx_vec", arma::hdf5_opts::append));
  A_vec.save(arma::hdf5_name(resume_file, "A_vec", arma::hdf5_opts::append));

  previous_x_.save(arma::hdf5_name(resume_file, "previous_x", arma::hdf5_opts::append));
  previous_dx_.save(arma::hdf5_name(resume_file, "previous_dx", arma::hdf5_opts::append));
  previous_A_.save(arma::hdf5_name(resume_file, "previous_A", arma::hdf5_opts::append));

  for (IndexType i = 0U; i < all_interfaces.n_elem; ++i) {
    all_interfaces[i].save(
      arma::hdf5_name(resume_file, "all_interfaces_" + std::to_string(i),
      arma::hdf5_opts::append));
  }

  deltas.save(arma::hdf5_name(resume_file, "deltas", arma::hdf5_opts::append));
  capillaries_with_interfaces.save(arma::hdf5_name(resume_file, "capillaries_with_interfaces",
                                                   arma::hdf5_opts::append));
  plugs.save(arma::hdf5_name(resume_file, "plugs", arma::hdf5_opts::append));
  last_plug_info.save(arma::hdf5_name(resume_file, "last_plug_info", arma::hdf5_opts::append));
  last_plug_inlet_info.save(arma::hdf5_name(resume_file,
                                            "last_plug_inlet_info",
                                            arma::hdf5_opts::append));

  const arma::uvec last_jump_info_size = {last_jump_info.GetNumberOfJumps()};
  last_jump_info_size.save(arma::hdf5_name(resume_file,
                                           "last_jump_info_size",
                                           arma::hdf5_opts::append));

  arma::umat last_jump_info_data(last_jump_info_size[0], 6);

  for (auto it = last_jump_info.cbegin(); it != last_jump_info.cend(); ++it) {
    try {
      auto i = static_cast<unsigned int>(std::distance(last_jump_info.cbegin(), it));
      last_jump_info_data(i, 0) = (*it).GetOrigin();
      last_jump_info_data(i, 1) = static_cast<unsigned int>((*it).GetSide());
      last_jump_info_data(i, 2) = (*it).GetNode();
      last_jump_info_data(i, 3) = (*it).IsCrossJump();
      last_jump_info_data(i, 4) = (*it).IsFirstInCross();

      IndexType destinations_size = (*it).GetNumberOfDestinations();
      last_jump_info_data(i, 5) = destinations_size;

      arma::umat destinations_data(destinations_size, 3);
      for (IndexType d = 0U; d < destinations_size; ++d)  {
        destinations_data(d, 0) =  (*it).GetDestinationCapillary(d);
        destinations_data(d, 1) =  static_cast<unsigned int>((*it).GetDestinationSide(d));
        destinations_data(d, 2) =  static_cast<unsigned int>((*it).GetDestinationType(d));
      }

      destinations_data.save(arma::hdf5_name(resume_file,
                                            "destinations_data_" + std::to_string(i),
                                            arma::hdf5_opts::append));
    } catch(...) {
      LOG(FATAL) << "Error saving Resume file!";
    }
  }

  last_jump_info_data.save(arma::hdf5_name(resume_file,
                                           "last_jump_info_data",
                                           arma::hdf5_opts::append));

  fluid_at_source.save(arma::hdf5_name(resume_file, "fluid_at_source", arma::hdf5_opts::append));
  pressures_.save(arma::hdf5_name(resume_file, "pressures", arma::hdf5_opts::append));
  caps_with_interfaces_not_plugged_.save(arma::hdf5_name(resume_file,
                                                         "caps_with_interfaces_not_plugged",
                                                         arma::hdf5_opts::append));;
  caps_residual_limits_.save(arma::hdf5_name(resume_file,
                                             "caps_residual_limits",
                                             arma::hdf5_opts::append));;
  roots_found_.save(arma::hdf5_name(resume_file, "roots_found", arma::hdf5_opts::append));

  LOG(INFO) << "Resume file saved for time " << current_time_;
}  // end of DynamicCapillaryNetworkAlgorithm::StoreResumeInfo() method



void DynamicCapillaryNetworkAlgorithm::LoadResumeInfo(const N_Vector x,
                                                       const N_Vector dx,
                                                       const SUNMatrix A,
                                                       IndexType *output_time_index,
                                                       double *next_target_time) {
/**
 * The \c LoadResumeInfo method restore the variables state at last time a root was found defined
 * by the at resume file. The following datasets are restored:
 * - \c current_time
 * - \c next_target_time
 * - \c output_time_index
 * - \c resume_next_target_time
 * - \c last_jump_time
 * - \c x
 * - \c dx
 * - \c A
 * - \c previous_x
 * - \c previous_dx
 * - \c previous_A
 * - \c all_interfaces
 * - \c deltas
 * - \c capillaries_with_interfaces
 * - \c plugs
 * - \c last_plug_info
 * - \c last_plug_inlet_info
 * - \c last_jump_info_size
 * - \c last_jump_info_data
 * - \c caps_with_interfaces_not_plugged
 * - \c caps_with_residual_limits
 * - \c roots_found
 */

  // Resume file path
  std::string resume_file = GetFolder() + "/resume.h5";

  if (!std::filesystem::exists(resume_file + ".tgz") && !std::filesystem::exists(resume_file))
    return;

  // Extracting resume file
  if (std::filesystem::exists(resume_file + ".tgz")) {
    std::string extract_command = "tar -xzf " + resume_file + ".tgz -C " + GetFolder() +
                                  " &>/dev/null";
    system(extract_command.c_str());
    std::filesystem::remove(resume_file + ".tgz");

    // Corrupted tgz file
    if (!std::filesystem::exists(resume_file)) return;
  }

  // time variables
  arma::vec time_info;
  arma::uvec output_time_idx;

  arma::Col<IndexType> fluid_at_source;
  arma::field<arma::Col<double>> all_interfaces;
  arma::Col<double> deltas;
  arma::Col<double> pressures;
  arma::Col<IndexType> plugs;
  arma::Col<IndexType> capillaries_with_interfaces;
  arma::Col<double> caps_residual_limits;
  arma::Col<int> roots_found;

  arma::uvec last_jump_info_size;
  arma::umat last_jump_info_data;
  arma::Col<IndexType> last_plug_info;
  arma::Col<IndexType> last_plug_inlet_info;

  arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);
  arma::vec dx_vec(NV_DATA_S(dx), static_cast<arma::uword>(NV_LENGTH_S(dx)), false, true);
  arma::vec A_vec(NV_DATA_S(A), static_cast<arma::uword>(SM_LDATA_D(A)), false, true);

  InterfaceJumpEventInfo last_jump_info;

  const IndexType n_links = context_->number_of_links_;
  all_interfaces.set_size(n_links);

  try {
    time_info.load(arma::hdf5_name(resume_file, "time_info"));
    output_time_idx.load(arma::hdf5_name(resume_file, "output_time_index"));

    last_jump_info_size.load(arma::hdf5_name(resume_file, "last_jump_info_size"));
    last_jump_info_data.load(arma::hdf5_name(resume_file, "last_jump_info_data"));

    x_vec.load(arma::hdf5_name(resume_file, "x_vec"));
    dx_vec.load(arma::hdf5_name(resume_file, "dx_vec"));
    A_vec.load(arma::hdf5_name(resume_file, "A_vec"));

    previous_x_.load(arma::hdf5_name(resume_file, "previous_x"));
    previous_dx_.load(arma::hdf5_name(resume_file, "previous_dx"));
    previous_A_.load(arma::hdf5_name(resume_file, "previous_A"));

    fluid_at_source.load(arma::hdf5_name(resume_file, "fluid_at_source"));
    deltas.load(arma::hdf5_name(resume_file, "deltas"));
    plugs.load(arma::hdf5_name(resume_file, "plugs"));
    pressures.load(arma::hdf5_name(resume_file, "pressures"));
    last_plug_info.load(arma::hdf5_name(resume_file, "last_plug_info"));
    last_plug_inlet_info.load(arma::hdf5_name(resume_file, "last_plug_inlet_info"));
    capillaries_with_interfaces.load(arma::hdf5_name(resume_file, "capillaries_with_interfaces"));
    caps_with_interfaces_not_plugged_.load(arma::hdf5_name(resume_file,
                                                          "caps_with_interfaces_not_plugged"));
    caps_residual_limits_.load(arma::hdf5_name(resume_file, "caps_residual_limits"));
    roots_found_.load(arma::hdf5_name(resume_file, "roots_found"));

    for (IndexType i = 0; i < last_jump_info_size[0]; ++i) {
      IndexType capillary = last_jump_info_data.at(i, 0);
      Side jump_side = static_cast<Side>(last_jump_info_data.at(i, 1));
      IndexType node = last_jump_info_data.at(i, 2);
      bool is_cross = last_jump_info_data.at(i, 3);
      bool first_in_cross = last_jump_info_data.at(i, 4);
      IndexType destinations_size = last_jump_info_data.at(i, 5);

      last_jump_info.AddJumpToJumpVector(capillary, jump_side, node, is_cross, first_in_cross);

      arma::umat destinations_data;
      destinations_data.load(arma::hdf5_name(resume_file,
                                             "destinations_data_" + std::to_string(i)));
      for (IndexType j = 0; j < destinations_size; ++j) {
        IndexType destination_capillary = destinations_data.at(j, 0);
        Side arriving_side = static_cast<Side>(destinations_data.at(j, 1));
        JumpType jump_type = static_cast<JumpType>(destinations_data.at(j, 2));

        last_jump_info.AddDestination(i, destination_capillary, arriving_side, jump_type);
      }
    }
    arma::Col<double> interface;
    for (IndexType i = 0U; i < n_links; ++i) {
      interface.load(arma::hdf5_name(resume_file, "all_interfaces_" + std::to_string(i)));
      all_interfaces[i] = interface;
    }

    current_time_ = time_info(0);
    *next_target_time = time_info(1);
    last_jump_time_ = time_info(2);
  } catch (...) {
    LOG(FATAL) << "Resume file corrupted. Please restart the analysis.";
  }

  *output_time_index = arma::conv_to<IndexType>::from(output_time_idx);

  interface_->SetLastJumpInfo(last_jump_info);
  interface_->SetInterfacePositions(all_interfaces);
  interface_->SetFluidAtSource(fluid_at_source);
  interface_->SetDeltas(deltas);
  interface_->SetPlugs(plugs);
  interface_->SetLastPlugInfo(last_plug_info);
  interface_->SetLastPlugInletInfo(last_plug_inlet_info);
  interface_->SetCapillariesWithInterfaces(capillaries_with_interfaces);
  SetPressures(pressures);
}  // end of DynamicCapillaryNetworkAlgorithm::LoadResumeInfo() method



void DynamicCapillaryNetworkAlgorithm::InitialiseIDASolver(void **mem,
                                                           N_Vector *x,
                                                           N_Vector *dx,
                                                           SUNMatrix *A,
                                                           SUNLinearSolver *LS,
                                                           const double t0) {
/**
 * Initialise the IDA solver using input parameters.
 * Creates the IDA memory block, sets the userdata, identifies the differential
 * and the algebraic equations, allocates the ida vectors, set the residual
 * function, the initial time and the ida vectors.
 *
 * \param[in,out] mem      Standard pointer used by \c IDASolver.
 * \param[in,out] x        SUNDIALS Vector with interface positions.
 * \param[in,out] dx       SUNDIALS Vector with interface movement speed.
 * \param[in] A            SUNDIALS matrix used by the linear solver.
 * \param[in] LS           SUNDIALS linear solver configuration.
 * \param[in] t0           Initial time of the simulation.
 */
  // Casting \c number_of_equations_ to int32_t as SUNDIALS requires this type as inputs
  int32_t n_equations = static_cast<int32_t>(number_of_equations_);

  // Hold IDA return values;
  int retval;

  // Allocating int32_t dynamic memory array as required by SUNDIALS root function
  roots_found_.set_size(0);
  caps_with_interfaces_not_plugged_.set_size(0);
  caps_residual_limits_.set_size(0);

  // Initialise position and velocity vectors
  *x = N_VNew_Serial(n_equations);
  if (check_retval(reinterpret_cast<void *>(x), "N_VNew_Serial", 0))
    LOG(FATAL) << "Error in creation of the `x` vector.";

  *dx = N_VNew_Serial(n_equations);
  if (check_retval(reinterpret_cast<void *>(dx), "N_VNew_Serial", 0))
    LOG(FATAL) << "Error in creation of the `dx` vector.";

  ComputeInitialValues(*x, *dx);

  *A = SUNDenseMatrix(n_equations, n_equations);
  if (check_retval(reinterpret_cast<void *>(A), "SUNDenseMatrix", 0))
    LOG(FATAL) << "Error in creation of the `A` Sundials Dense Matrix.";

  /* Create dense SUNLinearSolver object */
  *LS = SUNLinSol_Dense(*x, *A);
  if (check_retval(reinterpret_cast<void *>(*LS), "SUNDenseLinearSolver", 0))
    LOG(FATAL) << "Error in creation of the `LS` Sundials Dense Linear Solver.";

  double *atval;
  current_time_ = t0;
  last_jump_time_ = t0;

  N_Vector avtol = nullptr;

  // Set link_tol and node_tol to the SUNDIALS equations
  avtol = N_VNew_Serial(n_equations);
  if (check_retval(reinterpret_cast<void *>(avtol), "N_VNew_Serial", 0))
    LOG(FATAL) << "Error in creation of the `avtol` vector of error tolerances.";
  atval = N_VGetArrayPointer(avtol);
  arma::vec eq_vec(atval, static_cast<arma::uword>(NV_LENGTH_S(avtol)), false, true);
  eq_vec.head(context_->number_of_links_).fill(absolute_link_tolerance_);
  eq_vec.tail(context_->interior_nodes_.n_elem).fill(absolute_node_tolerance_);

  // Set id and userdata values
  N_Vector id = N_VNew_Serial(n_equations);
  SetIDValues(id);

  // Create IDA memory block
  *mem = IDACreate();

  // Set userdata(container to forward data through the solver)
  IDASetUserData(*mem, this);

  // identifies the differential and the algebraic equations
  IDASetId(*mem, id);

  // Set the residual function(Residual)
  IDAInit(*mem, Residual, t0, *x, *dx);

  // Set error tolerances
  retval = IDASVtolerances(*mem, relative_tolerance_, avtol);
  if (check_retval(&retval, "IDASVtolerances", 1)) exit(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(*mem, *LS, *A);
  if (check_retval(&retval, "IDADlsSetLinearSolver", 1))
    LOG(FATAL) << "Error in attach the `A` matrix and `LS` Sundials Dense Linear Solver.";
}  // end of DynamicCapillaryNetworkAlgorithm::InitialiseIDASolver() method



void DynamicCapillaryNetworkAlgorithm::ComputeInitialValues(const N_Vector x, const N_Vector dx) {
/**
 * This method computes initial \f$ \dot{x_i} \f$ values based on the initial state of the system.
 * Initial values are computed as described in ComputeDifferentialEquation().
 *
 * \param[in] x       SUNDIALS Vector with effective interface positions.
 * \param[in] dx      SUNDIALS Vector with effective interface speeds.
 */

  arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);
  arma::vec dx_vec(NV_DATA_S(dx), static_cast<arma::uword>(NV_LENGTH_S(dx)), false, true);

  // Fill first elements related to interface positions
  for (IndexType i = 0U; i < context_->number_of_links_; ++i) {
    if (interface_->IsNumberOfInterfacesAtCapillaryOdd(i)) {
      x_vec(i) = interface_->CalculateEffectiveInterfacePositionAtCapillary(i);
    } else {
      x_vec(i) = interface_->GetInterfacePositionAtCapillaryNearSide(i, Side::source);
    }
  }

  // Fill last elements related to interior pressures
  x_vec.tail(context_->interior_nodes_.n_elem).zeros();

  // Populate pressures_ with the values from x
  SetPressures(context_->interior_nodes_, x_vec.tail(context_->interior_nodes_.n_elem));

  // Fill last elements related to interior pressures
  x_vec.tail(context_->interior_nodes_.n_elem).zeros();

  // For each differential equation, use the equation according to the number of interfaces
  ComputeDifferentialEquationForAllCapillaries(dx_vec, x_vec);

  // For each algebraic equation, impose mass conservation
  dx_vec.tail(context_->interior_nodes_.n_elem).zeros();
}  // end of DynamicCapillaryNetworkAlgorithm::ComputeInitialValues() method



void DynamicCapillaryNetworkAlgorithm::ComputeDifferentialEquationForAllCapillaries(
  arma::Col<double> &differential,
  const arma::Col<double> &x_array) {
/**
 * This method computes the effective interface speed \f$ \dot{x} \f$ for each capillary in the
 * system.
 *
 * This method call \c ComputeDifferentialEquation() for each capillary.
 *
 * \param[out] differential   Effective interface speed of all capillaries.
 * \param[in]  x_array        Double array with effective interface positions. When the number of
 *                            interfaces is odd, x_array also represents the saturation in the
 *                            capillary. When the number of interfaces is even, the saturation is
 *                            constant and given by CalculateEffectiveInterfacePositionAtCapillary
 *                            and \c x_array vary by the rate of the velocity of the fluid.
 */

  #pragma omp parallel for
  for (IndexType i = 0U; i < context_->number_of_links_; ++i) {
    differential(i) = ComputeDifferentialEquation(i, x_array);
  }
}  // end of DynamicCapillaryNetworkAlgorithm::ComputeDifferentialEquationForAllCapillaries() method



void DynamicCapillaryNetworkAlgorithm::ComputeDifferentialEquationForSelectedCapillaries(
  arma::Col<double> &differential,
  const arma::Col<double> &x_array,
  const arma::Col<IndexType> &capillaries) {
/**
 * This method computes the effective interface speed \f$ \dot{x} \f$ for each selected capillary.
 *
 * This method call \c ComputeDifferentialEquation() for each selected capillary.
 *
 * \param[out] differential   Effective interface speed of selected capillaries.
 * \param[in]  x_array        Double array with effective interface positions. When the number of
 *                            interfaces is odd, x_array also represents the saturation in the
 *                            capillary. When the number of interfaces is even, the saturation is
 *                            constant and given by CalculateEffectiveInterfacePositionAtCapillary
 *                            and \c x_array vary by the rate of the velocity of the fluid.
 * \param[out] capillaries    Selected capillaries used to compute the effective interface speed.
 */
  #pragma omp parallel for
  for (IndexType i = 0U; i < capillaries.n_elem; ++i) {
    IndexType capillary = capillaries(i);
    differential(capillary) = ComputeDifferentialEquation(capillary, x_array);
  }
}  // end of DynamicCapillaryNetworkAlgorithm::ComputeDifferentialEquationForSelectedCapillaries()



double DynamicCapillaryNetworkAlgorithm::ComputeDifferentialEquation(
  const IndexType link_index,
  const arma::Col<double> &x_array) {
/**
 * This method computes the effective interface speed \f$ \dot{x} \f$ inside a given capillary.
 *
 * If the capillary is plugged
 *
 * \f[
 *    \dot{x} = 0
 * \f]
 *
 * Else,
 *
 * &nbsp;&nbsp; if the capillary has no interfaces
 *
 *   \f[
 *      \dot{x} = \frac{D^2 \Delta P}{32 L^2 \eta_s}
 *   \f]
 *
 * &nbsp;&nbsp; if the capillary has an odd number of interfaces
 *
 *   \f[
 *      \dot{x} = \frac{D [4 \sigma \cos \theta + D \Delta P]}{4 L [C_1 C_2 D
 *                + 8 L (\eta_t + (\eta_s - \eta_t) x)]}
 *   \f]
 *
 * &nbsp;&nbsp; if the capillary has an even number of interfaces
 *
 *   \f[
 *      \dot{x} = \frac{D^2 \Delta P}{32 L^2 (\eta_t + (\eta_s - \eta_t) x)}
 *   \f]
 *
 * Where:
 *  - \f$ x \f$ is the effective interface position at capillary \c link_index.
 *  - \f$ \dot{x} \f$ is the effective interface speed at capillary \c link_index.
 *  - \f$ D \f$ is the diameter of capillary \c link_index.
 *  - \f$ \sigma \f$ is the interfacial tension between the two immiscible fluids.
 *  - \f$ \theta \f$ is the contact angle at capillary \c link_index.
 *  - \f$ \Delta P \f$ is the pressure difference between the nodes in both ends of capillary \c link_index.
 *  - \f$ L \f$ is the length of capillary \c link_index.
 *  - \f$ C_1 \f$ and \f$ C_2 \f$ are fitting parameters of the MK model considering the fluids and the rock.
 *  - \f$ \eta_s \f$ and \f$ \eta_t \f$ are the fluid viscosity in the source and target nodes of the capillary, respectively.
 *
 * \param[out] link_index     Effective interface speed of capillary.
 * \param[in]  x_array        Double array with effective interface positions. When the number of
 *                            interfaces is odd, x_array also represents the saturation in the
 *                            capillary. When the number of interfaces is even, the saturation is
 *                            constant and given by CalculateEffectiveInterfacePositionAtCapillary
 *                            and \c x_array vary by the rate of the velocity of the fluid.
 * \retval     differential   Effective interface speed of capillary.
 */
  LinearMolecularKineticsPhysics *lmkp = static_cast<LinearMolecularKineticsPhysics *>(
    model_->GetPhysics().get());
  const double THETA = lmkp->GetWettabilityContactAngleInRadians();
  const double SIGMA = lmkp->GetFluidInterfaceInterfacialTension();
  const double MK = lmkp->GetWettabilityLinearMK();

  const double Di = 2.0 * std::sqrt(context_->link_squared_radius_(link_index));

  const double source_node_pressure = pressures_(context_->linked_nodes_(0, 2 * link_index));
  const double target_node_pressure = pressures_(context_->linked_nodes_(0, 2 * link_index + 1));

  const double DeltaP = source_node_pressure - target_node_pressure;

  const double avg_pressure = (source_node_pressure + target_node_pressure) / 2.0;

  const double fluid_viscosities[] = { lmkp->GetResidentFluidViscosity(avg_pressure),
                                       lmkp->GetInjectedFluidViscosity(avg_pressure) };

  const double Li = context_->link_length_(link_index);
  const IndexType fluid_at_source = interface_->GetFluidAtCapillaryAtSide(link_index, Side::source);
  const double Eta_s = fluid_viscosities[fluid_at_source];
  const double Eta_t = fluid_viscosities[!fluid_at_source];
  const double cosTHETA = fluid_at_source == 1 ? std::cos(THETA) : -std::cos(THETA);
  double Eta_eff;  // effective viscosity

  double differential;

  // TODO(rneumann): Explore branch optimisation by putting the most likely condition at the top
  // Having 0 interfaces is more probable than having non-zero interfaces
  if (interface_->IsCapillaryPlugged(link_index)) {
    differential = 0.0;
  } else {
    if (interface_->IsNumberOfInterfacesAtCapillaryOdd(link_index)) {
      Eta_eff = (Eta_t + (Eta_s * x_array(link_index))) - Eta_t * x_array(link_index);
      differential = (Di * (4.0 * SIGMA * cosTHETA + Di * DeltaP))
                   / (4.0 * Li * (MK * Di + 8.0 * Li * Eta_eff));
    } else {
      if (interface_->GetNumberOfInterfacesAtCapillary(link_index) == 0U) {
        Eta_eff = Eta_s;
      } else {
        const double eff_int_pos =
          interface_->CalculateEffectiveInterfacePositionAtCapillary(link_index);
        Eta_eff = (Eta_t + (Eta_s * eff_int_pos)) - (Eta_t * eff_int_pos);
      }
      differential = (Di * Di * DeltaP) / (32.0 * Li * Li * Eta_eff);
    }
  }

  // Log verbose information
  VLOG(1) << link_index << ":Di[" << Di << "]:DeltaP[" << DeltaP << "]:Li[" << Li
      << "]:fluid_at_source[" << fluid_at_source
      << "]:Eta_s[" << Eta_s << "]:Eta_t[" << Eta_t << "]:Eta_eff[" << Eta_eff
      << "]:cosTHETA[" << cosTHETA << "]:differential[" << differential
      << "]:x_array[" << x_array(link_index)
      << "]:IsCapillaryPlugged[" << interface_->IsCapillaryPlugged(link_index)
      << "]:NumberOfInterfaces[" << interface_->GetNumberOfInterfacesAtCapillary(link_index)
      << "]";
  return differential;
}  // end of DynamicCapillaryNetworkAlgorithm::ComputeDifferentialEquation() method



void DynamicCapillaryNetworkAlgorithm::SetIDValues(const N_Vector id) {
/**
 * Initialises IDA array that identifies the differential and the algebraic
 * equations. The differential equations are associated to the capillary and
 * are numbered first, following the indexing of the \c linked_nodes_.
 * The \c id array receives 1 for differential equations and 0 for algebraic equations
 *
 * \param[in] id       SUNDIALS Vector to receive the values.
 */
  arma::vec id_vec(NV_DATA_S(id), static_cast<arma::uword>(NV_LENGTH_S(id)), false, true);

  // Set first differential equations to 1 and last algebraic equations to 0
  id_vec.head(context_->number_of_links_).ones();
  id_vec.tail(context_->interior_nodes_.n_elem).zeros();
}  // end of DynamicCapillaryNetworkAlgorithm::SetIDValues() method



void DynamicCapillaryNetworkAlgorithm::RemoveBubblesFromOrigin(
  arma::Col<IndexType> &caps_with_roots,
  arma::Col<IndexType> &caps_roots_directions,
  arma::Col<IndexType> &changed_caps) {
/**
 * A "bubble" is a pair interfaces, very close to each other, that jump together. This method removes those
 * bubbles. Because our simulator deals with one interface per jump per capillary, this method
 * applies the required elimination of bubbles to return to the situation handled by our simulator
 *
 * There are two possibilities when removing bubbles. One is that it has no more interfaces to be
 * processed, so this method removes those capillaries from the list of capillaries with roots. The
 * other is that there is a single interface left in the capillary, and nothing more needs to be done.
 *
 * \param[in, out] caps_with_roots        Capillaries that were caught at the jump event.
 * \param[in, out] caps_roots_directions  Inform the direction of the jump event. 0 for source,
 *                                        1 for target.
 * \param[out]     changed_caps           Capillaries that were removed from \c caps_with_roots.
 */
  arma::Col<IndexType> caps_to_remove;

  // Certify that the capillaries that jumped have only one interface beyond their edges
  for (IndexType i = 0U; i < caps_with_roots.n_elem; ++i) {
    IndexType cap = caps_with_roots(i);
    Side side = (caps_roots_directions(i) == 0U ? Side::source : Side::target);
    arma::vec interface_positions_at_cap = interface_->GetInterfacePositionsAtCapillary(cap);

    // Calculate number of interfaces outside the [0,1] interval
    arma::uword number_of_interfaces_beyond_edge =
      (side == Side::source) ? arma::find(interface_positions_at_cap <= 0.0).eval().n_elem
                             : arma::find(interface_positions_at_cap >= 1.0).eval().n_elem;

    // If more than one interface crossed the edge at the same time, remove that bubble
    if (number_of_interfaces_beyond_edge > 1U) {
      LOG(INFO) << "Removing bubbles from " << side;
      LOG(INFO) << "number_of_interfaces_beyond_edge[" << number_of_interfaces_beyond_edge << "]";

      // Remove all interfaces that make up the bubbles. Leave 1 if number of interfaces is odd.
      arma::uword is_number_odd = number_of_interfaces_beyond_edge % 2U;
      for (IndexType j = 0U; j < number_of_interfaces_beyond_edge - is_number_odd; ++j) {
        interface_->RemoveInterface(cap, side);
      }

      // Also remove capillary from processing if number of interfaces is even
      if (!is_number_odd) {  // even number of interfaces beyond edge
        LOG(INFO) << "Removing capillary[" << caps_with_roots(i) << "] from processing";
        caps_to_remove.insert_rows(caps_to_remove.n_elem, arma::Col<IndexType>(
          { caps_with_roots(i) }));
      }
    }
  }

  // caps to remove
  for (const auto &capillary : caps_to_remove) {
    LOG(INFO) << "capillary: " << capillary << std::flush;
    LOG(INFO) << "caps_with_roots: " << caps_with_roots.t() << std::flush;
    LOG(INFO) << "caps_roots_directions: " << caps_roots_directions.t() << std::flush;
    IndexType index = arma::as_scalar(arma::find(caps_with_roots == capillary));
    LOG(INFO) << "index: " << index << std::flush;
    caps_with_roots.shed_row(index);
    caps_roots_directions.shed_row(index);
    changed_caps.insert_rows(changed_caps.n_elem, arma::Col<IndexType>({ capillary }));
  }
}  // end DynamicCapillaryNetworkAlgorithm::RemoveBubblesFromOrigin() method



void DynamicCapillaryNetworkAlgorithm::PlugCapillariesWhoseInterfacesJumpedThroughInletNodes(
  arma::Col<IndexType> &caps_with_roots,
  arma::Col<IndexType> &caps_roots_directions,
  arma::Col<IndexType> &changed_caps) {
/**
 * An inlet-node must inject fluid other than the resident fluid. When an interface jumps over an
 * input node, that node can inject the resident fluid. This method plugs capillaries whose
 * interfaces jump over inlet-nodes, preventing these nodes from being transformed into a resident
 * fluid inlet-node.
 *
 * \param[in, out] caps_with_roots        Capillaries that were caught at the jump event.
 * \param[in, out] caps_roots_directions  Inform the direction of the jump event. 0 for source,
 *                                        1 for target.
 * \param[out]     changed_caps           Capillaries that were removed from \c caps_with_roots.
 */
  // Remove capillaries plugged because of jumps through inlet nodes
  arma::Col<IndexType> caps_to_remove = interface_->PlugCapillariesJumpAtInletNode(
    caps_with_roots, caps_roots_directions);
  for (const auto &capillary : caps_to_remove) {
    LOG(INFO) << "capillary: " << capillary;
    LOG(INFO) << "caps_with_roots: " << caps_with_roots.t();
    LOG(INFO) << "caps_roots_directions: " << caps_roots_directions.t();
    IndexType index = arma::as_scalar(arma::find(caps_with_roots == capillary));
    LOG(INFO) << "index: " << index;
    caps_with_roots.shed_row(index);
    caps_roots_directions.shed_row(index);
    changed_caps.insert_rows(changed_caps.n_elem, arma::Col<IndexType>({ capillary }));
  }
}  // end of PlugCapillariesWhoseInterfacesJumpedThroughInletNodes() method



void DynamicCapillaryNetworkAlgorithm::InitialiseIDAState(void *mem,
                                                          const N_Vector x,
                                                          const N_Vector dx,
                                                          const SUNMatrix A,
                                                          const SUNLinearSolver LS,
                                                          const double t0,
                                                          const double next_target_time,
                                                          const IndexType output_time_index) {
/**
 * Sets the values of the current state to start the solver.
 * Updates the interface information, plugs and unplugs capillaries, updates
 * fluid at source of capillaries, computes new initial values, reinitialises
 * the solver for current state and corrects the initial values.
 *
 * \param[in] mem                       Standard pointer used by \c IDASolver.
 * \param[in] x                         SUNDIALS Vector with interface positions.
 * \param[in] dx                        SUNDIALS Vector with interface movement speed.
 * \param[in] A                         SUNDIALS matrix used by the linear solver.
 * \param[in] LS                        SUNDIALS linear solver configuration.
 * \param[in] t0                        Initial time of the simulation.
 * \param[in] next_target_time          Actual target time of the simulation.
 * \param[in] output_time_index         Output time index of the simulation.
 */
  // Updates the interface information based on movement of effective interface
  arma::vec x_vec(NV_DATA_S(x), static_cast<arma::uword>(NV_LENGTH_S(x)), false, true);
  arma::vec dx_vec(NV_DATA_S(dx), static_cast<arma::uword>(NV_LENGTH_S(dx)), false, true);
  arma::vec A_vec(NV_DATA_S(A), static_cast<arma::uword>(SM_LDATA_D(A)), false, true);
  if (current_time_ != t0) {
    IDAGetRootInfo(mem, roots_found_.memptr());
    const arma::uvec roots_idx(arma::find(roots_found_));
    arma::Col<IndexType> caps_with_roots;
    arma::Col<IndexType> caps_roots_directions;
    arma::Col<IndexType> changed_caps;

     // Selecting capillaries whose interfaces have jumped
    caps_with_roots.set_size(0);
    caps_roots_directions.set_size(0);
    for (const auto &id : roots_idx) {
      caps_with_roots.insert_rows(
        caps_with_roots.n_elem,
        arma::Col<IndexType>({ caps_with_interfaces_not_plugged_(id / 2U) }));

      // Insert 1 in case of interface jumping at target, 0 in case of interface jumping at source
      caps_roots_directions.insert_rows(caps_roots_directions.n_elem,
                                        arma::Col<IndexType>({ id % 2U }));
    }

    RemoveBubblesFromOrigin(caps_with_roots, caps_roots_directions, changed_caps);

    PlugCapillariesWhoseInterfacesJumpedThroughInletNodes(caps_with_roots,
                                                          caps_roots_directions,
                                                          changed_caps);

    // Check if plugging occurred in capillaries in which roots were found
    bool plugged = false;
    double plug_time_tol = 1e-6;
    if (std::abs(current_time_ - last_jump_time_) < plug_time_tol) {
      plugged = interface_->CheckPlugging(caps_with_roots, caps_roots_directions);
    }

    // Log relevant information
    std::stringstream sstream;
    sstream << "Roots found : " << caps_with_roots.t();
    LOG(INFO) << sstream.str();
    LOG(INFO) << "current_time_[" << current_time_ << "]:last_jump_time_[" << last_jump_time_
              << "]:plug_time_tol[" << plug_time_tol << "].";

    // Act upon plugging event
    if (plugged) {
      // Undo latest changes applied to SUNDIALS state variables
      RevertState(x_vec, dx_vec, A_vec);
      interface_->RevertUpdateInterfaceInformation();
      interface_->RevertJump(changed_caps);
      interface_->RevertPluggedCapillariesJumpAtInletNode(changed_caps);
      // Reverts Sundial variables for before jump
      UpdateSundialsVariables(x_vec, dx_vec, changed_caps);
      current_time_ = last_jump_time_;
    } else {
      // Apply changes to SUNDIALS state variables
      interface_->Jump(caps_with_roots, caps_roots_directions, changed_caps);
      // Updates Sundial variables with the new changes
      UpdateSundialsVariables(x_vec, dx_vec, changed_caps);
      StoreState(x_vec, dx_vec, A_vec);
      last_jump_time_ = current_time_;
    }

    number_of_roots_++;
    if (number_of_roots_ % resume_save_frequency_ == 0) {
      StoreResumeInfo(x, dx, A, output_time_index, next_target_time);
    }
  }


  // Reinitialise SUNDIALS IDA solver with the new state
  IDAReInit(mem, current_time_, x, dx);

  // Rebuild list of capillaries with interface that are not plugged
  const arma::Col<IndexType> capillaries = interface_->GetIndexOfCapillariesWithInterfaces();
  caps_with_interfaces_not_plugged_.set_size(0);
  caps_residual_limits_.set_size(0);
  for (const auto &capillary : capillaries) {
    if (!interface_->IsCapillaryPlugged(capillary)) {
      caps_with_interfaces_not_plugged_.insert_rows(
        caps_with_interfaces_not_plugged_.n_elem, arma::Col<IndexType>({ capillary }));

      // Get the effective, minimum and maximum interface values
      const double eff_x = interface_->CalculateEffectiveInterfacePositionAtCapillary(capillary);
      const double min_x = interface_->GetInterfacePositionAtCapillaryNearSide(capillary,
                                                                               Side::source);
      const double max_x = interface_->GetInterfacePositionAtCapillaryNearSide(capillary,
                                                                               Side::target);

      double lower_limit, upper_limit;

      if (interface_->IsNumberOfInterfacesAtCapillaryOdd(capillary)) {
        lower_limit = x_vec(capillary) - min_x;
        upper_limit = x_vec(capillary) + (1.0 - max_x);
      } else {
        lower_limit = 0.0;
        upper_limit = x_vec(capillary) + (1.0 - max_x);
      }

      caps_residual_limits_.insert_rows(
        caps_residual_limits_.n_elem, arma::Col<double>({ lower_limit }));

      caps_residual_limits_.insert_rows(
        caps_residual_limits_.n_elem, arma::Col<double>({ upper_limit }));

      if (x_vec(capillary) > 1.0 ||
          x_vec(capillary) < 0.0 ||
          min_x > 1.0 ||
          min_x < 0.0 ||
          max_x > 1.0 ||
          max_x < 0.0 ||
          (interface_->IsNumberOfInterfacesAtCapillaryOdd(capillary) &&
            eff_x != x_vec(capillary)) ||
          (!interface_->IsNumberOfInterfacesAtCapillaryOdd(capillary) &&
            x_vec(capillary) != min_x)) {
        LOG(INFO) << "Capillary[" << capillary << "]:eff_x[" << eff_x << "]:min_x[" << min_x
          << "]:max_x[" << max_x << "]:x_vec[" << x_vec[capillary] << "]:dx_vec["
          << dx_vec(capillary)
          << "]:NumberOfInterfaces[" << interface_->GetNumberOfInterfacesAtCapillary(capillary)
          << "]:lower_limit[" << lower_limit << "]:upper_limit[" << upper_limit
          << "]";
      }
    }
  }
  roots_found_.set_size(caps_residual_limits_.n_elem);
  IDARootInit(mem, roots_found_.n_elem, DetectJumps);

  IDASetLinearSolver(mem, LS, A);
  IDACalcIC(mem, IDA_YA_YDP_INIT, next_target_time);
  IDAGetConsistentIC(mem, x, dx);

  for (IndexType i = 0U; i < caps_with_interfaces_not_plugged_.n_elem; ++i) {
    if (x_vec(caps_with_interfaces_not_plugged_(i)) <= caps_residual_limits_(2U * i)) {
      LOG(INFO) << "Capillary: " << caps_with_interfaces_not_plugged_(i);
      LOG(INFO) << "x_vec(caps_with_interfaces_not_plugged_(i)) pass through lower limit";
      LOG(INFO) << "x_vec(..): " << x_vec(caps_with_interfaces_not_plugged_(i));
      LOG(INFO) << "caps_residual_limits_(2U * i): " << caps_residual_limits_(2U * i);
    }
    if (x_vec(caps_with_interfaces_not_plugged_(i)) >= caps_residual_limits_(2U * i + 1U)) {
      LOG(INFO) << "Capillary: " << caps_with_interfaces_not_plugged_(i);
      LOG(INFO) << "x_vec(caps_with_interfaces_not_plugged_(i)) pass through upper limit";
      LOG(INFO) << "x_vec(..): " << x_vec(caps_with_interfaces_not_plugged_(i));
      LOG(INFO) << "caps_residual_limits_(2U * i + 1U): " << caps_residual_limits_(2U * i + 1U);
    }
  }
}  // end of DynamicCapillaryNetworkAlgorithm::InitialiseIDAState() method



void DynamicCapillaryNetworkAlgorithm::UpdateSundialsVariables(
  arma::Col<double> &x,
  arma::Col<double> &dx,
  const arma::Col<IndexType> &capillaries) {
/**
 * Update the interface positions vector \c x with the appropriate values.
 *
 * \param[in] x             Double vector with interface positions and interior pressures
 * \param[in] dx            Double vector with interface velocity
 * \param[in] capillaries   Capillaries that have some modification
 */
  for (auto capillary : capillaries) {
    LOG(INFO) << "capillary: " << capillary << std::flush;
    LOG(INFO) << "all_interfaces: "
      << interface_->GetInterfacePositionsAtCapillary(capillary).t() << std::flush;
    if (interface_->IsNumberOfInterfacesAtCapillaryOdd(capillary)) {
      x(capillary) = interface_->CalculateEffectiveInterfacePositionAtCapillary(capillary);
      LOG(INFO) << "x(" << capillary << "): " << x(capillary) << std::flush;
    } else {
      x(capillary) = interface_->GetInterfacePositionAtCapillaryNearSide(capillary, Side::source);
      LOG(INFO) << "x(" << capillary << "): " << x(capillary) << std::flush;
    }
  }
  ComputeDifferentialEquationForSelectedCapillaries(dx, x, capillaries);
}  // end of DynamicCapillaryNetworkAlgorithm::UpdateSundialsVariables() method



void DynamicCapillaryNetworkAlgorithm::RevertState(arma::Col<double> &x,
                                                   arma::Col<double> &dx,
                                                   arma::Col<double> &A) {
/**
 * Revert the pressure of interior nodes in \c pressures_ to a previous state.
 *
 * \param[in] x   Vector with interface positions and interior pressures.
 * \param[in] dx  Vector with interface positions ratio and interior pressures ratio.
 * \param[in] A   Matrix with coefficients of the linear equations.
 */
  x = previous_x_;
  dx = previous_dx_;
  A = previous_A_;
}  // end of DynamicCapillaryNetworkAlgorithm::RevertPressures() method



void DynamicCapillaryNetworkAlgorithm::RevertPressures(arma::Col<double> &x) {
/**
 * Revert the pressure of interior nodes in \c pressures_ to a previous state.
 *
 * \param[in] x   Double vector with interface positions and interior pressures
 */
  x.tail(context_->interior_nodes_.n_elem) = previous_interior_pressures_;
}  // end of DynamicCapillaryNetworkAlgorithm::RevertPressures() method


void DynamicCapillaryNetworkAlgorithm::StoreState(const arma::Col<double> &x,
                                                  const arma::Col<double> &dx,
                                                  const arma::Col<double> &A) {
/**
 * Stores the pressure in the interiors with current values to be used in next state change
 * if a capillary becomes plugged
 *
 * \param[in] x   Double vector with interface positions and interior pressures
 * \param[in] dx  Vector with interface positions ratio and interior pressures ratio.
 * \param[in] A   Matrix with coefficients of the linear equations.
 */
  previous_x_ = x;
  previous_dx_ = dx;
  previous_A_ = A;
}  // end of DynamicCapillaryNetworkAlgorithm::StoreInteriorPressures() method



void DynamicCapillaryNetworkAlgorithm::StoreInteriorPressures(const arma::Col<double> &x) {
/**
 * Stores the pressure in the interiors with current values to be used in next state change
 * if a capillary becomes plugged
 *
 * \param[in] x   Double vector with interface positions and interior pressures
 */
  previous_interior_pressures_ = x.tail(context_->interior_nodes_.n_elem);
}  // end of DynamicCapillaryNetworkAlgorithm::StoreInteriorPressures() method



bool DynamicCapillaryNetworkAlgorithm::StopCriteriaMet(const double tFinal) {
/**
 * Checks if any stop criteria is met, so the solver is finished.
 *
 * \param[in] tFinal          Final simulation time expected by the user.
 * \retval    criterion_met   Boolean value representing whether any stop criterion was met.
 */
  return current_time_ >= tFinal;
}  // end of DynamicCapillaryNetworkAlgorithm::StopCriteriaMet() method



bool DynamicCapillaryNetworkAlgorithm::CheckIDAState(const int state) {
/**
 * Reads the \c IDASolve output and determine if there was a state change and stores it.
 *
 * \param[in] state         DAE state information.
 * \retval    state_changed Boolean value representing whether a state change occurred.
 */
  if (state < 0) {  // fail
    // something went wrong
    std::fprintf(stderr, "Error in IDASOLVE.\n");
    std::exit(-1);
  }

  bool state_changed = false;
  switch (state) {
    case IDA_SUCCESS:     state_changed = false;  break;
    case IDA_ROOT_RETURN: state_changed = true;   break;
  }

  return state_changed;
}  // end of DynamicCapillaryNetworkAlgorithm::CheckIDAState() method



int DynamicCapillaryNetworkAlgorithm::check_retval(void *returnvalue,
                                                   const char *funcname,
                                                   int opt) {
  /**
   * Check function return value.
   *
   *   opt == 0 means SUNDIALS function allocates memory so check if
   *            returned NULL pointer
   *   opt == 1 means SUNDIALS function returns an integer value so check if
   *            retval < 0
   *   opt == 2 means function allocates memory so check if returned
   *            NULL pointer
   *
   * \param[in] returnvalue return value.
   * \param[in] funcname    function name.
   * \param[in] opt         0 for SUNDIALS memory allocation.
   *                        1 for the guarantee that return value is positive
   *                        2 for function memory allocation.
   * \retval    retval      0 if no ERROR, 1 o/w.
   */
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    LOG(INFO) << "\nSUNDIALS_ERROR: " << funcname << "() failed - returned NULL pointer\n\n";
    return(1);
  } else if (opt == 1) {
    /* Check if retval < 0 */
    retval = reinterpret_cast<int *>(returnvalue);
    if (*retval < 0) {
      LOG(INFO) << "\nSUNDIALS_ERROR: " << funcname << " failed with retval = " << *retval
          << "\n\n";
      return(1);
    }
  } else if (opt == 2 && returnvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    LOG(INFO) << "\nMEMORY_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return(1);
  }
  return(0);
}

}  // namespace simulator
