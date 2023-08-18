/**
 * \file src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.h
 * \brief Contains the \c DynamicCapillaryNetwork class.
 *
 * \authors Giulia Duncan Coutinho \<coutingi@br.ibm.com\>
 * \authors Italo Cristiano Nievinski Lima \<italon@br.ibm.com\>
 * \authors Rodrigo Neumann Barros Ferreira \<rneumann@br.ibm.com\>
 * \authors Rodrigo Alves Prado da Silva \<rodrigo.alves@ibm.com\>
 * \copyright © IBM Corp.
 * \date 2017
 *
 * This header file contains the \c DynamicCapillaryNetworkAlgorithm class that derives from \c IAlgorithm.
 */


#ifndef SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_DYNAMIC_CAPILLARY_NETWORK_H_
#define SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_DYNAMIC_CAPILLARY_NETWORK_H_

#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <armadillo>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include "src/flow_simulator/algorithms/i_algorithm.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/models/dynamic_model_base.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/capillary/capillary_interface.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network_context.h"
#include "src/flow_simulator/algorithms/dynamic_capillary_network/experiments/dynamic_capillary_experiment_base.h"
#include "src/exec_manager/simulation_config.h"

namespace simulator {

/**
 * \class DynamicCapillaryNetworkAlgorithm dynamic_capillary_network.h "src/flow_simulator/algorithms/dynamic_capillary_network/dynamic_capillary_network.h"
 * \brief Simulator algorithm for Dynamic Capillary Network using centreline geometries
 *
 * The Dynamic Capillary Network Algorithm performs fluid flow simulations of two fluids as a
 * function of time using a set of connected centrelines as representation of the porous rock
 * geometry.
 */

class DynamicCapillaryNetworkAlgorithm : public IAlgorithm {
 public:
  /// Build and configure all resources needed for the algorithm.
  std::pair<int, std::string> Initialise(SimulationConfig &simulation_cfg);

  /// Builds the required geometric representation depending on the simulation.
  void BuildGeometricalRepresentation(void);

  /// Builds the physical equations of the Capillary Network Algorithm.
  void BuildPhysicalEquations(void);

  /// Solves the physical equation of the Capillary Network Algorithm.
  void SolvePhysicalEquations(void);

  /// Calculate derived quantities.
  void CalculateDerivedQuantities(void);

  /// Saves the resulting quantities to disk.
  void SaveResultsToDisk(void);

  /// Close all resources and log summary informations.
  void Finalize(void);

  /// Setter for \c pressures_(index)
  void SetPressure(const IndexType index, const double value) { pressures_(index) = value; }

  /// Setter for \c pressures_
  void SetPressures(const arma::Col<double> &pressures) { pressures_ = pressures; }

  /// Setter for \c pressures_(indexes, updated pressures)
  void SetPressures(const arma::Col<arma::uword> &indexes,
                    const arma::Col<double> &pressures) {
    pressures_(indexes) = pressures;
  }

  /// Getter for \c pressures_.
  const arma::Col<double> &GetPressures(void) const { return pressures_; }

  /// Getter for \c results_folder_.
  const std::string &GetFolder(void) const { return model_->GetGeometry()->folder_; }

  /// Getter for \c interface_.
  std::shared_ptr<CapillaryInterface> GetInterface(void) const { return interface_; }

  /// Getter for \c number of equations.
  IndexType GetNumberOfEquations(void) const { return number_of_equations_; }

  /// Getter for \c context.
  std::shared_ptr<DynamicCapillaryNetworkContext> GetContext(void) const { return context_; }

  /// Function called by IDASolve that find function roots that indicate capillary jumps.
  static int DetectJumps(const double t,
                         const N_Vector x,
                         const N_Vector dx,
                         double *gout,
                         void *user_data);

  /// Function called by IDASolve that compute the DAE residual.
  static int Residual(const double t,
                      const N_Vector x,
                      const N_Vector dx,
                      N_Vector resval,
                      void *user_data);

  /// Function called to check IDA return values.
  static int check_retval(void *returnvalue, const char *funcname, int opt);

  /**
   * \brief Stores the name of the folder where the snapshots will be saved
   *
   * The \c snapshots_folder stores the name of the folder where the snapshots will be saved during
   * the simulation.
   */
  const std::string snapshots_folder_ = "/snapshots";

 private:
  /// Calculates the number of equations in the DAE and saves it to \c number_of_equations_
  void SetNumberOfEquations(void) { number_of_equations_ = context_->number_of_links_
                                                         + context_->interior_nodes_.n_elem;
  }

  /// Sets up the save folder for the storage of results of the dynamic simulation
  void CreateSnapshotFolder(void);

  /// Store the partial results.
  void StoreStepResults(const N_Vector x,
                        const N_Vector dx,
                        IndexType output_time_index);

  /// Store partial results for resuming.
  void StoreResumeInfo(const N_Vector x,
                        const N_Vector dx,
                        const SUNMatrix A,
                        IndexType output_time_index,
                        double next_target_time);

  /// Loads last last saved state to resume the analysis
  void LoadResumeInfo(const N_Vector x,
                       const N_Vector dx,
                       const SUNMatrix A,
                       IndexType *output_time_index,
                       double *next_target_time);

  /// Removes pairs of very close interfaces that are beyond the limits of capillaries with roots.
  void RemoveBubblesFromOrigin(arma::Col<IndexType> &caps_with_roots,
                               arma::Col<IndexType> &caps_roots_directions,
                               arma::Col<IndexType> &changed_caps);

  /// Plug Capillaries witch interfaces jumped through inlet-nodes
  void PlugCapillariesWhoseInterfacesJumpedThroughInletNodes(
    arma::Col<IndexType> &caps_with_roots,
    arma::Col<IndexType> &caps_roots_directions,
    arma::Col<IndexType> &changed_caps);

  /// Initialise the IDA solver parameters and data structures.
  void InitialiseIDASolver(void **mem,
                           N_Vector *x,
                           N_Vector *dx,
                           SUNMatrix *A,
                           SUNLinearSolver *LS,
                           const double t0);

  /// Initialise the IDA state to solve until a change of state.
  void InitialiseIDAState(void *mem,
                          const N_Vector x,
                          const N_Vector dx,
                          const SUNMatrix A,
                          const SUNLinearSolver LS,
                          const double t0,
                          const double next_target_time,
                          const IndexType output_time_index);

  /// Computes the initial values of the DAE system.
  void ComputeInitialValues(const N_Vector x, const N_Vector dx);

  /// Computes the differential equation value of all links.
  void ComputeDifferentialEquationForAllCapillaries(arma::Col<double> &differential,
                                                    const arma::Col<double> &x_array);

  /// Computes the differential equation value of selected links.
  void ComputeDifferentialEquationForSelectedCapillaries(arma::Col<double> &differential,
                                                         const arma::Col<double> &x_array,
                                                         const arma::Col<IndexType> &capillaries);

  /// Computes the differential equation value of one link.
  double ComputeDifferentialEquation(const IndexType link_index, const arma::Col<double> &x_array);

  /// Set the identification of differential and algebraic equations.
  void SetIDValues(const N_Vector id);

  /// Returns true if the stop criteria is met.
  bool StopCriteriaMet(const double tFinal);

  /// Checks and stores changes of state returned by the IDA solver.
  bool CheckIDAState(const int state);

  /// Update the interface positions in array \c x with the appropriate values.
  void UpdateSundialsVariables(arma::Col<double> &x,
                               arma::Col<double> &dx,
                               const arma::Col<IndexType> &capillaries);

  /// Revert the previous state in \c previous_x_, \c previous_dx_, and \c previous_A_.
  void RevertState(arma::Col<double> &x, arma::Col<double> &dx, arma::Col<double> &A);

  /// Revert the pressures in \c pressures_ to a previous state.
  void RevertPressures(arma::Col<double> &x);

  /// Store the current state in \c previous_x_, \c previous_dx_, and \c previous_A_.
  void StoreState(const arma::Col<double> &x,
                  const arma::Col<double> &dx,
                  const arma::Col<double> &A);

  /// Store the current pressures in \c previous_interior_pressures_.
  void StoreInteriorPressures(const arma::Col<double> &x);

  /**
   * \brief Stores the number of equations of the system of DAEs
   *
   * The \c number_of_equations_ member stores the number of equations in the system of DAEs.
   * The number of equations is given by the number of capillaries in the capillary network
   * plus the number of non-boundary nodes (interior nodes).
   */
  IndexType number_of_equations_;

  /**
   * \brief Stores the current time of simulation
   *
   * The \c current_time member stores the time that is currently being computed by the simulator
   */
  double current_time_;

  /**
   * \brief Stores the resume flag of the simulation
   *
   * The \c resume member stores the boolean that indicates wheter or not to resume a simulation
   */
  bool resume_;

  /**
   * \brief Stores the time of the last interface jump in the simulation
   *
   * The \c last_jump_time_ member stores the time when the last interface jump occurred in the simulation
   */
  double last_jump_time_;

  /**
   * \brief This vector contains the pressures \f$ P_i \f$ at each node \f$ i \f$.
   *
   * The node pressures (in \c [Pa]) in the boundary is given by the boundary conditions
   * and the in the interiors are calculated by solving the transient model.
   */
  arma::Col<double> pressures_;

  /**
   * \brief This vector contains all x from IDA Sundials from the last state change.
   */
  arma::Col<double> previous_x_;

  /**
   * \brief This vector contains all dx from IDA Sundials from the last state change.
   */
  arma::Col<double> previous_dx_;

  /**
   * \brief This vector contains all A from IDA Sundials from the last state change.
   */
  arma::Col<double> previous_A_;

  /**
   * \brief This vector contains the pressures \f$ P_i \f$ at each interior node \f$ i \f$ at the
   * last state change.
   *
   * The node pressures (in \c [Pa]) in the in the interiors are calculated by solving the transient
   * model.
   */
  arma::Col<double> previous_interior_pressures_;

  /**
   * \brief Unique pointer to \c ModelBase instance defined by \c Initialise().
   *
   * Contains getters for \c DynamicCapillaryPhysicsBase and \c DynamicCapillaryGeometryBase objects.
   */
  std::unique_ptr<DynamicModelBase> model_;

  /**
   * \brief Unique pointer to \c StaticCapillaryExperimentBase instance defined by \c Initialise().
   *
   * This member variable stores a unique pointer to an object derived from
   * \c StaticCapillaryExperimentBase as set by the \c public method \c Initialise() according to
   * the \c name argument.
   */
  std::unique_ptr<DynamicCapillaryExperimentBase> experiment_;

  /**
   * \brief Shared pointer to \c CapillaryInterface instance defined by \c
   * SetCapillaryInterface().
   *
   * Contains methods needed by the algorithm
   */
  std::shared_ptr<CapillaryInterface> interface_;

  /**
   * \brief Unique pointer to \c DynamicCapillaryNetworkContext instance defined by \c
   * Initialise().
   *
   * Contains context variables needed by the classes in the Dynamic Capillary Network
   */
  std::shared_ptr<DynamicCapillaryNetworkContext> context_;

  /**
   * \brief Integer pointer used to recover the roots found in the \c IDASolve.
   */
  arma::Col<int> roots_found_;

  /**
   * \brief Vector with index of capillaries with interface but that are not plugged.
   */
  arma::Col<IndexType> caps_with_interfaces_not_plugged_;

  /**
   * \brief Vector with the limit values of capillaries that have interfaces and are not plugged.
   */
  arma::Col<double> caps_residual_limits_;

  /**
   * \brief Pointer to the Sundial memory space.
   */
  void *mem_;

  /**
   * \brief The initial simulation time.
   *
   * This floating-point variable stores the initial time (in seconds) at which the fluid flow
   * simulation starts.
   */
  double initial_time_;

  /**
   * \brief The final simulation time.
   *
   * This floating-point variable stores the final time (in seconds) at which the fluid flow
   * simulation ends.
   */
  double final_time_;

  /**
   * \brief Stores the simulation time steps size.
   */
  double time_step_size_;

  /**
   * \brief Stores the interval between saves of the resume file
   *u
   * The \c resume_save_frequency_  Stores the frequency between saves of the resume file. A resume
   * file wil be saved every resume_save_frequency_ roots found
   */
  int resume_save_frequency_ = 10;

  /**
   * \brief The relative tolerance used by IDA SUNDIALS as one of the stopping criteria in all
   * equations.
   */
  double relative_tolerance_;

  /**
   * \brief The absolute tolerance used by IDA SUNDIALS as one of the stopping criteria in partial
   * differential equations.
   */
  double absolute_link_tolerance_;

  /**
   * \brief The absolute tolerance used by IDA SUNDIALS as one of the stopping criteria in algebraic
   * equations.
   */
  double absolute_node_tolerance_;

  /**
   * \brief The number of roots found. Used to check wheter to save or not a resume file.
   * equations.
   */
  int number_of_roots_ = 0;
};  // end of class DynamicCapillaryNetworkAlgorithm

}  // namespace simulator

#endif  // SRC_FLOW_SIMULATOR_ALGORITHMS_DYNAMIC_CAPILLARY_NETWORK_DYNAMIC_CAPILLARY_NETWORK_H_
