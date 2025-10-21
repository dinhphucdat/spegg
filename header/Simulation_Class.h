#ifndef SIM_H
#define SIM_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include <libconfig.h++>

/**
 * @brief Abstract base class for simulations.
 *
 * Simulation provides common configuration fields and helper methods used by
 * concrete simulation subclasses. Subclasses must implement the pure virtual
 * methods to provide simulation-specific initialization and execution.
 *
 * Note: method and field names must not be changed — this header only adds
 * documentation for extraction by Doxygen.
 */
class Simulation
{
public:
	/**
	 * @brief Construct a Simulation object.
	 *
	 * Concrete initialisation is typically performed in subclass constructors
	 * and `initialize_classes()`.
	 */
	Simulation();

	/**
	 * @brief Virtual destructor.
	 */
	~Simulation();

	/**
	 * @brief Run the simulation.
	 *
	 * This pure virtual method must be implemented by derived classes and
	 * should contain the primary execution loop or orchestration logic.
	 */
	virtual void run() = 0;

protected:
	// ---------------------------------------------------------------------
	// Input data fields
	// ---------------------------------------------------------------------

	/**
	 * @brief Identifier for the thread running this simulation (if multithreaded).
	 */
	int threadID;

	/**
	 * @brief Number of time steps to run the simulation for.
	 */
	int nsteps;

	/**
	 * @brief Number of demes (sub-populations) in the simulation.
	 */
	int demes;

	/**
	 * @brief Initial population size.
	 */
	int initpop;

	/**
	 * @brief Maximum allowed population size.
	 */
	int maxpop;

	/**
	 * @brief Number of intra-step time steps (sub-steps within each main step).
	 */
	int intra_step_time_steps;

	/**
	 * @brief Number of genetic loci modelled.
	 */
	int nloci;

	/**
	 * @brief Number of phenotypes modelled.
	 */
	int nphenotypes;

	/**
	 * @brief RNG seed used to initialise random number generators.
	 */
	int seed;

	/**
	 * @brief Number of biotic variables tracked per individual or deme.
	 */
	int num_biotic_variables;

	/**
	 * @brief Number of abiotic variables tracked in the environment.
	 */
	int num_abiotic_variables;

	/**
	 * @brief Steps at which per-individual CSV output should be written.
	 *
	 * Each element is a time-step index where `individuals.csv` (or similar)
	 * output will be produced.
	 */
	std::vector<int> steps_to_output_individuals_csv;

	// ---------------------------------------------------------------------
	// Methods
	// ---------------------------------------------------------------------

	/**
	 * @brief Read simulation settings from configuration (e.g., libconfig).
	 *
	 * This method reads configuration values and populates the protected
	 * fields above. It is intended to be called by subclasses during
	 * initialization.
	 */
	void read_simulation_settings();

	/**
	 * @brief Initialize simulation-specific classes and data structures.
	 *
	 * Derived classes must implement this to allocate/setup species, demes,
	 * environment objects and any other simulation-specific resources.
	 */
	virtual void initialize_classes() = 0;
};

#endif
