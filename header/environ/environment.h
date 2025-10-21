#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <util/thrust_functors.h>
#include <map>

#include <libconfig.h++>

using namespace libconfig;

/**
 * @brief Simulates a habitat in which species go through their evolutional processes
 */
class environment
	{
	public:
		/**
		 * @brief Construct a new environment object
		 * 
		 * @param seed_val seed value for the random number generator
		 * @param num_biotic_variables number of biotic variables
		 * @param num_abiotic_variables number of abiotic variables
		 * @param num_demes number of demes
		 */
		environment(int seed_val, int num_biotic_variables, int num_abiotic_variables, int num_demes);
		/**
		 * Initializes the abiotic variables from a configuration file
		 * 
		 * @param filename configuration file name
		 */
		void initialize_abiotic_variables(const char *filename);
		/**
		 * @brief Destroy the environment object
		 * 
		 */
		~environment();

		/**
		 * @brief Random number generator for the habitat
		 * 
		 */
		curandGenerator_t gen;

		/// @brief Seed value for the random number generator
		int seed;
		/// @brief Number of demes
		int ndemes;
		/// @brief Number of biotic variables
		int nbiotic_vars;
		/// @brief Number of abiotic variables
		int nabiotic_vars;

		//void update(int num_biotic_variables, int month);
		/**
		 * Records the habitat details after the simulation into a file specified by its name
		 * 
		 * @note This function is meant to be provided with user-defined implementation.
		 * 
		 * @param output_file_environment output file name into which habitat details will be recorded
		 */
		void record_habitat_details(const char *output_file_environment);
		thrust::device_ptr<float> get_abiotic_vector_ptr(const char *abiotic_variable_name);


		// Data vectors - effect_of_inds_on_biotic_variables[prey][deme]
		/**
		 * Vector of lists of abiotic values, with the inner lists whose order relative to other inner lists corresponds to the abiotic names and whose size is of deme size
		 */
		thrust::device_vector<float> *biotic_variables;
		/**
		 * Vector of the lists of values which demonstrate the degrees to which each individual impact the habitat.
		 * 
		 * <p>The order of the inner lists is according to the order of the vector of abiotic factors' names, and the size of each inner list is the deme size, assuming that every individual of a deme exerts the same degree of impace on the habitat.</p>
		 */
		thrust::device_vector<float> *effect_of_inds_on_biotic_variable;

		/**
		 * Updates the habitat through each time step of the simulation
		 * 
		 * @note This is an abstract method which requires users to provide with its specific implementation.
		 */
		virtual void update() = 0;
		// Gradually move more of the other members into protected
	protected:
		/// @brief Vectors of names of abiotic variables
		std::vector<std::string> abiotic_variable_names;
		/// @brief Vector of lists of abiotic values, with the inner lists whose order relative to other inner lists corresponds to the abiotic names and whose size is of deme size
		thrust::device_vector<float> *abiotic_variables;
		/// @brief Map of abiotic variable names to their corresponding indices in the abiotic_variables vector
		std::map<std::string, int> abiotic_variable_indices;
	};
#endif
