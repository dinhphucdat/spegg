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
 * Simulates a habitat in which species go through their evolutional processes
 */
class environment
	{
	public:
		environment(int seed_val, int num_biotic_variables, int num_abiotic_variables, int num_demes);
		void initialize_abiotic_variables(const char *filename);
		~environment();

		//Random number generator
		curandGenerator_t gen;

		//Misc data ints
		int seed;
		int ndemes;
		int nbiotic_vars;
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
		std::vector<std::string> abiotic_variable_names;
		thrust::device_vector<float> *abiotic_variables;

		std::map<std::string, int> abiotic_variable_indices;
	};
#endif
