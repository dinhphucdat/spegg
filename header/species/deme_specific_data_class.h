#ifndef DEME_SPECIFIC_DATA_CLASS_H
#define DEME_SPECIFIC_DATA_CLASS_H

#include <species/add_kids/genetic_deme_specific_data.h>

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <libconfig.h++>
#include <string>
#include <curand.h>
#include <map>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <algorithm>

using namespace libconfig;

/**
 * @brief A class to read in and store deme-specific parameters
 * 
 */
class DemeSettings
	{
	public: 
		/**
		 * @brief Construct a new Deme Settings object
		 * 
		 * @param filename @c deme_config.txt
		 * @param species_ID ID of species, as specified in your @c inds object
		 */
		DemeSettings(const char *filename, int species_ID);
		/// @brief vector of parameters, with each parameter having a smaller array containing deme-specific information by deme
		thrust::device_vector<float> *deme_wide_parameters;
		/// @brief parameter name to index map
		std::map<std::string, int> parameter_index;
		/// @brief species-specific value name to value map
		std::map<std::string, float> species_specific_values;
		/// @brief pointer to genetic architecture object
		DemeGeneticsSettings *GeneticArchitecture;
		/// @brief recombination rates by loci, deme-specific
		thrust::device_vector<float> *deme_specific_recombination_rates;
		/**
		 * @brief Get the vector ptr object
		 * 
		 * @param parameter_name name of parameter
		 * @return thrust::device_ptr<float> pointer of subarray for that parameter's information for all demes
		 */
		thrust::device_ptr<float> get_vector_ptr(const char *parameter_name);
		/**
		 * @brief Get the allelic effects ptr object
		 * 
		 * @param locus_number locus number
		 * @return thrust::device_ptr<float> 
		 */
		thrust::device_ptr<float> get_allelic_effects_ptr(int locus_number);
		/**
		 * @brief Query if a parameter exists
		 * 
		 * @param parameter_name name of parameter being queried
		 * @return true if the parameter exists
		 * @return false if not
		 */
		bool does_parameter_exist(const char *parameter_name);
		/**
		 * @brief Returns the number of demes specified in the configuration file
		 * 
		 * @return int number of demes
		 */
		int check_number_of_demes();
	protected:
		/**
		 * @brief Specifies the parameter index for quick look-up
		 * 
		 */
		void specify_parameter_index();
		/**
		 * @brief Reads in the deme-specific parameters from @c deme_config.txt
		 * 
		 * @param filename @c deme_config.txt
		 * @param species_ID ID of species, as specified in your @c inds object
		 */
		void read_in_parameters(const char *filename, int species_ID);
		/// @brief Number of demes
		int Number_of_Demes;
		/// @brief Number of parameters
		int Number_of_Parameters;
		/// @brief Number of species-specific values
		int Number_of_Species_Specific_Values;
		/// @brief names of parameters
		std::vector<std::string> parameter_names;	
		/// @brief names of species-specific values
		std::vector<std::string> species_specific_values_names;
	};

#endif
