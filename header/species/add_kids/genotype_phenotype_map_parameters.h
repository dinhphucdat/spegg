#ifndef GENOTYPE_PHENOTYPE_MAP_PARAMETERS_H
#define GENOTYPE_PHENOTYPE_MAP_PARAMETERS_H

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
 * @brief Collects parameters about genotypes and phenotypes from @c deme_config.txt
 * @author xxx
 */
class GenotypePhenotypeMapParameters
	{
	friend class DemeGeneticsSettings;
	/**
	 * @brief Construct a new Genotype Phenotype Map Parameters object
	 * 
	 * This constructor reads in the data from the @c deme_config.txt in terms of genotypes 
	 * and phenotypes by deme and stores it into its underlying data structures
	 * 
	 * @param filename config file name
	 * @param species_ID species ID, defined by users' @c inds object
	 * @param phenotype_index phenotype index, as specified in the config file
	 * @param loci_names names of the loci
	 */
	GenotypePhenotypeMapParameters(const char *filename, int species_ID, int phenotype_index, std::vector<std::string> &loci_names);

	public: 
		/// @brief vector of parameters, with each parameter having a smaller array containing genotype-phenotype information by deme
		thrust::device_vector<float> *deme_specific_parameters;
		/**
		 * @brief Get the subarray containing genotype-phenotype for that specified parameter for all demes
		 * 
		 * @param parameter_name name of parameter
		 * @return @c thrust::device_ptr<float> pointer of subarray for that parameter's information for all demes
		 */
		thrust::device_ptr<float> get_vector_ptr(const char *parameter_name);
	
	protected:
		/// @brief number of parameters
		int Number_of_Parameters; 
		/// @brief phenotype index
		int phenotype_index;
		/// @brief list of parameter names
		std::vector<std::string> Names_of_Genotype_Phenotype_Map_Parameters;
		/// @brief number of demes
		int Number_of_Demes;
		/**
		 * @brief Reads in the genotype-phenotype parameters from @c deme_config.txt
		 * 
		 * @param filename @c deme_config.txt
		 * @param species_ID ID of species, as specified in your @c inds object
		 */
		void read_in_data(const char *filename, int species_ID);
		/**
		 * @brief Adds the look-up table for parameter names' indices
		 * 
		 */
		void specify_parameter_index();
		/// @brief A look-up table for parameters' indices
		std::map<std::string, int> parameter_index;
	};

#endif
