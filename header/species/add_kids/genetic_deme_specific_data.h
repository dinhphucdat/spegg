#ifndef GENETIC_DEME_SPECIFIC_DATA_CLASS_H
#define GENETIC_DEME_SPECIFIC_DATA_CLASS_H

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <libconfig.h++>
#include <string>
#include <curand.h>
#include <map>
#include <thrust/device_vector.h>
#include <thrust/functional.h>

#include <species/add_kids/genotype_phenotype_map_parameters.h>

using namespace libconfig;

/**
 * @brief A class to read in and store genetic parameters that are deme-specific
 * 
 */
class DemeGeneticsSettings
	{
	public: 
		/**
		 * @brief Construct a new Deme Genetics Settings object
		 * 
		 * This constructor reads in the data from the @c deme_config.txt in terms of genetic parameters by deme and stores it into its underlying data structures
		 * 
		 * @param filename config file name
		 * @param species_ID species ID, defined by users' @c inds object
		 */
		DemeGeneticsSettings(const char *filename, int species_ID);
		/// @brief recombination rates by loci
		thrust::device_vector<float> recombination_rates;
		/// @brief genphenotype map parameters by phenotype
		GenotypePhenotypeMapParameters **phen_gen_map_parm;
		/**
		 * @brief Get the vector ptr object
		 * 
		 * @param locus_number locus number
		 * @return thrust::device_ptr<float> 
		 */
		thrust::device_ptr<float> get_vector_ptr(int locus_number);
		/// @brief phenotype names vector
		std::vector<std::string> phenotype_names;
		/**
		 * @brief Get the mutation rates ptr object
		 * 
		 * @param locus_index locus index
		 * @return thrust::device_ptr<float> 
		 */
		thrust::device_ptr<float> get_mutation_rates_ptr(int locus_index);
		/**
		 * @brief Get the mutation magnitudes ptr object
		 * 
		 * @param locus_index locus index
		 * @return thrust::device_ptr<float> 
		 */
		thrust::device_ptr<float> get_mutation_magnitudes_ptr(int locus_index);
		/// @brief number of loci
		int Number_of_Loci;
		/// @brief number of phenotypes
		int Number_of_Phenotypes;

	protected:
		/**
		 * @brief reads in the deme-specific genetic data from @c deme_config.txt
		 * 
		 * @param filename @c deme_config.txt
		 * @param species_ID id of the species, as specified in your @c inds object
		 */
		void read_in_data(const char *filename, int species_ID);
		/// @brief number of demes
		int Number_of_Demes;
		/// @brief mutation rates by loci
		thrust::device_vector<float> *deme_specific_mutation_rates;
		/// @brief mutation magnitudes by loci
		thrust::device_vector<float> *deme_specific_mutation_magnitudes;
		/// @brief loci names vector
		std::vector<std::string> loci_names;
	};

#endif
