#ifndef GENETIC_DEME_SPECIFIC_DATA_CLASS_H
#define GENETIC_DEME_SPECIFIC_DATA_CLASS_H

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <libconfig.h++>
#include <string>
#include <curand.h>
#include <map>
#include <vector>
#include <thrust/device_vector.h>
#include <thrust/functional.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <util/python_thrust_api.h>

#include <species/add_kids/genotype_phenotype_map_parameters.h>

using namespace libconfig;
namespace py = pybind11;
using StringVector = std::vector<std::string>;
using String2DVector = std::vector<std::vector<std::string>>;
using FloatArrayVector = std::vector<py::array_t<float>>;

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
		/// @brief Destroys the object and its associated dynamically allocated objects
		~DemeGeneticsSettings();

		// -------------- NEW FUNCTIONALITY - READ FROM PYTHON ------------- //

		/**
		 * @brief Construct a new Deme Genetics Settings object.
		 * 
		 * This operation is designed to let users directly pass the parameters in 
		 * without having to create any intermediary file (i.e. @c deme_config.txt ), 
		 * thus reduces the overhead of making too many IO request if this project is 
		 * going to be optimized with any machine learning method.
		 * 
		 * @param speciesID species ID
		 * @param phenotypeNames a Python list of phenotype names
		 * @param genPhenParameterNamesAllPhenotypes phenotype parameters for every phenotype. This is a 2D list, with the first dimension being the number of phenotypes, the second one the number of subparameters for each phenotype
		 * @param demeSpecificPhenParametersAllPhenotypes a Python list of 2D numpy arrays. The list should have the size of the number of phenotypes, and for each numpy array, dimension 1 should be subparameters, and dimension 2 should be the value of that specific subparameter for every deme
		 * @param lociNames a Python list of loci names
		 * @param recombinationRates a Python numpy array of recombination rate. This is a 1D array, its size of the number of loci
		 * @param demeSpecificMutationRates a numpy array of deme-specific mutation rates. Dimension 1 is number of loci, dimension 2 is that locus's mutation rate for every deme
		 * @param demeSpecificMutationMagnitudes a numpy array of deme-specific mutation magnitudes. Dimension 1 is number of loci, dimension 2 is that locus's mutation magnitude for every deme
		 */
		DemeGeneticsSettings(
			const int&									 speciesID, 
			const StringVector& 			 		 	 phenotypeNames, 
			const String2DVector& 						 genPhenParameterNamesAllPhenotypes, 
			const FloatArrayVector& 		 			 demeSpecificPhenParametersAllPhenotypes,
			const StringVector& 			 			 lociNames, 
			const py::array_t<float>& 					 recombinationRates, 
			const py::array_t<float>&  					 demeSpecificMutationRates, 
			const py::array_t<float>&					 demeSpecificMutationMagnitudes
		);

		// -------------- END OF THIS NEW FUNCTIONALITY -------------------- //

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

		// -------------- NEW FUNCTIONALITY - READ FROM PYTHON ------------- //

		/**
		 * @brief Process the genotypic information.
		 * 
		 * A part of the new functionality, that is directly read parameters from Python interface
		 * 
		 * @param lociNames a Python list of loci names
		 * @param recombinationRates a Python numpy array of recombination rate. This is a 1D array, its size of the number of loci
		 * @param demeSpecificMutationRates a numpy array of deme-specific mutation rates. Dimension 1 is number of loci, dimension 2 is that locus's mutation rate for every deme
		 * @param demeSpecificMutationMagnitudes a numpy array of deme-specific mutation magnitudes. Dimension 1 is number of loci, dimension 2 is that locus's mutation magnitude for every deme
		 */
		void processGenotypicInfo(
			const StringVector&  			 lociNames, 
			const py::array_t<float>&  	 	 recombinationRates, 
			const py::array_t<float>& 	 	 demeSpecificMutationRates, 
			const py::array_t<float>&		 demeSpecificMutationMagnitudes
		);

		/**
		 * @brief Process the phenotypic information. 
		 * 
		 * A part of the new functionality, that is directly read parameters from Python interface
		 * 
		 * @param speciesID species ID
		 * @param phenotypeNames a Python list of phenotype names
		 * @param genPhenParameterNamesAllPhenotypes phenotype parameters for every phenotype. This is a 2D list, with the first dimension being the number of phenotypes, the second one the number of subparameters for each phenotype
		 * @param demeSpecificPhenParametersAllPhenotypes a Python list of 2D numpy arrays. The list should have the size of the number of phenotypes, and for each numpy array, dimension 1 should be subparameters, and dimension 2 should be the value of that specific subparameter for every deme
		 */
		void processPhenotypicInfo(
			const int& 				speciesID, 
			const StringVector& 	phenotypeNames, 
			const String2DVector& 	genPhenParameterNamesAllPhenotypes, 
			const FloatArrayVector& demeSpecificPhenParametersAllPhenotypes
		);

		// -------------- END OF THIS NEW FUNCTIONALITY -------------------- //
	};

#endif
