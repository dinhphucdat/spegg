#ifndef DEME_SPECIFIC_DATA_CLASS_H
#define DEME_SPECIFIC_DATA_CLASS_H

#include <species/add_kids/genetic_deme_specific_data.h>
#include <vector>
#include <map>
#include <string>

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <libconfig.h++>
#include <string>
#include <curand.h>
#include <map>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <util/python_thrust_api.h>

using namespace libconfig;
namespace py = pybind11;
namespace StringVector = std::vector<std::string>;
namespace String2DVector = std::vector<std::vector<std::string>>;
namespace StringFloatMap = std::map<std::string, float>;
namespace FloatArrayVector = std::vector<py::array_t<float>>;

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

		// ----------------- NEW FUNCTIONALITY - PYTHON PARAM PASSING ----------------- //

		/**
		 * @brief Construct a new Deme Settings object. 
		 * 
		 * This constructor is used for a Python-binding program for direct data 
		 * transferring and will bypass the use of @c deme_config.txt
		 * 
		 * @param numDemes number of demes
		 * @param speciesID species ID
		 * @param parameterNames a Python list of parameter names
		 * @param demeWideParameters a numpy array storing float values of deme-wise parameters. Should be corresponding to the order of the list of parameter names
		 * @param speciesSpecificValues a Python dictionary of species specific values. Keys must be value names and values should be the corresponding float values
		 * @param phenotypeNames a Python list of phenotype names
		 * @param genPhenParameterNamesAllPhenotypes phenotype parameters for every phenotype. This is a 2D list, with the first dimension being the number of phenotypes, the second one the number of subparameters for each phenotype
		 * @param demeSpecificPhenParametersAllPhenotypes this should be a list of numpy's 2d arrays. 
		 * The outer dimension should have the size of number of phenotypes. The first dimension 
		 * of the inner numpy arrays should have the size of phenotype's specific parameters, 
		 * and the second dimension of the numpy arrays should have the size of number of demes.
		 * @param lociNames a Python list of loci names
		 * @param recombinationRates a Python numpy array of recombination rate. This is a 1D array, its size of the number of loci
		 * @param demeSpecificMutationRates a numpy array of deme-specific mutation rates. Dimension 1 is number of loci, dimension 2 is that locus's mutation rate for every deme
		 * @param demeSpecificMutationMagnitudes a numpy array of deme-specific mutation magnitudes. Dimension 1 is number of loci, dimension 2 is that locus's mutation magnitude for every deme
		 */
		DemeSettings(
			const int&									 numDemes,
			const int&									 speciesID, 
			const StringVector&    			 			 parameterNames, 
			const py::array_t<float>& 		 			 demeWideParameters, 
			const StringFloatMap& 				 		 speciesSpecificValues, 
			const StringVector& 			 			 phenotypeNames, 
			const String2DVector& 						 genPhenParameterNamesAllPhenotypes, 
			const FloatArrayVector& 		 			 demeSpecificPhenParametersAllPhenotypes
			const StringVector& 			 			 lociNames, 
			const py::array_t<float>& 					 recombinationRates, 
			const py::array_t<float>&  					 demeSpecificMutationRates, 
			const py::array_t<float>&					 demeSpecificMutationMagnitudes
		);

		// ----------------- END OF THIS NEW FUNCTIONALITY ---------------------------- //

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

		// ----------------- NEW FUNCTIONALITY - PYTHON PARAM PASSING ----------------- //

		/**
		 * @brief Process the deme-wise parameters
		 * 
		 * @param numDemes number of demes
		 * @param parameterNames a Python list of parameter names
		 * @param demeWideParameters the corresponding array of parameter values
		 */
		void processParameters(
			const int&						 numDemes,
			const StringVector&    			 parameterNames, 
			const py::array_t<float>& 		 demeWideParameters
		);

		/**
		 * @brief Process the deme specific values
		 * 
		 * @param speciesSpecificValues a Python dictionary of species specific values. Keys must be value names and values should be the corresponding float values
		 */
		void processDemeSpecificValues(const StringFloatMap& speciesSpecificValues);

		// ----------------- END OF THIS NEW FUNCTIONALITY ---------------------------- //
	};

#endif
