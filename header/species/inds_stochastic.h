#ifndef STOCHASTIC_INDS_H
#define STOCHASTIC_INDS_H

#include <string>
#include <vector>
#include <map>

#include <species/inds.h>
#include <math/random_variables_functions.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <util/python_thrust_api.h>

namespace py = pybind11;
namespace StringVector = std::vector<std::string>;
namespace String2DVector = std::vector<std::vector<std::string>>;
namespace StringFloatMap = std::map<std::string, float>;
namespace FloatArrayVector = std::vector<py::array_t<float>>;

/**
 * This class is a subset of @link inds inds @endlink, used when the initial population's attributes are 
 * desired to be normally distributed across the individuals.
 */
class inds_stochastic : public inds
	{
	public:
		/**
		 * @c inds_stochastic class, a subset of @c inds class, if the randomness of the initial population is desired.
		 * 
		 * <p>This constructor will prime up the random generator that will be used throughout the program for 
		 * multiple purposes such as generate random initial pop information and other things as well.
		 * </p>
		 * @param size_val the initial size of the population of the species specified
		 * @param maxsize_val the maximal capacity that can hold individuals of that species
		 * @param seed_val the seed value of the random generator
		 * @param num_demes the number of demes into which individuals of that species are split
		 * @param species_ID_val the ID with which that specific species is tagged
		 * @sa inds::inds for additional notes.
		 */
		inds_stochastic(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);

		/**
		 * @brief Construct a new inds stochastic object
		 * 
		 * This operation is designed to let users directly pass the parameters in 
		 * without having to create any intermediary file (i.e. @c deme_config.txt ), 
		 * thus reduces the overhead of making too many IO request if this project is 
		 * going to be optimized with any machine learning method.
		 * 
		 * @param size_val the initial size of the population of the species specified
		 * @param maxsize_val the maximal capacity that can hold individuals of that species
		 * @param seed_val the seed value of the random generator
		 * @param num_demes the number of demes into which individuals of that species are split
		 * @param species_ID_val the ID with which that specific species is tagged@param size_val 
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
		inds_stochastic(
			int size_val, 
			int maxsize_val, 
			int seed_val, 
			int num_demes, 
			int species_ID_val, 
			const StringVector&    		parameterNames, 
			const py::array_t<float>& 	demeWideParameters, 
			const StringFloatMap& 		speciesSpecificValues, 
			const StringVector& 		phenotypeNames, 
			const String2DVector& 		genPhenParameterNamesAllPhenotypes, 
			const FloatArrayVector& 	demeSpecificPhenParametersAllPhenotypes, 
			const StringVector& 		lociNames, 
			const py::array_t<float>& 	recombinationRates, 
			const py::array_t<float>&  	demeSpecificMutationRates, 
			const py::array_t<float>&	demeSpecificMutationMagnitudes
		);

		/**
		 * Main generator of the class
		 */
		curandGenerator_t gen;
		/**
		 * How to update the state of the population. <strong>User-defined function</strong>.
		 */
		virtual void update(inds_stochastic **species) {};
		/**
		 * How to update the state of the population. <strong>User-defined function</strong>.
		 */
		virtual void update(inds_stochastic **species, environment *habitat) {};
		/**
		 * How to update the state of the population. <strong>User-defined function</strong>.
		 */
		virtual void update(inds_stochastic **species, environment *habitat, int intra_step_time_steps) {};
		/**
		 * How to update the state of the population. <strong>User-defined function</strong>.
		 */
		virtual void update(inds_stochastic **species, environment *habitat, int intra_step_time_steps, int current_time_step) {};
		/**
		 * How to add offspring into the population. <strong>User-defined function</strong>.
		 */
		virtual void addKids() {};
		/**
		 * How to add offspring into the population. <strong>User-defined function</strong>.
		 */
		virtual void addKids(environment *habitat) {};
		/// @brief seed value for the random number generator
		int seed;
		/// @brief phenotype index that determines mortality
		int MORTALITY_PHENOTYPE_INDEX;
		/// @brief phenotype index that determines fecundity
		int FECUNDITY_PHENOTYPE_INDEX;
	};

#endif
