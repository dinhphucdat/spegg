#ifndef INDS_H
#define INDS_H

#include <vector>
#include <map>
#include <string>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/gather.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <util/python_thrust_api.h>

#include <species/deme_specific_data_class.h>
#include <environ/environment.h>

namespace py = pybind11;
using StringVector = std::vector<std::string>;
using String2DVector = std::vector<std::vector<std::string>>;
using StringFloatMap = std::map<std::string, float>;
using FloatArrayVector = std::vector<py::array_t<float>>;

/*!
 *
 * A base class for storing the individuals belonging to a specific species, and performing 
 * three operations that are common to all species: removing the dead individuals (which is not 
 * implemented as a stream compaction for performanc reasons), sorting individuals according to 
 * their deme, and writing the output to a CSV file stored on the hard drive. The methods and 
 * data structures found in this class are those which apply to all species sPEGG can simulate.
 *
 * The basic philosophy behind inds is that the collection of individuals simulated for each species 
 * can be thought of as an N x M matrix, with N representing the number of individuals and M representing 
 * the number of attributes of these individuals (e.g., their genotypes, their ID, their demes, whether 
 * they are dead or alive, etc...). (MOVE TO QUICKSTART) By convention, every individual data point that 
 * is not a genotype, an ID, whether they are dead or alive, their sex, their age, and the IDs of their 
 * parents is designated as a "phenotype".
 *
 */
class inds
	{
	public:
		friend class Statistics;
		friend class Parents;
		friend class EggsNeonates;
		/**
		 * Instantiates an @c inds object with its related information specified in the @c deme_config.txt file.
		 * <br>
		 * <p>See <strong>deme_config.txt guide</strong> for more information.</p>
		 * <p>
		 * The program will terminate if:
		 * <ul>
		 *   <li>@c num_demes defined in @c Simulation.conf (essentially the value passed into this constructor) 
		 * exceeds the value declared in @c deme_config.txt.</li>
		 *   <li>@c size_val or @c maxsize_val are non-positive.</li>
		 *   <li>@c maxsize_val is less than @c size_val.</li>
		 * </ul>
		 * Ensure these constraints are satisfied before using this function.
		 * </p>
		 * @param size_val the initial size of the population of the species specified
		 * @param maxsize_val the maximal capacity that can hold individuals of that species
		 * @param num_demes the number of demes into which individuals of that species are split
		 * @param species_ID_val the ID with which that specific species is tagged
		 */
		inds(int size_val, int maxsize_val, int num_demes, int species_ID_val);

		/**
		 * @brief Construct a new inds object.
		 * 
		 * This operation is designed to let users directly pass the parameters in 
		 * without having to create any intermediary file (i.e. @c deme_config.txt ), 
		 * thus reduces the overhead of making too many IO request if this project is 
		 * going to be optimized with any machine learning method.
		 * 
		 * @param size_val the initial size of the population of the species specified
		 * @param maxsize_val the maximal capacity that can hold individuals of that species
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
		inds(
			int size_val, 
			int maxsize_val, 
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
		 * @brief Destroy the inds object
		 * 
		 */
		~inds();
		/**
		 * A function to reorganize your data points so that inds (0, 1, 2, ..., size) consist only of individuals whose vital status = 1.
		 * 
		 * <p>
		 * Rearranges the data structures of @c inds so that all the data points represented in each data structure 
		 * from 0 to number_of_individuals represent individuals that are alive. The data points past 
		 * number_of_individuals are garbaged. All operations on inds or its derived classes should, 
		 * therefore, operate only on @c number_of_individuals data points. 
		 * As a rule, @c number_of_individuals <= @c max_number_of_individuals; when the two are equal, 
		 * the data from dead individuals should simply be overwritten 
		 * until @c number_of_individuals < @c max_number_of_individuals again.
		 * </p>
		 * 
		 * Ideally this would involve stream compaction rather than relegating the dead individuals to occupy empty spaces, but some preliminary experiments suggested stream compaction entails substantial performance costs in thrust compared to the approach using sorting and gathering. This solution is subject to change; should performance improvements in Thrust or CUDA allow it, we will return to stream compaction routines here instead.
		 *
		 */
		void removeDead();
		/**
		 * @brief Set the Max Size object
		 * 
		 * @param n new maximum capacity of the population
		 */
		void setMaxSize(int n);
		/**
		 * @brief Sort individuals by their deme
		 * 
		 */
		void sortByDeme();
		/**
		 * Fills in the vectors with specific data about every individual from a CSV file.
		 * @note
		 * <ul>
		 * 	<li>If this object is not instantiated from @c inds_stochastic , a child of @c inds class, you may 
		 * be required to call this function in order to provide the program with individual-oriented information.</li>
		 * 	<li>If the row count, understood as initial population size, exceeds @c maxsize_val specified in the 
		 * constructor, this program will be terminated.</li>
		 * </ul>
		 * 
		 * <p>
		 * The CSV file may look like:
		 * <br>
		 * 
		 * <blockquote>
		 * @code{.csv}
		 * Index,Id,Status,Sex,Age,Deme,fgene0,fgene1,fgene2,mgene0,mgene1,mgene2,phen0,phen1,phen2
		 * 0,0,1,0,0,0,-0.0752043,-0.0743299,0.195556,-0.0759161,0.0367624,-0.0496664,9.2444,0.5,0.468834
		 * 1,1,1,1,0,0,0.000129433,-0.0521384,-0.228273,0.0881785,-0.00500755,-0.0091536,10.4415,0.5,-0.10614
		 * 2,2,1,0,0,0,0.153293,-0.111343,0.0840508,-0.0710807,0.157572,0.132522,10.4111,0.5,0.574859
		 * 3,3,1,1,0,0,0.0588492,0.0373746,0.0388453,-0.0822977,0.108823,-0.156742,9.88276,0.5,0.0731549
		 * 4,4,1,0,0,0,-0.0719927,-0.0156143,-0.0481384,-0.176604,-0.00593068,0.0552063,8.75702,0.5,0.260602
		 * 5,5,1,1,0,0,0.0448649,-0.0873122,-0.150351,0.167369,-0.144792,-0.140238,11.0612,0.5,-0.185884
		 * ...
		 * @endcode
		 * </blockquote>
		 * 
		 * @c maternal_ID and @c paternal_ID for each individual can be provided as well, these are optional.
		 * 
		 * <br>
		 * 
		 * See more at the <strong>CSV guide</strong>.
		 * </p>
		 * 
		 * @param filename name of the CSV file where initial values for specific individuals are stored
		 */
		virtual void initialize_from_CSV(const char *filename);
		/**
		 * Initializes the vectors storing information about @c fgenotype, @c mgenotype, and @c phenotype.
		 * 
		 * <p>This function allocates memory for the vectors with sizes determined by the specified 
		 * number of loci and phenotypes. It does not populate the vectors with data, serving only to 
		 * prepare them for subsequent use.</p>
		 * 
		 * @note
		 * This function is intended as a helper method and should not be overridden.
		 * 
		 * @param nloci_val number of loci for which genotypic data are stored, already provided in @c deme_config.txt
		 * @param nphen_val number of phenotypes, also already provided in @c deme_config.txt
		 */
		void initialize_individuals(int nloci_val, int nphen_val);
		/**
		 * @brief Exports the individuals' data to a CSV file
		 * 
		 */
		void exportCsv();
		/**
		 * @brief Exports the individuals' data to a CSV file with a specified filename
		 * 
		 * @param filename name of the CSV file to export to
		 */
		void exportCsv(const char *filename);
		/**
		 * @brief Exports the individuals' data to a CSV file with a specified timestep
		 * 
		 * @param timestep timestep to include in the filename
		 */
		void exportCsv(int timestep);
		/**
		 * @brief exports the individuals' data to a CSV file with a specified filename and timestep
		 * 
		 * @param filename name of the CSV file to export to
		 * @param timestep timestep to include in the filename
		 */
		void exportCsv(const char *filename, int timestep);
		/**
		 * @brief exports the individuals' data to a CSV file with a specified filename and timestep range
		 * 
		 * @param filename name of the CSV file to export to
		 * @param timestep1 starting timestep
		 * @param timestep2 ending timestep
		 */
		void exportCsv(const char *filename, int timestep1, int timestep2);

		// Input parameters
		DemeSettings *demeParameters;

		//Misc data ints
		/**
		 * size of the population of the species
		 */
		int size;
		/**
		 * maximum capacity of the population
		 */
		int maxsize;
		/**
		 * number of phenotypes defined in the species
		 */
		int nphen;
		/**
		 * number of loci (genes) in the species' genomes
		 */
		int nloci;
		/**
		 * next id pointing to the next individual, this is for iterative purpose
		 */
		int nextid;
		/**
		 * number of demes (clusters) into which the population is split
		 */
		int Num_Demes;
		/**
		 * ID for this species
		 */
		int species_ID;

		//Data vectors
		/**
		 * list of IDs for all individuals
		 */
		thrust::device_vector<int> id;
		/**
		 * list of statuses
		 */
		thrust::device_vector<int> status;
		/**
		 * list of sexes
		 */
		thrust::device_vector<int> sex;
		/**
		 * list of ages
		 */
		thrust::device_vector<int> age;
		/**
		 * list of the deme number that each individual belongs to
		 */
		thrust::device_vector<int> deme;
		/**
		 * list of the sets of all alleles inherited from each individual's dad.
		 * This is a list of lists, and is dynamically allocated.
		 */
		thrust::device_vector<float> *fgenotype;
		/**
		 * list of the sets of all alleles inherited from each individual's mom.
		 * This is a list of lists, and is dynamically allocated.
		 */
		thrust::device_vector<float> *mgenotype;
		/**
		 * list of the sets of all phenotypes expressed in each individual.
		 * This is a list of lists, and is dynamically allocated.
		 */
		thrust::device_vector<float> *phenotype;
		/**
		 * ids of each individual's mom
		 */
		thrust::device_vector<int> maternal_id;
		/**
		 * ids of each individual's dad
		 */
		thrust::device_vector<int> paternal_id;
		/**
		 * list of deme sizes. Demes are listed in order (ie. deme 0 will be at position 0, and so on) 
		 * with their corresponding number of individuals in that deme.
		 * 
		 * @note This does not belong to the population matrix, only describes the state of each deme.
		 */
		thrust::device_vector<int> deme_sizes;
		/**
		 * list of maximum sizes for all demes. Also not belongs to the population matrix as @link deme_sizes deme_sizes @endlink.
		 */
		thrust::device_vector<int> max_deme_sizes;

		// Calculate the number of individuals in each deme
		void demeCalculations();
	};

/**
 * @brief A functor to reassign the deme of dead individuals
 * 
 */
struct reassign_dead_deme_functor
	{
	/*!
	*
	* A functor to be called from Thrust to reassign the deme of individuals with status = 0 (i.e., dead individuals) to be equal to the maximum number of demes + 1. This is used because the way in which the dead are removed in sPEGG is actually through a sorting regime, rather than stream compaction for performance reasons. For more details, see also \link ../../src/cuda/species/inds.cu.
	*
	*/

	/// @brief number of demes
	int num_demes;
	/**
	 * @brief Construct a new reassign dead deme functor object
	 * 
	 * @param number_demes number of demes
	 */
	reassign_dead_deme_functor(int number_demes) : num_demes(number_demes)
	{};
	/*
		Elements in the tuple.
		----------------------
		0: status
		1: deme
	*/

	/**
	 * @brief Call operator to reassign the deme of dead individuals
	 * 
	 * @tparam tuple 
	 * @param t tuples containing: 0: status, 1: deme
	 * @return __host__ 
	 */
	template <typename tuple>
	__host__ __device__ 
	void operator() ( tuple t ) {
		if (thrust::get<0>(t)==0) /* if the individual is dead */
			{
			thrust::get<1>(t) = num_demes;
			}
		}
	};


#endif
