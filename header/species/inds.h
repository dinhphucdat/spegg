#ifndef INDS_H
#define INDS_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/gather.h>

#include <species/deme_specific_data_class.h>
#include <environ/environment.h>
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

		inds(int size_val, int maxsize_val, int num_demes, int species_ID_val);
		~inds();
	
		void removeDead();
		void setMaxSize(int n);
		void sortByDeme();
		virtual void initialize_from_CSV(const char *filename);
		void initialize_individuals(int nloci_val, int nphen_val);
		void exportCsv();
		void exportCsv(const char *filename);
		void exportCsv(int timestep);
		void exportCsv(const char *filename, int timestep);
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

struct reassign_dead_deme_functor
	{
	/*!
	*
	* A functor to be called from Thrust to reassign the deme of individuals with status = 0 (i.e., dead individuals) to be equal to the maximum number of demes + 1. This is used because the way in which the dead are removed in sPEGG is actually through a sorting regime, rather than stream compaction for performance reasons. For more details, see also \link ../../src/cuda/species/inds.cu.
	*
	*/

	int num_demes;
	reassign_dead_deme_functor(int number_demes) : num_demes(number_demes)
	{};
	/*
		Elements in the tuple.
		----------------------
		0: status
		1: deme
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
