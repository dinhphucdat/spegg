#ifndef SAMPLING_INPUT_H
#define SAMPLING_INPUT_H

#include <curand.h>
#include <iostream>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>
#include <thrust/remove.h>

/* The class sampling input calculates and provides:
	an int vector of individuals potentially doing the sampling. 
	an int vector of individuals that will potentially be subject to sampling.
	an int vector of the number of times an individual doing the sampling will sample other individuals.
*/

// TODO: Reduce promiscuity in data availability; there should be a get_vector_ptr() operation to access/interface with the data members

/**
 * @brief An interface for sampling input
 * 
 * The class sampling input calculates and provides:
 * an int vector of individuals potentially doing the sampling. 
 * an int vector of individuals that will potentially be subject to sampling.
 * an int vector of the number of times an individual doing the sampling will sample other individuals.
 * 
 */
class SamplingInput{
/* 
*The class sampling input calculates and provides:
* an int vector of individuals potentially doing the sampling. 
* an int vector of individuals that will potentially be subject to sampling.
* an int vector of the number of times an individual doing the sampling will sample other individuals.
*/
	public:	
	// outputs
		/// @brief list of individuals potentially conducting sampling
		thrust::device_vector<int> list_of_individuals_potentially_conducting_sampling;
		/// @brief list of individuals potentially subject to sampling
		thrust::device_vector<int> list_of_individuals_potentially_subject_to_sampling;
		/// @brief number of other individuals sampled by each sampling individual
		thrust::device_vector<int> number_of_other_individuals_sampled;
		/// @brief deme affiliation of sampling individuals
		thrust::device_vector<int> deme_affiliation_of_sampling_individuals;
		/// @brief sampleable individuals per deme
		thrust::device_vector<int> sampleable_individuals_per_deme;
		/// @brief focal species ID
		int focal_species_ID;
		/// @brief target species ID
		int target_species_ID;
		/**
		 * @brief Sets up the sampling input by specifying the individuals conducting sampling, 
		 * individuals being sampled, and the number of individuals being sampled
		 * 
		 * @param Num_Demes number of demes
		 * @param sampling_individuals_demes deme affiliation of sampling individuals
		 */
		void setup_sampling_individuals_demes(int Num_Demes, thrust::device_vector<int> &sampling_individuals_demes);
		/**
		 * @brief Sets up the sampleable individuals per deme
		 * 
		 * @param num_sampleable_individuals_in_demes number of sampleable individuals in demes
		 */
		void setup_sampleable_individuals_per_deme(thrust::device_vector<int> &num_sampleable_individuals_in_demes);
		/// @brief sampling scheme
		std::string sampling_scheme;
		/// @brief number of demes
		int Num_Demes;
	};
/**
 * @brief A derived class of @link SamplingInput @endlink that provides sampling input for 
 * inter-species interactions
 * 
 */
class SamplingInput_Interactions : public SamplingInput
	{
	/* 
	* Sampling input for sampling individuals from more than one species. This is mothballed as this solution doesn't have attractive performance characteristics.
	*/	
	public:
		// Factory to create appropriate subclass
		/**
		 * @brief Create a SamplingInput object. This serves as a factory method 
		 * to create derived classes of @link SamplingInput @endlink
		 * 
		 * @param species double pointer representing all species
		 * @param focal_species focal species
		 * @param alternative_species alternative species
		 * @return SamplingInput_Interactions* 
		 */
		static SamplingInput_Interactions *create_SamplingInput(class Species **species, int focal_species, int alternative_species);	
		/**
		 * @brief Construct a new SamplingInput_Interactions object
		 * 
		 */
		void prepare_sampling_input();
		/// @brief total population size of the focal species
		int focal_species_total_pop_size;

	protected:
		class Species **species;

		// data members

		int number_of_individuals_in_sampling_species;
		int number_of_individuals_in_species_being_sampled;
		int number_sampled;

		void set_species_IDs( int focal_species, int alternative_species);
		void set_parameters();

		// methods that vary across different sampling specifications
		virtual void specify_individuals_conducting_sampling() = 0;
		virtual void specify_individuals_being_sampled() = 0;
		virtual void determine_number_of_individuals_sampled() = 0;

		// Stock algorithms for determining the number of individuals sampled
		void determine_number_of_individuals_to_be_sampled_fixed(); 
		void determine_number_of_individuals_to_be_sampled_poisson(); 
	};		
#endif
