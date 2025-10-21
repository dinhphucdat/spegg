#ifndef SAMPLING_INPUT_MATING_H
#define SAMPLING_INPUT_MATING_H

#include <iostream>
#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>
#include <util/sampling_input.h>
#include <species/add_kids/assortative_mating_parents_class.h>

/**
 * @brief A derived class of @link SamplingInput @endlink that provides sampling input for mating purposes
 * 
 */
class SamplingInput_Mating : public SamplingInput
	{
	/*
	* A class whereby females or males of a given species will "sample" individuals of the opposite sex and evaluate their suitability as mates. A common example for a situation where this class would be necessary would be when female preference varies across individual females, who must evaluate a sequence of males and select a mate according to their innate preference function.
	*/
	public:
		/**
		 * @brief Construct a new SamplingInput_Mating object
		 * 
		 * @param mating_parents @c Assortative_mating_parents object representing the parents of a species
		 * @param Sampling_Parent 0: females sample males; 1: males sample females
		 */
		SamplingInput_Mating(Assortative_mating_parents *mating_parents, int Sampling_Parent);
		/**
		 * @brief Determines the number of individuals that will be sampled by each sampling individual
		 * 
		 * @param mating_parents @c Assortative_mating_parents object representing the parents of a species
		 */
		void determine_number_of_individuals_sampled(Assortative_mating_parents *mating_parents);
		/// @brief mating scheme 1: fixed number of individuals sampled; 0: poisson number of individuals sampled
		float mating_scheme;
	protected:
		/**
		 * @brief Stock algorithms for determining the number of individuals sampled - fixed number
		 * 
		 * @param mating_parents @c Assortative_mating_parents object representing the parents of a species
		 */
		void determine_number_of_individuals_to_be_sampled_fixed(Assortative_mating_parents *mating_parents); 
		/**
		 * @brief Stock algorithms for determining the number of individuals sampled - poisson number
		 * 
		 * @param mating_parents @c Assortative_mating_parents object representing the parents of a species
		 */
		void determine_number_of_individuals_to_be_sampled_poisson(Assortative_mating_parents *mating_parents); 
		/**
		 * @brief Determines the mating sampling scheme
		 * 
		 * @param species_index species index
		 */
		void determine_mate_sampling_scheme(int species_index);
	};
#endif
