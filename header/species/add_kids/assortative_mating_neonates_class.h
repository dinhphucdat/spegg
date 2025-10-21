#ifndef ASSORTATIVE_MATING_NEONATES_H
#define ASSORTATIVE_MATING_NEONATES_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <util/footimer2.h>
#include <species/inds.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <species/deme_specific_data_class.h>
#include <util/amplify.h>
#include <math/mating_thrust_prob_table_demes.h>
#include <math/random_variables_functions.h>
#include <species/add_kids/neonates_class.h>
#include <species/add_kids/assortative_mating_parents_class.h>
#include <util/one_dim_two_dim.h>

#include <thrust/adjacent_difference.h>

/**
 * @brief An extension to @link EggsNeonates @endlink where mating pairs are selectively chosen rather than randomly
 * 
 */
class Assortative_mating_neonates :  public EggsNeonates 
	{
	public:
		/**
		 * @brief Construct a new Assortative_mating_neonates object
		 * 
		 * @param species an @link inds_stochastic @endlink pointer representing the species
		 * @param kids_per_mom kids born per mother, meant to use the @c kids_per_mom vector from the @link Parent @endlink class
		 */
		Assortative_mating_neonates(
			inds_stochastic* species, 
			thrust::device_vector<int> kids_per_mom);
		
		/**
		 * @brief Simulates the a reproduction process which leads to changes in @link inds::mgenotype @endlink and @link inds::fgenotype @endlink
		 * 
		 * @param probability_pair_becomes_parents probability pair would become parents
		 * @param fathers_list list of fathers
		 * @param mothers_list list of mothers
		 * @param fgenotype maternal genotypes, possibly taken from @link inds::fgenotype @endlink
		 * @param mgenotype paternal genotypes, possibly taken from @link inds::mgenotype @endlink
		 */
		void inherit_genotypes_by_pair(
			thrust::device_vector<float> probability_pair_becomes_parents, 
			thrust::device_vector<int> fathers_in_pairs, 
			thrust::device_vector<int> mothers_in_pairs, 
			thrust::device_vector<int> pairs_demes);
		/**
		 * @brief get the mating pairs
		 * 
		 * This function performs algorithm to select the mating pairs and store into @link Assortative_mating_neonates::fathers_chosen @endlink and 
		 * @link Assortative_mating_neonates::mothers_chosen @endlink.
		 * 
		 * @param probability_pair_becomes_parents probability pair would become parents
		 * @param fathers_list list of fathers
		 * @param mothers_list list of mothers
		 */
		void get_mating_pair(
			thrust::device_vector<float> probability_pair_becomes_parents, 
			thrust::device_vector<int> fathers_in_pairs, 
			thrust::device_vector<int> mothers_in_pairs, 
			thrust::device_vector<int> pairs_demes
		);

		/**
		 * @brief Maps the id to the corresponding mother chosen for mating, get it stored in @c maternal_id and @c paternal_id
		 * 
		 * @param maternal_id vector storing the ids of mothers
		 * @param paternal_id vector storing the ids of fathers
		 * @param ids ids of individuals in the population
		 */
		void record_parents(thrust::device_vector<int> &maternal_id, 
				     thrust::device_vector<int> &paternal_id,
				     thrust::device_vector<int> &ids);

	protected:
		/// @brief list of indices of chosen mothers
		thrust::device_vector<int> mothers_chosen;
		/// @brief list of indices of chosen fathers
		thrust::device_vector<int> fathers_chosen;
		/// @brief number of subpopulations or demes
		int Num_Subpopulations;
		/// @brief current population size
		int current_pop_size;
		/**
		 * @brief Copy the chosen maternal genotype to the child's @c mgenotype
		 * 
		 * @param mother_index list of mothers' indices
		 * @param mgenotype paternal genotypes, possibly taken from @link inds::mgenotype @endlink
		 * @param fgenotype maternal genotypes, possibly taken from @link inds::fgenotype @endlink
		 */
		void get_maternally_derived_genotype_deterministic(
						     thrust::device_vector<float> *&mgenotype,
						     thrust::device_vector<float> *&fgenotype);
		/**
		 * @brief Copy the chosen maternal genotype to the child's @c fgenotype
		 * 
		 * @param mother_index list of mothers' indices
		 * @param mgenotype paternal genotypes, possibly taken from @link inds::mgenotype @endlink
		 * @param fgenotype maternal genotypes, possibly taken from @link inds::fgenotype @endlink
		 */
		void get_paternally_derived_genotype_deterministic(
						     thrust::device_vector<float> *&mgenotype,
						     thrust::device_vector<float> *&fgenotype);
		/// @brief list of demewise number of pairs
		thrust::device_vector<int> pairs_per_deme;
	};
#endif
