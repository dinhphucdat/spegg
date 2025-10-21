#ifndef PARENTS_CLASS_H
#define PARENTS_CLASS_H

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <curand.h>
#include <species/inds_stochastic.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <util/thrust_functors.h>

/**
 * @brief The class representing information about parents
 */
class Parents
	{
	friend class SamplingInput_Mating;

	public:
		/**
		 * @brief Constructs a Parents object from an @link inds_stochastic @endlink species instance.
		 *
		 * This constructor initializes the Parents class by copying relevant data from an
		 * `inds_stochastic` object. It sets up phenotype data, population size, deme structure,
		 * and prepares internal vectors for tracking reproductive behavior and probabilities.
		 *
		 * @param species Pointer to an `inds_stochastic` object containing the source data.
		 */
		Parents(inds_stochastic *species);
		/**
		 * @brief Sets up the reproductive eligibility and probabilities for the population.
	     *
	     * This function orchestrates the setup by calling:
	     * - `determine_parental_reproductive_potential()`
	     * - `finalize_parental_reproductive_probabilities()`
	     */
		virtual void setup_parents();
		/// @brief the vector storing the indices of the individuals that are chosen to be female parents
		thrust::device_vector<int> female_parents;
		/// @brief the vector storing the indices of the individuals that are chosen to be male parents
		thrust::device_vector<int> male_parents;
		/// @brief potential number of kids
		int Potential_Number_of_Kids;
		// TODO: maybe wipe it out? because potential should be represented by float, but this vector is not being used
		thrust::device_vector<int> reproductive_potential_per_deme;
		/// @brief reproductive males per deme
		thrust::device_vector<int> reproductive_males_per_deme;
		/// @brief reproductive males per deme
		thrust::device_vector<int> reproductive_females_per_deme;
		/// @brief number of kids per mom
		thrust::device_vector<int> kids_per_mom;

		/* stuff related to selection as parents */
		
		/// @brief probability of each individual to become a female parent
		thrust::device_vector<float> probability_individual_becomes_female_parent;
		/// @brief probability of each individual to become a male parent
		thrust::device_vector<float> probability_individual_becomes_male_parent;

	protected:
		/// @brief random generator
		curandGenerator_t gen;
		/// @brief object storing deme parameters
		DemeSettings *demeParameters;
		/// @brief population size
		int size;
		/// @brief number of demes
		int Num_Demes;
		/// @brief total number of parents
		int Total_Number_of_Parents;
		/// @brief list of phenotypes
		thrust::device_vector<float> *phenotype;
		/// @brief deme of each individual
		thrust::device_vector<int> deme;
		/// @brief individual's sex
		thrust::device_vector<int> sex;

		/* stuff related to eligibility */

		/// @brief if the individuals are females and will reproduce
		thrust::device_vector<int> will_reproduceF;
		/// @brief if the individuals are males and will reproduce
		thrust::device_vector<int> will_reproduceM;

		/**
	     * @brief Finalizes the probability that each individual becomes a parent.
	     *
	     * Calls virtual functions that can be overridden in derived classes:
	     * - `determine_probability_individual_becomes_female_parent()`
	     * - `determine_probability_individual_becomes_male_parent()`
	     */
		virtual void finalize_parental_reproductive_probabilities();

		/**
		 * @brief Determines which individuals are eligible to reproduce and their potential.
		 *
		 * This includes:
		 * - Determining male and female eligibility.
		 * - Reducing male counts per deme.
		 * - Calculating female fecundity.
		 * - Summing potential number of offspring.
		 */
		void determine_parental_reproductive_potential();

		/**
	     * @brief Determines which females are eligible to reproduce.
	     *
	     * Uses a `parental_eligibility_functor` with sex = 0 (female).
	     */
		virtual void determine_female_parent_eligibility();
		/**
	     * @brief Determines which females are eligible to reproduce.
	     *
	     * Uses a `parental_eligibility_functor` with sex = 1 (male).
	     */
		virtual void determine_male_parent_eligibility();
		
		/**
	     * @brief Virtual function to determine female parental probabilities.
	     *
	     * Intended to be overridden by species-specific implementations.
	     */
		virtual void determine_probability_individual_becomes_female_parent();
		/**
	     * @brief Virtual function to determine male parental probabilities.
	     *
	     * Intended to be overridden by species-specific implementations.
	     */
		virtual void determine_probability_individual_becomes_male_parent();

		/**
	     * @brief Calculates the fecundity of eligible females based on phenotype.
	     *
	     * Uses a phenotype index defined in `demeParameters` to compute fecundity
	     * and stores the result in `kids_per_mom`.
	     */
		virtual void female_fecundity();
	};

#endif
