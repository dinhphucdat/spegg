#ifndef STOCHASTIC_INDS_H
#define STOCHASTIC_INDS_H

#include <species/inds.h>
#include <math/random_variables_functions.h>

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
		 * @param num_demes the number of demes into which individuals of that species are split
		 * @param species_ID_val the ID with which that specific species is tagged
		 * @sa inds::inds for additional notes.
		 */
		inds_stochastic(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
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
