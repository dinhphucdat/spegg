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
		int seed;
		int MORTALITY_PHENOTYPE_INDEX;
		int FECUNDITY_PHENOTYPE_INDEX;
	};

#endif
