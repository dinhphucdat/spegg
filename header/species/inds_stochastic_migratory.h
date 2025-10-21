#ifndef STOCHASTIC_INDS_MIGRATORY_H
#define STOCHASTIC_INDS_MIGRATORY_H

#include <species/inds_stochastic.h>
#include <math/random_variables_functions.h>
#include <math/thrust_prob_table_demes.h>

/**
 * This class is a subset of @link inds inds @endlink class, and it also includes the functionality allowing 
 * the simulation of migration.
 */
class inds_stochastic_migratory : public inds_stochastic
	{
	public:
		/**
		 * <p>This constructor essentially initializes the <a href="">migration matrix</a></p>
		 * 
		 * <p>The object of this class expects a migration file required to be named as @c Migration_Probabilities.txt , 
		 * with each individual's migration rate in each line. <strong>See more the guide</strong> on @c Migration_Probabilities.txt .</p>
		 */
		inds_stochastic_migratory(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
		/**
		 * A user-defined function for migration behavior
		 */
		virtual void migrate() {};
		/**
		 * A user-defined function for migration behavior
		 */
		virtual void migrate(inds_stochastic_migratory **species) {};
		/**
		 * A user-defined function for migration behavior
		 */
		virtual void migrate(inds_stochastic_migratory **species, environment *habitat) {};
		

		/**
		 * The baseline migration rates are set by migration matrix object
		 */
		ThrustProbTable_demes Migration_Matrix;

		/**
		 * The vector individual_migration_rate stores the probability an individual will migrate on a given 
		 * time step; most likely this should be calculated as an update_behavior since migratory potential 
		 * is a phenotype and can be calculated using some functor like @code{.cpp} struct migration_probability_functor{} @endcode
		 * with code that calculates the probability that an individual migrates during each time step.
		 */
		thrust::device_vector<float> individual_migration_rate;
	};

#endif
