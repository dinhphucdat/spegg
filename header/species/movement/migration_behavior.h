#ifndef MIGRATION_BEHAVIOR_H
#define MIGRATION_BEHAVIOR_H

#include <species/inds_stochastic_migratory.h>
#include <math/thrust_prob_table_demes.h>
#include <species/movement/migration_kernel_functors.h>

/**
 * A user-defined class to define the migration behaviors of the species. The template is provided here.
 */
class MigrationBehavior
	{
	public:
		/**
		 * The factory function, which is static but its role is to create a new @c MigrationBehavior object 
		 * of your choice.
		 * 
		 * <p>Usually, users may have to create multiple child classes extending from this parent class and then 
		 * set up multiple conditions to which a corresponding object from a subset of @c MigrationBehavior class will be instantiated</p>
		 * 
		 * @code{.unparsed}
		 * MigrationBehavior* MigrationBehavior::create_migrationBehavior(inds_stochastic_migratory *species)
		 * {
		 * 		if (some condition)
		 * 			return new Subset1MigrationBehavior();
		 * 		else if (another condition)
		 * 			return new Subset2MigrationBehavior();
		 * 		...
		 * }
		 * @endcode
		 */
		static MigrationBehavior *create_migrationBehavior(inds_stochastic_migratory *species);
		
		/**
		 * Operates the migrating process through each time step of the simulation
		 * 
		 * <p>This is an <strong>abstract</strong> user-defined class.</p>
		 */
		virtual void migrate()=0;

		/**
		 * Destructor
		 */
		virtual ~MigrationBehavior() {};

		/**
		 * Variables storing the size of the population and the number of demes, aka. clusters of individuals
		 */
		int size, Number_of_Demes;
		/**
		 * Random generator, please prime this up before using it to guarantee its randomness!
		 */
		curandGenerator_t random_gen;

		/**
		 * Device vector of int type demonstrating how far each individual migrates from their initial location
		 * 
		 */
		thrust::device_vector<int> migration_offsets;
		/**
		 * Device vector of int type illustrating the destination to which individuals migrate.
		 * 
		 */
		thrust::device_vector<int> migrant_destinations;
		/**
		 * Device vector of float telling if individuals would migrate
		 * 
		 * <p>It could potentially be possibilities value from 0.0 to 1.0</p>
		 * 
		 */
		thrust::device_vector<float> will_migrate;

		/**
		 * Defined and used after the three vectors @link MigrationBehavior::migration_offsets @endlink , @link MigrationBehavior::migrant_destinations @endlink , and @link MigrationBehavior::will_migrate @endlink are filled up with values.
		 * 
		 * @sa determine_destination(inds_stochastic_migratory *species) , @sa determine_if_individuals_migrate(inds_stochastic_migratory *species)
		 */
		void move_individuals(inds_stochastic_migratory *species);
		/**
		 * User-defined function to fill up @link MigrationBehavior::migration_offsets @endlink , @link MigrationBehavior::migrant_destinations @endlink.
		 */
		void determine_destination(inds_stochastic_migratory *species);
		/**
		 * User-defined function to fill up @link MigrationBehavior::will_migrate @endlink.
		 */
		void determine_if_individuals_migrate(inds_stochastic_migratory *species);
	};

#endif
