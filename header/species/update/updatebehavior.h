#ifndef UPDATE_BEHAVIOR_H
#define UPDATE_BEHAVIOR_H

#include <species/inds_stochastic.h>
#include <species/update/survivorship_kernel_functors.h>
#include <environ/environment.h>
#include <thrust/sequence.h>

//virtual interface; 

/**
 * @brief An interface (pure virtual class) to update the species each generation
 * 
 */
class UpdateBehavior
	{
	public:
		/* the factory */

		/**
		 * @brief Create a updateBehavior object. This serves as a factory method to create derived classes of @link UpdateBehavior @endlink
		 * 
		 * @param species @link inds_stochastic @endlink double pointer representing all species
		 * @param habitat @link environment @endlink pointer representing the habitat
		 * @param species_ID ID of the species to be updated
		 * @return UpdateBehavior* 
		 */
		static UpdateBehavior *create_updateBehavior(inds_stochastic **species, environment *habitat, int species_ID);

		/* the actual updating */

		/**
		 * @brief How to update the state of the population. <strong>User-defined function</strong>.
		 * 
		 */
		virtual void update()=0;
		/**
		 * @brief Destroy the Update Behavior object
		 * 
		 */
		virtual ~UpdateBehavior() {};

		/* Generic functionality for simulating mortality by changing the vital state variable according to a bernoulli RV */

		/**
		 * @brief Determines the mortality of individuals in a species based on their survivorship probabilities. 
		 * 
		 * This simulates mortality by changing the vital state variable according to a Bernoulli random variable.
		 * 
		 * @param species an @link inds_stochastic @endlink pointer representing the species
		 */
		void determine_mortality(inds_stochastic *species);
	};

#endif
