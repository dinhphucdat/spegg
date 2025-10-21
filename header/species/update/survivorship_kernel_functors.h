#ifndef SURVIVORSHIP_KERNEL_FUNCTORS_H
#define SURVIVORSHIP_KERNEL_FUNCTORS_H

#include <species/inds.h>
#include "math.h"
#include <util/thrust_functors.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>

/**
 * @brief A functor to simulate mortality
 * 
 */
struct simulate_mortality
	{
	/* 
	* Functor related to survivorship: to be invoked in the update_mySpecies class's determine_mortality() function, specified in \loc cuda/species/update/survivorship_kernel.cu.
	*/
	/// @brief pointer to uniformly distributed random variables
	float *uniform_rv;
	/**
	 * @brief Construct a new simulate mortality object
	 * 
	 * @param uniform_random pointer to uniformly distributed random variables
	 */
	simulate_mortality(float* uniform_random) : uniform_rv(uniform_random)
	{};

	/* 
		Elements in the tuple.


		----------------------
		0: individual index
		1: individual's vital state
		2: probability of survivorship

	*/ 

	/**
	 * @brief Call operator to simulate mortality
	 * 
	 * @tparam tuple 
	 * @param t tuples containing: 0: individual index, 1: individual's vital state, 2: probability of survivorship
	 * @return __host__ 
	 */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int vital_state = thrust::get<1>(t);

		if (vital_state==1) /* If the individual is alive */
			{
			float unif = uniform_rv[thrust::get<0>(t)];
			float surv_val = thrust::get<2>(t);
			
			// Actually kill them
			if (unif < surv_val)
				vital_state = 1;
			else
				vital_state = 0;
			thrust::get<1>(t) = vital_state;
			}
		}	
};

#endif
