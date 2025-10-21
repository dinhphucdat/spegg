#ifndef SAMPLING_EVENT_H
#define SAMPLING_EVENT_H

#include <iostream>
#include <stdio.h>
#include <util/sampling_input.h>

/**
 * @brief An interface (pure virtual class) for sampling events
 * 
 */
class SamplingEvent
	{
	public:
	/* use a factory method */

		/**
		 * @brief Create a SamplingEvent object. This serves as a factory method to create derived classes of @link SamplingEvent @endlink
		 * 
		 * @param sampling_input @link SamplingInput @endlink object
		 * @param gen random number generator
		 * @return SamplingEvent* 
		 */
		static SamplingEvent *create_SamplingEvent(SamplingInput *sampling_input, curandGenerator_t gen);
		/// @brief Sampling individuals indices
		thrust::device_vector<int> sampling_individuals_indices;
		/// @brief Sampled individuals indices
		thrust::device_vector<int> sampled_individuals_indices;
		/**
		 * @brief Perform the sampling event. This is a pure virtual function that must be implemented by derived classes.
		 * 
		 */
		virtual void sample() = 0;
	
	protected:
		/// @brief Pointer to the SamplingInput object
		SamplingInput *sampling_input;
		/// @brief Random number generator
		curandGenerator_t gen;
	};

#endif
