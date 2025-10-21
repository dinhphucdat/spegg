#ifndef SAMPLE_WITH_REPLACEMENT_H
#define SAMPLE_WITH_REPLACEMENT_H

#include <iostream>
#include <stdio.h>
#include <util/sampling_event.h>

/**
 * @brief class that inherits from SamplingEvent and represents a sampling operation with replacement. 
 * 
 * It is constructed using a SamplingInput pointer and a curandGenerator_t generator, 
 * and it overrides the sample method to perform its specific sampling behavior.
 * 
 * @see SamplingEvent
 * 
 */
class Sample_With_Replacement : public SamplingEvent
	{
	public:
		/**
		 * @brief Construct a new Sample With Replacement object
		 * 
		 * @param sampling_input pointer to a SamplingInput object
		 * @param gen random number generator
		 */
		Sample_With_Replacement(SamplingInput *sampling_input, curandGenerator_t gen )
			{
			this->sampling_input = sampling_input;
			this->gen = gen;
			}
	/**
	 * @brief Perform the sampling event. This method implements the specific sampling logic for sampling with replacement.
	 * 
	 */
	void sample() override;
	};

#endif
