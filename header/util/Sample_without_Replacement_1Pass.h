#ifndef SAMPLE_WITHOUT_REPLACEMENT_1PASS_H
#define SAMPLE_WITHOUT_REPLACEMENT_1PASS_H

#include <iostream>
#include <stdio.h>
#include <util/sampling_event.h>

/**
 * @brief class that inherits from SamplingEvent and represents a sampling operation without replacement using a one-pass algorithm.
 * 
 * @see SamplingEvent
 * 
 */
class Sample_without_Replacement_1Pass : public SamplingEvent
	{
	public:
		/**
		 * @brief Construct a new Sample_without_Replacement_1Pass object
		 * 
		 * @param sampling_input pointer to a SamplingInput object
		 * @param gen random number generator
		 */
		Sample_without_Replacement_1Pass(SamplingInput *sampling_input, curandGenerator_t gen )
			{
			this->sampling_input = sampling_input;
			this->gen = gen;
			}
	/**
	 * @brief Perform the sampling event. This method implements the specific sampling 
	 * logic for sampling without replacement using a one-pass algorithm.
	 * 
	 */
	void sample() override;
	};

#endif
