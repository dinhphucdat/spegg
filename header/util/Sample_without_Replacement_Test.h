#ifndef SAMPLE_WITHOUT_REPLACEMENT_Test_H
#define SAMPLE_WITHOUT_REPLACEMENT_Test_H

#include <thrust/sort.h>
#include <thrust/copy.h>
#include <util/sampling_event.h>

/**
 * @brief class that inherits from SamplingEvent and represents a sampling operation without replacement, using a test method.
 * 
 * @see SamplingEvent
 * 
 */
class Sample_without_Replacement_Test : public SamplingEvent
	{
	public:
		/**
		 * @brief Construct a new Sample_without_Replacement_Test object
		 * 
		 * @param sampling_input pointer to a SamplingInput object
		 * @param gen random number generator
		 */
		Sample_without_Replacement_Test(SamplingInput *sampling_input, curandGenerator_t gen )
			{
			this->sampling_input = sampling_input;
			this->gen = gen;

			number_of_individuals_subject_to_sampling = sampling_input->list_of_individuals_potentially_subject_to_sampling.size();
			}
		/**
		 * @brief Perform the sampling event. This method implements the 
		 * specific sampling logic for sampling without replacement using a test method.
		 * 
		 */
		void sample() override;

	protected:
		/// @brief number of individuals to be sampled
		int number_of_individuals_subject_to_sampling;
		/// @brief unique uniform random variables for sampling individuals being sampled
		thrust::device_vector<double> unique_uniform_rvs;
		/**
		 * @brief Draws unique double random numbers for sampling individuals being sampled.
		 * 
		 */
		void draw_unique_double_randoms();
		/**
		 * @brief Shuffles the sampled individuals at random.
		 * 
		 */
		void shuffle_sampled_individuals_at_random();	
	};

#endif
