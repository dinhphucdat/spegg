#ifndef THRUST_PROB_TABLE_DEMES_H
#define THRUST_PROB_TABLE_DEMES_H

#include <math/thrust_prob_table.h>

/**
 * @brief A derived version of the @link ThrustProbTable @endlink capable of adjusting the look-up probabilities to fall within ranges specified by the deme.
 * 
 * For instance, if you have individuals with demes @c [0,0,1,1] and look up probabilities @code{.unparsed} [0.5,0.5, 0.5, 0.5] @endcode, a ThrustProbTable_demes object transforms this into @c [0.25,0.5,0.75,1.0] , and the random draws in deme 0 will be from @c U(0,0.5) and for deme 1 from @c U(0.5,1) so that deme-specific random draws can be made in parallel.
 * 
 * @see ThrustProbTable
 */
class ThrustProbTable_demes : public ThrustProbTable
	{
	public:
		
		/**
		 * @brief Adjusts random numbers to fall within the correct intervals for each deme.
		 *
		 * This function uses Thrust library functions to map uniformly distributed random numbers to the correct intervals
		 * based on the cumulative probability table and deme offsets.
		 *
		 * @param uniform_begin Iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end Iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param deme_offsets_begin Iterator pointing to the beginning of the range of deme offsets.
		 * @param deme_offsets_end Iterator pointing to the end of the range of deme offsets.
		 * @param inds_deme_begin Iterator pointing to the beginning of the range of indices for demes.
		 * @param inds_deme_end Iterator pointing to the end of the range of indices for demes.
		 *
		 * Example:
		 * @code
		 * thrust::device_vector<float> uniform = {0.15, 0.5, 0.85};
		 * thrust::device_vector<int> deme_offsets = {0, 2, 4};
		 * thrust::device_vector<int> inds_deme(3);
		 * ThrustProbTable_demes table_demes;
		 * table_demes.setup(probabilities.begin(), probabilities.end());
		 * table_demes.adjust_randoms(uniform.begin(), uniform.end(), deme_offsets.begin(), deme_offsets.end(), inds_deme.begin(), inds_deme.end());
		 * // inds_deme will contain the adjusted indices for the random numbers
		 * @endcode
		 * 
		 * @note Please use @link ThrustProbTable::setup @endlink before calling this function
		 * 
		 * @see ThrustProbTable::setup
		 */
		void adjust_randoms(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end,
		thrust::device_vector<int>::iterator deme_offsets_begin, thrust::device_vector<int>::iterator deme_offsets_end,
		thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
		/**
		 * @brief Adjusts random numbers based on fixed offsets so that they can be used for drawing individuals for mating
		 * 
		 * @param uniform_begin iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param inds_deme_begin iterator pointing to the beginning of the range of deme sizes in the inds object.
		 * @param inds_deme_end iterator pointing to the end of the range of deme sizes in the inds object.
		 */
		void adjust_randoms_fixed_offsets(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
	};

/**
 * @brief This class behaves similarly to the @link ThrustProbTable_demes @endlink class, but with arguments of type double
 * @see ThrustProbTable_demes
 */
class ThrustProbTable_demes_Double : public ThrustProbTableDouble
	{
	public:
		/**
		 * @brief Adjusts random numbers to fall within the correct intervals for each deme.
		 * 
		 * @param uniform_begin Iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end Iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param deme_offsets_begin Iterator pointing to the beginning of the range of deme offsets.
		 * @param deme_offsets_end Iterator pointing to the end of the range of deme offsets.
		 * @param inds_deme_begin Iterator pointing to the beginning of the range of indices for demes.
		 * @param inds_deme_end Iterator pointing to the end of the range of indices for demes.
		 * 
		 * @see ThrustProbTable_demes::adjust_randoms
		 */
		void adjust_randoms(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end,
		thrust::device_vector<int>::iterator deme_offsets_begin, thrust::device_vector<int>::iterator deme_offsets_end,
		thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
		/**
		 * @brief Adjusts random numbers based on fixed offsets so that they can be used for drawing individuals for mating
		 * 
		 * @param uniform_begin iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param inds_deme_begin iterator pointing to the beginning of the range of deme sizes in the inds object.
		 * @param inds_deme_end iterator pointing to the end of the range of deme sizes in the inds object.
		 */
		void adjust_randoms_fixed_offsets(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
	};

#endif
