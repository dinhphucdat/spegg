#ifndef THRUST_PROB_TABLE_H
#define THRUST_PROB_TABLE_H

#include <thrust/device_vector.h>

/**
 * @brief A probability table for selecting random individuals based on arbitrary discrete distributions
 */
class ThrustProbTable
	{
	public:
		/**
		 * @brief Initializes the cumulative probability table based on the input probability values.
		 * 
		 * This function computes the cumulative probabilities from the given probabilities using Thrust library functions.
		 * 
		 * Example:
		 * @code
		 * thrust::device_vector<float> probabilities = {0.1, 0.2, 0.3, 0.4};
		 * ThrustProbTable table;
		 * table.setup(probabilities.begin(), probabilities.end());
		 * // cumulative_prob will be {0.1, 0.3, 0.6, 1.0}
		 * @endcode
		 * 
		 * 
		 * @param prob_begin Iterator pointing to the beginning of the range of input probabilities.
 		 * @param prob_end Iterator pointing to the end of the range of input probabilities.
		 */
		void setup(thrust::device_vector<float>::iterator prob_begin, thrust::device_vector<float>::iterator prob_end);

		/**
		 * @brief Maps uniformly distributed random numbers to cumulative probabilities.
		 *
		 * This function uses `thrust::lower_bound` to find the positions in the cumulative probability table
		 * where the values from the uniform distribution should be inserted to maintain sorted order.
		 *
		 * @param uniform_begin Iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end Iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param result Iterator pointing to the beginning of the range where the resulting indices will be stored.
		 *
		 * Example:
		 * 
		 * @code{.cpp}
		 * thrust::device_vector<float> uniform = {0.15, 0.5, 0.85};
		 * thrust::device_vector<int> result(3);
		 * ThrustProbTable table;
		 * table.setup(probabilities.begin(), probabilities.end());
		 * table.draw(uniform.begin(), uniform.end(), result.begin());
		 * // result will contain the indices corresponding to the uniform values in the cumulative_prob vector
		 * @endcode
		 * 
		 */
		void draw(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator result);
		/**
		 * @brief Gets a pointer to the cumulative probabilities.
		 * 
		 * This function returns a device pointer to the first element of the cumulative probabilities vector.
		 * 
		 * @return thrust::device_ptr<float> - A device pointer to the first element of the cumulative probabilities vector.
		 */
		thrust::device_ptr<float> get_cumulative_probabilities_ptr();

	protected:
		/// @brief A look-up table to be used for simulating from arbitrary discrete distributions.
		thrust::device_vector<float> cumulative_prob;
	};

/**
 * @brief This class behaves similarly to @link ThrustProbTable @endlink class, but receives argument of type doubles
 * @see ThrustProbTable
 */
class ThrustProbTableDouble
	{
	public:
		/**
		 * @brief Initializes the cumulative probability table based on the input probability values.
		 * 
		 * @see ThrustProbTable::setup
		 * 
		 * @param prob_begin Iterator pointing to the beginning of the range of input probabilities.
		 * @param prob_end Iterator pointing to the end of the range of input probabilities.
		 */
		void setup(thrust::device_vector<double>::iterator prob_begin, thrust::device_vector<double>::iterator prob_end);
		/**
		 * @brief Maps uniformly distributed random numbers to cumulative probabilities.
		 * 
		 * @see ThrustProbTable::draw
		 * 
		 * @param uniform_begin Iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end Iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param result Iterator pointing to the beginning of the range where the resulting indices will be stored.
		 */
		void draw(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end, thrust::device_vector<int>::iterator result);
		/**
		 * @brief Gets a pointer to the cumulative probabilities.
		 * 
		 * @see ThrustProbTable::get_cumulative_probabilities_ptr
		 * 
		 * @return thrust::device_ptr<double> - A device pointer to the first element of the cumulative probabilities vector.
		 */
		thrust::device_ptr<double> get_cumulative_probabilities_ptr();

	protected:
		/// @brief A look-up table to be used for simulating from arbitrary discrete distributions. @see ThrustProbTable::cumulative_prob
		thrust::device_vector<double> cumulative_prob;
	};
#endif
