#ifndef REDUCE_BY_KEY_WITH_ZEROS_H
#define REDUCE_BY_KEY_WITH_ZEROS_H

#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/unique.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/extrema.h>

/** @{ @name Custom cumulative sum */

/**
 * @brief Reduces values by keys, ensuring that keys with zero values are represented in the output
 * 
 * @code{.cpp}
 * max_possible_keys = 4
 * thrust::device_vector<int> key(5);
 * thrust::device_vector<int> val(5);
 * key[0] = 1;
 * key[1] = 1;
 * key[2] = 1;
 * key[3] = 2;
 * key[4] = 2;
 * for (int i=0; i < 5; i++)
 *      val[i] = i;
 * thrust::device_vector<int> key_out(max_possible_keys);
 * thrust::device_vector<int> output(max_possible_keys);
 * thrust::reduce_by_key(key.begin(), key.end(), val.begin(), key_out.begin(), output.begin());	
 * 
 * returns, in the vector output = [3, 7] and key_out=[1,2]. However, it does not give us 
 * output_complete = [0, 3, 7, 0] for keys =[0,3]. This function returns output_complete.
 * @endcode
 * 
 * @param key_vector vector of keys
 * @param values vector of values
 * @param values_output output vector to store the reduced values by keys, including zeroes for keys with no associated values
 * @param number_of_values number of values in the input vectors
 * @param number_of_possible_keys number of possible keys
 */
void reduce_by_key_with_zeros(thrust::device_vector<int> &key_vector, thrust::device_vector<int> &values, thrust::device_vector<int> &values_output, int number_of_values, int number_of_possible_keys);

/**
 * @brief Reduces floating point values by keys, ensuring that keys with zero values are represented in the output
 * 
 * @see reduce_by_key_with_zeros(thrust::device_vector<int> &key_vector, thrust::device_vector<int> &values, thrust::device_vector<int> &values_output, int number_of_values, int number_of_possible_keys)
 * 
 * @param key_vector vector of keys
 * @param values vector of floating point values
 * @param values_output output vector to store the reduced floating point values by keys, including zeroes for keys with no associated values
 * @param number_of_values number of values in the input vectors
 * @param number_of_possible_keys number of possible keys
 */
void reduce_by_key_with_zeros(thrust::device_vector<int> &key_vector, thrust::device_vector<float> &values, thrust::device_vector<float> &values_output, int number_of_values, int number_of_possible_keys);

/** @} */

#endif
