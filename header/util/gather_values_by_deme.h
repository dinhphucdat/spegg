#ifndef VALUES_BY_DEME
#define VALUES_BY_DEME

#include <util/thrust_functors.h>

#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/scan.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/remove.h>

/* 

the function gather_values_by_deme works like this:

Suppose you have a vector A of indices

A = [A1, A2, A3, ..., An]

and a corresponding vector of the demes to which the individual Ai belongs

D = [D1, D2, D3 ,..., Dn]

If each deme has some deme-specific value stored in a vector X

B = [X1, X2, ..., Xm]

where m is the number of demes. Then gather_values_by_deme creates a resulting vector C such that:

C = [X[D1], X[D2], ..., X[Dm]]


*/

/** @{ @name @c gather operations acted on demes */

/**
 * @brief Gathers values based on demes
 * 
 * Suppose you have a vector A of indices
 * @code
 * A = [A1, A2, A3, ..., An]
 * @endcode
 * and a corresponding vector of the demes to which the individual Ai belongs
 * @code
 * D = [D1, D2, D3 ,..., Dn]
 * @endcode
 * If each deme has some deme-specific value stored in a vector X
 * @code
 * B = [X1, X2, ..., Xm]
 * @endcode
 * where m is the number of demes. Then gather_values_by_deme creates a resulting
 * vector C such that:
 * @code
 * C = [X[D1], X[D2], ..., X[Dm]]
 * @endcode
 * 
 * @param indices indices of individuals
 * @param demes demes of individuals
 * @param deme_specific_value deme-specific values
 * @param ans output vector containing the gathered values based on demes
 */
void gather_values_by_deme(thrust::device_vector<int> &indices, 
			   thrust::device_vector<int> &demes, 
			   thrust::device_vector<int> &deme_specific_value, 
			   thrust::device_vector<int> &ans);

/**
 * @brief Gathers values based on demes - floating point deme-specific values
 * 
 * @see gather_values_by_deme(thrust::device_vector<int> &indices, thrust::device_vector<int> &demes, thrust::device_vector<int> &deme_specific_value, thrust::device_vector<int> &ans))
 * 
 * @param indices indices of individuals
 * @param demes demes of individuals
 * @param deme_specific_value floating point deme-specific values
 * @param ans output vector containing the gathered values based on demes
 */
void gather_values_by_deme(thrust::device_vector<int> &indices, 
			   thrust::device_vector<int> &demes, 
			   thrust::device_vector<float> &deme_specific_value, 
			   thrust::device_vector<float> &ans);

/**
 * @brief Gathers values based on demes - floating point deme-specific values, with the deme-specific values provided as an iterator
 * 
 * @param indices indices of individuals
 * @param demes demes of individuals
 * @param deme_specific_values_begin iterator pointing to the beginning of the range of floating point deme-specific values
 * @param ans output vector containing the gathered values based on demes
 */
void gather_values_by_deme(thrust::device_vector<int> &indices,
			   thrust::device_vector<int> &demes,
                           thrust::device_vector<float>::iterator deme_specific_values_begin,
                           thrust::device_vector<float> &ans);

/** @} */

#endif
