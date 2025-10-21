#ifndef REMOVE_DUPLICATE_PAIRS_H
#define REMOVE_DUPLICATE_PAIRS_H

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

/*************************
//
//
// util function for identifying and removing duplicates in pairs made up of two vectors; e.g.
// vectorA = [0,0,0, 1,1,1, 2,2,2, 0,0]
// vectorB = [3,4,5, 6,6,7, 8,9,9, 6,3]
// 
// Remove pairs so that:
//
// vectorA = [0,0,0,0, 1,1, 2,2]
// vectorB = [3,4,5,6, 6,7, 8,9]
//
// Note that vectorA and vectorB must be of the same size.
//
*************************/

/** @{ @name Duplication Removals */

/**
 * @brief Removes duplicate pairs made up of two vectors
 * 
 * Suppose you have two vectors, vectorA and vectorB, representing pairs of values. 
 * This function identifies and removes duplicate pairs, ensuring that each unique pair 
 * is represented only once in the output vectors.
 * 
 * E.g.,
 * @code
 * vectorA = [0,0,0, 1,1,1, 2,2,2, 0,0]
 * vectorB = [3,4,5, 6,6,7, 8,9,9, 6,3]
 * @endcode
 * After calling this function, the output will be:
 * @code
 * vectorA = [0,0,0,0, 1,1, 2,2]
 * vectorB = [3,4,5,6, 6,7, 8,9]
 * @endcode
 * 
 * @note The input vectors vectorA and vectorB must be of the same size.
 * 
 * @param vectorA The first vector of the pairs
 * @param vectorB The second vector of the pairs
 */
void remove_duplicate_pairs(thrust::device_vector<int> &vectorA, 
			    thrust::device_vector<int> &vectorB );

/** @} */

#endif
