#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <curand.h>
#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/distance.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/replace.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/gather.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <thrust/unique.h>
#include <thrust/set_operations.h>
#include <thrust/extrema.h>

#include <util/thrust_functors.h>
#include <math/thrust_probabilities.h>

/**
 * @{ @name Histogram's functions
 */

/**
 * @brief Calculate the histograms from the @c data and output into @c counts 
 * 
 * @param data vector of int to calculate histograms
 * @param counts vector to buffer histograms into
 * @param counts_size number of bins for histograms
 */
void calculate_histogram(thrust::device_vector<int> &data, thrust::device_vector<int> &counts, 		 int counts_size);
/**
 * @brief Calculate the histograms from the @c data and output into @c counts , deal with float data
 * 
 * @param data vector of float to calculate histograms
 * @param counts vector to buffer histograms into
 * @param counts_size number of bins for histograms
 * 
 * @see calculate_histogram
 */
void calculate_histogram(thrust::device_vector<float> &data, thrust::device_vector<int> &counts, 		 int counts_size);
/**
 * @brief Extract the counts for a subset of the histogram elements pertaining to specific values
 * 
 * Example: @c data=(0,0,0,1,1,2,2,2,2,3) with histogram @c counts=(3,2,4,1), extract only counts for element 1 and 3 so end up with @c histogram=(2,1)
 * 
 * @param data vector of int to calculate histograms
 * @param counts vector to buffer histograms into
 * @param counts_size number of bins for histograms
 * @param subset_indices indices chosen to output the histogram. Histograms of other value indices will be truncated
 */
void calculate_histogram_subset(thrust::device_vector<int> &data, thrust::device_vector<int> &counts, int counts_size, thrust::device_vector<int> &subset_indices);

/** @} */

/**
 * @brief Convert values to fall inside bin groups
 * 
 * <strongElements in the tuple:</strong>
 * <br>
 * - 0: original data value
 * - 1: minimum data value
 * - 2: maximum data value
 * - 3: number of bins
 * - 4: new data value
 * 
 */
struct adjust_histogram_data_values
	{
	/* 
		Elements in the tuple.
		----------------------
		0: original data value
		1: minimum data value
		2: maximum data value
		3: number of bins
		4: new data value
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		thrust::get<4>(t) = thrust::get<0>(t) - thrust::get<1>(t);
		thrust::get<4>(t) = thrust::get<4>(t)/(thrust::get<2>(t) - thrust::get<1>(t));
		thrust::get<4>(t) *= (float) (thrust::get<3>(t) - 1);
		}
	
	};

#endif
