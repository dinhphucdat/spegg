#ifndef MIGRATION_KERNEL_FUNCTORS_H
#define MIGRATION_KERNEL_FUNCTORS_H

#include <species/inds.h>
#include "math.h"
#include <util/thrust_functors.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>

/**
 * A thrust functor that draws a Bernoulli random variable which determines whether an individual migrates during each time step. The migration probability is presumably calculated separately, as is the uniform random number. 
 * 
 * 
 * <ol>Elements in the tuple.
 *
 *		<li> 0 - whether the individual migrates, an output vector</li>
 *		<li> 1 - their migration probability, a predefined vector</li>
 *		<li> 2 - the associated random number, drawn by Bernoulli function</li>
 * </ol>
 */
struct determine_if_migratory
	{
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<2>(t) < thrust::get<1>(t))
			thrust::get<0>(t) = 1;
		}
	};

/**
 * A thrust functor that simulates migration by reassigning the deme of the migrating individual.
 * 
 * <ol> Elements in the tuple.
 *		<li> 0: individual's deme, this will be the one that is modified during this functor's execution</li>
 *		<li> 1: whether the individual migrates</li>
 *		<li> 2: individual's destination</li>
 * </ol>
 */
struct move_individuals_functor
	{
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		thrust::get<0>(t) = thrust::get<1>(t)*thrust::get<2>(t) + (1 - thrust::get<1>(t))*thrust::get<0>(t);
		}
	};
#endif
