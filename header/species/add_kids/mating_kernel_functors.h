#ifndef REPRODUCTION_KERNEL_FUNCTORS_H
#define REPRODUCTION_KERNEL_FUNCTORS_H

#include <species/inds.h>
#include "math.h"
#include <util/thrust_functors.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>

/**
 * @brief A functor to determine whether an individual is eligible to be a parent based on its
 * 
 */
struct parental_eligibility_functor
	{
	/// @brief the sex of the parent to be eligible
	int parental_sex;
	/// @brief Construct a new parental eligibility functor object
	parental_eligibility_functor(int sex) : parental_sex(sex)
	{};

	/* 
		Elements in the tuple.

		----------------------

		0: whether the individual will reproduce
		1: the individual's sex

	*/ 

	/**
	 * @brief Call operator to determine whether an individual is eligible to be a parent based
	 * 
	 * @tparam tuple 
	 * @param t tuples containing: 0: whether the individual will reproduce, 1: the individual's sex
	 * @return __host__ 
	 */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<1>(t)==parental_sex)
			{
			thrust::get<0>(t) = 1;
			}
		}
	};

/**
 * @brief A functor to determine the fecundity
 * 
 */
struct female_fecundity_functor
	{
	/* 
		Elements in the tuple.

		----------------------
		0: whether the individual is a reproductive female
		1: the individual's fecundity phenotype
		2: the individual's final fecundity score
	*/ 

	/**
	 * @brief Call operator to determine the fecundity
	 * 
	 * @tparam tuple 
	 * @param t tuples containing: 0: whether the individual is a reproductive female, 1: the individual's fecundity phenotype, 2: the individual's final fecundity score
	 * @return __host__ 
	 */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<0>(t) > 0)
			{
			thrust::get<2>(t) = thrust::get<1>(t) ;
			}

		}
	};

#endif
