#ifndef THRUST_FUNCTORS_H
#define THRUST_FUNCTORS_H

#include <thrust/tuple.h>
#include <cmath>

/** @{ @name MAX_RV */
#define MAX_RV 0.999991f /**< Used because floats occasionally cause indexing problems in reassign_functor(). */
/** @} */

/**
 * @brief Unary functor that adds a constant RHS value to its argument.
 * @tparam T numeric type
 */
template <typename T>
struct unary_plus
{
	T rhs; /**< Right-hand-side additive constant */
	unary_plus(const T &r) : rhs(r) {}

	__host__ __device__
	T operator()(const T &lhs) {
		return lhs + rhs;
	}
};

/**
 * @brief Unary functor that subtracts a constant RHS value from its argument.
 * @tparam T numeric type
 */
template <typename T>
struct unary_minus
{
	T rhs; /**< Right-hand-side subtractive constant */
	unary_minus(const T &r) : rhs(r) {}

	__host__ __device__
	T operator()(const T &lhs) {
		return lhs - rhs;
	}
};

/**
 * @brief Unary functor that computes modulus with a constant RHS.
 * @tparam T integer type
 */
template <typename T>
struct unary_modulus
{
	T rhs; /**< Right-hand-side modulus constant */
	unary_modulus(const T &r) : rhs(r) {}

	__host__ __device__
	T operator()(const T &lhs) {
		return lhs % rhs;
	}
};

/// Comparison functors -----------------------------------------------------
template <typename T>
struct unary_less
{
	T rhs;
	unary_less(const T &r) : rhs(r) {}

	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs < rhs;
	}
};

template <typename T>
struct unary_less_equal
{
	T rhs;
	unary_less_equal(const T &r) : rhs(r) {}

	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs <= rhs;
	}
};

template <typename T>
struct unary_greater
{
	T rhs;
	unary_greater(const T &r) : rhs(r) {}

	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs > rhs;
	}
};

template <typename T>
struct unary_not_equal
{
	T rhs;
	unary_not_equal(const T &r) : rhs(r) {}

	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs != rhs;
	}
};

template <typename T>
struct unary_greater_equal
{
	T rhs;
	unary_greater_equal(const T &r) : rhs(r) {}

	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs >= rhs;
	}
};

/**
 * @brief Predicate: true if integer is less than zero.
 */
struct is_less_than_zero
{
	__host__ __device__
	bool operator()(int x)
	{
		return x < 0;
	}
};

/**
 * @brief Predicate: true if float is less than zero.
 */
struct is_less_than_zero_f
{
	__host__ __device__
	bool operator()(float x)
	{
		return x < 0;
	}
};

/**
 * @brief Functor that writes the pairwise minimum of two inputs into a third tuple element.
 *
 * Tuple layout expected:
 *  - 0: left-hand value
 *  - 1: right-hand value
 *  - 2: output value (pairwise min)
 */
struct pairwise_min
{
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<0>(t) <= thrust::get<1>(t))
			thrust::get<2>(t) = thrust::get<0>(t);
		else
			thrust::get<2>(t) = thrust::get<1>(t);
	}
};

/**
 * @brief Reassign a sampled floating-point uniform [0,1) value into an integer index range per-population.
 *
 * Expected tuple layout:
 *  - 0: sampled uniform float
 *  - 1: source population index (sampler_from)
 *  - 2: output integer index (written by this functor)
 *
 * The functor uses an externally provided array `sampled_offset` that stores cumulative offsets
 * for each source population. This mirrors the original implementation's behavior on GPU.
 */
struct reassign_functor
{
	int *sampled_offset; /**< pointer to array of cumulative offsets */
	reassign_functor(int *d_offset) : sampled_offset(d_offset) {}

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int sampler_from = thrust::get<1>(t);

		if (sampler_from != 0)
		{
			thrust::get<2>(t) = (int) (thrust::get<0>(t) * (sampled_offset[sampler_from] - sampled_offset[sampler_from-1]) + (sampled_offset[sampler_from-1]));

			// Correction for floating point edge-case that would produce an out-of-range index.
			if (thrust::get<2>(t) == (sampled_offset[sampler_from] - sampled_offset[sampler_from-1]) + (sampled_offset[sampler_from-1]))
			{
				thrust::get<2>(t) = thrust::get<2>(t) - 1;
			}
		}
		else
		{
			thrust::get<2>(t) = (int) (thrust::get<0>(t) * (sampled_offset[sampler_from]));
			if (thrust::get<2>(t) == sampled_offset[sampler_from])
			{
				thrust::get<2>(t) = sampled_offset[sampler_from] - 1;
			}
		}
	}
};

/**
 * @brief Simplified reassign functor that maps a float to an integer in [lower, upper).
 *
 * Expected tuple layout:
 *  - 0: floating-point number in [0,1)
 *  - 1: lower bound (int)
 *  - 2: upper bound (int)
 *  - 3: output integer
 */
struct reassign_functor_simplified
{
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<0>(t) > MAX_RV)
		{
			thrust::get<0>(t) = MAX_RV;
		}

		thrust::get<3>(t) = (int) (thrust::get<0>(t) * (thrust::get<2>(t) - thrust::get<1>(t)));
	}
};

/**
 * @brief Functor that converts a float in tuple[0] to nearest int in tuple[1].
 *
 * Tuple: [float_input, int_output]
 */
struct float_to_int
{
	template <typename tuple>
	__host__ __device__
	void operator() (tuple t) {
		thrust::get<1>(t) = (int) std::floor(thrust::get<0>(t) + 0.5f);
	}
};

/**
 * @brief Converts an int (tuple[0]) to float (tuple[1]).
 */
struct int_to_float
{
	template <typename tuple>
	__host__ __device__
	void operator() (tuple t) {
		thrust::get<1>(t) = (float) thrust::get<0>(t);
	}
};

/**
 * @brief Converts an int (tuple[0]) to double (tuple[1]).
 */
struct int_to_double
{
	template <typename tuple>
	__host__ __device__
	void operator() (tuple t) {
		thrust::get<1>(t) = (double) thrust::get<0>(t);
	}
};

/**
 * @brief Computes absolute value of a float stored in tuple[0] in-place.
 */
struct fabs_functor
{
	template <typename tuple>
	__host__ __device__
	void operator() (tuple t) {
		if (thrust::get<0>(t) < 0)
		{
			thrust::get<0>(t) = -1 * thrust::get<0>(t);
		}
	}
};

#endif
