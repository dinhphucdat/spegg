#ifndef THRUST_PROBABILITIES_H
#define THRUST_PROBABILITIES_H

/**
 * @{ \name Essential constants for random generators
 */
/// @brief Pi value
#define PI 3.14159265358979f
/// @brief The maximum value for Poisson random variables
#define POISSON_MAX_RV 0.99999991f
/// @brief The maximum number of iterations for rejection sampling
#define MAXIMUM_ITERATIONS_REJECTION_SAMPLING 1000000
/**
 * @}
 */

#include "math.h"
#include <thrust/functional.h>
#include <stdio.h>

/**
 * @{ \name Random functors
 * 
 */
/**
 * @brief Converts a uniform random variable to a discrete value between 0 and maxSize-1.
 */
struct discrete_uniform
	{
	// Takes the float value given by curandGenerateUniform() and converts it into a discrete value between 0, n-1
	const int maxSize;
	/**
	 * @brief Construct a new discrete uniform object
	 * 
	 * @param _maxSize The maximum size of the discrete range.
	 */
	discrete_uniform(int _maxSize) : maxSize(_maxSize) {};
	

	/**
	 * @brief Operator to convert a uniform random variable to a discrete value.
	 *
	 * @param rv The uniform random variable.
	 * @return __host__ The discrete value.
	 */
	__host__ __device__
	int operator()(const float& rv) const {
		int ans = 0;
		float ansf = ((float) maxSize)*rv;
		ans = (int) ansf;
		return ans;
		}
	};

/**
 * @brief Generates Bernoulli-distributed random variates.
 */
struct bernoulli_rv
	{
	const float probability;
	/**
	 * @brief Construct a new bernoulli rv object
	 * 
	 * @param _probability The probability parameter for the Bernoulli distribution.
	 */
	bernoulli_rv(float _probability) : probability(_probability) {};
	
	/**
	 * @brief Operator to generate a Bernoulli random variate.
	 *
	 * @param rv The uniform random variable.
	 * @return __host__ The Bernoulli random variate (0 or 1).
	 */
	__host__ __device__
	int operator()(const float& rv) const {
		int ans = 0;
		if (rv <= probability)
			ans = 1;
		if (rv > probability)
			ans = 0;
		return ans;
		}
	};

/**
 * @brief Generates Poisson-distributed random variates.
 */
struct poisson_rv
	{
	const float lambda;
	/**
	 * @brief Construct a new poisson rv object
	 * 
	 * @param _lambda The lambda parameter for the Poisson distribution.
	 */
	poisson_rv(float _lambda) : lambda(_lambda) {};
	
	/**
	 * @brief Operator to generate a Poisson random variate.
	 *
	 * @param rv The uniform random variable.
	 * @return __host__ The Poisson random variate.
	 */
	__host__ __device__
	int operator()(const float& rv) const {
		float use_rv = rv;
		if (use_rv > POISSON_MAX_RV)
			{
			use_rv = POISSON_MAX_RV;
			}
		int ans = 0;
		float p = exp(-lambda);
		float F = p;
		int maxits = 0;
		while ((use_rv > F) && (maxits < MAXIMUM_ITERATIONS_REJECTION_SAMPLING))
			{
			p = lambda*p/(ans + 1);
			F = F + p;
			ans = ans + 1;
			maxits++;
			}
		return ans;
		}
	};

/**
 * @brief Generates discrete normal-distributed random variates.
 * 
 * <strong>Elements in the tuple.</strong>
 * <br>
 * - 0: uniform rv 1
 * - 1: uniform rv 2
 * - 2: ans
 */
struct discrete_normal_rv
	{
	const float mean;
	const float sd;
	/**
	 * @brief Construct a new discrete normal rv object
	 * 
	 * @param _mean The mean of the normal distribution
	 * @param _sd The standard deviation of the normal distribution
	 */
	discrete_normal_rv(float _mean, float _sd) : mean(_mean), sd(_sd) {};
	/**
	 * @brief Operator to generate a discrete normal random variate.
	 * 
	 * @tparam tuple tuple type
	 * @param t Tuple containing uniform random variables and the output variate.
	 * @return __host__ 
	 */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<2>(t) = (int) sd * r * cosf(phi) + mean;
		}
	};
/**
 * @brief Generates normal-distributed random variates.
 * 
 * <strong>Elements in the tuple.</strong>
 * 
 * <br>
 * 
 * - 0: uniform rv 1
 * - 1: uniform rv 2
 * - 2: ans
 */
struct normal_rv
	{
	const float mean;
	const float sd;
	/**
	 * @brief Construct a new normal rv object
	 * 
	 * @param _mean The mean of the normal distribution.
	 * @param _sd The standard deviation of the normal distribution.
	 */
	normal_rv(float _mean, float _sd) : mean(_mean), sd(_sd) {};

	/**
	 * @brief Operator to generate a normal random variate.
	 * 
	 * @tparam tuple argument tuples
	 * @param t Tuple containing uniform random variables and the output variate.
	 * @return __host__ 
	 */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<2>(t) = sd * r * cosf(phi) + mean;
		}
	};
/**
 * @brief Generates normal-distributed random variates with different parameters.
 * 
 * 
 * This functor uses the Box-Muller transform to convert uniform random variables into Gaussian-distributed random variables.
 * It allows for different means and standard deviations for each generated variate.
 * 
 * Elements in the tuple:
 * 
 * - 0: Uniform random variable 1
 * - 1: Uniform random variable 2
 * - 2: Mean of the Gaussian distribution
 * - 3: Standard deviation of the Gaussian distribution
 * - 4: Output Gaussian random variate
 * 
 */
struct normal_rv_different_parameters
	{	
	/**
	 * @brief Generates normal-distributed random variates with different parameters.
	 * 
	 * @tparam tuple tuple type
	 * @param t Tuple containing uniform random variables and the output variate.
	 * @return __host__ 
	 */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<4>(t) = thrust::get<3>(t) * r * cosf(phi) + thrust::get<2>(t);
		}
	};

/**
 * @brief Generates discrete normal-distributed random variates with different parameters.
 * @see normal_rv_different_parameters
 */
struct discrete_normal_rv_different_parameters
	{
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<4>(t) = (int) (thrust::get<3>(t) * r * cosf(phi) + thrust::get<2>(t));
		}
	};

/**
 * @brief Generates Poisson-distributed random variates with different parameters.
 *
 * This functor uses rejection sampling to convert uniform random variables into Poisson-distributed random variables.
 * It allows for different lambda parameters for each generated variate.
 *
 * Elements in the tuple:
 * - 0: Uniform random variable
 * - 1: Lambda parameter for the Poisson distribution
 * - 2: Output Poisson random variate
 * 
 * @see poisson_rv
 */
struct poisson_rv_different_parameters
	{
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ans = 0;
		float lambda = thrust::get<1>(t);
		float rv = thrust::get<0>(t);
		if (rv > POISSON_MAX_RV)
			{
			rv = POISSON_MAX_RV;
			}

		float p = exp(-lambda);
		float F = p;
		int maxits = 0;
		while ((rv > F) && (maxits < MAXIMUM_ITERATIONS_REJECTION_SAMPLING))
			{
			maxits++;
			p = lambda*p/(ans + 1);
			F = F + p;
			ans = ans + 1;
			}
		thrust::get<2>(t) = ans;
		}
	};

/**
 * @brief Generates Bernoulli-distributed random variates with different probability parameters.
 *
 * Elements in the tuple:
 * - 0: Uniform random variable
 * - 1: Probability of success
 * - 2: Output Bernoulli random variate
 * 
 * @see poisson_rv
 */
struct bernoulli_rv_different_parameters
	{
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ans = 0;
		float prob = thrust::get<1>(t);
		float rv = thrust::get<0>(t);
		if (rv <= prob)
			{
			ans = 1;
			}
		if (rv > prob)
			{
			ans = 0; 
			}
		thrust::get<2>(t) = ans;
		}
	};
/**
 * @}
 * 
 */
#endif
