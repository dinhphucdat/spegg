#ifndef RANDOM_VARIABLES_H
#define RANDOM_VARIABLES_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <util/thrust_functors.h>

/**
 * @{ @name Random variables functions 
 */
/**
 * @brief Prime up the @c gen generator object to guarantee the natural randomness
 * 
 * @param gen generator object
 * @param seed Initial seed used to prime the @c gen up
 */
void prime_random_number_generator(curandGenerator_t gen, int seed);
/**
 * @brief Draw normal distribution into the @c random_variates float vector
 * 
 * @param samples_needed sample size
 * @param mean mean
 * @param stddev standard deviation
 * @param random_variates output float vector for the generated normal distribution
 * @param gen @c curandGenerator_t object
 */
void draw_gaussian(int samples_needed, float mean, float stddev, thrust::device_vector<float> &random_variates, curandGenerator_t gen);
/**
 * @brief Draw normal distribution into the @c random_variates int vector
 * 
 * This function performs the same way as @link draw_gaussian @endlink but it outputs the result as a sequence of integers instead of floating-point numbers
 * 
 * @param samples_needed sample size
 * @param mean mean
 * @param stddev standard deviation
 * @param random_variates output int vector for the generated normal distribution
 * @param gen @c curandGenerator_t object 
 * 
 * @see draw_gaussian
 */
void draw_discrete_gaussian(int samples_needed, float mean, float stddev, thrust::device_vector<int> &random_variates, curandGenerator_t gen);
/**
 * @brief Generates Poisson-distributed random variates.
 *
 * This function generates Poisson-distributed random variates using the specified lambda parameter and outputs the sequence into @c random_variates 
 * Note: CUDA 5 also implements native Poisson and Gaussian simulators, which might be used in the future.
 *
 * @param samples_needed Number of random variates to generate.
 * @param lambda The lambda parameter for the Poisson distribution.
 * @param random_variates Device vector to store the generated random variates.
 * @param gen CURAND generator.
 */
void draw_poisson(int samples_needed, float lambda, thrust::device_vector<int> &random_variates, curandGenerator_t gen);
/**
 * @brief Generates Gaussian-distributed random variates with different means and standard deviations
 * 
 * This function generates Gaussian-distributed random variates using specified means and standard deviations for each sample. It uses two sets of uniform random numbers to generate the Gaussian variates.
 * 
 * This function uses the Box-Muller transform to convert uniform random variables into Gaussian-distributed random variables.
 * 
 * @param samples_needed Number of random variates to generate
 * @param Means Device vector containing the means for the Gaussian distribution
 * @param Standard_deviations Device vector containing the standard deviations for the Gaussian distribution
 * @param random_variates Device vector to store the generated random variates
 * @param gen @c curandGenerator_t object
 */
void draw_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen);
/**
 * @brief Generates Gaussian-distributed random variates with the same mean everywhere and standard deviations
 * 
 * This function generates Gaussian-distributed random variates using specified mean and standard deviations for each sample. It uses two sets of uniform random numbers to generate the Gaussian variates.
 * 
 * @param samples_needed Number of random variates to generate
 * @param mean_value the mean value used to generate gaussian all the way across
 * @param Standard_deviations Device vector containing the standard deviations for the Gaussian distribution
 * @param random_variates Device vector to store the generated random variates
 * @param gen @c curandGenerator_t object
 * 
 * @see void draw_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen)
 */
void draw_gaussian_different_parameters(int samples_needed, float mean_value, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen);
/**
 * @brief Generates discrete Gaussian-distributed random variates with different means and standard deviations
 * 
 * @param samples_needed Number of random variates to generate
 * @param Means Device vector containing the means for the Gaussian distribution
 * @param Standard_deviations Device vector containing the standard deviations for the Gaussian distribution
 * @param random_variates Integer device vector to store the generated random variates
 * @param gen @c curandGenerator_t object
 * 
 * @see void draw_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen)
 * @sa void draw_discrete_gaussian(int samples_needed, float mean, float stddev, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
 */
void draw_discrete_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<int> &random_variates, curandGenerator_t gen);
/**
 * @brief Generates Poisson-distributed random variates with different lambda parameters
 * 
 * This function generates Poisson-distributed random variates using specified lambda parameters for each sample. It uses uniform random numbers to generate the Poisson variates
 * 
 * @param samples_needed Number of random variates to generate
 * @param lambdas Device vector containing the lambda parameters for the Poisson distribution
 * @param random_variates Device vector to store the generated random variates
 * @param gen @c curandGenerator_t object
 * 
 * @see void draw_poisson(int samples_needed, float lambda, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
 */
void  draw_poisson_different_parameters(int samples_needed, thrust::device_vector<float> &lambdas, thrust::device_vector<int> &random_variates, curandGenerator_t gen);
/**
 * @brief Generates Bernoulli-distributed random variates, with integers.
 *
 * This function generates Bernoulli-distributed random variates using the specified probability parameter.
 * It uses uniform random numbers to generate the Bernoulli variates.
 *
 * @param samples_needed Number of random variates to generate.
 * @param probability The probability parameter for the Bernoulli distribution.
 * @param random_variates Device vector to store the generated random variates.
 * @param gen @c curandGenerator_t object
 */
void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<int> &random_variates, curandGenerator_t gen);
/**
 * @brief Generates Bernoulli-distributed random variate, with floats.
 *
 * @see void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
 *
 * @param samples_needed Number of random variates to generate.
 * @param probability The probability parameter for the Bernoulli distribution.
 * @param random_variates Device vector to store the generated random variates.
 * @param gen @c curandGenerator_t object
 */
void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<float> &random_variates, curandGenerator_t gen);

/**
 * @brief Generates Bernoulli-distributed random variates with different probability parameters.
 *
 * This function generates Bernoulli-distributed random variates using specified probability parameters for each sample.
 * It uses uniform random numbers to generate the Bernoulli variates.
 *
 * @param samples_needed Number of random variates to generate.
 * @param probabilities Device vector containing the probability parameters for the Bernoulli distribution.
 * @param random_variates Device vector to store the generated random variates.
 * @param gen @c curandGenerator_t object.
 * 
 * @see void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<float> &random_variates, curandGenerator_t gen)
 */
void draw_bernoulli_different_parameters(int samples_needed, thrust::device_vector<float> &probabilities, thrust::device_vector<int> &random_variates, curandGenerator_t gen);

/**
 * @}
 */

#endif

