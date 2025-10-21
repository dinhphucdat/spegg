#ifndef MATING_THRUST_PROB_TABLE_DEMES_H
#define MATING_THRUST_PROB_TABLE_DEMES_H

#include <math/thrust_prob_table.h>

/**
 * @brief A derived class of @link ThrustProbTable @endlink that adjusts random numbers based on deme sizes for mating purposes
 * 
 */
class mating_ThrustProbTable_demes : public ThrustProbTable
	{
	public:
		/**
		 * @brief Adjusts random numbers based on deme sizes so that they can be used for drawing individuals for mating
		 * 
		 * @param uniform_begin iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param inds_deme_begin iterator pointing to the beginning of the range of deme sizes in the inds object.
		 * @param inds_deme_end iterator pointing to the end of the range of deme sizes in the inds object.
		 */
		void adjust_randoms(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
		/**
		 * @brief Adjusts random numbers based on fixed offsets so that they can be used for drawing individuals for mating
		 * 
		 * @param uniform_begin iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param inds_deme_begin iterator pointing to the beginning of the range of deme sizes in the inds object.
		 * @param inds_deme_end iterator pointing to the end of the range of deme sizes in the inds object.
		 */
		void adjust_randoms_fixed_offsets(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
		/**
		 * @brief Determines the key offsets based on the histogram of key types
		 * 
		 * This function computes the starting index (offset) for each key type based on the provided histogram of key types.
		 * The offsets are stored in the `key_offsets` member variable.
		 * 
		 * @param number_of_key_types The total number of distinct key types.
		 * @param key_histogram_vector A device vector containing the histogram of key types, where each element represents the count of individuals for that key type.
		 */
		void determine_key_offsets(int number_of_key_types, thrust::device_vector<int> &key_histogram_vector);
		/// @brief A device vector storing the starting index (offset) for each key type
		thrust::device_vector<int> key_offsets;
	};

/**
 * @brief A derived class of @link ThrustProbTableDouble @endlink that adjusts random numbers based on deme sizes for mating purposes
 * 
 */
class mating_ThrustProbTable_demes_Double : public ThrustProbTableDouble
	{
	public:
		/**
		 * @brief Adjusts random numbers based on deme sizes so that they can be used for drawing individuals for mating
		 * 
		 * @param uniform_begin iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param inds_deme_begin iterator pointing to the beginning of the range of deme sizes in the inds object.
		 * @param inds_deme_end iterator pointing to the end of the range of deme sizes in the inds object.
		 */
		void adjust_randoms(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
		/**
		 * @brief Adjusts random numbers based on fixed offsets so that they can be used for drawing individuals for mating
		 * 
		 * @param uniform_begin iterator pointing to the beginning of the range of uniformly distributed random numbers.
		 * @param uniform_end iterator pointing to the end of the range of uniformly distributed random numbers.
		 * @param inds_deme_begin iterator pointing to the beginning of the range of deme sizes in the inds object.
		 * @param inds_deme_end iterator pointing to the end of the range of deme sizes in the inds object.
		 */
		void adjust_randoms_fixed_offsets(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
		/**
		 * @brief Determines the key offsets based on the histogram of key types
		 * 
		 * @param number_of_key_types The total number of distinct key types.
		 * @param key_histogram_vector A device vector containing the histogram of key types, where each element represents the count of individuals for that key type.
		 */
		void determine_key_offsets(int number_of_key_types, thrust::device_vector<int> &key_histogram_vector);
		/// @brief A device vector storing the starting index (offset) for each key type
		thrust::device_vector<int> key_offsets;
	};

#endif
