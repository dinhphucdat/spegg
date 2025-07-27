#ifndef EGGS_NEONATES_H
#define EGGS_NEONATES_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <species/inds_stochastic.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <species/deme_specific_data_class.h>
#include <util/amplify.h>
#include <math/mating_thrust_prob_table_demes.h>
#include <math/random_variables_functions.h>

/**
 * @brief Responsible for inheritting the generic materials across generations
 * @author xxx
 */
class EggsNeonates 
	{
	public:
		/**
		 * @brief Construct a new @c EggsNeonates object
		 * 
		 * This constructor initializes the EggsNeonates object with the given species and kids_per_mom
		 * 
		 * @param species pointer to the object of type @link inds @endlink, representing species
		 * @param kids_per_mom how many kids per female
		 * 
		 */
		EggsNeonates(inds_stochastic *species, thrust::device_vector<int> &kids_per_mom);
		/**
		 * @brief Inherits the genotypes to the offspring
		 * 
		 * This function takes in the probability of each individual in the @c species being either mom 
		 * or dad and manages to inherit the genotypes so that the offspring can have the traits reflective 
		 * of the traits of their parents. Last but not least, it performs mutation as specified in the 
		 * deme configuration file and finally copies back the offspring's genotypes to 
		 * the @c mgenotype and @c fgenotype of the @c species object
		 * 
		 * @param probability_individuals_become_mothers probability of each individual being a mother
		 * @param probability_individuals_become_fathers probability of each individual being a father
		 */
		void inherit_genotypes(thrust::device_vector<float> &probability_individuals_become_mothers,
				       thrust::device_vector<float> &probability_individuals_become_fathers);
		/// @brief the size of the population, reflected in @link EggsNeonates::species @endlink before genetic inheritance
		int previous_pop_size;
		/// @brief total number of children being born after the inheritance process
		int Total_Number_of_Neonates;
		/// @brief number of newborns per deme
		thrust::device_vector<int> Neonates_per_Deme;
		/// @brief the list of demes that any newborn belongs to
		thrust::device_vector<int> kids_deme; 

	protected:
		/// @brief the @link inds @endlink object representing individuals of a species
		class inds *species;
		/// @brief random generator, not meant for modification. This is the fixed object
		curandGenerator_t gen;
		/// @brief number of loci, number of phenotypes, and number of demes in the population
		int nloci, nphen, Num_Demes;
		/// @brief Magnitudes of mutation, specified in @c deme_config.txt
		thrust::device_vector<float> mutation_magnitude;
		/// @brief mutation rates, specfied in @c deme_config.txt
		thrust::device_vector<float> mutation_rate;
		/**
		 * @brief Filters the number of offspring per deme to fit in the carrying capacity of the population
		 * 
		 * @param subpopParameters the object of @link DemeSettings @endlink in which species parameters are stored
		 * @param everybodys_deme individual's demes
		 * @param kids_per_mom number of kids per mom
		 * @param current_deme_sizes current sizes of demes
		 * @param maximum_deme_sizes carrying capacities of demes
		 * 
		 * @note @c protected function, meant for extending this base class or override
		 */
		void Determine_Neonate_Population_Sizes(DemeSettings *subpopParameters,
							thrust::device_vector<int> &everybodys_deme,
							thrust::device_vector<int> &kids_per_mom,			
							thrust::device_vector<int> &current_deme_sizes,
							thrust::device_vector<int> &maximum_deme_sizes);
		
		/**
		 * @brief Calculate the offspring's genotypes at the maternally inherited loci.
		 *
		 * This function determines the mother of each offspring using <b>Alias sampling method</b> and generates a haploid gamete from the mother to be grafted onto the offspring. 
		 * Note that the result of this function will be copied back to the @c mgenotype and @c fgenotype 
		 * as a continuation to the genotypes of individuals already existed in the population
		 *
		 * @param probability_individuals_become_mothers Reference to a @c thrust::device_vector<float> containing the probability of individuals becoming mothers.
		 * @param mgenotype Pointer to a @c thrust::device_vector<float> representing the maternal genotype.
		 * @param fgenotype Pointer to a @c thrust::device_vector<float> representing the paternal genotype.
		 *
		 * @note Please also call @link EggsNeonates::get_paternally_derived_genotype @endlink along with this maternal inheritance
		 * @see EggsNeonates::get_paternally_derived_genotype
		 */
		void get_maternally_derived_genotype(thrust::device_vector<float> &probability_individuals_become_mothers, thrust::device_vector<float> *&mgenotype, thrust::device_vector<float> *&fgenotype);
		/**
		 * @brief Calculate the offspring's genotypes at the paternally inherited loci.
		 *
		 * This function determines the father of each offspring using <b>Alias sampling method</b> and generates a haploid gamete from the father to be grafted onto the offspring. 
		 * Note that the result of this function will be copied back to the @c mgenotype and @c fgenotype 
		 * as a continuation to the genotypes of individuals already existed in the population
		 *
		 * @param probability_individuals_become_fathers Reference to a @c thrust::device_vector<float> containing the probability of individuals becoming fathers.
		 * @param mgenotype Pointer to a @c thrust::device_vector<float> representing the maternal genotype.
		 * @param fgenotype Pointer to a @c thrust::device_vector<float> representing the paternal genotype.
		 *
		 * @note Please also call @link EggsNeonates::get_maternally_derived_genotype @endlink along with this paternal inheritance
		 * @see EggsNeonates::get_maternally_derived_genotype
		 */
		void get_paternally_derived_genotype(thrust::device_vector<float> &probability_individuals_become_fathers, thrust::device_vector<float> *&mgenotype, thrust::device_vector<float> *&fgenotype);

		
		void egg_mortality(DemeSettings *subpopParameters);
		/**
		 * @brief Copies either paternal or maternal genomes to each newborn based on random value generated compared to the recombination rate. See @link EggsNeonates::recomb_rate @endlink
		 * 
		 * @param rand random values in float
		 * @param parent parent indices
		 * @param parity 0 or 1, the toggle mechanism to determine if paternal or maternal chromosomes would be compied into the newborns' ones
		 * @param parents_fgenotype maternal genotypes, selected from @link inds::fgenotype @endlink during marental selection
		 * @param parents_mgenotype paternal genotypes, selected from @link inds::mgenotype @endlink during parental selection
		 * @param kids_genotype reference to kids genotype
		 * @param locus_ID the id of locus of interest
		 * 
		 * @see EggsNeonates::recomb_rate
		 */
		void recombine(thrust::device_vector<float> &rand,
				     thrust::device_vector<int> &parent,
				     thrust::device_vector<int> &parity,
				     thrust::device_vector<float> *&parents_fgenotype,
				     thrust::device_vector<float> *&parents_mgenotype,
				     thrust::device_vector<float> &kids_genotype,
				     int locus_ID);
		/// @brief recombination rates
		thrust::host_vector<float> recomb_rate;
		/**
		 * @brief Performs mutation to each locus based on mutation rates
		 * 
		 * Determine, for each offspring's locus, whether there will be a mutation at that locus, and 
		 * if so what the magnitude of that mutation will be and how that changes the offspring's allelic 
		 * value. The current behavior assumes mutations are gaussian about the parental allelic value, and 
		 * that the mutation parameters (mutation rate and sd of gaussian) vary by deme. Possible expansions 
		 * include allowing the mutation rate itself to be an individual-specific phenotype, or alternative 
		 * mutational models (e.g., point mutations that have a categorical rather than quantitative effect.
		 * 
		 * @param parents_fgenotype genotypes of mothers
		 * @param parents_mgenotype genotypes of fathers
		 * 
		 */
		void mutate(thrust::device_vector<float> *&parents_fgenotype,thrust::device_vector<float> *&parents_mgenotype);
		/**
		 * @brief Superseded by the @link GenotypePhenotypeMap::create_genotype_phenotype_map @endlink
		 * @deprecated
		 */
		void prepare_genotype_phenotype_map();
		/**
		 * @brief Copies processed information back to the @c species object
		 * 
		 */
		void integrate_kids();
	};

// make sure there are no more kids than spaces available


/**
 * @struct adjust_kids_functor
 * @brief Ensures that the number of kids in a subpopulation does not exceed the available capacity.
 *
 * This functor is designed to be used with Thrust algorithms. It adjusts the number of kids
 * in a subpopulation to ensure that the total population (existing + new kids) does not exceed
 * the carrying capacity of the subpopulation.
 *
 * The input is a tuple with the following elements:
 * - `0`: Number of kids to be added to the subpopulation.
 * - `1`: Current size of the subpopulation.
 * - `2`: Carrying capacity of the subpopulation.
 *
 * If the sum of the current population and the number of kids exceeds the carrying capacity,
 * the number of kids is reduced accordingly. If the subpopulation is already over capacity,
 * the number of kids is set to zero.
 */
struct adjust_kids_functor
	{
	
	/**
     * @brief Adjusts the number of kids to fit within the subpopulation's carrying capacity.
     *
     * @tparam tuple A Thrust tuple containing three elements:
     *         - `0`: int& kids (modifiable)
     *         - `1`: const int& current_population
     *         - `2`: const int& capacity
     *
     * @param t The tuple containing the kids count, current population, and capacity.
     */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
			if (thrust::get<0>(t) + thrust::get<1>(t) > thrust::get<2>(t)) 
			{
				if (thrust::get<2>(t) - thrust::get<1>(t) >= 0 )
					{
					thrust::get<0>(t) = thrust::get<2>(t) - thrust::get<1>(t);
					}
				if (thrust::get<2>(t) - thrust::get<1>(t) < 0 ) // if too crowded, no kids.
					{
					thrust::get<0>(t) = 0;
					}
			}
		}
	};


/**
 * @struct recombination_functor
 * @brief Simulates genetic recombination between two parent genotypes.
 *
 * This functor is intended for use with Thrust algorithms. It models the recombination
 * process by selecting a genotype from either the female or male parent based on a 
 * recombination rate and a parity value.
 *
 * The functor uses a tuple with the following elements:
 * - `0`: Index of the parent individual.
 * - `1`: A random float value used to determine if recombination occurs.
 * - `2`: Parity value (used to alternate between parents).
 * - `3`: Output genotype value (to be set based on recombination).
 *
 * If the random number is less than the recombination rate, the parity is toggled.
 * The genotype is then selected from the female or male parent based on the parity.
 */
struct recombination_functor
{
	float *fgenotype, *mgenotype;
	float recomb_rate;

	/**
     * @brief Constructor to initialize genotype pointers and recombination rate.
     * @param fgene Pointer to female genotype array.
     * @param mgene Pointer to male genotype array.
     * @param rate Recombination rate (probability).
     */
	recombination_functor(float *fgene, float *mgene, float rate) : fgenotype(fgene), mgenotype(mgene), recomb_rate(rate)
	{};
	
	
	/**
     * @brief Applies the recombination logic to a tuple of parent data.
     *
     * @tparam tuple A Thrust tuple containing:
     *         - `0`: int parent_index
     *         - `1`: float random_value
     *         - `2`: int& parity (modifiable)
     *         - `3`: float& genotype (modifiable)
     *
     * @param t The tuple containing parent index, random value, parity, and genotype.
     */
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int parent_index = thrust::get<0>(t);

		if (thrust::get<1>(t) < recomb_rate) 
			{
			thrust::get<2>(t) = (thrust::get<2>(t) + 1)%2;
			}
		if (thrust::get<2>(t) == 0)
			thrust::get<3>(t) = fgenotype[parent_index];
		if (thrust::get<2>(t) == 1)
			thrust::get<3>(t) = mgenotype[parent_index];

		float answer = thrust::get<3>(t);
		}
};
#endif
