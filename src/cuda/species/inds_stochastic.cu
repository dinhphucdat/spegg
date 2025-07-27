#include <species/inds_stochastic.h>

/**
 * Use @c inds_stochastic class, a subset of @c inds class, if the randomness of the initial population is desired.
 * 
 * <p>This constructor will prime up the random generator that will be used throughout the program for 
 * multiple purposes such as generate random initial pop information and other things as well.
 * </p>
 * @param size_val the initial size of the population of the species specified
 * @param maxsize_val the maximal capacity that can hold individuals of that species
 * @param num_demes the number of demes into which individuals of that species are split
 * @param species_ID_val the ID with which that specific species is tagged
 * @sa inds::inds for additional notes.
 */
inds_stochastic::inds_stochastic(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : inds(size_val, maxsize_val, ndemes, species_ID_val)
	{
/*
*
* Initialize the curand generator gen based on the seed argument, using CURAND_RNG_PSEUDO_DEFAULT. Draw 100 random uniform variables and store them in the vector rand; then discard the contents of the rand vector. This is all done to initialize the CUDA random number generator object. Note the rand vector will then be deallocated once the prime_random_number_generator scope is ended. For reasons not entirely clear, this can't seem to be done inside an external function prime_random_number_generator
*
*/
	seed = seed_val;

	int size = 100;
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gen, seed);

	//curand declarations
	thrust::device_vector<float> rand(size);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(gen, rand_ptr, size); // priming up the random number generator takes some time, get it done early.
	rand.clear();

	//Specify the indices among the phenotypes for the fitness components
	MORTALITY_PHENOTYPE_INDEX = (int) demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"];
	FECUNDITY_PHENOTYPE_INDEX = (int) demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
	}

