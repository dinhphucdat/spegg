#include <species/add_kids/assortative_mating_neonates_class.h>

/* One problem with a separate neonates class is that all this stuff has to get copied back into example species. It might just not be worth it depending on the performance cost. */

#define NLOCI 50

float recomb_array_assort[NLOCI] = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,0.5f, 0.5f, 0.5f, 0.5f, 0.5f};

// TODO: The parent class has up to 2 constructors, but here I only see one which is a longer parameter list one. So think in order to 
// ensure that we can initialize the parent class properly, we need to create another two-parameter constructor as the parent class did 
// which takes in the species object and kids_per_mom vector!
Assortative_mating_neonates::Assortative_mating_neonates(inds_stochastic* species, thrust::device_vector<int> &kids_per_mom) : EggsNeonates(species, kids_per_mom) 
	{
	// how many subpopulations depends on how many demes there are
	
	// the number of the current population is the sum of all the deme sizes
	this->Num_Subpopulations = species->demeParameters->check_number_of_demes();
	current_pop_size = species->size;
	thrust::device_vector<int> pair_populations = species->deme;
	// the pairs per deme cannot be the deme vector from species object!!!!
	pairs_per_deme.resize(Num_Subpopulations);
	
	// Determine how many pairs are in each subpopulation
	thrust::counting_iterator<int> search_begin(0);
	thrust::device_vector<int> temp_Subpop_sizes;
	// TODO: Note that the inds_stochastic object has `deme_size` vector, so is it a need to 
	// recalculate the number of individuals in each deme?
	temp_Subpop_sizes.resize(Num_Subpopulations);

	thrust::upper_bound(pair_populations.begin(), pair_populations.end(),
                      search_begin, search_begin + Num_Subpopulations,
                      temp_Subpop_sizes.begin());

	thrust::adjacent_difference(temp_Subpop_sizes.begin(), temp_Subpop_sizes.end(),
                              pairs_per_deme.begin());
		
	}

void Assortative_mating_neonates::inherit_genotypes_by_pair(thrust::device_vector<float> &probability_pair_becomes_parents,
				thrust::device_vector<int> &fathers_list,
				thrust::device_vector<int> &mothers_list,
				thrust::device_vector<float> *&fgenotype,
				thrust::device_vector<float> *&mgenotype)
	{

	mothers_chosen.resize(Total_Number_of_Neonates);
	fathers_chosen.resize(Total_Number_of_Neonates);

	get_mating_pair(probability_pair_becomes_parents, fathers_list, mothers_list);

	get_maternally_derived_genotype_deterministic(mothers_chosen, mgenotype, fgenotype);

	get_paternally_derived_genotype_deterministic(fathers_chosen, mgenotype, fgenotype);

	mutate(mgenotype, fgenotype);
	}

void Assortative_mating_neonates::get_mating_pair(thrust::device_vector<float> &probability_pair_becomes_parents,
						   thrust::device_vector<int> &fathers_list,
						   thrust::device_vector<int> &mothers_list
						  )
	{
	mating_ThrustProbTable_demes at;
	thrust::device_vector<int> pair_index(Total_Number_of_Neonates);
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
/*
	Feed reproductive probablity into the setup of the alias table.
	Draw from the alias table to determine mothers.
*/	
	at.setup(probability_pair_becomes_parents.begin(), probability_pair_becomes_parents.begin() + previous_pop_size);

	curandGenerateUniform(this->gen, rand_ptr, Total_Number_of_Neonates);
 
	at.determine_key_offsets( Num_Subpopulations, species->deme_sizes );
 
	at.adjust_randoms(rand.begin(), rand.end(), kids_deme.begin(), kids_deme.end());
 
	at.draw(rand.begin(), rand.end(), pair_index.begin());

 		
	// thrust::gather(pair_index.begin(), pair_index.end(), fathers_list.begin(), fathers_chosen.begin());

 
	thrust::gather(pair_index.begin(), pair_index.end(), species->deme.begin(), mothers_chosen.begin());

	}

/* use the functions blah_blah_deterministic() if the parents have already been chosen and you just need to copy genotypes */
	
void Assortative_mating_neonates::get_maternally_derived_genotype_deterministic(thrust::device_vector<int> &mother_index,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype)
	{
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(this->gen, rand_ptr, Total_Number_of_Neonates);
	
	thrust::device_vector<int> parity(Total_Number_of_Neonates);
	thrust::fill(parity.begin(), parity.end(), 0);

	for (int i = 0 ; i < nloci ; i++) 
		{
		curandGenerateUniform(this->gen, rand_ptr, Total_Number_of_Neonates);
        recombine(rand, mother_index, parity, fgenotype, mgenotype, fgenotype[i], i);
		}
	}

void Assortative_mating_neonates::get_paternally_derived_genotype_deterministic(thrust::device_vector<int> &father_index,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype)
	{
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(this->gen, rand_ptr, Total_Number_of_Neonates);
	
	//Reset parity to zeroes
	thrust::device_vector<int> parity(Total_Number_of_Neonates);
	thrust::fill(parity.begin(), parity.end(), 0);

	//Recombination for mgenotype
	for (int i = 0 ; i < nloci ; i++) 
		{
		curandGenerateUniform(this->gen, rand_ptr,Total_Number_of_Neonates);
		recombine(rand, father_index, parity, fgenotype, mgenotype, mgenotype[i], i);
		}
	}

void Assortative_mating_neonates::record_parents(thrust::device_vector<int> &maternal_id, thrust::device_vector<int> &paternal_id, thrust::device_vector<int> &ids)
	{
	thrust::gather(mothers_chosen.begin(), mothers_chosen.begin() + Total_Number_of_Neonates, ids.begin(), maternal_id.begin() + current_pop_size);
	thrust::gather(fathers_chosen.begin(), fathers_chosen.begin() + Total_Number_of_Neonates, ids.begin(), paternal_id.begin() + current_pop_size);
	}
