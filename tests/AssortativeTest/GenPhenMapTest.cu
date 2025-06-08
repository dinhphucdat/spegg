// #include <species/add_kids/genotype_phenotype_map.h>
#include <GenPhenMapTest.h>


//
// Include your mySpecies_genotype_phenotype_maps header file here
//

GenotypePhenotypeMap *GenotypePhenotypeMap::create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids)
{
	if (phenotype_index == species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]) {
		// std::cout << "start creating factory fecundity" << std::endl;
		return new FecundityPhenotype(species, phenotype_index, index_case, num_kids);
	}
	
	if (phenotype_index == species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"]) {
		// std::cout << "start creating factory mortality" << std::endl;
		return new MortalityPhenotype(species, phenotype_index, index_case, num_kids);
	}
	
	if (phenotype_index == species->demeParameters->species_specific_values["WING_COLOR_INDEX"]) {
		// std::cout << "start creating factory wing" << std::endl;
		return new WingColorPhenotype(species, phenotype_index, index_case, num_kids);
	}
	// std::cout << "end creating factory gen phen map" << std::endl;
}

void FecundityPhenotype::calculate_phenotype(inds* species)
{
    // wrap the parameters from vectors to arrays:
    float *constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
    float *coefficient0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);
    float *coefficient1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
    //Instantiate functor
    FecundityCalculator fecundityCal_functor(constants, coefficient0, coefficient1);
    //Perform genotype-phenotype map operation with for_each.

    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[0].begin() + index_case, species->mgenotype[0].begin() + index_case, species->fgenotype[1].begin() + index_case, species->mgenotype[1].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case)), thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[0].begin() + index_case + num_kids, species->mgenotype[0].begin() + index_case + num_kids, species->fgenotype[1].begin() + index_case + num_kids, species->mgenotype[1].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids)), fecundityCal_functor);

}

void MortalityPhenotype::calculate_phenotype(inds* species)
{
    // wrap the parameters from vectors to arrays:
    float* constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
    float* coefficient = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);

    // initializes functor
    MortalityCalculator mortalityCalc_functor(constants, coefficient);

    // perform calculation using functor
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[0].begin() + index_case, species->mgenotype[0].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case)), thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[0].begin() + index_case + num_kids, species->mgenotype[0].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids)), mortalityCalc_functor);

    
}

void WingColorPhenotype::calculate_phenotype(inds* species)
{
    float* constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
    float* coefficients = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);

    WingColorCalculator wingColorCalc_functor(constants, coefficients);

    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[2].begin() + index_case, species->mgenotype[2].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case)), thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[2].begin() + index_case + num_kids, species->mgenotype[2].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids)), wingColorCalc_functor);
}