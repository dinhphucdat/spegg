#include <species/add_kids/genetic_deme_specific_data.h>

DemeGeneticsSettings::DemeGeneticsSettings(const char *filename, int species_ID)
	{
	read_in_data(filename, species_ID);
	}


void DemeGeneticsSettings::read_in_data(const char *filename, int species_ID)
	{
	Config cfg;
	try
		{
		cfg.readFile(filename);
		}
	catch(const FileIOException &fioex)
		{
		std::cerr << "No genetic config file." << std::endl;
		}
	
	try
		{
		Number_of_Demes = cfg.lookup("number_of_demes");
		}
 	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'number_of_demes' setting in configuration file." << std::endl;
		}

	const Setting& root = cfg.getRoot();
	const Setting &species_genetic_specification = root["species_genetics"];

	try
		{
		species_genetic_specification[species_ID].lookupValue("number_of_loci", Number_of_Loci);

		loci_names.resize(Number_of_Loci);
		if (species_genetic_specification[species_ID]["loci_names"].getLength() != Number_of_Loci)
			std::cout << "Length of loci names differs from number of loci." << std::endl;

		for (int i=0; i < species_genetic_specification[species_ID]["loci_names"].getLength(); i++)
			{
			loci_names[i] = species_genetic_specification[species_ID]["loci_names"][i].c_str();
			}
		}

	catch(const SettingNotFoundException &nfex)
		{
		// Toggle if specifying the names of loci becomes important
		std::cout << "warning: No loci name. Use int indices to pass in genetics details " << std::endl;
		}

	const Setting &recombination_rates_config = species_genetic_specification[species_ID]["recombination_rates"];

	if (species_genetic_specification[species_ID]["recombination_rates"].getLength() != Number_of_Loci)
		std::cout << "Length of recombination_rates differs from number of loci.\n" << std::endl;

	recombination_rates.resize(Number_of_Loci);

	for (int i=0; i < Number_of_Loci; i++)
		{
		float val = 0;
		val = recombination_rates_config[i];
		recombination_rates[i] = val;
		}

	// Read in the locus-specific information (for now, just mutation rates)
	const Setting &locus_specification = species_genetic_specification[species_ID]["locus_specifications"];
	
	if (Number_of_Loci > 0)
		{
		deme_specific_mutation_rates = new thrust::device_vector<float>[Number_of_Loci];
		deme_specific_mutation_magnitudes = new thrust::device_vector<float>[Number_of_Loci];

		for (int i=0; i < Number_of_Loci; i++)
			{
			const Setting &mutation_settings = locus_specification[i]["mutation_parameters"];
			deme_specific_mutation_rates[i].resize(Number_of_Demes);
			deme_specific_mutation_magnitudes[i].resize(Number_of_Demes);

			for (int j=0; j < Number_of_Demes; j++)
				{
				const Setting &deme_values = mutation_settings[j]; 
				float val = 0;
				deme_values.lookupValue("MUTATION_RATE",val);
				deme_specific_mutation_rates[i][j] = val;
				deme_values.lookupValue("MUTATION_MAGNITUDE",val);
				deme_specific_mutation_magnitudes[i][j] = val;
				}
			}
		}
	
	else
		{
		std::cout << "warning: your species has no loci associated with it" << std::endl;
		}


	// Read in the genotype phenotype maps
	species_genetic_specification[species_ID].lookupValue("number_of_phenotypes", Number_of_Phenotypes);

	if (Number_of_Phenotypes > 0)
		{
		phenotype_names.resize(Number_of_Phenotypes);

		for (int i = 0; i < species_genetic_specification[species_ID]["phenotype_names"].getLength(); i++)
			{
			phenotype_names[i] = species_genetic_specification[species_ID]["phenotype_names"][i].c_str();
			}

		phen_gen_map_parm = new GenotypePhenotypeMapParameters *[Number_of_Phenotypes];

		for (int i=0; i < Number_of_Phenotypes; i++)
			phen_gen_map_parm[i] = new GenotypePhenotypeMapParameters(filename, species_ID, i, loci_names);
		}
	
	else
		{
		std::cout << "Warning: your species has no phenotypes associated with it" << std::endl;
		}
	}

thrust::device_ptr<float> DemeGeneticsSettings::get_mutation_rates_ptr(int locus_index)
	{
	return(&deme_specific_mutation_rates[locus_index][0]);
	}

thrust::device_ptr<float> DemeGeneticsSettings::get_mutation_magnitudes_ptr(int locus_index)
	{
	return(&deme_specific_mutation_magnitudes[locus_index][0]);
	}

// --------------------- DESTRUCTOR ----------------------- //

DemeGeneticsSettings::~DemeGeneticsSettings() {
	// clean up the mutation magnitudes and rates
	delete[] deme_specific_mutation_rates;
	delete[] deme_specific_mutation_magnitudes;
	// clean up GenotypePhenotypeMapParameters**
	for (int i = 0; i < Number_of_Phenotypes; i++) {
		delete phen_gen_map_parm[i];
		phen_gen_map_parm[i] = nullptr;
	}
	delete[] phen_gen_map_parm;
}

// --------------------- NEW FUNCTIONALITY - READ FROM PYTHON ------------------------ //
DemeGeneticsSettings::DemeGeneticsSettings(
	const int&									 speciesID, 
	const std::vector<std::string>& 			 phenotypeNames, 
	const std::vector<std::vector<std::string>>& genPhenParameterNamesAllPhenotypes, 
	const std::vector<py::array_t<float>>& 		 demeSpecificPhenParametersAllPhenotypes
	const std::vector<std::string>& 			 lociNames, 
	const py::array_t<float>& 					 recombinationRates, 
	const py::array_t<float>&  					 demeSpecificMutationRates, 
	const py::array_t<float>&					 demeSpecificMutationMagnitudes
) {
	processGenotypicInfo(
		lociNames, 
		recombinationRates, 
		demeSpecificMutationRates, 
		demeSpecificMutationMagnitudes
	);

	processPhenotypicInfo(
		speciesID, 
		phenotypeNames, 
		genPhenParameterNamesAllPhenotypes, 
		demeSpecificPhenParametersAllPhenotypes
	);
}

void DemeGeneticsSettings::processGenotypicInfo(
	const std::vector<std::string>&  lociNames, 
	const py::array_t<float>&  	 	 recombinationRates, 
	const py::array_t<float>& 	 	 demeSpecificMutationRates, 
	const py::array_t<float>&		 demeSpecificMutationMagnitudes
) {
	// loci names
	this->Number_of_Loci = lociNames.size();
	loci_names.resize(Number_of_Loci);
	std::copy(lociNames.begin(), lociNames.end(), loci_names.begin());

	// recombination rates
	py::array_t<float, py::array::c_style> recomRatesArr = (
		recombinationRates.cast<py::array_t<float, py::array::c_style>>()
	);
	const py::buffer_info buf = recomRatesArr.request();
	if (buf.shape.size() != 1) {
		throw std::runtime_error("The recombination rates array should only have one dimension");
	}
	if (buf.shape[0] != Number_of_Loci) {
		throw std::runtime_error("The recombination rates array should have the size of Number of Loci");
	}
	numpy_array_to_thrust_vector<float>(
		recombination_rates, 
		recomRatesArr
	);

	// deme-specific mutation rates
	py::array_t<float, py::array::c_style> demeSpecMutRatesArr = (
		demeSpecificMutationRates.cast<py::array_t<float, py::array::c_style>>()
	);
	const py::buffer_info bufRate = demeSpecMutRatesArr.request();
	if (bufRate.shape.size() != 2) {
		throw std::runtime_error("The mutation rates array should have two dimensions: " + 
			"(Number of Loci, Number of Demes)"
		);
	}
	if (bufRate.shape[0] != Number_of_Loci) {
		throw std::runtime_error("The mutation rates outer dimension must match number of Loci");
	}
	Number_of_Demes = static_cast<int>(buf.shape[1]);
	deme_specific_mutation_rates = new thrust::device_vector<float>[Number_of_Loci];
	numpy_array_to_thrust_vector<float>(
		deme_specific_mutation_rates, 
		demeSpecMutRatesArr
	);

	// deme-specific mutation magnitudes
	py::array_t<float, py::array::c_style> demeSpecMutMag = (
		demeSpecificMutationMagnitudes.cast<py::array_t<float, py::array::c_style>>()
	);
	const py::buffer_info bufMag = demeSpecMutMag.request();
	if (bufMag.shape.size() != 2) {
		throw std::runtime_error(
			"The mutation magnitude array should only have two dimensions: " + 
			"(Number_of_Loci, Number_of_Demes)"
		);
	}
	if (bufRate.shape[0] != Number_of_Loci) {
		throw std::runtime_error("The mutation magnitude outer dimension must match number of loci");
	}
	if (bufRate.shape[1] != Number_of_Demes) {
		throw std::runtime_error("The mutation magnitude inner dimension must match number of demes");
	}
	deme_specific_mutation_magnitudes = new thrust::device_vector<float>[Number_of_Loci];
	numpy_array_to_thrust_vector<float>(
		deme_specific_mutation_magnitudes, 
		demeSpecMutMag
	);
}

void DemeGeneticsSettings::processPhenotypicInfo(
	const int& 									 speciesID, 
	const std::vector<std::string>& 			 phenotypeNames, 
	const std::vector<std::vector<std::string>>& genPhenParameterNamesAllPhenotypes, 
	const std::vector<py::array_t<float>>& 		 demeSpecificPhenParametersAllPhenotypes
) {
	// number of phenotypes, phenotype names
	this->Number_of_Phenotypes = phenotypeNames.size();
	phenotype_names.resize(Number_of_Phenotypes);
	std::copy(phenotypeNames.begin(), phenotypeNames.end(), phenotype_names.begin());

	// initiate phenotypic information
	phen_gen_map_parm = new GenotypePhenotypeMapParameters*[Number_of_Phenotypes];
	for (int i = 0; i < Number_of_Phenotypes; i++) {
		phen_gen_map_parm[i] = new GenotypePhenotypeMapParameters(
			speciesID, 
			i, 
			genPhenParameterNamesAllPhenotypes, 
			demeSpecificPhenParametersAllPhenotypes
		);
	}
}

// --------------------- END OF THIS NEW FUNCTIONALITY ------------------------------- //