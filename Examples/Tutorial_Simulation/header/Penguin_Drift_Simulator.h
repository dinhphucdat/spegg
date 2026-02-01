#ifndef PENGUIN_SIMULATOR_H
#define PENGUIN_SIMULATOR_H

#include <Simulation_Class.h>
#include "Penguins.h"
#include <math/statistics_class.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <util/python_thrust_api.h>

namespace py = pybind11;
using StringVector = std::vector<std::string>;
using String2DVector = std::vector<std::vector<std::string>>;
using StringFloatMap = std::map<std::string, float>;
using FloatArrayVector = std::vector<py::array_t<float>>;

class Penguin_Drift_Simulator : public Simulation
	{
	public:
		Penguin_Drift_Simulator();
		~Penguin_Drift_Simulator();
		void run();
	private:
		inds_stochastic **array;
		Statistics *stats_penguins;
		void initialize_classes();

		int nspecies;
	};
#endif
