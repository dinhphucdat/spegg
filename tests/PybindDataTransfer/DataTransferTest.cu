#include <util/python_thrust_api.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <iostream>

namespace py = pybind11;

TEST(NumpyToDevice, Integer1D) {
    // start python interpreter

    py::module mymod = py::module::import("TestDataTransfer");
    const int size = 10;
    auto result = mymod.attr("getInt1D")();
    thrust::device_vector<int> d_vec_int(size);
    py::array_t<int, py::array::c_style> result_arr = result.cast<py::array_t<int, py::array::c_style>>();
    numpy_array_to_thrust_vector<int>(d_vec_int, result_arr);
    thrust::host_vector<int> h_vec_int = d_vec_int;
    // testing
    thrust::host_vector<int> expect_vtr = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    EXPECT_THAT(h_vec_int, ::testing::ContainerEq(expect_vtr));
}

TEST(NumpyToDevice, Integer2D) {
    // start python interpreter

    py::module mymod = py::module::import("TestDataTransfer");
    const int rows = 3;
    const int cols = 5;
    auto result = mymod.attr("getInt2D")();
    thrust::device_vector<int>* d_vec_int_2d = new thrust::device_vector<int>[rows];
    for (int i = 0; i < rows; i++) {
        d_vec_int_2d[i].resize(cols);
    }
    py::array_t<int, py::array::c_style> result_arr = result.cast<py::array_t<int, py::array::c_style>>();
    numpy_array_to_thrust_vector<int>(d_vec_int_2d, result_arr);
    // testing
    thrust::host_vector<int> expected_row0 = { 0, 1, 2, 3, 4 };
    thrust::host_vector<int> expected_row1 = { 5, 6, 7, 8, 9 };
    thrust::host_vector<int> expected_row2 = { 10, 11, 12, 13, 14 };
    // observed vectors
    thrust::host_vector<int> obs_row1 = d_vec_int_2d[0];
    thrust::host_vector<int> obs_row2 = d_vec_int_2d[1];
    thrust::host_vector<int> obs_row3 = d_vec_int_2d[2];
    EXPECT_THAT(obs_row1, ::testing::ContainerEq(expected_row0)) << "Int2D - Row 0 mismatch!";
    EXPECT_THAT(obs_row2, ::testing::ContainerEq(expected_row1)) << "Int2D - Row 1 mismatch!";
    EXPECT_THAT(obs_row3, ::testing::ContainerEq(expected_row2)) << "Int2D - Row 2 mismatch!";
}

TEST(NumpyToDevice, Float1D) {
    // start python interpreter

    py::module my_mod = py::module::import("TestDataTransfer");
    const int size = 10;
    auto result = my_mod.attr("getFloat1D")();
    thrust::device_vector<float> d_vec_float(size);
    py::array_t<float, py::array::c_style> result_arr = result.cast<py::array_t<float, py::array::c_style>>();
    numpy_array_to_thrust_vector<float>(d_vec_float, result_arr);
    thrust::host_vector<float> h_vec_float = d_vec_float;
    // testing
    thrust::host_vector<float> expect_vtr = { 0.0f, 1.1f, 2.2f, 3.3f, 4.4f, 5.5f, 6.6f, 7.7f, 8.8f, 9.9f };
    EXPECT_THAT(h_vec_float, ::testing::ContainerEq(expect_vtr));
}

TEST(NumpyToDevice, Float2d) {
    // start python interpreter

    py::module mymod = py::module::import("TestDataTransfer");
    const int rows = 3;
    const int cols = 5;
    auto result = mymod.attr("getFloat2D")();
    thrust::device_vector<float>* d_vec_float_2d = new thrust::device_vector<float>[rows];
    for (int i = 0; i < rows; i++) {
        d_vec_float_2d[i].resize(cols);
    }
    py::array_t<float, py::array::c_style> result_arr = result.cast<py::array_t<float, py::array::c_style>>();
    numpy_array_to_thrust_vector<float>(d_vec_float_2d, result_arr);
    // testing
    thrust::host_vector<float> expected_row0 = {  0.0,   1.1,   2.2,   3.3,   4.4 };
    thrust::host_vector<float> expected_row1 = {  5.5,   6.6,   7.7,   8.8,   9.9 };
    thrust::host_vector<float> expected_row2 = { 10.1, 11.11, 12.12, 13.13, 14.14 };
    // observed vectors
    thrust::host_vector<float> obs_row1 = d_vec_float_2d[0];
    thrust::host_vector<float> obs_row2 = d_vec_float_2d[1];
    thrust::host_vector<float> obs_row3 = d_vec_float_2d[2];
    EXPECT_THAT(obs_row1, ::testing::ContainerEq(expected_row0)) << "Float2D - Row 0 mismatch!";
    EXPECT_THAT(obs_row2, ::testing::ContainerEq(expected_row1)) << "Float2D - Row 1 mismatch!";
    EXPECT_THAT(obs_row3, ::testing::ContainerEq(expected_row2)) << "Float2D - Row 2 mismatch!";
}

int main(int argc, char** argv) {
    py::scoped_interpreter guard{};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
