#ifndef PYTHON_THRUST_API_H
#define PYTHON_THRUST_API_H

#include <thrust/device_vector.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <thrust/copy.h>
#include <cuda_runtime.h>
#include <vector>
#include <stdexcept> // for std::runtime_error


namespace py = pybind11;

/**
 * @brief Converts a 2D numpy array to a pointer to thrust::device_vector<T>'s.
 * 
 * @tparam T 
 * @param target the thrust::device_vector<T>* to populate. Remember, the dimension size must match
 * @param py_array numpy array to convert
 */
template <typename T>
void numpy_array_to_thrust_vector(
    thrust::device_vector<T>*& target, 
    const py::array_t<T, py::array::c_style>& py_array
) {
    // request the python's ndarray's data buffer
    py::buffer_info buf = py_array.request();
    // metadata: size, first dimension, second dimension
    const size_t size = buf.size;
    if (size == 0) {
        throw std::runtime_error("Numpy array is empty!");
    }
    // remember that buf.size == buf.shape[0] * buf.shape[1] * ... * buf.shape[n]
    // because ndarray is represented as contiguous array in memory
    std::vector<ssize_t> shape = buf.shape;
    if (shape.size() != 2) {
        throw std::runtime_error("Only 2D numpy arrays are supported!");
    }
    const size_t num_rows = shape[0];
    const size_t num_cols = shape[1];
    // get the pointer of the buffer
    T* buf_ptr = static_cast<T*>(buf.ptr);
    // intialize an empty thrust::device_vector*
    // if the driver decides to initialize a thrust::device_vector pointer themselves, 
    // they need to make sure that it is of the correct size
    if (!target) {
        target = new thrust::device_vector<T>[size];
    }
    // get the raw pointer of the device_vector's data
    
    // copy data from host to device
    // remember, the size is counted in bytes
    {
        py::gil_scoped_release release;
        for (size_t r = 0; r < num_rows; r++) {
            // resize the rth device_vector if necessary
            target[r].resize(num_cols);
            // raw pointer for the rth row
            T* d_ptr_row = thrust::raw_pointer_cast(target[r].data());
            // increment the base address of buf_ptr
            T* buf_ptr_row = buf_ptr + r * num_cols;
            // copy the rth row
            cudaError_t err = cudaMemcpy(
                d_ptr_row, buf_ptr_row, num_cols * sizeof(T), cudaMemcpyHostToDevice
            );

            if (err != cudaSuccess) {
                throw std::runtime_error(
                    std::string("cudaMemcpy failed in transfering data at row ") + 
                    std::to_string(r + 1) + ": " + 
                    cudaGetErrorString(err)
                );
            }
        }
    }
}

/**
 * @brief Converts a 1D numpy array to a thrust::device_vector<T>.
 * 
 * @tparam T 
 * @param target the thrust::device_vector<T> to populate. Remember, the dimension size must match
 * @param py_array the numpy array to convert
 */
template <typename T>
void numpy_array_to_thrust_vector(
    thrust::device_vector<T>& target, 
    const py::array_t<T, py::array::c_style>& py_array
) {
    // request the python's array's data buffer
    py::buffer_info buf = py_array.request();
    // metadata: dimension - ONLY ACCEPT 1D
    const std::vector<ssize_t> shape = buf.shape;
    if (shape.size() != 1) {
        throw std::runtime_error("Only 1D numpy arrays are supported!");
    }
    const size_t size = buf.size;
    if (size <= 0) {
        throw std::runtime_error("Numpy array is empty!");
    }
    // get the pointer of buffer
    T* buf_ptr = static_cast<T*>(buf.ptr);
    // resize the target vector if necessary
    if (target.size() != size) {
        target.resize(size);
    }
    // get the raw pointer of the device_vector's data
    T* d_ptr = thrust::raw_pointer_cast(target.data());
    // copy data from host to device
    // remember, the size is counted in bytes
    {
        py::gil_scoped_release release;
        cudaError_t err = cudaMemcpy(
            d_ptr, buf_ptr, size * sizeof(T), cudaMemcpyHostToDevice
        );
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("cudaMemcpy failed in transfering data!") + cudaGetErrorString(err));
        }
    }
}

#endif