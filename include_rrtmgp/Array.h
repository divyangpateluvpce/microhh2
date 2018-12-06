#ifndef ARRAY_H
#define ARRAY_H

#include <array>
#include <vector>
#include <iostream>

template<int N>
inline std::array<int, N> calc_strides(const std::array<int, N>& dims)
{
    std::array<int, N> strides;
    strides[N-1] = 1;
    for (int i=N-2; i>=0; --i)
        strides[i] = strides[i+1]*dims[i+1];

    return strides;
}

template<int N>
inline int c_index(const std::array<int, N>& index, const std::array<int, N>& stride)
{
    int sum = 0;
    for (int i=0; i<N; ++i)
        sum += index[i]*stride[i];

    return sum;
}

template<int N>
inline int fortran_index(const std::array<int, N>& index, const std::array<int, N>& stride)
{
    int sum = 0;
    for (int i=0; i<N; ++i)
        sum += (index[i]-1)*stride[N-1-i];

    return sum;
}

template<int N>
inline int product(const std::array<int, N>& array)
{
    int product = array[0];
    for (int i=1; i<N; ++i)
        product *= array[i];

    return product;
}

template<typename T>
class Array_iterator
{
    public:
        Array_iterator(const std::vector<T>& data, const int n) : data(data), n(n) {}
        Array_iterator operator++() { ++n; }
        std::pair<int, T> operator*() const { return std::make_pair(n, data[n]); }

    private:
        const std::vector<T>& data;
        int n;

    friend bool operator!=(const Array_iterator<T>& left, const Array_iterator<T>& right) { return left.n != right.n; }
};


template<typename T, int N>
class Array
{
    public:
        Array(const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(ncells),
            strides(calc_strides<N>(dims))
        {}

        Array(std::vector<T>&& data, const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(data),
            strides(calc_strides<N>(dims))
        {} // CvH Do we need to size check data?

        inline std::vector<T>& v() { return data; }

        inline void operator=(std::vector<T>&& data)
        {
            // CvH check size.
            this->data = data;
        }

        // C++-style indexing with []
        inline T& operator[](const std::array<int, N>& indices)
        {
            const int index = c_index<N>(indices, strides);
            return data[index];
        }

        inline T operator[](const std::array<int, N>& index) const
        {
            const int i = c_index<N>(index, strides);
            return data[i];
        }

        // Fortran-style indexing with ()
        inline T& operator()(const std::array<int, N>& indices)
        {
            const int index = fortran_index<N>(indices, strides);
            return data[index];
        }

        inline T operator()(const std::array<int, N>& index) const
        {
            const int i = fortran_index<N>(index, strides);
            return data[i];
        }

        inline Array_iterator<T> begin() { return Array_iterator<T>(data, 0); }
        inline Array_iterator<T> end()   { return Array_iterator<T>(data, ncells); }

        const std::array<int, N> dims;

    private:
        const int ncells;
        std::vector<T> data;
        const std::array<int, N> strides;
};
#endif
