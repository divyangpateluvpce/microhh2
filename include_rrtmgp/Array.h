#ifndef ARRAY_H
#define ARRAY_H

#include <array>
#include <vector>
#include <iostream>

template<int N>
inline std::array<int, N> calc_strides(const std::array<int, N>& dims)
{
    std::array<int, N> strides;
    strides[0] = 1;
    for(int i=1; i<N; ++i)
        strides[i] = strides[i-1]*dims[i-1];

    return strides;
}

template<int N>
inline int dot(const std::array<int, N>& left, const std::array<int, N>& right)
{
    int sum = 0;
    for (int i=0; i<N; ++i)
        sum += left[i]*right[i];

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

template<typename T, int N>
struct Array
{
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

    void operator=(std::vector<T>&& data)
    {
        // CvH check size.
        this->data = data;
    }

    inline T& operator()(const std::array<int, N>& indices)
    {
        const int index = dot<N>(indices, strides);
        return data[index];
    }

    inline T operator()(const std::array<int, N>& index) const
    {
        const int i = dot<N>(index, strides);
        return data[i];
    }

    const std::array<int, N> dims;
    const int ncells;
    std::vector<T> data;
    const std::array<int, N> strides;
};
#endif
