#ifndef ARRAY_H
#define ARRAY_H

#include <vector>
#include <algorithm>

template<typename T>
class Array_1d
{
    public:
        Array_1d() : itot_(0) {};

        Array_1d(const int itot) :
            itot_(itot),
            data_(itot) {}

        Array_1d(std::vector<T>&& data, const int itot) :
            itot_(itot),
            data_(data) {}

        Array_1d(Array_1d&&) = default;
        Array_1d& operator=(Array_1d&&) = default;

        // Fortran indexing.
        T& f(const int i) { return data_[i-1]; }
        const T f(const int i) const { return data_[i-1]; }

        T dim1() const { return itot_; }

        void set_zero()
        {
            std::fill(data_.begin(), data_.end(), 0.);
        }

        void resize(const int itot_new)
        { 
            data_.resize(itot_new);
            itot_ = itot_new;
        }

        T max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        }

        T min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        }

        T* get_ptr()
        {
            return data_.data();
        }

    private:
        int itot_;
        std::vector<T> data_;
};

template<typename T>
class Array_2d
{
    public:
        Array_2d() : itot_(0), jtot_(0) {};

        // Initialization goes in C order, fastest goes last!
        Array_2d(const int jtot, const int itot) :
            itot_(itot),
            jtot_(jtot),
            data_(itot*jtot) {}

        Array_2d(std::vector<T>&& data, const int jtot, const int itot) :
            itot_(itot),
            jtot_(jtot),
            data_(data) {}

        Array_2d(Array_2d&&) = default;
        Array_2d& operator=(Array_2d&&) = default;

        // C-style indexing.
        T& operator()(const int j, const int i)
        {
            const int ij = i + j*itot_;
            return data_[ij];
        }

        const T operator()(const int j, const int i) const
        {
            const int ij = i + j*itot_;
            return data_[ij];
        }

        // Fortran indexing.
        T& f(const int i, const int j)
        {
            const int ij = (i-1) + (j-1)*itot_;
            return data_[ij];
        }

        const T f(const int i, const int j) const
        {
            const int ij = (i-1) + (j-1)*itot_;
            return data_[ij];
        }

        T dim1() const { return itot_; }
        T dim2() const { return jtot_; }

        void set_zero()
        {
            std::fill(data_.begin(), data_.end(), 0.);
        }

        T max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        }

        T min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        }

        T* get_ptr()
        {
            return data_.data();
        }

    private:
        int itot_;
        int jtot_;
        std::vector<T> data_;
};

template<typename T>
class Array_3d
{
    public:
        Array_3d() : itot_(0), jtot_(0), ktot_(0) {};

        Array_3d(const int ktot, const int jtot, const int itot) :
            itot_(itot),
            jtot_(jtot),
            ktot_(ktot),
            data_(itot*jtot*ktot) {}

        Array_3d(std::vector<T>&& data, const int ktot, const int jtot, const int itot) :
            itot_(itot),
            jtot_(jtot),
            ktot_(ktot),
            data_(data) {}

        Array_3d(Array_3d&&) = default;
        Array_3d& operator=(Array_3d&&) = default;

        // C-indexing.
        T& operator()(const int k, const int j, const int i)
        {
            const int ij = i + j*itot_ + k*itot_*jtot_;
            return data_[ij];
        }

        const T operator()(const int k, const int j, const int i) const
        {
            const int ij = i + j*itot_ + k*itot_*jtot_;
            return data_[ij];
        }

        // Fortran indexing.
        T& f(const int i, const int j, const int k)
        {
            const int ij = (i-1) + (j-1)*itot_ + (k-1)*itot_*jtot_;
            return data_[ij];
        }

        const T f(const int i, const int j, const int k) const
        {
            const int ij = (i-1) + (j-1)*itot_ + (k-1)*itot_*jtot_;
            return data_[ij];
        }

        T dim1() const { return itot_; }
        T dim2() const { return jtot_; }
        T dim3() const { return ktot_; }

        void set_zero()
        {
            std::fill(data_.begin(), data_.end(), 0.);
        }

        T max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        }

        T min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        }

        T* get_ptr()
        {
            return data_.data();
        }

    private:
        int itot_;
        int jtot_;
        int ktot_;
        std::vector<T> data_;
};

template<typename T_a, typename T_b>
bool has_same_dims(const Array_1d<T_a>& a, const Array_1d<T_b>& b)
{
    return ( a.dim1() == b.dim1() );
}

template<typename T_a, typename T_b>
bool has_same_dims(const Array_2d<T_a>& a, const Array_2d<T_b>& b)
{
    return ( ( a.dim1() == b.dim1() ) &&
             ( a.dim2() == b.dim2() ) );
}

template<typename T_a, typename T_b>
bool has_same_dims(const Array_3d<T_a>& a, const Array_3d<T_b>& b)
{
    return ( ( a.dim1() == b.dim1() ) &&
             ( a.dim2() == b.dim2() ) &&
             ( a.dim3() == b.dim3() ) );
}
#endif
