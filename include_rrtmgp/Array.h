#ifndef ARRAY_H
#define ARRAY_H

#include <vector>
#include <algorithm>

template<typename TF>
class Array_1d
{
    public:
        Array_1d() : itot_(0) {};

        Array_1d(const int itot) :
            itot_(itot),
            data_(itot) {}

        Array_1d(std::vector<TF>&& data, const int itot) :
            data_(data),
            itot_(itot) {}

        // Array_1d(Array_1d&&) = default;
        // Array_1d& operator=(Array_1d&&) = default;

        TF& operator()(const int i)
        {
            return data_[i];
        }

        const TF operator()(const int i) const
        {
            return data_[i];
        }

        TF& operator[](const int i) { return data_[i]; }
        const TF operator[](const int i) const { return data_[i]; }

        TF dim1() const { return itot_; }

        void set_zero()
        {
            std::fill(data_.begin(), data_.end(), 0.);
        }

        void resize(const int itot_new)
        { 
            data_.resize(itot_new);
            itot_ = itot_new;
        }

        TF max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        }

        TF min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        }

        TF* get_ptr()
        {
            return data_.data();
        }

    private:
        int itot_;
        std::vector<TF> data_;
};

template<typename TF>
class Array_2d
{
    public:
        Array_2d() : itot_(0), jtot_(0) {};

        Array_2d(const int itot, const int jtot) :
            itot_(itot),
            jtot_(jtot),
            data_(itot*jtot) {}

        Array_2d(std::vector<TF>&& data, const int itot, const int jtot) :
            itot_(itot),
            jtot_(jtot),
            data_(data) {}

        // Array_2d(Array_2d&&) = default;
        // Array_2d& operator=(Array_2d&&) = default;

        TF& operator()(const int j, const int i)
        {
            const int ij = i + j*itot_;
            return data_[ij];
        }

        const TF operator()(const int j, const int i) const
        {
            const int ij = i + j*itot_;
            return data_[ij];
        }

        TF& operator[](const int i) { return data_[i]; }
        const TF operator[](const int i) const { return data_[i]; }

        TF dim0() const { return jtot_; }
        TF dim1() const { return itot_; }

        void set_zero()
        {
            std::fill(data_.begin(), data_.end(), 0.);
        }

        TF max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        }

        TF min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        }

        TF* get_ptr()
        {
            return data_.data();
        }

    private:
        int itot_;
        int jtot_;
        std::vector<TF> data_;
};

template<typename TF>
class Array_3d
{
    public:
        Array_3d() : itot_(0), jtot_(0), ktot_(0) {};

        Array_3d(const int ktot, const int jtot, const int itot) :
            itot_(itot),
            jtot_(jtot),
            ktot_(ktot),
            data_(itot*jtot*ktot) {}

        // Array_3d(Array_3d&&) = default;
        // Array_3d& operator=(Array_3d&&) = default;

        TF& operator()(const int k, const int j, const int i)
        {
            const int ij = i + j*itot_ + k*itot_*jtot_;
            return data_[ij];
        }

        const TF operator()(const int k, const int j, const int i) const
        {
            const int ij = i + j*itot_ + k*itot_*jtot_;
            return data_[ij];
        }

        TF& operator[](const int i) { return data_[i]; }
        const TF operator[](const int i) const { return data_[i]; }

        TF dim0() const { return ktot_; }
        TF dim1() const { return jtot_; }
        TF dim2() const { return itot_; }

        void set_zero()
        {
            std::fill(data_.begin(), data_.end(), 0.);
        }

        TF max() const
        {
            return *std::max_element(std::begin(data_), std::end(data_));
        }

        TF min() const
        {
            return *std::min_element(std::begin(data_), std::end(data_));
        }

        TF* get_ptr()
        {
            return data_.data();
        }

    private:
        int itot_;
        int jtot_;
        int ktot_;
        std::vector<TF> data_;
};

template<typename TF_a, typename TF_b>
bool has_same_dims(const Array_1d<TF_a>& a, const Array_1d<TF_b>& b)
{
    return ( a.dim0() == b.dim0() );
}

template<typename TF_a, typename TF_b>
bool has_same_dims(const Array_2d<TF_a>& a, const Array_2d<TF_b>& b)
{
    return ( ( a.dim0() == b.dim0() ) &&
             ( a.dim1() == b.dim1() ) );
}

template<typename TF_a, typename TF_b>
bool has_same_dims(const Array_3d<TF_a>& a, const Array_3d<TF_b>& b)
{
    return ( ( a.dim0() == b.dim0() ) &&
             ( a.dim1() == b.dim1() ) &&
             ( a.dim2() == b.dim2() ) );
}
#endif
