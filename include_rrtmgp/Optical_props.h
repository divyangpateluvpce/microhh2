#ifndef OPTICAL_PROPS_H
#define OPTICAL_PROPS_H

#include "Array.h"

class Optical_props
{
    public:
        Optical_props() {};
        virtual ~Optical_props() {};
        virtual int get_ncol() const = 0;
        virtual int get_nlay() const = 0;
        virtual int get_ngpt() const = 0;
        virtual std::string get_name() const = 0;

    private:
        // increment_bygpt
        // increment_byband
};

template<typename TF>
class Optical_props_arry : public Optical_props
{
    public:
        Optical_props_arry() {};
        virtual ~Optical_props_arry() {};
        virtual Array<TF,3>& get_tau() = 0;
        virtual Array<TF,3>& get_ssa() = 0;
        virtual Array<TF,3>& get_g() = 0;

        // virtual void get_subset(const int, const int, std::unique_ptr<Optical_props_arry>&) = 0;
};

/*
template<typename TF>
class Optical_props_2str : public Optical_props
{
    public:
        Optical_props_2str(
                const int ncol, const int nlay, const int ngpt,
                Array_3d<TF>&& tau, Array_3d<TF>&& ssa, Array_3d<TF>&& g,
                const std::string name="") :
            tau_(std::move(tau)), ssa_(std::move(ssa)), g_(std::move(g)), name_(name)
        {};
        ~Optical_props_2str() {};

        void get_subset(const int col_start, const int block_size, std::unique_ptr<Optical_props>& optical_props)
        {
        }

    private:
        Array_3d<TF> tau_;
        Array_3d<TF> ssa_;
        Array_3d<TF> g_;
        const std::string name_;
};
*/

template<typename TF>
class Optical_props_1scl : public Optical_props_arry<TF>
{
    public:
        // Initializer constructor.
        Optical_props_1scl(
                const int ncol, const int nlay, const int ngpt,
                Array<TF,3>&& tau,
                const std::string name="") :
            tau_(std::move(tau)),
            name_(name)
        {};

        // Subset constructor.
        Optical_props_1scl(
                std::unique_ptr<Optical_props_arry<TF>>& full,
                const int icol_start,
                const int ncol) :
            name_(full->get_name())
        {
            const int ncol_full = full->get_ncol();
            const int nlay = full->get_nlay();
            const int ngpt = full->get_ngpt();

            if ( (icol_start < 0) ||
                 ( (icol_start + ncol) > ncol_full) )
            {
                throw std::runtime_error("Optical_props::get_subset: out of range");
            }

            // Implement the class switching. Design is a little problematic at the moment.
            tau_.resize(ncol, nlay, ngpt);

            // Assign the sub-column.
            for (int igpt=0; igpt<ngpt; ++igpt)
                for (int ilay=0; ilay<nlay; ++ilay)
                    for (int icol=0; icol<ncol; ++icol)
                        tau_(icol, ilay, igpt) = full->get_tau()(icol+icol_start, ilay, igpt);
        }

        ~Optical_props_1scl() {};

        int get_ncol() const { return tau_.dim1(); }
        int get_nlay() const { return tau_.dim2(); }
        int get_ngpt() const { return tau_.dim3(); }
        std::string get_name() const { return name_; }

        Array<TF,3>& get_tau() { return tau_; }
        Array<TF,3>& get_ssa() { throw std::runtime_error("Not available in this class"); }
        Array<TF,3>& get_g  () { throw std::runtime_error("Not available in this class"); }

    private:
        Array<TF,3> tau_;
        const std::string name_;
};
#endif
