#ifndef GAS_CONCS_H
#define GAS_CONCS_H

#include <vector>

template<typename TF>
class Gas_concs
{
    public:
        Gas_concs(const std::string& name, std::vector<TF> values, const int nlay, const int ncol) :
                name(name), ncol(ncol), nlay(nlay), w(values) {}
        Gas_concs(const std::string& name, const TF value) :
                name(name), nlay(1), ncol(1), w(1, value) {}
        void print_w() const
        {
            std::cout << name << std::endl;
            for (TF v : w)
                std::cout << v << std::endl;
        }

    private:
        std::string name;
        const int ncol;
        const int nlay;
        std::vector<TF> w; 
};
#endif
