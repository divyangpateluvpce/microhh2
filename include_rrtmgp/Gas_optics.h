#ifndef GAS_OPTICS_H
#define GAS_OPTICS_H

#include "Optical_props.h"
#include "Array.h"

template<typename TF>
class Gas_optics
{
    public:
        Gas_optics(
            const std::vector<Gas_concs<TF>>& available_gases,
            Array<std::string,1> gas_names,
            Array<int,3> key_species,
            Array<int,2> band2gpt,
            Array<TF,2> band_lims_wavenum,
            Array<TF,1> press_ref,
            TF press_ref_trop,
            Array<TF,1> temp_ref,
            TF temp_ref_p,
            TF temp_ref_t,
            Array<TF,3> vmr_ref,
            Array<TF,4> kmajor,
            Array<TF,3> kminor_lower,
            Array<TF,3> kminor_upper,
            Array<std::string,1> gas_minor,
            Array<std::string,1> identifier_minor,
            Array<std::string,1> minor_gases_lower,
            Array<std::string,1> minor_gases_upper,
            Array<int,2> minor_limits_gpt_lower,
            Array<int,2> minor_limits_gpt_upper,
            Array<int,1> minor_scales_with_density_lower,
            Array<int,1> minor_scales_with_density_upper,
            Array<std::string,1> scaling_gas_lower,
            Array<std::string,1> scaling_gas_upper,
            Array<int,1> scale_by_complement_lower,
            Array<int,1> scale_by_complement_upper,
            Array<int,1> kminor_start_lower,
            Array<int,1> kminor_start_upper,
            Array<TF,2> totplnk,
            Array<TF,4> planck_frac,
            Array<TF,3> rayl_lower,
            Array<TF,3> rayl_upper)
        {}
};
#endif
