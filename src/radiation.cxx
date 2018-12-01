/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <string>
#include <iostream>
#include <boost/algorithm/string.hpp>

#include "radiation.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "netcdf_interface.h"

#include "Gas_concs.h"

namespace
{
}

template<typename TF>
Radiation<TF>::Radiation(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields)
{
    // Read the switches from the input
    std::string swradiation_in = inputin.get_item<std::string>("radiation", "swradiation", "", "0");

    if (swradiation_in == "0")
        swradiation = Radiation_type::Disabled;
    else if (swradiation_in == "1")
        swradiation = Radiation_type::Enabled;
    else
        throw std::runtime_error("Invalid option for \"swradiation\"");
}

template<typename TF>
Radiation<TF>::~Radiation()
{
}

template<typename TF>
void Radiation<TF>::init()
{
    if (swradiation == Radiation_type::Disabled)
        return;
}

template<typename TF>
void Radiation<TF>::create(Thermo<TF>& thermo, Netcdf_handle& input_nc)
{
    if (swradiation == Radiation_type::Disabled)
        return;

    Netcdf_group group_nc = input_nc.get_group("radiation");

    Netcdf_file coef_lw_nc(master, "coefficients_lw.nc", Netcdf_mode::Read);

    int layer = group_nc.get_variable_dimensions("pres_layer").at("layer");
    int level = group_nc.get_variable_dimensions("pres_level").at("level");

    // Download pressure and temperature data.
    pres_layer.resize(layer);
    pres_level.resize(level);
    temp_layer.resize(layer);
    temp_level.resize(level);

    group_nc.get_variable(pres_layer, "pres_layer", {0}, {layer});
    group_nc.get_variable(pres_level, "pres_level", {0}, {level});
    group_nc.get_variable(temp_layer, "temp_layer", {0}, {layer});
    group_nc.get_variable(temp_level, "temp_level", {0}, {level});

    // Read the gas concentrations.
    std::vector<Gas_concs<TF>> gas_conc_array;

    // Download surface boundary conditions for long wave.
    surface_emissivity.resize(1);
    surface_temperature.resize(1);

    group_nc.get_variable(surface_emissivity, "surface_emissivity", {0}, {1});
    group_nc.get_variable(surface_temperature, "surface_temperature", {0}, {1});

    // READ K-DISTRIBUTION MOVE TO SEPARATE FUNCTION LATER...
    // Read k-distribution information.
    int n_temps          = coef_lw_nc.get_dimension_size("temperature");
    int n_press          = coef_lw_nc.get_dimension_size("pressure");
    int n_absorbers      = coef_lw_nc.get_dimension_size("absorber");
    int n_char           = coef_lw_nc.get_dimension_size("string_len");
    int n_minorabsorbers = coef_lw_nc.get_dimension_size("minor_absorber");
    int n_extabsorbers   = coef_lw_nc.get_dimension_size("absorber_ext");
    int n_mixingfracs    = coef_lw_nc.get_dimension_size("mixing_fraction");
    int n_layers         = coef_lw_nc.get_dimension_size("atmos_layer");
    int n_bnds           = coef_lw_nc.get_dimension_size("bnd");
    int n_gpts           = coef_lw_nc.get_dimension_size("gpt");
    int n_pairs          = coef_lw_nc.get_dimension_size("pair");
    int n_minor_absorber_intervals_lower = coef_lw_nc.get_dimension_size("minor_absorber_intervals_lower");
    int n_minor_absorber_intervals_upper = coef_lw_nc.get_dimension_size("minor_absorber_intervals_upper");
    int n_internal_sourcetemps = coef_lw_nc.get_dimension_size("temperature_Planck");
    int n_contributors_lower = coef_lw_nc.get_dimension_size("contributors_lower");
    int n_contributors_upper = coef_lw_nc.get_dimension_size("contributors_upper");

    // Read gas names.
    std::vector<std::string> gas_names;
    for (int n=0; n<n_absorbers; ++n)
    {
        std::vector<char> gas_name_char(n_char);
        coef_lw_nc.get_variable(gas_name_char, "gas_names", {n,0}, {1,n_char});
        std::string gas_name(gas_name_char.begin(), gas_name_char.end());
        boost::trim(gas_name);
        gas_names.push_back(gas_name);
    }

    std::vector<int> key_species;
    // coef_lw_nc.get_variable(key_species, "key_species", {2,n_layers,n_bnds});
    coef_lw_nc.get_variable(key_species, "key_species", {n_bnds,n_layers,2});

    std::vector<TF> band_lims;
    coef_lw_nc.get_variable(band_lims, "bnd_limits_wavenumber", {n_bnds,2});

    std::vector<int> band2gpt;
    coef_lw_nc.get_variable(band2gpt, "bnd_limits_gpt", {n_bnds,2});

    std::vector<TF> press_ref;
    coef_lw_nc.get_variable(press_ref, "press_ref", {n_press});

    std::vector<TF> temp_ref;
    coef_lw_nc.get_variable(temp_ref, "temp_ref", {n_temps});

    /*
    temp_ref_p        = read_field(ncid, 'absorption_coefficient_ref_P')
    temp_ref_t        = read_field(ncid, 'absorption_coefficient_ref_T')
    press_ref_trop    = read_field(ncid, 'press_ref_trop')
    kminor_lower      = read_field(ncid, 'kminor_lower', &
        ncontributors_lower, nmixingfracs, ntemps)
    kminor_upper      = read_field(ncid, 'kminor_upper', &
        ncontributors_upper, nmixingfracs, ntemps)
    gas_minor = read_char_vec(ncid, 'gas_minor', nminorabsorbers)
    identifier_minor = read_char_vec(ncid, 'identifier_minor', nminorabsorbers)
    minor_gases_lower = read_char_vec(ncid, 'minor_gases_lower', nminor_absorber_intervals_lower)
    minor_gases_upper = read_char_vec(ncid, 'minor_gases_upper', nminor_absorber_intervals_upper)
    minor_limits_gpt_lower &
                      = int(read_field(ncid, 'minor_limits_gpt_lower', npairs,nminor_absorber_intervals_lower))
    minor_limits_gpt_upper &
                      = int(read_field(ncid, 'minor_limits_gpt_upper', npairs,nminor_absorber_intervals_upper))
    minor_scales_with_density_lower &
                      = read_logical_vec(ncid, 'minor_scales_with_density_lower', nminor_absorber_intervals_lower)
    minor_scales_with_density_upper &
                      = read_logical_vec(ncid, 'minor_scales_with_density_upper', nminor_absorber_intervals_upper)
    scale_by_complement_lower &
                      = read_logical_vec(ncid, 'scale_by_complement_lower', nminor_absorber_intervals_lower)
    scale_by_complement_upper &
                      = read_logical_vec(ncid, 'scale_by_complement_upper', nminor_absorber_intervals_upper)
    scaling_gas_lower &
                      = read_char_vec(ncid, 'scaling_gas_lower', nminor_absorber_intervals_lower)
    scaling_gas_upper &
                      = read_char_vec(ncid, 'scaling_gas_upper', nminor_absorber_intervals_upper)
    kminor_start_lower &
                      = read_field(ncid, 'kminor_start_lower', nminor_absorber_intervals_lower)
    kminor_start_upper &
                      = read_field(ncid, 'kminor_start_upper', nminor_absorber_intervals_upper)
    vmr_ref           = read_field(ncid, 'vmr_ref', nlayers, nextabsorbers, ntemps)

    kmajor            = read_field(ncid, 'kmajor',  ngpts, nmixingfracs,  npress+1, ntemps)
    if(var_exists(ncid, 'rayl_lower')) then
      rayl_lower = read_field(ncid, 'rayl_lower',   ngpts, nmixingfracs,            ntemps)
      rayl_upper = read_field(ncid, 'rayl_upper',   ngpts, nmixingfracs,            ntemps)
    end if
    */
    // END READ K-DISTRIBUTION

    throw 666;
}

template<typename TF>
void Radiation<TF>::exec(Thermo<TF>& thermo)
{
    if (swradiation == Radiation_type::Disabled)
        return;
}

template class Radiation<double>;
template class Radiation<float>;
