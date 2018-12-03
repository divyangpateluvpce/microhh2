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
#include <numeric>
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
#include "Gas_optics.h"

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

namespace
{
    void get_variable_string(
            std::vector<std::string>& var,
            const std::string& var_name,
            std::vector<int> i_count,
            Netcdf_handle& input_nc,
            const int string_len,
            bool trim=false)
    {
        // Multiply all elements in i_count.
        int total_count = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());

        // Add the string length as the rightmost dimension.
        i_count.push_back(string_len);

        // Multiply all elements in i_count.
        int total_count_char = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());

        // Read the entire char array;
        std::vector<char> var_char;
        input_nc.get_variable(var_char, var_name, i_count);

        for (int n=0; n<total_count; ++n)
        {
            std::string s(var_char.begin()+n*string_len, var_char.begin()+(n+1)*string_len);
            if (trim)
                boost::trim(s);
            var.push_back(s);
        }
    }
}

template<typename TF>
void Radiation<TF>::create(Thermo<TF>& thermo, Netcdf_handle& input_nc)
{
    if (swradiation == Radiation_type::Disabled)
        return;

    Netcdf_group group_nc = input_nc.get_group("radiation");

    Netcdf_file coef_lw_nc(master, "coefficients_lw.nc", Netcdf_mode::Read);

    int n_lay = group_nc.get_variable_dimensions("pres_layer").at("layer");
    int n_lev = group_nc.get_variable_dimensions("pres_level").at("level");
    int n_col = 1;

    // Download pressure and temperature data.
    pres_layer.resize(n_lay);
    pres_level.resize(n_lev);
    temp_layer.resize(n_lay);
    temp_level.resize(n_lev);

    // CvH This is only functional for 1D tests.... ncol needs to be included.
    group_nc.get_variable(pres_layer, "pres_layer", {0}, {n_lay});
    group_nc.get_variable(pres_level, "pres_level", {0}, {n_lev});
    group_nc.get_variable(temp_layer, "temp_layer", {0}, {n_lay});
    group_nc.get_variable(temp_level, "temp_level", {0}, {n_lev});

    int top_at_1 = pres_layer[0] < pres_layer[n_lay-1];

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
    get_variable_string(gas_names, "gas_names", {n_absorbers}, coef_lw_nc, n_char, true);

    std::vector<int> key_species;
    coef_lw_nc.get_variable(key_species, "key_species", {n_bnds,n_layers, 2});

    std::vector<double> band_lims;
    coef_lw_nc.get_variable(band_lims, "bnd_limits_wavenumber", {n_bnds, 2});

    std::vector<int> band2gpt;
    coef_lw_nc.get_variable(band2gpt, "bnd_limits_gpt", {n_bnds, 2});

    std::vector<double> press_ref;
    coef_lw_nc.get_variable(press_ref, "press_ref", {n_press});

    std::vector<double> temp_ref;
    coef_lw_nc.get_variable(temp_ref, "temp_ref", {n_temps});

    double temp_ref_p;
    coef_lw_nc.get_variable(temp_ref_p, "absorption_coefficient_ref_P");

    double temp_ref_t;
    coef_lw_nc.get_variable(temp_ref_t, "absorption_coefficient_ref_T");

    double press_ref_trop;
    coef_lw_nc.get_variable(press_ref_trop, "press_ref_trop");

    std::vector<double> kminor_lower, kminor_upper;
    coef_lw_nc.get_variable(kminor_lower, "kminor_lower", {n_temps, n_mixingfracs, n_contributors_lower});
    coef_lw_nc.get_variable(kminor_upper, "kminor_upper", {n_temps, n_mixingfracs, n_contributors_upper});

    std::vector<std::string> gas_minor, identifier_minor;
    get_variable_string(gas_minor, "gas_minor", {n_minorabsorbers}, coef_lw_nc, n_char, false);
    get_variable_string(identifier_minor, "identifier_minor", {n_minorabsorbers}, coef_lw_nc, n_char, false);

    std::vector<std::string> minor_gases_lower, minor_gases_upper;
    get_variable_string(minor_gases_lower, "minor_gases_lower", {n_minor_absorber_intervals_lower}, coef_lw_nc, n_char, false);
    get_variable_string(minor_gases_upper, "minor_gases_upper", {n_minor_absorber_intervals_upper}, coef_lw_nc, n_char, false);

    std::vector<int> minor_limits_gpt_lower, minor_limits_gpt_upper;
    coef_lw_nc.get_variable(minor_limits_gpt_lower, "minor_limits_gpt_lower", {n_minor_absorber_intervals_lower, n_pairs});
    coef_lw_nc.get_variable(minor_limits_gpt_upper, "minor_limits_gpt_upper", {n_minor_absorber_intervals_upper, n_pairs});

    std::vector<int> minor_scales_with_density_lower, minor_scales_with_density_upper;
    coef_lw_nc.get_variable(minor_scales_with_density_lower, "minor_scales_with_density_lower", {n_minor_absorber_intervals_lower});
    coef_lw_nc.get_variable(minor_scales_with_density_upper, "minor_scales_with_density_upper", {n_minor_absorber_intervals_upper});

    std::vector<int> scale_by_complement_lower, scale_by_complement_upper;
    coef_lw_nc.get_variable(scale_by_complement_lower, "scale_by_complement_lower", {n_minor_absorber_intervals_lower});
    coef_lw_nc.get_variable(scale_by_complement_upper, "scale_by_complement_upper", {n_minor_absorber_intervals_upper});

    std::vector<std::string> scaling_gas_lower, scaling_gas_upper;
    get_variable_string(scaling_gas_lower, "scaling_gas_lower", {n_minor_absorber_intervals_lower}, coef_lw_nc, n_char, false);
    get_variable_string(scaling_gas_upper, "scaling_gas_upper", {n_minor_absorber_intervals_upper}, coef_lw_nc, n_char, false);

    std::vector<int> kminor_start_lower, kminor_start_upper;
    coef_lw_nc.get_variable(kminor_start_lower, "kminor_start_lower", {n_minor_absorber_intervals_lower});
    coef_lw_nc.get_variable(kminor_start_upper, "kminor_start_upper", {n_minor_absorber_intervals_upper});

    std::vector<double> vmr_ref;
    coef_lw_nc.get_variable(vmr_ref, "vmr_ref", {n_temps, n_extabsorbers, n_layers});

    std::vector<double> kmajor;
    coef_lw_nc.get_variable(kmajor, "kmajor", {n_temps, n_press+1, n_mixingfracs, n_gpts});

    std::vector<double> rayl_lower, rayl_upper;
    if (coef_lw_nc.variable_exists("rayl_lower"))
    {
        // rayl_lower = read_field(ncid, 'rayl_lower',   ngpts, nmixingfracs,            ntemps)
        // rayl_upper = read_field(ncid, 'rayl_upper',   ngpts, nmixingfracs,            ntemps)
    }

    // Is it really LW if so read these variables as well.
    std::vector<double> totplnk, planck_frac;
    if (coef_lw_nc.variable_exists("totplnk"))
    {
        coef_lw_nc.get_variable(totplnk, "totplnk", {n_bnds, n_internal_sourcetemps});
        coef_lw_nc.get_variable(planck_frac, "plank_fraction", {n_temps, n_press+1, n_mixingfracs, n_gpts});
    }
    // END READ K-DISTRIBUTION

    // Read the gas concentrations.
    std::vector<Gas_concs<TF>> available_gases;
    for (const std::string& gas_name : gas_names)
    {
        if (gas_name == "h2o" || gas_name == "o3")
        {
            std::vector<TF> conc;
            group_nc.get_variable(conc, gas_name, {n_lay, n_col});
            available_gases.emplace_back(gas_name, conc, n_lay, n_col);
        }
        else
        {
            TF conc;
            group_nc.get_variable(conc, gas_name);
            available_gases.emplace_back(gas_name, conc);
        }
    }

    // Construct the k-distribution.
    Gas_optics<TF> kdist(
            available_gases, gas_names, key_species,
            band2gpt, band_lims,
            press_ref, press_ref_trop, temp_ref,
            temp_ref_p, temp_ref_t, vmr_ref,
            kmajor, kminor_lower, kminor_upper,
            gas_minor, identifier_minor,
            minor_gases_lower, minor_gases_upper,
            minor_limits_gpt_lower, minor_limits_gpt_upper,
            minor_scales_with_density_lower,
            minor_scales_with_density_upper,
            scaling_gas_lower, scaling_gas_upper,
            scale_by_complement_lower,
            scale_by_complement_upper,
            kminor_start_lower,
            kminor_start_upper,
            totplnk, planck_frac, rayl_lower, rayl_upper);

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
