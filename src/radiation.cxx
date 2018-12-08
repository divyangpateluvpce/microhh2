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

#include "Array.h"
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
    std::vector<std::string> get_variable_string(
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
        var_char = input_nc.get_variable<char>(var_name, i_count);

        std::vector<std::string> var;

        for (int n=0; n<total_count; ++n)
        {
            std::string s(var_char.begin()+n*string_len, var_char.begin()+(n+1)*string_len);
            if (trim)
                boost::trim(s);
            var.push_back(s);
        }

        return var;
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

    Array<double,2> pres_layer(group_nc.get_variable<double>("pres_layer", {n_lay, n_col}), {n_lay, n_col});
    Array<double,2> pres_level(group_nc.get_variable<double>("pres_level", {n_lev, n_col}), {n_lay, n_col});
    Array<double,2> temp_layer(group_nc.get_variable<double>("temp_layer", {n_lay, n_col}), {n_lay, n_col});
    Array<double,2> temp_level(group_nc.get_variable<double>("temp_level", {n_lev, n_col}), {n_lay, n_col});

    const int top_at_1 = pres_layer({0, 0}) < pres_layer({n_lay-1, 0});

    // Download surface boundary conditions for long wave.
    Array<double,1> surface_emissivity (group_nc.get_variable<double>("surface_emissivity" , {n_col}), {n_col});
    Array<double,1> surface_temperature(group_nc.get_variable<double>("surface_temperature", {n_col}), {n_col});

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
    Array<std::string,1> gas_names(get_variable_string("gas_names", {n_absorbers}, coef_lw_nc, n_char, true), {n_absorbers});

    Array<int,3> key_species(coef_lw_nc.get_variable<int>("key_species", {n_bnds, n_layers, 2}), {n_bnds, n_layers, 2});
    Array<double,2> band_lims(coef_lw_nc.get_variable<double>("bnd_limits_wavenumber", {n_bnds, 2}), {n_bnds, 2});
    Array<int,2> band2gpt(coef_lw_nc.get_variable<int>("bnd_limits_gpt", {n_bnds, 2}), {n_bnds, 2});
    Array<double,1> press_ref(coef_lw_nc.get_variable<double>("press_ref", {n_press}), {n_press});
    Array<double,1> temp_ref(coef_lw_nc.get_variable<double>("temp_ref", {n_temps}), {n_temps});

    double temp_ref_p = coef_lw_nc.get_variable<double>("absorption_coefficient_ref_P");
    double temp_ref_t = coef_lw_nc.get_variable<double>("absorption_coefficient_ref_T");
    double press_ref_trop = coef_lw_nc.get_variable<double>("press_ref_trop");

    Array<double,3> kminor_lower(coef_lw_nc.get_variable<double>("kminor_lower", {n_temps, n_mixingfracs, n_contributors_lower}), {n_temps, n_mixingfracs, n_contributors_lower});
    Array<double,3> kminor_upper(coef_lw_nc.get_variable<double>("kminor_upper", {n_temps, n_mixingfracs, n_contributors_upper}), {n_temps, n_mixingfracs, n_contributors_lower});

    Array<std::string,1> gas_minor(get_variable_string("gas_minor", {n_minorabsorbers}, coef_lw_nc, n_char, false), {n_minorabsorbers});
    Array<std::string,1> identifier_minor(get_variable_string("identifier_minor", {n_minorabsorbers}, coef_lw_nc, n_char, false), {n_minorabsorbers});

    Array<std::string,1> minor_gases_lower(get_variable_string("minor_gases_lower", {n_minor_absorber_intervals_lower}, coef_lw_nc, n_char, false), {n_minor_absorber_intervals_lower});
    Array<std::string,1> minor_gases_upper(get_variable_string("minor_gases_upper", {n_minor_absorber_intervals_upper}, coef_lw_nc, n_char, false), {n_minor_absorber_intervals_upper});

    Array<int,2> minor_limits_gpt_lower(coef_lw_nc.get_variable<int>("minor_limits_gpt_lower", {n_minor_absorber_intervals_lower, n_pairs}), {n_minor_absorber_intervals_lower, n_pairs});
    Array<int,2> minor_limits_gpt_upper(coef_lw_nc.get_variable<int>("minor_limits_gpt_upper", {n_minor_absorber_intervals_upper, n_pairs}), {n_minor_absorber_intervals_upper, n_pairs});

    Array<int,1> minor_scales_with_density_lower(coef_lw_nc.get_variable<int>("minor_scales_with_density_lower", {n_minor_absorber_intervals_lower}), {n_minor_absorber_intervals_lower});
    Array<int,1> minor_scales_with_density_upper(coef_lw_nc.get_variable<int>("minor_scales_with_density_upper", {n_minor_absorber_intervals_upper}), {n_minor_absorber_intervals_upper});

    Array<int,1> scale_by_complement_lower(coef_lw_nc.get_variable<int>("scale_by_complement_lower", {n_minor_absorber_intervals_lower}), {n_minor_absorber_intervals_lower});
    Array<int,1> scale_by_complement_upper(coef_lw_nc.get_variable<int>("scale_by_complement_upper", {n_minor_absorber_intervals_upper}), {n_minor_absorber_intervals_upper});

    Array<std::string,1> scaling_gas_lower(get_variable_string("scaling_gas_lower", {n_minor_absorber_intervals_lower}, coef_lw_nc, n_char, false), {n_minor_absorber_intervals_lower});
    Array<std::string,1> scaling_gas_upper(get_variable_string("scaling_gas_upper", {n_minor_absorber_intervals_upper}, coef_lw_nc, n_char, false), {n_minor_absorber_intervals_upper});

    Array<int,1> kminor_start_lower(coef_lw_nc.get_variable<int>("kminor_start_lower", {n_minor_absorber_intervals_lower}), {n_minor_absorber_intervals_lower});
    Array<int,1> kminor_start_upper(coef_lw_nc.get_variable<int>("kminor_start_upper", {n_minor_absorber_intervals_upper}), {n_minor_absorber_intervals_upper});

    Array<double,3> vmr_ref(coef_lw_nc.get_variable<double>("vmr_ref", {n_temps, n_extabsorbers, n_layers}), {n_temps, n_extabsorbers, n_layers});

    Array<double,4> kmajor(coef_lw_nc.get_variable<double>("kmajor", {n_temps, n_press+1, n_mixingfracs, n_gpts}), {n_temps, n_press+1, n_mixingfracs, n_gpts});

    Array<double,3> rayl_lower({n_temps, n_mixingfracs, n_gpts});
    Array<double,3> rayl_upper({n_temps, n_mixingfracs, n_gpts});

    if (coef_lw_nc.variable_exists("rayl_lower"))
    {
        // rayl_lower = read_field(ncid, 'rayl_lower',   ngpts, nmixingfracs,            ntemps)
        // rayl_upper = read_field(ncid, 'rayl_upper',   ngpts, nmixingfracs,            ntemps)
    }

    // Is it really LW if so read these variables as well.
    Array<double,2> totplnk({n_bnds, n_internal_sourcetemps});
    Array<double,4> planck_frac({n_temps, n_press+1, n_mixingfracs, n_gpts});

    if (coef_lw_nc.variable_exists("totplnk"))
    {
        totplnk = coef_lw_nc.get_variable<double>("totplnk", {n_bnds, n_internal_sourcetemps});
        planck_frac = coef_lw_nc.get_variable<double>("plank_fraction", {n_temps, n_press+1, n_mixingfracs, n_gpts});
    }
    // END READ K-DISTRIBUTION

    // Read the gas concentrations.
    std::vector<Gas_concs<double>> available_gases;

    for (int i=0; i<gas_names.get_dims()[0]; ++i)
    {
        const std::string& gas_name = gas_names[{i}];
        if (gas_name == "h2o" || gas_name == "o3")
        {
            Array<double,2> conc(group_nc.get_variable<double>(gas_name, {n_lay, n_col}), {n_lay, n_col});
            available_gases.emplace_back(gas_name, conc.v(), n_lay, n_col);
        }
        else
        {
            double conc = group_nc.get_variable<double>(gas_name);
            available_gases.emplace_back(gas_name, conc);
        }
    }

    for (auto g : gas_names)
        std::cout << "CvH: " << g.first << " = " << g.second << std::endl;

    // Construct the k-distribution.
    Gas_optics<TF> kdist(
            available_gases,
            gas_names,
            key_species,
            band2gpt,
            band_lims,
            press_ref,
            press_ref_trop,
            temp_ref,
            temp_ref_p,
            temp_ref_t,
            vmr_ref,
            kmajor,
            kminor_lower,
            kminor_upper,
            gas_minor,
            identifier_minor,
            minor_gases_lower,
            minor_gases_upper,
            minor_limits_gpt_lower,
            minor_limits_gpt_upper,
            minor_scales_with_density_lower,
            minor_scales_with_density_upper,
            scaling_gas_lower,
            scaling_gas_upper,
            scale_by_complement_lower,
            scale_by_complement_upper,
            kminor_start_lower,
            kminor_start_upper,
            totplnk,
            planck_frac,
            rayl_lower,
            rayl_upper);

    for (auto g : gas_names)
        std::cout << "CvH: " << g.first << " = " << g.second << std::endl;

    throw 666;
}

template<typename TF>
void Radiation<TF>::exec(Thermo<TF>& thermo)
{
    if (swradiation == Radiation_type::Disabled)
        return;
}

template class Radiation<double>;
// template class Radiation<float>;
