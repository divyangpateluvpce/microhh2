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

#ifndef RADIATION_H
#define RADIATION_H

#include <vector>
#include "field3d_operators.h"

enum class Radiation_type {Enabled, Disabled};

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Thermo;
class Netcdf_handle;

template<typename TF>
class Radiation
{
    public:
        Radiation(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Radiation();
        void init();
        void create(Thermo<TF>&, Netcdf_handle&);
        void exec(Thermo<TF>&);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_operators<TF> field3d_operators;

        Radiation_type swradiation;

        std::vector<TF> pres_layer;
        std::vector<TF> temp_layer;
        std::vector<TF> pres_level;
        std::vector<TF> temp_level;
};
#endif
