/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <cuda_runtime_api.h>
#include "master.h"
#include "tools.h"
#define ENV_LOCAL_RANK	"MV2_COMM_WORLD_LOCAL_RANK"

void Master::SetDeviceBeforeInit() // It's tested with single node multi-GPU.
{
  #ifdef USEMPI
	int devCount = 0;
  cuda_safe_call(cudaGetDeviceCount(&devCount));
  cuda_safe_call(cudaSetDevice(md.mpiid % devCount));
  cuda_check_error();
  #endif
}
