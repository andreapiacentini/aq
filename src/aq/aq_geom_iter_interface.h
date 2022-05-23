/*
 * This file is part of the Air Quality Ensemble Data Assimilation suite AQ.
 *
 * (C) Copyright 2022 CERFACS
 *
 * AQ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * AQ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * A copy of the GNU Lesser General Public License is distributed
 * along with AQ (files LICENSE.md, COPYING and COPYING.LESSER).
 */

#ifndef AQ_AQ_GEOM_ITER_INTERFACE_H_
#define AQ_AQ_GEOM_ITER_INTERFACE_H_

#include <memory>
#include <string>
#include <vector>

#include "aq/interface.h"

namespace aq {
extern "C" {
  void aq_geom_iter_setup_f90(F90iter &, const F90geom &, const int &);
  void aq_geom_iter_clone_f90(F90iter &, const F90iter &);
  void aq_geom_iter_delete_f90(F90iter &);
  void aq_geom_iter_equals_f90(const F90iter &, const F90iter&, int &);
  void aq_geom_iter_current_f90(const F90iter &, double &, double &);
  void aq_geom_iter_next_f90(const F90iter &);
}
}  // namespace aq
#endif  // AQ_AQ_GEOM_ITER_INTERFACE_H_
