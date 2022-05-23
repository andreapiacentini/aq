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

#ifndef AQ_AQ_GEOM_INTERFACE_H_
#define AQ_AQ_GEOM_INTERFACE_H_

#define AQ_STRLEN 256

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "aq/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
extern "C" {
  void aq_geom_setup_f90(F90geom &, const eckit::Configuration &, const eckit::mpi::Comm *,
                         atlas::grid::GridImpl *,
                         atlas::functionspace::FunctionSpaceImpl *);
  void aq_geom_fill_atlas_fieldset_f90(const F90geom &, atlas::field::FieldSetImpl *);
  void aq_geom_clone_f90(F90geom &, const F90geom &);
  void aq_geom_info_f90(const F90geom &, int &, int &, int &,
                        double &, double &, int &, char *, char *, char *);
  void aq_geom_delete_f90(F90geom &);
}
}  // namespace aq
#endif  // AQ_AQ_GEOM_INTERFACE_H_
