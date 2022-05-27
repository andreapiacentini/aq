/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
