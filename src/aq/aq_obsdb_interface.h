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

#ifndef AQ_AQ_OBSDB_INTERFACE_H_
#define AQ_AQ_OBSDB_INTERFACE_H_

#include <memory>
#include <string>
#include <vector>

#include "aq/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace aq {
extern "C" {
  void aq_obsdb_setup_f90(F90odb &, const eckit::Configuration &,
                          const util::DateTime &, const util::DateTime &,
                          const bool &, F90odb &);
  void aq_obsdb_delete_f90(F90odb &, const bool &);
  void aq_obsdb_read_f90(F90odb &);
  void aq_obsdb_get_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void aq_obsdb_put_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void aq_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              atlas::field::FieldSetImpl *, std::vector<util::DateTime> &);
  void aq_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration &, const util::DateTime &,
                             const util::Duration &, const int &, int &);
  void aq_obsdb_nobs_f90(const F90odb &, const int &, const char *, int &);
}

}  // namespace aq
#endif  // AQ_AQ_OBSDB_INTERFACE_H_
