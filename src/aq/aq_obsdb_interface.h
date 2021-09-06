/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
                          const eckit::mpi::Comm *);
  void aq_obsdb_delete_f90(F90odb &);
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
