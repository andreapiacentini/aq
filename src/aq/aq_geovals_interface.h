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

#ifndef AQ_AQ_GEOVALS_INTERFACE_H_
#define AQ_AQ_GEOVALS_INTERFACE_H_

#include <memory>
#include <string>
#include <vector>

#include "aq/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace aq {
  class Locations;

extern "C" {
  void aq_geovals_setup_f90(F90geovals &, const int &, const oops::Variables &,
                            const eckit::mpi::Comm *);
  void aq_geovals_create_f90(F90geovals &, const oops::Variables &, const eckit::mpi::Comm *);
  void aq_geovals_delete_f90(F90geovals &);
  void aq_geovals_copy_f90(const F90geovals &, const F90geovals &);
  void aq_geovals_fill_f90(const F90geovals &, const int &, const char *, const int &,
                           const int *, const int &, const double *);
  void aq_geovals_fillad_f90(const F90geovals &, const int &, const char *, const int &,
                             const int *, const int &, double *);
  void aq_geovals_zero_f90(const F90geovals &);
  void aq_geovals_abs_f90(const F90geovals &);
  void aq_geovals_random_f90(const F90geovals &);
  void aq_geovals_mult_f90(const F90geovals &, const double &);
  void aq_geovals_add_f90(const F90geovals &, const F90geovals &);
  void aq_geovals_diff_f90(const F90geovals &, const F90geovals &);
  void aq_geovals_schurmult_f90(const F90geovals &, const F90geovals &);
  void aq_geovals_divide_f90(const F90geovals &, const F90geovals &);
  void aq_geovals_rms_f90(const F90geovals &, double &);
  void aq_geovals_dotprod_f90(const F90geovals &, const F90geovals &, double &);
  void aq_geovals_stats_f90(const F90geovals &, int &, double &, double &, double &, double &);
  void aq_geovals_maxloc_f90(const F90geovals &, double &, int &, const oops::Variables &);
  void aq_geovals_read_file_f90(const F90geovals &, const eckit::Configuration &);
  void aq_geovals_write_file_f90(const F90geovals &, const eckit::Configuration &);
  void aq_geovals_analytic_init_f90(const F90geovals &, const Locations &,
                                    const eckit::Configuration &);
}
}  // namespace aq
#endif  // AQ_AQ_GEOVALS_INTERFACE_H_
