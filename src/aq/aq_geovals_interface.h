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
  void aq_geovals_setup_f90(F90geovals &, const Locations &, const oops::Variables &,
                            const eckit::mpi::Comm *);
  void aq_geovals_create_f90(F90geovals &, const oops::Variables &, const eckit::mpi::Comm *);
  void aq_geovals_delete_f90(F90geovals &);
  void aq_geovals_copy_f90(const F90geovals &, const F90geovals &);
  void aq_geovals_fill_f90(const F90geovals &, const int &, const int &, const int &,
                           const double &);
  void aq_geovals_fillad_f90(const F90geovals &, const int &, const int &, const int &, double &);
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
