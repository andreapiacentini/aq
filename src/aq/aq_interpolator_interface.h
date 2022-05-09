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

#ifndef AQ_AQ_INTERPOLATOR_INTERFACE_H_
#define AQ_AQ_INTERPOLATOR_INTERFACE_H_

#include "aq/interface.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace aq {

extern "C" {
  void aq_interpolator_create_f90(F90interp &, const F90geom &,
                                  const int &, const double &, const double &);
  void aq_interpolator_delete_f90(F90interp &);
  void aq_interpolator_apply_f90(const F90interp &, const F90flds &,
                                 const oops::Variables &,
                                 const int &, const bool &, double &);
  void aq_interpolator_applyAD_f90(const F90interp &, const F90flds &,
                                   const oops::Variables &,
                                   const int &, const bool &, const double &);
}

}  // namespace aq
#endif  // AQ_AQ_INTERPOLATOR_INTERFACE_H_
