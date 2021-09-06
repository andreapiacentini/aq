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

#ifndef AQ_AQ_GETVALUES_INTERFACE_H_
#define AQ_AQ_GETVALUES_INTERFACE_H_

#include <memory>
#include <string>
#include <vector>

#include "aq/interface.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace aq {
  class Locations;

extern "C" {
  void aq_getvalues_setup_f90(const F90hmat &);
  void aq_getvalues_delete_f90(const F90hmat &);
  void aq_getvalues_interp_f90(const Locations &, const F90flds &,
                               const util::DateTime &,
                               const util::DateTime &, const F90geovals &);
  void aq_getvalues_build_f90(const Locations &, const F90flds &,
                              const util::DateTime &,
                              const util::DateTime &, const F90hmat &);
  void aq_getvalues_interp_tl_f90(const Locations &, const F90flds &,
                                  const util::DateTime &,
                                  const util::DateTime &, const F90hmat &, const F90geovals &);
  void aq_getvalues_interp_ad_f90(const Locations &, const F90flds &,
                                  const util::DateTime &,
                                  const util::DateTime &, const F90hmat &, const F90geovals &);
}

}  // namespace aq
#endif  // AQ_AQ_GETVALUES_INTERFACE_H_
