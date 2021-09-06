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

#ifndef AQ_AQ_INSITU_INTERFACE_H_
#define AQ_AQ_INSITU_INTERFACE_H_

#include "aq/interface.h"

namespace aq {
extern "C" {
  void aq_insitu_equiv_f90(const F90geovals &, const F90ovec &, const double &);
  void aq_insitu_equiv_tl_f90(const F90geovals &, const F90ovec &, const double &);
  void aq_insitu_equiv_ad_f90(const F90geovals &, const F90ovec &, double &);
}

}  // namespace aq
#endif  // AQ_AQ_INSITU_INTERFACE_H_
