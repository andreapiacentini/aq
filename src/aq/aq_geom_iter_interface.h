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

#ifndef AQ_AQ_GEOM_ITER_INTERFACE_H_
#define AQ_AQ_GEOM_ITER_INTERFACE_H_

#include <memory>
#include <string>
#include <vector>

#include "aq/interface.h"

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
