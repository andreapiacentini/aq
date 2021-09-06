/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_AQ_LOCS_F_H_
#define AQ_AQ_LOCS_F_H_

#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"

#include "aq/Locations.h"
#include "oops/util/DateTime.h"

// ------------------------------------------------------------------------------
// These functions provide tools for interfacing Fortran and C++ string objects
// ------------------------------------------------------------------------------

namespace aq {

extern "C" {
  int aq_locs_nlocs_f90(aq::Locations*);
  atlas::field::FieldImpl* aq_locs_lonlat_f90(aq::Locations*);
  atlas::field::FieldImpl* aq_locs_altitude_f90(aq::Locations*);
  util::DateTime& aq_locs_times_f90(aq::Locations*, size_t &);
}

}  // namespace aq

#endif  // AQ_AQ_LOCS_F_H_
