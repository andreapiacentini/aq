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
