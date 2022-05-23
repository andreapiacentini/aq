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

#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"

#include "aq/aq_locs_f.h"
#include "aq/Locations.h"
#include "oops/util/DateTime.h"

namespace aq {

// -----------------------------------------------------------------------------
int aq_locs_nlocs_f90(aq::Locations* locs) {
    return locs->size();
}
atlas::field::FieldImpl* aq_locs_lonlat_f90(aq::Locations* locs) {
    return locs->lonlat().get();
}
atlas::field::FieldImpl* aq_locs_altitude_f90(aq::Locations* locs) {
    return locs->altitude().get();
}
util::DateTime& aq_locs_times_f90(aq::Locations* locs, size_t & idx) {
    return locs->times(idx);
}
// -----------------------------------------------------------------------------

}  // namespace aq
