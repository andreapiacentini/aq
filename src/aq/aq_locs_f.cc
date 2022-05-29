/*
 * (C) Copyright 2020 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
