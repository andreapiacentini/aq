/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "aq/ObsIteratorAQ.h"

// -----------------------------------------------------------------------------
namespace aq {

ObsIteratorAQ::ObsIteratorAQ(const ObsIteratorAQ & other): index_(other.index_),
      locslonlat_(other.locslonlat_) {}

// -----------------------------------------------------------------------------
ObsIteratorAQ::ObsIteratorAQ(const LocationsAQ & locations, int index):
      index_(index), locslonlat_(locations.lonlat()) {}

// -----------------------------------------------------------------------------
bool ObsIteratorAQ::operator==(const ObsIteratorAQ & other) const {
  return (index_ == other.index_);
}

// -----------------------------------------------------------------------------
bool ObsIteratorAQ::operator!=(const ObsIteratorAQ & other) const {
  return (index_!= other.index_);
}

// -----------------------------------------------------------------------------
eckit::geometry::Point2 ObsIteratorAQ::operator*() const {
  auto lonlat = atlas::array::make_view<double, 2>(locslonlat_);
  return eckit::geometry::Point2(lonlat(index_, 0), lonlat(index_, 1));
}

// -----------------------------------------------------------------------------
ObsIteratorAQ& ObsIteratorAQ::operator++() {
  index_++;
  return *this;
}

// -----------------------------------------------------------------------------

}  // namespace aq
