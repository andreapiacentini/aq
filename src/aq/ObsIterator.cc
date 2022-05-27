/*
 * (C) Copyright 2021 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "aq/ObsIterator.h"

// -----------------------------------------------------------------------------
namespace aq {

ObsIterator::ObsIterator(const ObsIterator & other): index_(other.index_),
      locslonlat_(other.locslonlat_) {}

// -----------------------------------------------------------------------------
ObsIterator::ObsIterator(const Locations & locations, int index):
      index_(index), locslonlat_(locations.lonlat()) {}

// -----------------------------------------------------------------------------
bool ObsIterator::operator==(const ObsIterator & other) const {
  return (index_ == other.index_);
}

// -----------------------------------------------------------------------------
bool ObsIterator::operator!=(const ObsIterator & other) const {
  return (index_!= other.index_);
}

// -----------------------------------------------------------------------------
eckit::geometry::Point2 ObsIterator::operator*() const {
  auto lonlat = atlas::array::make_view<double, 2>(locslonlat_);
  return eckit::geometry::Point2(lonlat(index_, 0), lonlat(index_, 1));
}

// -----------------------------------------------------------------------------
ObsIterator& ObsIterator::operator++() {
  index_++;
  return *this;
}

// -----------------------------------------------------------------------------

}  // namespace aq
