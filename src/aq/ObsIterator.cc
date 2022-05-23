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
