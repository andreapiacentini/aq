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

#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"

#include "aq/aq_geom_iter_interface.h"
#include "aq/GeometryIterator.h"
#include "aq/interface.h"

// -----------------------------------------------------------------------------

namespace aq {

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator& iter) {
  aq_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry& geom, const int & index) {
  aq_geom_iter_setup_f90(keyIter_, geom.toFortran(), index);
}


// -----------------------------------------------------------------------------

GeometryIterator::~GeometryIterator() {
  aq_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator==(const GeometryIterator & other) const {
  int equals = 0;
  aq_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 1);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator!=(const GeometryIterator & other) const {
  int equals = 0;
  aq_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 0);
}

// -----------------------------------------------------------------------------

eckit::geometry::Point2 GeometryIterator::operator*() const {
  double lat, lon;
  aq_geom_iter_current_f90(keyIter_, lat, lon);
  return eckit::geometry::Point2(lat, lon);
}

// -----------------------------------------------------------------------------

GeometryIterator& GeometryIterator::operator++() {
  aq_geom_iter_next_f90(keyIter_);
  return *this;
}

// -----------------------------------------------------------------------------

void GeometryIterator::print(std::ostream & os) const {
  os << "GeometryIterator, key: " <<  keyIter_;
}

// -----------------------------------------------------------------------------

}  // namespace aq
