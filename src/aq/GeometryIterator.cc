/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
