/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/config/Configuration.h"
#include "model/AqFortran.h"
#include "model/GeometryAQIterator.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace aq {

// -----------------------------------------------------------------------------

GeometryAQIterator::GeometryAQIterator(const GeometryAQIterator& iter) {
  aq_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryAQIterator::GeometryAQIterator(const GeometryAQ& geom, const int & index) {
  aq_geom_iter_setup_f90(keyIter_, geom.toFortran(), index);
}


// -----------------------------------------------------------------------------

GeometryAQIterator::~GeometryAQIterator() {
  aq_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometryAQIterator::operator==(const GeometryAQIterator & other) const {
  int equals = 0;
  aq_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 1);
}

// -----------------------------------------------------------------------------

bool GeometryAQIterator::operator!=(const GeometryAQIterator & other) const {
  int equals = 0;
  aq_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 0);
}

// -----------------------------------------------------------------------------

eckit::geometry::Point2 GeometryAQIterator::operator*() const {
  double lat, lon;
  aq_geom_iter_current_f90(keyIter_, lat, lon);
  return eckit::geometry::Point2(lat, lon);
}

// -----------------------------------------------------------------------------

GeometryAQIterator& GeometryAQIterator::operator++() {
  aq_geom_iter_next_f90(keyIter_);
  return *this;
}

// -----------------------------------------------------------------------------

void GeometryAQIterator::print(std::ostream & os) const {
  os << "GeometryAQIterator, key: " <<  keyIter_;
}

// -----------------------------------------------------------------------------

}  // namespace aq
