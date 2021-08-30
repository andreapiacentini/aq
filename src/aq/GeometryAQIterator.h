/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_GEOMETRYAQITERATOR_H_
#define AQ_GEOMETRYAQITERATOR_H_

#include <iterator>
#include <string>

#include "eckit/geometry/Point2.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/AqFortran.h"
#include "aq/GeometryAQ.h"

namespace aq {

class GeometryAQ;

// -----------------------------------------------------------------------------
class GeometryAQIterator: public std::iterator<std::forward_iterator_tag,
                                               eckit::geometry::Point2>,
                          public util::Printable,
                          private util::ObjectCounter<GeometryAQIterator> {
 public:
  static const std::string classname() {return "aq::GeometryAQIterator";}

  GeometryAQIterator(const GeometryAQIterator &);
  explicit GeometryAQIterator(const GeometryAQ & geom, const int & index = 1);
  ~GeometryAQIterator();

  bool operator==(const GeometryAQIterator &) const;
  bool operator!=(const GeometryAQIterator &) const;
  eckit::geometry::Point2 operator*() const;
  GeometryAQIterator& operator++();

  const F90iter & toFortran() const {return keyIter_;}

 private:
  void print(std::ostream &) const;
  F90iter keyIter_;
};

}  // namespace aq

#endif  // AQ_GEOMETRYAQITERATOR_H_
