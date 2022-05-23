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

#ifndef AQ_GEOMETRYITERATOR_H_
#define AQ_GEOMETRYITERATOR_H_

#include <iterator>
#include <string>

#include "eckit/geometry/Point2.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/Geometry.h"
#include "aq/interface.h"

namespace aq {

class Geometry;

// -----------------------------------------------------------------------------
class GeometryIterator: public std::iterator<std::forward_iterator_tag,
                                               eckit::geometry::Point2>,
                          public util::Printable,
                          private util::ObjectCounter<GeometryIterator> {
 public:
  static const std::string classname() {return "aq::GeometryIterator";}

  GeometryIterator(const GeometryIterator &);
  explicit GeometryIterator(const Geometry & geom, const int & index = 1);
  ~GeometryIterator();

  bool operator==(const GeometryIterator &) const;
  bool operator!=(const GeometryIterator &) const;
  eckit::geometry::Point2 operator*() const;
  GeometryIterator& operator++();

  const F90iter & toFortran() const {return keyIter_;}

 private:
  void print(std::ostream &) const;
  F90iter keyIter_;
};

}  // namespace aq

#endif  // AQ_GEOMETRYITERATOR_H_
