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

#ifndef AQ_OBSITERATOR_H_
#define AQ_OBSITERATOR_H_

#include <iterator>
#include <memory>
#include <string>

#include "eckit/geometry/Point2.h"

#include "aq/Locations.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace aq {

/// Iterator over all observations
class ObsIterator: public std::iterator<std::forward_iterator_tag,
                                       eckit::geometry::Point2>,
                   public util::Printable,
                   private util::ObjectCounter<ObsIterator> {
 public:
  static const std::string classname() {return "aq::ObsIterator";}

  ObsIterator(const ObsIterator &);
  ObsIterator(const Locations &, int);

  bool operator==(const ObsIterator &) const;
  bool operator!=(const ObsIterator &) const;

  /// return location of current observation
  eckit::geometry::Point2 operator*() const;

  ObsIterator& operator++();

 private:
  void print(std::ostream & os) const override {os << index_;}

  /// index of a current observation
  int index_;
  /// atlas field of the lons and lats
  atlas::Field locslonlat_;
};

}  // namespace aq

#endif  // AQ_OBSITERATOR_H_
