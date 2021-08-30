/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_MODEL_OBSITERATORAQ_H_
#define AQ_MODEL_OBSITERATORAQ_H_

#include <iterator>
#include <memory>
#include <string>

#include "eckit/geometry/Point2.h"

#include "aq/LocationsAQ.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace aq {

/// Iterator over all observations
class ObsIteratorAQ: public std::iterator<std::forward_iterator_tag,
                                       eckit::geometry::Point2>,
                   public util::Printable,
                   private util::ObjectCounter<ObsIteratorAQ> {
 public:
  static const std::string classname() {return "aq::ObsIteratorAQ";}

  ObsIteratorAQ(const ObsIteratorAQ &);
  ObsIteratorAQ(const LocationsAQ &, int);

  bool operator==(const ObsIteratorAQ &) const;
  bool operator!=(const ObsIteratorAQ &) const;

  /// return location of current observation
  eckit::geometry::Point2 operator*() const;

  ObsIteratorAQ& operator++();

 private:
  void print(std::ostream & os) const override {os << index_;}

  /// index of a current observation
  int index_;
  /// atlas field of the lons and lats
  atlas::Field locslonlat_;
};

}  // namespace aq

#endif  // AQ_MODEL_OBSITERATORAQ_H_
