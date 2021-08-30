/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_LOCATIONSAQ_H_
#define AQ_MODEL_LOCATIONSAQ_H_

#include <iomanip>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/aq/AqFortran.h"

namespace aq {

// -----------------------------------------------------------------------------
/// LocationsAQ class to handle locations for AQ model.
class LocationsAQ : public util::Printable,
                    private util::ObjectCounter<LocationsAQ> {
 public:
  static const std::string classname() {return "aq::LocationsAQ";}

// Constructors and basic operators
  LocationsAQ(atlas::FieldSet &, std::vector<util::DateTime> &&);
  LocationsAQ(const eckit::Configuration &, const eckit::mpi::Comm &);
  LocationsAQ(const LocationsAQ &);
  ~LocationsAQ() {}

// Utilities
  int size() const {return pointcloud_->size();}
  atlas::functionspace::PointCloud & pointcloud() {return *pointcloud_;}
  atlas::Field lonlat() const {return pointcloud_->lonlat();}
  atlas::Field & altitude() {ASSERT(altitude_); return *altitude_;}
  util::DateTime & times(size_t idx) {return times_[idx];}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<atlas::functionspace::PointCloud> pointcloud_;
  std::unique_ptr<atlas::Field> altitude_;
  std::vector<util::DateTime> times_;
};

}  // namespace aq

#endif  // AQ_MODEL_LOCATIONSAQ_H_
