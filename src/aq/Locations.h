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

#ifndef AQ_LOCATIONS_H_
#define AQ_LOCATIONS_H_

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

namespace aq {

// -----------------------------------------------------------------------------
/// Locations class to handle locations for AQ model.
class Locations : public util::Printable,
                    private util::ObjectCounter<Locations> {
 public:
  static const std::string classname() {return "aq::Locations";}

// Constructors and basic operators
  Locations(atlas::FieldSet &, std::vector<util::DateTime> &&, const eckit::mpi::Comm &);
  Locations(const eckit::Configuration &, const eckit::mpi::Comm &);
  Locations(const Locations &);
  ~Locations() {}

// Utilities
  int size() const {return pointcloud_->size();}
  atlas::functionspace::PointCloud & pointcloud() {return *pointcloud_;}
  atlas::Field lonlat() const {return pointcloud_->lonlat();}
  atlas::Field & altitude() {ASSERT(altitude_); return *altitude_;}
  const std::vector<double> & latitudes() const {return lats_;}
  const std::vector<double> & longitudes() const {return lons_;}
  const std::vector<util::DateTime> & times() const {return times_;}
  util::DateTime & times(size_t idx) {return times_[idx];}

  /// communicator
  const eckit::mpi::Comm & comm() const {return comm_;}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<atlas::functionspace::PointCloud> pointcloud_;
  std::unique_ptr<atlas::Field> altitude_;
  std::vector<util::DateTime> times_;
  std::vector<double> lats_;
  std::vector<double> lons_;
  const eckit::mpi::Comm & comm_;
};

}  // namespace aq

#endif  // AQ_LOCATIONS_H_
