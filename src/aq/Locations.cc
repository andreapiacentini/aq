/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "aq/Locations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Random.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;
using atlas::option::name;

namespace aq {
// -------------------------------------------------------------------------
/*! This constructor creates test locations based on settings in the config file */
Locations::Locations(const eckit::Configuration & config, const eckit::mpi::Comm & comm)
: comm_(comm)
{
  if (comm_.rank() == 0) {
    std::vector<double> z;
    unsigned int nrandom;
    if (config.has("lats") || config.has("nrandom")) {
      /*! These are optional lists of locations that the user can
      * specify in the config file for testing
      */
      lons_ = config.getDoubleVector("lons");
      lats_ = config.getDoubleVector("lats");

      ASSERT(lons_.size() == lats_.size());
      for (std::size_t jj=0; jj < lons_.size(); ++jj) {
        z.push_back(0.);
      }

      /*! Instead of or in addition to the specific locations in the config file
      * let the user specify a number of random locations
      */
      nrandom = config.getInt("nrandom", 0);

    } else {
      /*! If the config file does not specify otherwise, default to 4 random locations */
      nrandom = 4;
    }

    if (nrandom > 0) {
      /*! To create random locations, we need to know the dimensions of the
      * computational domain
      */
      double lonmin = config.getDouble("random lon min");
      double lonmax = config.getDouble("random lon max");
      double latmin = config.getDouble("random lat min");
      double latmax = config.getDouble("random lat max");

      /*! Now create random locations
      * use a specified random seed for reproducibility
      */

      std::uint32_t rseed = 11;
      util::UniformDistribution<double> randlons(nrandom, lonmin, lonmax, rseed);
      util::UniformDistribution<double> randlats(nrandom, latmin, latmax);
      randlons.sort();

      for (std::size_t jj=0; jj < nrandom; ++jj) {
        lons_.push_back(randlons[jj]);
        lats_.push_back(randlats[jj]);
        z.push_back(0.);
      }
    }

    const unsigned int nlocs = lons_.size();
    ASSERT(nlocs > 0);

    /*! render lat, lon as an atlas functionspace */
    atlas::Field field_lonlat("lonlat", make_datatype<double>(), make_shape(nlocs, 2));
    auto lonlat = make_view<double, 2>(field_lonlat);
    for ( unsigned int j = 0; j < nlocs; ++j ) {
        lonlat(j, 0) = lons_[j];
        lonlat(j, 1) = lats_[j];
    }
    pointcloud_.reset(new atlas::functionspace::PointCloud(field_lonlat));

    /*! render levels as an atlas field */
    altitude_.reset(new atlas::Field("altitude", make_datatype<double>(),
                     make_shape(nlocs)));
    auto zpc = make_view<double, 1>(*altitude_);
    for (unsigned int j = 0; j < nlocs; ++j) {
        zpc(j) = z[j];
    }

    /*! Now get times */
    if (config.has("window end")) {
      util::DateTime winbgn(config.getString("window begin"));
      util::DateTime winend(config.getString("window end"));
      util::Duration window_len(winend - winbgn);
      util::Duration dt = window_len / nlocs;

      times_.clear();
      for (unsigned int j=0; j < nlocs; ++j) {
        times_.push_back(winbgn + (j+1)*dt);
      }
    }
  } else {
    const unsigned int nlocs = 0;
    atlas::Field field_lonlat("lonlat", make_datatype<double>(), make_shape(nlocs, 2));
    pointcloud_.reset(new atlas::functionspace::PointCloud(field_lonlat));
    altitude_.reset(new atlas::Field("altitude", make_datatype<double>(),
                     make_shape(nlocs)));
  }
}

// -------------------------------------------------------------------------
/*! Deep copy constructor */
Locations::Locations(const Locations & other)
: comm_(other.comm_)
{
  pointcloud_.reset(new atlas::functionspace::PointCloud(other.lonlat()));
  altitude_.reset(new atlas::Field(*other.altitude_));
  times_ = other.times_;
  lats_ = other.lats_;
  lons_ = other.lons_;
}
// -------------------------------------------------------------------------
/*! Constructor from fields and times.  These may be obtained from the obsdb,
 * passed to C++ from Fortran
 */
Locations::Locations(atlas::FieldSet & fields,
                         std::vector<util::DateTime> && times, const eckit::mpi::Comm & comm)
: times_(times), comm_(comm)
{
  if (comm_.rank() == 0) {
    pointcloud_.reset(new atlas::functionspace::PointCloud(fields.field("lonlat")));
    if (fields.has_field("altitude")) {
      altitude_.reset(new atlas::Field(fields.field("altitude")));
    } else {
      altitude_.reset(new atlas::Field(atlas::Field()));
    }
  } else {
    const unsigned int nlocs = 0;
    atlas::Field field_lonlat("lonlat", make_datatype<double>(), make_shape(nlocs, 2));
    pointcloud_.reset(new atlas::functionspace::PointCloud(field_lonlat));
    altitude_.reset(new atlas::Field("altitude", make_datatype<double>(),
                     make_shape(nlocs)));
  }
  lats_.resize(pointcloud_->size());
  lons_.resize(pointcloud_->size());
  atlas::Field field_lonlat = pointcloud_->lonlat();
  auto lonlat = make_view<double, 2>(field_lonlat);
  for (size_t jloc = 0; jloc < pointcloud_->size(); ++jloc) {
    lats_[jloc] = lonlat(jloc, 1);
    lons_[jloc] = lonlat(jloc, 0);
  }
}
// -------------------------------------------------------------------------
void Locations::localCoords(const util::DateTime & t1, const util::DateTime & t2,
                            std::vector<double> & lats, std::vector<double> & lons,
                            std::vector<size_t> & indx) const {
  lats.clear();
  lons.clear();
  indx.clear();
  int nobs = pointcloud_->size();
  atlas::Field field_lonlat = pointcloud_->lonlat();
  auto lonlat = make_view<double, 2>(field_lonlat);
  for (size_t jj=0; jj < static_cast<size_t>(nobs); ++jj) {
    if (times_[jj] >= t1 && times_[jj] <= t2) {
      lats.push_back(lonlat(jj, 1));
      lons.push_back(lonlat(jj, 0));
      indx.push_back(jj);
    }
  }
}
// -------------------------------------------------------------------------
void Locations::print(std::ostream & os) const {
  int nobs = pointcloud_->size();
  atlas::Field field_lonlat = pointcloud_->lonlat();
  auto lonlat = make_view<double, 2>(field_lonlat);
  auto z = make_view<double, 1>(*altitude_);
  for (size_t jj=0; jj < static_cast<size_t>(nobs); ++jj) {
    os << "location " << jj << std::setprecision(2) << ": lon = " << lonlat(jj, 0)
       << ", lat = " << lonlat(jj, 1) << ", z = " << z(jj) << std::endl;
  }
}
// -------------------------------------------------------------------------
}  // namespace aq
