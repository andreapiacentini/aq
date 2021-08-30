/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/ObsSpaceAQ.h"

#include <map>
#include <string>
#include <utility>

#include "atlas/array.h"
#include "atlas/field.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Sphere.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace aq {
// -----------------------------------------------------------------------------
// initialization for the static map
std::map < std::string, F90odb > ObsSpaceAQ::theObsFileRegister_;
int ObsSpaceAQ::theObsFileCount_ = 0;

// -----------------------------------------------------------------------------

ObsSpaceAQ::ObsSpaceAQ(const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                       const util::DateTime & bgn, const util::DateTime & end,
                       const eckit::mpi::Comm & timeComm)
  : oops::ObsSpaceBase(config, comm, bgn, end), obsname_(config.getString("obs type")),
  winbgn_(bgn), winend_(end), obsvars_(), comm_(comm)
{
  typedef std::map< std::string, F90odb >::iterator otiter;

  eckit::LocalConfiguration fileconf(config);
  std::string ofin("-");
  if (config.has("obsdatain")) {
    ofin = config.getString("obsdatain.obsfile");
  }
  std::string ofout("-");
  if (config.has("obsdataout")) {
    ofout = config.getString("obsdataout.obsfile");
    if (timeComm.size() > 1) {
      std::ostringstream ss;
      ss << "_" << timeComm.rank();
      std::size_t found = ofout.find_last_of(".");
      if (found == std::string::npos) found = ofout.length();
      std::string fileout = ofout.insert(found, ss.str());
      fileconf.set("obsdataout.obsfile", fileout);
    }
  }
  oops::Log::trace() << "ObsSpaceAQ: Obs files are: " << ofin << " and " << ofout << std::endl;
  std::string ref = ofin + ofout;
  if (ref == "--") {
    ABORT("Underspecified observation files.");
  }

  ref = ref + bgn.toString() + end.toString();
  // otiter it = theObsFileRegister_.find(ref);
  // if ( it == theObsFileRegister_.end() ) {
    // Open new file
    oops::Log::trace() << "ObsSpaceAQ::getHelper: " << "Opening " << ref << std::endl;
    aq_obsdb_setup_f90(key_, fileconf, bgn, end, &comm_);
    // theObsFileRegister_[ref] = key_;
  // } else {
  //   // File already open
  //   oops::Log::trace() << "ObsSpaceAQ::getHelper: " << ref << " already opened." << std::endl;
  //   key_ = it->second;
  // }
  // theObsFileCount_++;

  // Set variables simulated for different obstypes
  if (obsname_ == "Stream") obsvars_.push_back("Stream");
  if (obsname_ == "O3") obsvars_.push_back("O3");
  if (obsname_ == "CO") obsvars_.push_back("CO");
  if (obsname_ == "WSpeed") obsvars_.push_back("WSpeed");
  if (obsname_ == "Wind") {
    obsvars_.push_back("Uwind");
    obsvars_.push_back("Vwind");
  }

  //  Generate locations etc... if required
  if (config.has("generate")) {
    const eckit::LocalConfiguration gconf(config, "generate");
    const util::Duration first(gconf.getString("begin"));
    const util::DateTime start(winbgn_ + first);
    const util::Duration freq(gconf.getString("obs_period"));
    int nobstimes = 0;
    util::DateTime now(start);
    while (now <= winend_) {
      ++nobstimes;
      now += freq;
    }
    int iobs;
    aq_obsdb_generate_f90(key_, obsname_.size(), obsname_.c_str(), gconf,
                          start, freq, nobstimes, iobs);
  }
}

// -----------------------------------------------------------------------------

ObsSpaceAQ::~ObsSpaceAQ() {}

// -----------------------------------------------------------------------------

void ObsSpaceAQ::save() const {
  // ASSERT(theObsFileCount_ > 0);
  // theObsFileCount_--;
  // if (theObsFileCount_ == 0) {
  //  theObsFileRegister_.clear();
    aq_obsdb_delete_f90(key_);
  // }
}

// -----------------------------------------------------------------------------

void ObsSpaceAQ::getdb(const std::string & col, int & keyData) const {
  aq_obsdb_get_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
}

// -----------------------------------------------------------------------------

void ObsSpaceAQ::putdb(const std::string & col, const int & keyData) const {
  aq_obsdb_put_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
}

// -----------------------------------------------------------------------------

std::unique_ptr<LocationsAQ> ObsSpaceAQ::locations() const {
  atlas::FieldSet fields;
  std::vector<util::DateTime> times;
  aq_obsdb_locations_f90(key_, obsname_.size(), obsname_.c_str(), fields.get(), times);
  return std::unique_ptr<LocationsAQ>(new LocationsAQ(fields, std::move(times)));
}

// -----------------------------------------------------------------------------

int ObsSpaceAQ::nobs() const {
  int iobs;
  aq_obsdb_nobs_f90(key_, obsname_.size(), obsname_.c_str(), iobs);
  return iobs;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
ObsIteratorAQ ObsSpaceAQ::begin() const {
  return ObsIteratorAQ(*this->locations(), 0);
}
// -----------------------------------------------------------------------------
ObsIteratorAQ ObsSpaceAQ::end() const {
  return ObsIteratorAQ(*this->locations(), this->nobs());
}
// -----------------------------------------------------------------------------

void ObsSpaceAQ::print(std::ostream & os) const {
  os << "ObsSpace for " << obsname_ << ", " << this->nobs() << " obs";
}

// -----------------------------------------------------------------------------

}  // namespace aq
