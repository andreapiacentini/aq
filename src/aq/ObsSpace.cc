/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/ObsSpace.h"

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

#include "aq/aq_obsdb_interface.h"

using atlas::array::make_view;

namespace aq {
// -----------------------------------------------------------------------------
// initialization for the static map
std::map < std::string, F90odb > ObsSpace::theObsFileRegister_;
std::map < std::string, F90odb > ObsSpace::theObsDbRegister_;
std::map < std::string, int > ObsSpace::theObsDbPerFile_;
int ObsSpace::theObsFileCount_ = 0;
int ObsSpace::theObsDbCount_ = 0;

// -----------------------------------------------------------------------------

ObsSpace::ObsSpace(const Parameters_ & params, const eckit::mpi::Comm & comm,
                       const util::DateTime & bgn, const util::DateTime & end,
                       const eckit::mpi::Comm & timeComm)
  : oops::ObsSpaceBase(params, comm, bgn, end), obsname_(params.obsType),
  winbgn_(bgn), winend_(end), obsvars_(), comm_(comm)
{
  if (comm_.rank() == 0) {
    typedef std::map< std::string, F90odb >::iterator otiter;
    typedef std::map< std::string, F90odb >::iterator fiter;

    eckit::LocalConfiguration fileconf = params.toConfiguration();
    std::string ofin("-");
    if (params.obsdatain.value() != boost::none) {
      ofin = params.obsdatain.value()->obsfile;
    }
    std::string ofout("-");
    if (params.obsdataout.value() != boost::none) {
      ofout = params.obsdataout.value()->obsfile;
      if (timeComm.size() > 1) {
        std::ostringstream ss;
        ss << "_" << timeComm.rank();
        std::size_t found = ofout.find_last_of(".");
        if (found == std::string::npos) found = ofout.length();
        std::string fileout = ofout.insert(found, ss.str());
        fileconf.set("obsdataout.obsfile", fileout);
      }
    }
    oops::Log::trace() << "ObsSpace: Obs files are: " << ofin << " and " << ofout << std::endl;
    fileref_ = ofin + ofout;
    if (fileref_ == "--") {
      ABORT("Underspecified observation files.");
    }
    std::string dbref = fileref_ + obsname_ + bgn.toString() + end.toString();

    otiter it = theObsDbRegister_.find(dbref);
    if ( it == theObsDbRegister_.end() ) {
      // Setup new db
      oops::Log::trace() << "ObsSpace::getHelper: " << "Setting up " << dbref << std::endl;
      fiter jt = theObsFileRegister_.find(fileref_);
      if ( jt == theObsFileRegister_.end() ) {
        aq_obsdb_setup_f90(key_, fileconf, bgn, end, true, key_);
        theObsFileRegister_[fileref_] = key_;
        theObsDbPerFile_[fileref_] = 1;
        theObsFileCount_++;
      } else {
        F90odb keyold = jt->second;
        aq_obsdb_setup_f90(key_, fileconf, bgn, end, false, keyold);
        theObsDbPerFile_[fileref_]++;
      }
      theObsDbRegister_[dbref] = key_;
      theObsDbCount_++;
    } else {
      // File already open
      oops::Log::trace() << "ObsSpace::getHelper: " << dbref << " already setup." << std::endl;
      key_ = it->second;
    }

    aq_obsdb_read_f90(key_);

    // Set variables simulated for different obstypes
    obsvars_.push_back(obsname_);

    // For this model the processed varaibles are the same as the simulated variables.
    assimvars_ = obsvars_;

    //  Generate locations etc... if required
    if (params.generate.value() != boost::none) {
      const ObsGenerateParameters &gParams = *params.generate.value();
      const util::Duration first(gParams.begin);
      const util::DateTime start(winbgn_ + first);
      const util::Duration freq(gParams.obsPeriod);
      int nobstimes = 0;
      util::DateTime now(start);
      while (now <= winend_) {
        ++nobstimes;
        now += freq;
      }
      int iobs;
      aq_obsdb_generate_f90(key_, obsname_.size(), obsname_.c_str(), gParams.toConfiguration(),
                            start, freq, nobstimes, iobs);
    }
  }
}

// -----------------------------------------------------------------------------

ObsSpace::~ObsSpace() {}

// -----------------------------------------------------------------------------

void ObsSpace::save() const {
  if (comm_.rank() == 0) {
    typedef std::map< std::string, int >::iterator fiter;
    fiter it = theObsDbPerFile_.find(fileref_);
    int dbperfilecount = it->second;
    ASSERT(dbperfilecount > 0);
    theObsDbPerFile_[fileref_]--;
    dbperfilecount--;
    if (dbperfilecount == 0) {
      aq_obsdb_delete_f90(key_, true);
      theObsFileCount_--;
      if (theObsFileCount_ == 0) {
        theObsFileRegister_.clear();
        theObsDbPerFile_.clear();
        theObsDbRegister_.clear();
      }
    } else {
      aq_obsdb_delete_f90(key_, false);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsSpace::getdb(const std::string & col, int & keyData) const {
  if (comm_.rank() == 0) {
    aq_obsdb_get_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
  }
}

// -----------------------------------------------------------------------------

void ObsSpace::putdb(const std::string & col, const int & keyData) const {
  if (comm_.rank() == 0) {
    aq_obsdb_put_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
  }
}

// -----------------------------------------------------------------------------

std::unique_ptr<Locations> ObsSpace::locations() const {
  atlas::FieldSet fields;
  std::vector<util::DateTime> times;
  if (comm_.rank() == 0) {
    aq_obsdb_locations_f90(key_, obsname_.size(), obsname_.c_str(), fields.get(), times);
  }
  return std::unique_ptr<Locations>(new Locations(fields, std::move(times), comm_));
}

// -----------------------------------------------------------------------------

int ObsSpace::nobs() const {
  int iobs = 0;
  if (comm_.rank() == 0) {
    aq_obsdb_nobs_f90(key_, obsname_.size(), obsname_.c_str(), iobs);
  }
  return iobs;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
ObsIterator ObsSpace::begin() const {
  return ObsIterator(*this->locations(), 0);
}
// -----------------------------------------------------------------------------
ObsIterator ObsSpace::end() const {
  return ObsIterator(*this->locations(), this->nobs());
}
// -----------------------------------------------------------------------------

void ObsSpace::print(std::ostream & os) const {
  if (comm_.rank() == 0) {
    os << "ObsSpace for " << obsname_ << ", " << this->nobs() << " obs";
  }
}

// -----------------------------------------------------------------------------

}  // namespace aq
