/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_OBSSPACEAQ_H_
#define AQ_MODEL_OBSSPACEAQ_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/geometry/Point2.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/ObsSpaceBase.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "oops/aq/AqFortran.h"
#include "oops/aq/LocationsAQ.h"
#include "oops/aq/ObsIteratorAQ.h"

namespace eckit {
  class Configuration;
}

namespace aq {
  class ObsIteratorAQ;

/// \brief ObsSpace for AQ model
//  \details ObsSpaceAQ is created for each obs type. The underlying Fortran
//  structure (key_) is created for each matching input-output filename pair
//  (i.e. different obstypes can be stored in the same Fortran structure).
//  For mapping between ObsSpaceAQ and Fortran structures,
//  ObsSpaceAQ::theObsFileRegister_ map is used
class ObsSpaceAQ : public oops::ObsSpaceBase {
 public:
  /// create full ObsSpace (read or generate data)
  ObsSpaceAQ(const eckit::Configuration &, const eckit::mpi::Comm &,
             const util::DateTime &, const util::DateTime &, const eckit::mpi::Comm &);
  ~ObsSpaceAQ();

  /// save and close file
  void save() const;

  /// read data or metadata
  void getdb(const std::string &, int &) const;
  /// save data or metadata
  void putdb(const std::string &, const int &) const;

  /// create locations for the whole time window
  std::unique_ptr<LocationsAQ> locations() const;

  /// return number of observations (unique locations)
  int nobs() const;

  /// return variables simulated by ObsOperators
  const oops::Variables & obsvariables() const { return obsvars_; }

  /// observation type
  const std::string & obsname() const {return obsname_;}

  /// iterator to the first observation
  ObsIteratorAQ begin() const;
  /// iterator to the observation past-the-last
  ObsIteratorAQ end() const;

  /// interface with Fortran
  const F90odb & toFortran() const {return key_;}

 private:
  void print(std::ostream &) const;

  mutable F90odb key_;               // pointer to Fortran structure
  const std::string obsname_;        // corresponds with obstype
  const util::DateTime winbgn_;      // window for the observations
  const util::DateTime winend_;
  oops::Variables obsvars_;          // variables simulated by ObsOperators
  const eckit::mpi::Comm & comm_;

  // defines mapping for Fortran structures
  static std::map < std::string, F90odb > theObsFileRegister_;
  static int theObsFileCount_;  // number of files used
};

}  // namespace aq

#endif  // AQ_MODEL_OBSSPACEAQ_H_
