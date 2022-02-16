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

#ifndef AQ_OBSSPACE_H_
#define AQ_OBSSPACE_H_

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
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "aq/interface.h"
#include "aq/Locations.h"
#include "aq/ObsIterator.h"

namespace aq {
  class ObsIterator;

/// Contents of the `obsdatain` or `obsdataout` YAML section.
class ObsDataParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsDataParameters, Parameters)

 public:
  /// File path.
  oops::RequiredParameter<std::string> obsfile{"obsfile", this};
};

/// Contents of the `obserrors` YAML section.
class ObsErrorsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorsParameters, Parameters)

 public:
  /// File path.
  oops::RequiredParameter<std::string> type{"type", this};
  oops::OptionalParameter<double> value{"value", this};
};

/// Options controlling generation of artificial observations.
class ObsGenerateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsGenerateParameters, Parameters)

 public:
  oops::RequiredParameter<util::Duration> begin{"begin", this};
  oops::RequiredParameter<util::Duration> obsPeriod{"obs_period", this};
  /// Number of observations to generate in each time slot.
  oops::RequiredParameter<int> obsDensity{"obs_density", this};
  oops::RequiredParameter<int> nval{"nval", this};
  oops::RequiredParameter<double> obsError{"obs_error", this};
};

/// Options controlling the variable tranformations.
class ObsTransfvarParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsTransfvarParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> method{"method", this};
  oops::RequiredParameter<std::vector<double>> parameters{"parameters", this};
};

/// \brief Configuration parameters for the AQ application's ObsSpace.
class ObsSpaceParameters : public oops::ObsSpaceParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsSpaceParameters, ObsSpaceParametersBase)

 public:
  /// Type of observations.
  oops::RequiredParameter<std::string> obsType{"obs type", this};
  /// Instrument Name.
  oops::RequiredParameter<std::string> instName{"instr name", this};
  /// File from which to load observations.
  oops::OptionalParameter<ObsDataParameters> obsdatain{"obsdatain", this};
  /// File to which to save observations and analysis.
  oops::OptionalParameter<ObsDataParameters> obsdataout{"obsdataout", this};
  /// Specification of observation errors.
  oops::OptionalParameter<ObsErrorsParameters> obserrors{"obserrors", this};
  /// Options controlling generation of artificial observations.
  oops::OptionalParameter<ObsGenerateParameters> generate{"generate", this};
  /// Options controlling the variable transformations.
  oops::OptionalParameter<ObsTransfvarParameters> transfvar{"transfvar", this};
};

/// \brief ObsSpace for AQ model
/// \details ObsSpace is created for each obs type. The underlying Fortran
/// structure (key_) is created for each matching input-output filename pair
/// (i.e. different obstypes can be stored in the same Fortran structure).
/// For mapping between ObsSpace and Fortran structures,
/// ObsSpace::theObsFileRegister_ map is used
class ObsSpace : public oops::ObsSpaceBase {
 public:
  typedef ObsSpaceParameters Parameters_;

  /// create full ObsSpace (read or generate data)
  ObsSpace(const Parameters_ &, const eckit::mpi::Comm &,
             const util::DateTime &, const util::DateTime &, const eckit::mpi::Comm &);
  ~ObsSpace();

  /// save and close file
  void save() const;

  /// read data or metadata
  void getdb(const std::string &, int &) const;
  /// save data or metadata
  void putdb(const std::string &, const int &) const;

  /// create locations for the whole time window
  std::unique_ptr<Locations> locations() const;

  /// return number of observations (unique locations)
  int nobs() const;

  /// return variables simulated by ObsOperators
  const oops::Variables & obsvariables() const { return obsvars_; }

  /// observation type
  const std::string & obsname() const {return obsname_;}

  /// communicator
  const eckit::mpi::Comm & comm() const {return comm_;}

  /// iterator to the first observation
  ObsIterator begin() const;
  /// iterator to the observation past-the-last
  ObsIterator end() const;

  /// interface with Fortran
  const F90odb & toFortran() const {return key_;}

 private:
  void print(std::ostream &) const;

  mutable F90odb key_;               // pointer to Fortran structure
  const std::string obsname_;        // corresponds with obstype
  const util::DateTime winbgn_;      // window for the observations
  const util::DateTime winend_;
  std::string fileref_;              // Reference to the observation file
  oops::Variables obsvars_;          // variables simulated by ObsOperators
  const eckit::mpi::Comm & comm_;

  // defines mapping for Fortran structures
  static std::map < std::string, F90odb > theObsFileRegister_;
  static std::map < std::string, int > theObsDbPerFile_;  // number of db per file
  static std::map < std::string, F90odb > theObsDbRegister_;
  static int theObsDbCount_;  // number of db used
  static int theObsFileCount_;  // number of files used
};

}  // namespace aq

#endif  // AQ_OBSSPACE_H_
