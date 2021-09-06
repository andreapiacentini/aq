/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_LINEAROBSOPBASE_H_
#define AQ_LINEAROBSOPBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

#include "aq/ObsSpace.h"

namespace aq {
class GeoVals;
class ObsAuxControl;
class ObsAuxIncrement;
class ObsVec;

// -----------------------------------------------------------------------------
/// Base class for observation operators

class LinearObsOpBase : public util::Printable,
                        private boost::noncopyable {
 public:
  LinearObsOpBase() = default;

/// Obs Operator
  virtual void setTrajectory(const GeoVals &, const ObsAuxControl &) = 0;
  virtual void simulateObsTL(const GeoVals &, ObsVec &, const ObsAuxIncrement &) const = 0;
  virtual void simulateObsAD(GeoVals &, const ObsVec &, ObsAuxIncrement &) const = 0;

/// Other
  virtual const oops::Variables & requiredVars() const = 0;  // Required from Model

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class LinearObsOpFactory {
 public:
  static LinearObsOpBase * create(const ObsSpace &, const eckit::Configuration &);
  virtual ~LinearObsOpFactory() = default;
 protected:
  explicit LinearObsOpFactory(const std::string &);
 private:
  virtual LinearObsOpBase * make(const ObsSpace &, const eckit::Configuration &) = 0;
  static std::map < std::string, LinearObsOpFactory * > & getMakers() {
    static std::map < std::string, LinearObsOpFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class LinearObsOpMaker : public LinearObsOpFactory {
  virtual LinearObsOpBase * make(const ObsSpace & odb, const eckit::Configuration & conf)
    { return new T(odb, conf); }
 public:
  explicit LinearObsOpMaker(const std::string & name) : LinearObsOpFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_LINEAROBSOPBASE_H_
