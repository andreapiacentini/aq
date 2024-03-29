/*
 * (C) Copyright 2017-2018 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_OBSOPBASE_H_
#define AQ_OBSOPBASE_H_

#include <map>
#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

#include "aq/ObsSpace.h"

namespace oops {
template <typename OBS> class Locations;
}

namespace aq {
class GeoVals;
class ObsAuxControl;
class ObsVec;
struct ObsTraits;

// -----------------------------------------------------------------------------
/// Base class for observation operators

class ObsOpBase : public util::Printable,
                    private boost::noncopyable {
 public:
  typedef oops::Locations<ObsTraits> Locations_;

  ObsOpBase() = default;

/// Obs Operator
  virtual void simulateObs(const GeoVals &, ObsVec &, const ObsAuxControl &) const = 0;

/// Other
  virtual const oops::Variables & requiredVars() const = 0;  // Required from Model
  virtual Locations_ locations() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class ObsOpFactory {
 public:
  static ObsOpBase * create(const ObsSpace &, const eckit::Configuration &);
  virtual ~ObsOpFactory() = default;
 protected:
  explicit ObsOpFactory(const std::string &);
 private:
  virtual ObsOpBase * make(const ObsSpace &, const eckit::Configuration &) = 0;
  static std::map < std::string, ObsOpFactory * > & getMakers() {
    static std::map < std::string, ObsOpFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsOpMaker : public ObsOpFactory {
  virtual ObsOpBase * make(const ObsSpace & odb, const eckit::Configuration & conf)
    { return new T(odb, conf); }
 public:
  explicit ObsOpMaker(const std::string & name) : ObsOpFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSOPBASE_H_
