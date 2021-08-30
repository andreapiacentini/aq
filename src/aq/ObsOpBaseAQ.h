/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_OBSOPBASEAQ_H_
#define AQ_OBSOPBASEAQ_H_

#include <map>
#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

#include "aq/ObsSpaceAQ.h"

namespace aq {
class GomAQ;
class LocationsAQ;
class ObsBias;
class ObsVecAQ;

// -----------------------------------------------------------------------------
/// Base class for observation operators

class ObsOpBaseAQ : public util::Printable,
                    private boost::noncopyable {
 public:
  ObsOpBaseAQ() = default;

/// Obs Operator
  virtual void simulateObs(const GomAQ &, ObsVecAQ &, const ObsBias &) const = 0;

/// Other
  virtual const oops::Variables & requiredVars() const = 0;  // Required from Model
  virtual std::unique_ptr<LocationsAQ> locations() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class ObsOpFactory {
 public:
  static ObsOpBaseAQ * create(const ObsSpaceAQ &, const eckit::Configuration &);
  virtual ~ObsOpFactory() = default;
 protected:
  explicit ObsOpFactory(const std::string &);
 private:
  virtual ObsOpBaseAQ * make(const ObsSpaceAQ &, const eckit::Configuration &) = 0;
  static std::map < std::string, ObsOpFactory * > & getMakers() {
    static std::map < std::string, ObsOpFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsOpMaker : public ObsOpFactory {
  virtual ObsOpBaseAQ * make(const ObsSpaceAQ & odb, const eckit::Configuration & conf)
    { return new T(odb, conf); }
 public:
  explicit ObsOpMaker(const std::string & name) : ObsOpFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSOPBASEAQ_H_
