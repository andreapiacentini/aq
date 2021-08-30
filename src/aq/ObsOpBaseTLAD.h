/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_MODEL_OBSOPBASETLAD_H_
#define AQ_MODEL_OBSOPBASETLAD_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

#include "oops/aq/ObsSpaceAQ.h"

namespace aq {
class GomAQ;
class ObsBias;
class ObsBiasIncrement;
class ObsVecAQ;

// -----------------------------------------------------------------------------
/// Base class for observation operators

class ObsOpBaseTLAD : public util::Printable,
                      private boost::noncopyable {
 public:
  ObsOpBaseTLAD() = default;

/// Obs Operator
  virtual void setTrajectory(const GomAQ &, const ObsBias &) = 0;
  virtual void simulateObsTL(const GomAQ &, ObsVecAQ &, const ObsBiasIncrement &) const = 0;
  virtual void simulateObsAD(GomAQ &, const ObsVecAQ &, ObsBiasIncrement &) const = 0;

/// Other
  virtual const oops::Variables & requiredVars() const = 0;  // Required from Model

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class ObsOpTLADFactory {
 public:
  static ObsOpBaseTLAD * create(const ObsSpaceAQ &, const eckit::Configuration &);
  virtual ~ObsOpTLADFactory() = default;
 protected:
  explicit ObsOpTLADFactory(const std::string &);
 private:
  virtual ObsOpBaseTLAD * make(const ObsSpaceAQ &, const eckit::Configuration &) = 0;
  static std::map < std::string, ObsOpTLADFactory * > & getMakers() {
    static std::map < std::string, ObsOpTLADFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsOpTLADMaker : public ObsOpTLADFactory {
  virtual ObsOpBaseTLAD * make(const ObsSpaceAQ & odb, const eckit::Configuration & conf)
    { return new T(odb, conf); }
 public:
  explicit ObsOpTLADMaker(const std::string & name) : ObsOpTLADFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODEL_OBSOPBASETLAD_H_
