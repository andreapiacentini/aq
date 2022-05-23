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

namespace aq {
class GeoVals;
class Locations;
class ObsAuxControl;
class ObsVec;

// -----------------------------------------------------------------------------
/// Base class for observation operators

class ObsOpBase : public util::Printable,
                    private boost::noncopyable {
 public:
  ObsOpBase() = default;

/// Obs Operator
  virtual void simulateObs(const GeoVals &, ObsVec &, const ObsAuxControl &) const = 0;

/// Other
  virtual const oops::Variables & requiredVars() const = 0;  // Required from Model
  virtual std::unique_ptr<Locations> locations() const = 0;

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
