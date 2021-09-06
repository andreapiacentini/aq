/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODELAUXCONTROL_H_
#define AQ_MODELAUXCONTROL_H_

#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace aq {
  class Geometry;
  class ModelAuxIncrement;

/// Model error for the AQ model.
/*!
 * This class is used to manipulate parameters of the model that
 * can be estimated in the assimilation. This includes model bias for
 * example but could be used for other parameters to be estimated.
 * This is sometimes referred to as augmented state or augmented
 * control variable in the litterature.
 * The augmented state is understood here as an augmented 4D state.
 */

// -----------------------------------------------------------------------------

class ModelAuxControl : public util::Printable,
                  private boost::noncopyable,
                  private util::ObjectCounter<ModelAuxControl> {
 public:
  static const std::string classname() {return "aq::ModelAuxControl";}

  ModelAuxControl(const Geometry &, const eckit::Configuration &) {}
  ModelAuxControl(const Geometry &, const ModelAuxControl &) {}
  ModelAuxControl(const ModelAuxControl &, const bool) {}
  ~ModelAuxControl() {}

  ModelAuxControl & operator+=(const ModelAuxIncrement &) {return *this;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODELAUXCONTROL_H_
