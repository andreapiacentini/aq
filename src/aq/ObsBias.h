/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_OBSBIAS_H_
#define AQ_MODEL_OBSBIAS_H_

#include <array>
#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "model/ObsBiasParameters.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace aq {
  class ObsBiasIncrement;
  class ObsSpaceAQ;

/// Class to handle observation bias parameters.

// -----------------------------------------------------------------------------

class ObsBias : public util::Printable,
                private boost::noncopyable,
                private util::ObjectCounter<ObsBias> {
 public:
  typedef ObsBiasParameters Parameters_;

  static const unsigned int ntypes = 4;
  static const std::string classname() {return "aq::ObsBias";}

  ObsBias(const ObsSpaceAQ &, const Parameters_ &);
  ObsBias(const ObsBias &, const bool);
  ~ObsBias() {}

  ObsBias & operator+=(const ObsBiasIncrement &);
  ObsBias & operator=(const ObsBias &);

/// I/O and diagnostics
  void read(const Parameters_ &) {}
  void write(const Parameters_ &) const {}
  double norm() const;

  const double & operator[](const unsigned int ii) const {return bias_[ii];}

  /// Other
  const oops::Variables & requiredVars() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}

  const double & stream() const {return bias_[0];}
  const double & wind() const {return bias_[1];}
  const double & wspd() const {return bias_[3];}

 private:
  void print(std::ostream &) const;
  std::array<double, ntypes> bias_;
  bool active_;
  const oops::Variables geovars_;
  const oops::Variables hdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODEL_OBSBIAS_H_
