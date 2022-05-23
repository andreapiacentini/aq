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

#ifndef AQ_OBSAUXCONTROL_H_
#define AQ_OBSAUXCONTROL_H_

#include <array>
#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/ObsAuxParameters.h"

namespace aq {
  class ObsAuxIncrement;
  class ObsSpace;

/// Class to handle observation bias parameters.

// -----------------------------------------------------------------------------

class ObsAuxControl : public util::Printable,
                private boost::noncopyable,
                private util::ObjectCounter<ObsAuxControl> {
 public:
  typedef ObsAuxControlParameters Parameters_;

  static const unsigned int ntypes = 4;
  static const std::string classname() {return "aq::ObsAuxControl";}

  ObsAuxControl(const ObsSpace &, const Parameters_ &);
  ObsAuxControl(const ObsAuxControl &, const bool);
  ~ObsAuxControl() {}

  ObsAuxControl & operator+=(const ObsAuxIncrement &);
  ObsAuxControl & operator=(const ObsAuxControl &);

/// I/O and diagnostics
  void read(const Parameters_ &) {}
  void write(const Parameters_ &) const {}
  double norm() const;

  const double & operator[](const unsigned int ii) const {return bias_[ii];}

  /// Other
  const oops::Variables & requiredVars() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}

  const double & insitu() const {return bias_[0];}

 private:
  void print(std::ostream &) const;
  std::array<double, ntypes> bias_;
  bool active_;
  const oops::Variables geovars_;
  const oops::Variables hdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSAUXCONTROL_H_
