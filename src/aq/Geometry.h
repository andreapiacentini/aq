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

#ifndef AQ_GEOMETRY_H_
#define AQ_GEOMETRY_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "aq/GeometryIterator.h"
#include "aq/interface.h"

namespace oops {
  class Variables;
}

namespace aq {

class GeometryAqParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryAqParameters, Parameters)

 public:
  /// Domain size
  oops::RequiredParameter<int> nx{"nx", this};
  oops::RequiredParameter<int> ny{"ny", this};
  oops::RequiredParameter<int> levels{"levels", this};
  oops::RequiredParameter<double> dx{"dx", this};
  oops::RequiredParameter<double> dy{"dy", this};
  oops::RequiredParameter<double> xmin{"xmin", this};
  oops::RequiredParameter<double> ymin{"ymin", this};
  oops::RequiredParameter<std::string> domname{"domname", this};
  oops::Parameter<std::string> orientation{"orientation", "up", this};
  oops::Parameter<std::string> model{"model", "MOCAGE", this};
  oops::Parameter<int> halo{"halo", 0, this};
  oops::Parameter<int> mod_levels{"model levels", -1, this};
};

class GeometryIterator;

// -----------------------------------------------------------------------------
/// Geometry handles geometry for AQ model.

class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
 public:
  typedef GeometryAqParameters Parameters_;

  static const std::string classname() {return "aq::Geometry";}

  Geometry(const GeometryAqParameters &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  ~Geometry();

  const F90geom & toFortran() const {return keyGeom_;}

  GeometryIterator begin() const;
  GeometryIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::Grid * atlasGrid() const {return atlasGrid_.get();}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}

  std::vector<size_t> variableSizes(const oops::Variables & vars) const;

  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  int halo_ = 0;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::StructuredGrid> atlasGrid_;
  std::unique_ptr<atlas::functionspace::StructuredColumns> atlasFunctionSpace_;
  std::unique_ptr<atlas::functionspace::StructuredColumns> atlasFunctionSpaceNoHalo_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_GEOMETRY_H_
