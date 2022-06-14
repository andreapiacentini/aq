/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
  atlas::Grid grid() const {return grid_;}
  const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
  atlas::FunctionSpace & functionSpace() {return functionSpace_;}
  const atlas::FieldSet & extraFields() const {return extraFields_;}
  atlas::FieldSet & extraFields() {return extraFields_;}

  std::vector<size_t> variableSizes(const oops::Variables & vars) const;

  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  int halo_ = 0;
  const eckit::mpi::Comm & comm_;
  eckit::LocalConfiguration gridConfig_;
  atlas::StructuredGrid grid_;
  atlas::functionspace::StructuredColumns functionSpace_;
  atlas::functionspace::StructuredColumns functionSpaceSurf_;
  atlas::functionspace::StructuredColumns functionSpaceNoHalo_;
  atlas::FieldSet extraFields_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_GEOMETRY_H_
