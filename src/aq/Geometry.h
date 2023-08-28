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
#include "oops/util/Printable.h"

#include "aq/GeometryIterator.h"
#include "aq/interface.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace aq {

class GeometryIterator;

// -----------------------------------------------------------------------------
/// Geometry handles geometry for AQ model.

class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname() {return "aq::Geometry";}

  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
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
  bool levelsAreTopDown() const {return true;}

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
  atlas::functionspace::StructuredColumns functionSpaceNoHalo_;
  atlas::FieldSet extraFields_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_GEOMETRY_H_
