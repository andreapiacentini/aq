/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_GEOMETRYAQ_H_
#define AQ_MODEL_GEOMETRYAQ_H_

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

#include "oops/aq/AqFortran.h"
#include "oops/aq/GeometryAQIterator.h"

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

class GeometryAQIterator;

// -----------------------------------------------------------------------------
/// GeometryAQ handles geometry for AQ model.

class GeometryAQ : public util::Printable,
                   private util::ObjectCounter<GeometryAQ> {
 public:
  typedef GeometryAqParameters Parameters_;

  static const std::string classname() {return "aq::GeometryAQ";}

  GeometryAQ(const GeometryAqParameters &, const eckit::mpi::Comm &);
  GeometryAQ(const GeometryAQ &);
  ~GeometryAQ();

  const F90geom & toFortran() const {return keyGeom_;}

  GeometryAQIterator begin() const;
  GeometryAQIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::Grid * atlasGrid() const {return atlasGrid_.get();}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}

  std::vector<size_t> variableSizes(const oops::Variables & vars) const;

 private:
  GeometryAQ & operator=(const GeometryAQ &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::StructuredGrid> atlasGrid_;
  std::unique_ptr<atlas::functionspace::StructuredColumns> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};
// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODEL_GEOMETRYAQ_H_
