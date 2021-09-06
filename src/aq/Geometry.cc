/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "aq/aq_geom_interface.h"
#include "aq/Geometry.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
Geometry::Geometry(const GeometryAqParameters & params,
                       const eckit::mpi::Comm & comm) : comm_(comm) {
  // Get geometry subconfiguration
  eckit::LocalConfiguration geomConfig(params.toConfiguration());
  geomConfig.set("type", "regional");

  // Setup regional grid
  atlasGrid_.reset(new atlas::StructuredGrid(geomConfig));

  // Setup partitioner
  atlas::grid::Partitioner partitioner("checkerboard");

  // Setup distribution
  atlas::grid::Distribution distribution(*atlasGrid_, partitioner);

  // Setup function space
  atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(*atlasGrid_, distribution,
  geomConfig));

  // Setup Fortran geometry
  aq_geom_setup_f90(keyGeom_, geomConfig,  &comm_, atlasGrid_->get(), atlasFunctionSpace_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  aq_geom_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  // Copy ATLAS grid
  atlasGrid_.reset(new atlas::StructuredGrid(*(other.atlasGrid_)));

  // Copy ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(
                            *(other.atlasFunctionSpace_)));

  // Copy Fortran geometry
  aq_geom_clone_f90(keyGeom_, other.keyGeom_);

  // Copy ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }
}
// -----------------------------------------------------------------------------
Geometry::~Geometry() {
  aq_geom_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
GeometryIterator Geometry::begin() const {
  return GeometryIterator(*this);
}
// -----------------------------------------------------------------------------
GeometryIterator Geometry::end() const {
  int nx = 0;
  int ny = 0;
  int nz;
  double deltax;
  double deltay;
  int mod_levels;
  char orientation[AQ_STRLEN];
  char domname[AQ_STRLEN];
  char model[AQ_STRLEN];
  aq_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay, mod_levels, orientation, domname, model);
  return GeometryIterator(*this, nx*ny+1);
}
// -------------------------------------------------------------------------------------------------
std::vector<double> Geometry::verticalCoord(std::string & vcUnits) const {
  // returns vertical coordinate in untis of vcUnits
  int nx = 0;
  int ny = 0;
  int nz;
  double deltax;
  double deltay;
  int mod_levels;
  char orientation[AQ_STRLEN];
  char domname[AQ_STRLEN];
  char model[AQ_STRLEN];
  aq_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay, mod_levels, orientation, domname, model);
  std::vector<double> vc(nz);
  if (vcUnits == "levels") {
    for (int i=0; i < nz; ++i) {vc[i]=i+1;}
  } else {
    std::stringstream errorMsg;
    errorMsg << "Uknown vertical coordinate unit " << vcUnits << std::endl;
    ABORT(errorMsg.str());
  }
  oops::Log::debug() << "AQ vert coord: " << vc << std::endl;
  return vc;
}
// -------------------------------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  // Note: in aq we always do trilinear interpolation, so GeoVaLs are always
  // size 1.
  std::vector<size_t> sizes(vars.size(), 1);
  return sizes;
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  int nx;
  int ny;
  int nz;
  double deltax;
  double deltay;
  int mod_levels;
  char orientation[AQ_STRLEN];
  char domname[AQ_STRLEN];
  char model[AQ_STRLEN];
  aq_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay, mod_levels, orientation, domname, model);
  os << "Geometry:" << std::endl;
  os << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
  os << "vert coord orientation = " << orientation << std::endl;
  os << "deltax = " << deltax << ", deltay = " << deltay << std::endl;
  os << "model = " << model << ", domain = " << domname << ", levels = " << mod_levels;
}
// -----------------------------------------------------------------------------
}  // namespace aq
