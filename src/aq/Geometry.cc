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

#include <algorithm>
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
  gridConfig_ = params.toConfiguration();
  gridConfig_.set("type", "regional");
  if (gridConfig_.has("halo")) {
    halo_ = gridConfig_.getInt("halo");
  }

  // Set default communicator
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // Setup regional grid
  grid_ = atlas::StructuredGrid(gridConfig_);

  // Setup partitioner
  atlas::grid::Partitioner partitioner = atlas::grid::Partitioner("checkerboard");

  // Setup function space
  functionSpace_ = atlas::functionspace::StructuredColumns(grid_, partitioner, gridConfig_);
  // Setup surf function space (always with halo for parallel interpolation)
  eckit::LocalConfiguration gridConfigSurf(params.toConfiguration());
  gridConfigSurf.set("type", "regional");
  gridConfigSurf.set("halo", 1);
  functionSpaceSurf_ = atlas::functionspace::StructuredColumns(
                               grid_, partitioner, gridConfigSurf));
  // Extra function space without halo (coincident with the previous if halo is zero
  if (halo_ > 0) {
    eckit::LocalConfiguration gridConfigNoHalo(params.toConfiguration());
    gridConfigNoHalo.set("type", "regional");
    gridConfigNoHalo.set("halo", 0);
    functionSpaceNoHalo_ = atlas::functionspace::StructuredColumns(
                                  grid_, partitioner, gridConfigNoHalo);
  }

  // Setup Fortran geometry
  aq_geom_setup_f90(keyGeom_, gridConfig_,  &comm_, grid_.get(),
                    functionSpace_.get(), functionSpaceSurf_.get());

  // Fill ATLAS fieldset
  extraFields_ = atlas::FieldSet();
  aq_geom_fill_extra_fields_f90(keyGeom_, extraFields_.get());
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  // Copy ATLAS grid
  gridConfig_ = other.gridConfig_;
  grid_ = atlas::StructuredGrid(gridConfig_);

  // Copy ATLAS function space
  functionSpace_ = atlas::functionspace::StructuredColumns(other.functionSpace_);
  halo_ = other.halo_;
  if (halo_ > 0) {
    functionSpaceNoHalo_ = atlas::functionspace::StructuredColumns(other.functionSpaceNoHalo_);
  }

  // Copy Fortran geometry
  aq_geom_clone_f90(keyGeom_, other.keyGeom_);

  // Copy ATLAS fieldset
  extraFields_ = atlas::FieldSet();
  for (auto & field : other.extraFields_) {
    extraFields_.add(field);
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
void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
                      const bool halo) const {
  const atlas::functionspace::StructuredColumns * fspace;
  if (halo_ > 0) {
    if (halo) {
      fspace = &(functionSpace_);
    } else {
      fspace = &(functionSpaceNoHalo_);
    }
    const auto lonlat = atlas::array::make_view<double, 2>(fspace->lonlat());
    const size_t npts = fspace->size();
    lats.resize(npts);
    lons.resize(npts);
    for (size_t jj = 0; jj < npts; ++jj) {
      lats[jj] = lonlat(jj, 1);
      lons[jj] = lonlat(jj, 0);
    }
  } else {
    /* If the fields have no halo, the locations have to be only on proc 0, therefore latlon
       is provided only there */
    if (comm_.rank() == 0) {
      const size_t npts = grid_.size();
      lats.resize(npts);
      lons.resize(npts);
      atlas::Grid::PointLonLat lonlat;
      size_t jj = 0;
      for (atlas::idx_t jy = 0; jy < grid_.ny(); ++jy) {
        for (atlas::idx_t jx = 0; jx < grid_.nx(jy); ++jx) {
          lonlat = grid_.lonlat(jx, jy);
          lats[jj] = lonlat[1];
          lons[jj] = lonlat[0];
          jj++;
        }
      }
    } else {
      const size_t npts = 0;
      lats.resize(npts);
      lons.resize(npts);
    }
  }
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
