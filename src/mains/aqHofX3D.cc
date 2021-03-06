/*
 * (C) Copyright 2018-2020 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/Traits.h"
#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::HofX3D<aq::Traits, aq::ObsTraits> hofx;
  return run.execute(hofx);
}
