/*
 * Copyright 2016 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <string>
#include <pcu_util.h>
#include <cstdlib>

#include <reel.h>
#include <pcu_io.h>
#include "apfFile.h"
#include "apf.h"
#include "apfShape.h"
#include "apfTagData.h"
#include "apfNumbering.h"

namespace apf {

static void save_string(pcu_file* file, const char* s) {
  pcu_write_string(file, s);
}

static std::string restore_string(pcu_file* file) {
  char* s;
  pcu_read_string(file, &s);
  std::string s_cpp(s);
  free(s);
  return s_cpp;
}

static void save_int(pcu_file* file, int x) {
  unsigned tmp = static_cast<unsigned>(x);
  PCU_WRITE_UNSIGNED(file, tmp);
}

static int restore_int(pcu_file* file) {
  unsigned tmp;
  PCU_READ_UNSIGNED(file, tmp);
  return static_cast<int>(tmp);
}

static void save_field_meta(pcu_file* file, apf::Field* field) {
  save_string(file, getName(field));
  save_int(file, getValueType(field));
  save_int(file, countComponents(field));
  save_string(file, getShape(field)->getName());
}

static void restore_field_meta(pcu_file* file, apf::Mesh* mesh) {
  std::string field_name = restore_string(file);
  int value_type = restore_int(file);
  int ncomps = restore_int(file);
  std::string shape_name = restore_string(file);
  apf::FieldShape* shape = getShapeByName(shape_name.c_str());
  if (shape == 0)
    reel_fail("field shape \"%s\" could not be found\n", shape_name.c_str());
  makeField(mesh, field_name.c_str(), value_type, ncomps,
      shape, new TagDataOf<double>);
}

static void save_numbering_meta(pcu_file* file, apf::Numbering* numbering) {
  save_string(file, getName(numbering));
  save_int(file, countComponents(numbering));
  save_string(file, getShape(numbering)->getName());
}

static void restore_numbering_meta(pcu_file* file, apf::Mesh* mesh) {
  std::string numbering_name = restore_string(file);
  int ncomps = restore_int(file);
  std::string shape_name = restore_string(file);
  apf::FieldShape* shape = getShapeByName(shape_name.c_str());
  if (shape == 0)
    reel_fail("numbering shape \"%s\" could not be found\n", shape_name.c_str());
  createNumbering(mesh, numbering_name.c_str(), shape, ncomps);
}

static int latest_version_number = 3;

void save_meta(pcu_file* file, apf::Mesh* mesh) {
  save_string(file, mesh->getShape()->getName());
/* the first version of this system made the mistake of
 * not including a version number.
 * to introduce a backwards-compatible version number,
 * what we do here is to write a negative version
 * number where the old code would have expected the
 * number of fields.
 */
  save_int(file, -latest_version_number);
  save_int(file, mesh->countFields());
  for (int i = 0; i < mesh->countFields(); ++i) {
    save_field_meta(file, mesh->getField(i));
  }
  save_int(file, mesh->countNumberings());
  for (int i = 0; i < mesh->countNumberings(); ++i) {
    save_numbering_meta(file, mesh->getNumbering(i));
  }
}

void restore_meta(pcu_file* file, apf::Mesh* mesh) {
  std::string shape_name = restore_string(file);
  apf::FieldShape* shape = getShapeByName(shape_name.c_str());
  PCU_ALWAYS_ASSERT(shape != 0);
  if (shape != mesh->getShape()) mesh->changeShape(shape, false);
  int nfields_or_version = restore_int(file);
  int nfields;
  int version;
  if (nfields_or_version >= 0) {
    nfields = nfields_or_version;
    version = 1;
  } else {
    nfields = restore_int(file);
    version = -nfields_or_version;
  }
  PCU_ALWAYS_ASSERT(version <= latest_version_number);
  PCU_ALWAYS_ASSERT(nfields >= 0);
  PCU_ALWAYS_ASSERT(nfields < 256);
  for (int i = 0; i < nfields; ++i) {
    restore_field_meta(file, mesh);
  }
  if (version >= 3) {
    int nnumberings = restore_int(file);
    PCU_ALWAYS_ASSERT(nnumberings >= 0);
    PCU_ALWAYS_ASSERT(nnumberings < 256);
    for (int i = 0; i < nnumberings; ++i) {
      restore_numbering_meta(file, mesh);
    }
  }
}

} // end namespace apf
