/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crvSnap.h"

#include <gmi.h>
#include <maSnap.h>

namespace crv {

/** \brief returns true if there is a geometric
      degeneracy in that direction at that parameter
     \param p is the parameters of the original coordinate
     \details degeneracy is when at a point, changing
     one parameter does not change the point location.
     This is tested numerically, by moving a point along
     the other parameter and checking if its location has
     changed. If a surface has degenerate coordinates,
     every point on it needs to be checked for degeneracy
     for edge splitting, or similar processes. */
  // TODO: Remove this if not used anywhere else since a better
  // version is now implemented in ma/maSnap.cc
static bool checkIsDegenerate(apf::Mesh* m, apf::ModelEntity* g,
    apf::Vector3 const& p, int axis)
{
  gmi_ent* e = (gmi_ent*)g;
  gmi_model* gm = m->getModel();
  int md = gmi_dim(gm, e);
  if (md != 2)
    return 0;
  apf::Vector3 x,y,q;
  gmi_eval(gm, e, &p[0], &x[0]); // original point
  int other = axis ? 0 : 1; // move along other direction
  double range[2]; // compute the range to determine how to perturb
  gmi_range(gm, e, other, range);
  if (range[0] > range[1])
    std::swap(range[0],range[1]);
  // this just guarantees we are checking a point sufficiently far
  if (range[1] - p[other] > p[other]-range[0])
    q[other] = 0.25*(range[1]-range[0])+p[other];
  else
    q[other] = p[other]-0.25*(range[1]-range[0]);
  q[axis] = p[axis];
  gmi_eval(gm, e, &q[0], &y[0]); // point along that axis
  return ((x-y).getLength() < 1e-13);
}

static void interpolateParametricCoordinates(
    apf::Mesh* m,
    apf::ModelEntity* g,
    double t,
    apf::Vector3 const& a,
    apf::Vector3 const& b,
    apf::Vector3& p)
{
  // the corresponding function in ma takes care of degenerate
  // cases, now. So it is enough to just call that!
  ma::interpolateParametricCoordinates(m, g, t, a, b, p);
}

static void transferParametricBetween(
    apf::Mesh* m,
    apf::ModelEntity* g,
    apf::MeshEntity* v[2],
    double t,
    apf::Vector3& p)
{
  apf::Vector3 ep[2];
  for (int i=0; i < 2; ++i)
    m->getParamOn(g,v[i],ep[i]);
  crv::interpolateParametricCoordinates(m,g,t,ep[0],ep[1],p);
}

void transferParametricOnEdgeSplit(
    apf::Mesh* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p)
{
  apf::ModelEntity* g = m->toModel(e);
  apf::MeshEntity* ev[2];
  m->getDownward(e,0,ev);
  crv::transferParametricBetween(m, g, ev, t, p);
}

/* see bezier.tex in SCOREC/docs repo for more info on this,
   but its straightforward, two line splits */
void transferParametricOnTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p)
{
  apf::ModelEntity* g = m->toModel(e);
  apf::MeshEntity* ev[3];
  m->getDownward(e,0,ev); // pick two points, split on edge
  // adjust so the degenerate point is in ev[2]
  m->getParamOn(g,ev[0],p);
  if(checkIsDegenerate(m,g,p,0) || checkIsDegenerate(m,g,p,1)){
    apf::MeshEntity* v = ev[0];
    ev[0] = ev[1]; ev[1] = ev[2];ev[2] = v;
    t[2] = 1-t[0]-t[1]; t[0] = t[1]; t[1] = t[2];
  } else {
    m->getParamOn(g,ev[1],p);
    if(checkIsDegenerate(m,g,p,0) || checkIsDegenerate(m,g,p,1)){
      apf::MeshEntity* v = ev[1];
      ev[1] = ev[0]; ev[0] = ev[2]; ev[2] = v;
      t[2] = 1-t[0]-t[1]; t[1] = t[0];t[0] = t[2];
    }
  }
  apf::Vector3 p1,p2;
  m->getParamOn(g,ev[2],p2);
  // two linear splits
  crv::transferParametricBetween(m, g, ev, t[0]/(1.-t[1]), p1);
  crv::interpolateParametricCoordinates(m,g,t[1],p1,p2,p);
}

void transferParametricOnGeometricEdgeSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p)
{
  apf::ModelEntity* g = m->toModel(e);
  apf::MeshEntity* ev[2];
  m->getDownward(e,0,ev);
  apf::Vector3 p0,p1,cpt;
  m->getPoint(ev[0],0,p0);
  m->getPoint(ev[1],0,p1);
  apf::Vector3 pt = p0*(1.-t)+p1*t;
  m->getClosestPoint(g,pt,cpt,p);
}

void transferParametricOnGeometricTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p)
{
  apf::ModelEntity* g = m->toModel(e);
  apf::MeshEntity* ev[3];
  m->getDownward(e,0,ev);
  // split in physical space, project
  apf::Vector3 p0,p1,p2,cpt;
  m->getPoint(ev[0],0,p0);
  m->getPoint(ev[1],0,p1);
  m->getPoint(ev[2],0,p2);
  apf::Vector3 pt = p0*(1.-t[0]-t[1])+p1*t[0]+p2*t[1];
  m->getClosestPoint(g,pt,cpt,p);
}

} // namespace crv

