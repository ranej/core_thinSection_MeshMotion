#include <apfIntegrate.h>
#include <apfMesh.h>
#include <pcu_util.h>
#include <cstdlib>
#include <iostream>

static void testType(int type, double expectedSum)
{
  apf::EntityIntegration const* eg = apf::getIntegration(type);
  PCU_ALWAYS_ASSERT(eg);
  for (int i = 0; i < eg->countIntegrations(); ++i) {
    apf::Integration const* g = eg->getIntegration(i);
    PCU_ALWAYS_ASSERT(g);
    double sum = 0;
    for (int j = 0; j < g->countPoints(); ++j) {
      apf::IntegrationPoint const* p = g->getPoint(j);
      PCU_ALWAYS_ASSERT(p);
      sum += p->weight;
    }
    if (fabs(sum - expectedSum) > 1e-10) {
      std::cerr << apf::Mesh::typeName[type] << " rule #" << i
        << " (order " << g->getAccuracy() << ") is wrong:\n";
      std::cerr << "adds to " << sum << ", should be " << expectedSum << '\n';
      abort();
    }
  }
}

int main()
{
  testType(apf::Mesh::EDGE,     2.0);
  testType(apf::Mesh::TRIANGLE, 1.0/2.0);
  testType(apf::Mesh::QUAD,     4.0);
  testType(apf::Mesh::TET,      1.0/6.0);
  testType(apf::Mesh::HEX,      8.0);
  return 0;
}
