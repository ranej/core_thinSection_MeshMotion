#include <PCU.h>
#include <lionPrint.h>
#include "phThinSectionStack.h"
#include "phOutput.h"
namespace ph {
  void getThinSectionStack(Output& o) {
    o.nTS = 0;
    o.nTSMeshVertices = 0;
    if(PCU_Comm_Self() == 0)
      lion_oprint(1,"warning! \'%s\' requires the Simmetrix SimAdvMeshing library\n",__func__);
    return;
  }
}
