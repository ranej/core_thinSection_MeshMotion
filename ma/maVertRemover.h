/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_VERTREMOVER_H
#define MA_VERTREMOVER_H

#include "maCollapse.h"

namespace ma {

/* tries to remove a vertex by collapsing
   one of the adjacent edges */

class VertRemover
{
  public:
/* Init instead of constructor because C++ doesn't allow
   creation of arrays of objects that have no default constructors */
    void Init(Adapt* a);
    void setVert(Entity* v);
    Entity* getVert() {return vert;}
    apf::Up& getEdges() {return edges;}
    void findEdges();
    bool tryToCollapse(Entity* e);
    bool run();
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* vert;
    apf::Up edges;
    Collapse collapse;
};

}

#endif
