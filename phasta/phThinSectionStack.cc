#include <PCU.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include "phOutput.h"
#include "phThinSectionStack.h"
#include "phLinks.h"
#include "phAdjacent.h"
#include "apfSIM.h"
#include "gmi_sim.h"
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimAdvMeshing.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

namespace ph {
void getThinSectionStack(Output& o)
{
  o.nThinSectionStacks = 0;
  o.nThinSectionStackMeshVertices = 0;

  Input& in = *o.in;
  if (in.simmetrixMesh == 1) {
    if (in.writeSimLog)
      Sim_logOn("getThinStack.log");
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    // get simmetrix mesh
    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(o.mesh);
    pParMesh parMesh = apf_msim->getMesh();
    pMesh mesh = PM_mesh(parMesh,0);

    // get simmetrix model
    gmi_model* gmiModel = apf_msim->getModel();
    pGModel model = gmi_export_sim(gmiModel);

//  Algorithm: Get growth curve info

//    pGEntity gEntity;
    pGFace gFace;
    pVertex vertex;
    pFace face;
//    pGFace destFace;
    pRegion region;
    pVertex vert_temp;

    pPList gEdges = PList_new();
    pPList gVertices = PList_new();
    pPList gSrcFaces = PList_new();
    pPList regions_stack = PList_new();
    pPList faces_stack = PList_new();
//    pPList gDestFaces =PList_new();

//    pFace destFace_temp;
//    int list_size = 0;
    int iThinEntity = 0;
    int count = 0;
//    pPList Fvertices = PList_new();
//    int vtype = 0;
//    int stack_size;

    GFIter gFIter = GM_faceIter(model);
    while((gFace = GFIter_next(gFIter))){	//Loop over gFaces
//    IF gFace has extrusion/thin stack mesh attribute
//	cout << "Working on GFace: " << GEN_tag((pGEntity)gFace) << endl;
//        bool isThinStackFace = false;
//	bool isTest = true;
      FIter fIter = M_classifiedFaceIter(mesh, gFace, 0);	
      while((face = FIter_next(fIter))){	//Loop over mFaces on gFace
//	cout << "Working on mface: " << EN_id(face) << endl;
        if(EN_isExtrusionEntity(face) == 1){
	  count = count+1;
//	  cout << "potential extrusion entity" << endl;
	  for(int j = 0; j < 2; j++){
            if(EN_isExtrusionEntity(F_region(face, j)) == 1){
//	      cout << "Extrusion Entity Found" <<endl;
              region = F_region(face, j);
//	      pPList Fvertices = PList_new();
//	      Fvertices = F_vertices(face, 1);
              iThinEntity = Extrusion_3DRegionsAndLayerFaces(region, regions_stack, faces_stack, 0);
//	      stack_size = PList_size(faces_stack);
	      if(iThinEntity == 1){
		PList_appUnique(gSrcFaces, F_whatIn((pFace)PList_item(faces_stack, 0)));
		break;
	      }
            }  
	  }  
	}
      }
      
    FIter_delete(fIter);
    }


    GFIter_delete(gFIter);
//    FIter_delete(fIter);

//    std::cout<< PList_size(gSrcFaces) << " : gSrcFaces list size \n" ;

//    for(int i = 0; i < PList_size(gSrcFaces); i++ ){
//      cout<< GEN_tag((pGEntity)PList_item(gSrcFaces,i)) << ": source face tag" <<endl;
//    }

// debug part 2
//EVERYTHING UPTO THIS POINT IS GOOD//
//
    //  get seeds of all growth curves
    pPList allSeeds = PList_new();
    pPList gFaces = PList_new();
    pPList seeds = PList_new();
    pPList base_vertices = PList_new();
    pPList fvertices = PList_new();
    pPList fvertices2 = PList_new();
    pPList allStackVertices = PList_new();
    pPList StackVertices = PList_new();
    pPList vvertices = PList_new();
//    pPList iTSnv_list = PList_new();

//   std::list<int> iTSnv_list;
    int iTSnv_listsize = 0;
//      pRegion region;

//   cout << PList_size(gSrcFaces) << ":List_size" << endl;

//  FOR each gEntity in gEntities
    for(int i = 0; i < 1; i++){
      if(PList_size(gSrcFaces)<1)
	break;

      gFace = (pGFace)PList_item(gSrcFaces,i);
//      cout << "Working on GFace : " << GEN_tag(gFace)<< endl;

//    Get mesh faces classified on gEntity
      FIter fIter2 = M_classifiedFaceIter(mesh, gFace, 0);
//      cout << "Looping over mesh faces" <<endl;

//    FOR each mesh face
      while((face = FIter_next(fIter2))){
//      Create an empty list (seeds) for storing potential seed edges of vertex
        PList_clear(seeds);
        PList_clear(regions_stack);
        PList_clear(faces_stack);

        for(int j = 0; j < 2; j++){
          if(EN_isExtrusionEntity(F_region(face, j)) == 1){
//	    cout << "Extrusion Entity found" << endl;
            region = F_region(face, j);
            PList_appUnique(seeds,region);                                          //add region to seeds List

            //Get region stack and faces_stack
            iThinEntity = Extrusion_3DRegionsAndLayerFaces(region, regions_stack, faces_stack, 0);
            PList_clear(fvertices);
	    fvertices = F_vertices((pFace)PList_item(faces_stack,0),1);                                     //get mesh vertices for base face

            for(int k = 0; k < PList_size(fvertices); k++){                                                         //Loop over vertices of base face
	      if(PList_contains(base_vertices, PList_item(fvertices,k))==0){                //Check if vertex is already in list. If not
//		cout << "Unique base_vertex found: " << EN_id((pVertex)PList_item(fvertices,k))<<endl;	
                PList_appUnique(base_vertices, PList_item(fvertices,k));                //add vertex to list
//                cout << "base_vertices_size" << PList_size(base_vertices) << endl;
                PList_clear(StackVertices);
		iTSnv_listsize = PList_size(faces_stack);
		PList_append(StackVertices, PList_item(fvertices,k));
		vert_temp = (pVertex)PList_item(fvertices,k);		
		double xyz[3];
		V_coord(vert_temp,xyz);
//		cout << "base coords: " <<  xyz[0] << " " << xyz[1] << " " << xyz[2] <<endl;

                for(int faces_stack_iter = 0; faces_stack_iter < PList_size(faces_stack)-1; faces_stack_iter++){                                 //Loop over faces in faces_stack
		  PList_clear(vvertices);
		  vvertices =  V_vertices(vert_temp);	
                  for(int vvertices_iter = 0; vvertices_iter < PList_size(vvertices); vvertices_iter++){                                 //Loop over faces in faces_stack
                    PList_clear(fvertices2);
                    fvertices2 = F_vertices((pFace)PList_item(faces_stack,faces_stack_iter+1),1);   //Get Stack vertices
		    for (int temp=0; temp<PList_size(fvertices2) ; temp++){
//			cout<< "face vertex ids:"<< EN_id((pVertex)PList_item(fvertices2,temp)) <<endl;
		    }
//		    cout<< "Vertex id in consideration:" << EN_id((pVertex)PList_item(vvertices,vvertices_iter)) <<endl;
                    if(PList_contains(fvertices2, PList_item(vvertices,vvertices_iter))==1){
//		      cout<< "Stack Vertex found" << endl;
		      vert_temp = (pVertex)PList_item(vvertices,vvertices_iter);
		      PList_append(StackVertices, PList_item(vvertices,vvertices_iter));
		      V_coord((pVertex)PList_item(vvertices,vvertices_iter),xyz);
//		      cout << "Stack Verticex id: "<< EN_id((pVertex)PList_item(vvertices,vvertices_iter)) << endl;
//		      cout << "coords: " <<  xyz[0] << " " << xyz[1] << " " << xyz[2] <<endl;
		  }
          //For face kk in connectivity stack, add vertex k to allStackVertices
                }
              }
              PList_appPList(allStackVertices, StackVertices);
            }
          }
        }
      }
     }
    FIter_delete(fIter2);
    }

    
// cout << "List size iTSnv:" << iTSnv_listsize <<endl;


  int nTS = PList_size(base_vertices);
  o.nThinSectionStacks = nTS;
	
  int nTSMeshVertices = PList_size(allStackVertices);						
  o.nThinSectionStackMeshVertices = nTSMeshVertices;

  o.arrays.iTSnv = new int[nTS];

	for(int iTSnv_iter = 0; iTSnv_iter < nTS; iTSnv_iter++){
		 o.arrays.iTSnv[iTSnv_iter] = iTSnv_listsize;
	}

//  count = 0;
//  for(int i : iTSnv_list){
//    o.arrays.iTSnv[count] = i;
//    count = count+1;	
//  }

//  iTSnv_list.clear();
  o.arrays.iTSlv = new apf::MeshEntity*[nTSMeshVertices];

  for(int i = 0; i < PList_size(allStackVertices); i++){
    vertex = (pVertex)PList_item(allStackVertices,i);

    apf::MeshEntity* me = reinterpret_cast<apf::MeshEntity*> (vertex);
    o.arrays.iTSlv[i] = me;
  }
  
  cout << "phThinSection" << endl;
  lion_oprint(1,"%s: rank %d, nTS, nv: %d, %d\n", __func__, PCU_Comm_Self(), nTS, nTSMeshVertices);

  cout << "Debug 1" << endl;
    PCU_Add_Ints(&nTS,sizeof(nTS));
  cout << "Debug 2" << endl;
    PCU_Add_Ints(&nTSMeshVertices,sizeof(nTSMeshVertices));

    if(PCU_Comm_Self() == 0)
      lion_oprint(1,"%s: total nTS, nv: %d, %d\n", __func__, nTS, nTSMeshVertices);
      
	cout << "after PCU_COMM_self" << endl;
	//Delete all lists
	PList_delete(gEdges);
	PList_delete(gVertices);
	PList_delete(allSeeds);
	PList_delete(gFaces);
	PList_delete(seeds);
	PList_delete(regions_stack);
	PList_delete(faces_stack);
	PList_delete(base_vertices);
	PList_delete(fvertices);
	PList_delete(fvertices2);
	PList_delete(allStackVertices);
	PList_delete(StackVertices);
	PList_delete(gSrcFaces);
	PList_delete(vvertices);
	
    cout << "done" <<endl;

    //clean up utility
    Progress_delete(progress);
    if (in.writeSimLog)
      Sim_logOff();
  }
  else {
    if(PCU_Comm_Self() == 0)
      lion_oprint(1,"%s: warning! not implemented for MDS mesh\n",__func__);
  }
  return;
}
} //end namespace ph
