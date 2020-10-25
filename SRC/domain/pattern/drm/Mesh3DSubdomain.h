#ifndef Mesh3DSubdomain_h
#define Mesh3DSubdomain_h

<<<<<<< HEAD
#include <stdlib.h>
#include <iostream>
#include <OPS_Globals.h>
#include <StandardStream.h>
#include <Timer.h>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <NDMaterial.h>

#include <ArrayOfTaggedObjects.h>

// includes for the domain classes
#include <Domain.h>
#include <ShadowSubdomain.h>
#include <ActorSubdomain.h>
#include <Node.h>
#include <Element.h>
#include <ZeroLength.h>
#include <ElasticMaterial.h>
#include <ElasticIsotropicThreeDimensional.h>

#include <SP_Constraint.h>

#include <SingleDomAllSP_Iter.h>
#include <NodeIter.h>

#include <Vertex.h>
#include <VertexIter.h>
#include <Graph.h>

#include <Vector.h>
#include <Matrix.h>
#include "GeometricBrickDecorator.h"

#include <map>
#include <set>
#include <math.h>

using namespace std;
=======
#include <Element.h>
#include <map>
#include <set>

class Domain;

>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

class Mesh3DSubdomain {

  public :
    
    //general constructor
    Mesh3DSubdomain(Domain * inpDomain) ;


    Mesh3DSubdomain() ;


    //destructor
    virtual ~Mesh3DSubdomain();
    
    int lastEle();

    int lastNode();
    
    int lastStartNode();
  
    void allocate_e_Nodes(double xMin, double xMax, 
			  double yMin, double yMax,
			  double zMin, double zMax,
			  std::map<int,int>& eNodes);

    void allocateBoundaryLayerElements(double xMin, double xMax,
				       double yMin, double yMax,
				       double zMin, double zMax,
				       std::map<int,Element*>& elements,
				       std::map<int,Vector*>& storage,
				       std::map<int,int>& storage2);

    private :
     
    Domain * myDomain;
    int myLastStartNodeTag;
    int myLastEle;
    int myLastNode;
<<<<<<< HEAD



=======
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
};

#endif






