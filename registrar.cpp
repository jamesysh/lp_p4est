#include "registrar.h"
#include "geometry.h"
#include "geometry_pellet.h"
#include "state_pellet.h"
#include "boundary_pellet.h"
namespace{

    GeometryRegistrar<PelletLayer> g1("pelletlayer");
    StateRegistrar<PelletState> s1("pelletstate");
    BoundaryRegistrar<PelletInflowBoundary> b1("inflowboundary");

}
