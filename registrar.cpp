#include "registrar.h"
#include "geometry.h"
#include "geometry_pellet.h"
#include "state_pellet.h"
#include "geometry_disk.h"
#include "boundary_pellet.h"
#include "state_gresho.h"

namespace{

    GeometryRegistrar<PelletLayer> g1("pelletlayer");
    StateRegistrar<PelletState> s1("pelletstate");
    BoundaryRegistrar<PelletInflowBoundary> b1("inflowboundary");
    GeometryRegistrar<Disk> g2("disk");
    StateRegistrar<Gresho2DState> s2("gresho2dstate");

}
