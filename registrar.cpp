#include "registrar.h"
namespace{

    GeometryRegistrar<PelletLayer> g1("pelletlayer");
    StateRegistrar<PelletState> s1("pelletstate");
    BoundaryRegistrar<PelletInflowBoundary> b1("inflowboundary");
    GeometryRegistrar<MultiPelletLayer> g2("multipellet");
}
