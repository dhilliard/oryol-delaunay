#include "Mesh.h"
#include "Geo2D.h"
constexpr double EPSILON = 0.01;
constexpr double EPSILON_SQUARED = 0.0001;
using namespace Geo2D;

void Delaunay::Mesh::Setup(double width, double height)
{
	//Clear everything
	vertices.Clear();
	freeVertices.Clear();
	faces.Clear();
	freeFaces.Clear();
	//Infinite vertex

	//Other vertices in CW order
}
