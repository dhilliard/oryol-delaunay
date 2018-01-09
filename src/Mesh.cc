#include "Mesh.h"
#include "Geo2D.h"
constexpr double EPSILON = 0.01;
constexpr double EPSILON_SQUARED = 0.0001;
using namespace Geo2D;

int getDirection(double x1, double y1, double x2, double y2, double x3, double y3) {
	// dot product with the orthogonal vector pointing left vector of eUp:
	double dot = (x3 - x1) * (y2 - y1) + (y3 - y1) * (-x2 + x1);
	return (dot == 0) ? 0 : ((dot > 0) ? 1 : -1);
}

int getDirectionSafe(double x1, double y1, double x2, double y2, double x3, double y3) {
	// dot product with the orthogonal vector pointing left vector of eUp:
	double dot = (x3 - x1) * (y2 - y1) + (y3 - y1) * (-x2 + x1);

	// check sign
	if (dot == 0)
	{
		return 0;
	}
	else if (dot > 0)
	{
		if (DistanceSquaredPointToLine(x3, y3, x1, y1, x2, y2) <= EPSILON_SQUARED)
			return 0;
		else
			return 1;
	}
	else
	{
		if (DistanceSquaredPointToLine(x3, y3, x1, y1, x2, y2) <= EPSILON_SQUARED)
			return 0;
		else
			return -1;
	}
}

