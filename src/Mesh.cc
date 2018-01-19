#include "Mesh.h"
#include "Geo2D.h"
#include <cmath>
#include <limits>
#include "Core/Containers/Set.h"
#include "glm/geometric.hpp"

constexpr double EPSILON = 0.01;
constexpr double EPSILON_SQUARED = 0.0001;
using namespace Geo2D;
using namespace Delaunay;

int rand_range(int min, int max) {
	return min + (float)(std::rand() / RAND_MAX)*(max - min);
}
//TODO: Investigate a better approach for initialising everything.
void Delaunay::Mesh::Setup(double width, double height)
{
	//Clear everything
	vertices.Clear();
	faces.Clear();
	//Infinite vertex
	Vertex & vInf = requestVertex(width*0.5, height * 0.5, false);
	//Other vertices in CCW order
	Vertex & vBL = requestVertex(0, 0);
	Vertex & vBR = requestVertex(width, 0);
	Vertex & vTR = requestVertex(width, height);
	Vertex & vTL = requestVertex(0, height);
	vBL.constraintCount = 4;
	vBR.constraintCount = 4;
	vTR.constraintCount = 4;
	vTL.constraintCount = 4;
	
	//Request faces
	Face & fTL_BL_BR = requestFace(); //eBR_TL & eTL_BL & eBL_BR [0]
	Face & fBR_TR_TL = requestFace(); //eTL_BR & eBR_TR & eTR_TL [1]
	Face & fTL_TR_vInf = requestFace(); //eInf_TL & eTL_TR & eTR_Inf [2]
	Face & fBL_TL_vInf = requestFace(); //eInf_BL & eBL_TL & eTL_Inf [3]
	Face & fBR_BL_vInf = requestFace(); //eInf_BR & eBR_BL & eBL_Inf [4]
	Face & fTR_BR_vInf = requestFace(); //eInf_TR & eTR_BR & eBR_Inf [5]
	//Setup vertex data
	fTL_BL_BR.edges[0].destinationVertex = indexFor(vTL);
	fTL_BL_BR.edges[1].destinationVertex = indexFor(vBL);
	fTL_BL_BR.edges[2].destinationVertex = indexFor(vBR);
	fTL_BL_BR.flags = 0b1100;
	
	fBR_TR_TL.edges[0].destinationVertex = indexFor(vBR);
	fBR_TR_TL.edges[1].destinationVertex = indexFor(vTR);
	fBR_TR_TL.edges[2].destinationVertex = indexFor(vTL);
	fBR_TR_TL.flags = 0b1100;
	
	fTL_TR_vInf.edges[0].destinationVertex = indexFor(vTL);
	fTL_TR_vInf.edges[1].destinationVertex = indexFor(vTR);
	fTL_TR_vInf.edges[2].destinationVertex = indexFor(vInf);
	fTL_TR_vInf.flags = 0b0101;

	fBL_TL_vInf.edges[0].destinationVertex = indexFor(vBL);
	fBL_TL_vInf.edges[1].destinationVertex = indexFor(vTL);
	fBL_TL_vInf.edges[2].destinationVertex = indexFor(vInf);
	fBL_TL_vInf.flags = 0b0101;

	fBR_BL_vInf.edges[0].destinationVertex = indexFor(vBR);
	fBR_BL_vInf.edges[1].destinationVertex = indexFor(vBL);
	fBR_BL_vInf.edges[2].destinationVertex = indexFor(vInf);
	fBR_BL_vInf.flags = 0b0101;

	fTR_BR_vInf.edges[0].destinationVertex = indexFor(vTR);
	fTR_BR_vInf.edges[1].destinationVertex = indexFor(vBR);
	fTR_BR_vInf.edges[2].destinationVertex = indexFor(vInf);
	fTR_BR_vInf.flags = 0b0101;
	//Setup shared edges
	//Weld bottom of pyramid together
	fTL_BL_BR.edges[0].oppositeHalfEdge = indexFor(fBR_TR_TL, 0); // eBR_TL.opposite = eTL_BR
	fBR_TR_TL.edges[0].oppositeHalfEdge = indexFor(fTL_BL_BR, 0); // eTL_BR.opposite = eBR_TL
	
	//Weld outer triangles to inner square
	fTL_TR_vInf.edges[1].oppositeHalfEdge = indexFor(fBR_TR_TL, 2); //eTL_TR.opposite = eTR_TL
	fBR_TR_TL.edges[2].oppositeHalfEdge = indexFor(fTL_TR_vInf, 1); //eTR_TL.opposite = eTL_TR

	fBL_TL_vInf.edges[1].oppositeHalfEdge = indexFor(fTL_BL_BR, 1); //eBL_TL.opposite = eTL_BL
	fTL_BL_BR.edges[1].oppositeHalfEdge = indexFor(fBL_TL_vInf, 1); //eTL_BL.opposite = eBL_TL

	fBR_BL_vInf.edges[1].oppositeHalfEdge = indexFor(fTL_BL_BR, 2); //eBR_BL.opposite = eBL_BR
	fTL_BL_BR.edges[2].oppositeHalfEdge = indexFor(fBR_BL_vInf, 1); //eBL_BR.opposite = eBR_BL
	
	fTR_BR_vInf.edges[1].oppositeHalfEdge = indexFor(fBR_TR_TL, 1); //eTR_BR.opposite = eBR_TR
	fBR_TR_TL.edges[1].oppositeHalfEdge = indexFor(fTR_BR_vInf, 1); //eBR_TR.opposite = eTR_BR

	//Weld neighbouring edges of pyramid together
	fTL_TR_vInf.edges[0].oppositeHalfEdge = indexFor(fBL_TL_vInf, 2); //eInf_TL.opposite = eTL_Inf
	fBL_TL_vInf.edges[2].oppositeHalfEdge = indexFor(fTL_TR_vInf, 0); //eTL_Inf.opposite = eInf_TL
	
	fBL_TL_vInf.edges[0].oppositeHalfEdge = indexFor(fBR_BL_vInf, 2); //eInf_BL.opposite = eBL_Inf
	fBR_BL_vInf.edges[2].oppositeHalfEdge = indexFor(fBL_TL_vInf, 0); //eBL_Inf.opposite = eInf_BL
	
	fBR_BL_vInf.edges[0].oppositeHalfEdge = indexFor(fTR_BR_vInf, 2); // eInf_BR.opposite = eBR_Inf
	fTR_BR_vInf.edges[2].oppositeHalfEdge = indexFor(fBR_BL_vInf, 0); // eBR_Inf.opposite = eInf_BR

	fTR_BR_vInf.edges[0].oppositeHalfEdge = indexFor(fTL_TR_vInf, 2); // eInf_TR.opposite = eTR_Inf
	fTL_TR_vInf.edges[2].oppositeHalfEdge = indexFor(fTR_BR_vInf, 0); // eTR_Inf.opposite = eInf_TR

	//Set edge on each vertex
	vInf.edge = indexFor(fTL_TR_vInf, 2);
	vBL.edge = indexFor(fTL_BL_BR, 1);
	vBR.edge = indexFor(fTL_BL_BR, 2);
	vTR.edge = indexFor(fBR_TR_TL, 1);
	vTL.edge = indexFor(fBR_TR_TL, 2);
}

Delaunay::Index Delaunay::Mesh::InsertVertex(double x, double y)
{
	vertices.Reserve(1);
	Index vertex = Vertex::InvalidIndex;
	//Make sure the vertex is inside the bounding box the mesh was initialised with.
	//Locate the primitive the vertex falls on
	LocateResult result = this->Locate(x, y);
	switch (result.type) {
	case LocateResult::Vertex:
		vertex = result.object;
		break;
	case LocateResult::Edge:
		vertex = splitEdge(result.object, x, y);
		break;
	case LocateResult::Face:
		vertex = splitFace(result.object, x, y);
		break;
	}
	//Restore delaunay condition
	//TODO: In rare cases this code can fail and cause a crash because edge references are not stable and may occasionally hit an invalid half-edge index
	while (!edgesToCheck.Empty()) {
		Index h = edgesToCheck.PopFront();
		
		if (!isHalfEdgeConstrained(h) && !isDelaunay(h)) {
			h = flipEdge(h);
			HalfEdge & current = edgeAt(h);
			HalfEdge & opposite = edgeAt(current.oppositeHalfEdge);
			//If the newly flipped edge is now pointing towards the center, it was not before
			if (opposite.destinationVertex == centerVertex || current.destinationVertex == centerVertex) { 
				edgesToCheck.Add(Face::prevHalfEdge(indexFor(opposite)));
				edgesToCheck.Add(Face::nextHalfEdge(indexFor(current)));
			}
			else {
				edgesToCheck.Add(Face::prevHalfEdge(h));
				edgesToCheck.Add(Face::nextHalfEdge(h));
			}
			
		}
	}
	

	return vertex;
}

bool Delaunay::Mesh::DeleteVertex(Index index)
{
	// Few typical cases
	//If vertex.constraintCount == 0 it is completely free of constraint --> Remove immediately
	//If vertex.constraintCount == 2 it is the end point of a constraint --> Cannot be removed until the constraint is broken first
	//If vertex.constraintCount >= 4 it is either a midpoint of a constraint and/or belongs to multiple constraints --> Ignore these and return false.
	Vertex & vertex = this->vertices[index];
	Oryol::Array<Index> bound; //Outer bound is required by untriangulate
	Oryol::Array<Index> taggedFaces; //Contains the face indices that will be removed from the mesh
	if (vertex.constraintCount == 0) {
		for (Index h : vertex.OutgoingEdges(*this)) {
			bound.Add(edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge);
			taggedFaces.Add(h / 4);
		}
	}
	else if (vertex.constraintCount == 4) {

	}

	return false;
}

Delaunay::Index Delaunay::Mesh::InsertConstraintSegment(double x1, double y1, double x2, double y2)
{
	return Index();
}

void Delaunay::Mesh::DeleteConstraintSegment(Index index)
{

}

Delaunay::Index Delaunay::Mesh::splitEdge(Index h, double x, double y)
{
	centerVertex = Vertex::InvalidIndex;
	edgesToCheck.Clear();
	//Have to specify this otherwise code that uses the indexFor() shortcut for generating indices will fail.
	//The ideal approach would be to use indices the entire way through and would avoid this particular "hack"
	faces.Reserve(4);
	HalfEdge & eUp_Down = edgeAt(h);
	HalfEdge & eDown_Up = edgeAt(eUp_Down.oppositeHalfEdge);

	Vertex & vDown = this->vertices[eUp_Down.destinationVertex];
	Vertex & vUp = this->vertices[eDown_Up.destinationVertex];
	//Check to see if the point is close enough to a vertex and return the corresponding halfedge
	if (Geo2D::DistanceSquared(vUp.position - glm::dvec2{x, y}) <= EPSILON_SQUARED)
		return h;
	if (Geo2D::DistanceSquared(vDown.position - glm::dvec2{ x,y }) <= EPSILON_SQUARED)
		return eUp_Down.oppositeHalfEdge;
	//These are the edges which surround the pair of faces we will be removing
	HalfEdge & eDown_Right = edgeAt(Face::nextHalfEdge(indexFor(eUp_Down)));
	//Oryol::Log::Info("eDown_Right: %d Opposite: %d\n", indexFor(eDown_Right),eDown_Right.oppositeHalfEdge);
	HalfEdge & eUp_Left = edgeAt(Face::nextHalfEdge(indexFor(eDown_Up)));
	//Oryol::Log::Info("eUp_Left: %d Opposite: %d\n", indexFor(eUp_Left), eUp_Left.oppositeHalfEdge);
	HalfEdge & eRight_Up = edgeAt(Face::prevHalfEdge(indexFor(eUp_Down)));
	//Oryol::Log::Info("eRight_Up: %d Opposite: %d\n", indexFor(eRight_Up), eRight_Up.oppositeHalfEdge);
	HalfEdge & eLeft_Down = edgeAt(Face::prevHalfEdge(indexFor(eDown_Up)));

	//Make sure we store references to all these vertices because it will make things much easier.
	Vertex & vLeft = this->vertices[eUp_Left.destinationVertex];
	Vertex & vRight = this->vertices[eDown_Right.destinationVertex];

	//Create the 4 new faces that we need
	Face & fUp_Left_Center = requestFace();
	Face & fLeft_Down_Center = requestFace();
	Face & fDown_Right_Center = requestFace();
	Face & fRight_Up_Center = requestFace();
	//Create the new vertex we need
	//As the point may not be exactly on the line it must be projected on the line segment
	glm::dvec2 p = Geo2D::OrthogonallyProjectPointOnLineSegment(vDown.position, vUp.position, { x,y });
	Vertex & vCenter = requestVertex(p.x, p.y);
	vCenter.edge = indexFor(fUp_Left_Center, 2);
	//Set vertex data on the face edges
	fUp_Left_Center.edges[0].destinationVertex = indexFor(vUp);
	fUp_Left_Center.edges[1].destinationVertex = indexFor(vLeft);
	fUp_Left_Center.edges[2].destinationVertex = indexFor(vCenter);

	fLeft_Down_Center.edges[0].destinationVertex = indexFor(vLeft);
	fLeft_Down_Center.edges[1].destinationVertex = indexFor(vDown);
	fLeft_Down_Center.edges[2].destinationVertex = indexFor(vCenter);

	fDown_Right_Center.edges[0].destinationVertex = indexFor(vDown);
	fDown_Right_Center.edges[1].destinationVertex = indexFor(vRight);
	fDown_Right_Center.edges[2].destinationVertex = indexFor(vCenter);

	fRight_Up_Center.edges[0].destinationVertex = indexFor(vRight);
	fRight_Up_Center.edges[1].destinationVertex = indexFor(vUp);
	fRight_Up_Center.edges[2].destinationVertex = indexFor(vCenter);

	//Weld the new faces together
	fUp_Left_Center.edges[2].oppositeHalfEdge = indexFor(fLeft_Down_Center, 0); //eLeft_Center.opposite = eCenter_Left
	fLeft_Down_Center.edges[0].oppositeHalfEdge = indexFor(fUp_Left_Center, 2); //eCenter_Left.opposite = eLeft_Center

	fLeft_Down_Center.edges[2].oppositeHalfEdge = indexFor(fDown_Right_Center, 0); //eDown_Center.opposite = eCenter_Down
	fDown_Right_Center.edges[0].oppositeHalfEdge = indexFor(fLeft_Down_Center, 2); //eCenter_Down.opposite = eDown_Center

	fDown_Right_Center.edges[2].oppositeHalfEdge = indexFor(fRight_Up_Center, 0); //eRight_Center.opposite = eCenter_Right
	fRight_Up_Center.edges[0].oppositeHalfEdge = indexFor(fDown_Right_Center, 2); //eCenter_Right.opposite = eRight_Center

	fRight_Up_Center.edges[2].oppositeHalfEdge = indexFor(fUp_Left_Center, 0); //eUp_Center.opposite = eCenter_Up
	fUp_Left_Center.edges[0].oppositeHalfEdge = indexFor(fRight_Up_Center, 2); //eCenter_Up.opposite = eUp_Center

	//Weld the new faces to the existing faces
	fLeft_Down_Center.edges[1].oppositeHalfEdge = eLeft_Down.oppositeHalfEdge;
	edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = indexFor(fLeft_Down_Center, 1);

	fRight_Up_Center.edges[1].oppositeHalfEdge = eRight_Up.oppositeHalfEdge;
	edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = indexFor(fRight_Up_Center, 1);

	fUp_Left_Center.edges[1].oppositeHalfEdge = eUp_Left.oppositeHalfEdge;
	edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = indexFor(fUp_Left_Center, 1);

	fDown_Right_Center.edges[1].oppositeHalfEdge = eDown_Right.oppositeHalfEdge;
	edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = indexFor(fDown_Right_Center, 1);

	//Patch edge references in each vertex to refer to the newly created faces (if needed)
	if (indexFor(eDown_Right) == vRight.edge)
		vRight.edge = indexFor(fDown_Right_Center, 1);
	if (indexFor(eLeft_Down) == vDown.edge || indexFor(eUp_Down) == vDown.edge)
		vDown.edge = indexFor(fLeft_Down_Center, 1);
	if (indexFor(eUp_Left) == vLeft.edge)
		vLeft.edge = indexFor(fUp_Left_Center, 1);
	if (indexFor(eRight_Up) == vUp.edge || indexFor(eDown_Up) == vUp.edge)
		vUp.edge = indexFor(fRight_Up_Center, 1);

	//Copy constraint state to the newly constructed faces
	if (isHalfEdgeConstrained(indexFor(eDown_Right)))
		fDown_Right_Center.flags |= (1 << 2);
	if (isHalfEdgeConstrained(indexFor(eLeft_Down)))
		fLeft_Down_Center.flags |= (1 << 2);
	if (isHalfEdgeConstrained(indexFor(eUp_Left)))
		fUp_Left_Center.flags |= (1 << 2);
	if (isHalfEdgeConstrained(indexFor(eRight_Up)))
		fRight_Up_Center.flags |= (1 << 2);

	if (!isHalfEdgeReal(eDown_Up.oppositeHalfEdge)) {
		//eUp_Down is not real therefore tag the appropriate new faces as not real
		fRight_Up_Center.flags |= 0x1;
		fDown_Right_Center.flags |= 0x1;
	}
	else if (!isHalfEdgeReal(eUp_Down.oppositeHalfEdge)) {
		//eDown_Up is not real therefore tag the new faces as not real 
		fUp_Left_Center.flags |= 0x1;
		fLeft_Down_Center.flags |= 0x1;
	}
	//Handle the fact that the split edge is constrained
	if (isHalfEdgeConstrained(h)) {
		//Tag the half-edges as constrained
		fRight_Up_Center.flags |= (1 << 1);
		fDown_Right_Center.flags |= (1 << 1);
		fUp_Left_Center.flags |= (1 << 1);
		fDown_Right_Center.flags |= (1 << 1);

		//Find the constraint segments associated with the original half-edge and insert the new vertex between the two existing vertices
		//A typical approach would be to have a list on each half-edge but since that wouldn't work with the Face == 4xHalf-Edge constraint
		//We check the Map for vertices to constraint segments
	}
	//Recycle the old faces
	recycleFace(eUp_Down.oppositeHalfEdge / 4);
	recycleFace(h / 4);
	//Tag the edges for delaunay condition checking
	centerVertex = indexFor(vCenter);

	edgesToCheck.Add(indexFor(fRight_Up_Center,1));
	edgesToCheck.Add(indexFor(fUp_Left_Center,1));
	edgesToCheck.Add(indexFor(fDown_Right_Center,1));
	edgesToCheck.Add(indexFor(fRight_Up_Center,1));

	return centerVertex;
}

//TODO: Replace any vertex references with their corresponding index to avoid unnecessary calculations
Delaunay::Index Delaunay::Mesh::flipEdge(Index h)
{
	if (h % 4 == 0)
		o_error("Not a valid half-edge index");
	faces.Reserve(2);
	//The triangles owning these half-edges will be replaced with a triangle pair with a common edge running left to right
	HalfEdge & eUp_Down = edgeAt(h);
	HalfEdge & eDown_Up = edgeAt(eUp_Down.oppositeHalfEdge);

	//These belong to the two faces which will be removed so their opposites must be patched.
	HalfEdge & eDown_Right = edgeAt(Face::nextHalfEdge(indexFor(eUp_Down)));
	HalfEdge & eUp_Left = edgeAt(Face::nextHalfEdge(indexFor(eDown_Up)));
	HalfEdge & eRight_Up = edgeAt(Face::prevHalfEdge(indexFor(eUp_Down)));
	HalfEdge & eLeft_Down = edgeAt(Face::prevHalfEdge(indexFor(eDown_Up)));

	//Make sure we store references to all these vertices because it will make things much easier.
	Vertex & vDown = this->vertices[eUp_Down.destinationVertex];
	Vertex & vUp = this->vertices[eDown_Up.destinationVertex];
	Vertex & vLeft = this->vertices[eUp_Left.destinationVertex];
	Vertex & vRight = this->vertices[eDown_Right.destinationVertex];

	//Request new faces; Ideally you should be able to request a pair/triplet of faces at the same time
	Face & fLeft_Right_Up = requestFace();
	Face & fRight_Left_Down = requestFace();

	//Set up vertex references for both new faces
	fLeft_Right_Up.edges[0].destinationVertex = indexFor(vLeft);
	fLeft_Right_Up.edges[1].destinationVertex = indexFor(vRight);
	fLeft_Right_Up.edges[2].destinationVertex = indexFor(vUp);

	fRight_Left_Down.edges[0].destinationVertex = indexFor(vRight);
	fRight_Left_Down.edges[1].destinationVertex = indexFor(vLeft);
	fRight_Left_Down.edges[2].destinationVertex = indexFor(vDown);

	//Weld the two new faces together along eLeft_Right and eRight_Left
	fLeft_Right_Up.edges[1].oppositeHalfEdge = indexFor(fRight_Left_Down, 1);
	fRight_Left_Down.edges[1].oppositeHalfEdge = indexFor(fLeft_Right_Up, 1);

	//Patch the new faces with the correct opposite references
	fLeft_Right_Up.edges[0].oppositeHalfEdge = eUp_Left.oppositeHalfEdge;
	edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = indexFor(fLeft_Right_Up, 0);

	fLeft_Right_Up.edges[2].oppositeHalfEdge = eRight_Up.oppositeHalfEdge;
	edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = indexFor(fLeft_Right_Up, 2);

	fRight_Left_Down.edges[0].oppositeHalfEdge = eDown_Right.oppositeHalfEdge;
	edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = indexFor(fRight_Left_Down, 0);

	fRight_Left_Down.edges[2].oppositeHalfEdge = eLeft_Down.oppositeHalfEdge;
	edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = indexFor(fRight_Left_Down, 2);

	//Patch edge references in each vertex to refer to the newly created faces (if needed)
	if (indexFor(eDown_Right) == vRight.edge)
		vRight.edge = indexFor(fRight_Left_Down, 0);
	if (indexFor(eLeft_Down) == vDown.edge || indexFor(eUp_Down) == vDown.edge)
		vDown.edge = indexFor(fRight_Left_Down, 2);
	if (indexFor(eUp_Left) == vLeft.edge)
		vLeft.edge = indexFor(fLeft_Right_Up, 0);
	if (indexFor(eRight_Up) == vUp.edge || indexFor(eDown_Up) == vUp.edge);
		vUp.edge = indexFor(fLeft_Right_Up, 2);
	//Copy constraint state to the newly constructed faces
	if (isHalfEdgeConstrained(indexFor(eDown_Right)))
		fRight_Left_Down.flags |= (1 << 1);
	if (isHalfEdgeConstrained(indexFor(eLeft_Down)))
		fRight_Left_Down.flags |= (1 << 3);
	if (isHalfEdgeConstrained(indexFor(eUp_Left)))
		fLeft_Right_Up.flags |= (1 << 1);
	if (isHalfEdgeConstrained(indexFor(eRight_Up)))
		fLeft_Right_Up.flags |= (1 << 3);

	//Release the old faces
	recycleFace(indexFor(eDown_Up) / 4); //Left face
	recycleFace(indexFor(eUp_Down) / 4); //Right face;
	
	return indexFor(fLeft_Right_Up, 1); //Returns the Half-Edge running left to right for the newly created triangle pair
}

//TODO: Replace any vertex references with their corresponding index to avoid unnecessary calculations
Delaunay::Index Delaunay::Mesh::splitFace(Index f, double x, double y)
{
	edgesToCheck.Clear();
	faces.Reserve(3);
	Face & fA_B_C = this->faces[f];
	Vertex & vA = this->vertices[fA_B_C.edges[0].destinationVertex];
	Vertex & vB = this->vertices[fA_B_C.edges[1].destinationVertex];
	Vertex & vC = this->vertices[fA_B_C.edges[2].destinationVertex];

	Face & fA_B_New = requestFace();
	Face & fB_C_New = requestFace();
	Face & fC_A_New = requestFace();

	Vertex & vNew = requestVertex(x, y);
	vNew.edge = indexFor(fA_B_New, 2);

	fA_B_New.edges[0].destinationVertex = indexFor(vA);
	fA_B_New.edges[1].destinationVertex = indexFor(vB);
	fA_B_New.edges[2].destinationVertex = indexFor(vNew);

	fB_C_New.edges[0].destinationVertex = indexFor(vB);
	fB_C_New.edges[1].destinationVertex = indexFor(vC);
	fB_C_New.edges[2].destinationVertex = indexFor(vNew);

	fC_A_New.edges[0].destinationVertex = indexFor(vC);
	fC_A_New.edges[1].destinationVertex = indexFor(vA);
	fC_A_New.edges[2].destinationVertex = indexFor(vNew);

	//Copy eC_A from old face
	fC_A_New.edges[1].oppositeHalfEdge = fA_B_C.edges[0].oppositeHalfEdge;
	edgeAt(fA_B_C.edges[0].oppositeHalfEdge).oppositeHalfEdge = indexFor(fC_A_New, 1);
	
	//Copy eA_B from old face
	fA_B_New.edges[1].oppositeHalfEdge = fA_B_C.edges[1].oppositeHalfEdge;
	edgeAt(fA_B_C.edges[1].oppositeHalfEdge).oppositeHalfEdge = indexFor(fA_B_New, 1);
	
	//Copy eB_C from old face
	fB_C_New.edges[1].oppositeHalfEdge = fA_B_C.edges[2].oppositeHalfEdge;
	edgeAt(fA_B_C.edges[2].oppositeHalfEdge).oppositeHalfEdge = indexFor(fB_C_New, 1);


	fA_B_New.edges[0].oppositeHalfEdge = indexFor(fC_A_New, 2); //eNew_A.opposite = eA_New
	fC_A_New.edges[2].oppositeHalfEdge = indexFor(fA_B_New, 0); //eA_New.opposite = eNew_A

	fB_C_New.edges[0].oppositeHalfEdge = indexFor(fA_B_New, 2); //eNew_B.opposite = eB_New
	fA_B_New.edges[2].oppositeHalfEdge = indexFor(fB_C_New, 0); //eB_New.opposite = eNew_B

	fC_A_New.edges[0].oppositeHalfEdge = indexFor(fB_C_New, 2); //eNew_C.opposite = eC_New
	fB_C_New.edges[2].oppositeHalfEdge = indexFor(fC_A_New, 0); //eC_New.opposite = eNew_C

	//Patch edge references in each vertex to refer to the newly created faces (if needed)
	if (indexFor(fA_B_C, 0) == vA.edge)
		vA.edge = indexFor(fC_A_New, 1);
	if (indexFor(fA_B_C, 1) == vB.edge)
		vB.edge = indexFor(fA_B_New, 1);
	if (indexFor(fA_B_C, 2) == vC.edge)
		vC.edge = indexFor(fB_C_New, 1);

	//Copy across constraint state to the new faces
	if (isHalfEdgeConstrained(indexFor(fA_B_C, 0)))
		fC_A_New.flags |= (1 << 2);
	if (isHalfEdgeConstrained(indexFor(fA_B_C, 1)))
		fA_B_New.flags |= (1 << 2);
	if (isHalfEdgeConstrained(indexFor(fA_B_C, 2)))
		fB_C_New.flags |= (1 << 2);

	//Free the old face
	recycleFace(indexFor(fA_B_C));

	//Tag the modified edges for delaunay condition checking
	edgesToCheck.Add(indexFor(fA_B_New, 1));
	edgesToCheck.Add(indexFor(fB_C_New, 1));
	edgesToCheck.Add(indexFor(fC_A_New, 1));
	return indexFor(vNew);
}

Delaunay::Mesh::LocateResult Delaunay::Mesh::Locate(double x, double y)
{
	LocateResult result{ Face::InvalidIndex, LocateResult::None };
	Index bestVertex = Vertex::InvalidIndex;
	{
		//Seed the random generator
		std::srand(x * 10 + y * 4);
		//Find closest vertex by randomly sampling the vertices (excluding the infinite vertex)
		int vertexSampleCount = std::pow(this->cachedVertices.Size(), 1 / 3.);
		double minDistanceSquared = std::numeric_limits<double>::infinity();
		for (int i = 0; i < vertexSampleCount; i++) {
			Index cacheIndex = rand_range(1, this->cachedVertices.Size() - 1);
			Index index = this->cachedVertices.ValueAtIndex(cacheIndex);
			Vertex & vertex = this->vertices[index];
			double distanceSquared = Geo2D::DistanceSquared(glm::dvec2(x,y)-vertex.position);
			if (distanceSquared < minDistanceSquared) {
				minDistanceSquared = distanceSquared;
				bestVertex = index;
			}
		}
	}
	//Start jump-and-walk search with the first face associated with the best vertex
	Index currentEdge = this->vertices[bestVertex].edge;
	Index currentFace = currentEdge / 4;
	{
		//TODO: Extract this functionality into its own iterator
		const Index firstFace = currentFace;
		//Ensure the selected face is real and that we havent looped back to where we were.
		while (this->faces[currentFace].flags & 0x01 > 0) {
			Index nextEdge = edgeAt(currentEdge).oppositeHalfEdge;
			currentFace = nextEdge / 4;
			currentEdge = Face::nextHalfEdge(nextEdge);
			if (currentFace == firstFace) return result;
		}
	}
	
	Oryol::Set<Index> visitedFaces;
	int iterations = 0;
	while (visitedFaces.Contains(currentFace) || !(result = isInFace(x, y, this->faces[currentFace]))) {
		visitedFaces.Add(currentFace);
		iterations++;
		if (iterations == 50) {
			//Log this as it is taking longer than expected
		}
		else if (iterations > 1000) {
			//Bail out if too many iterations have elapsed
			result.type = LocateResult::None;
			break;
		}
		Index nextFace = Face::InvalidIndex;
		//Find the best direction to look in.
		for (HalfEdge & h : this->faces[currentFace].edges) {
			//Determine if the position falls to the right of the current half edge (thus outside of the current face)
			Vertex & originVertex = this->vertices[edgeAt(h.oppositeHalfEdge).destinationVertex];
			Vertex & destinationVertex = this->vertices[h.destinationVertex];
			if (Geo2D::DetermineSide(originVertex.position, destinationVertex.position, { x,y }) == -1) {
				nextFace = h.oppositeHalfEdge / 4;
				break;
			}
		}
		if (nextFace != Face::InvalidIndex) {
			currentFace = nextFace;
		}
		else {
			result.type = LocateResult::None;
			break; //Something has gone wrong so log it and bail
		}
	}
	
	return result;
}

void Delaunay::Mesh::SetDebugDraw(DebugDraw * debug)
{
	this->debugDraw = debug;
}

void Delaunay::Mesh::DrawDebugData()
{
	if (debugDraw == nullptr) return;
	Oryol::Set<Index> visitedFaces;
	for (Index vIndex : cachedVertices) {
		Vertex & vertex = this->vertices[vIndex];
		debugDraw->DrawVertex(vertex.position);
		const Index first = vertex.edge;
		Index h = first;
		do {
			Face & face = this->faces[h / 4];
			HalfEdge & current = edgeAt(h);
			//Render Edges
			if ((h < current.oppositeHalfEdge) && (current.destinationVertex != 0) && (edgeAt(current.oppositeHalfEdge).destinationVertex != 0)) {
				Vertex & origin = this->vertices[current.destinationVertex];
				Vertex & destination = this->vertices[edgeAt(current.oppositeHalfEdge).destinationVertex];
				debugDraw->DrawEdge(origin.position, destination.position, isHalfEdgeConstrained(h));
			}
			//TODO: Render Faces
			if(!visitedFaces.Contains(h/4)){
				visitedFaces.Add(h / 4);
			}
			h = edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge;
		} while (h != first);
		
	}
}


Delaunay::Mesh::LocateResult Delaunay::Mesh::isInFace(double x, double y, Face & face)
{
	LocateResult result{ Vertex::InvalidIndex, LocateResult::None };
	glm::dvec2 p = { x,y };

	glm::dvec2 v1 = this->vertices[face.edges[0].destinationVertex].position;
	glm::dvec2 v2 = this->vertices[face.edges[1].destinationVertex].position;
	glm::dvec2 v3 = this->vertices[face.edges[2].destinationVertex].position;
	
	if (Geo2D::Sign(v3, v1, p) >= 0.0 && Geo2D::Sign(v1, v2, p) >= 0.0 && Geo2D::Sign(v2, v3, p) >= 0.0) {
		//Cache the proximity info rather than calculate it more than necessary
		bool proximity[3] = {
			Geo2D::DistanceSquaredPointToLineSegment(v3,v1,p) <= EPSILON_SQUARED,
			Geo2D::DistanceSquaredPointToLineSegment(v1,v2,p) <= EPSILON_SQUARED,
			Geo2D::DistanceSquaredPointToLineSegment(v2,v3,p) <= EPSILON_SQUARED
		};

		if (proximity[0]) {
			if (proximity[1]) {
				result.object = face.edges[0].destinationVertex;
				result.type = LocateResult::Vertex;
			}
			else if (proximity[2]) {
				result.object = face.edges[2].destinationVertex;
				result.type = LocateResult::Vertex;
			}
			else {
				result.object = indexFor(face, 0); //eV3_V1
				result.type = LocateResult::Edge;
			}
		}
		else if (proximity[1]) {
			if (proximity[2]) {
				result.object = face.edges[1].destinationVertex;
				result.type = LocateResult::Vertex;
			}
			else {
				result.object = indexFor(face, 1); //eV1_V2
				result.type = LocateResult::Edge;
			}
		}
		else if (proximity[2]) {
			result.object = indexFor(face, 2); //eV2_V3
			result.type = LocateResult::Edge;
		}
		else {
			result.object = indexFor(face);
			result.type = LocateResult::Face;
		}
	}
	return result;
}

bool Delaunay::Mesh::isDelaunay(Index h)
{
	Face & face = this->faces[h / 4];
	Vertex & vA = this->vertices[face.edges[0].destinationVertex];
	Vertex & vB = this->vertices[face.edges[1].destinationVertex];
	Vertex & vC = this->vertices[face.edges[2].destinationVertex];

	Vertex & vOpp = this->vertices[edgeAt(Face::nextHalfEdge(edgeAt(h).oppositeHalfEdge)).destinationVertex];


	glm::dvec2 circumcenter = Geo2D::ComputeCircumcenter(vA.position, vB.position, vC.position);
	double squaredRadius = Geo2D::DistanceSquared(vA.position - circumcenter);
	double squaredDistance = Geo2D::DistanceSquared(vOpp.position - circumcenter);
	return squaredDistance >= squaredRadius;
}

inline Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) {
	return *(reinterpret_cast<HalfEdge*>(&faces[0])+ index);
}
inline const Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) const {
	return *(reinterpret_cast<const HalfEdge*>(&faces[0]) + index);
}

inline Delaunay::Index Delaunay::Mesh::Face::nextHalfEdge(Index h) {
	++h;
	return (h & 3) != 0 ? h : h - 3;
}

inline Delaunay::Index Delaunay::Mesh::Face::prevHalfEdge(Index h) {
	--h;
	return (h & 3) != 0 ? h : h + 3;
}

inline Delaunay::Index Delaunay::Mesh::indexFor(Face & face) {
	return faces.Distance(face);
}

inline Delaunay::Index Delaunay::Mesh::indexFor(Vertex & vertex) {
	return vertices.Distance(vertex);
}

inline Delaunay::Index Delaunay::Mesh::indexFor(HalfEdge & edge) {
	return faces.Distance(edge);
}

//Compute index for edge within a face.
inline Delaunay::Index Delaunay::Mesh::indexFor(Face & face, Index edge) {
	return indexFor(face) * 4 + edge + 1;
}

inline Delaunay::Mesh::Face & Delaunay::Mesh::requestFace() {
	return faces.Request();
}

void Delaunay::Mesh::recycleFace(Index f)
{
	this->faces[f] = Face(); //This is done such that any attempt to use recycled face will become apparent
	faces.Release(this->faces[f]);
}

inline Delaunay::Mesh::Vertex & Delaunay::Mesh::requestVertex(double x, double y, bool cache) {
	Vertex * vertex = &vertices.Request(x,y,-1);
	if (cache)
		cachedVertices.Add(indexFor(*vertex));
	return *vertex;
}

Delaunay::Mesh::HalfEdge::HalfEdge()
	: destinationVertex(Vertex::InvalidIndex), oppositeHalfEdge(HalfEdge::InvalidIndex){}

Delaunay::Mesh::Face::Face() 
	: flags(0), matId(-1) {}

Delaunay::Index Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::operator*()
{
	return current;
}

void Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::operator++()
{	
	current = mesh.edgeAt(Mesh::Face::nextHalfEdge(current)).oppositeHalfEdge;
	if (current == first)
		current = Mesh::HalfEdge::InvalidIndex;
}

bool Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::operator!=(const IncomingHalfEdgeIterator & rhs)
{
	return current != Mesh::HalfEdge::InvalidIndex;
}

Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::begin()
{
	return *this;
}
Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::end()
{
	return *this;
}

Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::IncomingHalfEdgeIterator(Mesh & mesh, Index first)
	: mesh(mesh), current(first), first(first)
{
	/*
	current = mesh.vertices[vertex].edge;
	if (mesh.edgeAt(current).destinationVertex != vertex)
	current = mesh.edgeAt(current).oppositeHalfEdge;
	first = current;
	*/
	
}

Delaunay::Index Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::operator*()
{
	return current;
}

Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::OutgoingHalfEdgeIterator(Mesh & mesh, Index first)
	: mesh(mesh), current(first), first(first)
{
	/*
	current = mesh.vertices[vertex].edge;
	if (mesh.edgeAt(current).destinationVertex == vertex)
	current = mesh.edgeAt(current).oppositeHalfEdge;
	first = current;
	*/
}
void Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::operator++() {
	current = Mesh::Face::nextHalfEdge(mesh.edgeAt(current).oppositeHalfEdge);
	if (current == first)
		current = Mesh::HalfEdge::InvalidIndex;
}
bool Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::operator!=(const OutgoingHalfEdgeIterator & rhs)
{
	return current != Mesh::HalfEdge::InvalidIndex;
}
Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::begin()
{
	return *this;
}
Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator & Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::end()
{
	return *this;
}

inline Delaunay::Mesh::Vertex::Vertex(double x, double y, Index e) : position(x, y), edge(e), constraintCount(0), endPointCount(0) {}

Mesh::HalfEdge::IncomingHalfEdgeIterator Mesh::Vertex::IncomingEdges(Mesh & mesh) {
	return HalfEdge::IncomingHalfEdgeIterator(mesh, edge);
}
Mesh::HalfEdge::OutgoingHalfEdgeIterator Mesh::Vertex::OutgoingEdges(Mesh & mesh) {
	return HalfEdge::OutgoingHalfEdgeIterator(mesh, edge);
}