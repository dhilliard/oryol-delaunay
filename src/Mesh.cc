#include "Mesh.h"
#include "Geo2D.h"
#include <cmath>
#include <limits>
#include "Core/Containers/Set.h"
#include "Core/Containers/Queue.h"
#include "Core/Assertion.h"
#include "glm/geometric.hpp"

constexpr double EPSILON = 0.01;
constexpr double EPSILON_SQUARED = 0.0001;
using namespace Geo2D;
using namespace Delaunay;
using Index = Mesh::HalfEdge::Index;

int rand_range(int min, int max) {
	return min + (float)(rand() / RAND_MAX)*(max - min);
}

struct Mesh::Impl {
public:
    static void LogHalfEdge(Mesh & mesh, Index h){
        HalfEdge & edge = mesh.edgeAt(h);
        HalfEdge & next = mesh.edgeAt(Face::nextHalfEdge(h));
        HalfEdge & prev = mesh.edgeAt(Face::prevHalfEdge(h));
        //These are a very useful debugging aid
        Oryol::Log::Info("\nEdge: (id:%u, origin: %u, destination: %u, constrained: %c)\n", h, prev.destinationVertex, edge.destinationVertex, edge.constrained ? 'Y' : 'N');
        Oryol::Log::Info("Next: (id:%u, origin: %u, destination: %u, constrained: %c)\n", Face::nextHalfEdge(h), edge.destinationVertex, next.destinationVertex, next.constrained ? 'Y' : 'N');
        Oryol::Log::Info("Prev: (id:%u, origin: %u, destination: %u, constrained: %c)\n", Face::prevHalfEdge(h), next.destinationVertex, prev.destinationVertex, prev.constrained ? 'Y' : 'N');
    }
    static Index GetOriginVertex(Mesh & mesh, Index h){
        return mesh.edgeAt(Mesh::Face::prevHalfEdge(h)).destinationVertex;
    }
	static bool CheckFaceIsCounterClockwise(Mesh & mesh, Index a, Index b, Index c) {
		Vertex & vA = mesh.vertices[a];
		Vertex & vB = mesh.vertices[b];
		Vertex & vC = mesh.vertices[c];
		return Geo2D::CounterClockwise(vA.position, vB.position, vC.position);
	}
    static Mesh::LocateRef IsInFace(Mesh & mesh, Index faceIndex, const glm::dvec2 & p)
    {
        LocateRef result{ HalfEdge::InvalidIndex, LocateRef::None};
        Face & face = mesh.faces[faceIndex];
        
        glm::dvec2 v1 = mesh.vertices[face.edges[0].destinationVertex].position;
        glm::dvec2 v2 = mesh.vertices[face.edges[1].destinationVertex].position;
        glm::dvec2 v3 = mesh.vertices[face.edges[2].destinationVertex].position;
        
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
                    result.type = LocateRef::Vertex;
                }
                else if (proximity[2]) {
                    result.object = face.edges[2].destinationVertex;
                    result.type = LocateRef::Vertex;
                }
                else {
                    result.object = faceIndex * 4 + 1; //eV3_V1
                    result.type = LocateRef::Edge;
                }
            }
            else if (proximity[1]) {
                if (proximity[2]) {
                    result.object = face.edges[1].destinationVertex;
                    result.type = LocateRef::Vertex;
                }
                else {
                    result.object = faceIndex * 4 + 2; //eV1_V2
                    result.type = LocateRef::Edge;
                }
            }
            else if (proximity[2]) {
                result.object = faceIndex * 4 + 3; //eV2_V3
                result.type = LocateRef::Edge;
            }
            else {
                result.object = faceIndex;
                result.type = LocateRef::Face;
            }
            o_assert(result.type != LocateRef::None);
        }
        return result;
    }

	static bool IsDelaunay(Mesh & mesh, Index h)
	{
        const HalfEdge & edge = mesh.edgeAt(h);
        Index cw = Face::prevHalfEdge(h);
        Index ccw = Face::nextHalfEdge(h);
        
		Index ivA = edge.destinationVertex;
		Index ivB = mesh.edgeAt(cw).destinationVertex;
		Index ivC = mesh.edgeAt(ccw).destinationVertex;
        Index ivD = mesh.edgeAt(Face::nextHalfEdge(edge.oppositeHalfEdge)).destinationVertex;
        
        const glm::dvec2 & pA = mesh.vertices[ivA].position;
        const glm::dvec2 & pB = mesh.vertices[ivB].position;
        const glm::dvec2 & pC = mesh.vertices[ivC].position;
        const glm::dvec2 & pD = mesh.vertices[ivD].position;

		glm::dvec2 circumcenter = Geo2D::ComputeCircumcenter(pA, pB, pC);
		double squaredRadius = Geo2D::DistanceSquared(pB - circumcenter);
		double squaredDistance = Geo2D::DistanceSquared(pD - circumcenter);
		return squaredDistance >= squaredRadius;
	}
    
    //Returns halfedge from new face pair created as a result of the flip
    static Index FlipEdge(Mesh & mesh, Index h) {
        //o_error("fix_me");
        const HalfEdge eUp_Down = mesh.edgeAt(h);

        const Index iLRU = mesh.faces.Add();
        const Index iRLD = mesh.faces.Add();
        const Index ipL_R = mesh.edgeInfo.Add({iLRU * 4 + 2,{}});
        
        //Recycle the old faces to use them as the new faces
        Face & fLeft_Right_Up = mesh.faces[iLRU];
        Face & fRight_Left_Down = mesh.faces[iRLD];
        {
            //Check preconditions first
			o_assert2(!eUp_Down.constrained,"A constrained edge cannot be flipped");
            //o_assert2(!fLeft_Right_Up.isReal() || !fRight_Left_Down.isReal(), "Both faces must be real in order to be flipped");
            //o_assert2(fLeft_Right_Up.matId == fRight_Left_Down.matId, "The matId should be the same across both faces");
        }
        
        const Index iRight_Up = Face::prevHalfEdge(h);
        const Index iDown_Right = Face::nextHalfEdge(h);
        const Index iUp_Left = Face::nextHalfEdge(eUp_Down.oppositeHalfEdge);
        const Index iLeft_Down = Face::prevHalfEdge(eUp_Down.oppositeHalfEdge);
        
        const HalfEdge eRight_Up = mesh.edgeAt(iRight_Up);
        const HalfEdge eDown_Right = mesh.edgeAt(iDown_Right);
        const HalfEdge eUp_Left = mesh.edgeAt(iUp_Left);
        const HalfEdge eLeft_Down = mesh.edgeAt(iLeft_Down);
        
        
        
        const Index iUp = eRight_Up.destinationVertex;
        const Index iDown = eLeft_Down.destinationVertex;
        const Index iLeft = eUp_Left.destinationVertex;
        const Index iRight = eDown_Right.destinationVertex;

		o_assert(Geo2D::CounterClockwise(mesh.vertices[iLeft].position, mesh.vertices[iRight].position, mesh.vertices[iUp].position));
		o_assert(Geo2D::CounterClockwise(mesh.vertices[iRight].position, mesh.vertices[iLeft].position, mesh.vertices[iDown].position));
    
		//Construct the new faces in place of the old ones
        fLeft_Right_Up.edges[0] = {iLeft, eUp_Left.oppositeHalfEdge, eUp_Left.constrained, eUp_Left.edgePair}; //eUp_Left
        fLeft_Right_Up.edges[1] = {iRight, iRLD * 4 + 2, false, ipL_R}; //eLeft_Right
        fLeft_Right_Up.edges[2] = {iUp, eRight_Up.oppositeHalfEdge, eRight_Up.constrained, eRight_Up.edgePair}; //eRight_Up
        
        fRight_Left_Down.edges[0] = {iRight, eDown_Right.oppositeHalfEdge, eDown_Right.constrained, eDown_Right.edgePair}; //eDown_Right
        fRight_Left_Down.edges[1] = {iLeft, iLRU * 4 + 2, false, ipL_R}; //eRight_Left
        fRight_Left_Down.edges[2] = {iDown, eLeft_Down.oppositeHalfEdge, eLeft_Down.constrained, eLeft_Down.edgePair}; //eLeft_Down
        
		//Patch opposite half edge references
        mesh.edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = iLRU * 4 + 1;
        mesh.edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = iLRU * 4 + 3;
        
        mesh.edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = iRLD * 4 + 1;
        mesh.edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = iRLD * 4 + 3;
        
        mesh.edgeInfo[eUp_Left.edgePair].edge = iLRU * 4 + 1;
        mesh.edgeInfo[eRight_Up.edgePair].edge = iLRU * 4 + 3;
        mesh.edgeInfo[eDown_Right.edgePair].edge = iRLD * 4 + 1;
        mesh.edgeInfo[eLeft_Down.edgePair].edge = iRLD * 4 + 3;
        

		//Patch vertices to refer to the new half edges
		mesh.vertices[iUp].edge = iLRU * 4 + 3;
		mesh.vertices[iDown].edge = iRLD * 4 + 3;
		mesh.vertices[iLeft].edge = iRLD * 4 + 2;
		mesh.vertices[iRight].edge = iLRU * 4 + 2;
        
        mesh.edgeInfo.Erase(eUp_Down.edgePair);
        mesh.faces.Erase(h/4);
        mesh.faces.Erase(eUp_Down.oppositeHalfEdge / 4);

        return iLRU * 4 + 2; //eLeft_Right
		
    }
    
    
    //Recycles one face + creates 2 new faces
    //Returns new vertex that is created as a result of splitting the face
    static Index SplitFace(Mesh & mesh, Index oldFace, const glm::dvec2 & p, Oryol::Array<Index> * edgesToCheck = nullptr){

		const Face fA_B_C = mesh.faces[oldFace];

		const HalfEdge & eC_A = fA_B_C.edges[0];
		const HalfEdge & eA_B = fA_B_C.edges[1];
		const HalfEdge & eB_C = fA_B_C.edges[2];

		const Index vA = eC_A.destinationVertex;
		const Index vB = eA_B.destinationVertex;
		const Index vC = eB_C.destinationVertex;
        
        mesh.faces.Reserve(3);
        Index iC_A_Center = mesh.faces.Add();
        Index iA_B_Center = mesh.faces.Add();
        Index iB_C_Center = mesh.faces.Add();
        
		//Create the new vertex
		const Index vCenter = mesh.vertices.Add({p,iA_B_Center * 4 + 3, 0, 0 });

		o_assert(Geo2D::CounterClockwise(mesh.vertices[vC].position, mesh.vertices[vA].position, mesh.vertices[vCenter].position));
		o_assert(Geo2D::CounterClockwise(mesh.vertices[vA].position, mesh.vertices[vB].position, mesh.vertices[vCenter].position));
		o_assert(Geo2D::CounterClockwise(mesh.vertices[vB].position, mesh.vertices[vC].position, mesh.vertices[vCenter].position));
        
		//Create the new edge pair info structs
        Index ipCenter_C = mesh.edgeInfo.Add({iC_A_Center*4 + 1, {}});
        Index ipCenter_A = mesh.edgeInfo.Add({iA_B_Center * 4 + 1, {}});
        Index ipCenter_B = mesh.edgeInfo.Add({iB_C_Center * 4 + 1, {}});

		
        Face & fC_A_Center = mesh.faces[iC_A_Center];
		fC_A_Center.edges[0] = {vC, iB_C_Center * 4 + 3, false, ipCenter_C}; //eCenter_C
		fC_A_Center.edges[1] = {vA, eC_A.oppositeHalfEdge, eC_A.constrained, eC_A.edgePair}; //eC_A
		fC_A_Center.edges[2] = {vCenter, iA_B_Center * 4 + 1, false, ipCenter_A}; //eA_Center

        Face & fA_B_Center = mesh.faces[iA_B_Center];
		fA_B_Center.edges[0] = {vA, iC_A_Center * 4 + 3, false, ipCenter_A}; //eCenter_A
		fA_B_Center.edges[1] = {vB, eA_B.oppositeHalfEdge, eA_B.constrained, eA_B.edgePair}; //eA_B
		fA_B_Center.edges[2] = {vCenter, iB_C_Center * 4 + 1, false, ipCenter_B}; //eB_Center

        Face & fB_C_Center = mesh.faces[iB_C_Center];
		fB_C_Center.edges[0] = {vB, iA_B_Center * 4 + 3, false, ipCenter_B}; //eCenter_B
		fB_C_Center.edges[1] = {vC, eB_C.oppositeHalfEdge, eB_C.constrained, eB_C.edgePair}; //eB_C
		fB_C_Center.edges[2] = {vCenter, iC_A_Center * 4 + 1, false, ipCenter_C}; //eC_Center
		
		//Patch opposite half edge references
		mesh.edgeAt(eC_A.oppositeHalfEdge).oppositeHalfEdge = iC_A_Center * 4 + 2;
		mesh.edgeAt(eA_B.oppositeHalfEdge).oppositeHalfEdge = iA_B_Center * 4 + 2;
		mesh.edgeAt(eB_C.oppositeHalfEdge).oppositeHalfEdge = iB_C_Center * 4 + 2;
        
        mesh.edgeInfo[eC_A.edgePair].edge = iC_A_Center * 4 + 2;
        mesh.edgeInfo[eA_B.edgePair].edge = iA_B_Center * 4 + 2;
        mesh.edgeInfo[eB_C.edgePair].edge = iB_C_Center * 4 + 2;

		//Patch existing vertices to refer to the new half edges
		mesh.vertices[vA].edge = iC_A_Center * 4 + 2;
		mesh.vertices[vB].edge = iA_B_Center * 4 + 2;
		mesh.vertices[vC].edge = iB_C_Center * 4 + 2;
		
        if(edgesToCheck){
            edgesToCheck->Add(iC_A_Center * 4 + 2); //eC_A
            edgesToCheck->Add(iA_B_Center * 4 + 2); //eA_B
            edgesToCheck->Add(iB_C_Center * 4 + 2); //eB_C
        }
        mesh.faces.Erase(oldFace);
		return vCenter;
    }
    
    //Deletes 2 faces + creates 4 faces
    //Returns new vertex that is created as a result of splitting the edge
    static Index SplitEdge(Mesh & mesh, Index h, const glm::dvec2 & p, Index * centerVertex = nullptr, Oryol::Array<Index> * edgesToCheck = nullptr) {
        //o_error("fix_me");
		
		const HalfEdge eDown_Up = mesh.edgeAt(h);
		const HalfEdge eUp_Down = mesh.edgeAt(eDown_Up.oppositeHalfEdge);

		//Check to see if the point is close enough to a vertex and return the corresponding vertex
		if (Geo2D::DistanceSquared(mesh.vertices[eDown_Up.destinationVertex].position - p) <= EPSILON_SQUARED)
			return eDown_Up.destinationVertex;
		if (Geo2D::DistanceSquared(mesh.vertices[eUp_Down.destinationVertex].position - p) <= EPSILON_SQUARED)
			return eUp_Down.destinationVertex;

		const Index iUp_Left = Face::nextHalfEdge(h);
		const Index iLeft_Down = Face::prevHalfEdge(h);
		const Index iRight_Up = Face::prevHalfEdge(eDown_Up.oppositeHalfEdge);
		const Index iDown_Right = Face::nextHalfEdge(eDown_Up.oppositeHalfEdge);

		const HalfEdge eRight_Up = mesh.edgeAt(iRight_Up);
		const HalfEdge eDown_Right = mesh.edgeAt(iDown_Right);
		const HalfEdge eUp_Left = mesh.edgeAt(iUp_Left);
		const HalfEdge eLeft_Down = mesh.edgeAt(iLeft_Down);
		

		const Index iUp = eRight_Up.destinationVertex;
		const Index iDown = eLeft_Down.destinationVertex;
		const Index iLeft = eUp_Left.destinationVertex;
		const Index iRight = eDown_Right.destinationVertex;

		

        mesh.faces.Reserve(4);
        const Index iULC = mesh.faces.Add();
        const Index iLDC = mesh.faces.Add();
        const Index iDRC = mesh.faces.Add();
        const Index iRUC = mesh.faces.Add();

		const Index iCenter = mesh.vertices.Add({
            Geo2D::OrthogonallyProjectPointOnLineSegment(mesh.vertices[iDown].position, mesh.vertices[iUp].position,p),
            iDRC * 4 + 1, 0, 0
        });
		
        Index ipCenter_Up = mesh.edgeInfo.Add({ iULC * 4 + 1, {} });
        Index ipCenter_Left = mesh.edgeInfo.Add({ iLDC * 4 + 1, {} });
        Index ipCenter_Right = mesh.edgeInfo.Add({ iRUC * 4 + 1, {} });
        Index ipCenter_Down = mesh.edgeInfo.Add({ iDRC * 4 + 1, {} });

        Face & fUp_Left_Center = mesh.faces[iULC];
		fUp_Left_Center.edges[0] = { iUp, iRUC * 4 + 3, eUp_Down.constrained, ipCenter_Up }; //eCenter_Up
		fUp_Left_Center.edges[1] = eUp_Left;
		fUp_Left_Center.edges[2] = {iCenter, iLDC * 4 + 1, false, ipCenter_Left}; //eLeft_Center

        Face & fLeft_Down_Center = mesh.faces[iLDC];
		fLeft_Down_Center.edges[0] = {iLeft, iULC * 4 + 3, false, ipCenter_Left}; //eCenter_Left
		fLeft_Down_Center.edges[1] = eLeft_Down;
		fLeft_Down_Center.edges[2] = {iCenter, iDRC * 4 + 1, eUp_Down.constrained, ipCenter_Down};

        Face & fDown_Right_Center = mesh.faces[iDRC];
		fDown_Right_Center.edges[0] = {iDown, iLDC * 4 + 3, eUp_Down.constrained, ipCenter_Down};
		fDown_Right_Center.edges[1] = eDown_Right;
		fDown_Right_Center.edges[2] = {iCenter, iRUC * 4 + 1, false, ipCenter_Right };

        Face & fRight_Up_Center = mesh.faces[iRUC];
		fRight_Up_Center.edges[0] = {iRight, iDRC * 4 + 3, false, ipCenter_Right };
		fRight_Up_Center.edges[1] = eRight_Up;
		fRight_Up_Center.edges[2] = {iCenter, iULC * 4 + 1, eUp_Down.constrained, ipCenter_Up};

		//Repair the opposite face references
		mesh.edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = iULC * 4 + 2;
		mesh.edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = iDRC * 4 + 2;
		mesh.edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = iLDC * 4 + 2;
		mesh.edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = iRUC * 4 + 2;

        mesh.edgeInfo[eUp_Left.edgePair].edge = iULC * 4 + 2;
        mesh.edgeInfo[eDown_Right.edgePair].edge = iDRC * 4 + 2;
        mesh.edgeInfo[eLeft_Down.edgePair].edge = iLDC * 4 + 2;
        mesh.edgeInfo[eRight_Up.edgePair].edge = iRUC * 4 + 2;
        
		mesh.vertices[iUp].edge = iULC * 4 + 1;
		mesh.vertices[iDown].edge = iDRC * 4 + 1;
		mesh.vertices[iLeft].edge = iLDC * 4 + 1;
		mesh.vertices[iRight].edge = iRUC * 4 + 1;

		if (eUp_Down.constrained) {
			mesh.vertices[iCenter].constraintCount += 2;
			//If eUp_Down is constrained we need to split the edge and add the edgePair belonging to eCenter_Up to each constraint
			EdgeInfo & info = mesh.edgeInfo[eUp_Down.edgePair];
			for (auto c : info.constraints) {
				ConstraintSegment & segment = mesh.constraints[c];
				Index index = segment.edgePairs.FindIndexLinear(eUp_Down.edgePair);
				segment.edgePairs[index] = ipCenter_Up;
                const auto & up = mesh.vertices[iUp].position;
                const auto & start = mesh.vertices[segment.startVertex].position;
                //We want to preserve the relative ordering of the edgePair array
                if(Geo2D::DistanceSquared(start - p) > Geo2D::DistanceSquared(start - up)){
                    segment.edgePairs.Insert(index + 1, ipCenter_Down);
                } else {
                    segment.edgePairs.Insert(index, ipCenter_Down);
                }
                
			}
            mesh.edgeInfo[ipCenter_Up].constraints = info.constraints;
            mesh.edgeInfo[ipCenter_Down].constraints = info.constraints;
		}
        if(centerVertex)
            *centerVertex = iCenter;
        if(edgesToCheck){
            edgesToCheck->Add(iRUC * 4 + 2); //eRight_Up
            edgesToCheck->Add(iULC * 4 + 2); //eUp_Left
            edgesToCheck->Add(iLDC * 4 + 2); //eLeft_Down
            edgesToCheck->Add(iDRC * 4 + 2); //eDown_Right
        }

		mesh.edgeInfo.Erase(eDown_Up.edgePair);
        
        mesh.faces.Erase(h / 4);
        mesh.faces.Erase(eDown_Up.oppositeHalfEdge / 4);
        
        
        
		return iCenter;
    }

	//Recycles faces in intersectedEdges before creating new faces
	static Index createConstrainedEdge(Mesh & mesh, Index segmentID, const Oryol::Array<Index> & intersectedEdges, Oryol::Array<Index> & leftBound, Oryol::Array<Index> & rightBound) {
		//As we have to create the constrained edge ourselves we have to add one to both left and right bound
		o_assert(leftBound.Size() + 1 >= 3);
		o_assert(rightBound.Size() + 1 >= 3);
		o_assert(intersectedEdges.Size() > 0);
        for(int i = 0; i < leftBound.Size(); i++){
            mesh.edgeAt(leftBound[i]).oppositeHalfEdge = -1;
        }
        for(int i = 0; i < rightBound.Size(); i++){
            mesh.edgeAt(rightBound[i]).oppositeHalfEdge = -1;
        }
        untriangulate(mesh, intersectedEdges);
		Index h = triangulate(mesh, leftBound, true);
		rightBound.Add(h);
		triangulate(mesh, rightBound, false);
        
        return TagEdgeAsConstrained(mesh, h, segmentID);
	}
    static void untriangulate(Mesh & mesh, const Oryol::Array<Index> & intersectedEdges, bool loop = false){
        //Frees faces/edge pairs associated with the edges in question
        //The number of faces that will be freed is intersectedEdges + 1
        //But number of edge pairs freed will be equal to intersected edges
		if(!loop)
			mesh.faces.Erase(mesh.edgeAt(intersectedEdges.Back()).oppositeHalfEdge / 4);

        for(Index h : intersectedEdges){
            HalfEdge & e = mesh.edgeAt(h);
            mesh.edgeInfo.Erase(e.edgePair);
            mesh.faces.Erase(h / 4);
        }
    }
	//Expects bound to be a CW list of outer edges surrounding the hole to be triangulated.
	//Triangulate handles both closed and open edge contours
	//Open contours occur when triangulating the first side of an edge pair.
    //Triangulate is fairly simplistic but works;
    //  * Picks the first edge in bound or if dealing with an open contour creates a "virtual" halfedge.
    //  * The virtual half edge doesn't have an associated opposite halfedge nor does it have an edge pair already
    //      associated with it.
    static Index triangulate(Mesh & mesh, Oryol::Array<Index> & bound, bool open) {
        
        const unsigned int edgeCount = bound.Size();
        const unsigned int firstEdge = 0;
        const unsigned int lastEdge = edgeCount - (open ? 1 : 2);
        //This is for the purposes of debugging; In order to ensure that we have a valid contour the destination
        //vertex of each contour should be the same as the origin vertex of the previous contour (CW sequence of outer edges)
        for(unsigned int i = 1; i < edgeCount; i++){
            o_assert(mesh.edgeAt(bound[i]).destinationVertex == Impl::GetOriginVertex(mesh, bound[i-1]));
        }
        if(!open){
            o_assert(mesh.edgeAt(bound.Front()).destinationVertex == Impl::GetOriginVertex(mesh, bound.Back()));
        }
        
        //Dealing with an open contour is different compared to a a closed contour
        //instead of simply using the last edge in bound we need to manually obtain the first and last vertices
        Index ivA, ivB;
        if(open){
            //We're dealing with an open contour so we have to go through more work to find our first 2 vertices
            //ivA = mesh.edgeAt(bound.Front()).destinationVertex;
			ivA = GetOriginVertex(mesh,bound.Back());
            ivB = mesh.edgeAt(bound.Front()).destinationVertex;
            o_assert(edgeCount >= 2);
        } else {

            ivA = mesh.edgeAt(bound.Back()).destinationVertex;
            ivB = GetOriginVertex(mesh, bound.Back());
            o_assert(edgeCount >= 3);
        }
        o_assert(ivA != ivB);
        
        //Most straight forward case is a triangle hole which needs to be filled
        if((open && edgeCount == 2) || (!open && edgeCount == 3)){
            Index ivC = mesh.edgeAt(bound[1]).destinationVertex;
            //else ivC = mesh.edgeAt(bound[1]).destinationVertex;
            
            //Initialise the face data
            Index ieC_B = bound[0];
            Index ieA_C = bound[1];
			Index ieB_A = open ? -1 : bound[2];
            
            HalfEdge & eA_C = mesh.edgeAt(ieA_C);
            HalfEdge & eC_B = mesh.edgeAt(ieC_B);
            
            o_assert(ivC != ivA);
            o_assert(ivC != ivB);
            o_assert(eA_C.destinationVertex == ivC);
            o_assert(eC_B.destinationVertex == ivB);
            o_assert(open || mesh.edgeAt(ieB_A).destinationVertex == ivA);
            o_assert(GetOriginVertex(mesh, ieC_B) == ivC);
            o_assert(GetOriginVertex(mesh, ieA_C) == ivA);
            o_assert(CheckFaceIsCounterClockwise(mesh,ivA,ivB,ivC));
            
            Index iA_B_C = mesh.faces.Add();
            Index ipA_B = open ? mesh.edgeInfo.Add() : mesh.edgeAt(ieB_A).edgePair;
			bool eAB_Constrained = open ? false : mesh.edgeAt(ieB_A).constrained;
            Face & fA_B_C = mesh.faces[iA_B_C];
			fA_B_C.edges[0] = { ivA, ieA_C, eA_C.constrained, eA_C.edgePair };
			fA_B_C.edges[1] = { ivB, ieB_A, eAB_Constrained, ipA_B};
			fA_B_C.edges[2] = { ivC, ieC_B, eC_B.constrained, eC_B.edgePair };
			
			//Fix up opposite half edges
			eA_C.oppositeHalfEdge = iA_B_C * 4 + 1;
			eC_B.oppositeHalfEdge = iA_B_C * 4 + 3;
			if (!open)
				mesh.edgeAt(ieB_A).oppositeHalfEdge = iA_B_C * 4 + 2;
			//Fix up edge pair data
			mesh.edgeInfo[eA_C.edgePair].edge = iA_B_C * 4 + 1;
			mesh.edgeInfo[ipA_B].edge = iA_B_C * 4 + 2;
			mesh.edgeInfo[eC_B.edgePair].edge = iA_B_C * 4 + 3;

			//Fix up vertices
			mesh.vertices[ivA].edge = iA_B_C * 4 + 1;
			mesh.vertices[ivB].edge = iA_B_C * 4 + 2;
			mesh.vertices[ivC].edge = iA_B_C * 4 + 3;
            
			return iA_B_C * 4 + 2;
          
        } else {
            Index index; //Indicates which halfedge will be used to construct the face with
            bool delaunay = true;
            const glm::dvec2 & pA = mesh.vertices[ivA].position;
            const glm::dvec2 & pB = mesh.vertices[ivB].position;
            
            //In the more complex case, a face only satisfies delaunay condition if all other vertices are outside of the face's circumcircle
            for(unsigned int i = firstEdge; i < (lastEdge-1); i++){
                delaunay = true; //Before we've checked it against the other vertices each vertex is "potentially" delaunay compliant
                index = i;
                const Index halfEdge = bound[i];
                const Index ivC = mesh.edgeAt(halfEdge).destinationVertex;
                const glm::dvec2 & pC = mesh.vertices[ivC].position;
                //As some holes can be odd shapes we only care to check vertices above the edge in question
                //Also because otherwise we'd overlap already existing faces.
                if(Geo2D::Sign(pA,pB,pC) > 0.0){
                    glm::dvec2 circumcenter = Geo2D::ComputeCircumcenter(pA, pB, pC);
                    double radiusSquared = Geo2D::DistanceSquared(circumcenter-pC) - EPSILON_SQUARED;
                    for(unsigned int j = open ? 0 : 1; j < edgeCount; j++ ){
                        Index ivD = mesh.edgeAt(bound[j]).destinationVertex;
                        const glm::dvec2 & pD = mesh.vertices[ivD].position;
                        double distanceSquared = Geo2D::DistanceSquared(pD-circumcenter);
                        if(distanceSquared < radiusSquared){
                            delaunay = false;
                            break;
                        }
                    }
                    if(delaunay) break; //Only want the first delaunay satisfying vertex.
                }
            }
            //If no triangle satisfying delaunay was found, just create a face out of the last 2 edges
            //This will usually occur with perfect n-sided polygons
            if(!delaunay){
                index = lastEdge;
            }
            //o_error("Check me");
            Index edgeA = -1, edgeB = -1;
            //Recurse into the left hole
            if(index >= (firstEdge+1)) {
                Oryol::Array<Index> boundA;
                for(Index h : bound.MakeSlice(firstEdge, index+1)){
                    boundA.Add(h);
                }
                edgeA = triangulate(mesh, boundA, true);
            }
            //Recurse into the right hole
            if(index < (lastEdge-1)){
                
                Oryol::Array<Index> boundB;
                for(Index h : bound.MakeSlice(index+1,lastEdge)){
                    boundB.Add(h);
                }
                edgeB = triangulate(mesh, boundB, true);
            }
            
            
            //Build the middle triangle -> the returned half edge needs to be trampolined up to the caller.
            Oryol::Array<Index> middleBound;
            if(index == firstEdge) {
                //No hole was found to the right so the remaining hole is only 3 sides.
                if(open){
                    middleBound = { bound.Front(), edgeB };
                } else {
                    middleBound = { bound.Front(), edgeB, bound.Back() };
                }
            } else if(index == (lastEdge-1)){
                //No hole was found to the left so the remaining hole is only 3 sides
                if(open){
                    middleBound = { edgeA, bound.Back() };
                } else {
                    middleBound = { edgeA, bound[edgeCount - 2], bound.Back() };
                }
            } else {
                //Vertex for middle triangle occurred in the middle of the edge array
                //Therefore we need to consume the halfedges returned from left hole and right hole triangulation calls.
                if(open){
                    middleBound = {edgeA,edgeB};
                } else {
                    middleBound = {edgeA,edgeB,bound.Back()};
                }
            }
            //o_error("Check me");
            return triangulate(mesh, middleBound, open);
        }
	}
    static Index TagEdgeAsConstrained(Mesh & mesh, Index h, Index segmentID){
        HalfEdge & edge = mesh.edgeAt(h);
        EdgeInfo & edgePair = mesh.edgeInfo[edge.edgePair];
        if(!edge.constrained){
            HalfEdge & opposite = mesh.edgeAt(edge.oppositeHalfEdge);
            edge.constrained = true;
            opposite.constrained = true;
            mesh.vertices[edge.destinationVertex].constraintCount += 1;
            mesh.vertices[opposite.destinationVertex].constraintCount += 1;
        }
        edgePair.constraints.Add(segmentID);
        return edge.edgePair;
    }
    
};
void Delaunay::Mesh::Setup(double width, double height)
{
    enum vIndices : Index {
        vInfinite, vBottomLeft, vBottomRight, vTopRight, vTopLeft
    };
    enum eIndices : Index {
        eBR_TL = 1, eTL_BL, eBL_BR,
        eTL_BR = 5, eBR_TR, eTR_TL,
        eInf_TL = 9, eTL_TR, eTR_Inf,
        eInf_BL = 13, eBL_TL, eTL_Inf,
        eInf_BR = 17, eBR_BL, eBL_Inf,
        eInf_TR = 21, eTR_BR, eBR_Inf
    };
    enum pIndices : Index {
      pTL_BR, pBL_TL, pBR_BL, pBR_TR, pTL_TR, pTR_Inf, pTL_Inf, pBR_Inf, pBL_Inf
    };
    enum cIndices : Index {
      cTL_BL, cBL_BR, cBR_TR, cTR_TL
    };

	boundingBox.min = { 0.0,0.0 };
	boundingBox.max = { width,height };

	//Clear everything
	vertices.Clear();
	faces.Clear();
    constraints.Clear();
    
    //Add vertices
    vertices.Add({ {width * 0.5f, height * 0.5}, eTR_Inf,0, 0});
    vertices.Add({ {0,0}, eTL_BL, 2,2 });
    vertices.Add({ {width, 0}, eBL_BR, 2, 2 });
    vertices.Add({ {width, height}, eBR_TR, 2, 2 });
    vertices.Add({ {0, height}, eTR_TL, 2, 2 });
    
    //Add faces
    faces.Add({ 0, 0, 0, {{vTopLeft,eTL_BR, false, pTL_BR}, {vBottomLeft,eBL_TL, true, pBL_TL}, {vBottomRight,eBR_BL, true, pBR_BL}}}); //fTL_BL_BR
    faces.Add({ 0, 0, 0, {{vBottomRight,eBR_TL, false, pTL_BR},{vTopRight,eTR_BR,true,pBR_TR},{ vTopLeft,eTL_TR, true, pTL_TR}}}); //fBR_TR_TL
    faces.Add({ 0, 0, 0, {{vTopLeft, eTL_Inf, false, pBL_Inf},{vTopRight,eTR_TL, true, pTL_TR}, {vInfinite, eInf_TR, false, pTR_Inf}}}); //fTL_TR_vInf
    faces.Add({ 0, 0, 0, {{vBottomLeft, eBL_Inf, false, pBL_Inf},{vTopLeft, eTL_BL, true, pBL_TL},{vInfinite, eInf_TL, false, pTL_Inf}}}); //fBL_TL_vInf
    faces.Add({ 0, 0, 0, {{vBottomRight, eBR_Inf, false, pBR_Inf},{vBottomLeft,eBL_BR, true, pBR_BL},{vInfinite, eInf_BL, false, pBL_Inf}}}); //fBR_BL_vInf
    faces.Add({ 0, 0, 0, {{vTopRight, eTR_Inf, false, pTR_Inf},{vBottomRight, eBR_TR, true, pBR_TR},{vInfinite, eInf_BR, false, pBR_Inf}}}); //fTR_BR_vInf

    //Add edge information
    edgeInfo.Add({eTL_BR, {}});
    Index ipBL_TL = edgeInfo.Add({eBL_TL, {}});
    edgeInfo[ipBL_TL].constraints.Add(cTL_BL);
    Index ipBR_BL = edgeInfo.Add({eBR_BL, {}});
    edgeInfo[ipBR_BL].constraints.Add(cBL_BR);
    Index ipBR_TR = edgeInfo.Add({eBR_TR, {}});
    edgeInfo[ipBR_TR].constraints.Add(cBR_TR);
    Index ipTL_TR = edgeInfo.Add({eTL_TR, {}});
    edgeInfo[ipTL_TR].constraints.Add(cTR_TL);
    
    edgeInfo.Add({eTR_Inf, {}});
    edgeInfo.Add({eTL_Inf, {}});
    edgeInfo.Add({eBR_Inf, {}});
    edgeInfo.Add({eBL_Inf, {}});
    
    //Add edge constraints
    constraints.Add({vTopLeft,vBottomLeft,{pBL_TL}});
    constraints.Add({vBottomLeft,vBottomRight,{pBR_BL}});
    constraints.Add({vBottomRight,vTopRight,{pBR_TR}});
    constraints.Add({vTopRight,vTopLeft,{pTL_TR}});
    
}


uint32_t Delaunay::Mesh::InsertVertex(const glm::dvec2 & p)
{
	Index centerVertex;
	Oryol::Array<Index> edgesToCheck;
    HalfEdge::Index vertex = HalfEdge::InvalidIndex;
	//Make sure the vertex is inside the bounding box the mesh was initialised with.
	//Locate the primitive the vertex falls on
	LocateRef result = this->Locate(p);
	switch (result.type) {
	case LocateRef::Vertex:
		vertex = result.object;
		break;
	case LocateRef::Edge:
		vertex = Impl::SplitEdge(*this,result.object, p, &centerVertex, &edgesToCheck);
		break;
	case LocateRef::Face:
		vertex = Impl::SplitFace(*this,result.object, p, &edgesToCheck);
		break;
    default:
            o_error("Mesh::Locate() couldn't find a primitive \n");
            
	}
	//Restore delaunay condition
	
	while (!edgesToCheck.Empty()) {
		Index h = edgesToCheck.PopFront();
		if (!edgeAt(h).constrained && !Impl::IsDelaunay(*this, h)) {
            h = Impl::FlipEdge(*this, h);
            const HalfEdge current = edgeAt(h);
            //const HalfEdge opposite = edgeAt(current.oppositeHalfEdge);
            if(current.destinationVertex == centerVertex){
                edgesToCheck.Add(Face::prevHalfEdge(h));
                edgesToCheck.Add(Face::nextHalfEdge(current.oppositeHalfEdge));
            } else {
                edgesToCheck.Add(Face::nextHalfEdge(h));
                edgesToCheck.Add(Face::prevHalfEdge(h));
            }
		}
	}
	
	return vertex;
}

uint32_t Delaunay::Mesh::InsertConstraintSegment(const glm::dvec2 & p1, const glm::dvec2 & p2){
    //Clip the vertices against the mesh's AABB
	auto clipped = Geo2D::ClipSegment(p1, p2, boundingBox);
	//Check to see if the segment is inside the bounding box and the segment has adequate length.
	if (!clipped.success || DistanceSquared(clipped.a - clipped.b) < EPSILON_SQUARED)
		return -1;
    
    Index iSegment = constraints.Add({});
    
    //Insert the first and last vertices
    ConstraintSegment & segment = constraints[iSegment];
    segment.startVertex = InsertVertex(clipped.a);
    segment.endVertex = InsertVertex(clipped.b);
    //Oryol::Log::Info("\nInserting Segment with id: %u\n",iSegment);
    const glm::dvec2 tangent = { -(clipped.b.y - clipped.a.y), clipped.b.x - clipped.a.x };

    //o_assert(Geo2D::Sign(a,b,this->vertices[segment.endVertex].position) < 0.0);
    vertices[segment.startVertex].endPointCount += 1;
    vertices[segment.endVertex].endPointCount += 1;
    
    //Sweep from the start vertex to the final vertex recording all edges we intersect
	Oryol::Array<Index> intersectedEdges, leftBound, rightBound;
    Oryol::Set<Index> visitedVertices;
    Index currentEdge, currentVertex = segment.startVertex;
	LocateRef::Code currentType = LocateRef::Vertex;
    bool done = false;
    while(true){
		done = false;
        const glm::dvec2 currentPosition = vertices[currentVertex].position;
        const glm::dvec2 tangentSegmentA = currentPosition + 0.5 * tangent;
        const glm::dvec2 tangentSegmentB = currentPosition - 0.5 * tangent;
        
		if (currentType == LocateRef::Vertex) {
            //Process vertex index
            //Oryol::Log::Info("Processing Vertex with id: %u\n",currentVertex);
            o_assert(!visitedVertices.Contains(currentVertex));
            o_assert(currentVertex != 0);
            visitedVertices.Add(currentVertex);
            const Index first = this->GetOutgoingEdgeFor(currentVertex);
            {
                Index h = first;
                do {
                    HalfEdge & edge = edgeAt(h);
                    o_assert(edge.destinationVertex != currentVertex);
                    const glm::dvec2 & vertexPosition = this->vertices[edge.destinationVertex].position;
                    //First case we check for is when the current edge is directly connected to the final vertex
                    if (edge.destinationVertex == segment.endVertex) {
                        //Oryol::Log::Info("Found final vertex\n");
                        Index pair = Impl::TagEdgeAsConstrained(*this, h, iSegment);
                        segment.edgePairs.Add(pair);

                        return iSegment; //TODO: Don't just return here because probably have to do some clean up by this point
                    }
                    
                    //Next we check if we've hit a vertex which is in approximately in line with our target vertex
                    //Also make sure we're heading in the right direction
                    if (Geo2D::DistanceSquaredPointToLineSegment(clipped.a, clipped.b, vertexPosition) <= EPSILON_SQUARED
                        /*&& Geo2D::Sign(tangentSegmentA,tangentSegmentB,vertexPosition) > 0.0*/ && edge.destinationVertex != 0) {
                        
                        o_assert(!visitedVertices.Contains(edge.destinationVertex));
                        //Oryol::Log::Info("Advance to next vertex\n");
                        Index pair = Impl::TagEdgeAsConstrained(*this, h, iSegment);
                        segment.edgePairs.Add(pair);

                        currentVertex = edge.destinationVertex;
                        currentType = LocateRef::Vertex;
                            
                        done = true;
                        break;
                        
                    }
                } while((h = this->GetNextOutgoingEdge(h)) != first);
            }
			if (done) 
				continue; //Common cases are handled so there's no sense in waiting around.
            {
                //Oryol::Log::Info("Processing adjacent intersections for vertex\n");
                Index h = first;
                //Process adjacent edge intersections
                do {
                    const Index iAdj = Face::nextHalfEdge(h);
                    HalfEdge & adjacent = edgeAt(iAdj);
                    
                    const glm::dvec2 pA = this->vertices[edgeAt(h).destinationVertex].position;
                    const glm::dvec2 pB = this->vertices[adjacent.destinationVertex].position;
                    
                    glm::dvec2 intersection;
                    if (Geo2D::ComputeIntersection(pA, pB, clipped.a, clipped.b, &intersection)
                        /*&& Geo2D::Sign(tangentSegmentA,tangentSegmentB,intersection) > 0.0*/) {
                        //Oryol::Log::Info("Intersected Segment: ")
                        //Ensure that the adjacent edge is on the correct side of the current vertex
                        if (adjacent.constrained) {
                            //o_error("Check me");
                            //If the edge we've hit is constrained we will need to split the edge and we advance the search from the newly created vertex.
                            Index newVertex = Impl::SplitEdge(*this, iAdj, intersection);
                            //Oryol::Log::Info("Created new vertex at intersection with id %u\n",newVertex);
                            o_assert(!visitedVertices.Contains(newVertex));
                            const Index first = this->GetOutgoingEdgeFor(newVertex);
                            Index h = first;
                            do {
                                HalfEdge & edge = edgeAt(h);
                                if(edge.destinationVertex == currentVertex){
                                    Index pair = Impl::TagEdgeAsConstrained(*this, h, iSegment);
                                    segment.edgePairs.Add(pair);
                                    break;
                                }
                            } while((h = this->GetNextOutgoingEdge(h)) != first);

                            currentVertex = newVertex;
                            currentType = LocateRef::Vertex;
                        }
                        else {
                            
                            intersectedEdges.Add(iAdj);
                            Index ccw = Face::nextHalfEdge(iAdj);
                            Index cw = Face::prevHalfEdge(iAdj);
                            
                            rightBound.Insert(0,edgeAt(ccw).oppositeHalfEdge);
                            leftBound.Add(edgeAt(cw).oppositeHalfEdge);
                            
                            currentEdge = adjacent.oppositeHalfEdge;
                            currentType = LocateRef::Edge;
                        }
                        done = true;
                        break;
                    }
                } while((h = this->GetNextOutgoingEdge(h)) != first);
            }
            o_assert2(done,"By this point we should have found an adjacent edge to hit\n");
		}
		else if (currentType == LocateRef::Edge) {
            //Oryol::Log::Info("Processing Edge with id: %u\n",currentEdge);
           // o_error("Implement me\n");
            //Process Edge Index
            HalfEdge & nextEdge = edgeAt(Face::nextHalfEdge(currentEdge));
            if(nextEdge.destinationVertex == segment.endVertex){
                //We've found our final vertex -> trigger triangulation
                Index cw = Face::prevHalfEdge(currentEdge);
                Index ccw = Face::nextHalfEdge(currentEdge);
				leftBound.Add( edgeAt(ccw).oppositeHalfEdge);
				rightBound.Insert(0,edgeAt(cw).oppositeHalfEdge);
                //o_error("Check me");
				Index pair = Impl::createConstrainedEdge(*this, iSegment, intersectedEdges, leftBound, rightBound);
				segment.edgePairs.Add(pair);
				return iSegment;
            } else if(Geo2D::DistanceSquaredPointToLineSegment(clipped.a, clipped.b, this->vertices[nextEdge.destinationVertex].position) <= EPSILON_SQUARED) {
                //We've hit a vertex -> trigger triangulation
				Index newVertex = nextEdge.destinationVertex;
                o_assert(!visitedVertices.Contains(newVertex));
                //o_error("Check me");
				leftBound.Add( edgeAt(Face::nextHalfEdge(currentEdge)).oppositeHalfEdge);
				rightBound.Insert(0,edgeAt(Face::prevHalfEdge(currentEdge)).oppositeHalfEdge);
				Index pair = Impl::createConstrainedEdge(*this, iSegment, intersectedEdges, leftBound, rightBound);
				segment.edgePairs.Add(pair);

				intersectedEdges.Clear();
				leftBound.Clear();
				rightBound.Clear();
				
				currentVertex = newVertex;
				currentType = LocateRef::Vertex;
            } else {
                //Oryol::Log::Info("Advancing through adjacent edges\n");
                Index cw = Face::prevHalfEdge(currentEdge);
                Index ccw = Face::nextHalfEdge(currentEdge);
                const glm::dvec2 & pA = vertices[edgeAt(currentEdge).destinationVertex].position;
                const glm::dvec2 & pB = vertices[edgeAt(cw).destinationVertex].position;
                const glm::dvec2 & pC = vertices[edgeAt(ccw).destinationVertex].position;
                o_assert(Geo2D::ComputeIntersection(clipped.a, clipped.b, pA, pB));
                //o_assert(Geo2D::Sign(pA,pB,pC) < 0);
                glm::dvec2 intersection;

                //First Test the ccw segment defined by A-C
                if (Geo2D::ComputeIntersection(clipped.a, clipped.b, pA, pC, &intersection)) {
                    //Oryol::Log::Info("Hit CCW segment: %u\n",ccw);
					HalfEdge & edge = edgeAt(ccw);
					if (edge.constrained) {
						// We've hit a constrained edge -> trigger triangulation
						Index newVertex = Impl::SplitEdge(*this, ccw, intersection);
                        //Oryol::Log::Info("Created new vertex at intersection with id %u\n",newVertex);
						o_assert(currentVertex != newVertex);
                        o_assert(!visitedVertices.Contains(newVertex));

                        const Index first = this->GetOutgoingEdgeFor(newVertex);
                        Index h = first;
                        do {
                            HalfEdge & outgoing = edgeAt(h);
                            if(outgoing.destinationVertex == Impl::GetOriginVertex(*this,leftBound.Back())){
                                //We've hit our only edge on the left hand side
                                leftBound.Add(h);
                            }
                            else if(outgoing.destinationVertex == edgeAt(rightBound.Front()).destinationVertex){
                                //o_error("Check me");
                                if(outgoing.constrained){
                                    //We've hit our final edge on the right hand side
                                    rightBound.Insert(0,outgoing.oppositeHalfEdge);
                                } else {
                                    //Hit an intermediate edge so add it to the right bound
                                    Index adj = Face::prevHalfEdge(outgoing.oppositeHalfEdge);
                                    rightBound.Insert(0,edgeAt(adj).oppositeHalfEdge);
                                    intersectedEdges.Add(outgoing.oppositeHalfEdge);
                                }
                            }
                        } while((h = this->GetNextOutgoingEdge(h)) != first);
                        //o_error("Check me");
						Index pair = Impl::createConstrainedEdge(*this, iSegment, intersectedEdges, leftBound, rightBound);
						segment.edgePairs.Add(pair);

						intersectedEdges.Clear();
						leftBound.Clear();
						rightBound.Clear();

						currentVertex = newVertex;
						currentType = LocateRef::Vertex;
					}
					else {
                        
						intersectedEdges.Add(ccw);
						rightBound.Insert(0, edgeAt(cw).oppositeHalfEdge);
						currentEdge = edgeAt(ccw).oppositeHalfEdge;
						currentType = LocateRef::Edge;

					}
                }
                else if(Geo2D::ComputeIntersection(clipped.a, clipped.b, pB, pC, &intersection)) {
                    //By the process of elimination it can only be the CW segment defined by C-B
                    //Oryol::Log::Info("Hit CW segment: %u\n",cw);
					HalfEdge & edge = edgeAt(cw);
					if (edge.constrained) {
						// We've hit a constrained edge -> trigger triangulation
                        //o_error("Check me");
						Index newVertex = Impl::SplitEdge(*this, edgeAt(cw).oppositeHalfEdge, intersection);
                        //Oryol::Log::Info("Created new vertex at intersection with id %u\n",newVertex);
                        o_assert(currentVertex != newVertex);
                        o_assert(!visitedVertices.Contains(newVertex));

                        const Index first = this->GetOutgoingEdgeFor(newVertex);
                        Index h = first;

						do {
                            
                            HalfEdge & outgoing = edgeAt(h);
                            if(outgoing.destinationVertex == edgeAt(rightBound.Front()).destinationVertex){
                                rightBound.Insert(0,outgoing.oppositeHalfEdge);
                            }
                            if(outgoing.destinationVertex == Impl::GetOriginVertex(*this, leftBound.Back())){
                                if(outgoing.constrained){
                                    leftBound.Add(h);
                                } else {
                                    HalfEdge & adj = edgeAt(Face::nextHalfEdge(h));
                                    leftBound.Add(adj.oppositeHalfEdge);
                                    intersectedEdges.Add(h);
                                }
                            }

                        } while((h = this->GetPrevOutgoingEdge(h)) != first);
                        //o_error("Check me\n");
						Index pair = Impl::createConstrainedEdge(*this, iSegment, intersectedEdges, leftBound, rightBound);
                        
						segment.edgePairs.Add(pair);

						intersectedEdges.Clear();
						leftBound.Clear();
						rightBound.Clear();

						currentVertex = newVertex;
						currentType = LocateRef::Vertex;
					}
					else {
                        //o_error("Check me");
                        

                        intersectedEdges.Add(cw);
                        leftBound.Add(edgeAt(ccw).oppositeHalfEdge);
                        
						currentEdge = edgeAt(cw).oppositeHalfEdge;
						currentType = LocateRef::Edge;
					}
					
                } else {
                    //If we didn't hit either segment this means that something has corrupted the face/halfedge data
                    o_error("Mesh::InsertConstraintSegment(): Didn't hit CW or CCW segment\n");
                }
                
            }
        } else {
            o_error("Illegal state\n");
        }
        
    }
    o_error("Something has gone wrong so eventually we should just clean up instead of asserting\n");
    return -1;
    
}


bool Delaunay::Mesh::RemoveVertex(const uint32_t vertexID)
{
	//This function handles the following cases for "permissible" vertex removal
	//vertexID must not be an end point and must either have zero or two constrained edges originating from it.
	Vertex & vertex = this->vertices[vertexID];
	
	if (vertex.endPointCount == 0) {
		if (vertex.constraintCount == 0) {
			//This case handles the completely unconstrained vertex. So we figure out the outer bounds
			//And remove the faces surrounding the vertex + finally the vertex.
			Oryol::Array<Index> bound, intersectedEdges;
			const Index first = this->GetOutgoingEdgeFor(vertexID);
			Index h = first;
			do {
				o_assert(!edgeAt(h).constrained);
				intersectedEdges.Add(h);
				Index adj = edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge;
				bound.Insert(0,adj);
			} while ((h = this->GetNextOutgoingEdge(h)) != first);
			Impl::untriangulate(*this, intersectedEdges, true);
			this->vertices.Erase(vertexID);
			Impl::triangulate(*this, bound, false);
			return true;
		}
		else if (vertex.constraintCount == 2) {
			//For this case the vertex has two constrained edges coming off it so we first have to determine those.
			//We take the naive approach and scan through all outgoing half-edges to find our two constrained half-edges.
			Oryol::Array<Index> leftBound, rightBound, intersectedEdges;
			Index hCenterUp = -1, hCenterDown = -1;
			{
				const Index first = this->GetOutgoingEdgeFor(vertexID);
				Index h = first;
				do {
					HalfEdge & outgoing = edgeAt(h);
					if (outgoing.constrained) {
						if (hCenterUp == (Index)-1) {
							hCenterUp = h;
							//rightBound.Add(edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge);
						}
						else if (hCenterDown == (Index)-1) {
							hCenterDown = h;
							//leftBound.Insert(0, edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge);
						}
						else
							o_error("The vertex has more than two constrained edges\n");
					}
				} while ((h = this->GetNextOutgoingEdge(h)) != first);
				o_assert(hCenterUp != (Index)-1 && hCenterDown != (Index)-1);
			}
			//TODO: Add an assert that verifies that the up, down and center vertices are colinear
			//Once we've identified up and down constraint edges then we can loop through the outgoing edges again, building our left and right bounds
			{
				const Index last = hCenterDown;
				Index h = hCenterUp;
				do {
					Index iAdj = Face::nextHalfEdge(h);
					intersectedEdges.Add(h);
					leftBound.Add(edgeAt(iAdj).oppositeHalfEdge);
				} while ((h = this->GetPrevOutgoingEdge(h)) != last);
			}
			{
				const Index last = edgeAt(hCenterDown).oppositeHalfEdge;
				Index h = edgeAt(hCenterUp).oppositeHalfEdge;
				do {
					HalfEdge & adjacent = edgeAt(Face::prevHalfEdge(h));
                    intersectedEdges.Add(Face::nextHalfEdge(h));
					rightBound.Insert(0,adjacent.oppositeHalfEdge);
				} while ((h = this->GetNextIncomingEdge(h)) != last);
			}
            //Before we can triangulate the hole we need to retain some information about the old edgePair
            Index ipCenterUp = edgeAt(hCenterUp).edgePair;
            Index ipCenterDown = edgeAt(hCenterDown).edgePair;
            //Naively we can assume the constraints on ipCenterUp are the same as ipCenterDown
            Oryol::Set<Index> edgeConstraints = edgeInfo[ipCenterUp].constraints;
			//Clean up our mess
			Impl::untriangulate(*this, intersectedEdges, true);
            vertices.Erase(vertexID);
            //Then we triangulate our left and right bounds, retaining a half edge from the call to triangulate.
			Index hUp_Down = Impl::triangulate(*this, leftBound, true);
			rightBound.Add(hUp_Down);
			Impl::triangulate(*this, rightBound, false);
			//Once done we set the new edge to be constrained and modify all constraints using this vertex to replace the two old edgePairs with our single new edge pair
            HalfEdge & eUp_Down = edgeAt(hUp_Down);
            eUp_Down.constrained = true;
            edgeAt(eUp_Down.oppositeHalfEdge).constrained = true;
            edgeInfo[eUp_Down.edgePair].constraints = edgeConstraints;
            for(Index cIndex : edgeConstraints){
                ConstraintSegment & segment = this->constraints[cIndex];
                int pairIndex = segment.edgePairs.FindIndexLinear(ipCenterDown);
                segment.edgePairs[pairIndex] = eUp_Down.edgePair;
                segment.edgePairs.Erase(segment.edgePairs.FindIndexLinear(ipCenterUp));
            }
            return true;
		}
	}
	return false;
}

void Delaunay::Mesh::RemoveConstraintSegment(const uint32_t constraintID){
    ConstraintSegment & segment = this->constraints[constraintID];
    //First things first; clean edge pairs associated with the constraint segment
    Oryol::Array<Index> segmentVertices {segment.startVertex};
    for(Index pairIndex : segment.edgePairs){
        EdgeInfo & edgePair = this->edgeInfo[pairIndex];
        HalfEdge & edge = edgeAt(edgePair.edge);
        HalfEdge & opposite = edgeAt(edge.oppositeHalfEdge);
        
        edgePair.constraints.Erase(constraintID);
        if(edgePair.constraints.Size() == 0){
            edge.constrained = false;
            opposite.constrained = false;
            vertices[edge.destinationVertex].constraintCount -= 1;
            vertices[opposite.destinationVertex].constraintCount -= 1;
        }
        if(segmentVertices.Back() == edge.destinationVertex){
            segmentVertices.Add(opposite.destinationVertex);
        } else {
            segmentVertices.Add(edge.destinationVertex);
        }
    }
    //Then clean up our vertices
    vertices[segment.startVertex].endPointCount -= 1;
    vertices[segment.endVertex].endPointCount -= 1;
    for(Index vIndex : segmentVertices){
        RemoveVertex(vIndex);
    }
    this->constraints.Erase(constraintID);
    
}


Delaunay::Mesh::LocateRef Delaunay::Mesh::Locate(const glm::dvec2 & p)
{
    LocateRef result { Index(-1), LocateRef::None};
	Index currentFace = -1;

	{
		Index bestVertex = HalfEdge::InvalidIndex;
		//Seed the random generator
		srand(p.x * 10 + p.y * 4);
		//Find closest vertex by randomly sampling the vertices (excluding the infinite vertex)
		int vertexSampleCount = std::pow(this->vertices.Size(), 1 / 3.);
		double minDistanceSquared = std::numeric_limits<double>::infinity();
		for (int i = 0; i < vertexSampleCount; i++) {
			Index vIndex = rand_range(1, this->vertices.Size() - 1);
            Index index = this->vertices.ActiveIndexAtIndex(vIndex);
			const Vertex & vertex = this->vertices[index];
			double distanceSquared = Geo2D::DistanceSquared(p-vertex.position);
			if (distanceSquared < minDistanceSquared) {
				minDistanceSquared = distanceSquared;
				bestVertex =index;
			}
		}
		//Start jump-and-walk search with the first face associated with the best vertex
        const Index first = this->GetOutgoingEdgeFor(bestVertex);
        Index h = first;
		do {
			if (this->faces[h / 4].isReal()) {
				currentFace = h / 4;
				break;
			}
        } while((h = this->GetNextOutgoingEdge(h)) != first);
	}
	
	Oryol::Set<Index> visitedFaces;
	int iterations = 0;
    while (!visitedFaces.Contains(currentFace) && !(result = Impl::IsInFace(*this,currentFace,p))) {
		visitedFaces.Add(currentFace);
		iterations++;
		if (iterations == 50) {
			//Log this as it is taking longer than expected
		}
		else if (iterations > 1000) {
			//Bail out if too many iterations have elapsed
            Oryol::Log::Info("Mesh::Locate({%f,%f}) has taken 1000 iterations to locate the closest primitive", p.x, p.y);
			result.type = LocateRef::None;
			break;
		}
		Index nextFace = HalfEdge::InvalidIndex;
		//Find the best direction to look in.
        for (int i = 1; i < 4; i++) {
			//Determine if the position falls to the right of the current half edge (thus outside of the current face)
            Index h = currentFace * 4 + i;
            Vertex & originVertex = this->vertices[Impl::GetOriginVertex(*this, h)];
			Vertex & destinationVertex = this->vertices[edgeAt(h).destinationVertex];
			if (Geo2D::Sign(originVertex.position, destinationVertex.position, p) < 0.0) {
				nextFace = edgeAt(h).oppositeHalfEdge / 4;
				break;
			}
		}
		if (nextFace != HalfEdge::InvalidIndex) {
			currentFace = nextFace;
		}
		else {
            Oryol::Log::Info("Something has gone wrong inside Mesh::Locate()");
			result.type = LocateRef::None;
			break; //Something has gone wrong so log it and bail
		}
	}
	
	return result;
}

inline Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) {
	return faces[index / 4].edges[(index & 3) - 1];
}





