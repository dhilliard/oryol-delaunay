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
    static Index GetOriginVertex(Mesh & mesh, Index h){
        return mesh.edgeAt(Mesh::Face::prevHalfEdge(h)).destinationVertex;
    }
    //Deletes face at specified index and EraseSwaps last face into the vacated index + patches references.
    static void DeleteFace(Mesh & mesh, Face & f){
        Index id = &f - mesh.faces.begin();
        mesh.faces.EraseSwap(id);
        Face & face = mesh.faces[id];
        for(int i = 0; i < 4; i++){
            mesh.edgeAt(face.edges[0].oppositeHalfEdge).oppositeHalfEdge = id * 4 + 1 + i;
            mesh.edgeInfo[face.edges[i].edgePair].edge = id * 4 + 1 + i;
        }
    }
    static Mesh::ObjectRef IsInFace(Mesh & mesh, const glm::dvec2 & p, Face & face)
    {
        ObjectRef result{ HalfEdge::InvalidIndex, ObjectRef::None, (size_t)-1 };
        size_t faceIndex = &face - mesh.faces.begin();
        
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
                    result.type = ObjectRef::Vertex;
                }
                else if (proximity[2]) {
                    result.object = face.edges[2].destinationVertex;
                    result.type = ObjectRef::Vertex;
                }
                else {
                    result.object = faceIndex * 4 + 1; //eV3_V1
                    result.type = ObjectRef::Edge;
                }
            }
            else if (proximity[1]) {
                if (proximity[2]) {
                    result.object = face.edges[1].destinationVertex;
                    result.type = ObjectRef::Vertex;
                }
                else {
                    result.object = faceIndex * 4 + 2; //eV1_V2
                    result.type = ObjectRef::Edge;
                }
            }
            else if (proximity[2]) {
                result.object = faceIndex * 4 + 3; //eV2_V3
                result.type = ObjectRef::Edge;
            }
            else {
                result.object = faceIndex;
                result.type = ObjectRef::Face;
            }
        }
        switch(result.type){
            case ObjectRef::Vertex:
                result.generation = mesh.vertices[result.object].generation;
                break;
            case ObjectRef::Face:
            case ObjectRef::Edge:
                result.generation = mesh.faces[faceIndex].generation;
                break;
            default:
                result.generation = (size_t)-1;
        }
        return result;
    }

	static bool IsDelaunay(Mesh & mesh, Index h)
	{
		Face & face = mesh.faces[h / 4];
		Vertex & vA = mesh.vertices[face.edges[0].destinationVertex];
		Vertex & vB = mesh.vertices[face.edges[1].destinationVertex];
		Vertex & vC = mesh.vertices[face.edges[2].destinationVertex];

		Vertex & vOpp = mesh.vertices[mesh.edgeAt(Face::nextHalfEdge(mesh.edgeAt(h).oppositeHalfEdge)).destinationVertex];


		glm::dvec2 circumcenter = Geo2D::ComputeCircumcenter(vA.position, vB.position, vC.position);
		double squaredRadius = Geo2D::DistanceSquared(vA.position - circumcenter);
		double squaredDistance = Geo2D::DistanceSquared(vOpp.position - circumcenter);
		return squaredDistance >= squaredRadius;
	}
    
    //Returns halfedge from new face pair created as a result of the flip
    static Index FlipEdge(Mesh & mesh, Index h) {
        const HalfEdge eUp_Down = mesh.edgeAt(h);

        const Index iRightFace = h / 4;
        const Index iLeftFace = eUp_Down.oppositeHalfEdge / 4;
        
        //Recycle the old faces to use them as the new faces
        Face & fLeft_Right_Up = mesh.faces[iRightFace];
        Face & fRight_Left_Down = mesh.faces[iLeftFace];
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
    
		//Construct the new faces in place of the old ones
        fLeft_Right_Up.edges[0] = {iLeft, eUp_Left.oppositeHalfEdge, eUp_Left.constrained, eUp_Left.edgePair}; //eUp_Left
        fLeft_Right_Up.edges[1] = {iRight, iLeftFace * 4 + 2, false, eUp_Down.edgePair}; //eLeft_Right
        fLeft_Right_Up.edges[2] = {iUp, eRight_Up.oppositeHalfEdge, eRight_Up.constrained, eRight_Up.edgePair}; //eRight_Up
        fLeft_Right_Up.generation++;
        
        fRight_Left_Down.edges[0] = {iRight, eDown_Right.oppositeHalfEdge, eDown_Right.constrained, eDown_Right.edgePair}; //eDown_Right
        fRight_Left_Down.edges[1] = {iLeft, iRightFace * 4 + 2, false, eUp_Down.edgePair}; //eRight_Left
        fRight_Left_Down.edges[2] = {iDown, eLeft_Down.oppositeHalfEdge, eLeft_Down.constrained, eLeft_Down.edgePair}; //eLeft_Down
        fRight_Left_Down.generation++;
        
		//Patch opposite half edge references
        mesh.edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = iRightFace * 4 + 1;
        mesh.edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = iRightFace * 4 + 3;
        
        mesh.edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = iLeftFace * 4 + 1;
        mesh.edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = iLeftFace * 4 + 3;
        
        mesh.edgeInfo[eUp_Left.edgePair].edge = iRightFace * 4 + 1;
        mesh.edgeInfo[eRight_Up.edgePair].edge = iRightFace * 4 + 3;
        mesh.edgeInfo[eDown_Right.edgePair].edge = iLeftFace * 4 + 1;
        mesh.edgeInfo[eLeft_Down.edgePair].edge = iLeftFace * 4 + 3;
        

		//Patch vertices to refer to the new half edges
		mesh.vertices[iUp].edge = iRightFace * 4 + 3;
		mesh.vertices[iDown].edge = iLeftFace * 4 + 3;
		mesh.vertices[iLeft].edge = iLeftFace * 4 + 2;
		mesh.vertices[iRight].edge = iRightFace * 4 + 2;
        
        return iRightFace * 4 + 2; //eLeft_Right
		
    }
    
    
    //Recycles one face + creates 2 new faces
    //Returns new vertex that is created as a result of splitting the face
    static Index SplitFace(Mesh & mesh, Index f, const glm::dvec2 & p, Oryol::Array<Index> * edgesToCheck = nullptr){
		const Index faceOffset = mesh.faces.Size();
		const Index pairOffset = mesh.edgeInfo.Size();

		const Face fA_B_C = mesh.faces[f];
		
		const HalfEdge & eC_A = fA_B_C.edges[0];
		const HalfEdge & eA_B = fA_B_C.edges[1];
		const HalfEdge & eB_C = fA_B_C.edges[2];

		const Index vA = eC_A.destinationVertex;
		const Index vB = eA_B.destinationVertex;
		const Index vC = eB_C.destinationVertex;
		//Create the new vertex
		const Index vCenter = mesh.vertices.Size();
		mesh.vertices.Add({ p,faceOffset * 4 + 3, 0, 0, 0 });

		//Create the new edge pair info structs
        mesh.edgeInfo.Add({ f*4 + 1, {} });
        mesh.edgeInfo.Add({ faceOffset * 4 + 1, {} });
        mesh.edgeInfo.Add({ faceOffset * 4 + 5, {} });

		mesh.faces.Reserve(2);
		Face & fC_A_Center = mesh.faces[f];
		Face & fA_B_Center = mesh.faces.Add();
		Face & fB_C_Center = mesh.faces.Add();

		fC_A_Center.edges[0] = {vC, faceOffset * 4 + 7, false, pairOffset + 0}; //eCenter_C
		fC_A_Center.edges[1] = {vA, eC_A.oppositeHalfEdge, eC_A.constrained, eC_A.edgePair}; //eC_A
		fC_A_Center.edges[2] = {vCenter, faceOffset * 4 + 1, false, pairOffset + 1}; //eA_Center
		fC_A_Center.generation++;

		fA_B_Center.matId = fA_B_C.matId;
		fA_B_Center.flags = fA_B_C.flags;
		fA_B_Center.edges[0] = {vA, f * 4 + 3, false, pairOffset + 1}; //eCenter_A
		fA_B_Center.edges[1] = {vB, eA_B.oppositeHalfEdge, eA_B.constrained, eA_B.edgePair}; //eA_B
		fA_B_Center.edges[2] = {vCenter, faceOffset * 4 + 5, false, pairOffset + 2}; //eB_Center

		fB_C_Center.matId = fA_B_C.matId;
		fB_C_Center.flags = fA_B_C.flags;
		fB_C_Center.edges[0] = {vB, faceOffset * 4 + 3, false, pairOffset + 2}; //eCenter_B
		fB_C_Center.edges[1] = {vC, eB_C.oppositeHalfEdge, eB_C.constrained, eB_C.edgePair}; //eB_C
		fB_C_Center.edges[2] = {vCenter, f * 4 + 1, false, pairOffset + 0}; //eC_Center
		
		//Patch opposite half edge references
		mesh.edgeAt(eC_A.oppositeHalfEdge).oppositeHalfEdge = f * 4 + 2;
		mesh.edgeAt(eA_B.oppositeHalfEdge).oppositeHalfEdge = faceOffset * 4 + 2;
		mesh.edgeAt(eB_C.oppositeHalfEdge).oppositeHalfEdge = faceOffset * 4 + 6;
        
        mesh.edgeInfo[eC_A.edgePair].edge = f * 4 + 2;
        mesh.edgeInfo[eA_B.edgePair].edge = faceOffset * 4 + 2;
        mesh.edgeInfo[eB_C.edgePair].edge = faceOffset * 4 + 6;

		//Patch existing vertices to refer to the new half edges
		mesh.vertices[vA].edge = f * 4 + 2;
		mesh.vertices[vB].edge = faceOffset * 4 + 2;
		mesh.vertices[vC].edge = faceOffset * 4 + 6;
		
        if(edgesToCheck){
            edgesToCheck->Add(f * 4 + 2); //eC_A
            edgesToCheck->Add(faceOffset * 4 + 2); //eA_B
            edgesToCheck->Add(faceOffset * 4 + 6); //eB_C
        }

		return vCenter;
    }
    
    //Deletes 2 faces + creates 4 faces
    //Returns new vertex that is created as a result of splitting the edge
    static Index SplitEdge(Mesh & mesh, Index h, const glm::dvec2 & p, Index * centerVertex = nullptr, Oryol::Array<Index> * edgesToCheck = nullptr) {
		const HalfEdge eUp_Down = mesh.edgeAt(h);
		const HalfEdge eDown_Up = mesh.edgeAt(eUp_Down.oppositeHalfEdge);

		//Check to see if the point is close enough to a vertex and return the corresponding halfedge
		if (Geo2D::DistanceSquared(mesh.vertices[eDown_Up.destinationVertex].position - p) <= EPSILON_SQUARED)
			return h;
		if (Geo2D::DistanceSquared(mesh.vertices[eUp_Down.destinationVertex].position - p) <= EPSILON_SQUARED)
			return eUp_Down.oppositeHalfEdge;

		const Index faceOffset = mesh.faces.Size();
		const Index pairOffset = mesh.edgeInfo.Size();

		const Index iRightFace = h / 4;
		const Index iLeftFace = eUp_Down.oppositeHalfEdge / 4;

		const Index iDown_Right = Face::nextHalfEdge(h);
		const Index iRight_Up = Face::prevHalfEdge(h);

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

		const Index iULC = iRightFace;
		const Index iLDC = iLeftFace;
		const Index iDRC = faceOffset;
		const Index iRUC = faceOffset + 1;

		const Index iCenter = mesh.vertices.Size();
		mesh.vertices.Add({ Geo2D::OrthogonallyProjectPointOnLineSegment(mesh.vertices[iDown].position, mesh.vertices[iUp].position,p), iULC * 4 + 3, 0, 0, 0});
		
        mesh.edgeInfo.Add({ iULC * 4 + 1, {} });
        mesh.edgeInfo.Add({ iLDC * 4 + 1, {} });
        mesh.edgeInfo.Add({ iRUC * 4 + 1, {} });
		mesh.edgeInfo[eUp_Down.edgePair].edge = iDRC * 4 + 1;

		mesh.faces.Reserve(2);
		Face & fUp_Left_Center = mesh.faces[iRightFace];
		Face & fLeft_Down_Center = mesh.faces[iLeftFace];
		Face & fDown_Right_Center = mesh.faces.Add();
		Face & fRight_Up_Center = mesh.faces.Add();

		fUp_Left_Center.edges[0] = { iUp, iRUC * 4 + 3, eUp_Down.constrained, pairOffset + 0 }; //eCenter_Up
		fUp_Left_Center.edges[1] = {iLeft, eUp_Left.oppositeHalfEdge, eUp_Left.constrained, eUp_Left.edgePair};
		fUp_Left_Center.edges[2] = {iCenter, iLDC * 4 + 1, false, pairOffset + 1}; //eLeft_Center
		fUp_Left_Center.generation++;

		fLeft_Down_Center.edges[0] = {iLeft, iULC * 4 + 3, false, pairOffset + 1}; //eCenter_Left
		fLeft_Down_Center.edges[1] = {iDown, eLeft_Down.oppositeHalfEdge, eLeft_Down.constrained, eLeft_Down.edgePair};
		fLeft_Down_Center.edges[2] = {iCenter, iDRC * 4 + 1, eUp_Down.constrained, eUp_Down.edgePair};
		fLeft_Down_Center.generation++;

		fDown_Right_Center.edges[0] = {iDown, iLDC * 4 + 3, eUp_Down.constrained, eUp_Down.edgePair};
		fDown_Right_Center.edges[1] = {iRight, eDown_Right.oppositeHalfEdge, eDown_Right.constrained, eDown_Right.edgePair};
		fDown_Right_Center.edges[2] = {iCenter, iRUC * 4 + 1, false, pairOffset + 2 };

		fRight_Up_Center.edges[0] = {iRight, iDRC * 4 + 3, false, pairOffset + 2 };
		fRight_Up_Center.edges[1] = {iUp, eRight_Up.oppositeHalfEdge, eRight_Up.constrained, eRight_Up.edgePair};
		fRight_Up_Center.edges[2] = {iCenter, iULC * 4 + 1, eUp_Down.constrained, pairOffset + 0};

		//Repair the opposite face references
		mesh.edgeAt(eUp_Left.oppositeHalfEdge).oppositeHalfEdge = iULC * 4 + 2;
		mesh.edgeAt(eDown_Right.oppositeHalfEdge).oppositeHalfEdge = iDRC * 4 + 2;
		mesh.edgeAt(eLeft_Down.oppositeHalfEdge).oppositeHalfEdge = iLDC * 4 + 2;
		mesh.edgeAt(eRight_Up.oppositeHalfEdge).oppositeHalfEdge = iRUC * 4 + 2;

        mesh.edgeInfo[eUp_Left.edgePair].edge = iULC * 4 + 2;
        mesh.edgeInfo[eDown_Right.edgePair].edge = iDRC * 4 + 2;
        mesh.edgeInfo[eLeft_Down.edgePair].edge = iLDC * 4 + 2;
        mesh.edgeInfo[eRight_Up.edgePair].edge = iRUC * 4 + 2;
        
		mesh.vertices[iUp].edge = eUp_Left.oppositeHalfEdge;
		mesh.vertices[iDown].edge = eDown_Right.oppositeHalfEdge;
		mesh.vertices[iLeft].edge = eLeft_Down.oppositeHalfEdge;
		mesh.vertices[iRight].edge = eRight_Up.oppositeHalfEdge;

		if (eUp_Down.constrained) {
			mesh.vertices[iCenter].constraintCount += 2;
			//If eUp_Down is constrained we need to split the edge and add the edgePair belonging to eCenter_Up to each constraint
			EdgeInfo & info = mesh.edgeInfo[eUp_Down.edgePair];
			for (auto c : info.constraints) {
				ConstraintSegment & segment = mesh.constraints[c];
				Index index = segment.edgePairs.FindIndexLinear(eUp_Down.edgePair);
				segment.edgePairs.Insert(index+1, pairOffset);
			}
            mesh.edgeInfo[pairOffset].constraints = info.constraints;
		}
        if(centerVertex)
            *centerVertex = iCenter;
        if(edgesToCheck){
            edgesToCheck->Add(iRUC * 4 + 2); //eRight_Up
            edgesToCheck->Add(iULC * 4 + 2); //eUp_Left
            edgesToCheck->Add(iLDC * 4 + 2); //eLeft_Down
            edgesToCheck->Add(iDRC * 4 + 2); //eDown_Right
        }
		return iCenter;
    }
	//Returns the edge pair index associated with the new constrained edge
	//Recycles faces in intersectedEdges before creating new faces
	static Index createConstrainedEdge(Mesh & mesh, const Oryol::Array<Index> & intersectedEdges, Oryol::Array<Index> & leftBound, Oryol::Array<Index> & rightBound) {
		//As we have to create the constrained edge ourselves we have to add one to both left and right bound
		o_assert(leftBound.Size() + 1 >= 3);
		o_assert(rightBound.Size() + 1 >= 3);
		o_assert(intersectedEdges.Size() > 0);
        
		Index h = triangulate(mesh, leftBound, true);
		rightBound.Add(h);
		triangulate(mesh, rightBound, false);
        untriangulate(mesh, intersectedEdges);
		HalfEdge & edge = mesh.edgeAt(h);
		edge.constrained = true;
		mesh.edgeAt(edge.oppositeHalfEdge).constrained = true;
		return edge.edgePair;
	}
    static void untriangulate(Mesh & mesh, const Oryol::Array<Index> & intersectedEdges){
        //Frees faces/edge pairs associated with the edges in question
        //The number of faces that will be freed is intersectedEdges + 1
        //But number of edge pairs freed will be equal to intersected edges
        //TODO: Implement me.
        
    }
	//Expects bound to be a CW list of outer edges surrounding the hole to be triangulated.
	//Triangulate handles both closed and open edge contours
	//Open contours occur when triangulating the first side of an edge pair.
    //Triangulate is fairly simplistic but works;
    //  * Picks the first edge in bound or if dealing with an open contour creates a "virtual" halfedge.
    //  * The virtual half edge doesn't have an associated opposite halfedge nor does it have an edge pair already
    //      associated with it.
    static Index triangulate(Mesh & mesh, Oryol::Array<Index> & bound, bool open) {
        //Dealing with an open contour is different compared to a a closed contour
        //instead of simply using the first edge in bound we need to manually obtain the first and last vertices
        const unsigned int edgeCount = bound.Size();
        Index ivA, ivB;
        if(open){
            //We're dealing with an open contour so we have to go through more work to find our first 2 vertices
            ivA = GetOriginVertex(mesh,bound.Back());
            ivB = mesh.edgeAt(bound.Front()).destinationVertex;
        } else {
            ivA = mesh.edgeAt(bound.Back()).destinationVertex;
            ivB = mesh.edgeAt(bound.Front()).destinationVertex;
        }
        o_assert(ivA != ivB);
        
        //Most straight forward case is a triangle hole which needs to be filled
        if((open && edgeCount == 2) || (!open && edgeCount == 3)){
            
            Index ivC = mesh.edgeAt(bound[1]).destinationVertex;
            o_assert(ivC != ivA);
            o_assert(ivC != ivB);
            o_assert(Geo2D::CounterClockwise(mesh.vertices[ivA].position, mesh.vertices[ivB].position, mesh.vertices[ivC].position));
            //Initialise the face data
            HalfEdge & eA_C = mesh.edgeAt(bound[1]);
            HalfEdge & eC_B = mesh.edgeAt(bound[0]);
            
            const Index offset = mesh.faces.Size();
            
            Index pair;
            if(open){
                pair = mesh.edgeInfo.Size();
                mesh.edgeInfo.Add({offset * 4 + 2,{}});
            } else {
                pair = mesh.edgeAt(bound.Back()).edgePair;
                mesh.edgeInfo[pair].edge = offset * 4 + 2;
            }
            
            HalfEdge eC_A = {ivA, bound[1], eA_C.constrained, eA_C.edgePair};
            HalfEdge eA_B = {ivB, open ? Index(-1) : bound[2], false, pair}; //Possibly open edge of the contour
            HalfEdge eB_C = {ivC, bound[0], eC_B.constrained, eC_B.edgePair};
            
            mesh.faces.Add({0,0,0,{eC_A, eA_B, eB_C}});
            
            mesh.edgeAt(bound[0]).oppositeHalfEdge = offset * 4 + 3;
            mesh.edgeAt(bound[1]).oppositeHalfEdge = offset * 4 + 1;
            if(!open){
                mesh.edgeAt(bound[2]).oppositeHalfEdge = offset * 4 + 2;
            }
            
            mesh.vertices[ivA].edge = bound[1];
            mesh.vertices[ivB].edge = open ? -1 : bound[2];
            mesh.vertices[ivC].edge = bound[0];
            
            return offset * 4 + 2;
        
        } else {
            Index index; //Indicates which halfedge will be used to construct the face with
            bool delaunay = true;
            const glm::dvec2 & pA = mesh.vertices[ivA].position;
            const glm::dvec2 & pB = mesh.vertices[ivB].position;
            
            //In the more complex case, a face only satisfies delaunay condition if all other vertices are outside of the face's circumcircle
            for(unsigned int i = open ? 1 : 2; i < edgeCount; i++){
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
            //If no triangle satisfying delaunay was found, just create a face out of the first 2 edges
            //This will usually occur with perfect n-sided polygons
            if(!delaunay){
                index = open ? 1 : 2;
            }
            Index left = -1, right = -1;
            //Recurse into the left hole
            if(index > 1) {
                o_error("Implement me");
            }
            //Recurse into the right hole
            if(index < edgeCount - 1){
                o_error("Implement me");
            }
            //Build the middle triangle -> the returned half edge needs to be trampolined up to the caller.
            Oryol::Array<Index> middleBound;
            if((open && index == 1) || (open && index == 2)){
                o_error("Implement me");
            } else if(index == (edgeCount-1)){
                o_error("Implement me");
            } else {
                //Vertex for middle triangle occurred in the middle of the edge array
                //Therefore we need to consume the halfedges returned from left hole and right hole triangulation calls.
                o_error("Implement me");
            }
            return triangulate(mesh, middleBound, open);
        }
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
    vertices.Add({ {width * 0.5f, height * 0.5}, eTR_Inf,0,0, 0 });
    vertices.Add({ {0,0}, eTL_BL, 2,2, 0 });
    vertices.Add({ {width, 0}, eBL_BR, 2, 2, 0 });
    vertices.Add({ {width, height}, eBR_TR, 2, 2, 0 });
    vertices.Add({ {0, height}, eTR_TL, 2, 2, 0 });
    
    //Add faces
    faces.Add({ 0, 0, 0, {{vTopLeft,eTL_BR, false, pTL_BR}, {vBottomLeft,eBL_TL, true, pBL_TL}, {vBottomRight,eBR_BL, true, pBR_BL}}}); //fTL_BL_BR
    faces.Add({ 0, 0, 0, {{vBottomRight,eBR_TL, false, pTL_BR},{vTopRight,eTR_BR,true,pBR_TR},{ vTopLeft,eTL_TR, true, pTL_TR}}}); //fBR_TR_TL
    faces.Add({ 0, 0, 0, {{vTopLeft, eTL_Inf, false, pBL_Inf},{vTopRight,eTR_TL, true, pTL_TR}, {vInfinite, eInf_TR, false, pTR_Inf}}}); //fTL_TR_vInf
    faces.Add({ 0, 0, 0, {{vBottomLeft, eBL_Inf, false, pBL_Inf},{vTopLeft, eTL_BL, true, pBL_TL},{vInfinite, eInf_TL, false, pTL_Inf}}}); //fBL_TL_vInf
    faces.Add({ 0, 0, 0, {{vBottomRight, eBR_Inf, false, pBR_Inf},{vBottomLeft,eBL_BR, true, pBR_BL},{vInfinite, eInf_BL, false, pBL_Inf}}}); //fBR_BL_vInf
    faces.Add({ 0, 0, 0, {{vTopRight, eTR_Inf, false, pTR_Inf},{vBottomRight, eBR_TR, true, pBR_TR},{vInfinite, eInf_BR, false, pBR_Inf}}}); //fTR_BR_vInf

    //Add edge information
    edgeInfo.Add({eTL_BR, {}});
    edgeInfo.Add({eBL_TL, {cTL_BL}});
    edgeInfo.Add({eBR_BL, {cBL_BR}});
    edgeInfo.Add({eBR_TR, {cBR_TR}});
    edgeInfo.Add({eTL_TR, {cTR_TL}});
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


size_t Delaunay::Mesh::InsertVertex(const glm::dvec2 & p)
{
	Index centerVertex;
	Oryol::Array<Index> edgesToCheck;
	vertices.Reserve(1);
    HalfEdge::Index vertex = HalfEdge::InvalidIndex;
	//Make sure the vertex is inside the bounding box the mesh was initialised with.
	//Locate the primitive the vertex falls on
	ObjectRef result = this->Locate(p);
	switch (result.type) {
	case ObjectRef::Vertex:
		vertex = result.object;
		break;
	case ObjectRef::Edge:
		vertex = Impl::SplitEdge(*this,result.object, p, &centerVertex, &edgesToCheck);
		break;
	case ObjectRef::Face:
		vertex = Impl::SplitFace(*this,result.object, p, &edgesToCheck);
		break;
	}
	//Restore delaunay condition
	
	while (!edgesToCheck.Empty()) {
		Index h = edgesToCheck.PopFront();
		if (!edgeAt(h).constrained && !Impl::IsDelaunay(*this,h)) {
            h = Impl::FlipEdge(*this, h);
            const HalfEdge current = edgeAt(h);
            const HalfEdge opposite = edgeAt(current.oppositeHalfEdge);
            if(opposite.destinationVertex == centerVertex){
                edgesToCheck.Add(Face::nextHalfEdge(h));
                edgesToCheck.Add(Face::prevHalfEdge(h));
            } else {
                edgesToCheck.Add(Face::nextHalfEdge(current.oppositeHalfEdge));
                edgesToCheck.Add(Face::prevHalfEdge(current.oppositeHalfEdge));
            }
		}
	}
	
	return vertex;
}

size_t Delaunay::Mesh::InsertConstraintSegment(const glm::dvec2 & p1, const glm::dvec2 & p2){
    //Clip the vertices against the mesh's AABB
	auto clipped = Geo2D::ClipSegment(p1, p2, boundingBox);
	//Check to see if the segment is inside the bounding box and the segment has adequate length.
	if (!clipped.success || DistanceSquared(clipped.a - clipped.b) < EPSILON_SQUARED)
		return -1;

	size_t iSegment = constraints.Size();
    
    //Insert the first and last vertices
    ConstraintSegment & segment = constraints.Add();
    segment.startVertex = InsertVertex(clipped.a);
    segment.endVertex = InsertVertex(clipped.b);
    
    vertices[segment.startVertex].constraintCount += 1;
    vertices[segment.startVertex].endPointCount += 1;
    
    vertices[segment.endVertex].constraintCount += 1;
    vertices[segment.endVertex].endPointCount += 1;
    
    //Sweep from the start vertex to the final vertex recording all edges we intersect
	Oryol::Array<Index> intersectedEdges, leftBound, rightBound, segmentVertices = { segment.startVertex };
    Index currentEdge, currentVertex = segment.startVertex;
	ObjectRef::Code currentType = ObjectRef::Vertex;
    bool done = false;
    while(true){
		done = false;
		if (currentType == ObjectRef::Vertex) {
            //Process vertex index
			for (Index h : this->vertices[currentVertex].OutgoingEdges(*this)) {
				HalfEdge & edge = edgeAt(h);
				const glm::dvec2 & vertexPosition = this->vertices[edge.destinationVertex].position;
				//First case we check for is when the current vertex is directly connected to the final vertex
				if (edge.destinationVertex == segment.endVertex) {
					//Oryol::Log::Info("Found final vertex\n");
					segment.edgePairs.Add((Index)edge.edgePair);
                    edgeInfo[edge.edgePair].constraints.Add(iSegment);
					edge.constrained = true;
					edgeAt(edge.oppositeHalfEdge).constrained = true;
					return iSegment; //TODO: Don't just return here because probably have to do some clean up by this point
				}
				//Next we check if we've hit a vertex which is in approximately in line with our target vertex
				//Also make sure we haven't hit the previous vertex otherwise we'll get stuck in an infinite loop
				if (segmentVertices.Back() != edge.destinationVertex && Geo2D::DistanceSquaredPointToLineSegment(clipped.a, clipped.b, vertexPosition) <= EPSILON_SQUARED) {
					//Oryol::Log::Info("Advance to next vertex\n");
					segmentVertices.Add(edge.destinationVertex);
					segment.edgePairs.Add((Index)edge.edgePair);
                    edgeInfo[edge.edgePair].constraints.Add(iSegment);
                    
					edge.constrained = true;
					edgeAt(edge.oppositeHalfEdge).constrained = true;
                    
					this->vertices[edge.destinationVertex].constraintCount += 1;
                    
					currentVertex = edge.destinationVertex;
                    currentType = ObjectRef::Vertex;
                    
					done = true;
					break;
				}
			}
			if (done) continue; //Common cases are handled so there's no sense in waiting around.
			//Process adjacent edge intersections
			for (Index h : this->vertices[currentVertex].OutgoingEdges(*this)) {
				const Index iAdj = Face::nextHalfEdge(h);
				HalfEdge & adjacent = edgeAt(iAdj);
				const glm::dvec2 pA = this->vertices[edgeAt(h).destinationVertex].position;
				const glm::dvec2 pB = this->vertices[adjacent.destinationVertex].position;
				glm::dvec2 intersection;
				if (Geo2D::ComputeIntersection(pA, pB, clipped.a, clipped.b, &intersection)) {
					if (adjacent.constrained) {
						//If the edge we've hit is constrained we will need to split the edge and we advance the search from the newly created vertex.
                        Index newVertex = Impl::SplitEdge(*this, iAdj, intersection);
                        for(Index h : this->vertices[newVertex].OutgoingEdges(*this)){
                            HalfEdge & edge = edgeAt(h);
                            if(edge.destinationVertex == currentVertex){
                                
                                segment.edgePairs.Add((Index)edge.edgePair);
                                edgeInfo[edge.edgePair].constraints.Add(iSegment);
                                
                                edge.constrained = true;
                                edgeAt(edge.oppositeHalfEdge).constrained = true;
                                
                                this->vertices[edge.destinationVertex].constraintCount += 1;
                                
                                break;
                            }
                        }
                        currentVertex = newVertex;
                        currentType = ObjectRef::Vertex;
					}
					else {                        
						intersectedEdges.Add(iAdj);
						leftBound.Insert(0,edgeAt(Face::nextHalfEdge(iAdj)).oppositeHalfEdge);
						rightBound.Add(edgeAt(Face::prevHalfEdge(iAdj)).oppositeHalfEdge);
						currentEdge = adjacent.oppositeHalfEdge;
						currentType = ObjectRef::Edge;
					}
					done = true;
					break;
				}
			}
		}
		else if (currentType == ObjectRef::Edge) {
            //Process Edge Index
            HalfEdge & nextEdge = edgeAt(Face::nextHalfEdge(currentEdge));
            if(nextEdge.destinationVertex == segment.endVertex){
                //We've found our final vertex -> trigger triangulation
				leftBound.Insert(0, edgeAt(Face::prevHalfEdge(currentEdge)).oppositeHalfEdge);
				rightBound.Add(edgeAt(Face::nextHalfEdge(currentEdge)).oppositeHalfEdge);
				Index newSegment = Impl::createConstrainedEdge(*this, intersectedEdges, leftBound, rightBound);
				segment.edgePairs.Add(newSegment);
				edgeInfo[newSegment].constraints.Add(iSegment);
				return iSegment;
            } else if(Geo2D::DistanceSquaredPointToLineSegment(clipped.a, clipped.b, this->vertices[nextEdge.destinationVertex].position) <= EPSILON_SQUARED) {
                //We've hit a vertex -> trigger triangulation
				Index newVertex = nextEdge.destinationVertex;
				leftBound.Insert(0, edgeAt(Face::prevHalfEdge(currentEdge)).oppositeHalfEdge);
				rightBound.Add(edgeAt(Face::nextHalfEdge(currentEdge)).oppositeHalfEdge);
				Index newSegment = Impl::createConstrainedEdge(*this, intersectedEdges, leftBound, rightBound);
				segment.edgePairs.Add(newSegment);
				edgeInfo[newSegment].constraints.Add(iSegment);
				
				intersectedEdges.Clear();
				leftBound.Clear();
				rightBound.Clear();
				
				segmentVertices.Add(currentVertex);
				currentVertex = newVertex;
				currentType = ObjectRef::Vertex;
            } else {
                Index cw = Face::prevHalfEdge(currentEdge);
                Index ccw = Face::nextHalfEdge(currentEdge);
                const glm::dvec2 & pA = vertices[edgeAt(cw).destinationVertex].position;
                const glm::dvec2 & pB = vertices[edgeAt(ccw).destinationVertex].position;
                const glm::dvec2 & pC = vertices[edgeAt(currentEdge).destinationVertex].position;
                glm::dvec2 intersection;

                //First Test the cw segment defined by B-A
                if (Geo2D::ComputeIntersection(clipped.a, clipped.b, pB, pA, &intersection)) {
					HalfEdge & edge = edgeAt(cw);
					if (edge.constrained) {
						// We've hit a constrained edge -> trigger triangulation
						Index newVertex = Impl::SplitEdge(*this, cw, intersection);
						for (Index h : this->vertices[newVertex].OutgoingEdges(*this)) {
							HalfEdge & edge = edgeAt(h);
							HalfEdge & opposite = edgeAt(edge.oppositeHalfEdge);
							if (edge.destinationVertex == edgeAt(leftBound.Front()).destinationVertex) {
								leftBound.Insert(0, edge.oppositeHalfEdge);
							}
							if (edge.destinationVertex == edgeAt(rightBound.Back()).destinationVertex) {
								rightBound.Add(edge.oppositeHalfEdge);
							}
						}
						Index newSegment = Impl::createConstrainedEdge(*this, intersectedEdges, leftBound, rightBound);
						segment.edgePairs.Add(newSegment);
						edgeInfo[newSegment].constraints.Add(iSegment);
						intersectedEdges.Clear();
						leftBound.Clear();
						rightBound.Clear();

						segmentVertices.Add(currentVertex);
						currentVertex = newVertex;
						currentType = ObjectRef::Vertex;
					}
					else {
						intersectedEdges.Add(cw);
						leftBound.Insert(0, edgeAt(Face::nextHalfEdge(cw)).oppositeHalfEdge);
						currentEdge = edge.oppositeHalfEdge;
						currentType = ObjectRef::Edge;
					}
                }
                else {
                    //By the process of elimination it can only be the CCW segment defined by C-B
   
					HalfEdge & edge = edgeAt(ccw);
					if (edge.constrained) {
						// We've hit a constrained edge -> trigger triangulation
						Geo2D::ComputeIntersection(clipped.a, clipped.b, pC, pB, &intersection);
						Index newVertex = Impl::SplitEdge(*this, ccw, intersection);
						for (Index h : this->vertices[newVertex].OutgoingEdges(*this)) {
							HalfEdge & edge = edgeAt(h);
							HalfEdge & opposite = edgeAt(edge.oppositeHalfEdge);
							if (edge.destinationVertex == edgeAt(leftBound.Front()).destinationVertex) {
								leftBound.Insert(0, edge.oppositeHalfEdge);
							}
							if (edge.destinationVertex == edgeAt(rightBound.Back()).destinationVertex) {
								rightBound.Add(edge.oppositeHalfEdge);
							}
						}
						Index newSegment = Impl::createConstrainedEdge(*this, intersectedEdges, leftBound, rightBound);
						segment.edgePairs.Add(newSegment);
						edgeInfo[newSegment].constraints.Add(iSegment);

						intersectedEdges.Clear();
						leftBound.Clear();
						rightBound.Clear();

						segmentVertices.Add(currentVertex);
						currentVertex = newVertex;
						currentType = ObjectRef::Vertex;
					}
					else {
						intersectedEdges.Add(ccw);
						rightBound.Add(edgeAt(Face::prevHalfEdge(ccw)).oppositeHalfEdge);
						currentEdge = edge.oppositeHalfEdge;
						currentType = ObjectRef::Edge;
					}
					
                }
                
            }
        } else {
            o_error("Illegal state");
        }
        
    }
    
    return iSegment;
    
}

/*
bool Delaunay::Mesh::DeleteVertex(Index index)
{
	// Few typical cases
	//If vertex.constraintCount == 0 it is completely free of constraint --> Remove immediately
	//If vertex.constraintCount > 0 && vertex.endPointCount != 0 it is the end point of a constraint --> Cannot be removed until the constraint is broken first
	//If vertex.constraintCount == 2 && vertex.endPointCount == 0 it is the midpoint of a constraint and is no longer needed
	Vertex & vertex = this->vertices[index];
	Oryol::Array<Index> bound; //Outer bound is required by untriangulate
	Oryol::Array<Index> taggedFaces; //Contains the face indices that will be removed from the mesh
	if (vertex.constraintCount == 0 && vertex.endPointCount == 0) {
		for (Index h : vertex.OutgoingEdges(*this)) {
			bound.Add(edgeAt(Face::nextHalfEdge(h)).oppositeHalfEdge);
			taggedFaces.Add(h / 4);
		}
	}
	else if (vertex.constraintCount == 2 && vertex.endPointCount == 0) {

	}

	return false;
}

*/

Delaunay::Mesh::ObjectRef Delaunay::Mesh::Locate(const glm::dvec2 & p)
{
    ObjectRef result { size_t(-1), ObjectRef::None, size_t(-1) };
	Index currentFace = -1;

	{
		Index bestVertex = HalfEdge::InvalidIndex;
		//Seed the random generator
		srand(p.x * 10 + p.y * 4);
		//Find closest vertex by randomly sampling the vertices (excluding the infinite vertex)
		int vertexSampleCount = std::pow(this->vertices.Size(), 1 / 3.);
		double minDistanceSquared = std::numeric_limits<double>::infinity();
		for (int i = 0; i < vertexSampleCount; i++) {
			Index index = rand_range(1, this->vertices.Size() - 1);
			Vertex & vertex = this->vertices[index];
			double distanceSquared = Geo2D::DistanceSquared(p-vertex.position);
			if (distanceSquared < minDistanceSquared) {
				minDistanceSquared = distanceSquared;
				bestVertex = index;
			}
		}
		//Start jump-and-walk search with the first face associated with the best vertex
		for (Index h : this->vertices[bestVertex].OutgoingEdges(*this)) {
			if (this->faces[h / 4].isReal()) {
				currentFace = h / 4;
				break;
			}
		}
	}
	
	Oryol::Set<Index> visitedFaces;
	int iterations = 0;
	while (visitedFaces.Contains(currentFace) || !(result = Impl::IsInFace(*this,p, this->faces[currentFace]))) {
		visitedFaces.Add(currentFace);
		iterations++;
		if (iterations == 50) {
			//Log this as it is taking longer than expected
		}
		else if (iterations > 1000) {
			//Bail out if too many iterations have elapsed
			result.type = ObjectRef::None;
			break;
		}
		Index nextFace = HalfEdge::InvalidIndex;
		//Find the best direction to look in.
		for (HalfEdge & h : this->faces[currentFace].edges) {
			//Determine if the position falls to the right of the current half edge (thus outside of the current face)
			Vertex & originVertex = this->vertices[edgeAt(h.oppositeHalfEdge).destinationVertex];
			Vertex & destinationVertex = this->vertices[h.destinationVertex];
			if (Geo2D::Sign(originVertex.position, destinationVertex.position, p) < 0) {
				nextFace = h.oppositeHalfEdge / 4;
				break;
			}
		}
		if (nextFace != HalfEdge::InvalidIndex) {
			currentFace = nextFace;
		}
		else {
			result.type = ObjectRef::None;
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
    for (int i = 1; i < vertices.Size(); i++) {
		Vertex & vertex = this->vertices[i];
		debugDraw->DrawVertex(vertex.position);
		for (Index h : vertex.IncomingEdges(*this)) {
			HalfEdge & current = edgeAt(h);
            o_assert(current.destinationVertex == (Index)i);
			//Render Edges
			if ((h < current.oppositeHalfEdge) && (current.destinationVertex != 0) && (edgeAt(current.oppositeHalfEdge).destinationVertex != 0)) {
				Vertex & origin = this->vertices[current.destinationVertex];
				Vertex & destination = this->vertices[edgeAt(current.oppositeHalfEdge).destinationVertex];
				debugDraw->DrawEdge(origin.position, destination.position, current.constrained);
			}
		}
	}
}

inline Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) {
	return *(reinterpret_cast<HalfEdge*>(&faces[0])+ index);
}
inline const Delaunay::Mesh::HalfEdge & Delaunay::Mesh::edgeAt(Index index) const {
	return *(reinterpret_cast<const HalfEdge*>(&faces[0]) + index);
}

inline Index Delaunay::Mesh::Face::nextHalfEdge(Index h) {
	++h;
	return (h & 3) != 0 ? h : h - 3;
}

inline Index Delaunay::Mesh::Face::prevHalfEdge(Index h) {
	--h;
	return (h & 3) != 0 ? h : h + 3;
}


Index Delaunay::Mesh::HalfEdge::IncomingHalfEdgeIterator::operator*()
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
{}

Index Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::operator*()
{
	return current;
}

Delaunay::Mesh::HalfEdge::OutgoingHalfEdgeIterator::OutgoingHalfEdgeIterator(Mesh & mesh, Index first)
	: mesh(mesh), current(first), first(first)
{

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

Mesh::HalfEdge::IncomingHalfEdgeIterator Mesh::Vertex::IncomingEdges(Mesh & mesh) {
    Index id = this - mesh.vertices.begin();
	Index first = edge;
	if (mesh.edgeAt(edge).destinationVertex != id)
		first = mesh.edgeAt(edge).oppositeHalfEdge;
	return HalfEdge::IncomingHalfEdgeIterator(mesh, first);
}
Mesh::HalfEdge::OutgoingHalfEdgeIterator Mesh::Vertex::OutgoingEdges(Mesh & mesh) {
    Index id = this - mesh.vertices.begin();
	Index first = edge;
	if (mesh.edgeAt(edge).destinationVertex == id)
		first = mesh.edgeAt(edge).oppositeHalfEdge;
	return HalfEdge::OutgoingHalfEdgeIterator(mesh, first);
}
