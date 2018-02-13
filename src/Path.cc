
#include "Path.h"
#include "Mesh.h"
using namespace Delaunay;
//This function mainly ensures that there is sufficient space through the adjacent face to ensure the circle
//representing the agent can make it through.
bool IsEdgeWalkable(Mesh & mesh, uint32_t hFrom, uint32_t throughFace, uint32_t hTo, const double diameterSquared){
    const Mesh::HalfEdge & eTo = mesh.EdgeAt(hTo);
    const Mesh::HalfEdge & eToOpp = mesh.EdgeAt(eTo.oppositeHalfEdge);
    const Mesh::HalfEdge & eFrom = mesh.EdgeAt(hFrom);
    const Mesh::HalfEdge & eFromOpp = mesh.EdgeAt(eFrom.oppositeHalfEdge);
    uint32_t ivA, ivB, ivC, adjacent;
    //For the edges to be walkable both edges must have adequate clearance for the circle representing the agent to pass through
    //There are a couple typical cases. For acute or right angled triangles we can directly check the length of each edge is greater than diameterSquared
    //Obtuse triangles require more work to assess whether the agent circle can go through hFrom through to hTo
    if(eTo.destinationVertex == eFrom.destinationVertex){
        ivB = eToOpp.destinationVertex;
        ivA = eFromOpp.destinationVertex;
        ivC = eTo.destinationVertex;
        adjacent = Mesh::Face::prevHalfEdge(hFrom);
        
    } else if(eToOpp.destinationVertex == eFromOpp.destinationVertex){
        ivB = eTo.destinationVertex;
        ivA = eFrom.destinationVertex;
        ivC = eToOpp.destinationVertex;
        adjacent = Mesh::Face::nextHalfEdge(hFrom);
    } else {
        o_error("If neither of these are true we should revise this function\n");
    }
    const Mesh::Vertex & vertexA = mesh.VertexAt(ivA);
    const Mesh::Vertex & vertexB = mesh.VertexAt(ivB);
    const Mesh::Vertex & vertexC = mesh.VertexAt(ivC);
    if(glm::dot(vertexC.position - vertexA.position,vertexB.position - vertexA.position) <= 0){
        //AC
        if(Geo2D::DistanceSquared(vertexC.position - vertexA.position) >= diameterSquared)
            return true;
        else
            return false;
    }

    if(mesh.EdgeAt(adjacent).constrained){
        if(Geo2D::DistanceSquaredPointToLineSegment(vertexA.position, vertexB.position, vertexC.position) >= diameterSquared)
            return true;
        else
            return false;
    } else {
        //Check if neighbouring face(s) has/have enough clearance to allow the agent circle to navigate freely without hitting a constraint
        //First check if this face can be quickly discarded
        if(Geo2D::DistanceSquared(vertexC.position - vertexA.position) < diameterSquared || Geo2D::DistanceSquared(vertexC.position - vertexB.position) < diameterSquared)
            return false;
        else {
            
            Oryol::Set<uint32_t> checkedFaces;
            Oryol::Array<uint32_t> edgesToCheck;
            {
                const Mesh::HalfEdge & e = mesh.EdgeAt(adjacent);
                checkedFaces.Add(e.oppositeHalfEdge/4);
                edgesToCheck.Add(e.oppositeHalfEdge);
            };
            //o_error("Check me");
            //If not constrained there are two potential edges for each subsequent face we check.
            //Rather than using a stack to check for potentially constrained edges it may be possible to iterate through the neighboring halfedges directly
            while(!edgesToCheck.Empty()){
                uint32_t h = edgesToCheck.PopFront();
                
                const Mesh::HalfEdge & edge = mesh.EdgeAt(h);
                const Mesh::HalfEdge & next = mesh.EdgeAt(Mesh::Face::nextHalfEdge(h));
                const Mesh::HalfEdge & prev = mesh.EdgeAt(Mesh::Face::prevHalfEdge(h));
                const Mesh::Vertex & pivot = mesh.VertexAt(next.destinationVertex);
        
                if(!checkedFaces.Contains(next.oppositeHalfEdge/4)
                   && Geo2D::DistanceSquaredPointToLineSegment(pivot.position, mesh.VertexAt(edge.destinationVertex).position, vertexC.position) < diameterSquared){
                    if(next.constrained){
                        return false;
                    } else {
                        edgesToCheck.Add(next.oppositeHalfEdge);
                        checkedFaces.Add(next.oppositeHalfEdge/4);
                    }
                }
                
                
                if(!checkedFaces.Contains(prev.oppositeHalfEdge/4)
                   && Geo2D::DistanceSquaredPointToLineSegment(pivot.position, mesh.VertexAt(prev.destinationVertex).position, vertexC.position) < diameterSquared){
                    if(prev.constrained){
                        return false;
                    } else {
                        edgesToCheck.Add(prev.oppositeHalfEdge);
                        checkedFaces.Add(prev.oppositeHalfEdge/4);
                    }
                }
             
            }
            
        }
    }
    return true;
}

bool Path::FindPath(Mesh & mesh, const glm::dvec2 & start, const glm::dvec2 & end, const double radius, Oryol::Array<uint32_t> & pathFaces, Oryol::Array<uint32_t> & pathEdges){
    const double radiusSquared = radius * radius;
    const double diameterSquared = 4 * radiusSquared;
    uint32_t fromFace = -1, toFace = -1;
    
    //If either point falls on a vertex it's already considered constrained
    //Locate our first point on the mesh in question.
    {
        if(!mesh.GetBoundingBox().IsPointInside(start))
            return false;
        auto location = mesh.Locate(start);
        switch(location.type){
            case Mesh::LocateRef::None:
                return false;
            case Mesh::LocateRef::Vertex:
                fromFace = mesh.GetIncomingEdgeFor(location.object) / 4;
                break;
            case Mesh::LocateRef::Edge:
                fromFace = location.object / 4;
                break;
            case Mesh::LocateRef::Face:
                fromFace = location.object;
                break;
        }
    }
    //Locate our last point on the mesh in question.
    {
        if(!mesh.GetBoundingBox().IsPointInside(end))
            return false;
        auto location = mesh.Locate(end);
        switch(location.type){
            case Mesh::LocateRef::None:
                return false;
            case Mesh::LocateRef::Vertex:
                toFace = mesh.GetIncomingEdgeFor(location.object) / 4;
                break;
            case Mesh::LocateRef::Edge:
                toFace = location.object / 4;
                break;
            case Mesh::LocateRef::Face:
                toFace = location.object;
                break;
        }
    }
    o_assert(mesh.FaceAt(fromFace).isReal());
    o_assert(mesh.FaceAt(toFace).isReal());
    uint32_t currentFace;
    Oryol::Set<uint32_t> closed;
    Oryol::Set<uint32_t> openFaces;
    Oryol::Array<uint32_t> sortedOpenFaces;
    Oryol::Map<uint32_t, uint32_t> cameFrom;
    Oryol::Map<uint32_t, glm::dvec2> entryPositions;
    Oryol::Map<uint32_t,uint32_t> entryEdges;
    Oryol::Map<uint32_t, double> fScores; //F = G + H
    Oryol::Map<uint32_t, double> gScores; //Cost for face from start
    Oryol::Map<uint32_t, double> hScores; //Estimated Cost for face to end;
    
    hScores.AddUnique(fromFace,Geo2D::DistanceSquared(end-start));
    gScores.AddUnique(fromFace,0);
    fScores.AddUnique(fromFace,hScores[fromFace] + gScores[fromFace]);
    entryPositions.AddUnique(fromFace,start);
    entryEdges.AddUnique(fromFace,-1);
    openFaces.Add(fromFace);
    sortedOpenFaces.Add(fromFace);
    
    //The main criteria the A-Star search attempts to satisfy are;
    // * that we dont cross any constrained edges
    // * the circle representing the agent is able to pass from an edge through a face to the subsequent edge
    // * additional functionality would call into user code that is able to determine whether or not a constrained edge is passable or not.
    while(!openFaces.Empty()){
        
        currentFace = sortedOpenFaces.PopFront();
        openFaces.Erase(currentFace);
        
        if(!closed.Contains(currentFace)){
            if(currentFace == toFace)
                break;
            const Mesh::Face & face = mesh.FaceAt(currentFace);
            for(int i = 1; i < 4; i++){
                const Mesh::HalfEdge & e = face.edges[i-1];
                if(e.constrained) //TODO: Replace this condition with a callback
                    continue;
                uint32_t adjacentFace = e.oppositeHalfEdge/4;
                if(!closed.Contains(adjacentFace)){
                    o_assert(mesh.FaceAt(adjacentFace).isReal());
                    
                    //We have to validate that the face is passable
                    if(currentFace != fromFace && radius > 0 && !IsEdgeWalkable(mesh,entryEdges[currentFace],currentFace, e.oppositeHalfEdge, diameterSquared)){
                        continue;
                    }
                    const Mesh::Vertex & vA = mesh.VertexAt(e.destinationVertex);
                    const Mesh::Vertex & vB = mesh.VertexAt(mesh.EdgeAt(e.oppositeHalfEdge).destinationVertex);
                    
                    //TODO: Fix this metric because occasionally it can cause abnormally long paths
                    //A better way to calculate the entry position for very long edge lengths would be offset from the location of the pivot vertex
                    const auto entryPosition = (vA.position + vB.position) * 0.5;
                    
                    const double h = Geo2D::DistanceSquared(entryPosition - end);
                    const double g = gScores[currentFace] + Geo2D::DistanceSquared(entryPositions[currentFace] - entryPosition);
                    const double f = h + g;

                    if(!openFaces.Contains(adjacentFace)){
                        openFaces.Add(adjacentFace);
                        sortedOpenFaces.Add(adjacentFace);
                        entryPositions.AddUnique(adjacentFace,entryPosition);
                        entryEdges.AddUnique(adjacentFace,e.oppositeHalfEdge);
                        fScores.AddUnique(adjacentFace,f);
                        gScores.AddUnique(adjacentFace,g);
                        hScores.AddUnique(adjacentFace,h);
                        cameFrom.AddUnique(adjacentFace,currentFace);
                    } else if(fScores[adjacentFace] > f){
                        //We've found a better score for the adjacent face so rewrite those values.
                        entryPositions[adjacentFace] = entryPosition;
                        entryEdges[adjacentFace] = e.oppositeHalfEdge;
                        fScores[adjacentFace] = f;
                        gScores[adjacentFace] = g;
                        hScores[adjacentFace] = h;
                        cameFrom[adjacentFace] = currentFace;
                    }
                }
            }
            closed.Add(currentFace);
            //Only sort when the open edge list is modified
            std::sort(sortedOpenFaces.begin(), sortedOpenFaces.end(), [&fScores](uint32_t a, uint32_t b){ return fScores[a] < fScores[b]; });
            
        }
        currentFace = -1;
    }
    if(currentFace == (uint32_t)-1)
        return false;
    //Once the search is complete, reconstruct the sequence of faces in path.
    //path can then be fed into subsequent path refinement functions such as string pulling
    pathFaces.Add(currentFace);
    while(currentFace != fromFace){
        pathEdges.Insert(0,entryEdges[currentFace]);
        currentFace = cameFrom[currentFace];
        pathFaces.Insert(0,currentFace);
        
    }
    return true;
}



