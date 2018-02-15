//------------------------------------------------------------------------------
//  Triangle.cc
//------------------------------------------------------------------------------
#include "Pre.h"
#include "Core/Main.h"
#include "Gfx/Gfx.h"
#include "Input/Input.h"
#include "DebugBatch.h"
#include "glm/gtc/matrix_transform.hpp"
#include "Mesh.h"
#include "Geo2D.h"
#include "Path.h"
#include "Core/Time/Clock.h"
#include "IMUI/IMUI.h"
#include "imgui.h"

using namespace Oryol;
using namespace Delaunay;

class MeshDraw : private DebugBatch {
public:
    void Setup(const Oryol::GfxSetup & gfx){
        DebugBatch::Setup(gfx);
    }
    void Draw(const Mesh & mesh){
        for(uint32_t vIndex : mesh.ActiveVertexIndices()){
            if(vIndex != 0){
                const Delaunay::Mesh::Vertex & vertex = mesh.VertexAt(vIndex);
                this->DrawVertex(vertex.position);
                const uint32_t first = mesh.GetOutgoingEdgeFor(vIndex);
                uint32_t h = first;
                do {
                    const Delaunay::Mesh::HalfEdge & edge = mesh.EdgeAt(h);
                    const Delaunay::Mesh::HalfEdge & opposite = mesh.EdgeAt(edge.oppositeHalfEdge);
                    if((h < edge.oppositeHalfEdge) && edge.destinationVertex != 0 && opposite.destinationVertex != 0){
                        const Mesh::Vertex & destination = mesh.VertexAt(edge.destinationVertex);
                        this->DrawEdge(vertex,destination,edge.constrained);
                    }
                } while((h = mesh.GetNextOutgoingEdge(h)) != first);
            }
        }
        
    }
    void Submit(const glm::mat4 & mvp){
        DebugBatch::Draw(mvp);
    }
public:
    void DrawVertex(glm::vec2 position,float radius = 5, Color color = {1,1,1}) {
		this->Point(position.x, position.y, radius, color);
	}
    void DrawEdge(const Mesh::Vertex & origin, const Mesh::Vertex & destination, bool constrained) {
		if (constrained)
			this->Line(origin.position.x, origin.position.y, destination.position.x, destination.position.y, { 1,0,0,0.8f });
		else
			this->Line(origin.position.x, origin.position.y, destination.position.x, destination.position.y, { 1,1,1,0.8f });
	}
    void DrawFace(const Mesh & mesh,uint32_t fIndex,Color color){
        const Mesh::Face & face = mesh.FaceAt(fIndex);
        const Mesh::Vertex & vA = mesh.VertexAt(face.edges[0].destinationVertex);
        const Mesh::Vertex & vB = mesh.VertexAt(face.edges[1].destinationVertex);
        const Mesh::Vertex & vC = mesh.VertexAt(face.edges[2].destinationVertex);
        this->Triangle(vA.position.x, vA.position.y, vB.position.x, vB.position.y, vC.position.x, vC.position.y, color);
    }
};

class DelaunayApp : public App {
public:
    AppState::Code OnRunning();
    AppState::Code OnInit();
    AppState::Code OnCleanup();
	glm::mat4 projectionMatrix;
    int mouseMode = 0;
    bool showMenu = true;
	MeshDraw debug;
    TimePoint lastTimePoint;
	Delaunay::Mesh mesh;
    Delaunay::Mesh::LocateRef location;
    glm::vec2 locatedPosition;
    Array<uint32_t> pathFaces,pathEdges;

};
OryolMain(DelaunayApp);

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnInit() {
    // setup rendering system
    Gfx::Setup(GfxSetup::Window(800, 600, "Oryol Delaunay Sample"));
    Input::Setup();
    IMUI::Setup();
    
    debug.Setup(Gfx::GfxSetup());

	mesh.Setup(500, 500);
    //mesh.InsertConstraintSegment({0,100}, {175,100});
    //mesh.InsertConstraintSegment({400,100}, {225,100});
    //mesh.InsertConstraintSegment({100,300}, {300,300});
    //mesh.InsertConstraintSegment({0,200}, {300,200});
    //mesh.InsertConstraintSegment({300,200}, {500,400});
	mesh.InsertConstraintSegment({350,400}, {400,350});
	mesh.InsertConstraintSegment({ 150,400 }, { 100,350 });
	mesh.InsertConstraintSegment({150,400}, {450,400});
	mesh.InsertConstraintSegment({ 100,350 }, { 300,200 });

    Path::FindPath(mesh, {85,65}, {450,450}, 25, pathFaces, pathEdges);
    
	projectionMatrix = glm::ortho<float>(0, 800, 600, 0, -10, 10);
    projectionMatrix *= glm::translate(glm::mat4(1.0f), {50,50,0});
    lastTimePoint = Clock::Now();
    return App::OnInit();
}
static const char * names[5] = { "None","Vertex","Edge","Face" };
const char * mouseModeStrings[] = {"Locate Position","Add Constraint"};
static const float menuWidth = 200;
//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnRunning() {
    
    Gfx::BeginPass();
    
	IMUI::NewFrame(Clock::LapTime(this->lastTimePoint));
    ImGui::SetNextWindowPos(ImVec2((float)800 - menuWidth - 10, 10));
    ImGui::SetNextWindowSize(ImVec2((float)menuWidth, (float)600 - 20));
    ImGui::Begin("Testbed Controls", &showMenu, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
    ImGui::PushAllowKeyboardFocus(false); // Disable TAB
    
    ImGui::Text("Mouse Mode");
    if(ImGui::Combo("##MouseMode",&mouseMode,mouseModeStrings,IM_ARRAYSIZE(mouseModeStrings))){
        locatedPosition = {};
        location = {};
    }
    if (mouseMode == 0){
        if(!ImGui::GetIO().WantCaptureMouse && Input::MouseButtonDown(MouseButton::Left)) {
            locatedPosition = Input::MousePosition();
            locatedPosition.x = std::round(locatedPosition.x);
            locatedPosition.y = std::round(locatedPosition.y);
            location = mesh.Locate({locatedPosition.x - 50,locatedPosition.y - 50});
        }
        ImGui::Separator();
        ImGui::Text("Location: (%.1f, %.1f)",locatedPosition.x,locatedPosition.y);
        ImGui::Text("Current Object: %s", names[location.type]);
        ImGui::Text("Object Index: %u",location.object);
        ImGui::Separator();
        if(location.type == Mesh::LocateRef::Vertex){
            const auto & vertex = mesh.VertexAt(location.object);
            ImGui::Text("End Point Count: %zu",vertex.endPointCount);
            ImGui::Text("Constrained Edge Count: %zu",vertex.constraintCount);
            debug.DrawVertex(vertex.position, 25, {0,1,0,0.3});
        } else if(location.type == Mesh::LocateRef::Face){
            const auto & face = mesh.FaceAt(location.object);
            ImGui::Text("Vertex A: %u",face.edges[0].destinationVertex);
            ImGui::Text("Vertex B: %u",face.edges[1].destinationVertex);
            ImGui::Text("Vertex C: %u",face.edges[2].destinationVertex);
            debug.DrawFace(mesh, location.object, {0,1,0,0.3});
        }
    }

    ImGui::PopAllowKeyboardFocus();
    ImGui::End();
    
    ImGui::Render();
    for(auto f : pathFaces){
        debug.DrawFace(mesh, f, {0.7f,0,0.3f,0.4f});
    }
    debug.Draw(mesh);
    debug.Submit(projectionMatrix);

    Gfx::EndPass();
    Gfx::CommitFrame();
    
    // continue running or quit?
    return Gfx::QuitRequested() ? AppState::Cleanup : AppState::Running;
}

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnCleanup() {
	Input::Discard();
    Gfx::Discard();
	
    return App::OnCleanup();
}
