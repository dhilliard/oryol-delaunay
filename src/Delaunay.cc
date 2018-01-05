//------------------------------------------------------------------------------
//  Triangle.cc
//------------------------------------------------------------------------------
#include "Pre.h"
#include "Core/Main.h"
#include "Gfx/Gfx.h"
#include "DebugDraw.h"
#include "glm/gtc/matrix_transform.hpp"

using namespace Oryol;

class DelaunayApp : public App {
public:
    AppState::Code OnRunning();
    AppState::Code OnInit();
    AppState::Code OnCleanup();
	glm::mat4 projectionMatrix;
	DebugDraw debug;
};
OryolMain(DelaunayApp);

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnInit() {
    // setup rendering system
    Gfx::Setup(GfxSetup::Window(400, 400, "Oryol Delaunay Sample"));
    
	debug.Setup(GfxSetup());
	projectionMatrix = glm::ortho<float>(0, 400, 400, 0, -10, 10);
    return App::OnInit();
}

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnRunning() {
    
    Gfx::BeginPass();
	debug.Line(0, 0, 100, 100, { 1,1,1 });
	debug.Point(200, 200, 10, { 1,0,1 });
	debug.Draw(projectionMatrix);

    Gfx::EndPass();
    Gfx::CommitFrame();
    
    // continue running or quit?
    return Gfx::QuitRequested() ? AppState::Cleanup : AppState::Running;
}

//------------------------------------------------------------------------------
AppState::Code
DelaunayApp::OnCleanup() {
    Gfx::Discard();
    return App::OnCleanup();
}
