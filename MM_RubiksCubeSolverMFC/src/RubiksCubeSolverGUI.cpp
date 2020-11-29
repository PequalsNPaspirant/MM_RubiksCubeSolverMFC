//=======================================================================================================//
//   Copyright (c) 2018 Maruti Mhetre                                                                    //
//   All rights reserved.                                                                                //
//=======================================================================================================//
//   Redistribution and use of this software in source and binary forms, with or without modification,   //
//   are permitted for personal, educational or non-commercial purposes provided that the following      //
//   conditions are met:                                                                                 //
//   1. Redistributions of source code must retain the above copyright notice, this list of conditions   //
//      and the following disclaimer.                                                                    //
//   2. Redistributions in binary form must reproduce the above copyright notice, this list of           //
//      conditions and the following disclaimer in the documentation and/or other materials provided     //
//      with the distribution.                                                                           //
//   3. Neither the name of the copyright holder nor the names of its contributors may be used to        //
//      endorse or promote products derived from this software without specific prior written            //
//      permission.                                                                                      //
//=======================================================================================================//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR      //
//   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND    //
//   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR          //
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL   //
//   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   //
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER  //
//   IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT   //
//   OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     //
//=======================================================================================================//

// RubiksCube.cpp : Defines the entry point for the application.
//
#include "stdafx.h"

//#include "Resource.h"
#include <windows.h>
#include <windowsx.h>
#include <cassert>
//#include <thread>
//#include <chrono>
#include <experimental/filesystem> // or #include <filesystem>
using namespace std;

#include "RubiksCubeSolverGUI.h"
#include "RubiksCubeSolverUtils.h"
#include "RubiksCubeSolverTest.h"
#include "MM_RubiksCubeSolverMFCDlg.h"

namespace mm {

	RubiksCubeSolverGUI::RubiksCubeSolverGUI()
		: //WND_WIDTH(800),
		//WND_HEIGHT(800),
		//messageWndHeight(110),
		//currentModelName_("RubiksCubeModel_v7"),
		//currentModelName_("RubiksCubeModel_v8"),
		//currentModelName_("RubiksCubeModel_v9"),
		currentModelName_("RubiksCubeModel_v10"),
		//scene_(*this, "RubiksCubeModel_v1", 3),
		//scene_(*this, "RubiksCubeModel_v2", 3),
		//scene_(*this, "RubiksCubeModel_v3", 3),
		//scene_(*this, "RubiksCubeModel_v3", 2),
		//scene_(*this, "RubiksCubeModel_v4", 2),
		//scene_(*this, "RubiksCubeModel_v4", 3),
		//scene_(*this, "RubiksCubeModel_v4", 4),
		//scene_(*this, "RubiksCubeModel_v4", 10),
		//scene_(*this, "RubiksCubeModel_v4", 100),
		//scene_(*this, "RubiksCubeModel_v5", 4),
		//scene_(*this, "RubiksCubeModel_v5", 5),
		//scene_(*this, "RubiksCubeModel_v5", 6),
		//scene_(*this, "RubiksCubeModel_v5", 7),
		//scene_(*this, "RubiksCubeModel_v5", 8),
		//scene_(*this, "RubiksCubeModel_v6", 3),
		//scene_(*this, "RubiksCubeModel_v6", 5),
		//scene_(*this, "RubiksCubeModel_v7", 3),
		scene_(*this, currentModelName_, 3),
		framesPerRotation_(20), //moderate
		sleepTimeMilliSec_(20), //moderate
		rubiksCubeSize_(3),
		tester_(*this),
		solutionDirectory_{ "C:/RubiksCubeSolutions" }
		//selMenuAnimationSpeed(ID_ANIMATIONSPEED_MODERATE),
		//selMenuRubiksCubeSize(ID_RUBIK_3X3X3)
	{
		//Create a directory to store Rubik's Cube scramble and solutions
		namespace fs = std::experimental::filesystem;
		if (!fs::is_directory(solutionDirectory_) || !fs::exists(solutionDirectory_)) {
			fs::create_directory(solutionDirectory_);
		}
	}

	RubiksCubeSolverGUI::~RubiksCubeSolverGUI()
	{
		renderingThread_.join();
	}

	void RubiksCubeSolverGUI::initialize(HWND hWnd)
	{
		g_hWnd = hWnd;
		renderingThread_ = std::thread(&mm::RubiksCubeSolverGUI::render, this);
	}

	bool RubiksCubeSolverGUI::activateRubiksCube(int size)
	{
		rubikCubeSize_ = size;
		cubeType_ = cubeType::rubiksCube;
		firstGenCommand_ = firstGenerationCommands::eSetCubeType;
		return activateRenderingThread();
	}

	bool RubiksCubeSolverGUI::activateMirrorCube()
	{
		cubeType_ = cubeType::mirrorCube;
		firstGenCommand_ = firstGenerationCommands::eSetCubeType;
		return activateRenderingThread();
	}

	void RubiksCubeSolverGUI::exitUI()
	{
		interruptAnimation_ = true;
		keepRunning_ = false;
	}

	void RubiksCubeSolverGUI::Scramble(const string& scramblingAlgo, bool animateIn)
	{
		scramblingAlgo_ = scramblingAlgo;
		animate_ = animateIn;
		firstGenCommand_ = firstGenerationCommands::eScramble;
		activateRenderingThread();
	}
	void RubiksCubeSolverGUI::Solve(bool animateIn)
	{
		animate_ = animateIn;
		firstGenCommand_ = firstGenerationCommands::eSolve;
		activateRenderingThread();
	}
	void RubiksCubeSolverGUI::runTests(bool animateIn)
	{
		animate_ = animateIn;
		firstGenCommand_ = firstGenerationCommands::eRunTests;
		activateRenderingThread();
	}
	bool RubiksCubeSolverGUI::setRubiksCubeSize(unsigned int size)
	{
		if (0 < size && size < 101)
		{
			firstGenCommand_ = firstGenerationCommands::eResizeRubiksCube;
			rubikCubeSize_ = size;
			return activateRenderingThread();
		}
		else
		{
			RubiksCubeSolverUtils::CreateOkDialog("The Rubik's cube size must be in the range [1, 100]");
			return false;
		}
	}

	//second generation commands
	void RubiksCubeSolverGUI::resetRubiksCube()
	{
		if(firstGenCommand_ != firstGenerationCommands::eNoCommand) //There is first generation command active
			interruptAnimation_ = true;

		secondGenCommand_ = secondGenerationCommands::eResetRubiksCube;
		activateRenderingThread(true);
	}

	//Third generation commands
	void RubiksCubeSolverGUI::setAnimationSpeed(unsigned int speed)
	{
		//firstGenCommand_ = firstGenerationCommands::eNoCommand;
		//secondGenCommand_ = secondGenerationCommands::eSetAnimationSpeed;
		animationSpeed_ = speed;
		framesPerRotation_ = (106 - animationSpeed_) / 2;
		sleepTimeMilliSec_ = (106 - animationSpeed_) / 2;
		//activateRenderingThread();
	}
	void RubiksCubeSolverGUI::fitToScreen()
	{
		//firstGenCommand_ = firstGenerationCommands::eNoCommand;
		//secondGenCommand_ = secondGenerationCommands::eFitToScreen;
		//interruptAnimation_ = true;
		scene_.fitToScreen();
		thirdGenCommand_ = thirdGenerationCommands::eFitToScreen;
		activateRenderingThread();
	}

#if 0
	//Implementation using atomic variable and spin lock
	void RubiksCubeSolverGUI::render()
	{
		if (!graphicsAreaCreated_)
		{
			createGraphicsArea();
			graphicsAreaCreated_ = true;
		}

		while (keepRunning_)
		{
			//bool expected = true;
			//bool desired = false;
			//if (renderNow_.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
			if(renderNow_.load(std::memory_order_acquire))
			{
				try 
				{
					commandHandlerThirdGen(); //Always run these commands
					commandHandlerFirstGen();
				}
				catch (bool flag)
				{
					//The previous command is broken/interrupted, reset the rubik cube.
					//If we need to perform different actions on the type of interrupt, we can include that information in the exception object
					interruptAnimation_ = false; //reset the flag
					firstGenCommand_ = firstGenerationCommands::eNoCommand;
				}
				commandHandlerSecondGen();
				renderNow_.store(false, std::memory_order_release);

				//std::thread t1(&RubiksCubeSolverGUI::Scramble, this, animate);
				//t1.detach();
			}
		}
	}
#endif

	void RubiksCubeSolverGUI::waitOnConditionVariable()
	{
		//classic double (null) check case to avoid unnecessary locking and unlocking mutex (can happen if we dont have this check and ready_ is true)
		//Another benefit of this check is that, when ready_ is true, it is not reset to false at the end of this function
		if (!ready_) 
		{
			std::unique_lock<std::mutex> lock(mtx_);
			while (!ready_)
				cv_.wait(lock);
			//ready_ = false; //reset flag so that next time we come here, the thread will get stuck in wait()
		}
	}

	void RubiksCubeSolverGUI::setReadySynchronously(bool ready)
	{
		std::unique_lock<std::mutex> lock(mtx_);
		ready_ = ready;
	}

	void RubiksCubeSolverGUI::render()
	{
		if (!graphicsAreaCreated_)
		{
			createGraphicsArea();
			graphicsAreaCreated_ = true;
		}

		while (keepRunning_)
		{
			waitOnConditionVariable(); //blocking call...waits until condition variable is notified

			try
			{
				firstGenCommandInProgress_ = firstGenCommand_;
				commandHandlerThirdGen(); //Always run these commands
				commandHandlerFirstGen();
			}
			catch (bool)
			{
				//The previous command is broken/interrupted, reset the rubik cube.
				//If we need to perform different actions on the type of interrupt, we can include that information in the exception object
				interruptAnimation_ = false; //reset the flag
				firstGenCommand_ = firstGenerationCommands::eNoCommand;
			}

			commandHandlerSecondGen();

			firstGenCommandInProgress_ = firstGenerationCommands::eNoCommand;
			setReadySynchronously(false);
		}
	}

	void RubiksCubeSolverGUI::createGraphicsArea()
	{
		GetClientRect(g_hWnd, &g_rWnd);
		g_hDC = GetDC(g_hWnd);

		if (!setupPixelFormat(g_hDC))
			PostQuitMessage(-1);

		g_hRC = wglCreateContext(g_hDC);
		wglMakeCurrent(g_hDC, g_hRC);

		//scene_.initOpenGl(g_rWnd.right, g_rWnd.bottom - messageWndHeight);
		scene_.initOpenGl(g_rWnd.right, g_rWnd.bottom);
		scene_.initScene();

		redrawWindow();
	}

	bool RubiksCubeSolverGUI::setupPixelFormat(HDC hdc)
	{
		PIXELFORMATDESCRIPTOR pfd = { 0 };
		int pixelformat;

		pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR); // Set the size of the structure
		pfd.nVersion = 1; // Always set this to 1

		pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER; // Pass in the appropriate OpenGL flags
		pfd.dwLayerMask = PFD_MAIN_PLANE;					// We want the standard mask (this is ignored anyway)
		pfd.iPixelType = PFD_TYPE_RGBA;						// We want RGB and Alpha pixel type
		pfd.cColorBits = SCREEN_DEPTH;						// Here we use our #define for the color bits
		pfd.cDepthBits = SCREEN_DEPTH;						// Depthbits is ignored for RGBA, but we do it anyway
		pfd.cAccumBits = 0;									// No special bitplanes needed
		pfd.cStencilBits = 0;								// We desire no stencil bits

															// This gets us a pixel format that best matches the one passed in from the device
		if ((pixelformat = ChoosePixelFormat(hdc, &pfd)) == false)
			return false;

		// This sets the pixel format that we extracted from above
		if (SetPixelFormat(hdc, pixelformat, &pfd) == false)
			return false;

		return true;
	}

	void RubiksCubeSolverGUI::redrawWindow()
	{
		scene_.renderScene();
		bool result = SwapBuffers(g_hDC);
	}

	void RubiksCubeSolverGUI::applyAlgorithm(const string& algo, bool animate)
	{
		string invalidStep;
		scene_.scramble(algo, animate, invalidStep);
	}

	unique_ptr<RubiksCubeModel> RubiksCubeSolverGUI::replaceModelBy(const string& modelName, int size, bool animate)
	{
		unique_ptr<RubiksCubeModel> originalModel = scene_.replaceModelBy(modelName, size);
		if(animate)
			displayUpdatedStats();
		return std::move(originalModel);
	}

	unique_ptr<RubiksCubeModel> RubiksCubeSolverGUI::replaceModelBy(unique_ptr<RubiksCubeModel>&& newModel, bool animate)
	{
		unique_ptr<RubiksCubeModel> originalModel = scene_.replaceModelBy(std::move(newModel));
		return std::move(originalModel);
	}

	bool RubiksCubeSolverGUI::isSolved()
	{
		return scene_.isSolved();
	}

	string getCommaSeparatedTimeDuration(unsigned long long duration)
	{
		string durationStr = "000,000.000,000,000";
		int pos = static_cast<int>(durationStr.length());
		for (; pos > 0 && duration > 0; pos)
		{
			if (durationStr[--pos] == '0')
			{
				durationStr[pos] = '0' + duration % 10;
				duration /= 10;
			}
		}
		if (pos > 6)
			pos = 6;
		durationStr = durationStr.substr(pos);
		durationStr += " sec";

		return durationStr;
	}

	void RubiksCubeSolverGUI::displayUpdatedStats()
	{
		//Fill up message area with different color
		//RECT messageWndRect;
		//messageWndRect.left = 0;
		//messageWndRect.right = WND_WIDTH;
		//messageWndRect.top = 0;
		//messageWndRect.bottom = messageWndHeight;
		//HBRUSH brush = CreateSolidBrush(RGB(255, 255, 255));
		//FillRect(g_hDCMessage, &messageWndRect, brush);

		unsigned int size;
		unsigned int scramblingSteps;
		string scramblingAlgo;
		unsigned int solutionSteps;
		string solution;
		unsigned long long duration;
		string statusStr;
		scene_.getUpdatedStats(size, scramblingSteps, scramblingAlgo, solutionSteps, solution, duration, statusStr);

		string sizeStr{ to_string(scene_.getRubiksCubeSize()) };
		string rubikCubeSize("Rubik's Cube Size: " + sizeStr + "x" + sizeStr + "x" + sizeStr);
		//string scramblingStepsStr("Scrambling Steps: " + (scramblingSteps > 0 ? to_string(scramblingSteps) : ""));
		string scrambleMsg("Scrambling Algorithm (" + to_string(scramblingSteps) + " steps): " + scramblingAlgo);
		//string solutionStepsStr("Solution Steps: " + (solutionSteps > 0 ? to_string(solutionSteps) : ""));
		string solutionMsg("Solution Algorithm (" + to_string(solutionSteps) + " steps): " + solution);
		string durationStr("Time required to solve: " + getCommaSeparatedTimeDuration(duration));

		//Hourglass animations
		static const vector<string> hourglass{
			""
		};

		//static const vector<string> hourglass{
		//	"~",
		//	"~~",
		//	"~~~",
		//	"~~~~"
		//};

		//static const vector<string> hourglass{
		//	".",
		//	"..",
		//	"...",
		//	"...."
		//};

		//static const vector<string> hourglass{ "|", "/", "-", "\\" };

		//static const vector<string> hourglass{
		//	">-----",
		//	"->----",
		//	"-->---",
		//	"--->--",
		//	"---->-",
		//	"----->",
		//	"-----<",
		//	"----<-",
		//	"---<--",
		//	"--<---",
		//	"-<----",
		//	"<-----",
		//};

		//static const vector<string> hourglass{
		//	">    ---",
		//	"->    --",
		//	"-->    -",
		//	"--->    ",
		//	" --->   ",
		//	"  --->  ",
		//	"   ---> ",
		//	"    --->",
		//};

		//vector<string> hourglass{
		//	"|****       |",
		//	"|*** *      |",
		//	"|** * *     |",
		//	"|* * * *    |",
		//	"| * * * *   |",
		//	"|  * * * *  |",
		//	"|   * * * * |",
		//	"|    * * * *|",
		//	"|     * * **|",
		//	"|      * ***|",
		//	"|       ****|",
		//};
		//hourglass.insert(hourglass.end(), hourglass.rbegin(), hourglass.rend());

		//vector<string> hourglass{
		//	"|====       |",
		//	"|=== =      |",
		//	"|== = =     |",
		//	"|= = = =    |",
		//	"| = = = =   |",
		//	"|  = = = =  |",
		//	"|   = = = = |",
		//	"|    = = = =|",
		//	"|     = = ==|",
		//	"|      = ===|",
		//	"|       ====|",
		//};
		//hourglass.insert(hourglass.end(), hourglass.rbegin(), hourglass.rend());

		//vector<string> hourglass{
		//	"|====          |",
		//	"|=== =         |",
		//	"|== =  =       |",
		//	"|= =  =   =    |",
		//	"| =  =   =    =|",
		//	"|   =   =    ==|",
		//	"|      =    ===|",
		//	"|          ====|"
		//};
		//hourglass.insert(hourglass.end(), hourglass.rbegin(), hourglass.rend());

		//vector<string> hourglass{
		//	"|====             |",
		//	" |====            |",
		//	"   |====          |",
		//	"   |=== =         |",
		//	"   |== =  =       |",
		//	"   |= =  =   =    |",
		//	"   | =  =   =    =|",
		//	"   |   =   =    ==|",
		//	"   |      =    ===|",
		//	"   |          ====|",
		//	"   |            ====|",
		//	"   |             ====|"
		//};
		//hourglass.insert(hourglass.end(), hourglass.rbegin(), hourglass.rend());

		//vector<string> hourglass{
		//	"{====             }",
		//	" {====            }",
		//	"   {====          }",
		//	"   {=== =         }",
		//	"   {== =  =       }",
		//	"   {= =  =   =    }",
		//	"   { =  =   =    =}",
		//	"   {   =   =    ==}",
		//	"   {      =    ===}",
		//	"   {          ====}",
		//	"   {            ====}",
		//	"   {             ====}"
		//};
		//hourglass.insert(hourglass.end(), hourglass.rbegin(), hourglass.rend());

		//vector<string> hourglass{
		//	"|XX|====              | X X | ",
		//	" |XX|====             | X X | ",
		//	" | XX|=== =           | X X | ",
		//	" | X X|== =  =        | X X | ",
		//	" | X X |= =  =   =    | X X | ",
		//	" | X X | =  =   =    =| X X | ",
		//	" | X X |   =   =     ==|X X | ",
		//	" | X X |      =      ===|XX | ",
		//	" | X X |             ====|XX| ",
		//};
		//int count = hourglass.size();
		//for (int i = 0; i < count; ++i)
		//	hourglass.push_back(string{ hourglass[i].rbegin(), hourglass[i].rend() });

		static int hourglassIndex = -1;
		hourglassIndex = (++hourglassIndex) % hourglass.size();
		string hourglassStr{};
		//if (isScrambling_ || isSolving_)
		{
			hourglassStr = hourglass[hourglassIndex];
		}

		CMMRubiksCubeSolverMFCDlg::getMainDailog().displayMessage({
			rubikCubeSize,
			//scramblingStepsStr,
			scrambleMsg,
			//solutionStepsStr,
			solutionMsg,
			durationStr,
			statusStr,
			hourglassStr
			});
	}

	string RubiksCubeSolverGUI::generateScramblingAlgo(int length)
	{
		return scene_.generateScramblingAlgo(length);
	}

	void RubiksCubeSolverGUI::pauseAnimation(bool pause)
	{
		//scene_.pauseAnimation(pause);
		setReadySynchronously(!pause);
		if(pause == false)
			cv_.notify_all();
	}

	// ============================== Message Handling ================================= //

	void RubiksCubeSolverGUI::OnRubiksCubeChanged(HWND hWnd)
	{
		// check for solution
		//if (scene_.g_cCube.IsSolved())
		{
			//TCHAR solvedMsg[MAX_LOADSTRING];
			//LoadString(g_hInstance, IDS_SOLVED, solvedMsg, MAX_LOADSTRING);
			//MessageBox(g_hWnd, solvedMsg, g_szTitle, MB_OK);
		}
	}

	void RubiksCubeSolverGUI::OnPaint(HWND hWnd)
	{
		redrawWindow();
	}

	BOOL RubiksCubeSolverGUI::OnEraseBackground(HWND hwnd, HDC hdc)
	{
		//redrawWindow();
		return TRUE;
	}

	//  Process WM_LBUTTONDOWN message for window/dialog: 
	void RubiksCubeSolverGUI::OnLButtonDown(HWND hWnd, BOOL fDoubleClick, int x, int y, UINT keyFlags)
	{
		// set up tracking for when mouse leaves window
		TRACKMOUSEEVENT tme;
		tme.cbSize = sizeof(TRACKMOUSEEVENT);
		tme.dwFlags = TME_LEAVE;
		tme.hwndTrack = hWnd;
		TrackMouseEvent(&tme);

		g_bMouseDown = true;

		if ((g_nHitCount = scene_.getSelectedObjects(x, y, g_rWnd.right, g_rWnd.bottom)) > 0)
			g_vMouseDown = scene_.mapCoordinates(x, y);

		g_nPrevX = x;
		g_nPrevY = y;

		//redrawWindow();
	}

	//  Process WM_LBUTTONUP message for window/dialog: 
	void RubiksCubeSolverGUI::OnLButtonUp(HWND hWnd, int x, int y, UINT keyFlags)
	{
		g_bMouseDown = false;

		/*
		if (g_bRotating)
		{
		if (abs(g_nRotationAngle) >= 45)
		{
		int turns = 1;

		if (g_bFlipRotation)
		g_nRotationAngle *= -1;

		if (g_nRotationAngle < 0)
		turns = 3;
		else if (g_nRotationAngle > 0)
		turns = 1;

		if (g_vRotationAxis.x)
		g_cCube.Tilt(g_nRotatingSection, turns);
		else if (g_vRotationAxis.y)
		g_cCube.Rotate(g_nRotatingSection, turns);
		else if (g_vRotationAxis.z)
		g_cCube.Turn(g_nRotatingSection, turns);
		}

		g_bRotating = false;
		g_bFlipRotation = false;

		PostMessage(g_hWnd, RC_CHANGED, 0, 0);
		}
		*/

		//redrawWindow();
	}

	//  Process WM_DESTROY message for window/dialog: 
	void RubiksCubeSolverGUI::OnDestroy(HWND hWnd)
	{
		//deInit();
	}

	//  Process WM_MOUSEMOVE message for window/dialog: 
	void RubiksCubeSolverGUI::OnRotate(int rotate, int tilt)
	{
		thirdGenCommand_ = thirdGenerationCommands::eRotate;
		//if (!g_bMouseDown)
		//{
		//	/*
		//	if ((g_nHitCount = GetSelectedObjects(x, y, g_rWnd.right, g_rWnd.bottom)) > 0)
		//	SetCursor(g_hHand);
		//	else
		//	SetCursor(g_hArrow);
		//	*/
		//}

		// moving camera
		//else if (g_nHitCount == 0 && g_bMouseDown)
		//{
			//if (x < g_nPrevX)
			//	scene_.g_cCamera.Rotate(-5);
			//else if (x > g_nPrevX)
			//	scene_.g_cCamera.Rotate(5);

			//if (y < g_nPrevY)
			//	scene_.g_cCamera.Tilt(5);
			//else if (y > g_nPrevY)
			//	scene_.g_cCamera.Tilt(-5);

		if(rotate != 0)
			scene_.g_cCamera.Rotate(rotate);
		if(tilt != 0)
			scene_.g_cCamera.Tilt(tilt);
		//firstGenCommand_ = firstGenerationCommands::eNoCommand;
		activateRenderingThread(true);
			//redrawWindow();
		//}
		/*
		// rotating section
		else if (g_nHitCount > 0 && g_bMouseDown && !g_bRotating)
		{
		int deltaX = abs(x - g_nPrevX);
		int deltaY = abs(y - g_nPrevY);
		int i, j , k = 0;
		Face face;

		GetCubeSelection(&i, &j, &k, &face, g_nHitCount);

		if (deltaX > 3 || deltaY > 3)
		{
		if (deltaX >= deltaY)
		{
		if (face == Top || face == Bottom)
		{
		float phi = g_cCamera.GetPhi();

		if (phi < 45 || (phi >= 135 && phi <= 225))
		{
		g_vRotationAxis = CVector3(1, 0, 0);
		g_nRotatingSection = i;
		g_vMouseDown.x = 0;

		g_bFlipRotation = !g_bFlipRotation;
		}

		else
		{
		g_vRotationAxis = CVector3(0, 0, 1);
		g_nRotatingSection = k;
		g_vMouseDown.z = 0;
		}
		}
		else if (face == Left)
		{
		g_vRotationAxis = CVector3(0, 1, 0);
		g_nRotatingSection = j;
		g_vMouseDown.y = 0;
		}
		else if (face == Right)
		{
		g_vRotationAxis = CVector3(0, 1, 0);
		g_nRotatingSection = j;
		g_vMouseDown.y = 0;
		}
		else if (face == Front)
		{
		g_vRotationAxis = CVector3(0, 1, 0);
		g_nRotatingSection = j;
		g_vMouseDown.y = 0;
		}
		else if (face == Back)
		{
		g_vRotationAxis = CVector3(0, 1, 0);
		g_nRotatingSection = j;
		g_vMouseDown.y = 0;
		}

		if (x - g_nPrevX < 0)
		g_bFlipRotation = !g_bFlipRotation;
		}

		else
		{
		if (face == Top || face == Bottom)
		{
		float phi = g_cCamera.GetPhi();

		if (phi < 45 || (phi >= 135 && phi <= 225))
		{

		g_vRotationAxis = CVector3(0, 0, 1);
		g_nRotatingSection = k;
		g_vMouseDown.z = 0;
		}

		else
		{
		g_vRotationAxis = CVector3(1, 0, 0);
		g_nRotatingSection = i;
		g_vMouseDown.x = 0;
		}

		g_bFlipRotation = !g_bFlipRotation;
		}
		else if (face == Left)
		{
		g_vRotationAxis = CVector3(0, 0, 1);
		g_nRotatingSection = k;
		g_vMouseDown.z = 0;
		}
		else if (face == Right)
		{
		g_vRotationAxis = CVector3(0, 0, 1);
		g_nRotatingSection = k;
		g_vMouseDown.z = 0;
		g_bFlipRotation = !g_bFlipRotation;
		}
		else if (face == Front)
		{
		g_vRotationAxis = CVector3(1, 0, 0);
		g_nRotatingSection = i;
		g_vMouseDown.x = 0;
		}
		else if (face == Back)
		{
		g_vRotationAxis = CVector3(1, 0, 0);
		g_nRotatingSection = i;
		g_vMouseDown.x = 0;
		g_bFlipRotation = !g_bFlipRotation;
		}

		if (y - g_nPrevY < 0)
		g_bFlipRotation = !g_bFlipRotation;
		}

		if (g_cCamera.IsFlipped())
		g_bFlipRotation = !g_bFlipRotation;

		if (g_cCamera.GetPhi() >= 90 && g_cCamera.GetPhi() <= 270 && (face == Top || face == Bottom))
		g_bFlipRotation = !g_bFlipRotation;

		g_bRotating = true;
		}
		}

		if (g_bRotating)
		{
		CVector3 pos = MapCoordinates(x, y);

		if (g_vRotationAxis.x)
		pos.x = 0;
		else if (g_vRotationAxis.y)
		pos.y = 0;
		else if (g_vRotationAxis.z)
		pos.z = 0;

		double angle = g_vMouseDown.GetAngle(pos);

		g_nRotationAngle = (int)(angle * 180 / PI);

		if (g_nRotationAngle > 90)
		g_nRotationAngle = 90;
		else if (g_nRotationAngle < -90)
		g_nRotationAngle = -90;
		}
		*/
		//g_nPrevX = x;
		//g_nPrevY = y;

		//redrawWindow();
	}

	void RubiksCubeSolverGUI::OnPan(int horizontal, int vertical)
	{
		return;

		thirdGenCommand_ = thirdGenerationCommands::ePan;
		activateRenderingThread(true);

		CVector3 pos = scene_.g_cCamera.GetLookAt();

		//Pan vertically i.e. move pos in up direction by amount = vertical
		CVector3 up = scene_.g_cCamera.GetUp();
		pos.x += (up.x * vertical);
		pos.y += (up.y * vertical);
		pos.z += (up.z * vertical);

		//Pan horizontally i.e. move pos in horizontal direction by amount = horizontal
		CVector3 screenNormal = scene_.g_cCamera.GetScreenNormal().Unit();
		CVector3 panDir = up ^ screenNormal;
		pos.x += (panDir.x * horizontal);
		pos.y += (panDir.y * horizontal);
		pos.z += (panDir.z * horizontal);
		scene_.g_cCamera.SetLookAt(pos);

	}

	//  Process WM_MOUSEWHEEL message for window/dialog: 
	void RubiksCubeSolverGUI::OnZoom(float distance)
	{
		scene_.g_cCamera.Move(distance);
		thirdGenCommand_ = thirdGenerationCommands::eZoom;
		activateRenderingThread(true);
	}

	//  Process WM_SIZE message for window/dialog: 
	void RubiksCubeSolverGUI::OnSize(int cx, int cy)
	{
		g_rWnd.left = 0;
		g_rWnd.right = cx;
		g_rWnd.top = 0;
		g_rWnd.bottom = cy;
		thirdGenCommand_ = thirdGenerationCommands::eResizeWindow;
		activateRenderingThread(true);
	}

	void RubiksCubeSolverGUI::OnSizeImpl()
	{
		//if (!graphicsAreaCreated_)
		//	return;

		//////////GetClientRect(g_hWnd, &g_rWnd);
		//////////g_hDC = GetDC(g_hWnd);

		/////////if (!setupPixelFormat(g_hDC))
		////////	PostQuitMessage(-1);

		//g_hRC = wglCreateContext(g_hDC);
		//wglMakeCurrent(g_hDC, g_hRC);

		//scene_.initOpenGl(g_rWnd.right, g_rWnd.bottom);
		scene_.sizeOpenGlScreen(g_rWnd.right, g_rWnd.bottom);

		/////////scene_.initScene();

		//redrawWindow();

		///////activateRenderingThread(true);
	}

	void RubiksCubeSolverGUI::OnMouseLeave(HWND hWnd)
	{
		g_bMouseDown = false;
		//g_nRotationAngle = 0;
	}

#if 0
	//Implementation using atomic variable and spin lock
	bool RubiksCubeSolverGUI::activateRenderingThread(bool force /*= false*/)
	{
		if (!graphicsAreaCreated_)
			return false;

		if (force)
		{
			//Activate always
			renderNow_.store(true, std::memory_order_release);
		}
		else
		{
			bool expected = false;
			bool desired = true;
			if (!renderNow_.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
			{
				//Display a message box
				if(firstGenCommand_ != firstGenerationCommands::eNoCommand)
				{
					RubiksCubeSolverUtils::CreateOkDialog("Another command is in progress. Please have patience!");
					return false;
				}
			}
		}

		return true;
	}
#endif

	bool RubiksCubeSolverGUI::activateRenderingThread(bool force /*= false*/)
	{
		if (!graphicsAreaCreated_)
			return false;

		//Display a message box
		if (!force && firstGenCommandInProgress_ != firstGenerationCommands::eNoCommand)
		{
			RubiksCubeSolverUtils::CreateOkDialog("Another command is in progress. Please have patience!");
			return false;
		}

		setReadySynchronously(true);
		cv_.notify_all();

		return true;
	}

	void RubiksCubeSolverGUI::commandHandlerFirstGen()
	{
		switch (firstGenCommand_)
		{
		case firstGenerationCommands::eSetCubeType:
			SetCubeTypeImpl();
			break;
		case firstGenerationCommands::eScramble:
			ScrambleImpl();
			break;
		case firstGenerationCommands::eSolve:
			unsigned int solutionSteps;
			unsigned long long duration;
			Solve(solutionSteps, duration);
			break;
		case firstGenerationCommands::eRunTests:
			runRubiksCubeTests();
			break;
		case firstGenerationCommands::eResizeRubiksCube:
			replaceModelBy(currentModelName_, rubikCubeSize_, true);
			break;
		default:
			//do nothing
			break;
		}

		if (firstGenCommand_ != firstGenerationCommands::eNoCommand)
			redrawWindow();
		firstGenCommand_ = firstGenerationCommands::eNoCommand;
	}

	void RubiksCubeSolverGUI::commandHandlerSecondGen()
	{
		switch (secondGenCommand_)
		{
		case secondGenerationCommands::eResetRubiksCube:
		{
			bool animate = true;
			Reset(animate);
		}
		break;
		//case secondGenerationCommands::eFitToScreen:
		//	fitToScreenImpl();
		//	break;
		default:
			//do nothing
			break;
		}

		if (secondGenCommand_ != secondGenerationCommands::eNoCommand)
			redrawWindow();
		secondGenCommand_ = secondGenerationCommands::eNoCommand;
	}

	void RubiksCubeSolverGUI::commandHandlerThirdGen()
	{
		switch (thirdGenCommand_)
		{
		case thirdGenerationCommands::eResizeWindow:
		{
			OnSizeImpl();
		}
		break;

		default:
			//do nothing
			break;
		}

		if(thirdGenCommand_ != thirdGenerationCommands::eNoCommand)
			redrawWindow();
		thirdGenCommand_ = thirdGenerationCommands::eNoCommand;
	}

	void RubiksCubeSolverGUI::Reset(bool animate)
	{
		scene_.Reset(animate);
	}

	void RubiksCubeSolverGUI::SetCubeTypeImpl()
	{
		switch (cubeType_)
		{
		case cubeType::rubiksCube:
			scene_.activateRubiksCube(rubikCubeSize_);
			break;
		case cubeType::mirrorCube:
			scene_.activateMirrorCube();
			break;
		}
	}

	void RubiksCubeSolverGUI::ScrambleImpl()
	{
		string invalidStep;
		if(!scene_.scramble(scramblingAlgo_, animate_, invalidStep))
			RubiksCubeSolverUtils::CreateOkDialog("Invalid scrambling algo! The below step is not a valid step.\nInvalid step: " + invalidStep);

		//if (!animate_)
		//	redrawWindow();
		//if(!animate_)
		displayUpdatedStats();
	}

	string RubiksCubeSolverGUI::Solve(unsigned int& solutionSteps, unsigned long long& duration)
	{
		string solution = scene_.Solve(solutionSteps, duration, animate_);
		//displayUpdatedStats();
		if (!animate_ && !devTestingMode_)
		{
			displayUpdatedStats();
			bool userDecision = RubiksCubeSolverUtils::CreateYesNoDialog("Do you want to see the animation of the solution?");
			if (userDecision)
			{
				unsigned int size;
				unsigned int scramblingSteps;
				string scramblingAlgo;
				string solution2;
				unsigned long long duration;
				string status;
				scene_.getUpdatedStats(size, scramblingSteps, scramblingAlgo, solutionSteps, solution, duration, status);

				//Go back to original position
				bool animate = false;
				string invalidStep;
				scene_.scramble(scramblingAlgo, animate, invalidStep);
				animate = true;
				solution2 = scene_.Solve(solutionSteps, duration, animate);
				//TODO: Check if both solutions are same or not
				//assert(solution == solution2); //sometimes initial steps may not be same
			}
		}
		return solution;
	}

	void RubiksCubeSolverGUI::runRubiksCubeTests()
	{
		devTestingMode_ = true;
		tester_.testRubiksCube(animate_);
		devTestingMode_ = false;
	}

}

