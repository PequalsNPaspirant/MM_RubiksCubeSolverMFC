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
using namespace std;

#include "RubiksCubeSolverGUI.h"
#include "RubiksCubeSolverUtils.h"
#include "RubiksCubeSolverTest.h"

//int APIENTRY _tWinMain(HINSTANCE hInstance,
//                     HINSTANCE hPrevInstance,
//                     LPTSTR    lpCmdLine,
//                     int       nCmdShow)
//{
//	UNREFERENCED_PARAMETER(hPrevInstance);
//	UNREFERENCED_PARAMETER(lpCmdLine);
//
//	mm::RubiksCubeSolverGUI& mainWindow = mm::RubiksCubeSolverGUI::getInstance();
//	mainWindow.createMainWindow(hInstance);
//	//mainWindow.createGraphicsArea();
//	std::thread t1(&mm::RubiksCubeSolverGUI::render, &mainWindow);
//	//t1.detach();
//
//	mainWindow.createMessageWindow();
//	WPARAM w = mainWindow.enterMainLoop();
//	t1.join();
//	return int(w);
//}

namespace mm {

	//RubiksCubeSolverGUI& RubiksCubeSolverGUI::getInstance()
	//{
	//	static RubiksCubeSolverGUI mainWindow;
	//	return mainWindow;
	//}

	RubiksCubeSolverGUI::RubiksCubeSolverGUI()
		: //WND_WIDTH(800),
		//WND_HEIGHT(800),
		//messageWndHeight(110),
		currentModelName_("RubiksCubeModel_v7"),
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
		scene_(*this, "RubiksCubeModel_v7", 3),
		framesPerRotation_(20), //moderate
		sleepTimeMilliSec_(20), //moderate
		rubiksCubeSize_(3),
		tester_(*this)
		//selMenuAnimationSpeed(ID_ANIMATIONSPEED_MODERATE),
		//selMenuRubiksCubeSize(ID_RUBIK_3X3X3)
	{
		
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

	void RubiksCubeSolverGUI::exitUI()
	{
		interruptAnimation_ = true;
		keepRunning = false;
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
		activateRenderingThread();
	}

	

	//void RubiksCubeSolverGUI::createMainWindow(HINSTANCE hInstance)
	//{
	//	LoadString(NULL, IDS_APP_TITLE, g_szTitle, MAX_LOADSTRING);
	//	LoadString(NULL, IDC_RUBIKSCUBE, g_szWindowClass, MAX_LOADSTRING);
	//	g_hAccelTable = LoadAccelerators(g_hInstance, MAKEINTRESOURCE(IDC_RUBIKSCUBE));
	//	g_hInstance = hInstance;

	//	WNDCLASSEX wndClass;

	//	g_hArrow = LoadCursor(NULL, IDC_ARROW);
	//	g_hHand = LoadCursor(NULL, IDC_HAND);

	//	ZeroMemory(&wndClass, sizeof(WNDCLASSEX));
	//	wndClass.cbSize = sizeof(WNDCLASSEX);
	//	wndClass.style = CS_HREDRAW | CS_VREDRAW;
	//	wndClass.lpfnWndProc = RubiksCubeSolverGUI::WndProcCallback;
	//	wndClass.hInstance = hInstance;
	//	wndClass.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_RUBIKSCUBE));
	//	wndClass.hCursor = g_hArrow;
	//	wndClass.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	//	wndClass.lpszClassName = g_szWindowClass;
	//	wndClass.cbClsExtra = 0;
	//	wndClass.lpszMenuName = MAKEINTRESOURCE(IDC_RUBIKSCUBE);
	//	wndClass.hIconSm = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_RUBIKSCUBE));

	//	RegisterClassEx(&wndClass);

	//	DWORD dwStyle = 0;

	//	if (g_bFullScreen)
	//	{
	//		dwStyle = WS_POPUP | WS_CLIPSIBLINGS | WS_CLIPCHILDREN;
	//		changeToFullScreen();
	//	}
	//	else
	//	{
	//		//dwStyle = WS_OVERLAPPEDWINDOW | WS_CLIPSIBLINGS | WS_CLIPCHILDREN | WS_MAXIMIZE;
	//		dwStyle = WS_OVERLAPPEDWINDOW | WS_CLIPSIBLINGS | WS_CLIPCHILDREN;
	//	}

	//	RECT rWnd;
	//	rWnd.left = 0;
	//	rWnd.right = WND_WIDTH;
	//	rWnd.top = 0;
	//	rWnd.bottom = WND_HEIGHT;

	//	AdjustWindowRectEx(&rWnd, dwStyle, true, 0);

	//	//Create menu beore creating main window
	//	createMenu();

	//	g_hWnd = CreateWindowEx(0,
	//		g_szWindowClass,
	//		g_szTitle,
	//		dwStyle,
	//		CW_USEDEFAULT, CW_USEDEFAULT,
	//		WND_WIDTH,
	//		WND_HEIGHT, // - messageWndHeight,
	//		NULL, hMenu, g_hInstance, NULL);

	//	if (!g_hWnd) return;

	//	ShowWindow(g_hWnd, SW_SHOWNORMAL);
	//	//ShowWindow(hWnd, SW_MAXIMIZE);
	//	UpdateWindow(g_hWnd);

	//	SetFocus(g_hWnd);
	//}

	//bool RubiksCubeSolverGUI::changeToFullScreen()
	//{
	//	DEVMODE dmSettings;	// Device Mode variable

	//	memset(&dmSettings, 0, sizeof(dmSettings)); // Makes Sure Memory's Cleared

	//												// Get current settings -- This function fills our the settings
	//												// This makes sure NT and Win98 machines change correctly
	//	if (!EnumDisplaySettings(NULL, ENUM_CURRENT_SETTINGS, &dmSettings))
	//		return false;

	//	// selected screen width and height
	//	dmSettings.dmPelsWidth = WND_WIDTH;
	//	dmSettings.dmPelsHeight = WND_HEIGHT;

	//	// This function actually changes the screen to full screen
	//	// CDS_FULLSCREEN Gets Rid Of Start Bar.
	//	// We always want to get a result from this function to check if we failed
	//	int result = ChangeDisplaySettings(&dmSettings, CDS_FULLSCREEN);

	//	// Check if we didn't recieved a good return message From the function
	//	if (result != DISP_CHANGE_SUCCESSFUL)
	//		return false;

	//	return true;
	//}

	//void RubiksCubeSolverGUI::createMenu()
	//{
	//	hMenu = LoadMenu(g_hInstance, MAKEINTRESOURCE(IDC_RUBIKSCUBE));

	//	//Uncheck all menu items
	//	//HMENU hmenu = GetMenu(g_hWnd);
	//	CheckMenuItem(hMenu, ID_ANIMATIONSPEED_VERYSLOW, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_ANIMATIONSPEED_SLOW, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_ANIMATIONSPEED_MODERATE, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_ANIMATIONSPEED_FAST, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_ANIMATIONSPEED_VERYFAST, MF_UNCHECKED | MF_ENABLED);

	//	CheckMenuItem(hMenu, ID_RUBIK_1X1X1, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_2X2X2, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_3X3X3, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_4X4X4, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_5X5X5, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_6X6X6, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_7X7X7, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_8X8X8, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_9X9X9, MF_UNCHECKED | MF_ENABLED);
	//	CheckMenuItem(hMenu, ID_RUBIK_10X10X10, MF_UNCHECKED | MF_ENABLED);

	//	//Check defalt menu selection on start
	//	CheckMenuItem(hMenu, selMenuAnimationSpeed, MF_CHECKED | MF_UNHILITE);
	//	//CheckMenuRadioItem(hMenu, ID_ANIMATIONSPEED_VERYSLOW, ID_ANIMATIONSPEED_VERYFAST, selMenuAnimationSpeed, MF_BYCOMMAND);
	//	CheckMenuItem(hMenu, selMenuRubiksCubeSize, MF_CHECKED | MF_UNHILITE);
	//}

	void RubiksCubeSolverGUI::render()
	{
		createGraphicsArea();

		while (keepRunning)
		{
			//bool expected = true;
			//bool desired = false;
			//if (renderNow_.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
			if(renderNow_.load(std::memory_order_acquire))
			{
				try 
				{
					commandHandlerFirstGen();
				}
				catch (bool flag)
				{
					//The previous command is broken/interrupted, reset the rubik cube.
					//If we need to perform different actions on the type of interrupt, we can include that information in the exception object
					//bool animate = true;
					//Reset(animate);
					interruptAnimation_ = false; //reset the flag
					firstGenCommand_ = firstGenerationCommands::eNoCommand;
					//commandHandlerSecondGen();
				}
				commandHandlerSecondGen();
				renderNow_.store(false, std::memory_order_release);

				//GdiFlush();
				//redrawWindow();

				//bool animate = true;
				//Scramble(animate);

				//std::thread t1(&RubiksCubeSolverGUI::Scramble, this, animate);
				//t1.detach();
			}
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

		Reset(true); //there is some bug in initial display: the Rubiks cube is displayed scattered. Reset() is a workaround.
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

	//void RubiksCubeSolverGUI::createMessageWindow()
	//{
	//	//g_hWndList = GetDlgItem(g_hWnd, IDC_LIST2);
	//	//
	//	//ScrollBar_Show(g_hWndList, TRUE);
	//	//
	//	//ListBox_ResetContent(g_hWndList);

	//	int gapUnderScrollBar = 20;
	//	int top = 0;
	//	int left = 0;
	//	int width = WND_WIDTH;
	//	int height = messageWndHeight + gapUnderScrollBar;
	//	g_hWndMessage = CreateWindowEx(
	//		NULL,
	//		L"LISTBOX",
	//		NULL,
	//		WS_CHILD | WS_VISIBLE | WS_HSCROLL | WS_VSCROLL,
	//		left,
	//		top,
	//		width,
	//		height,
	//		g_hWnd,
	//		NULL,
	//		GetModuleHandle(NULL),
	//		NULL);
	//	//ListBox_AddString(g_hWndMessage, L"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
	//	ListBox_Enable(g_hWndMessage, TRUE);
	//	ListBox_SetHorizontalExtent(g_hWndMessage, width);
	//	g_hDCMessage = GetDC(g_hWndMessage);

	//	//https://docs.microsoft.com/en-us/windows/desktop/api/wingdi/nf-wingdi-createfonta
	//	//HFONT hFont = CreateFont(15, 0, 0, 0, /*FW_BOLD*/FW_LIGHT, 0, 0, 0, 0, 0, 0, 2, FF_SCRIPT, L"SYSTEM_FIXED_FONT");
	//	charHeight = 15;
	//	HFONT hFont = CreateFont(
	//		charHeight, //cHeight
	//		0, //cWidth
	//		0, //cEscapement
	//		0, //cOrientation
	//		300, //cWeight
	//		true, //bItalic
	//		false, //bUnderline
	//		false,  //bStrikeOut
	//		DEFAULT_CHARSET, //iCharSet
	//		OUT_DEFAULT_PRECIS, //iOutPrecision
	//		CLIP_DEFAULT_PRECIS, //iClipPrecision
	//		DEFAULT_QUALITY, //iQuality
	//		DEFAULT_PITCH, //iPitchAndFamily
	//		L"Arial" //pszFaceName
	//	);
	//	charWidth = 20;
	//	HFONT hTmp = (HFONT)SelectObject(g_hDCMessage, hFont);
	//	SendMessage(g_hWndMessage, WM_SETFONT, (WPARAM)hFont, 0);

	//	SetBkMode(g_hDCMessage, TRANSPARENT);
	//	//SetBkColor(g_hDCMessage, RGB(255, 0, 0));
	//	//SetTextColor(g_hDCMessage, RGB(50, 50, 50));
	//	SetTextColor(g_hDCMessage, RGB(0, 0, 255));

	//	displayMessage();
	//}

	void RubiksCubeSolverGUI::redrawWindow()
	{
		scene_.renderScene();
		bool result = SwapBuffers(g_hDC);
	}

	//bool RubiksCubeSolverGUI::CreateYesNoDialog(const string& message)
	//{
	//	wstring wMessage(message.begin(), message.end());
	//	if (MessageBox(g_hWnd, wMessage.c_str(), g_szTitle, MB_YESNO | MB_ICONQUESTION | MB_APPLMODAL) == IDYES)
	//		return true;
	//	else
	//		return false;
	//}

	//void RubiksCubeSolverGUI::CreateOkDialog(const string& message)
	//{
	//	wstring wMessage(message.begin(), message.end());
	//	MessageBox(g_hWnd, wMessage.c_str(), g_szTitle, MB_OK | MB_ICONINFORMATION | MB_APPLMODAL);
	//}

	void RubiksCubeSolverGUI::applyAlgorithm(const string& algo, bool animate)
	{
		string invalidStep;
		scene_.scramble(algo, animate, invalidStep);
	}

	unique_ptr<RubiksCubeModel> RubiksCubeSolverGUI::replaceModelBy(const string& modelName, int size, bool animate)
	{
		unique_ptr<RubiksCubeModel> originalModel = scene_.replaceModelBy(modelName, size);
		//if(animate)
		//	activateRenderingThread();
			//redrawWindow();
		displayUpdatedStats();
		return std::move(originalModel);
	}

	unique_ptr<RubiksCubeModel> RubiksCubeSolverGUI::replaceModelBy(unique_ptr<RubiksCubeModel>&& newModel, bool animate)
	{
		unique_ptr<RubiksCubeModel> originalModel = scene_.replaceModelBy(std::move(newModel));
		//if (animate)
		//	activateRenderingThread();
			//redrawWindow();
		return std::move(originalModel);
	}

	bool RubiksCubeSolverGUI::isSolved()
	{
		return scene_.isSolved();
	}

	
	void RubiksCubeSolverGUI::getUpdatedStats(unsigned int& size, unsigned int& scramblingSteps, string& scramblingAlgo, unsigned int& solutionSteps, string& solution, unsigned long long& duration)
	{
		scene_.getUpdatedStats(size, scramblingSteps, scramblingAlgo, solutionSteps, solution, duration);
	}

	string getCommaSeparatedTimeDuration(unsigned long long duration)
	{
		string durationStr = "000,000.000,000,000";
		int pos = durationStr.length();
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

		//if (message.empty())
		//{
		unsigned int size;
		unsigned int scramblingSteps;
			string scramblingAlgo;
			unsigned int solutionSteps;
			string solution;
			unsigned long long duration;
			scene_.getUpdatedStats(size, scramblingSteps, scramblingAlgo, solutionSteps, solution, duration);
			RubiksCubeSolverUtils::displayMessage(scene_.getRubiksCubeSize(), scramblingSteps, scramblingAlgo, solutionSteps, solution, 
				"Time required to solve: " + getCommaSeparatedTimeDuration(duration));
		//}
		//else
		//{
		//	RubiksCubeSolverMainWindow::getInstance().displayMessage();

			//wstring wStrMsg(message.begin(), message.end());

			//RECT c = { 0, 0, 0, 0 };
			//DrawText(g_hDCMessage, wStrMsg.c_str(), wStrMsg.length(), &c, DT_CALCRECT);

			//int clearance = 20;
			////if (WND_WIDTH < c.right + clearance)
			//ListBox_SetHorizontalExtent(g_hWndMessage, c.right + clearance);

			//ListBox_ResetContent(g_hWndMessage);
			//displayMessage_currentLine(6, 4, message);
			//BOOL result = UpdateWindow(g_hWndMessage);
			//wstring wText(message.begin(), message.end());
			//RECT textRect;
			//textRect.left = 6;
			//textRect.right = WND_WIDTH;
			//textRect.top = 4;
			//textRect.bottom = messageWndHeight;
			//
			//SetBkMode(g_hDC, TRANSPARENT);
			//SetTextColor(g_hDC, RGB(50, 50, 50));
			//DrawText(g_hDC, wText.c_str(), -1, &textRect, DT_SINGLELINE | DT_NOCLIP);

			//The following lines refreshes the screen and the message is displayed on screen (dont know the reason)
			//HDC wdc = GetWindowDC(g_hWnd);
			//DeleteDC(wdc);
		//}
	}

	string RubiksCubeSolverGUI::generateScramblingAlgo(int length)
	{
		return scene_.generateScramblingAlgo(length);
	}

	//string getCommaSeparatedTimeDuration(unsigned long long duration)
	//{
	//	string durationStr = "000,000.000,000,000";
	//	int pos = durationStr.length() - 1;
	//	for (; pos > 0 && duration > 0; --pos)
	//	{
	//		if (durationStr[pos] == '0')
	//		{
	//			durationStr[pos] = '0' + duration % 10;
	//			duration /= 10;
	//		}
	//	}
	//	if (pos > 6)
	//		pos = 6;
	//	durationStr = durationStr.substr(pos);
	//	durationStr += " sec";

	//	return durationStr;
	//}

	//void RubiksCubeSolverGUI::displayMessage(int scramblingSteps, const string& scramblingAlgo, 
	//	int solutionSteps, const string& solution, unsigned long long duration)
	//{
	//	int left = 6;
	//	int top = 4;
	//	int lineHeight = 18;

	//	string size = to_string(scene_.getRubiksCubeSize());
	//	string rubikCubeSize("Rubik's Cube Size: " + size + "x" + size + "x" + size);
	//	string scramblingStepsStr("Scrambling Steps: " + (scramblingSteps > 0 ? to_string(scramblingSteps) : ""));
	//	string scrambleMsg("Scrambling Algorithm: " + scramblingAlgo);
	//	string solutionStepsStr("Solution Steps: " + (solutionSteps > 0 ? to_string(solutionSteps) : ""));
	//	string solutionMsg("Solution Algorithm: " + solution);
	//	string durationStr;
	//	if(duration > 0)
	//	{
	//		durationStr = getCommaSeparatedTimeDuration(duration);
	//	}
	//	string timeToSolveMsg("Time required to solve: " + durationStr);

	//	//int horizontalExtent = ListBox_GetHorizontalExtent(g_hWndMessage);
	//	string* pStr = &rubikCubeSize;
	//	if (pStr->length() < scramblingStepsStr.length())
	//		pStr = &scramblingStepsStr;
	//	if (pStr->length() < scrambleMsg.length())
	//		pStr = &scrambleMsg;
	//	if (pStr->length() < solutionStepsStr.length())
	//		pStr = &solutionStepsStr;
	//	if (pStr->length() < solutionMsg.length())
	//		pStr = &solutionMsg;
	//	if (pStr->length() < timeToSolveMsg.length())
	//		pStr = &timeToSolveMsg;

	//	wstring wStrMsg(pStr->begin(), pStr->end());

	//	RECT c = { 0, 0, 0, 0 };
	//	DrawText(g_hDCMessage, wStrMsg.c_str(), wStrMsg.length(), &c, DT_CALCRECT);

	//	int clearance = 20;
	//	//if (WND_WIDTH < c.right + clearance)
	//	ListBox_SetHorizontalExtent(g_hWndMessage, c.right + clearance);

	//	ListBox_ResetContent(g_hWndMessage);
	//	int currentTop = top;
	//	displayMessage_currentLine(left, currentTop, std::move(rubikCubeSize));

	//	currentTop += lineHeight;
	//	displayMessage_currentLine(left, currentTop, std::move(scramblingStepsStr));

	//	currentTop += lineHeight;
	//	displayMessage_currentLine(left, currentTop, std::move(scrambleMsg));

	//	currentTop += lineHeight;
	//	displayMessage_currentLine(left, currentTop, std::move(solutionStepsStr));

	//	currentTop += lineHeight;
	//	displayMessage_currentLine(left, currentTop, std::move(solutionMsg));

	//	currentTop += lineHeight;
	//	displayMessage_currentLine(left, currentTop, std::move(timeToSolveMsg));

	//	//SendMessage(g_hWndMessage, WM_SETREDRAW, TRUE, 0);
	//	BOOL result = UpdateWindow(g_hWndMessage);
	//	//RECT d;
	//	//GetClientRect(g_hWndMessage, &d);
	//	//InvalidateRect(g_hWndMessage, &d, TRUE);
	//}

	//void RubiksCubeSolverGUI::displayMessage_currentLine(int left, int top, const string& line)
	//{
	//	wstring wText(line.begin(), line.end());
	//	RECT textRect;
	//	textRect.left = left;

	//	int WND_WIDTH = 100; //TODO: temp fix

	//	textRect.right = WND_WIDTH;
	//	textRect.top = top;
	//	textRect.bottom = messageWndHeight;

	//	//HDC hDC = GetWindowDC(g_hWnd);
	//	//HFONT hFont = CreateFont(15, 0, 0, 0, /*FW_BOLD*/FW_LIGHT, 0, 0, 0, 0, 0, 0, 2, FF_SCRIPT, L"SYSTEM_FIXED_FONT");
	//	//HFONT hTmp = (HFONT)SelectObject(g_hDC, hFont);
	//	//LOGFONT logfont;
	//	//GetObject(hFont, sizeof(LOGFONT), &logfont);
	//	//logfont.lfHeight = logfont.lfHeight / 2;

	//	//HFONT hNewFont = CreateFontIndirect(&logfont);
	//	//HFONT hOldFont = (HFONT)SelectObject(hDC, hNewFont);

	//	HDC hdc = GetDC(g_hWndMessage);
	//	SetTextColor(hdc, RGB(0, 0, 255));
	//	SetTextColor(g_hDCMessage, RGB(0, 0, 255));
	//	//HDC wdc = GetWindowDC(g_hWnd);
	//	//g_hDC = GetWindowDC(g_hWnd);
	//	//SetBkColor(wdc, RGB(255, 255, 255));
	//	//SetBkMode(g_hDC, TRANSPARENT);
	//	////SetTextColor(g_hDC, RGB(255, 255, 255));
	//	//SetTextColor(g_hDC, RGB(50, 50, 50));
	//	//https://docs.microsoft.com/en-us/windows/desktop/api/winuser/nf-winuser-drawtext
	//	//DrawText(g_hDCMessage, wText.c_str(), -1, &textRect, DT_SINGLELINE | DT_LEFT);
	//	ListBox_AddString(g_hWndMessage, wText.c_str());
	//	//DeleteDC(wdc);

	//	// Always select the old font back into the DC
	//	//SelectObject(hDC, hOldFont);
	//	//DeleteObject(hNewFont);
	//	DeleteDC(hdc);
	//}

	//WPARAM RubiksCubeSolverGUI::enterMainLoop()
	//{
	//	MSG msg;

	//	while (1)
	//	{
	//		// Check if there was a message
	//		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
	//		{
	//			if (!TranslateAccelerator(msg.hwnd, g_hAccelTable, &msg))
	//			{
	//				// if the message wasn't to quit
	//				if (msg.message == WM_QUIT)
	//					break;

	//				// render the scene first if the 
	//				// cube state was changed
	//				else if (msg.message == RC_CHANGED)
	//				{
	//					redrawWindow();
	//				}

	//				// find out what the message does
	//				TranslateMessage(&msg);

	//				// execute the message
	//				DispatchMessage(&msg);
	//			}
	//		}
	//		// if no message
	//		else
	//		{
	//			// Do computationally expensive things here.  We want to render the scene
	//			// every frame, so we call our rendering function here.  Even though the scene
	//			// doesn't change, it will bottle neck the message queue if we don't do something.
	//			// Usually WaitMessage() is used to make sure the app doesn't eat up the CPU.
	//			//redrawWindow();
	//			WaitMessage();
	//			//::Sleep(10);
	//		}
	//	}

	//	return msg.wParam;
	//}

	// ======================== Callback functions ======================== //

	//LRESULT CALLBACK RubiksCubeSolverGUI::WndProcCallback(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
	//{
	//	RubiksCubeSolverGUI& mainWindow = RubiksCubeSolverGUI::getInstance();

	//	switch (uMsg)
	//	{
	//		HANDLE_MSG(hWnd, WM_COMMAND, mainWindow.OnCommand);
	//		HANDLE_MSG(hWnd, RC_CHANGED, mainWindow.OnRubiksCubeChanged);
	//		//HANDLE_MSG(hWnd, WM_PAINT, mainWindow.OnPaint);
	//		HANDLE_MSG(hWnd, WM_LBUTTONDOWN, mainWindow.OnLButtonDown);
	//		HANDLE_MSG(hWnd, WM_LBUTTONUP, mainWindow.OnLButtonUp);
	//		HANDLE_MSG(hWnd, WM_MOUSEMOVE, mainWindow.OnMouseMove);
	//		HANDLE_MSG(hWnd, WM_SIZE, mainWindow.OnSize);
	//		HANDLE_MSG(hWnd, WM_MOUSELEAVE, mainWindow.OnMouseLeave);
	//		HANDLE_MSG(hWnd, WM_MOUSEWHEEL, mainWindow.OnMouseWheel);
	//		HANDLE_MSG(hWnd, WM_DESTROY, mainWindow.OnDestroy);
	//		//HANDLE_MSG(hWnd, WM_ERASEBKGND, mainWindow.OnEraseBackground);
	//	default:
	//		return DefWindowProc(hWnd, uMsg, wParam, lParam);
	//	}
	//}

	//LRESULT CALLBACK RubiksCubeSolverGUI::AboutProcCallback(HWND hwndDlg, UINT uMsg, WPARAM wParam, LPARAM lParam)
	//{
	//	switch (uMsg)
	//	{
	//	case WM_COMMAND:
	//		switch (wParam)
	//		{
	//		case IDOK:
	//			EndDialog(hwndDlg, 0);
	//			return true;
	//		}
	//		break;
	//	}

	//	return false;
	//}

	//void RubiksCubeSolverGUI::deInit()
	//{
	//	if (g_hRC)
	//	{
	//		wglMakeCurrent(NULL, NULL);
	//		wglDeleteContext(g_hRC);
	//	}

	//	if (g_hDC)
	//		ReleaseDC(g_hWnd, g_hDC);

	//	if (g_bFullScreen)
	//		ChangeDisplaySettings(NULL, 0);

	//	UnregisterClass(g_szWindowClass, g_hInstance);
	//	PostQuitMessage(0);
	//}

	// ============================== Message Handling ================================= //

	void RubiksCubeSolverGUI::OnRubiksCubeChanged(HWND hWnd)
	{
		// check for solution
		//if (scene_.g_cCube.IsSolved())
		{
			TCHAR solvedMsg[MAX_LOADSTRING];
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
	void RubiksCubeSolverGUI::OnMouseMove(int rotate, int tilt)
	{
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

	//  Process WM_MOUSEWHEEL message for window/dialog: 
	void RubiksCubeSolverGUI::OnMouseWheel(float distance)
	{
		//int rotations = zDelta / WHEEL_DELTA;
		//float distance = rotations * 1.0f;

		scene_.g_cCamera.Move(distance);
		//firstGenCommand_ = firstGenerationCommands::eNoCommand;
		activateRenderingThread(true);
		//redrawWindow();
	}

	//  Process WM_SIZE message for window/dialog: 
	void RubiksCubeSolverGUI::OnSize(int cx, int cy)
	{
		//if (!g_bFullScreen)
		//{
		//	WND_WIDTH = cx;
		//	WND_HEIGHT = cy;
			scene_.sizeOpenGlScreen(cx, cy);
		//	MoveWindow(g_hWndMessage, 0, 0, cx, messageWndHeight, false);
		//	GetClientRect(hWnd, &g_rWnd);
		//}
		//firstGenCommand_ = firstGenerationCommands::eNoCommand;
		//std::this_thread::sleep_for(std::chrono_duration<std::chrono::milliseconds>(10));
		activateRenderingThread(true);
		//redrawWindow();
	}

	void RubiksCubeSolverGUI::OnMouseLeave(HWND hWnd)
	{
		g_bMouseDown = false;
		//g_nRotationAngle = 0;
	}

	//  Process WM_COMMAND message for window/dialog: 
	//void RubiksCubeSolverGUI::OnCommand(HWND hwnd, int id, HWND hwndCtl, UINT codeNotify)
	//{
	//	commandId_ = id;
	//	commandGeneration_ = 1;
	//	//Uncheck last selected menu items
	//	HMENU hMenu = GetMenu(hwnd);
	//	if (commandId_ >= ID_RUBIK_1X1X1 && commandId_ <= ID_RUBIK_10X10X10 && commandId_ != selMenuRubiksCubeSize)
	//	{
	//		commandGeneration_ = 1;
	//		CheckMenuItem(hMenu, selMenuRubiksCubeSize, MF_UNCHECKED | MF_ENABLED);
	//		selMenuRubiksCubeSize = commandId_;
	//		CheckMenuItem(hMenu, selMenuRubiksCubeSize, MF_CHECKED | MFS_GRAYED);
	//	}
	//	else if (commandId_ >= ID_ANIMATIONSPEED_VERYSLOW && commandId_ <= ID_ANIMATIONSPEED_VERYFAST && commandId_ != selMenuAnimationSpeed)
	//	{
	//		commandGeneration_ = 2;
	//		CheckMenuItem(hMenu, selMenuAnimationSpeed, MF_UNCHECKED | MF_ENABLED);
	//		selMenuAnimationSpeed = commandId_;
	//		CheckMenuItem(hMenu, selMenuAnimationSpeed, MF_CHECKED | MFS_GRAYED);
	//		switch (commandId_)
	//		{
	//		case ID_ANIMATIONSPEED_VERYSLOW: 
	//			framesPerRotation_ = 40;
	//			sleepTimeMilliSec_ = 40;
	//			break;
	//		case ID_ANIMATIONSPEED_SLOW:
	//			framesPerRotation_ = 30;
	//			sleepTimeMilliSec_ = 30;
	//			break;
	//		case ID_ANIMATIONSPEED_MODERATE:
	//			framesPerRotation_ = 20;
	//			sleepTimeMilliSec_ = 20;
	//			break;
	//		case ID_ANIMATIONSPEED_FAST:
	//			framesPerRotation_ = 10;
	//			sleepTimeMilliSec_ = 10;
	//			break;
	//		case ID_ANIMATIONSPEED_VERYFAST:
	//			framesPerRotation_ = 3;
	//			sleepTimeMilliSec_ = 3;
	//			break;
	//		}
	//	}
	//	else if (commandId_ == IDM_ABOUT)
	//	{
	//		commandGeneration_ = 2;
	//		DialogBox(g_hInstance, MAKEINTRESOURCE(IDD_ABOUTBOX),
	//			g_hWnd, reinterpret_cast<DLGPROC>(RubiksCubeSolverGUI::AboutProcCallback));
	//	}
	//	else if (commandId_ == IDM_EXIT)
	//	{
	//		commandGeneration_ = 2;
	//		PostQuitMessage(0);
	//	}
	//	else if (commandId_ == IDM_RUBIKSCUBE_RESET)
	//	{
	//		commandGeneration_ = 2;
	//		//Reset should force reset the Rubik's cube, even though other operations are in progress
	//		setResetRubiksCube(true);
	//	}

	//	activateRenderingThread();
	//}

	bool RubiksCubeSolverGUI::activateRenderingThread(bool force /*= false*/)
	{
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
				//if (commandGeneration_ == 1)
					//CreateOkDialog("Another animation is in progress. Please have patience!");
				if(firstGenCommand_ != firstGenerationCommands::eNoCommand)
				{
					RubiksCubeSolverUtils::CreateOkDialog("Another animation is in progress. Please have patience!");
					return false;
				}
			}
		}

			return true;
	}

	void RubiksCubeSolverGUI::commandHandlerFirstGen()
	{
		switch (firstGenCommand_)
		{
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

		redrawWindow();
		firstGenCommand_ = firstGenerationCommands::eNoCommand;

		// menu
		//if (commandId_ == IDM_RUBIKSCUBE_SCRAMBLE)
		//{
		//	ScrambleImpl();
		//}
		//else if (commandId_ == IDM_RUBIKSCUBE_SCRAMBLE_ANIM)
		//{
		//	bool animate = true;
		//	Scramble(animate);

		//	//bool expected = false;
		//	//bool desired = true;
		//	//if (isScrambling.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
		//	//{
		//	//	//GdiFlush();
		//	//	redrawWindow();
		//	//	bool animate = true;
		//	//	std::thread t1(&RubiksCubeSolverGUI::Scramble, this, animate);
		//	//	t1.detach();

		//	//	isScrambling.store(false, std::memory_order_release);
		//	//}

		//	//bool expected = false;
		//	//bool desired = true;
		//	//isScrambling.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed);
		//}
		//else if (commandId_ == IDM_RUBIKSCUBE_SOLVE)
		//{
		//	unsigned int solutionSteps;
		//	unsigned long long duration;
		//	bool askForAnimation = true;
		//	string solution = SolveOnCopy(solutionSteps, duration, askForAnimation);
		//}
		//else if (commandId_ == IDM_RUBIKSCUBE_SOLVE_ANIM)
		//{
		//	unsigned int solutionSteps;
		//	unsigned long long duration;
		//	bool animate = true;
		//	string solution = Solve(solutionSteps, duration, animate);
		//}
		//else if (commandId_ == ID_RUBIKSCUBE_TEST)
		//{
		//	bool animate = false;
		//	testRubiksCube(animate);
		//}
		//else if(commandId_ == ID_RUBIKSCUBE_TEST_ANIM)
		//{
		//	bool animate = true;
		//	testRubiksCube(animate);
		//}
		//else if (commandId_ == ID_RUBIKSCUBE_FITTOSCREEN)
		//{
		//	fitToScreen();
		//}
		//else if (commandId_ == ID_RUBIK_1X1X1)
		//{
		//	replaceModelBy(currentModelName_, 1, true);
		//}
		//else if (commandId_ == ID_RUBIK_2X2X2)
		//{
		//	replaceModelBy(currentModelName_, 2, true);
		//}
		//else if (commandId_ == ID_RUBIK_3X3X3)
		//{
		//	replaceModelBy(currentModelName_, 3, true);
		//}
		//else if (commandId_ == ID_RUBIK_4X4X4)
		//{
		//	replaceModelBy(currentModelName_, 4, true);
		//}
		//else if (commandId_ == ID_RUBIK_5X5X5)
		//{
		//	replaceModelBy(currentModelName_, 5, true);
		//}
		//else if (commandId_ == ID_RUBIK_6X6X6)
		//{
		//	replaceModelBy(currentModelName_, 6, true);
		//}
		//else if (commandId_ == ID_RUBIK_7X7X7)
		//{
		//	replaceModelBy(currentModelName_, 7, true);
		//}
		//else if (commandId_ == ID_RUBIK_8X8X8)
		//{
		//	replaceModelBy(currentModelName_, 8, true);
		//}
		//else if (commandId_ == ID_RUBIK_9X9X9)
		//{
		//	replaceModelBy(currentModelName_, 9, true);
		//}
		//else if (commandId_ == ID_RUBIK_10X10X10)
		//{
		//	replaceModelBy(currentModelName_, 10, true);
		//}
		//else if (commandId_ == ID_RUBIK_INCREASEBYONE)
		//{
		//	replaceModelBy(currentModelName_, scene_.getRubiksCubeSize() + 1, true);
		//}
		//else if (commandId_ == ID_RUBIK_INCREASEBYFIVE)
		//{
		//	replaceModelBy(currentModelName_, scene_.getRubiksCubeSize() + 5, true);
		//}
		//else if (commandId_ == ID_RUBIK_INCREASEBYTEN)
		//{
		//	replaceModelBy(currentModelName_, scene_.getRubiksCubeSize() + 10, true);
		//}
		//else if (commandId_ == IDM_RUBIKSCUBE_RESET) //This reset is called when no other command is in progress
		//{
		//	bool animate = true;
		//	Reset(animate);
		//}
		//else //The command may be mouse move, mouse wheel scroll, resize window etc
		//	redrawWindow();

		//commandId_ = -1; //Remove this reset of flag...it's bad design
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

		redrawWindow();
		secondGenCommand_ = secondGenerationCommands::eNoCommand;
	}

	void RubiksCubeSolverGUI::Reset(bool animate)
	{
		scene_.Reset(animate);
	}

	void RubiksCubeSolverGUI::ScrambleImpl()
	{
		//g_hDC = GetDC(g_hWnd);
		//g_hRC = wglCreateContext(g_hDC);
		//wglMakeCurrent(g_hDC, g_hRC);

		//string algo = scene_.generateScramblingAlgo(25);
		//wstring wAlgo(algo.begin(), algo.end());
		//wstring wMessage = L"Scramble using following Algorithm?";
		//wMessage = wMessage
		//	+ L"\nAlgorithm      : " + wAlgo
		//	+ L"\nNumber of steps: " + to_wstring(25);
		//if (MessageBox(g_hWnd, wMessage.c_str(),
		//	g_szTitle, MB_YESNO | MB_ICONQUESTION | MB_APPLMODAL) == IDYES)
		
		string invalidStep;
		if(!scene_.scramble(scramblingAlgo_, animate_, invalidStep))
			RubiksCubeSolverUtils::CreateOkDialog("Invalid scrambling algo! The below step is not a valid step.\nInvalid step: " + invalidStep);

		//if (!animate_)
		//	redrawWindow();
		if(!animate_)
			displayUpdatedStats();
	}

	//string RubiksCubeSolverGUI::SolveOnCopy(unsigned int& solutionSteps, unsigned long long& duration, bool askForAnimation)
	//{
	//	string solution1 = scene_.SolveOnCopy(solutionSteps, duration);
	//	//redrawWindow();
	//	//TODO: Display updated stats and solution on screen here

	//	if(!askForAnimation)
	//		return solution1;

	//	wstring wMessage = L"\nDo you want to see the animation of the solution?";
	//	//if (MessageBox(g_hWnd, wMessage.c_str(),
	//		g_szTitle, MB_YESNO | MB_ICONQUESTION | MB_APPLMODAL) == IDYES)
	//	{
	//		string solution2 = scene_.Solve(solutionSteps, duration, true);

	//		if (solution1 != solution2)
	//			//CreateOkDialog("Both the solutions are not same!");
	//			return "";
	//	}

	//	return solution1;
	//}

	string RubiksCubeSolverGUI::Solve(unsigned int& solutionSteps, unsigned long long& duration)
	{
		string solution = scene_.Solve(solutionSteps, duration, animate_);
		displayUpdatedStats();
		if (!animate_)
		{
			//displayUpdatedStats();
			bool userDecision = RubiksCubeSolverUtils::CreateYesNoDialog("Do you want to see the animation of the solution?");
			if (userDecision)
			{
				unsigned int size;
				unsigned int scramblingSteps;
				string scramblingAlgo;
				string solution2;
				unsigned long long duration;
				scene_.getUpdatedStats(size, scramblingSteps, scramblingAlgo, solutionSteps, solution, duration);

				//Go back to original position
				bool animate = false;
				string invalidStep;
				scene_.scramble(scramblingAlgo, animate, invalidStep);
				animate = true;
				solution2 = scene_.Solve(solutionSteps, duration, animate);
				//TODO: Check if both solutions are same or not
				assert(solution == solution2);
			}
		}
		return solution;
	}

	void RubiksCubeSolverGUI::runRubiksCubeTests()
	{
		tester_.testRubiksCube(animate_);
	}

	//void RubiksCubeSolverGUI::fitToScreenImpl()
	//{
	//	scene_.fitToScreen();
	//	//redrawWindow();
	//}
}

