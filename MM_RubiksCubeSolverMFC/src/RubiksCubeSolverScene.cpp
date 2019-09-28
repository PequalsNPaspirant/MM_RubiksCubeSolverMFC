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
#include <gl\gl.h>
#include <gl\glu.h>
#include <cassert>
using namespace std;

//#include "Resource.h"
#include "RubiksCubeSolverScene.h"
#include "RubiksCubeSolverUtils.h"

namespace mm {

	GLuint RubiksCubeSolverScene::g_pSelectBuffer[RubiksCubeSolverScene::SELECT_BUFFER_SIZE];
	const int RUBIKS_CUBE_SIZE = 3;

	RubiksCubeSolverScene::RubiksCubeSolverScene(RubiksCubeSolverGUI& refUI, const string& modelName, int size)
		: 
		rubicCubeModel_(RubiksCubeModelFactory::getRubiksCubeModel(modelName, size)),
		refUI_(refUI),
		rubiksCubeSize_(size)
	{
	}

	RubiksCubeSolverScene::~RubiksCubeSolverScene()
	{
		Textures::unloadAllTextures();
	}

	unique_ptr<RubiksCubeModel> RubiksCubeSolverScene::replaceModelBy(const string& modelName, int size)
	{
		unique_ptr<RubiksCubeModel> originalModel = std::move(rubicCubeModel_);
		rubicCubeModel_ = RubiksCubeModelFactory::getRubiksCubeModel(modelName, size);
		rubiksCubeSize_ = size;
		return std::move(originalModel);
	}

	unique_ptr<RubiksCubeModel> RubiksCubeSolverScene::replaceModelBy(unique_ptr<RubiksCubeModel>&& newModel)
	{
		unique_ptr<RubiksCubeModel> originalModel = std::move(rubicCubeModel_);
		rubicCubeModel_ = std::move(newModel);
		return std::move(originalModel);
	}

	void RubiksCubeSolverScene::initOpenGl(int nWidth, int nHeight)
	{
		glPolygonMode(GL_FRONT, GL_FILL);
		glClearColor(0.1f, 0.2f, 0.3f, 0.0f);
		//glClearColor(0.1f, 0.1f, 0.1f, 0.0f);

		//glFrontFace(GL_CCW);
		glEnable(GL_CULL_FACE); //This hides the mirror images of visible faces as those mirror images have face normal pointing into screen

		//glEnable(GL_COLOR_MATERIAL);
		//glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

		glClearDepth(1.0f);											// Depth Buffer Setup
		glEnable(GL_DEPTH_TEST);									// Enables Depth Testing
		glDepthFunc(GL_LEQUAL);										// The Type Of Depth Testing To Do
		glLineWidth(LINE_WIDTH);									// Set outline width

		//Load all textures after OpenGL context is loaded
		Textures::loadAllTextures();

		sizeOpenGlScreen(nWidth, nHeight);
	}

	void RubiksCubeSolverScene::sizeOpenGlScreen(int nWidth, int nHeight)
	{
		glViewport(0, 0, nWidth, nHeight);
		glMatrixMode(GL_PROJECTION);
		//glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		setFrustum(nWidth, nHeight);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}

	void RubiksCubeSolverScene::setFrustum(int nWidth, int nHeight)
	{
		GLdouble left, right;
		GLdouble top, bottom;
		GLdouble nearDist = 2;
		GLdouble farDist = 100000.0;

		if (nWidth < nHeight)
		{
			left = -1.0;
			right = 1.0;
			if (nWidth == 0)
				nWidth = 1;
			top = (double)nHeight / (double)nWidth;
			bottom = -top;

		}
		else
		{
			top = 1.0;
			bottom = -1.0;
			if (nHeight == 0)
				nHeight = 1;
			right = (double)nWidth / (double)nHeight;
			left = -right;
		}

		//glFrustum(left, right, bottom, top, 5.0, 100.0);
		glFrustum(left, right, bottom, top, nearDist, farDist);
	}

	void RubiksCubeSolverScene::initScene()
	{
		g_cCamera.SetDistance(45.0);
		g_cCamera.SetPhi((float)(PI / 4));
		g_cCamera.SetTheta((float)(PI / 4));
		//g_cCamera.SetPhi((float)(0.0));
		//g_cCamera.SetTheta((float)(0.0));
		//g_cCube.Randomize();

		fitToScreen();
	}

	void RubiksCubeSolverScene::renderScene()
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear The Screen And The Depth Buffer
		glLoadIdentity();									// Reset The View

		CVector3 pos = g_cCamera.GetEyePosition();
		CVector3 lookAt = g_cCamera.GetLookAt();
		CVector3 up = g_cCamera.GetUp();

		gluLookAt(
			pos.x, pos.y, pos.z,
			lookAt.x, lookAt.y, lookAt.z,
			up.x, up.y, up.z
		);

		glEnable(GL_LIGHTING);
		float color[] = { 1.0f, 1.0f, 1.0f, 1.0f };
		float* position = new float[4];

		g_cCamera.GetEyePosition().ToFloatArray(position);

		glLightfv(GL_LIGHT1, GL_AMBIENT, color);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, color);
		glLightfv(GL_LIGHT1, GL_POSITION, position);
		glLightfv(GL_LIGHT1, GL_SPECULAR, color);
		glLightfv(GL_LIGHT1, GL_SHININESS, color);

		delete[] position;
		position = NULL;

		glEnable(GL_LIGHT1);

		float shininess = 5.0f;
		glMaterialfv(GL_FRONT, GL_SPECULAR, color);
		//glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

		rubicCubeModel_->render();

		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHTING);

		//SwapBuffers(g_hDC);
	}

	INT RubiksCubeSolverScene::getSelectedObjects(int x, int y, int nWidth, int nHeight)
	{
		GLint viewport[4];

		glSelectBuffer(SELECT_BUFFER_SIZE, g_pSelectBuffer);
		glRenderMode(GL_SELECT);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();

		glGetIntegerv(GL_VIEWPORT, viewport);
		gluPickMatrix(x, viewport[3] - y, 1, 1, viewport);

		setFrustum(nWidth, nHeight);

		glMatrixMode(GL_MODELVIEW);

		renderScene();

		// restoring the original projection matrix
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glFlush();

		return glRenderMode(GL_RENDER);
	}

	CVector3 RubiksCubeSolverScene::mapCoordinates(int x, int y)
	{
		GLint viewport[4];
		GLdouble modelview[16];
		GLdouble projection[16];
		GLfloat winX, winY, winZ;
		GLdouble posX, posY, posZ;

		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		glGetIntegerv(GL_VIEWPORT, viewport);

		winX = (float)x;
		winY = (float)viewport[3] - (float)y;
		glReadPixels(x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);

		gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);

		return CVector3(posX, posY, posZ);
	}

	//void RubiksCubeSolverScene::getCubeSelection(int *x, int *y, int *z, Face *face, int g_nHitCount)
	//{
	//	GLuint names, *ptr, minZ, *ptrNames, numberOfNames;
	//	ptr = (GLuint *)g_pSelectBuffer;
	//	minZ = 0xffffffff;

	//	for (int i = 0; i < g_nHitCount && i < SELECT_BUFFER_SIZE; i++)
	//	{
	//		names = *ptr;
	//		ptr++;
	//		if (*ptr < minZ)
	//		{
	//			numberOfNames = names;
	//			minZ = *ptr;
	//			ptrNames = ptr + 2;
	//		}

	//		ptr += names + 2;
	//	}

	//	*x = ptrNames[0];
	//	*y = ptrNames[1];
	//	*z = ptrNames[2];
	//	*face = (Face)ptrNames[3];
	//}

	void RubiksCubeSolverScene::setAnimate(bool animate)
	{
		rubicCubeModel_->setAnimate(animate);
	}

	void RubiksCubeSolverScene::Reset(bool animate)
	{
		rubicCubeModel_->ResetCube(animate, &refUI_);
	}

	string RubiksCubeSolverScene::generateScramblingAlgo(int length)
	{
		bool includeNonStandardRotations = true;
		return rubicCubeModel_->generateScramblingAlgo(length, includeNonStandardRotations);
	}

	bool RubiksCubeSolverScene::pauseAnimation(bool pause)
	{
		return rubicCubeModel_->pauseAnimation(pause);
	}

	bool RubiksCubeSolverScene::scramble(const string& algo, bool animate, string& invalidStep)
	{
		return rubicCubeModel_->scramble(algo, animate, refUI_, invalidStep);
	}

	string RubiksCubeSolverScene::Solve(unsigned int& solutionSteps, unsigned long long& duration, bool animate)
	{
		string solution = rubicCubeModel_->solve(solutionSteps, duration, animate, refUI_);
		RubiksCubeSolverUtils::RunTimeAssert(rubicCubeModel_->isSolved());

		return solution;
	}

	bool RubiksCubeSolverScene::isSolved()
	{
		return rubicCubeModel_->isSolved();
	}

	void RubiksCubeSolverScene::fitToScreen()
	{
		//Set appropriate zoom factor
		//g_cCamera.SetDistance(45 + (10 * 2));
		float approxDistToFitScreen = rubiksCubeSize_ * 5;
		g_cCamera.SetDistance(approxDistToFitScreen);
		//Dim	width	actual dist		diff
		//10	20		- 65			45
		//4		08		- 30			22	
		//3		06		- 22			16
		//2		04		- 17			13
	}

	bool RubiksCubeSolverScene::activateRubiksCube()
	{
		return rubicCubeModel_->activateRubiksCube();
	}

	bool RubiksCubeSolverScene::activateMirrorCube()
	{
		return rubicCubeModel_->activateMirrorCube();
	}

	void RubiksCubeSolverScene::getUpdatedStats(unsigned int& size, unsigned int& scramblingSteps, string& scramblingAlgo,
		unsigned int& solutionSteps, string& solution, unsigned long long& duration, string& status)
	{
		rubicCubeModel_->getUpdatedStats(size, scramblingSteps, scramblingAlgo, solutionSteps, solution, duration, status);
	}
}

