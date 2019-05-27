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

#pragma once

#include <string>
using namespace std;

#include <gl\gl.h>
#include "Camera.h"
#include "RubiksCubeModel.h"

namespace mm {

	//forward declaration
	class RubiksCubeSolverGUI;

	class RubiksCubeSolverScene
	{
	public:
		RubiksCubeSolverScene(RubiksCubeSolverGUI& refUI, const string& modelName, int size);
		~RubiksCubeSolverScene();

		unique_ptr<RubiksCubeModel> replaceModelBy(const string& modelName, int size);
		unique_ptr<RubiksCubeModel> replaceModelBy(unique_ptr<RubiksCubeModel>&& newModel);

		/**/void initOpenGl(int nWidth, int nHeight);
		/**//**/void sizeOpenGlScreen(int nWidth, int nHeight);
		/**//**//**/void setFrustum(int nWidth, int nHeight);
		/**/void initScene();
		/**/void renderScene();	
		INT getSelectedObjects(int x, int y, int nWidth, int nHeight);
		//void getCubeSelection(int *x, int *y, int *z, Face *face, int g_nHitCount);
		CVector3 mapCoordinates(int x, int y);

		void setAnimate(bool animate);
		void Reset(bool animate);
		string generateScramblingAlgo(int length);
		bool scramble(const string& algo, bool animate, string& invalidStep);
		string Solve(unsigned int& solutionSteps, unsigned long long& duration, bool animate);
		//string SolveOnCopy(unsigned int& solutionSteps, unsigned long long& duration);
		bool isSolved();
		void fitToScreen();
		int getRubiksCubeSize() { return rubiksCubeSize_; }
		void getUpdatedStats(unsigned int& size, unsigned int& scramblingSteps, string& scramblingAlgo, unsigned int& solutionSteps, string& solution, unsigned long long& duration, string& status);
		bool pauseAnimation(bool pause);

		//deleted functions
		RubiksCubeSolverScene(const RubiksCubeSolverScene&) = delete;
		RubiksCubeSolverScene& operator=(const RubiksCubeSolverScene&) = delete;
		RubiksCubeSolverScene(RubiksCubeSolverScene&&) = delete;
		RubiksCubeSolverScene& operator=(RubiksCubeSolverScene&&) = delete;

	private:
		//const double CUBE_SIZE = 2.0;
		const float LINE_WIDTH = 2.0f;
		static const GLsizei SELECT_BUFFER_SIZE = 128;
		static GLuint g_pSelectBuffer[SELECT_BUFFER_SIZE];
		unique_ptr<RubiksCubeModel> rubicCubeModel_;
		RubiksCubeSolverGUI& refUI_;
		int rubiksCubeSize_;

	public:
		CCamera g_cCamera;
	};

}