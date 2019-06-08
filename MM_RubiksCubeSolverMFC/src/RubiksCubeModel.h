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
#include <vector>
#include <unordered_map>
#include <memory>
using namespace std;

#include <gl\gl.h>
#include <gl\glu.h>

namespace mm {

	enum Color
	{
		Yellow = 0,
		White = 1,
		Orange = 2,
		Red = 3,
		Blue = 4,
		Green = 5,

		Black = 6
	};

	class Textures
	{
	public:
		static void loadAllTextures();
		static void unloadAllTextures();
		static GLuint Textures::getTextureID(Color color);

	private:
		static void loadTexture(int nId, GLuint* texture);

		static vector<GLuint> g_pTextures;
	};

	class RubiksCubeSolverGUI;

	class RubiksCubeModel
	{
	public:
		virtual void setAnimate(bool animate) {}
		virtual void ResetCube(bool animate, RubiksCubeSolverGUI* ui) = 0;
		virtual void ResetCube(int size, double xt, double yt, double zt, bool animate, RubiksCubeSolverGUI* ui) {}
		//virtual int applyAlgorithm(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui) = 0;
		virtual string generateScramblingAlgo(int length, bool includeNonStandardRotations) = 0;
		virtual void scramble(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui) = 0;
		virtual bool scramble(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui, string& invalidStep) { return false; }
		virtual string solve(unsigned int& solutionSteps, unsigned long long& duration, bool animate, RubiksCubeSolverGUI& ui) = 0;
		virtual void render() = 0;
		virtual bool isSolved() = 0;

		virtual unique_ptr<RubiksCubeModel> copy() = 0;
		virtual string getModelName() = 0;
		virtual int getDimension() = 0;
		virtual bool activateRubiksCube() { return true; }
		virtual bool activateMirrorCube() { return true; }

		virtual void getUpdatedStats(unsigned int& size, unsigned int& scramblingSteps, string& scramblingAlgo, unsigned int& solutionSteps, string& solution, unsigned long long& duration, string& status) {}
		//virtual void setDisplayParameters(int scramblingSteps, const string& scramblingAlgo, int solutionSteps, const string& solution, unsigned long long duration) {}

		virtual bool isAlgoValid(const string& algo, string& invalidStep) { return false; }
		virtual bool pauseAnimation(bool pause) { return false; }
			
		virtual ~RubiksCubeModel() = 0
		{
		}
	};

	typedef unique_ptr<RubiksCubeModel> (*fptrRubiksCubeModelCreator)(int size);

	class RubiksCubeModelFactory
	{
	public:
		static unique_ptr<RubiksCubeModel> getRubiksCubeModel(const string& modelName, int size);

	};

	//a singlton class
	class RubiksCubeFactoryMap
	{
	public:
		static RubiksCubeFactoryMap& getRubiksCubeFactoryMap();
		void addEntry(const string& modelName, fptrRubiksCubeModelCreator fptr);
		void removeEntry(const string& modelName);
		fptrRubiksCubeModelCreator getEntry(const string& modelName);

		~RubiksCubeFactoryMap();
		RubiksCubeFactoryMap(const RubiksCubeFactoryMap&) = delete;
		RubiksCubeFactoryMap(RubiksCubeFactoryMap&&) = delete;
		RubiksCubeFactoryMap& operator=(const RubiksCubeFactoryMap&) = delete;		
		RubiksCubeFactoryMap& operator=(RubiksCubeFactoryMap&&) = delete;

	private:
		unordered_map<string, fptrRubiksCubeModelCreator> rubiksCubeFactoryMap_;

		RubiksCubeFactoryMap();
	};

	class RegisterRubiksCubeFactoryFunction
	{
	public:
		RegisterRubiksCubeFactoryFunction(const string& modelName, fptrRubiksCubeModelCreator fptr);
	};
}