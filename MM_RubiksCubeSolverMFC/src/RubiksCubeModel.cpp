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

#include "stdafx.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
using namespace std;

#include "Resource.h"
#include "RubiksCubeModel.h"

namespace mm {

	vector<GLuint> Textures::g_pTextures(7);

	void Textures::loadAllTextures()
	{
		loadTexture(IDB_WHITE, &g_pTextures[Color::White]);
		loadTexture(IDB_BLUE, &g_pTextures[Color::Blue]);
		loadTexture(IDB_ORANGE, &g_pTextures[Color::Orange]);
		loadTexture(IDB_RED, &g_pTextures[Color::Red]);
		loadTexture(IDB_GREEN, &g_pTextures[Color::Green]);
		loadTexture(IDB_YELLOW, &g_pTextures[Color::Yellow]);
		loadTexture(IDB_BLACK, &g_pTextures[Color::Black]);
	}

	void Textures::loadTexture(int nId, GLuint* texture)
	{
		// bitmap handle
		HBITMAP hBMP;

		// bitmap struct
		BITMAP   bmp;

		glGenTextures(1, texture);    // Create The Texture 
		hBMP = (HBITMAP)LoadImage(
			GetModuleHandle(NULL),
			MAKEINTRESOURCE(nId),
			//L"C:\\@_Programming\\Everything\\MM_RubiksCubeSolverMFC\\MM_RubiksCubeSolverMFC\\res\\blue.bmp",
			IMAGE_BITMAP, 0, 0, LR_CREATEDIBSECTION);

		if(hBMP == NULL)
		{
			DWORD lastError = GetLastError();
			return;
		}

		GetObject(hBMP, sizeof(bmp), &bmp);

		// Pixel Storage Mode (Word Alignment / 4 Bytes) 
		glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

		// bind to the texture ID
		glBindTexture(GL_TEXTURE_2D, *texture);

		glTexImage2D(
			GL_TEXTURE_2D,
			0,
			3,
			bmp.bmWidth, bmp.bmHeight,
			0,
			GL_BGR_EXT,
			GL_UNSIGNED_BYTE,
			bmp.bmBits
		);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		DeleteObject(hBMP);
	}

	void Textures::unloadAllTextures()
	{
		if (g_pTextures[Color::White])
			glDeleteTextures(1, &g_pTextures[Color::White]);

		if (g_pTextures[Color::Blue])
			glDeleteTextures(1, &g_pTextures[Color::Blue]);

		if (g_pTextures[Color::Orange])
			glDeleteTextures(1, &g_pTextures[Color::Orange]);

		if (g_pTextures[Color::Red])
			glDeleteTextures(1, &g_pTextures[Color::Red]);

		if (g_pTextures[Color::Green])
			glDeleteTextures(1, &g_pTextures[Color::Green]);

		if (g_pTextures[Color::Yellow])
			glDeleteTextures(1, &g_pTextures[Color::Yellow]);

		if (g_pTextures[Color::Black])
			glDeleteTextures(1, &g_pTextures[Color::Black]);
	}

	GLuint Textures::getTextureID(Color color)
	{
		return g_pTextures[color];
	}

	////Factory function declarations
	//unique_ptr<RubiksCubeModel> createRubiksCubeModel_v1(int size);
	//unique_ptr<RubiksCubeModel> createRubiksCubeModel_v2(int size);
	//unique_ptr<RubiksCubeModel> createRubiksCubeModel_v3(int size);
	//unique_ptr<RubiksCubeModel> createRubiksCubeModel_v4(int size);
	//unique_ptr<RubiksCubeModel> createRubiksCubeModel_v5(int size);
	//unique_ptr<RubiksCubeModel> createRubiksCubeModel_v6(int size);
	//unique_ptr<RubiksCubeModel> createRubiksCubeModel_v7(int size);

	//unordered_map<string, fptrRubiksCubeModelCreator> RubiksCubeModelFactory::rubiksCubeFactoryMap_ = {
	//	unordered_map<string, fptrRubiksCubeModelCreator>::value_type("RubiksCubeModel_v1", createRubiksCubeModel_v1),
	//	unordered_map<string, fptrRubiksCubeModelCreator>::value_type("RubiksCubeModel_v2", createRubiksCubeModel_v2),
	//	unordered_map<string, fptrRubiksCubeModelCreator>::value_type("RubiksCubeModel_v3", createRubiksCubeModel_v3),
	//	unordered_map<string, fptrRubiksCubeModelCreator>::value_type("RubiksCubeModel_v4", createRubiksCubeModel_v4),
	//	unordered_map<string, fptrRubiksCubeModelCreator>::value_type("RubiksCubeModel_v5", createRubiksCubeModel_v5),
	//	unordered_map<string, fptrRubiksCubeModelCreator>::value_type("RubiksCubeModel_v6", createRubiksCubeModel_v6),
	//	unordered_map<string, fptrRubiksCubeModelCreator>::value_type("RubiksCubeModel_v7", createRubiksCubeModel_v7)
	//};

	unique_ptr<RubiksCubeModel> RubiksCubeModelFactory::getRubiksCubeModel(const string& modelName, int size)
	{
		RubiksCubeFactoryMap& rubiksCubeFactoryMap = RubiksCubeFactoryMap::getRubiksCubeFactoryMap();
		fptrRubiksCubeModelCreator fPtr = rubiksCubeFactoryMap.getEntry(modelName);
		return (*fPtr)(size);
	}

	RubiksCubeFactoryMap& RubiksCubeFactoryMap::getRubiksCubeFactoryMap()
	{
		static RubiksCubeFactoryMap rubiksCubeFactoryMapSingleObject;
		return rubiksCubeFactoryMapSingleObject;
	}

	void RubiksCubeFactoryMap::addEntry(const string& modelName, fptrRubiksCubeModelCreator fptr)
	{
		rubiksCubeFactoryMap_[modelName] = fptr;
	}

	void RubiksCubeFactoryMap::removeEntry(const string& modelName)
	{
		size_t numElementsRemoved = rubiksCubeFactoryMap_.erase(modelName);
	}

	fptrRubiksCubeModelCreator RubiksCubeFactoryMap::getEntry(const string& modelName)
	{
		return rubiksCubeFactoryMap_[modelName];
	}

	RubiksCubeFactoryMap::RubiksCubeFactoryMap()
	{}

	RubiksCubeFactoryMap::~RubiksCubeFactoryMap()
	{}

	RegisterRubiksCubeFactoryFunction::RegisterRubiksCubeFactoryFunction(const string& modelName, fptrRubiksCubeModelCreator fptr)
	{
		RubiksCubeFactoryMap& rubiksCubeFactoryMap = RubiksCubeFactoryMap::getRubiksCubeFactoryMap();
		rubiksCubeFactoryMap.addEntry(modelName, fptr);
	}
}