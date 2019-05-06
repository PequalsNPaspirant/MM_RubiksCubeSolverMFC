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

#include "StdAfx.h"

#include <time.h>
#include <cassert>
#include <memory>
#include <chrono>
using namespace std;

//#include "Resource.h"
#include "RubiksCubeModel_v2.h"
#include "RubiksCubeSolverGUI.h"
#include "RubiksCubeSolverUtils.h"

namespace mm {

	//Factory function definition
	unique_ptr<RubiksCubeModel> createRubiksCubeModel_v2(int size)
	{
		return make_unique<RubiksCubeModel_v2>(size);
	}

	//Create a global object, so that its constructor is called before main and the factory map is initialized before main
	static RegisterRubiksCubeFactoryFunction object("RubiksCubeModel_v2", createRubiksCubeModel_v2);

	//==================== RubiksCubeModel_v2::Cube =========================

	const int RubiksCubeModel_v2::Cube::FACE_COUNT /* = 6*/;
	const double RubiksCubeModel_v2::CUBE_SIZE = 2.0;

	const RubiksCubeModel_v2::ColorRGB RubiksCubeModel_v2::ColorRGB::RGBColors[7] = {
		ColorRGB{ 255, 255, 0 },
		ColorRGB{ 255, 0, 0 },
		ColorRGB{ 0, 0, 255 },
		ColorRGB{ 0, 255, 0 },
		ColorRGB{ 255, 165, 0 },
		ColorRGB{ 255, 255, 255 },
		ColorRGB{ 0, 0, 0 }
	};

	RubiksCubeModel_v2::Cube::Cube(Color cTop, Color cBottom, Color cLeft, Color cRight, Color cFront, Color cBack, const Location& location)
		: faces_(FACE_COUNT)
	{
		faces_[Up] = cTop;
		faces_[Down] = cBottom;
		faces_[Left] = cLeft;
		faces_[Right] = cRight;
		faces_[Front] = cFront;
		faces_[Back] = cBack;

		location_ = location;
		//group_ = group;
	}

	RubiksCubeModel_v2::Cube::~Cube(void)
	{
	}

	Color RubiksCubeModel_v2::Cube::GetFaceColor(Face eFace) const
	{
		return faces_[eFace];
	}

	/*
	Make sure Cube is positioned in right way

	// Z-axis : Back -> Front // green  -> blue
	// X-axis : Left -> Right // orange -> red
	// Y-axis : Down -> Up    // white  -> yellow

		  yellow
			Y
			|
			. --> X red
		  /
		Z
	 blue

	*/

	void RubiksCubeModel_v2::Cube::rotate(CVector3 rotationAxis, double rotationAngle)
	{
		if (fabs(rotationAngle) < 0.00001)
			return;

		location_.rotate(rotationAxis, rotationAngle * PI / 180.0); //Angle should be in radians

		int numRotations = fabs(rotationAngle) / 90;
		if(rotationAxis == CVector3::XAxis)
		{
			while (--numRotations > -1)
			{
				if (rotationAngle > 0)
					TiltDown();
				else
					TiltUp();					
			}
		}
		else if (rotationAxis == CVector3::YAxis)
		{
			while (--numRotations > -1)
			{
				if (rotationAngle > 0)
					TurnRight();
				else
					TurnLeft();
			}
		}
		else if (rotationAxis == CVector3::ZAxis)
		{
			while (--numRotations > -1)
			{
				if (rotationAngle > 0)
					TiltLeft();
				else
					TiltRight();
			}
		}
	}

	// Aound X axis
	void RubiksCubeModel_v2::Cube::TiltUp()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Front];
		faces_[Front] = faces_[Down];
		faces_[Down] = faces_[Back];
		faces_[Back] = temp1;
	}

	// Aound X axis
	void RubiksCubeModel_v2::Cube::TiltDown()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Back];
		faces_[Back] = faces_[Down];
		faces_[Down] = faces_[Front];
		faces_[Front] = temp1;
	}

	//Around Y axis
	void RubiksCubeModel_v2::Cube::TurnLeft()
	{
		Color temp1 = faces_[Front];
		faces_[Front] = faces_[Right];
		faces_[Right] = faces_[Back];
		faces_[Back] = faces_[Left];
		faces_[Left] = temp1;
	}

	//Around Y axis
	void RubiksCubeModel_v2::Cube::TurnRight()
	{
		Color temp1 = faces_[Front];
		faces_[Front] = faces_[Left];
		faces_[Left] = faces_[Back];
		faces_[Back] = faces_[Right];
		faces_[Right] = temp1;
	}

	//Around Z axis
	void RubiksCubeModel_v2::Cube::TiltLeft()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Right];
		faces_[Right] = faces_[Down];
		faces_[Down] = faces_[Left];
		faces_[Left] = temp1;
	}

	//Around Z axis
	void RubiksCubeModel_v2::Cube::TiltRight()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Left];
		faces_[Left] = faces_[Down];
		faces_[Down] = faces_[Right];
		faces_[Right] = temp1;
	}

	bool RubiksCubeModel_v2::Cube::belongsTo(Face rotatingSection, int layerIndex, int size) const
	{
		if (rotatingSection == All)
			return true;

		RubiksCubeSolverUtils::RunTimeAssert(layerIndex > 0);

		double extend = (size - 1) / 2.0;
		/*
		Up = 0,
		Down = 1,
		Left = 2,
		Right = 3,
		Front = 4,
		Back = 5,
		*/
		static int diffData[6] = {-1, 1, 1, -1, -1, 1};
		double diff = diffData[rotatingSection] * (layerIndex - 1);

		double x, y, z;
		x = y = z = -1;
		bool retVal = false;
		switch (rotatingSection)
		{
		case Left:
			x = -extend + diff;
			retVal = (fabs(x - location_.x_) < 0.0001);
			break;
		case Right:
			x = extend + diff;
			retVal = (fabs(x - location_.x_) < 0.0001);
			break;
		case Down:
			y = -extend + diff;
			retVal = (fabs(y - location_.y_) < 0.0001);
			break;
		case Up:
			y = extend + diff;
			retVal = (fabs(y - location_.y_) < 0.0001);
			break;
		case Back:
			z = -extend + diff;
			retVal = (fabs(z - location_.z_) < 0.0001);
			break;
		case Front:
			z = extend + diff;
			retVal = (fabs(z - location_.z_) < 0.0001);
			break;
		}
		
		return retVal;
	}

	//==================== RubiksCubeModel_v2 =========================

	RubiksCubeModel_v2::RubiksCubeModel_v2(int size)
		: //cubes_(vector< vector< vector<Cube> > > (size, vector< vector<Cube> >(size, vector<Cube>(size)) ) ),
		//layerF_(vector< vector<Cube*> >(size, vector<Cube*>(size, nullptr))),
		//layerB_(vector< vector<Cube*> >(size, vector<Cube*>(size, nullptr))),
		//layerL_(vector< vector<Cube*> >(size, vector<Cube*>(size, nullptr))),
		//layerR_(vector< vector<Cube*> >(size, vector<Cube*>(size, nullptr))),
		//layerU_(vector< vector<Cube*> >(size, vector<Cube*>(size, nullptr))),
		//layerD_(vector< vector<Cube*> >(size, vector<Cube*>(size, nullptr))),
		size_(size)
		//g_bRotating(false),
		//g_bFlipRotation(false),
		//g_vRotationAxis(0, 0, 0),
		//g_nRotatingSection(None),
		//g_nRotationAngle(0)
	{
		ResetCube(false, nullptr);

		//for(int i = 0; i < size; ++i)
		//	for (int j = 0; j < size; ++j)
		//		layerF_[i][j] = cubes_[Location(i, j, size - 1)].get();

		//for (int i = 0; i < size; ++i)
		//	for (int j = 0; j < size; ++j)
		//		layerB_[i][j] = cubes_[Location(i, j, 0)].get();

		//for (int i = 0; i < size; ++i)
		//	for (int j = 0; j < size; ++j)
		//		layerL_[i][j] = cubes_[Location(0, i, j)].get();

		//for (int i = 0; i < size; ++i)
		//	for (int j = 0; j < size; ++j)
		//		layerR_[i][j] = cubes_[Location(size - 1, i, j)].get();

		//for (int i = 0; i < size; ++i)
		//	for (int j = 0; j < size; ++j)
		//		layerU_[i][j] = cubes_[Location(i, size - 1, j)].get();

		//for (int i = 0; i < size; ++i)
		//	for (int j = 0; j < size; ++j)
		//		layerD_[i][j] = cubes_[Location(i, 0, j)].get();

	}

	RubiksCubeModel_v2::RubiksCubeModel_v2(const RubiksCubeModel_v2& copy)
		: //cubes_(copy.cubes_),
		size_(copy.size_),
		g_bRotating(copy.g_bRotating),
		g_bFlipRotation(copy.g_bFlipRotation),
		g_vRotationAxis(copy.g_vRotationAxis),
		g_nRotatingSection(copy.g_nRotatingSection),
		g_nLayerIndex(copy.g_nLayerIndex),
		g_nRotationAngle(copy.g_nRotationAngle)
	{
		for (auto& obj : copy.cubes_)
		{
			cubes_[obj.first] = make_unique<Cube>(*obj.second.get());
		}
	}

	RubiksCubeModel_v2::~RubiksCubeModel_v2()
	{
	}

	void RubiksCubeModel_v2::ResetCube(bool animate, RubiksCubeSolverGUI* ui)
	{
		g_bRotating = false;
		g_bFlipRotation = false;
		g_vRotationAxis = CVector3(0, 0, 0);
		g_nRotatingSection = None;
		g_nLayerIndex = 1;
		g_nRotationAngle = 0;

		double extend = (size_ - 1) / 2.0;
		double x = -extend;
		for (int i = 0; i < size_; i++, ++x)
		{
			double y = -extend;
			for (int j = 0; j < size_; j++, ++y)
			{
				double z = -extend;
				for (int k = 0; k < size_; k++, ++z)
				{
					//int group = 0;
					//if (i == 0)
					//	group |= Groups::L;
					//else if (i == size_ - 1)
					//	group |= Groups::R;
					//
					//if (j == 0)
					//	group |= Groups::D;
					//else if (j == size_ - 1)
					//	group |= Groups::U;

					//if (k == 0)
					//	group |= Groups::B;
					//else if (k == size_ - 1)
					//	group |= Groups::F;

					Location loc(x, y, z);

					if (i == 0 || i == size_ - 1
						|| j == 0 || j == size_ - 1
						|| k == 0 || k == size_ - 1)
						cubes_[loc] = CreateCube(i, j, k, loc);
				}
			}
		}

		if (animate)
		{
			ui->redrawWindow();
		}
	}

	string RubiksCubeModel_v2::solve(unsigned int& solutionSteps, unsigned long long& duration, bool animate, RubiksCubeSolverGUI& ui)
	{
		RubiksCubeSolver solver(*this, animate, ui);
		using HRClock = std::chrono::high_resolution_clock;
		HRClock::time_point start_time = HRClock::now();
		string solution = solver.solve(solutionSteps);
		HRClock::time_point end_time = HRClock::now();
		std::chrono::nanoseconds time_span = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
		duration = time_span.count();

		return solution;
	}

	void RubiksCubeModel_v2::render()
	{
#ifdef _DEBUG
		// Draw Axis
		glBegin(GL_LINES);
		// x
		glColor3f(1.0f, 0.6f, 0.0f); // orange
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(CUBE_SIZE * 3, 0.0, 0.0);
		glColor3f(1.0f, 0.0f, 0.0f); // red
		glVertex3d(CUBE_SIZE * 3, 0.0, 0.0);
		glVertex3d(CUBE_SIZE * 4.5f, 0.0, 0.0);

		// y
		//glColor3f(0.0f, 1.0f, 0.0f);  // green
		glColor3f(1.0f, 1.0f, 1.0f);  // white
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, CUBE_SIZE * 3, 0.0);
		glColor3f(1.0f, 1.0f, 0.0f);  // yellow
		glVertex3d(0.0, CUBE_SIZE * 3, 0.0);
		glVertex3d(0.0, CUBE_SIZE * 4.5f, 0.0);

		// z
		glColor3f(0.0f, 1.0f, 0.0f);  // green
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, 0.0, CUBE_SIZE * 3);
		glColor3f(0.0f, 0.0f, 1.0f); // blue
		glVertex3d(0.0, 0.0, CUBE_SIZE * 3);
		glVertex3d(0.0, 0.0, CUBE_SIZE * 4.5f);
		glEnd();
#endif

		glInitNames();

		for (auto& obj : cubes_)
		{
			const Location& loc = obj.first;
			const Cube& cube = *obj.second.get();

			glPushMatrix();

			if (g_bRotating)
			{
				if (cube.belongsTo(g_nRotatingSection, g_nLayerIndex, size_))
				{
					int angle = g_bFlipRotation ? -g_nRotationAngle : g_nRotationAngle;
					glRotated(angle, g_vRotationAxis.x, g_vRotationAxis.y, g_vRotationAxis.z);
				}
			}

			renderIndividualCube(cube, cube.getLocation());

			glPopMatrix();
		}
	}

	void RubiksCubeModel_v2::renderIndividualCube(const RubiksCubeModel_v2::Cube& pCube, const RubiksCubeModel_v2::Location& location)
	{
		double x = location.x_;
		double y = location.y_;
		double z = location.z_;
		double extend = (size_ - 1) / 2.0;
		bool mirrorVisibleFaces = false;
		int offsetDist = 3;

		glPushName(x);
		glPushName(y);
		glPushName(z);

		// scale to -1 to +1
		//x--;
		//y--;
		//z--;

		glPushMatrix();

		glTranslated(x * CUBE_SIZE, y * CUBE_SIZE, z * CUBE_SIZE);

		Color top = pCube.GetFaceColor(Up);
		Color bottom = pCube.GetFaceColor(Down);
		Color left = pCube.GetFaceColor(Left);
		Color right = pCube.GetFaceColor(Right);
		Color back = pCube.GetFaceColor(Back);
		Color front = pCube.GetFaceColor(Front);

		glEnable(GL_TEXTURE_2D);

		// Front Face
		glPushName((GLuint)Front);
		glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(front));
		glBegin(GL_QUADS);
		ColorRGB colRgb = ColorRGB::RGBColors[front];
		glColor3ub(colRgb.r, colRgb.g, colRgb.b);
		glNormal3f(0.0f, 0.0f, 1.0f);
		glTexCoord2d(0.0, 0.0); glVertex3f(-1.0f, -1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
		glTexCoord2d(1.0, 0.0); glVertex3f(1.0f, -1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
		glTexCoord2d(1.0, 1.0); glVertex3f(1.0f, 1.0f, 1.0f);	// Top Right Of The Texture and Quad
		glTexCoord2d(0.0, 1.0); glVertex3f(-1.0f, 1.0f, 1.0f);	// Top Left Of The Texture and Quad
		glEnd();
		glPopName();

		if (mirrorVisibleFaces && g_bRotating && fabs(z - extend) < 0.0001)
		{
			glPushMatrix();
			glTranslated(0, 0, offsetDist * CUBE_SIZE);

			// Mirror Front Face
			glPushName((GLuint)Front);
			glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(front));
			glBegin(GL_QUADS);
			ColorRGB colRgb = ColorRGB::RGBColors[front];
			glColor3ub(colRgb.r, colRgb.g, colRgb.b);
			glNormal3f(0.0f, 0.0f, -1.0f);
			glTexCoord2d(0.0, 0.0); glVertex3f(-1.0f, -1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
			glTexCoord2d(1.0, 0.0); glVertex3f(1.0f, -1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
			glTexCoord2d(1.0, 1.0); glVertex3f(1.0f, 1.0f, 1.0f);	// Top Right Of The Texture and Quad
			glTexCoord2d(0.0, 1.0); glVertex3f(-1.0f, 1.0f, 1.0f);	// Top Left Of The Texture and Quad
			glEnd();
			glPopName();

			glTranslated(0, 0, -offsetDist * CUBE_SIZE);
			glPushMatrix();
		}

		// Back Face
		glPushName((GLuint)Back);
		glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(back));
		glBegin(GL_QUADS);
		colRgb = ColorRGB::RGBColors[back];
		glColor3ub(colRgb.r, colRgb.g, colRgb.b);
		glNormal3f(0.0f, 0.0f, -1.0f);
		glTexCoord2d(1.0, 0.0); glVertex3f(-1.0f, -1.0f, -1.0f);	// Bottom Right Of The Texture and Quad
		glTexCoord2d(1.0, 1.0); glVertex3f(-1.0f, 1.0f, -1.0f);	// Top Right Of The Texture and Quad
		glTexCoord2d(0.0, 1.0); glVertex3f(1.0f, 1.0f, -1.0f);	// Top Left Of The Texture and Quad
		glTexCoord2d(0.0, 0.0); glVertex3f(1.0f, -1.0f, -1.0f);	// Bottom Left Of The Texture and Quad
		glEnd();
		glPopName();

		if (fabs(z - -extend) < 0.0001)
		{
			glPushMatrix();
			glTranslated(0, 0, -offsetDist * CUBE_SIZE);

			// Mirror Back Face
			glPushName((GLuint)Back);
			glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(back));
			glBegin(GL_QUADS);
			colRgb = ColorRGB::RGBColors[back];
			glColor3ub(colRgb.r, colRgb.g, colRgb.b);
			glNormal3f(0.0f, 0.0f, +1.0f);
			glTexCoord2d(1.0, 0.0); glVertex3f(-1.0f, -1.0f, -1.0f);	// Bottom Right Of The Texture and Quad
			glTexCoord2d(1.0, 1.0); glVertex3f(-1.0f, 1.0f, -1.0f);	// Top Right Of The Texture and Quad
			glTexCoord2d(0.0, 1.0); glVertex3f(1.0f, 1.0f, -1.0f);	// Top Left Of The Texture and Quad
			glTexCoord2d(0.0, 0.0); glVertex3f(1.0f, -1.0f, -1.0f);	// Bottom Left Of The Texture and Quad
			glEnd();
			glPopName();

			glTranslated(0, 0, offsetDist * CUBE_SIZE);
			glPushMatrix();
		}

		// Up Face
		glPushName((GLuint)Up);
		glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(top));
		glBegin(GL_QUADS);
		colRgb = ColorRGB::RGBColors[top];
		glColor3ub(colRgb.r, colRgb.g, colRgb.b);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glTexCoord2d(0.0, 1.0); glVertex3f(-1.0f, 1.0f, -1.0f);	// Top Left Of The Texture and Quad
		glTexCoord2d(0.0, 0.0); glVertex3f(-1.0f, 1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
		glTexCoord2d(1.0, 0.0); glVertex3f(1.0f, 1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
		glTexCoord2d(1.0, 1.0); glVertex3f(1.0f, 1.0f, -1.0f);	// Top Right Of The Texture and Quad
		glEnd();
		glPopName();

		if (mirrorVisibleFaces && g_bRotating && fabs(y - extend) < 0.0001)
		{
			glPushMatrix();
			glTranslated(0, (offsetDist + 1) * CUBE_SIZE, 0);

			// Mirror Up Face
			glPushName((GLuint)Up);
			glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(top));
			glBegin(GL_QUADS);
			colRgb = ColorRGB::RGBColors[top];
			glColor3ub(colRgb.r, colRgb.g, colRgb.b);
			glNormal3f(0.0f, -1.0f, 0.0f);
			glTexCoord2d(0.0, 1.0); glVertex3f(-1.0f, 1.0f, -1.0f);	// Top Left Of The Texture and Quad
			glTexCoord2d(0.0, 0.0); glVertex3f(-1.0f, 1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
			glTexCoord2d(1.0, 0.0); glVertex3f(1.0f, 1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
			glTexCoord2d(1.0, 1.0); glVertex3f(1.0f, 1.0f, -1.0f);	// Top Right Of The Texture and Quad
			glEnd();
			glPopName();

			glTranslated(0, -(offsetDist + 1) * CUBE_SIZE, 0);
			glPushMatrix();
		}

		// Down Face
		glPushName((GLuint)Down);
		glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(bottom));
		glBegin(GL_QUADS);
		colRgb = ColorRGB::RGBColors[bottom];
		glColor3ub(colRgb.r, colRgb.g, colRgb.b);
		glNormal3f(0.0f, -1.0f, 0.0f);
		glTexCoord2d(1.0, 1.0); glVertex3f(-1.0f, -1.0f, -1.0f);	// Top Right Of The Texture and Quad
		glTexCoord2d(0.0, 1.0); glVertex3f(1.0f, -1.0f, -1.0f);	// Top Left Of The Texture and Quad
		glTexCoord2d(0.0, 0.0); glVertex3f(1.0f, -1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
		glTexCoord2d(1.0, 0.0); glVertex3f(-1.0f, -1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
		glEnd();
		glPopName();

		if (fabs(y - -extend) < 0.0001)
		{
			glPushMatrix();
			glTranslated(0, -(offsetDist + 1) * CUBE_SIZE, 0);

			// Down Face
			glPushName((GLuint)Down);
			glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(bottom));
			glBegin(GL_QUADS);
			colRgb = ColorRGB::RGBColors[bottom];
			glColor3ub(colRgb.r, colRgb.g, colRgb.b);
			glNormal3f(0.0f, +1.0f, 0.0f);
			glTexCoord2d(1.0, 1.0); glVertex3f(-1.0f, -1.0f, -1.0f);	// Top Right Of The Texture and Quad
			glTexCoord2d(0.0, 1.0); glVertex3f(1.0f, -1.0f, -1.0f);	// Top Left Of The Texture and Quad
			glTexCoord2d(0.0, 0.0); glVertex3f(1.0f, -1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
			glTexCoord2d(1.0, 0.0); glVertex3f(-1.0f, -1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
			glEnd();
			glPopName();

			glTranslated(0, (offsetDist + 1) * CUBE_SIZE, 0);
			glPushMatrix();
		}

		// Right face
		glPushName((GLuint)Right);
		glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(right));
		glBegin(GL_QUADS);
		colRgb = ColorRGB::RGBColors[right];
		glColor3ub(colRgb.r, colRgb.g, colRgb.b);
		glNormal3f(1.0f, 0.0f, 0.0f);
		glTexCoord2d(1.0, 0.0); glVertex3f(1.0f, -1.0f, -1.0f);	// Bottom Right Of The Texture and Quad
		glTexCoord2d(1.0, 1.0); glVertex3f(1.0f, 1.0f, -1.0f);	// Top Right Of The Texture and Quad
		glTexCoord2d(0.0, 1.0); glVertex3f(1.0f, 1.0f, 1.0f);	// Top Left Of The Texture and Quad
		glTexCoord2d(0.0, 0.0); glVertex3f(1.0f, -1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
		glEnd();
		glPopName();

		if (mirrorVisibleFaces && g_bRotating && fabs(x - extend) < 0.0001)
		{
			glPushMatrix();
			glTranslated(offsetDist * CUBE_SIZE, 0, 0);

			// Mirror Right face
			glPushName((GLuint)Right);
			glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(right));
			glBegin(GL_QUADS);
			colRgb = ColorRGB::RGBColors[right];
			glColor3ub(colRgb.r, colRgb.g, colRgb.b);
			glNormal3f(-1.0f, 0.0f, 0.0f);
			glTexCoord2d(1.0, 0.0); glVertex3f(1.0f, -1.0f, -1.0f);	// Bottom Right Of The Texture and Quad
			glTexCoord2d(1.0, 1.0); glVertex3f(1.0f, 1.0f, -1.0f);	// Top Right Of The Texture and Quad
			glTexCoord2d(0.0, 1.0); glVertex3f(1.0f, 1.0f, 1.0f);	// Top Left Of The Texture and Quad
			glTexCoord2d(0.0, 0.0); glVertex3f(1.0f, -1.0f, 1.0f);	// Bottom Left Of The Texture and Quad
			glEnd();
			glPopName();

			glTranslated(-offsetDist * CUBE_SIZE, 0, 0);
			glPushMatrix();
		}

		// Left Face
		glPushName((GLuint)Left);
		glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(left));
		glBegin(GL_QUADS);
		colRgb = ColorRGB::RGBColors[left];
		glColor3ub(colRgb.r, colRgb.g, colRgb.b);
		glNormal3f(-1.0f, 0.0f, 0.0f);
		glTexCoord2d(0.0, 0.0); glVertex3f(-1.0f, -1.0f, -1.0f);	// Bottom Left Of The Texture and Quad
		glTexCoord2d(1.0, 0.0); glVertex3f(-1.0f, -1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
		glTexCoord2d(1.0, 1.0); glVertex3f(-1.0f, 1.0f, 1.0f);	// Top Right Of The Texture and Quad
		glTexCoord2d(0.0, 1.0); glVertex3f(-1.0f, 1.0f, -1.0f);	// Top Left Of The Texture and Quad
		glEnd();
		glPopName();

		if (fabs(x - -extend) < 0.0001)
		{
			glPushMatrix();
			glTranslated(-offsetDist * CUBE_SIZE, 0, 0);

			// Mirror Left Face
			glPushName((GLuint)Left);
			glBindTexture(GL_TEXTURE_2D, Textures::getTextureID(left));
			glBegin(GL_QUADS);
			colRgb = ColorRGB::RGBColors[left];
			glColor3ub(colRgb.r, colRgb.g, colRgb.b);
			glNormal3f(+1.0f, 0.0f, 0.0f);
			glTexCoord2d(0.0, 0.0); glVertex3f(-1.0f, -1.0f, -1.0f);	// Bottom Left Of The Texture and Quad
			glTexCoord2d(1.0, 0.0); glVertex3f(-1.0f, -1.0f, 1.0f);	// Bottom Right Of The Texture and Quad
			glTexCoord2d(1.0, 1.0); glVertex3f(-1.0f, 1.0f, 1.0f);	// Top Right Of The Texture and Quad
			glTexCoord2d(0.0, 1.0); glVertex3f(-1.0f, 1.0f, -1.0f);	// Top Left Of The Texture and Quad
			glEnd();
			glPopName();

			glTranslated(offsetDist * CUBE_SIZE, 0, 0);
			glPushMatrix();
		}

		glPopName();
		glPopName();
		glPopName();

		glDisable(GL_TEXTURE_2D);

		glPopMatrix();
	}

	unique_ptr<RubiksCubeModel_v2::Cube> RubiksCubeModel_v2::CreateCube(double x, double y, double z, const Location& location)
	{
		Color left, right, top, bottom, front, back;

		if (x == 0)
		{
			left = Orange;
			right = Black;
		}
		else if (x == size_ - 1)
		{
			left = Black;
			right = Red;
		}
		else
		{
			left = Black;
			right = Black;
		}

		if (y == 0)
		{
			bottom = White;
			top = Black;
		}
		else if (y == size_ - 1)
		{
			bottom = Black;
			top = Yellow;
		}
		else
		{
			bottom = Black;
			top = Black;
		}

		if (z == 0)
		{
			back = Green;
			front = Black;
		}
		else if (z == size_ - 1)
		{
			back = Black;
			front = Blue;
		}
		else
		{
			back = Black;
			front = Black;
		}

		return make_unique<Cube>(top, bottom, left, right, front, back, location);
	}

	const RubiksCubeModel_v2::Cube& RubiksCubeModel_v2::GetCube(double x, double y, double z)
	{
		//if (!IsValidCube(x, y, z))
		//	RubiksCubeSolverUtils::RunTimeAssert
		
		return *cubes_[Location(x - 1, y - 1, z - 1)];
	}

	/*
	const Cube& RubiksCubeModel_v2::GetCube(Face layer0, Face layer1, Face layer2)
	{
		RubiksCubeSolverUtils::RunTimeAssert(layer0 != layer1 && layer1 != layer2);

		double x, y, z;
		double extend = (size_ - 1) / 2.0;
		for (int i = 0; i < 3; ++i)
		{
			Face layer;
			if (i == 0)
				layer = layer0;
			else if (i == 1)
				layer = layer1;
			else if (i == 2)
				layer = layer2;

			switch (layer)
			{
			case Left:
				x = -extend;
				break;
			case Right:
				x = extend;
				break;
			case Down:
				y = -extend;
				break;
			case Up:
				y = extend;
				break;
			case Back:
				z = -extend;
				break;
			case Front:
				z = extend;
				break;
			}
		}

		return *cubes_[Location(x, y, z)];		
	}
	*/

	bool RubiksCubeModel_v2::IsValidCube(int x, int y, int z)
	{
		return (x >= 0 && x < size_) &&
			(y >= 0 && y < size_) &&
			(z >= 0 && z < size_);
	}

	bool RubiksCubeModel_v2::isSolved()
	{
		return IsFaceSolved(Up) &&
			IsFaceSolved(Down) &&
			IsFaceSolved(Left) &&
			IsFaceSolved(Right) &&
			IsFaceSolved(Front) &&
			IsFaceSolved(Back);
	}

	bool RubiksCubeModel_v2::IsFaceSolved(RubiksCubeModel_v2::Face face)
	{
		double extend = (size_ - 1) / 2.0;

		if (face == Left || face == Right)
		{
			double iStart = (face == Left) ? -extend : extend;

			Color color = cubes_[Location(iStart, -extend, -extend)]->GetFaceColor(face);

			double jStart = -extend;
			for (int j = 0; j < size_; j++, ++jStart)
			{
				double kStart = -extend;
				for (int k = 0; k < size_; k++, ++kStart)
				{
					if (cubes_[Location(iStart, jStart, kStart)]->GetFaceColor(face) != color)
						return false;
				}
			}
		}
		else if (face == Up || face == Down)
		{
			double jStart = (face == Down) ? -extend : extend;

			Color color = cubes_[Location(-extend, jStart, -extend)]->GetFaceColor(face);

			double iStart = -extend;
			for (int i = 0; i < size_; i++, ++iStart)
			{
				double kStart = -extend;
				for (int k = 0; k < size_; k++, ++kStart)
				{
					if (cubes_[Location(iStart, jStart, kStart)]->GetFaceColor(face) != color)
						return false;
				}
			}
		}
		else if (face == Front || face == Back)
		{
			double kStart = (face == Back) ? -extend : extend;

			Color color = cubes_[Location(-extend, -extend, kStart)]->GetFaceColor(face);

			double iStart = -extend;
			for (int i = 0; i < size_; i++, ++iStart)
			{
				double jStart = -extend;
				for (int j = 0; j < size_; j++, ++jStart)
				{
					if (cubes_[Location(iStart, jStart, kStart)]->GetFaceColor(face) != color)
						return false;
				}
			}
		}

		return true;
	}

	unique_ptr<RubiksCubeModel> RubiksCubeModel_v2::copy()
	{
		return make_unique<RubiksCubeModel_v2>(*this);
	}

	string RubiksCubeModel_v2::getModelName()
	{
		return "RubiksCubeModel_v1";
	}

	int RubiksCubeModel_v2::getDimension()
	{
		return size_;
	}

	void RubiksCubeModel_v2::scramble(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui)
	{
		applyAlgorithm(algorithm, animate, ui);
	}

	int RubiksCubeModel_v2::applyAlgorithm(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui)
	{
		int solutionSteps = 0;
		g_bFlipRotation = false;

		for (int i = 0; i < algorithm.length(); ++i)
		{
			char face = algorithm[i];
			if (face == ' ')
				continue;

			//if (face >= 'a')
			//	face = face - 32; // Convert into Upper case char

			// Check if prime operation
			bool isPrime = false;
			int nextCharIndex = i + 1;
			if (nextCharIndex < algorithm.length() && algorithm[nextCharIndex] == '\'')
			{
				isPrime = true;
				++i;
			}
			// Check if multiple rotations
			// Check the layer index
			nextCharIndex = i + 1;
			int numRotations = 1;
			//int layerIndex = 1;
			if (nextCharIndex < algorithm.length() && '0' <= algorithm[nextCharIndex] && algorithm[nextCharIndex] <= '9')
			{
				numRotations = algorithm[i + 1] - '0';
				//layerIndex = algorithm[i + 1] - '0';
				++i;
			}

			applyStep(face, isPrime, numRotations, animate, ui);
			//applyStep(face, isPrime, layerIndex, animate, steps, ui);
			++solutionSteps;
		}

		return solutionSteps;
	}

	//const CVector3& RubiksCubeModel_v2::getRotationAxis(Groups rotationSection)
	//{
	//	switch (rotationSection)
	//	{
	//	case F:
	//		g_vRotationAxis = CVector3(0, 0, 1);
	//		break;

	//	case Z:
	//		g_vRotationAxis = CVector3(0, 0, 1);
	//		break;

	//	case B:
	//		g_vRotationAxis = CVector3(0, 0, 1);
	//		break;

	//	case L:
	//		g_vRotationAxis = CVector3(1, 0, 0);
	//		break;

	//	case X:
	//		g_vRotationAxis = CVector3(1, 0, 0);
	//		break;

	//	case R:
	//		g_vRotationAxis = CVector3(1, 0, 0);
	//		break;

	//	case U:
	//		g_vRotationAxis = CVector3(0, 1, 0);
	//		break;

	//	case Y:
	//		g_vRotationAxis = CVector3(0, 1, 0);
	//		break;

	//	case D:
	//		g_vRotationAxis = CVector3(0, 1, 0);
	//		break;

	//	default:
	//		RubiksCubeSolverUtils::RunTimeAssert(false);
	//		break;
	//	}
	//}

	void RubiksCubeModel_v2::applyStep(const char& face, bool isPrime, int numRotations, bool animate, RubiksCubeSolverGUI& ui)
	//void RubiksCubeModel_v2::applyStep(const char& face, bool isPrime, int layerIndex, bool animate /*= false*/, int steps /*= 0*/, RubiksCubeSolverGUI* ui /*= nullptr*/)
	{
		//cout << "\nApplying move: " << face;
		//if(isPrime)
		//	cout << "\'";
		//if(numRotations > 1)
		//	cout << numRotations;

		// Z-axis : Back -> Front // green  -> blue
		// X-axis : Left -> Right // orange -> red
		// Y-axis : Down -> Up    // white  -> yellow
		/*
		yellow
		Y
		|
		. --> X red
		/
		Z
		blue
		*/

		//int angle = 90;

		switch (face)
		{
		case 'F':
			g_vRotationAxis = CVector3(0, 0, 1);
			g_nRotatingSection = Front;
			g_nRotationAngle = -90;
			break;

		case 'Z':
			g_vRotationAxis = CVector3(0, 0, 1);
			g_nRotatingSection = All;
			g_nRotationAngle = 90;
			break;

		case 'B':
			g_vRotationAxis = CVector3(0, 0, 1);
			g_nRotatingSection = Back;
			g_nRotationAngle = 90;
			break;

		case 'L':
			g_vRotationAxis = CVector3(1, 0, 0);
			g_nRotatingSection = Left;
			g_nRotationAngle = 90;
			break;

		case 'X':
			g_vRotationAxis = CVector3(1, 0, 0);
			g_nRotatingSection = All;
			g_nRotationAngle = 90;
			break;

		case 'R':
			g_vRotationAxis = CVector3(1, 0, 0);
			g_nRotatingSection = Right;
			g_nRotationAngle = -90;
			break;

		case 'U':
			g_vRotationAxis = CVector3(0, 1, 0);
			g_nRotatingSection = Up;
			g_nRotationAngle = -90;
			break;

		case 'Y':
			g_vRotationAxis = CVector3(0, 1, 0);
			g_nRotatingSection = All;
			g_nRotationAngle = 90;
			break;

		case 'D':
			g_vRotationAxis = CVector3(0, 1, 0);
			g_nRotatingSection = Down;
			g_nRotationAngle = 90;
			break;

		default:
			//RubiksCubeSolverUtils::RunTimeAssert
			break;
		}

		g_nRotationAngle = g_nRotationAngle * numRotations;
		//g_nLayerIndex = layerIndex;
		if (isPrime)
			g_nRotationAngle = -g_nRotationAngle;

		if (animate)
		{
			g_bRotating = true;
			int angle = g_nRotationAngle;
			g_nRotationAngle = 0;
			int step = (angle - g_nRotationAngle) / ui.getFramesPerRotation();
			for (int i = 0; i < ui.getFramesPerRotation(); ++i)
			{
				g_nRotationAngle += step;
				ui.redrawWindow();
				Sleep(ui.getSleepTimeMilliSec());
			}

			//If after above loop, the target angle is not achieved
			if (g_nRotationAngle != angle)
			{
				g_nRotationAngle = angle;
				ui.redrawWindow();
			}
			g_bRotating = false;
		}

		//Fix cube position and Reset all parameters
		fixRubiksCubeFaces();
		if (animate)
			ui.redrawWindow();
		g_vRotationAxis = CVector3{0.0, 0.0, 0.0};
		g_nRotationAngle = 0;
		g_nRotatingSection = None;
	}


	void RubiksCubeModel_v2::fixRubiksCubeFaces()
	{
		Rotate(g_vRotationAxis, g_nRotatingSection, g_nLayerIndex, g_nRotationAngle);
	}

	void RubiksCubeModel_v2::Rotate(CVector3 rotationAxis, RubiksCubeModel_v2::Face rotatingSection, int layerIndex, double rotationAngle)
	{
		for (auto& obj : cubes_)
		{
			const Location& loc = obj.first;
			Cube& cube = *obj.second.get();

			if (cube.belongsTo(rotatingSection, layerIndex, size_))
			{
				cube.rotate(rotationAxis, rotationAngle);
				//bool noChangeInGroupAfterRotation = (rotationAxis == CVector3::XAxis && (cube.group_ == L || cube.group_ == R))
				//	|| (rotationAxis == CVector3::YAxis && (cube.group_ == L || cube.group_ == R))
				//	|| (rotationAxis == CVector3::ZAxis && (cube.group_ == L || cube.group_ == R));
				//if (rotationAxis != getRotationAxis(cube.group_))
				//cube.group_ = cube.location_.recalcGroup(size_); //Need to recalculate group always in case of any type of rotation since the cube may belong to multiple groups
			}
		}

		for (auto& obj : cubes_)
		{
			const Location& loc = obj.first;
			unique_ptr<Cube>& current = obj.second;

			//unique_ptr<Cube> current = std::move(cube);
			while (loc != current->getLocation())
			{
				unique_ptr<Cube> temp = std::move(cubes_[current->getLocation()]);
				cubes_[current->getLocation()] = std::move(current);
				current = std::move(temp);
			}
			//obj.second = std::move(current);
		}
	}

	string RubiksCubeModel_v2::generateScramblingAlgo(int length, bool includeNonStandardRotations)
	{
		char charSet[9] = { 
			'F', //Front
			'B', //Back
			'L', //Left
			'R', //Right
			'U', //Up
			'D',  //Down
			'X', //whole cube in X Axis direction X = L + R + center
			'Y', //whole cube in Y Axis direction Y = U + D + center
			'Z'  //whole cube in Z Axis direction Z = F + B + center
		};

		int numNotations = sizeof(charSet) / sizeof(char);
		int wholeCubeRotateNotations = 3; // 'X', 'Y' and 'Z'
		int numSingleLayerRotateNotations = numNotations - wholeCubeRotateNotations;
		if (!includeNonStandardRotations)
			numNotations = numSingleLayerRotateNotations;
		string retVal;
		for (int i = 0; i < length; ++i)
		{
			int index = rand() % numNotations;
			retVal += charSet[index];

			//Generate 50% prime rotations
			bool isPrime = (rand() % 2) == 1;
			if (isPrime)
				retVal += '\'';
			
			//Generate double rotations 1/3 times
			if (index < numSingleLayerRotateNotations) //Avoid X, Y and Z since layerIndex is not applicable for it
			{
				int rotations = rand() % 3; 
				if(rotations > 1)
					retVal += to_string(rotations);
			}
		}

		return retVal;
	}


	//=======================================================================================================
	// Solver
	//=======================================================================================================

	RubiksCubeModel_v2::RubiksCubeSolver::RubiksCubeSolver(RubiksCubeModel_v2& rubiksCube, bool animate, RubiksCubeSolverGUI& ui)
		: rubiksCube_(rubiksCube),
		solutionSteps_(0),
		animate_(animate),
		ui_(ui)
	{
	}

	string RubiksCubeModel_v2::RubiksCubeSolver::solve(unsigned int& solutionSteps)
	{
		solutionSteps_ = 0;
		solution_ = "";

		RubiksCubeModel_v2::RubiksCubeSolver::positionTheCube();
		RubiksCubeModel_v2::RubiksCubeSolver::buildCross();
		RubiksCubeModel_v2::RubiksCubeSolver::buildF2L();
		RubiksCubeModel_v2::RubiksCubeSolver::buildOLL();
		RubiksCubeModel_v2::RubiksCubeSolver::buildPLL();

		//verify
		RubiksCubeSolverUtils::RunTimeAssert(rubiksCube_.isSolved());

		solutionSteps = solutionSteps_;
		return solution_;
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::applyAlgorithm(const string& step)
	{
		solution_ += step;
		solutionSteps_ += rubiksCube_.applyAlgorithm(step, animate_, ui_);
	}

	bool RubiksCubeModel_v2::RubiksCubeSolver::isEdgeCube(const Cube& currentCube, const Color& first, const Color& second)
	{
		int firstCount = 0;
		int secondCount = 0;
		for (int i = 0; i < 6; ++i)
		{
			if (currentCube.GetFaceColor(Face(i)) == first)
				++firstCount;
			else if (currentCube.GetFaceColor(Face(i)) == second)
				++secondCount;
			else if (currentCube.GetFaceColor(Face(i)) != Color::Black)
				return false;
		}

		return firstCount == 1 && secondCount == 1;
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::buildCross_PlaceEdgePiece(const Color& targetColorFront, const Color& targetColorBottom)
	{
		//Cube* currentCube = nullptr;

		// Bring it from bottom later (y = 0) to top layer
		Cube currentCube = rubiksCube_.GetCube(1, 0, 2);
		Color c1 = currentCube.GetFaceColor(Face::Front);
		Color c2 = currentCube.GetFaceColor(Face::Down);

		if (c1 == targetColorFront && c2 == targetColorBottom)
		{
			//Do nothing
		}
		if (c1 == targetColorBottom && c2 == targetColorFront)
		{
			applyAlgorithm("F2");
		}
		currentCube = rubiksCube_.GetCube(2, 0, 1);
		c1 = currentCube.GetFaceColor(Face::Right);
		c2 = currentCube.GetFaceColor(Face::Down);
		if ((c1 == targetColorFront && c2 == targetColorBottom) || (c1 == targetColorBottom && c2 == targetColorFront))
		{
			applyAlgorithm("R2");
		}
		currentCube = rubiksCube_.GetCube(1, 0, 0);
		c1 = currentCube.GetFaceColor(Face::Back);
		c2 = currentCube.GetFaceColor(Face::Down);
		if ((c1 == targetColorFront && c2 == targetColorBottom) || (c1 == targetColorBottom && c2 == targetColorFront))
		{
			applyAlgorithm("B2");
		}
		currentCube = rubiksCube_.GetCube(0, 0, 1);
		c1 = currentCube.GetFaceColor(Face::Left);
		c2 = currentCube.GetFaceColor(Face::Down);
		if ((c1 == targetColorFront && c2 == targetColorBottom) || (c1 == targetColorBottom && c2 == targetColorFront))
		{
			applyAlgorithm("L2");
			//RubiksCubeAlgoExecuter::executeAlgorithm("L'F'");
		}

		// Bring it from middle later (y = 1) to top layer
		currentCube = rubiksCube_.GetCube(0, 1, 0);
		c1 = currentCube.GetFaceColor(Face::Left);
		c2 = currentCube.GetFaceColor(Face::Back);
		if ((c1 == targetColorFront && c2 == targetColorBottom) || (c1 == targetColorBottom && c2 == targetColorFront))
		{
			applyAlgorithm("LU'L'");
		}
		currentCube = rubiksCube_.GetCube(0, 1, 2);
		c1 = currentCube.GetFaceColor(Face::Left);
		c2 = currentCube.GetFaceColor(Face::Front);
		if (c1 == targetColorBottom && c2 == targetColorFront)
		{
			applyAlgorithm("F'");
		}
		else if (c1 == targetColorFront && c2 == targetColorBottom)
		{
			applyAlgorithm("F");
		}
		currentCube = rubiksCube_.GetCube(2, 1, 2);
		c1 = currentCube.GetFaceColor(Face::Front);
		c2 = currentCube.GetFaceColor(Face::Right);
		if (c1 == targetColorFront && c2 == targetColorBottom)
		{
			applyAlgorithm("F");
		}
		else if (c1 == targetColorBottom && c2 == targetColorFront)
		{
			applyAlgorithm("F'");
		}
		currentCube = rubiksCube_.GetCube(2, 1, 0);
		c1 = currentCube.GetFaceColor(Face::Back);
		c2 = currentCube.GetFaceColor(Face::Right);
		if ((c1 == targetColorFront && c2 == targetColorBottom) || (c1 == targetColorBottom && c2 == targetColorFront))
		{
			applyAlgorithm("R'UR");
		}

		// Bring it from top later (y = 2) to bottom layer at appropriate position
		currentCube = rubiksCube_.GetCube(1, 2, 0);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Back);
		if (c1 == targetColorFront && c2 == targetColorBottom)
		{
			applyAlgorithm("B'R'URBF2");
		}
		else if (c1 == targetColorBottom && c2 == targetColorFront)
		{
			applyAlgorithm("U2F2");
		}

		currentCube = rubiksCube_.GetCube(0, 2, 1);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Left);
		if (c1 == targetColorFront && c2 == targetColorBottom)
		{
			applyAlgorithm("LF'L'");
		}
		else if (c1 == targetColorBottom && c2 == targetColorFront)
		{
			applyAlgorithm("U'F2");
		}

		currentCube = rubiksCube_.GetCube(1, 2, 2);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Front);
		if (c1 == targetColorFront && c2 == targetColorBottom)
		{
			applyAlgorithm("FRUR'F2");
		}
		else if (c1 == targetColorBottom && c2 == targetColorFront)
		{
			applyAlgorithm("F2");
		}
		currentCube = rubiksCube_.GetCube(2, 2, 1);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Right);
		if (c1 == targetColorFront && c2 == targetColorBottom)
		{
			applyAlgorithm("R'FR");
		}
		else if (c1 == targetColorBottom && c2 == targetColorFront)
		{
			applyAlgorithm("UF2");
		}
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::positionTheCube()
	{
		/*
		Make sure Cube is positioned in right way

		// Z-axis : Back -> Front // green  -> blue
		// X-axis : Left -> Right // orange -> red
		// Y-axis : Down -> Up    // white  -> yellow

		yellow
		Y
		|
		. --> X red
		/
		Z
		blue
		*/

		// At the most 3 out of below 3 if clauses are executed
		Cube currentCube;
		Color c;

		// Check front face has blue center cube
		currentCube = rubiksCube_.GetCube(1, 1, 2);
		c = currentCube.GetFaceColor(Face::Front);
		if (c != Color::Blue)
		{
			if (c == Color::Green)
				applyAlgorithm("Y2");
			else if (c == Color::Orange)
				applyAlgorithm("Y'");
			else if (c == Color::Red)
				applyAlgorithm("Y");
			else if (c == Color::White)
				applyAlgorithm("X");
			else if (c == Color::Yellow)
				applyAlgorithm("X'");
		}

		//Check right face
		// Do not disturb front face, so rotate around only z axis
		currentCube = rubiksCube_.GetCube(2, 1, 1);
		c = currentCube.GetFaceColor(Face::Right);
		if (c != Color::Red)
		{
			if (c == Color::Orange)
				applyAlgorithm("Z2");
			else if (c == Color::Green)
				applyAlgorithm("Y");
			else if (c == Color::Blue)
				applyAlgorithm("Y'");
			else if (c == Color::White)
				applyAlgorithm("Z'");
			else if (c == Color::Yellow)
				applyAlgorithm("Z");
		}

		//Check top face
		currentCube = rubiksCube_.GetCube(1, 2, 1);
		c = currentCube.GetFaceColor(Face::Up);
		if (c != Color::Yellow)
		{
			if (c == Color::White)
				applyAlgorithm("X2");
			else if (c == Color::Green)
				applyAlgorithm("X'");
			else if (c == Color::Blue)
				applyAlgorithm("X");
			else if (c == Color::Orange)
				applyAlgorithm("Z");
			else if (c == Color::Red)
				applyAlgorithm("Z'");
		}
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::buildCross()
	{
		// Place blue-white at right position
		buildCross_PlaceEdgePiece(Color::Blue, Color::White);

		// Place red at right position
		applyAlgorithm("Y'");
		buildCross_PlaceEdgePiece(Color::Red, Color::White);

		// Place green at right position
		applyAlgorithm("Y'");
		buildCross_PlaceEdgePiece(Color::Green, Color::White);

		// Place orange at right position
		applyAlgorithm("Y'");
		buildCross_PlaceEdgePiece(Color::Orange, Color::White);

		applyAlgorithm("Y'");
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::buildF2L_PositionCornerPieces(const Color& targetColorFront, const Color& targetColorRight, const Color& targetColorBottom /*= Color::White*/)
	{
		Cube currentCube;
		Color c1, c2, c3;

		// Check bottom layer and bring target to top layer at (2, 2, 2)
		currentCube = rubiksCube_.GetCube(0, 0, 0);
		c1 = currentCube.GetFaceColor(Face::Back);
		c2 = currentCube.GetFaceColor(Face::Left);
		c3 = currentCube.GetFaceColor(Face::Down);
		if ((c1 == targetColorFront || c2 == targetColorFront || c3 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight || c3 == targetColorRight)
			&& (c1 == targetColorBottom || c2 == targetColorBottom || c3 == targetColorBottom)
			)
			applyAlgorithm("LU2L'");


		currentCube = rubiksCube_.GetCube(0, 0, 2);
		c1 = currentCube.GetFaceColor(Face::Front);
		c2 = currentCube.GetFaceColor(Face::Left);
		c3 = currentCube.GetFaceColor(Face::Down);
		if ((c1 == targetColorFront || c2 == targetColorFront || c3 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight || c3 == targetColorRight)
			&& (c1 == targetColorBottom || c2 == targetColorBottom || c3 == targetColorBottom)
			)
			applyAlgorithm("FU'F'U'");

		currentCube = rubiksCube_.GetCube(2, 0, 2);
		c1 = currentCube.GetFaceColor(Face::Front);
		c2 = currentCube.GetFaceColor(Face::Right);
		c3 = currentCube.GetFaceColor(Face::Down);
		if (c1 == targetColorFront || c2 == targetColorRight || c3 == targetColorBottom)
		{
			//do nothing
		}
		else if ((c1 == targetColorFront || c2 == targetColorFront || c3 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight || c3 == targetColorRight)
			&& (c1 == targetColorBottom || c2 == targetColorBottom || c3 == targetColorBottom)
			)
			applyAlgorithm("F'UF");


		currentCube = rubiksCube_.GetCube(2, 0, 0);
		c1 = currentCube.GetFaceColor(Face::Back);
		c2 = currentCube.GetFaceColor(Face::Right);
		c3 = currentCube.GetFaceColor(Face::Down);
		if ((c1 == targetColorFront || c2 == targetColorFront || c3 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight || c3 == targetColorRight)
			&& (c1 == targetColorBottom || c2 == targetColorBottom || c3 == targetColorBottom)
			)
			applyAlgorithm("R'URU");

		// Check top layer and bring target to (2, 2, 2)
		currentCube = rubiksCube_.GetCube(0, 2, 0);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Left);
		c3 = currentCube.GetFaceColor(Face::Back);
		if ((c1 == targetColorFront || c2 == targetColorFront || c3 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight || c3 == targetColorRight)
			&& (c1 == targetColorBottom || c2 == targetColorBottom || c3 == targetColorBottom)
			)
			applyAlgorithm("U2");

		currentCube = rubiksCube_.GetCube(0, 2, 2);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Left);
		c3 = currentCube.GetFaceColor(Face::Front);
		if ((c1 == targetColorFront || c2 == targetColorFront || c3 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight || c3 == targetColorRight)
			&& (c1 == targetColorBottom || c2 == targetColorBottom || c3 == targetColorBottom)
			)
			applyAlgorithm("U'");

		currentCube = rubiksCube_.GetCube(2, 2, 0);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Back);
		c3 = currentCube.GetFaceColor(Face::Right);
		if ((c1 == targetColorFront || c2 == targetColorFront || c3 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight || c3 == targetColorRight)
			&& (c1 == targetColorBottom || c2 == targetColorBottom || c3 == targetColorBottom)
			)
			applyAlgorithm("U");

		// Target is now in top layer at (2, 2, 2)
		currentCube = rubiksCube_.GetCube(2, 2, 2);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Front);
		c3 = currentCube.GetFaceColor(Face::Right);

		if (c1 == targetColorFront && c2 == targetColorBottom && c3 == targetColorRight)
			applyAlgorithm("F'U'F");
		else if (c1 == targetColorRight && c2 == targetColorFront && c3 == targetColorBottom)
			applyAlgorithm("RUR'");
		else if (c1 == targetColorBottom && c2 == targetColorRight && c3 == targetColorFront)
			applyAlgorithm("RUUR'U'RUR'");
		else
		{
			//RubiksCubeSolverUtils::RunTimeAssert
			int i = 0;
			++i;
		}
	}

	bool RubiksCubeModel_v2::RubiksCubeSolver::buildF2L_PositionEdgePieces(const Color& targetColorFront, const Color& targetColorRight)
	{
		Cube currentCube;
		Color c1, c2;
		bool retVal = true;
		string algo1("URU'R'U'F'UF");
		string algo2("U'F'UFURU'R'");

		//Check if aleady at position
		currentCube = rubiksCube_.GetCube(2, 1, 2);
		c1 = currentCube.GetFaceColor(Face::Front);
		c2 = currentCube.GetFaceColor(Face::Right);
		if (c1 == targetColorFront && c2 == targetColorRight)
			return true;
		else if (c1 == targetColorRight && c2 == targetColorFront) // If piece is stuck at right position but in wrong orientation
			applyAlgorithm(algo1);

		// Check top layer
		currentCube = rubiksCube_.GetCube(1, 2, 0);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Back);
		if ((c1 == targetColorFront || c2 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight))
			applyAlgorithm("U");

		currentCube = rubiksCube_.GetCube(0, 2, 1);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Left);
		if ((c1 == targetColorFront || c2 == targetColorFront)
			&& (c1 == targetColorRight || c2 == targetColorRight))
			applyAlgorithm("U'");


		currentCube = rubiksCube_.GetCube(1, 2, 2);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Front);
		if (c1 == targetColorRight && c2 == targetColorFront)
			applyAlgorithm(algo1);

		else if (c1 == targetColorFront && c2 == targetColorRight)
			applyAlgorithm("U'" + algo2);
		else
			retVal = false;

		if (retVal)
			return retVal;

		retVal = true;
		currentCube = rubiksCube_.GetCube(2, 2, 1);
		c1 = currentCube.GetFaceColor(Face::Up);
		c2 = currentCube.GetFaceColor(Face::Right);
		if (c1 == targetColorFront && c2 == targetColorRight)
			applyAlgorithm(algo2);
		else if (c1 == targetColorRight && c2 == targetColorFront)
			applyAlgorithm("U" + algo1);
		else
			retVal = false;

		//If we fail, check if any edge piece is stuck in second layer
		if (!retVal)
		{
			currentCube = rubiksCube_.GetCube(2, 1, 2);
			c1 = currentCube.GetFaceColor(Face::Front);
			c2 = currentCube.GetFaceColor(Face::Right);
			if (c1 != Color::Yellow && c2 != Color::Yellow)
				applyAlgorithm(algo1);
		}

		return retVal;
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::buildF2L()
	{
		//position corner pieces
		buildF2L_PositionCornerPieces(Color::Blue, Color::Red, Color::White);

		applyAlgorithm("Y'");
		buildF2L_PositionCornerPieces(Color::Red, Color::Green, Color::White);

		applyAlgorithm("Y'");
		buildF2L_PositionCornerPieces(Color::Green, Color::Orange, Color::White);

		applyAlgorithm("Y'");
		buildF2L_PositionCornerPieces(Color::Orange, Color::Blue, Color::White);

		applyAlgorithm("Y'");

		//position edge pieces
		int numIterations = 0;
		int done = 0;
		while (done != 15)
		{
			++numIterations;

			if (buildF2L_PositionEdgePieces(Color::Blue, Color::Red))
				done |= 1;

			applyAlgorithm("Y'");
			if (buildF2L_PositionEdgePieces(Color::Red, Color::Green))
				done |= 2;

			applyAlgorithm("Y'");
			if (buildF2L_PositionEdgePieces(Color::Green, Color::Orange))
				done |= 4;

			applyAlgorithm("Y'");
			if (buildF2L_PositionEdgePieces(Color::Orange, Color::Blue))
				done |= 8;

			applyAlgorithm("Y'");
		}
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::buildOLL()
	{
		// Step 1: build yellow cross on top face

		while (true)
		{
			Cube currentCube;
			Color c, c1, c2, c3, c4;
			string algo("FRUR'U'F'");

			/*
			Top Face
			*   c1  *
			c4  c   c2
			*   c3  *

			*/

			//Check if aleady at position
			currentCube = rubiksCube_.GetCube(1, 2, 1);
			c = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(1, 2, 0);
			c1 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(2, 2, 1);
			c2 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(1, 2, 2);
			c3 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(0, 2, 1);
			c4 = currentCube.GetFaceColor(Face::Up);
			if (c1 == Color::Yellow
				&& c2 == Color::Yellow
				&& c3 == Color::Yellow
				&& c4 == Color::Yellow) // We are sure that c is Color::Yellow
				break;

			if (c1 == Color::Yellow && c3 == Color::Yellow)
				applyAlgorithm("U");

			if (c2 == Color::Yellow && c4 == Color::Yellow)
			{
				applyAlgorithm(algo);
				continue;
			}

			if (c1 == Color::Yellow && c2 == Color::Yellow)
				applyAlgorithm("U'");
			else if (c2 == Color::Yellow && c3 == Color::Yellow)
				applyAlgorithm("U2");
			else if (c3 == Color::Yellow && c4 == Color::Yellow)
				applyAlgorithm("U");

			if (c1 == Color::Yellow && c4 == Color::Yellow)
			{
				applyAlgorithm(algo + algo);
				continue;
			}

			// Do the sequence once if none of above was executed
			applyAlgorithm(algo);
		}

		// Step 2: get all yellow on top face
		while (true)
		{
			Cube currentCube;
			Color c1, c2, c3, c4, c5, c6, c7, c8, c9;
			Color s1, s2, s3, s4, s5, s6, s7, s8;
			string algo("RUR'URUUR'");

			/*
			Top Face
			   s2      s3
			s1 c1  c2  c3 s4
			   c4  c5  c6
			s8 c7  c8  c9 s5
			   s7      s6
			*/

			//Check if aleady at position
			currentCube = rubiksCube_.GetCube(0, 2, 0);
			c1 = currentCube.GetFaceColor(Face::Up);
			s1 = currentCube.GetFaceColor(Face::Left);
			s2 = currentCube.GetFaceColor(Face::Back);
			currentCube = rubiksCube_.GetCube(1, 2, 0);
			c2 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(2, 2, 0);
			c3 = currentCube.GetFaceColor(Face::Up);
			s3 = currentCube.GetFaceColor(Face::Back);
			s4 = currentCube.GetFaceColor(Face::Right);
			currentCube = rubiksCube_.GetCube(0, 2, 1);
			c4 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(1, 2, 1);
			c5 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(2, 2, 1);
			c6 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(0, 2, 2);
			c7 = currentCube.GetFaceColor(Face::Up);
			s7 = currentCube.GetFaceColor(Face::Front);
			s8 = currentCube.GetFaceColor(Face::Left);
			currentCube = rubiksCube_.GetCube(1, 2, 2);
			c8 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(2, 2, 2);
			c9 = currentCube.GetFaceColor(Face::Up);
			s5 = currentCube.GetFaceColor(Face::Right);
			s6 = currentCube.GetFaceColor(Face::Front);
			if (c1 == Color::Yellow
				&& c2 == Color::Yellow
				&& c3 == Color::Yellow
				&& c4 == Color::Yellow
				&& c5 == Color::Yellow
				&& c6 == Color::Yellow
				&& c7 == Color::Yellow
				&& c8 == Color::Yellow
				&& c9 == Color::Yellow)
				break;

			int numYellowCorners = 0;
			if (c1 == Color::Yellow)
				++numYellowCorners;
			if (c3 == Color::Yellow)
				++numYellowCorners;
			if (c7 == Color::Yellow)
				++numYellowCorners;
			if (c9 == Color::Yellow)
				++numYellowCorners;

			int* n1 = nullptr;
			switch (numYellowCorners)
			{
			case 0:
				if (s6 == Color::Yellow)
					applyAlgorithm("U");
				else if (s4 == Color::Yellow)
					applyAlgorithm("U2");
				else if (s2 == Color::Yellow)
					applyAlgorithm("U'");
				else if (s8 == Color::Yellow)
					applyAlgorithm(""); // do nothing
				else
				{
					int* n = nullptr;
					*n = 10; //we should not come here
				}
				break;
			case 1: // After this you will go to case 1 again or you would be done
				if (c9 == Color::Yellow)
					applyAlgorithm("U");
				else if (c3 == Color::Yellow)
					applyAlgorithm("U2");
				else if (c1 == Color::Yellow)
					applyAlgorithm("U'");
				else if (c7 == Color::Yellow)
					applyAlgorithm(""); // do nothing
				else
				{
					int* n = nullptr;
					*n = 10; //we should not come here
				}
				break;
			case 2: // After this you will go to case 1 
				if (s5 == Color::Yellow)
					applyAlgorithm("U");
				else if (s3 == Color::Yellow)
					applyAlgorithm("U2");
				else if (s1 == Color::Yellow)
					applyAlgorithm("U'");
				else if (s7 == Color::Yellow)
					applyAlgorithm(""); // do nothing
				else
				{
					int* n = nullptr;
					*n = 10; //we should not come here
				}
				break;
			case 3:
				*n1 = 10; //we should not come here
				break;
			default:
				*n1 = 10; //we should not come here
				break;
			}

			// Do the sequence once and continue
			applyAlgorithm(algo);
		}
	}

	void RubiksCubeModel_v2::RubiksCubeSolver::buildPLL()
	{
		//Step 1
		while (true)
		{
			Cube currentCube;
			//Color c1, c2, c3, c4, c5, c6, c7, c8, c9;
			Color s1, s2, s3, s4, s5, s6, s7, s8;
			string algo("RB'RFFR'BRFFRR");

			/*
			Top Face
			s2  o1  s3
			s1 c1  c2  c3 s4
			o4	c4  c5  c6 o2
			s8 c7  c8  c9 s5
			s7  o3  s6
			*/

			//Check if aleady at position
			currentCube = rubiksCube_.GetCube(0, 2, 0);
			//c1 = currentCube.GetFaceColor(Face::Up);
			s1 = currentCube.GetFaceColor(Face::Left);
			s2 = currentCube.GetFaceColor(Face::Back);
			//currentCube = Scene::getInstance().g_cCube.GetCube(1, 2, 0);
			//c2 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(2, 2, 0);
			//c3 = currentCube.GetFaceColor(Face::Up);
			s3 = currentCube.GetFaceColor(Face::Back);
			s4 = currentCube.GetFaceColor(Face::Right);
			//currentCube = Scene::getInstance().g_cCube.GetCube(0, 2, 1);
			//c4 = currentCube.GetFaceColor(Face::Up);
			//currentCube = Scene::getInstance().g_cCube.GetCube(1, 2, 1);
			//c5 = currentCube.GetFaceColor(Face::Up);
			//currentCube = Scene::getInstance().g_cCube.GetCube(2, 2, 1);
			//c6 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(0, 2, 2);
			//c7 = currentCube.GetFaceColor(Face::Up);
			s7 = currentCube.GetFaceColor(Face::Front);
			s8 = currentCube.GetFaceColor(Face::Left);
			//currentCube = Scene::getInstance().g_cCube.GetCube(1, 2, 2);
			//c8 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(2, 2, 2);
			//c9 = currentCube.GetFaceColor(Face::Up);
			s5 = currentCube.GetFaceColor(Face::Right);
			s6 = currentCube.GetFaceColor(Face::Front);
			if (s2 == s3
				&& s4 == s5
				&& s6 == s7
				&& s8 == s1)
			{
				//Match centers and corners before proceeding
				//Get centers
				Color o1, o2, o3, o4;
				currentCube = rubiksCube_.GetCube(1, 1, 0);
				o1 = currentCube.GetFaceColor(Face::Back);
				currentCube = rubiksCube_.GetCube(2, 1, 1);
				o2 = currentCube.GetFaceColor(Face::Right);
				currentCube = rubiksCube_.GetCube(1, 1, 2);
				o3 = currentCube.GetFaceColor(Face::Front);
				currentCube = rubiksCube_.GetCube(0, 1, 1);
				o4 = currentCube.GetFaceColor(Face::Left);

				if (o1 == s4)
					applyAlgorithm("U'");
				else if (o1 == s6)
					applyAlgorithm("U2");
				else if (o1 == s8)
					applyAlgorithm("U");

				break;
			}

			//Rotate the complete cube to set the face "having two corner piece color same" as front face
			if (s4 == s5)
				applyAlgorithm("Y'");
			else if (s2 == s3)
				applyAlgorithm("Y2");
			else if (s1 == s8)
				applyAlgorithm("Y");

			applyAlgorithm(algo);
		}

		//Step 2
		while (true)
		{
			Cube currentCube;
			//Color c1, c2, c3, c4, c5, c6, c7, c8, c9;
			//Color s1, s2, s3, s4, s5, s6, s7, s8;
			Color e1, e2, e3, e4;
			Color s4, s6, s8;
			string algo("RU'RURURU'R'U'RR");

			//Get centers
			/*
			o1
			e1
			c1  c2  c3
			o4  e4	c4  c5  c6 e2 o2
			c7  c8  c9
			e3
			o3
			*/
			Color o1, o2, o3, o4;
			currentCube = rubiksCube_.GetCube(1, 1, 0);
			o1 = currentCube.GetFaceColor(Face::Back);
			currentCube = rubiksCube_.GetCube(2, 1, 1);
			o2 = currentCube.GetFaceColor(Face::Right);
			currentCube = rubiksCube_.GetCube(1, 1, 2);
			o3 = currentCube.GetFaceColor(Face::Front);
			currentCube = rubiksCube_.GetCube(0, 1, 1);
			o4 = currentCube.GetFaceColor(Face::Left);

			/*
			Top Face
			o1
			e1
			c1  c2  c3
			o4  e4	c4  c5  c6 e2 o2
			c7  c8  c9
			e3
			o3
			*/

			//Check if aleady at position
			//currentCube = Scene::getInstance().g_cCube.GetCube(0, 2, 0);
			//c1 = currentCube.GetFaceColor(Face::Up);
			//s1 = currentCube.GetFaceColor(Face::Left);
			//s2 = currentCube.GetFaceColor(Face::Back);
			currentCube = rubiksCube_.GetCube(1, 2, 0);
			//c2 = currentCube.GetFaceColor(Face::Up);
			e1 = currentCube.GetFaceColor(Face::Back);
			//currentCube = Scene::getInstance().g_cCube.GetCube(2, 2, 0);
			//c3 = currentCube.GetFaceColor(Face::Up);
			//s3 = currentCube.GetFaceColor(Face::Back);
			//s4 = currentCube.GetFaceColor(Face::Right);
			currentCube = rubiksCube_.GetCube(0, 2, 1);
			//c4 = currentCube.GetFaceColor(Face::Up);
			e4 = currentCube.GetFaceColor(Face::Left);
			//currentCube = Scene::getInstance().g_cCube.GetCube(1, 2, 1);
			//c5 = currentCube.GetFaceColor(Face::Up);
			currentCube = rubiksCube_.GetCube(2, 2, 1);
			//c6 = currentCube.GetFaceColor(Face::Up);
			e2 = currentCube.GetFaceColor(Face::Right);
			//currentCube = Scene::getInstance().g_cCube.GetCube(0, 2, 2);
			//c7 = currentCube.GetFaceColor(Face::Up);
			//s7 = currentCube.GetFaceColor(Face::Front);
			//s8 = currentCube.GetFaceColor(Face::Left);
			currentCube = rubiksCube_.GetCube(1, 2, 2);
			//c8 = currentCube.GetFaceColor(Face::Up);
			e3 = currentCube.GetFaceColor(Face::Front);
			//currentCube = Scene::getInstance().g_cCube.GetCube(2, 2, 2);
			//c9 = currentCube.GetFaceColor(Face::Up);
			//s5 = currentCube.GetFaceColor(Face::Right);
			//s6 = currentCube.GetFaceColor(Face::Front);

			//Match centers with corners, they may be misaligned after few iterations of last algo
			//if(o1 == s4)
			//	RubiksCubeAlgoExecuter::executeAlgorithm("U'");
			//else if(o1 == s6)
			//	RubiksCubeAlgoExecuter::executeAlgorithm("U2");
			//else if(o1 == s8)
			//	RubiksCubeAlgoExecuter::executeAlgorithm("U");

			if (e1 == o1
				&& e2 == o2
				&& e3 == o3
				&& e4 == o4)
			{
				break;
			}

			if (e2 == o2)
				applyAlgorithm("Y");
			else if (e3 == o3)
				applyAlgorithm("Y2");
			else if (e4 == o4)
				applyAlgorithm("Y'");

			applyAlgorithm(algo);
		}
	}

}

