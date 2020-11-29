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
#include "RubiksCubeModel_v1.h"
#include "RubiksCubeSolverGUI.h"
#include "RubiksCubeSolverUtils.h"

namespace mm {

	//Factory function definition
	unique_ptr<RubiksCubeModel> createRubiksCubeModel_v1(int size)
	{
		return make_unique<RubiksCubeModel_v1>(size);
	}

	//Create a global object, so that its constructor is called before main and the factory map is initialized before main
	static RegisterRubiksCubeFactoryFunction object("RubiksCubeModel_v1", createRubiksCubeModel_v1);

	//==================== RubiksCubeModel_v1::Cube =========================

	const int RubiksCubeModel_v1::Cube::FACE_COUNT /*= 6*/;
	const double RubiksCubeModel_v1::CUBE_SIZE = 2.0;

	const RubiksCubeModel_v1::ColorRGB RubiksCubeModel_v1::ColorRGB::RGBColors[7] = {
		ColorRGB{ 255, 255, 0 },
		ColorRGB{ 255, 0, 0 },
		ColorRGB{ 0, 0, 255 },
		ColorRGB{ 0, 255, 0 },
		ColorRGB{ 255, 165, 0 },
		ColorRGB{ 255, 255, 255 },
		ColorRGB{ 0, 0, 0 }
	};

	RubiksCubeModel_v1::Cube::Cube(Color cTop, Color cBottom, Color cLeft, Color cRight, Color cFront, Color cBack)
		: faces_(FACE_COUNT)
	{
		faces_[Up] = cTop;
		faces_[Down] = cBottom;
		faces_[Left] = cLeft;
		faces_[Right] = cRight;
		faces_[Front] = cFront;
		faces_[Back] = cBack;
	}

	RubiksCubeModel_v1::Cube::~Cube(void)
	{
	}

	Color RubiksCubeModel_v1::Cube::GetFaceColor(Face eFace) const
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

	// Aound X axis
	void RubiksCubeModel_v1::Cube::TiltUp()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Front];
		faces_[Front] = faces_[Down];
		faces_[Down] = faces_[Back];
		faces_[Back] = temp1;
	}

	// Aound X axis
	void RubiksCubeModel_v1::Cube::TiltDown()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Back];
		faces_[Back] = faces_[Down];
		faces_[Down] = faces_[Front];
		faces_[Front] = temp1;
	}

	//Around Y axis
	void RubiksCubeModel_v1::Cube::TurnLeft()
	{
		Color temp1 = faces_[Front];
		faces_[Front] = faces_[Right];
		faces_[Right] = faces_[Back];
		faces_[Back] = faces_[Left];
		faces_[Left] = temp1;
	}

	//Around Y axis
	void RubiksCubeModel_v1::Cube::TurnRight()
	{
		Color temp1 = faces_[Front];
		faces_[Front] = faces_[Left];
		faces_[Left] = faces_[Back];
		faces_[Back] = faces_[Right];
		faces_[Right] = temp1;
	}

	//Around Z axis
	void RubiksCubeModel_v1::Cube::TiltLeft()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Right];
		faces_[Right] = faces_[Down];
		faces_[Down] = faces_[Left];
		faces_[Left] = temp1;
	}

	//Around Z axis
	void RubiksCubeModel_v1::Cube::TiltRight()
	{
		Color temp1 = faces_[Up];
		faces_[Up] = faces_[Left];
		faces_[Left] = faces_[Down];
		faces_[Down] = faces_[Right];
		faces_[Right] = temp1;
	}

	//==================== RubiksCubeModel_v1 =========================

	RubiksCubeModel_v1::RubiksCubeModel_v1(int size)
		: cubes_(vector< vector< vector<Cube> > > (size, vector< vector<Cube> >(size, vector<Cube>(size)) ) ),
		size_(size),
		g_bRotating(false),
		g_bFlipRotation(false),
		g_vRotationAxis(0, 0, 0),
		g_nRotatingSection(-1),
		g_nRotationAngle(0)
	{
		ResetCube(false, nullptr);
	}

	RubiksCubeModel_v1::RubiksCubeModel_v1(const RubiksCubeModel_v1& copy)
		: cubes_(copy.cubes_),
		size_(copy.size_),
		g_bRotating(copy.g_bRotating),
		g_bFlipRotation(copy.g_bFlipRotation),
		g_vRotationAxis(copy.g_vRotationAxis),
		g_nRotatingSection(copy.g_nRotatingSection),
		g_nRotationAngle(copy.g_nRotationAngle)
	{
	}

	RubiksCubeModel_v1::~RubiksCubeModel_v1()
	{
	}

	void RubiksCubeModel_v1::ResetCube(bool animate, RubiksCubeSolverGUI* ui)
	{
		g_bRotating = false;
		g_bFlipRotation = false;
		g_vRotationAxis = CVector3(0, 0, 0);
		g_nRotatingSection = -1;
		g_nRotationAngle = 0;

		for (int i = 0; i < size_; i++)
		{
			for (int j = 0; j < size_; j++)
			{
				for (int k = 0; k < size_; k++)
				{
					cubes_[i][j][k] = CreateCube(i, j, k);
				}
			}
		}

		if (animate)
		{
			ui->redrawWindow();
		}
	}

	string RubiksCubeModel_v1::solve(unsigned int& solutionSteps, unsigned long long& duration, bool animate, RubiksCubeSolverGUI& ui)
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

	void RubiksCubeModel_v1::render()
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

		for (int i = 0; i < getSize(); i++)
		{
			for (int j = 0; j < getSize(); j++)
			{
				for (int k = 0; k < getSize(); k++)
				{
					glPushMatrix();

					if (g_bRotating)
					{
						if (g_vRotationAxis.x && i == g_nRotatingSection ||
							g_vRotationAxis.y && j == g_nRotatingSection ||
							g_vRotationAxis.z && k == g_nRotatingSection)
						{
							int angle = g_bFlipRotation ? -g_nRotationAngle : g_nRotationAngle;
							glRotated(angle, g_vRotationAxis.x, g_vRotationAxis.y, g_vRotationAxis.z);
						}
					}

					renderIndividualCube(GetCube(i, j, k), i, j, k);

					glPopMatrix();
				}
			}
		}
	}

	void RubiksCubeModel_v1::renderIndividualCube(const Cube& pCube, int x, int y, int z)
	{
		glPushName(x);
		glPushName(y);
		glPushName(z);

		// scale to -1 to +1
		x--;
		y--;
		z--;

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

		glPopName();
		glPopName();
		glPopName();

		glDisable(GL_TEXTURE_2D);

		glPopMatrix();
	}

	RubiksCubeModel_v1::Cube RubiksCubeModel_v1::CreateCube(int x, int y, int z)
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

		return Cube(top, bottom, left, right, front, back);
	}

	const RubiksCubeModel_v1::Cube& RubiksCubeModel_v1::GetCube(int x, int y, int z)
	{
		//if (!IsValidCube(x, y, z))
		//	RubiksCubeSolverUtils::RunTimeAssert
		
		return cubes_[x][y][z];
	}

	bool RubiksCubeModel_v1::IsValidCube(int x, int y, int z)
	{
		return (x >= 0 && x < size_) &&
			(y >= 0 && y < size_) &&
			(z >= 0 && z < size_);
	}

	void RubiksCubeModel_v1::Rotate(int section, int turns)
	{
		if (section >= 0 && section < 3)
		{
			for (int i = 0; i < size_; i++)
			{
				for (int k = 0; k < size_; k++)
				{
					for (int l = 0; l < turns; l++)
					{
						cubes_[i][section][k].TurnRight();
					}
				}
			}

			for (int i = 0; i < turns; i++)
			{
				Cube temp1 = cubes_[0][section][0];
				Cube temp2 = cubes_[0][section][2];

				cubes_[0][section][2] = temp1;

				temp1 = cubes_[2][section][2];
				cubes_[2][section][2] = temp2;

				temp2 = cubes_[2][section][0];
				cubes_[2][section][0] = temp1;

				cubes_[0][section][0] = temp2;

				temp1 = cubes_[1][section][0];
				temp2 = cubes_[0][section][1];

				cubes_[0][section][1] = temp1;
				temp1 = cubes_[1][section][2];

				cubes_[1][section][2] = temp2;
				temp2 = cubes_[2][section][1];

				cubes_[2][section][1] = temp1;
				cubes_[1][section][0] = temp2;
			}
		}
	}

	void RubiksCubeModel_v1::Tilt(int section, int turns)
	{
		if (section >= 0 && section < 3)
		{
			for (int j = 0; j < size_; j++)
			{
				for (int k = 0; k < size_; k++)
				{
					for (int l = 0; l < turns; l++)
					{
						cubes_[section][j][k].TiltDown();
					}
				}
			}

			for (int i = 0; i < turns; i++)
			{
				Cube temp1 = cubes_[section][0][0];
				Cube temp2 = cubes_[section][2][0];

				cubes_[section][2][0] = temp1;

				temp1 = cubes_[section][2][2];
				cubes_[section][2][2] = temp2;

				temp2 = cubes_[section][0][2];
				cubes_[section][0][2] = temp1;

				cubes_[section][0][0] = temp2;

				temp1 = cubes_[section][1][0];
				temp2 = cubes_[section][2][1];

				cubes_[section][2][1] = temp1;
				temp1 = cubes_[section][1][2];

				cubes_[section][1][2] = temp2;
				temp2 = cubes_[section][0][1];

				cubes_[section][0][1] = temp1;
				cubes_[section][1][0] = temp2;
			}
		}
	}

	void RubiksCubeModel_v1::Turn(int section, int turns)
	{
		if (section >= 0 && section < 3)
		{
			// rotate each cube
			for (int i = 0; i < size_; i++)
			{
				for (int j = 0; j < size_; j++)
				{
					for (int l = 0; l < turns; l++)
					{
						cubes_[i][j][section].TiltLeft();
					}
				}
			}

			for (int i = 0; i < turns; i++)
			{
				Cube temp1 = cubes_[0][0][section];
				Cube temp2 = cubes_[2][0][section];

				cubes_[2][0][section] = temp1;

				temp1 = cubes_[2][2][section];
				cubes_[2][2][section] = temp2;

				temp2 = cubes_[0][2][section];
				cubes_[0][2][section] = temp1;

				cubes_[0][0][section] = temp2;

				temp1 = cubes_[0][1][section];
				temp2 = cubes_[1][0][section];

				cubes_[1][0][section] = temp1;
				temp1 = cubes_[2][1][section];

				cubes_[2][1][section] = temp2;
				temp2 = cubes_[1][2][section];

				cubes_[1][2][section] = temp1;
				cubes_[0][1][section] = temp2;
			}
		}
	}

	void RubiksCubeModel_v1::Randomize()
	{
		int count = 0;
		bool done = false;
		srand((unsigned)time(NULL));

		while (!done)
		{
			int turns, section, axis;
			turns = (int)((double)rand() / (RAND_MAX + 1) * (4));
			section = (int)((double)rand() / (RAND_MAX + 1) * (size_));
			axis = (int)((double)rand() / (RAND_MAX + 1) * (3));

			switch (axis)
			{
			case 0:
				this->Rotate(section, turns);
				break;
			case 1:
				this->Tilt(section, turns);
			case 2:
				this->Turn(section, turns);
			}

			count++;

			if (count >= 20)
			{
				int diff = count - 20;
				int probability = (int)((double)rand() / (RAND_MAX + 1) * (100 - diff) + diff);

				if (probability >= 75)
					done = true;
			}
		}
	}

	bool RubiksCubeModel_v1::isSolved()
	{
		return IsFaceSolved(Up) &&
			IsFaceSolved(Down) &&
			IsFaceSolved(Left) &&
			IsFaceSolved(Right) &&
			IsFaceSolved(Front) &&
			IsFaceSolved(Back);
	}

	bool RubiksCubeModel_v1::IsFaceSolved(Face face)
	{
		if (face == Up || face == Down)
		{
			int j = (face == Up) ? 2 : 0;

			Color color = cubes_[0][j][0].GetFaceColor(face);

			for (int i = 0; i < size_; i++)
			{
				for (int k = 0; k < size_; k++)
				{
					if (cubes_[i][j][k].GetFaceColor(face) != color)
						return false;
				}
			}
		}

		else if (face == Left || face == Right)
		{
			int i = (face == Left) ? 0 : 2;

			Color color = cubes_[i][0][0].GetFaceColor(face);

			for (int j = 0; j < size_; j++)
			{
				for (int k = 0; k < size_; k++)
				{
					if (cubes_[i][j][k].GetFaceColor(face) != color)
						return false;
				}
			}
		}

		else if (face == Front || face == Back)
		{
			int k = (face == Front) ? 2 : 0;

			Color color = cubes_[0][0][k].GetFaceColor(face);

			for (int i = 0; i < size_; i++)
			{
				for (int j = 0; j < size_; j++)
				{
					if (cubes_[i][j][k].GetFaceColor(face) != color)
						return false;
				}
			}
		}

		return true;
	}

	unique_ptr<RubiksCubeModel> RubiksCubeModel_v1::copy()
	{
		return make_unique<RubiksCubeModel_v1>(*this);
	}

	string RubiksCubeModel_v1::getModelName()
	{
		return "RubiksCubeModel_v1";
	}

	int RubiksCubeModel_v1::getDimension()
	{
		return size_;
	}

	void RubiksCubeModel_v1::scramble(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui)
	{
		applyAlgorithm(algorithm, animate, ui);
	}

	int RubiksCubeModel_v1::applyAlgorithm(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui)
	{
		int solutionSteps = 0;
		g_bFlipRotation = false;

		for (int i = 0; i < algorithm.length(); ++i)
		{
			char face = algorithm[i];
			if (face == ' ')
				continue;

			if (face >= 'a')
				face = face - 32; // Convert into Upper case char

								  // Check if prime operation
			bool isPrime = false;
			int nextCharIndex = i + 1;
			if (nextCharIndex < algorithm.length() && algorithm[nextCharIndex] == '\'')
			{
				isPrime = true;
				++i;
			}
			// Check if multiple rotations
			nextCharIndex = i + 1;
			int numRotations = 1;
			if (nextCharIndex < algorithm.length() && '0' <= algorithm[nextCharIndex] && algorithm[nextCharIndex] <= '9')
			{
				numRotations = algorithm[i + 1] - '0';
				++i;
			}

			applyStep(face, isPrime, numRotations, animate, ui);
			++solutionSteps;
		}

		return solutionSteps;
	}

	void RubiksCubeModel_v1::applyStep(const char& face, bool isPrime, int numRotations, bool animate, RubiksCubeSolverGUI& ui)
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
			g_nRotatingSection = 2;
			g_nRotationAngle = -90;
			break;

		case 'Z':
			g_vRotationAxis = CVector3(0, 0, 1);
			g_nRotatingSection = 1;
			g_nRotationAngle = 90;
			break;

		case 'B':
			g_vRotationAxis = CVector3(0, 0, 1);
			g_nRotatingSection = 0;
			g_nRotationAngle = 90;
			break;

		case 'L':
			g_vRotationAxis = CVector3(1, 0, 0);
			g_nRotatingSection = 0;
			g_nRotationAngle = 90;
			break;

		case 'X':
			g_vRotationAxis = CVector3(1, 0, 0);
			g_nRotatingSection = 1;
			g_nRotationAngle = 90;
			break;

		case 'R':
			g_vRotationAxis = CVector3(1, 0, 0);
			g_nRotatingSection = 2;
			g_nRotationAngle = -90;
			break;

		case 'U':
			g_vRotationAxis = CVector3(0, 1, 0);
			g_nRotatingSection = 2;
			g_nRotationAngle = -90;
			break;

		case 'Y':
			g_vRotationAxis = CVector3(0, 1, 0);
			g_nRotatingSection = 1;
			g_nRotationAngle = 90;
			break;

		case 'D':
			g_vRotationAxis = CVector3(0, 1, 0);
			g_nRotatingSection = 0;
			g_nRotationAngle = 90;
			break;

		default:
			//RubiksCubeSolverUtils::RunTimeAssert
			break;
		}

		g_nRotationAngle = g_nRotationAngle * numRotations;
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
		g_nRotationAngle = 0;
		g_nRotatingSection = -1;
	}


	void RubiksCubeModel_v1::fixRubiksCubeFaces()
	{
		int turns = 0;
		if (g_nRotationAngle == 0)
			turns = 0;
		else if (g_nRotationAngle == 90)
			turns = 1;
		else if (g_nRotationAngle == 180 || g_nRotationAngle == -180)
			turns = 2;
		else if (g_nRotationAngle == -90)
			turns = 3;

		if (g_vRotationAxis.x)
			Tilt(g_nRotatingSection, turns);
		else if (g_vRotationAxis.y)
			Rotate(g_nRotatingSection, turns);
		else if (g_vRotationAxis.z)
			Turn(g_nRotatingSection, turns);
	}

	string RubiksCubeModel_v1::generateScramblingAlgo(int length, bool includeNonStandardRotations)
	{
		char charSet[9] = { 
			'F', //Front
			'Z', //Center layer between F and B
			'B', //Back
			'L', //Left
			'X', //Center layer between L and R
			'R', //Right
			'U', //Up
			'Y', //Center layer between U and D
			'D'  //Down
		};
		string retVal;
		for (int i = 0; i < length; ++i)
		{
			int index = rand() % 9;
			retVal += charSet[index];

			if (rand() % 2 == 0)
				retVal += '\'';
			else if (rand() % 10 == 0)
				retVal += '2';
		}

		return retVal;
	}


	//=======================================================================================================
	// Solver
	//=======================================================================================================

	RubiksCubeModel_v1::RubiksCubeSolver::RubiksCubeSolver(RubiksCubeModel_v1& rubiksCube, bool animate, RubiksCubeSolverGUI& ui)
		: rubiksCube_(rubiksCube),
		solutionSteps_(0),
		animate_(animate),
		ui_(ui)
	{
	}

	string RubiksCubeModel_v1::RubiksCubeSolver::solve(unsigned int& solutionSteps)
	{
		solutionSteps_ = 0;
		solution_ = "";

		RubiksCubeModel_v1::RubiksCubeSolver::positionTheCube();
		RubiksCubeModel_v1::RubiksCubeSolver::buildCross();
		RubiksCubeModel_v1::RubiksCubeSolver::buildF2L();
		RubiksCubeModel_v1::RubiksCubeSolver::buildOLL();
		RubiksCubeModel_v1::RubiksCubeSolver::buildPLL();

		//verify
		RubiksCubeSolverUtils::RunTimeAssert(rubiksCube_.isSolved());

		solutionSteps = solutionSteps_;
		return solution_;
	}

	void RubiksCubeModel_v1::RubiksCubeSolver::applyAlgorithm(const string& step)
	{
		solution_ += step;
		solutionSteps_ += rubiksCube_.applyAlgorithm(step, animate_, ui_);
	}

	bool RubiksCubeModel_v1::RubiksCubeSolver::isEdgeCube(const Cube& currentCube, const Color& first, const Color& second)
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

	void RubiksCubeModel_v1::RubiksCubeSolver::buildCross_PlaceEdgePiece(const Color& targetColorFront, const Color& targetColorBottom)
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

	void RubiksCubeModel_v1::RubiksCubeSolver::positionTheCube()
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

	void RubiksCubeModel_v1::RubiksCubeSolver::buildCross()
	{
		// Place blue-white at right position
		buildCross_PlaceEdgePiece(Color::Blue, Color::White);

		// Place red at right position
		applyAlgorithm("UY'D'");
		buildCross_PlaceEdgePiece(Color::Red, Color::White);

		// Place green at right position
		applyAlgorithm("UY'D'");
		buildCross_PlaceEdgePiece(Color::Green, Color::White);

		// Place orange at right position
		applyAlgorithm("UY'D'");
		buildCross_PlaceEdgePiece(Color::Orange, Color::White);

		applyAlgorithm("UY'D'");
	}

	void RubiksCubeModel_v1::RubiksCubeSolver::buildF2L_PositionCornerPieces(const Color& targetColorFront, const Color& targetColorRight, const Color& targetColorBottom /*= Color::White*/)
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

	bool RubiksCubeModel_v1::RubiksCubeSolver::buildF2L_PositionEdgePieces(const Color& targetColorFront, const Color& targetColorRight)
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

	void RubiksCubeModel_v1::RubiksCubeSolver::buildF2L()
	{
		//position corner pieces
		buildF2L_PositionCornerPieces(Color::Blue, Color::Red, Color::White);

		applyAlgorithm("UY'D'");
		buildF2L_PositionCornerPieces(Color::Red, Color::Green, Color::White);

		applyAlgorithm("UY'D'");
		buildF2L_PositionCornerPieces(Color::Green, Color::Orange, Color::White);

		applyAlgorithm("UY'D'");
		buildF2L_PositionCornerPieces(Color::Orange, Color::Blue, Color::White);

		applyAlgorithm("UY'D'");

		//position edge pieces
		int numIterations = 0;
		int done = 0;
		while (done != 15)
		{
			++numIterations;

			if (buildF2L_PositionEdgePieces(Color::Blue, Color::Red))
				done |= 1;

			applyAlgorithm("UY'D'");
			if (buildF2L_PositionEdgePieces(Color::Red, Color::Green))
				done |= 2;

			applyAlgorithm("UY'D'");
			if (buildF2L_PositionEdgePieces(Color::Green, Color::Orange))
				done |= 4;

			applyAlgorithm("UY'D'");
			if (buildF2L_PositionEdgePieces(Color::Orange, Color::Blue))
				done |= 8;

			applyAlgorithm("UY'D'");
		}
	}

	void RubiksCubeModel_v1::RubiksCubeSolver::buildOLL()
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

	void RubiksCubeModel_v1::RubiksCubeSolver::buildPLL()
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
				applyAlgorithm("UY'D'");
			else if (s2 == s3)
				applyAlgorithm("U2Y2D2");
			else if (s1 == s8)
				applyAlgorithm("U'YD");

			applyAlgorithm(algo);
		}

		//Step 2
		while (true)
		{
			Cube currentCube;
			//Color c1, c2, c3, c4, c5, c6, c7, c8, c9;
			//Color s1, s2, s3, s4, s5, s6, s7, s8;
			Color e1, e2, e3, e4;
			//Color s4, s6, s8;
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
				applyAlgorithm("U'YD");
			else if (e3 == o3)
				applyAlgorithm("U2Y2D2");
			else if (e4 == o4)
				applyAlgorithm("UY'D'");

			applyAlgorithm(algo);
		}
	}

}

