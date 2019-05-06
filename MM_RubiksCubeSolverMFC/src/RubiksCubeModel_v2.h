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
#include "Vector3.h"
#include "RubiksCubeModel.h"

namespace mm {

	class RubiksCubeSolverGUI;

	class RubiksCubeModel_v2 : public RubiksCubeModel
	{
	public:
		typedef unsigned char byte;

		//TODO: Define enum Axis
		/*
		enum Axis
		{
		XAxis = 0,
		YAxis = 1,
		ZAxis = 2
		};
		*/

		//enum Groups
		//{
		//	None = 0,

		//	L = 1,
		//	R = 2,

		//	D = 4,
		//	U = 8,

		//	B = 16,
		//	F = 32,

		//	All = 64 - 1 //Stores All flags from 1 to 32
		//};

		enum Face
		{
			Up = 0,
			Down = 1,
			Left = 2,
			Right = 3,
			Front = 4,
			Back = 5,

			All = 6,
			None = 7
		};

		struct ColorRGB
		{
			ColorRGB(byte red, byte green, byte blue)
			{
				r = red;
				g = green;
				b = blue;
			}

			ColorRGB()
			{
				r = g = b = 0;
			}

			byte r;
			byte g;
			byte b;

			bool operator==(const ColorRGB& c)
			{
				return (r == c.r) && (g == c.g) && (b == c.b);
			}

			bool operator!=(const ColorRGB& c)
			{
				return !((r == c.r) && (g == c.g) && (b == c.b));
			}

			static const ColorRGB RGBColors[7];
		};

		class Utility
		{
		public:
			//static Groups getGroup(char layer, int index)
			//{
			//}

			//static const CVector3& getAxis(char layer, int index)
			//{
			//	static vector<CVector3*> data{};
			//}
		};

		struct Location
		{
			Location()
				: x_(0.0), y_(0.0), z_(0.0)
			{}
			Location(double x, double y, double z)
				: x_(x), y_(y), z_(z)
			{}
			Location(const Location&) = default;
			Location(Location&&) = default;
			Location& operator=(const Location&) = default;
			Location& operator=(Location&&) = default;

			bool operator==(const Location& rhs) const
			{
				return fabs(x_ - rhs.x_) < 0.0001
					&& fabs(y_ - rhs.y_) < 0.0001
					&& fabs(z_ - rhs.z_) < 0.0001;
			}

			bool operator!=(const Location& rhs) const
			{
				return !(*this == rhs);
			}

			void rotate(CVector3 rotationAxis, double rotationAngle)
			{
				static vector<vector<double>> rotationMatrix(4, vector<double>(4));

				// Initialize rotation matrix
				for (int Row = 0; Row < 4; Row++)
					for (int Column = 0; Column < 4; Column++)
						if (Row == Column)
							rotationMatrix[Row][Column] = 1.0;
						else
							rotationMatrix[Row][Column] = 0.0;

				if (rotationAxis == CVector3::XAxis)
				{
					rotationMatrix[1][1] = cos(rotationAngle);
					rotationMatrix[1][2] = sin(rotationAngle);
					rotationMatrix[2][1] = -sin(rotationAngle);
					rotationMatrix[2][2] = cos(rotationAngle);
				}
				else if (rotationAxis == CVector3::YAxis)
				{
					rotationMatrix[0][0] = cos(rotationAngle);
					rotationMatrix[0][2] = -sin(rotationAngle);
					rotationMatrix[2][0] = sin(rotationAngle);
					rotationMatrix[2][2] = cos(rotationAngle);
				}
				else if (rotationAxis == CVector3::ZAxis)
				{
					rotationMatrix[0][0] = cos(rotationAngle);
					rotationMatrix[0][1] = sin(rotationAngle);
					rotationMatrix[1][0] = -sin(rotationAngle);
					rotationMatrix[1][1] = cos(rotationAngle);
				}

				for (int Row = 0; Row < 4; Row++)
					for (int Column = 0; Column < 4; Column++)
						if (fabs(rotationMatrix[Row][Column]) < 0.000001)
							rotationMatrix[Row][Column] = 0.0;

				vector< vector<double>> geomVecMatrix(1, vector<double>(4, 1.0));
				geomVecMatrix[0][0] = x_;
				geomVecMatrix[0][1] = y_;
				geomVecMatrix[0][2] = z_;
				vector< vector<double>> result(1, vector<double>(4, 0.0));

				for (int i = 0; i < 1; i++)
					for (int j = 0; j < 4; j++)
						for (int k = 0; k < 4; k++) // OR use Matrix2.m_row. Both are equal.
							result[i][j] += geomVecMatrix[i][k] * rotationMatrix[k][j];

				x_ = result[0][0];
				y_ = result[0][1];
				z_ = result[0][2];
			}

			//int recalcGroup(int size)
			//{
			//	double extend = (size - 1) / 2.0;
			//	int group = 0;
			//	if (fabs(x_ - (-extend)) < 0.0001)
			//		group |= Groups::L;
			//	else if (fabs(x_ - extend) < 0.0001)
			//		group |= Groups::R;

			//	if (fabs(y_ - (-extend)) < 0.0001)
			//		group |= Groups::D;
			//	else if (fabs(y_ - extend) < 0.0001)
			//		group |= Groups::U;

			//	if (fabs(z_ - (-extend)) < 0.0001)
			//		group |= Groups::B;
			//	else if (fabs(z_ - extend) < 0.0001)
			//		group |= Groups::F;

			//	return group;
			//}

			double x_;
			double y_;
			double z_;
		};

		class Cube
		{
		public:
			Cube() {}
			Cube(const Cube& copy)
				: faces_(copy.faces_),
				location_(copy.location_)
				//group_(copy.group_)
			{}
			Cube(Cube&&) = default;
			Cube& operator=(const Cube&) = default;
			Cube& operator=(Cube&&) = default;
			Cube(Color cTop, Color cBottom, Color cLeft,
				Color cRight, Color cFront, Color cBack, const Location& location);
			~Cube();
			Color GetFaceColor(Face eFace) const;
			void TiltUp();
			void TiltDown();
			void TurnLeft();
			void TurnRight();
			void TiltLeft();
			void TiltRight();

			const Location& getLocation() const { return location_; }
			void rotate(CVector3 rotationAxis, double rotationAngle);
			bool belongsTo(Face rotatingSection, int layerIndex, int size) const;

		private:
			static const int FACE_COUNT = 6;
			vector<Color> faces_;
			Location location_;
			//int group_;
		};

	public:
		RubiksCubeModel_v2(int size);
		~RubiksCubeModel_v2();
		RubiksCubeModel_v2(const RubiksCubeModel_v2& copy);

		void ResetCube(bool animate, RubiksCubeSolverGUI* ui) override;
		void scramble(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui) override;
		int applyAlgorithm(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui);
		string generateScramblingAlgo(int length, bool includeNonStandardRotations) override;
		string solve(unsigned int& solutionSteps, unsigned long long& duration, bool animate, RubiksCubeSolverGUI& ui) override;
		void render() override;
		void renderIndividualCube(const Cube& pCube, const Location& location);
		bool isSolved() override;
		bool IsFaceSolved(Face face);

		unique_ptr<RubiksCubeModel> copy() override;
		string getModelName() override;
		int getDimension() override;

		const Cube& GetCube(double x, double y, double z);
		//const Cube& GetCube(Face layer0, Face layer1, Face layer2);
		void Rotate(CVector3 rotationAxis, Face rotatingSection, int layerIndex, double rotationAngle);
		int getSize() { return size_; }

	private:
		//const CVector3& getRotationAxis(Groups rotationSection); //TODO: add this to localise group <--> Axis relation

		template <class T>
		static inline void hash_combine(std::size_t& seed, const T& v)  // Similar to boost::hash_combine<T>
		{
			/*
			The magic number below is the reciprocal of the golden ratio.
			Reference:
			https://stackoverflow.com/questions/4948780/magic-number-in-boosthash-combine
			http://burtleburtle.net/bob/hash/doobs.html
			phi = (1 + sqrt(5)) / 2
			2^32 / phi = 0x9e3779b9
			*/
			std::hash<T> hasher{};
			seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}

		struct Hasher
		{
			size_t operator()(const Location& key) const noexcept
			{
				std::size_t seed = 0;
				hash_combine(seed, key.x_);
				hash_combine(seed, key.y_);
				hash_combine(seed, key.z_);
				return seed;
			}
		};

		bool IsValidCube(int x, int y, int z);
		unique_ptr<Cube> CreateCube(double x, double y, double z, const Location& location);
		void applyStep(const char& face, bool isPrime, int numRotations, bool animate, RubiksCubeSolverGUI& ui);
		void fixRubiksCubeFaces();

		//vector< vector< vector<Cube> > > cubes_; // Total elements = size_ * size_ * size_
		unordered_map<Location, unique_ptr<Cube>, Hasher> cubes_; // Total elements = size_ * size_ * size_ - ( (size_-2) * (size_-2) * (size_-2) )
		//vector< vector<Cube*> > layerF_; //Front layer //Total elements = size_ * size_
		//vector< vector<Cube*> > layerB_; //Back layer
		//vector< vector<Cube*> > layerL_; //Left layer
		//vector< vector<Cube*> > layerR_; //Right layer
		//vector< vector<Cube*> > layerU_; //Upper layer
		//vector< vector<Cube*> > layerD_; //Down layer
		int size_;
		bool g_bRotating;
		bool g_bFlipRotation;
		CVector3 g_vRotationAxis;
		Face g_nRotatingSection;
		int g_nLayerIndex;
		double g_nRotationAngle;
		
		static const double CUBE_SIZE;


		//=======================================================================================================
		// Solver
		//=======================================================================================================

	public:

		class RubiksCubeSolver
		{
		public:
			RubiksCubeSolver(RubiksCubeModel_v2& rubiksCube, bool animate, RubiksCubeSolverGUI& ui);
			string solve(unsigned int& solutionSteps);

		private:
			void positionTheCube();
			void buildCross();
			void buildF2L();
			void buildOLL();
			void buildPLL();

			void buildCross_PlaceEdgePiece(const Color& targetColorFront, const Color& targetColorBottom);
			void buildF2L_PositionCornerPieces(const Color& targetColorFront, const Color& targetColorRight, const Color& targetColorBottom = Color::White);
			bool buildF2L_PositionEdgePieces(const Color& targetColorFront, const Color& targetColorRight);

		private:
			void applyAlgorithm(const string& step);
			bool isEdgeCube(const Cube& currentCube, const Color& first, const Color& second);

			RubiksCubeModel_v2& rubiksCube_;
			string solution_;
			int solutionSteps_;

			bool animate_;
			RubiksCubeSolverGUI& ui_;
		};
	};

}