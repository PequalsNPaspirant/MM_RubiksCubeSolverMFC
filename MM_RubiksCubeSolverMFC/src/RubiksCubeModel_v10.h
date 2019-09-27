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
#include <unordered_set>
#include <chrono>
using namespace std;

#include <gl\gl.h>
#include <gl\glu.h>
//#include <glm/glm.hpp>
#include "Vector3.h"
#include "RubiksCubeModel.h"

namespace mm {

	class RubiksCubeSolverGUI;

	class RubiksCubeModel_v10 : public RubiksCubeModel
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

		enum class cubeType
		{
			rubiksCube,
			mirrorCube
		};

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

			const vector<vector<double>>& getRotationMatrix(const CVector3& rotationAxis, double rotationAngle)
			{
				static vector<vector<double>> matrix(4, vector<double>(4, 0.0));
				double angle = rotationAngle * (PI / 180.0); //Angle should be in radians
				// Initialize rotation matrix
				for (int Row = 0; Row < 4; Row++)
					for (int Column = 0; Column < 4; Column++)
						if (Row == Column)
							matrix[Row][Column] = 1.0;
						else
							matrix[Row][Column] = 0.0;

				if (rotationAxis == CVector3::XAxis)
				{
					matrix[1][1] = cos(angle);
					matrix[1][2] = sin(angle);
					matrix[2][1] = -sin(angle);
					matrix[2][2] = cos(angle);
				}
				else if (rotationAxis == CVector3::YAxis)
				{
					matrix[0][0] = cos(angle);
					matrix[0][2] = -sin(angle);
					matrix[2][0] = sin(angle);
					matrix[2][2] = cos(angle);
				}
				else if (rotationAxis == CVector3::ZAxis)
				{
					matrix[0][0] = cos(angle);
					matrix[0][1] = sin(angle);
					matrix[1][0] = -sin(angle);
					matrix[1][1] = cos(angle);
				}

				for (int Row = 0; Row < 4; Row++)
					for (int Column = 0; Column < 4; Column++)
						if (fabs(matrix[Row][Column]) < 0.000001)
							matrix[Row][Column] = 0.0;

				return matrix;

				using RotationMatrix = vector<vector<double>>;
				static vector<vector<RotationMatrix>> rotationMatrixSet(3, vector<RotationMatrix>(3, RotationMatrix(4, vector<double>(4, 0.0))));

				static bool firstTime = true;
				if (firstTime)
				{
					firstTime = false;					
					for (int i = 0; i < 3; ++i)
					{
						for (int j = 0; j < 3; ++j)
						{
							double angle = (j + 1) * 90 * (PI / 180.0); //Angle should be in radians

							RotationMatrix& matrix = rotationMatrixSet[i][j];

							// Initialize rotation matrix
							for (int Row = 0; Row < 4; Row++)
								for (int Column = 0; Column < 4; Column++)
									if (Row == Column)
										matrix[Row][Column] = 1.0;
									else
										matrix[Row][Column] = 0.0;

							if (i == 0)
							{
								matrix[1][1] = cos(angle);
								matrix[1][2] = sin(angle);
								matrix[2][1] = -sin(angle);
								matrix[2][2] = cos(angle);
							}
							else if (i == 1)
							{
								matrix[0][0] = cos(angle);
								matrix[0][2] = -sin(angle);
								matrix[2][0] = sin(angle);
								matrix[2][2] = cos(angle);
							}
							else if (i == 2)
							{
								matrix[0][0] = cos(angle);
								matrix[0][1] = sin(angle);
								matrix[1][0] = -sin(angle);
								matrix[1][1] = cos(angle);
							}

							for (int Row = 0; Row < 4; Row++)
								for (int Column = 0; Column < 4; Column++)
									if (fabs(matrix[Row][Column]) < 0.000001)
										matrix[Row][Column] = 0.0;
						}
					}
				}

				int i = 0;
				if (rotationAxis == CVector3::YAxis)
					i = 1;
				else if (rotationAxis == CVector3::ZAxis)
					i = 2;

				if (rotationAngle < 0)
					rotationAngle += 360;

				int j = round(rotationAngle / 90.0);
				--j;

				return rotationMatrixSet[i][j];
			}

			void rotate(const CVector3& rotationAxis, double rotationAngle)
			{
				vector< vector<double>> geomVecMatrix(1, vector<double>(4, 1.0));
				geomVecMatrix[0][0] = x_;
				geomVecMatrix[0][1] = y_;
				geomVecMatrix[0][2] = z_;
				vector< vector<double>> result(1, vector<double>(4, 0.0));

				const vector<vector<double>>& rotationMatrix = getRotationMatrix(rotationAxis, rotationAngle);

				for (int i = 0; i < 1; i++)
					for (int j = 0; j < 4; j++)
						for (int k = 0; k < 4; k++) // OR use Matrix2.m_row. Both are equal.
							result[i][j] += geomVecMatrix[i][k] * rotationMatrix[k][j];

				x_ = result[0][0];
				y_ = result[0][1];
				z_ = result[0][2];
			}

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
			{}
			Cube(Cube&&) = default;
			Cube& operator=(const Cube&) = default;
			Cube& operator=(Cube&&) = default;
			Cube(Color cTop, Color cBottom, Color cLeft,
				Color cRight, Color cFront, Color cBack, const Location& location, const Location& cubeCenter);
			~Cube();
			Color GetFaceColor(Face eFace) const;
			void TiltUp();
			void TiltDown();
			void TurnLeft();
			void TurnRight();
			void TiltLeft();
			void TiltRight();

			const Location& getLocation() const { return location_; }
			void initializeMatrix();
			void rotateCubeCenter(CVector3 rotationAxis, double rotationAngle);
			void rotateLocation(CVector3 rotationAxis, double rotationAngle);
			void rotateThickness(CVector3 rotationAxis, double rotationAngle);
			void fixRubiksCubeFaces(CVector3 rotationAxis, double rotationAngle);
			void setThickness(double x, double y, double z) { xt_ = x; yt_ = y; zt_ = z; }
			void getThickness(double& x, double& y, double& z) const { x = xt_; y = yt_; z = zt_; }

			GLfloat matrixf_[16];
			double matrix_[16];
			//glGetFloatv(GL_MODELVIEW_MATRIX, matrixf);
			//glLoadMatrixf(matrixProjection.get());

		private:
			static const int FACE_COUNT = 6;
			vector<Color> faces_;
			Location cubeCenter_; //Location of actual cube center
			Location location_; //Location of cube in NxNxN Rubik's Cube matrix of equal sized sub-cubes
			double xt_, yt_, zt_; //thicknesses in each direction
		};

	public:
		RubiksCubeModel_v10(int size, double xt, double yt, double zt);
		~RubiksCubeModel_v10();
		RubiksCubeModel_v10(const RubiksCubeModel_v10& copy);

		void setAnimate(bool animate) override { animate_ = animate; }
		bool getAnimate() { return animate_; }
		void ResetCube(bool animate, RubiksCubeSolverGUI* ui) override {}
		void ResetCube(int size, double xt, double yt, double zt, bool animate, RubiksCubeSolverGUI* ui) override;
		void scramble(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui) override {}
		bool scramble(const string& algorithm, bool animate, RubiksCubeSolverGUI& ui, string& invalidStep) override;
		int applyAlgorithm(const string& algorithm);
		string generateScramblingAlgo(int length, bool includeNonStandardRotations) override;
		string solve(unsigned int& solutionSteps, unsigned long long& duration, bool animate, RubiksCubeSolverGUI& ui) override;
		void render() override;
		void renderIndividualCube(const Cube& pCube, const Location& location);
		bool isSolved() override;
		bool IsFaceSolved(Face face);
		void getUpdatedStats(unsigned int& size, unsigned int& scramblingSteps, string& scramblingAlgo, unsigned int& solutionSteps, string& solution, unsigned long long& duration, string& status) override;

		bool activateRubiksCube() override;
		bool activateMirrorCube() override;

		unique_ptr<RubiksCubeModel> copy() override;
		string getModelName() override;
		int getDimension() override;

		bool pauseAnimation(bool pause) override;

		Cube& GetCube(Face layer1, int layerIndex1, Face layer2, int layerIndex2, Face layer3, int layerIndex3);
		void fixRubiksCubeFaces(CVector3 rotationAxis, Face rotatingSection, int layerIndexFrom, int layerIndexTo, double rotationAngle);
		int getSize() const { return size_; }

	private:
		bool extractSteps(const string& algo);
		string invalidStep_;
		struct AlgoStep
		{
			string thisStep;
			char face;
			unsigned int layerIndexFrom; 
			unsigned int layerIndexTo;
			bool isPrime;
			unsigned int numRotations;
		};
		vector<AlgoStep> algoSteps_;
		bool pauseAnimation_{ false };

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
				//return key.x_ * 9 + key.y_ * 3 + key.z_;
				std::size_t seed = 0;
				hash_combine(seed, key.x_);
				hash_combine(seed, key.y_);
				hash_combine(seed, key.z_);
				return seed;
			}
		};

		bool IsValidCube(int x, int y, int z);
		void applyStep(const char& face, int layerIndexFrom, int layerIndexTo, bool isPrime, int numRotations);
		void fixRubiksCubeFaces();

		//vector< vector< vector<Cube> > > cubes_; // Total elements = size_ * size_ * size_
		unordered_map<Location, unique_ptr<Cube>, Hasher> cubes_; // Total elements = size_ * size_ * size_ - ( (size_-2) * (size_-2) * (size_-2) )
		//vector<unique_ptr<Cube>> cubes_;
		//vector< vector<Cube*> > layerF_; //Front layer //Total elements = size_ * size_
		//vector< vector<Cube*> > layerB_; //Back layer
		//vector< vector<Cube*> > layerL_; //Left layer
		//vector< vector<Cube*> > layerR_; //Right layer
		//vector< vector<Cube*> > layerU_; //Upper layer
		//vector< vector<Cube*> > layerD_; //Down layer

		//using Iterator = vector<unique_ptr<Cube>>::iterator;
		//using Edge = vector<Iterator>; // max number of elelements = size_
		//using Loop = Edge[4];
		//using Layer = vector<Loop>;
		//vector<Layer> frontLayers{ size_ }; //size = size_
		//vector<Layer> leftLayers{ size_ };
		//vector<Layer> upLayers{ size_ };

		int size_{ 3 };
		bool g_bRotating;
		bool g_bFlipRotation;
		CVector3 g_vRotationAxis;
		Face g_nRotatingSection;
		int g_nLayerIndexFrom;
		int g_nLayerIndexTo;
		double g_nRotationAngle;
		const int subCubeSize_{ 2 };
		int extend_{ subCubeSize_ * (size_ - 1) / 2 };
		cubeType cubeType_{ cubeType::rubiksCube };

		int scramblingSteps_;
		string scramblingAlgo_;
		bool isScrambling_{ false };
		bool isSolving_{ false };
		int solutionSteps_;
		string solution_;
		string solutionFileFullPath_;
		ofstream* pSolutionFile_;
		unsigned long long duration_;
		using HRClock = std::chrono::high_resolution_clock;
		HRClock::time_point startTime_;
		
		static const double scale_;
		bool animate_;
		//RubiksCubeSolverGUI& ui_;
		RubiksCubeSolverGUI* pUi_;

		static const vector<string> statusStrings;

		enum class status
		{
			scrambled = 0,
			face_reduction_to_3x3x3,
			edge_reduction_to_3x3x3,
			cross,
			F2L,
			OLL,
			PLL,
			solved,

			eMaxStatus
		};

		status status_{ status::solved };

	public:
		int getScramblingSteps() { return scramblingSteps_; }
		const string& getScramblingAlgo() { return scramblingAlgo_; }
		int getSolutionSteps() { return solutionSteps_; }
		const string& getSolution() { return solution_; }

		//=======================================================================================================
		// Solver
		//=======================================================================================================

	public:
		class RubiksCubeSolver_NxNxN
		{
		public:
			//RubiksCubeSolver_NxNxN(RubiksCubeModel_v10& rubiksCube, bool animate, RubiksCubeSolverGUI& ui);
			RubiksCubeSolver_NxNxN(RubiksCubeModel_v10& rubiksCube);
			string solve(unsigned int& solutionSteps);

		private:
			void reduceTo3x3x3();
			/**/void fixCenterCubes_singleFace(Color targetColor);
			/**//**/bool fixCenterCubes_moveTargetCubeFromUpToRightFace(Face fromFace, const string& preMove,
				int targetLineIndexFromLeft, int targetIndexFromUp, Color targetColor, 
				int columnFromLeftToAvoid, int centerColumnToAvoid);
			/**//**/bool fixCenterCubes_moveTargetCubeFromFrontToRightFace(Face fromFace, const string& preMove,
				int targetLineIndexFromLeft, int targetIndexFromUp, Color targetColor,
				int columnFromLeftToAvoid, int centerColumnToAvoid);
			/**//**/bool fixCenterCubes_moveTargetCubeToFrontFace(Face fromFace, const string& frontFacePreMove,
				int targetLineIndexFromLeft, int targetIndexFromUp, Color targetColor);
			/**/void fixCenterCubes_twoFaces(Color targetColor1, Color targetColor2);
			/**/void fixCenterCubes_twoFaces2(Color targetColor1, Color targetColor2);
			/**//**/bool fixCenterCubes_twoFaces_CheckIfLineExist(Face faceFront, Face faceLeft, int targetLineIndexFromLeft, Face faceUp, Color targetColorRightFace, const string& algo);
			/**//**/bool fixCenterCubes_twoFaces_moveTargetCubeFromRightToFrontFace(int targetLineIndexFromLeft, int targetIndexFromUp, Color targetColor, int startRowFromTopToAvoid, int endRowFromTopToAvoid);
			/**//**/bool fixCenterCubes_twoFaces_moveTargetCubeFromFrontToRightFace(int targetLineIndexFromLeft, int targetIndexFromUp, Color targetColor);

			/**/void fixEdgeCubes(Color targetColorUp, Color targetColorFront);
			/**//**/bool fixEdgeCubes_checkIfEdgeIsAlreadySolved(Color targetColorUp, Color targetColorFront);
			/**//**/bool fixEdgeCubes_ensureUpRightEdgeUnsolved();
			/**//**/bool fixEdgeCubes_bringToUpBackEdge(int targetIndexLeft, Color targetColorUp, Color targetColorFront);
			/**//**//**/bool fixEdgeCubes_bringToUpBackEdge_searchEdge(int targetIndex, Color targetColorUp, Color targetColorFront,
				Face face1, Face face2, Face face3, const string& algo, bool sameOrientation);
			/**//**/bool fixEdgeCubes_lastTwoEdges(Color targetColorUp1, Color targetColorFront1, Color targetColorUp2, Color targetColorBack2);

			void positionTheCube();
			void buildCross();
			/**/void buildCross_PlaceEdgePiece(const Color& targetColorFront, const Color& targetColorBottom);
			void buildF2L();
			/**/bool buildF2L_PositionEdgeColumns(const Color& targetColorFront, const Color& targetColorRight);
			void buildOLL();
			void buildPLL();

		private:
			void applyAlgorithm(const string& step);
			//bool isEdgeCube(const Cube& currentCube, const Color& first, const Color& second);

			RubiksCubeModel_v10& rubiksCube_;
			unordered_set<unsigned int> solvedEdges_;
		};
	};

}