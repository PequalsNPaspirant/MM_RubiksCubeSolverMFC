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
#include <unordered_map>
using namespace std;

namespace mm {

	class RubiksCubeSolverGUI;

	struct testInfoAggregate
	{
		testInfoAggregate()
		{}

		testInfoAggregate(const string& modelName, unsigned int size)
			: modelName_(modelName), size_(size), numTestCases_(0), nsAggregateDuration_(0)
		{}
		string modelName_;
		unsigned int size_;
		unsigned int numTestCases_;
		unsigned long long nsAggregateDuration_;
	};

	struct testInfo
	{
		testInfo(const string& modelName, unsigned int size, const string& scrambleAlgo, const string& idealSolution)
			: modelName_(modelName),
			size_(size),
			numSteps_(0),
			nsDuration_(0),
			scrambleAlgo_(scrambleAlgo),
			idealSolution_(idealSolution)
		{}

		testInfo(const string& modelName, unsigned int size, unsigned int numSteps, unsigned long long nsDuration, const string& scrambleAlgo,
			const string& idealSolution, const string& actualSolution)
			: modelName_(modelName),
			size_(size),
			numSteps_(numSteps),
			nsDuration_(nsDuration),
			scrambleAlgo_(scrambleAlgo),
			idealSolution_(idealSolution),
			actualSolution_(actualSolution)
		{}

		string modelName_;
		unsigned int size_;
		unsigned int numSteps_;
		unsigned long long nsDuration_;
		string scrambleAlgo_;
		string idealSolution_;
		string actualSolution_;
	};

	struct AlgoPairs
	{
		string scramble;
		string solution;
	};

	class RubiksCubeSolverTest
	{
	public:
		RubiksCubeSolverTest(RubiksCubeSolverGUI& refUI)
			: refUI_(refUI)
		{}

		bool testRubiksCube(bool animate);
		static string getCurrentLocalTimeInNanoSeconds2();

	private:
		void executeTest(testInfo& info, bool animate, unsigned int testNum);

		void writeResultsToCSVFile(ofstream& testResultsFile, const vector<testInfoAggregate>& testInfoAggregateSet);
		void writeResultsToCSVFile(ofstream& testResultsFile, const vector<testInfo>& testInfoSet);

		RubiksCubeSolverGUI& refUI_;

		static const int minSize;
		static const int maxSize;
		static const int incrementSize;
		static vector<int> scramblingAlgoLengths;
		static const int numAlgoOfEachLength;
		static vector<string> genericModels;

		void generateHardcodedTestCases(unordered_map<int, vector<AlgoPairs>>& scrambleAlgosHardcoded, vector<testInfoAggregate>& testInfoAggregateSetHardcoded);
		void generateBasicTestCases(vector<AlgoPairs>& scrambleAlgosBasic, vector<testInfoAggregate>& testInfoAggregateSetBasic);
		void generateGenericTestCases(vector<AlgoPairs>& scrambleAlgosGeneric, vector<testInfoAggregate>& testInfoAggregateSetGeneric);

		void generateTestCases(unordered_map<int, vector<AlgoPairs>>& scrambleAlgos, vector<testInfoAggregate>& testInfoAggregateSet);

		bool executeAllTests(unordered_map<int, vector<AlgoPairs>>& scrambleAlgos, vector<testInfoAggregate>& testInfoAggregateSet,
			vector<testInfo>& testInfoSet, bool animate);
	};
	
}