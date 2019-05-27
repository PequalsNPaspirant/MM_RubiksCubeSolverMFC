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

#include <vector>
#include <string>
using namespace std;

//#include "RubiksCubeSolverMainWindow.h"
#include "MM_RubiksCubeSolverMFCDlg.h"

namespace mm {

	class RubiksCubeSolverUtils
	{
	public:

		static bool CreateYesNoDialog(const string& message)
		{
			wstring wMessage(message.begin(), message.end());
			if (MessageBox(NULL, wMessage.c_str(), L"Choose Option:", MB_YESNO | MB_ICONQUESTION | MB_APPLMODAL) == IDYES)
				return true;
			else
				return false;
		}	

		static void CreateOkDialog(const string& message)
		{
			wstring wMessage(message.begin(), message.end());
			MessageBox(NULL, wMessage.c_str(), L"Warning!", MB_OK | MB_ICONINFORMATION | MB_APPLMODAL);
		}

		static void RunTimeAssert(bool expression, const string& msg = "")
		{
			if (!expression)
			{
				//cout << msg;
				//RubiksCubeSolverMainWindow::getInstance().redrawWindow();
				bool bContinue = RubiksCubeSolverUtils::CreateYesNoDialog("Assertion failure: " + msg + "\nDo you want to continue?");
				if (!bContinue)
				{
					int* nullPointer = nullptr;
					*nullPointer = 0;
				}
			}
		}

		static void displayMessage(const string& message)
		{
			std::vector<string> vmsg{ message };
			CMMRubiksCubeSolverMFCDlg::getMainDailog().displayMessage(vmsg);
		}
	};


}