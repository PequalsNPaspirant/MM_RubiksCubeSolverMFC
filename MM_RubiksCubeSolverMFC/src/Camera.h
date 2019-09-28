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

#include "Vector3.h"

namespace mm {

	class CCamera
	{
	public:
		CCamera(void);
		~CCamera(void);
		void SetDistance(float fDistance);
		float GetDistance(void);
		void Move(float fDistance);
		CVector3 GetLookAt(void);
		void SetLookAt(CVector3 vLookAt);
		CVector3 GetUp(void);
		void SetUp(CVector3 vUp);
		void Tilt(int nDegrees);
		void Rotate(int nDegrees);
		CVector3 GetPosition();
		CVector3 GetEyePosition();
		CVector3 GetScreenNormal();
		BOOL IsFlipped();
		float GetPhi();
		float GetTheta();
		VOID SetPhi(float fPhi);
		VOID SetTheta(float fPhi);

	private:
		float m_fDistance;
		float m_fPhi; //rotation angle - angle around vertical i.e. Z-Axis
		float m_fTheta; //tilt angle - angle around screen normal i.e. X-Axis
		bool m_bFlipped;
		CVector3 m_vLookAt; //The point/location the camera is looking at
		CVector3 m_vUp; //The up direction
	};

}
