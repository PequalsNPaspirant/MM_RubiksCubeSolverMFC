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

#include <math.h>

namespace mm {

	const double PI = 3.1415926535897932384626433832795;

	class CVector3
	{
	public:
		CVector3();
		CVector3(double x, double y, double z);
		~CVector3();
		double x;
		double y;
		double z;
		CVector3 operator+(const CVector3& v);
		CVector3 operator-(const CVector3& v);
		CVector3 operator-();
		CVector3 operator*(double s);
		CVector3 operator/(double s);
		double operator*(const CVector3& v);
		CVector3 operator^(const CVector3& v);
		CVector3& operator+=(const CVector3& v);
		CVector3& operator-=(const CVector3& v);
		CVector3& operator/=(double s);
		CVector3& operator*=(double s);
		CVector3 Unit(void);
		double Length(void);
		CVector3 OrthogonalTo(void);
		double GetAngle(CVector3& v);
		VOID ToFloatArray(float* position);

		static const CVector3 XAxis;
		static const CVector3 YAxis;
		static const CVector3 ZAxis;

		bool operator==(const CVector3& rhs) const
		{
			return fabs(x - rhs.x) < 0.0001
				&& fabs(y - rhs.y) < 0.0001
				&& fabs(z - rhs.z) < 0.0001;
		}

		bool operator!=(const CVector3& rhs) const
		{
			return !(*this == rhs);
		}

	private:
		double LengthSquared(void);
	};

}