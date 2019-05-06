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
#include "Vector3.h"
#include <math.h>

namespace mm {

	const CVector3 CVector3::XAxis{ 1.0, 0.0, 0.0 };
	const CVector3 CVector3::YAxis{ 0.0, 1.0, 0.0 };
	const CVector3 CVector3::ZAxis{ 0.0, 0.0, 1.0 };

	CVector3::CVector3(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	CVector3::CVector3()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}

	CVector3::~CVector3(void)
	{

	}

	CVector3 CVector3::operator +(const CVector3 &v)
	{
		double x = this->x + v.x;
		double y = this->y + v.y;
		double z = this->z + v.z;

		return CVector3(x, y, z);
	}

	CVector3 CVector3::operator-(const CVector3 &v)
	{
		return CVector3(x - v.x, y - v.y, z - v.z);
	}

	CVector3 CVector3::operator-()
	{
		return CVector3(-x, -y, -z);
	}

	CVector3 CVector3::operator*(double s)
	{
		return CVector3(s * x, s * y, s * z);
	}

	CVector3 CVector3::operator/(double s)
	{
		return CVector3(x / s, y / s, z / s);
	}

	// dot product
	double CVector3::operator*(const CVector3 &v)
	{
		return (x * v.x) + (y * v.y) + (z * v.z);
	}

	// cross product
	CVector3 CVector3::operator^(const CVector3 &v)
	{
		return CVector3(
			y * v.z - z * v.y,
			z * v.x - x * v.z,
			x * v.y - y * v.x);
	}

	CVector3& CVector3::operator+=(const CVector3& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;

		return *this;
	}

	CVector3& CVector3::operator-=(const CVector3& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;

		return *this;
	}

	CVector3& CVector3::operator/=(double s)
	{
		x /= s;
		y /= s;
		z /= s;

		return *this;
	}

	CVector3& CVector3::operator*=(double s)
	{
		x *= s;
		y *= s;
		z *= s;

		return *this;
	}

	CVector3 CVector3::Unit(void)
	{
		double len = LengthSquared();
		return len > 0.0 ? *this / sqrt(len) : CVector3();
	}

	double CVector3::Length(void)
	{
		return sqrt(LengthSquared());
	}

	CVector3 CVector3::OrthogonalTo(void)
	{
		if (x == 0)
			return CVector3(1, 0, 0);
		else
			return CVector3(y, -x, 0);
	}

	double CVector3::LengthSquared(void)
	{
		return pow(x, 2) + pow(y, 2) + pow(z, 2);
	}

	double CVector3::GetAngle(CVector3 &v)
	{
		CVector3 me(x, y, z);

		double dotProduct = me * v;

		double cosRatio = dotProduct / (me.Length() * v.Length());
		return acos(cosRatio);
	}

	VOID CVector3::ToFloatArray(float* position)
	{
		position[0] = (float)x;
		position[1] = (float)y;
		position[2] = (float)z;
		position[3] = 1.0f;
	}

}