
// MM_RubiksCubeSolverMFC.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CMMRubiksCubeSolverMFCApp:
// See MM_RubiksCubeSolverMFC.cpp for the implementation of this class
//

class CMMRubiksCubeSolverMFCApp : public CWinApp
{
public:
	CMMRubiksCubeSolverMFCApp();

// Overrides
public:
	virtual BOOL InitInstance();

// Implementation

	DECLARE_MESSAGE_MAP()
};

extern CMMRubiksCubeSolverMFCApp theApp;
