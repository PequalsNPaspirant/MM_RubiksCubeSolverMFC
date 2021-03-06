
// MM_RubiksCubeSolverMFCDlg.h : header file
//

#pragma once

#include "RubiksCubeSolverGUI.h"

// CMMRubiksCubeSolverMFCDlg dialog
class CMMRubiksCubeSolverMFCDlg : public CDialogEx
{
// Construction
public:
	static CMMRubiksCubeSolverMFCDlg& getMainDailog();
	CMMRubiksCubeSolverMFCDlg(CWnd* pParent = nullptr);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_MM_RUBIKSCUBESOLVERMFC_DIALOG };
#endif

// Implementation

public:
	afx_msg void OnBnClickedReset();
	afx_msg void OnBnClickedFittowindow();
	afx_msg void OnBnClickedScramble();
	afx_msg void OnBnClickedSolve();
	afx_msg void OnBnClickedAnimate();
	afx_msg void OnEnChangeEditsize();

	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint point);
	afx_msg void OnSize(UINT state, int cx, int cy);
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnEnChangeEditScrambleSteps();
	afx_msg void OnBnClickedScrambleGenerate();
	afx_msg void OnEnChangeEditScrambleAlgo();

	//bool CreateYesNoDialog(const string& message);
	//void CreateOkDialog(const string& message);

	void handleRubiksCubeSizeChange();
	void displayMessage(const vector<string>& messages);

protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	void OnClose();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()

	CToolTipCtrl toolTipCtrl_;

private:
	enum UpdateDataDirection
	{
		FromVariablesToControls = FALSE,
		FromControlsToVariables = TRUE
	};

	mm::RubiksCubeSolverGUI rubiksCubeSolverGui_;
	HWND hWndOpenGl;
	bool g_bMouseDown;
	bool bRMouseDown;
	int g_nPrevX;
	int g_nPrevY;
	bool animate_;
	unsigned int rubiksCubeSize_ = 3;
	unsigned int rubiksCubeChangedSize_ = 3;
	unsigned int numScramblingSteps_ = 25;
	CString scramblingAlgo_;
	CString cubeType_;
	CListBox listBoxCtrl_;

	//Unused
	int g_nHitCount{ 0 };	
public:
	afx_msg void OnBnClickedPause();
private:
	CButton togglePauseResume_;
	bool paused_{ false };
public:
	afx_msg void OnBnClickedRunTests();
	CComboBox comboBoxCubeType;
	afx_msg void OnCbnSelchangeCombotype();
	virtual BOOL PreTranslateMessage(MSG* pMsg);
};
