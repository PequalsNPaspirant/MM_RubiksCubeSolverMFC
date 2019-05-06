
// MM_RubiksCubeSolverMFCDlg.h : header file
//

#pragma once

#include "RubiksCubeSolverGUI.h"

// CMMRubiksCubeSolverMFCDlg dialog
class CMMRubiksCubeSolverMFCDlg : public CDialogEx
{
// Construction
public:
	CMMRubiksCubeSolverMFCDlg(CWnd* pParent = nullptr);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_MM_RUBIKSCUBESOLVERMFC_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedReset();
	afx_msg void OnBnClickedFittowindow();
	afx_msg void OnBnClickedScramble();
	afx_msg void OnBnClickedSolve();
	afx_msg void OnBnClickedAnimate();
	afx_msg void OnNMCustomdrawSpeed(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnTRBNThumbPosChangingSpeed(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnCbnSelchangeCombosize();
	afx_msg void OnEnChangeEditsize();
	afx_msg void OnNMThemeChangedScrollbarsize(NMHDR *pNMHDR, LRESULT *pResult);

private:
	mm::RubiksCubeSolverGUI rubiksCubeSolverGui_;
};
