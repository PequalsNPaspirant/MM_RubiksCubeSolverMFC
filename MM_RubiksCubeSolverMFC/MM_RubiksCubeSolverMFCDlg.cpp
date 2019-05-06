
// MM_RubiksCubeSolverMFCDlg.cpp : implementation file
//

#include "stdafx.h"
#include "MM_RubiksCubeSolverMFC.h"
#include "MM_RubiksCubeSolverMFCDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CMMRubiksCubeSolverMFCDlg dialog



CMMRubiksCubeSolverMFCDlg::CMMRubiksCubeSolverMFCDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_MM_RUBIKSCUBESOLVERMFC_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMMRubiksCubeSolverMFCDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CMMRubiksCubeSolverMFCDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_RESET, &CMMRubiksCubeSolverMFCDlg::OnBnClickedReset)
	ON_BN_CLICKED(IDC_FITTOWINDOW, &CMMRubiksCubeSolverMFCDlg::OnBnClickedFittowindow)
	ON_BN_CLICKED(IDC_SCRAMBLE, &CMMRubiksCubeSolverMFCDlg::OnBnClickedScramble)
	ON_BN_CLICKED(IDC_SOLVE, &CMMRubiksCubeSolverMFCDlg::OnBnClickedSolve)
	ON_BN_CLICKED(IDC_ANIMATE, &CMMRubiksCubeSolverMFCDlg::OnBnClickedAnimate)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SPEED, &CMMRubiksCubeSolverMFCDlg::OnNMCustomdrawSpeed)
	ON_NOTIFY(TRBN_THUMBPOSCHANGING, IDC_SPEED, &CMMRubiksCubeSolverMFCDlg::OnTRBNThumbPosChangingSpeed)
	ON_CBN_SELCHANGE(IDC_COMBOSIZE, &CMMRubiksCubeSolverMFCDlg::OnCbnSelchangeCombosize)
	ON_EN_CHANGE(IDC_EDITSIZE, &CMMRubiksCubeSolverMFCDlg::OnEnChangeEditsize)
	ON_NOTIFY(NM_THEMECHANGED, IDC_SCROLLBARSIZE, &CMMRubiksCubeSolverMFCDlg::OnNMThemeChangedScrollbarsize)
END_MESSAGE_MAP()


// CMMRubiksCubeSolverMFCDlg message handlers

BOOL CMMRubiksCubeSolverMFCDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here
	HWND hWnd = NULL;
	GetDlgItem(IDC_STATIC_OPENGL, &hWnd);
	rubiksCubeSolverGui_.initialize(hWnd);
	//HWND hWnd = GetSafeHwnd();
	//rubiksCubeSolverGui_.initialize(hWnd);
	//rubiksCubeSolverGui_.redrawWindow();

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CMMRubiksCubeSolverMFCDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CMMRubiksCubeSolverMFCDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
		rubiksCubeSolverGui_.redrawWindow();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CMMRubiksCubeSolverMFCDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CMMRubiksCubeSolverMFCDlg::OnBnClickedReset()
{
	// TODO: Add your control notification handler code here
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedFittowindow()
{
	// TODO: Add your control notification handler code here
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedScramble()
{
	// TODO: Add your control notification handler code here
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedSolve()
{
	// TODO: Add your control notification handler code here
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedAnimate()
{
	// TODO: Add your control notification handler code here
}


void CMMRubiksCubeSolverMFCDlg::OnNMCustomdrawSpeed(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	*pResult = 0;
}


void CMMRubiksCubeSolverMFCDlg::OnTRBNThumbPosChangingSpeed(NMHDR *pNMHDR, LRESULT *pResult)
{
	// This feature requires Windows Vista or greater.
	// The symbol _WIN32_WINNT must be >= 0x0600.
	NMTRBTHUMBPOSCHANGING *pNMTPC = reinterpret_cast<NMTRBTHUMBPOSCHANGING *>(pNMHDR);
	// TODO: Add your control notification handler code here
	*pResult = 0;
}


void CMMRubiksCubeSolverMFCDlg::OnCbnSelchangeCombosize()
{
	// TODO: Add your control notification handler code here
}


void CMMRubiksCubeSolverMFCDlg::OnEnChangeEditsize()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialogEx::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
}


void CMMRubiksCubeSolverMFCDlg::OnNMThemeChangedScrollbarsize(NMHDR *pNMHDR, LRESULT *pResult)
{
	// This feature requires Windows XP or greater.
	// The symbol _WIN32_WINNT must be >= 0x0501.
	// TODO: Add your control notification handler code here
	*pResult = 0;
}
