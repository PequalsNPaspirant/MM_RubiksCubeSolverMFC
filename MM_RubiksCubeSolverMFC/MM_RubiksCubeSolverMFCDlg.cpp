
// MM_RubiksCubeSolverMFCDlg.cpp : implementation file
//

#include "stdafx.h"
#include "resource.h"		// main symbols
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


CMMRubiksCubeSolverMFCDlg& CMMRubiksCubeSolverMFCDlg::getMainDailog()
{
	static CMMRubiksCubeSolverMFCDlg dlg;
	return dlg;
}

CMMRubiksCubeSolverMFCDlg::CMMRubiksCubeSolverMFCDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_MM_RUBIKSCUBESOLVERMFC_DIALOG, pParent),
	animate_(true)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMMRubiksCubeSolverMFCDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);

	DDX_Text(pDX, IDC_EDITSIZE, rubiksCubeChangedSize_);
	DDX_Text(pDX, IDC_EDIT_SCRAMBLE_STEPS, numScramblingSteps_);
	DDX_Text(pDX, IDC_EDIT_SCRAMBLE_ALGO, scramblingAlgo_);
	DDX_Control(pDX, IDC_LIST, listBoxCtrl_);
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
	//ON_NOTIFY(NM_CUSTOMDRAW, IDC_SPEED, &CMMRubiksCubeSolverMFCDlg::OnNMCustomdrawSpeed)
	//ON_NOTIFY(TRBN_THUMBPOSCHANGING, IDC_SPEED, &CMMRubiksCubeSolverMFCDlg::OnTRBNThumbPosChangingSpeed)
	ON_CBN_SELCHANGE(IDC_COMBOSIZE, &CMMRubiksCubeSolverMFCDlg::OnCbnSelchangeCombosize)
	ON_EN_CHANGE(IDC_EDITSIZE, &CMMRubiksCubeSolverMFCDlg::OnEnChangeEditsize)
	ON_EN_CHANGE(IDC_EDIT_SCRAMBLE_STEPS, &CMMRubiksCubeSolverMFCDlg::OnEnChangeEditScrambleSteps)
	//ON_NOTIFY(NM_THEMECHANGED, IDC_SCROLLBARSIZE, &CMMRubiksCubeSolverMFCDlg::OnNMThemeChangedScrollbarsize)
	ON_WM_LBUTTONUP()
	ON_WM_LBUTTONDOWN()
	ON_WM_MOUSEMOVE()
	ON_WM_MOUSEWHEEL()
	ON_WM_SIZE()
	ON_WM_HSCROLL()
	ON_WM_VSCROLL()
	ON_WM_CLOSE()
	
	//unused
	ON_EN_KILLFOCUS(IDC_EDITSIZE, &CMMRubiksCubeSolverMFCDlg::OnEnKillfocusEditsize)
	
	ON_BN_CLICKED(IDC_SCRAMBLE_GENERATE, &CMMRubiksCubeSolverMFCDlg::OnBnClickedScrambleGenerate)
	ON_EN_CHANGE(IDC_EDIT_SCRAMBLE_ALGO, &CMMRubiksCubeSolverMFCDlg::OnEnChangeEditScrambleAlgo)
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
	//HWND hWndAnimate = NULL;
	//GetDlgItem(IDC_ANIMATE, &hWndAnimate);
	//if(animate_)
	//	((CButton*)hWndAnimate)->SetCheck(1);
	//else
	//	((CButton*)hWndAnimate)->SetCheck(0);
	CheckDlgButton(IDC_ANIMATE, animate_);

	//initial animation speed is 40
	int intialSpeed = 40;
	HWND hWndSpeed = NULL;
	GetDlgItem(IDC_SPEED, &hWndSpeed);
	//((CScrollBar*)GetDlgItem(IDC_SPEED))->SetScrollPos(intialSpeed);
	//SetScrollRange(IDC_SPEED, 1, 100);
	//SetScrollPos(IDC_SPEED, intialSpeed);
	//CScrollBar* p = GetScrollBarCtrl(IDC_SPEED);
	//p->SetScrollRange(15, 90);
	//p->SetScrollPos(intialSpeed);
	::SendMessage(hWndSpeed, TBM_SETPOS, TRUE, intialSpeed);
	rubiksCubeSolverGui_.setAnimationSpeed(intialSpeed);

	hWndOpenGl = NULL;
	GetDlgItem(IDC_STATIC_OPENGL, &hWndOpenGl);
	rubiksCubeSolverGui_.initialize(hWndOpenGl);
	//HWND hWnd = GetSafeHwnd();
	//rubiksCubeSolverGui_.initialize(hWnd);
	//rubiksCubeSolverGui_.redrawWindow();

	//HWND hWndListBox = NULL;
	//GetDlgItem(IDC_LIST, &hWndListBox);
	//listBox_.Attach(hWndListBox);

	//listBoxCtrl_.AddString(L"sample 1");

	//listBoxCtrl_.SetHorizontalExtent(10000);

	rubiksCubeSolverGui_.displayUpdatedStats();
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

void CMMRubiksCubeSolverMFCDlg::OnClose()
{
	rubiksCubeSolverGui_.exitUI();

	CDialogEx::OnClose();
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
	rubiksCubeSolverGui_.resetRubiksCube();
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedFittowindow()
{
	// TODO: Add your control notification handler code here
	rubiksCubeSolverGui_.fitToScreen();
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedScramble()
{
	// TODO: Add your control notification handler code here
	CT2CA pszConvertedAnsiString(scramblingAlgo_);
	rubiksCubeSolverGui_.Scramble(std::string(pszConvertedAnsiString), animate_);
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedSolve()
{
	// TODO: Add your control notification handler code here
	rubiksCubeSolverGui_.Solve(animate_);
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedAnimate()
{
	// TODO: Add your control notification handler code here
	animate_ = !animate_;
	rubiksCubeSolverGui_.setAnimate(animate_);
}

void CMMRubiksCubeSolverMFCDlg::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	//if (nSBCode == SB_THUMBPOSITION)
	if (nSBCode == SB_THUMBTRACK)
	{
		CSliderCtrl* pSlider = reinterpret_cast<CSliderCtrl*>(pScrollBar);
		int NewPos = pSlider->GetPos();
		rubiksCubeSolverGui_.setAnimationSpeed(NewPos);
	}

	CDialogEx::OnHScroll(nSBCode, nPos, pScrollBar);
}

void CMMRubiksCubeSolverMFCDlg::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	if (nSBCode == SB_LINEUP)
		++rubiksCubeChangedSize_;
	else if(nSBCode == SB_LINEDOWN)
		--rubiksCubeChangedSize_;

	if (rubiksCubeChangedSize_ != rubiksCubeSize_)
	{
		UpdateData(UpdateDataDirection::FromVariablesToControls);
		handleRubiksCubeSizeChange();
	}

	CDialogEx::OnVScroll(nSBCode, nPos, pScrollBar);
}

//void CMMRubiksCubeSolverMFCDlg::OnNMCustomdrawSpeed(NMHDR *pNMHDR, LRESULT *pResult)
//{
//	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
//	// TODO: Add your control notification handler code here
//	*pResult = 0;
//}


//void CMMRubiksCubeSolverMFCDlg::OnTRBNThumbPosChangingSpeed(NMHDR *pNMHDR, LRESULT *pResult)
//{
//	// This feature requires Windows Vista or greater.
//	// The symbol _WIN32_WINNT must be >= 0x0600.
//	NMTRBTHUMBPOSCHANGING *pNMTPC = reinterpret_cast<NMTRBTHUMBPOSCHANGING *>(pNMHDR);
//	// TODO: Add your control notification handler code here
//	*pResult = 0;
//}


void CMMRubiksCubeSolverMFCDlg::OnCbnSelchangeCombosize()
{
	// TODO: Add your control notification handler code here
}

//void CMMRubiksCubeSolverMFCDlg::OnNMThemeChangedScrollbarsize(NMHDR *pNMHDR, LRESULT *pResult)
//{
//	// This feature requires Windows XP or greater.
//	// The symbol _WIN32_WINNT must be >= 0x0501.
//	// TODO: Add your control notification handler code here
//	*pResult = 0;
//}

void CMMRubiksCubeSolverMFCDlg::OnLButtonUp(UINT nFlags, CPoint point)
{
	g_bMouseDown = false;

	CDialogEx::OnLButtonDown(nFlags, point);
}

void CMMRubiksCubeSolverMFCDlg::OnLButtonDown(UINT nFlags, CPoint point)
{
	// set up tracking for when mouse leaves window
	TRACKMOUSEEVENT tme;
	tme.cbSize = sizeof(TRACKMOUSEEVENT);
	tme.dwFlags = TME_LEAVE;
	tme.hwndTrack = hWndOpenGl;
	TrackMouseEvent(&tme);

	g_bMouseDown = true;
	g_nPrevX = point.x;
	g_nPrevY = point.y;

	CDialogEx::OnLButtonUp(nFlags, point);
}

void CMMRubiksCubeSolverMFCDlg::OnMouseMove(UINT nFlags, CPoint point)
{
	//CDC* dc;
	//dc = GetDC();
	//CString str;
	//str.Format(_T("x: %d  y: %d"), point.x, point.y);
	//dc->TextOut(10, 10, str);
	//ReleaseDC(dc);

	if (!g_bMouseDown)
	{
		/*
		if ((g_nHitCount = GetSelectedObjects(x, y, g_rWnd.right, g_rWnd.bottom)) > 0)
		SetCursor(g_hHand);
		else
		SetCursor(g_hArrow);
		*/
	}

	// moving camera
	else if (g_nHitCount == 0 && g_bMouseDown)
	{
		//if (x < g_nPrevX)
		//	scene_.g_cCamera.Rotate(-5);
		//else if (x > g_nPrevX)
		//	scene_.g_cCamera.Rotate(5);

		//if (y < g_nPrevY)
		//	scene_.g_cCamera.Tilt(5);
		//else if (y > g_nPrevY)
		//	scene_.g_cCamera.Tilt(-5);

		int rotate = 0;
		if (point.x < g_nPrevX)
			rotate = (-5);
		else if (point.x > g_nPrevX)
			rotate = (5);

		int tilt = 0;
		if (point.y < g_nPrevY)
			tilt = (-5);
		else if (point.y > g_nPrevY)
			tilt = (5);

		//activateRenderingThread();
		//redrawWindow();

		if (rotate != 0 || tilt != 0)
			rubiksCubeSolverGui_.OnMouseMove(rotate, tilt);

		g_nPrevX = point.x;
		g_nPrevY = point.y;
	}

	CDialogEx::OnMouseMove(nFlags, point);
}

BOOL CMMRubiksCubeSolverMFCDlg::OnMouseWheel(UINT nFlags, short zDelta, CPoint point)
{
	int rotations = zDelta / WHEEL_DELTA;
	float distance = rotations * 1.0f;

	//scene_.g_cCamera.Move(distance);

	//rubiksCubeSolverGui_.activateRenderingThread();
	//redrawWindow();
	rubiksCubeSolverGui_.OnMouseWheel(distance);

	return CDialogEx::OnMouseWheel(nFlags, zDelta, point);
}

void CMMRubiksCubeSolverMFCDlg::OnSize(UINT state, int cx, int cy)
{
	//if (!g_bFullScreen)
	{
		//WND_WIDTH = cx;
		//WND_HEIGHT = cy;
		////rubiksCubeSolverGui_.scene_.sizeOpenGlScreen(cx, cy - messageWndHeight);
		//rubiksCubeSolverGui_.OnSize(cx, cy - messageWndHeight);

		//rubiksCubeSolverGui_.OnSize(cx, cy);
		RECT g_rWnd;
		::GetClientRect(hWndOpenGl, &g_rWnd);
		rubiksCubeSolverGui_.OnSize(g_rWnd.right - g_rWnd.left, g_rWnd.bottom - g_rWnd.top);
		//RedrawWindow(&g_rWnd);
		//MoveWindow(g_hWndMessage, 0, 0, cx, messageWndHeight, false);
		//GetClientRect(hWnd, &g_rWnd);
	}
	//rubiksCubeSolverGui_.activateRenderingThread(true);
	//redrawWindow();
	//rubiksCubeSolverGui_.redrawWindow();

	CDialogEx::OnSize(state, cx, cy);
	//rubiksCubeSolverGui_.OnSize(cx, cy);
}

void CMMRubiksCubeSolverMFCDlg::OnEnKillfocusEditsize()
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
	UpdateData(UpdateDataDirection::FromControlsToVariables);
	handleRubiksCubeSizeChange();	
}

void CMMRubiksCubeSolverMFCDlg::handleRubiksCubeSizeChange()
{
	bool success = rubiksCubeSolverGui_.setRubiksCubeSize(rubiksCubeChangedSize_);
	if (success)
		rubiksCubeSize_ = rubiksCubeChangedSize_;
	else
	{
		rubiksCubeChangedSize_ = rubiksCubeSize_;
		UpdateData(UpdateDataDirection::FromVariablesToControls);
	}
}

//bool CMMRubiksCubeSolverMFCDlg::CreateYesNoDialog(const string& message)
//{
//	wstring wMessage(message.begin(), message.end());
//	int msgboxID = MessageBox(
//		(LPCWSTR)wMessage.c_str(),
//		//(LPCWSTR)L"Resource not available\nDo you want to try again?",
//		(LPCWSTR)L"Choose Option",
//		//MB_ICONWARNING | MB_CANCELTRYCONTINUE | MB_DEFBUTTON2
//		MB_YESNO
//		);
//	switch (msgboxID)
//	{
//	case IDYES:
//		// TODO: add code
//		break;
//	case IDNO:
//		// TODO: add code
//		break;
//	}
//
//	return msgboxID;
//}
//
//void CMMRubiksCubeSolverMFCDlg::CreateOkDialog(const string& message)
//{
//	wstring wMessage(message.begin(), message.end());
//	int msgboxID = MessageBox(
//		(LPCWSTR)wMessage.c_str(),
//		//(LPCWSTR)L"Resource not available\nDo you want to try again?",
//		(LPCWSTR)L"Warning",
//		//MB_ICONWARNING | MB_CANCELTRYCONTINUE | MB_DEFBUTTON2
//		MB_OK
//	);
//}


void CMMRubiksCubeSolverMFCDlg::OnEnChangeEditScrambleSteps()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialogEx::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData(UpdateDataDirection::FromControlsToVariables);
}


void CMMRubiksCubeSolverMFCDlg::OnBnClickedScrambleGenerate()
{
	// TODO: Add your control notification handler code here
	string scramblingAlgo = rubiksCubeSolverGui_.generateScramblingAlgo(numScramblingSteps_);
	scramblingAlgo_ = CString(scramblingAlgo.c_str());
	UpdateData(UpdateDataDirection::FromVariablesToControls);
}


void CMMRubiksCubeSolverMFCDlg::OnEnChangeEditScrambleAlgo()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialogEx::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData(UpdateDataDirection::FromControlsToVariables);
}

void CMMRubiksCubeSolverMFCDlg::displayMessage(const vector<string>& messages)
{
	listBoxCtrl_.ResetContent();

	int maxIndex = 0;
	TEXTMETRIC tm;
	CDC* pDC = listBoxCtrl_.GetDC();
	CFont* pFont = listBoxCtrl_.GetFont();
	// Select the listbox font, save the old font
	CFont* pOldFont = pDC->SelectObject(pFont);
	// Get the text metrics for avg char width
	pDC->GetTextMetrics(&tm);

	for(int i = 0; i < messages.size(); ++i)
	{
		const string& str = messages[i];
		CString cstr = CString(str.c_str());
		listBoxCtrl_.AddString(cstr);

		if (messages[maxIndex].length() < messages[i].length())
			maxIndex = i;
	}

	CString cstr = CString(messages[maxIndex].c_str());
	//listBoxCtrl_.GetText(i, cstr);
	CSize sz = pDC->GetTextExtent(cstr);

	// Add the avg width to prevent clipping
	sz.cx += tm.tmAveCharWidth;

	// Select the old font back into the DC
	pDC->SelectObject(pOldFont);
	listBoxCtrl_.ReleaseDC(pDC);

	int extraClearance = 10;
	int curr = listBoxCtrl_.GetHorizontalExtent();
	if (sz.cx > curr)
		listBoxCtrl_.SetHorizontalExtent(sz.cx + extraClearance);
}