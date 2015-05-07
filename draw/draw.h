#ifndef _DRAW_H
#define _DRAW_H (1)

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TImage.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TString.h"
#include "TStyle.h"


const Color_t kDefaultColor   = kRed-7;
const Float_t default_xoffset =   1.30;
const Float_t default_yoffset =   1.80;
const Float_t default_xlow    =  0.195;
const Float_t default_yup     =  0.905;


// Functions
//------------------------------------------------------------------------------
TH1F*         GetTH1F            (TString            filename,
				  TString            name);

TH2F*         GetTH2F            (TString            filename,
				  TString            name);

void          DrawTH1F           (TPad*              canvas,
				  TH1*               hist,
				  TString            option,
				  TString            title   = "",
				  TString            xtitle  = "",
				  TString            ytitle  = "",
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          DrawTH2F           (TPad*              canvas,
				  TH2F*              hist,
				  TString            option,
				  TString            title   = "",
				  TString            xtitle  = "",
				  TString            ytitle  = "",
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          Axis               (TAxis*             axis,
				  TString            title,
				  Float_t            offset,
				  Bool_t             center = kFALSE);

void          AxisTH1            (TH1*               hist,
				  TString            xtitle,
				  TString            ytitle,
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          AxisTHStack        (THStack*           hist,
				  TString            xtitle,
				  TString            ytitle,
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          FillOpt            (TH1F*              hist,
				  Color_t            fcolor = kDefaultColor,
				  Style_t            fstyle = 1001);

void          LineOpt            (TH1F*              hist,
				  Color_t            lcolor = kDefaultColor,
				  Width_t            lwidth = 0,
				  Style_t            lstyle = kSolid);

void          MarkerOpt          (TH1F*              hist,
				  Color_t            mcolor = kDefaultColor,
				  Size_t             msize  = 0,
				  Style_t            mstyle = kFullCircle);

void          ExitPrint          (TString            method);

void          MoveStats          (TCanvas*           canvas,
				  Float_t            xoffset = 2.0);

void          SetTF1             (TF1*               f,
				  TString            option,
				  Color_t            lcolor,
				  Style_t            lstyle);

void          SetTGraph          (TGraph*            g,
				  TString            option,
				  Color_t            mcolor,
				  Style_t            mstyle,
				  TString            title  = "",
				  TString            xtitle = "",
				  TString            ytitle = "");

TLegend*      DrawLegend         (Float_t            x1,
				  Float_t            y1,
				  TH1*               hist,
				  TString            label,
				  TString            option,
				  Float_t            tsize   = 0.04,
				  Float_t            xoffset = 0.25,
				  Float_t            yoffset = 0.07);

void          DrawLogo           (TCanvas*           canvas,
				  TString            image,
				  Float_t            xlow = default_xlow,
				  Float_t            yup  = default_yup);

TCanvas*      SetCanvas          (TString            title,
				  Long_t             iCanvas = -1);

void          SaveCanvas         (TCanvas*           canvas,
				  TString            title,
				  TString            extension = ".png",
				  Long_t             iCanvas   = -1);

Float_t       Ratio              (Float_t            a,
				  Float_t            b);

Float_t       RatioError         (Float_t            a,
				  Float_t            b,
				  Float_t            aErr = -999,
				  Float_t            bErr = -999);

TH1F*         GetRatio           (TH1F*              hnum,
				  TH1F*              hden,
				  Float_t            mean = -999);

void          DrawRatioBinContent(TH1F*              hnum,
				  TH1F*              hden);

void          SetLatex           (TLatex*            latex,
				  Short_t            align,
				  Float_t            tangle,
				  Color_t            tcolor,
				  Float_t            tsize);

void          ProjectY           (TH2F*              hist);

Double_t      DeltaPhi           (Double_t           phi1,
				  Double_t           phi2);

Double_t      DeltaR             (Double_t           eta1,
				  Double_t           eta2,
				  Double_t           phi1,
				  Double_t           phi2);

TPaveText*    SetTPaveText       (Double_t           x1,
				  Double_t           y1,
				  Double_t           x2,
				  Double_t           y2);

TH1F*         Chi2Rebin          (TH1F*              hist1,
				  TH1F*              hist2,
				  Double_t&          chi2,
				  Int_t&             NDF,
				  Double_t           low       =  0.0,
				  Double_t           high      =  0.0,
				  Int_t              verbose   =    0,
				  Int_t              threshold = -999,
				  Int_t              fitParams =    0);

void          fitGaus            (TH1*               hist,
				  TF1*               f,
				  Double_t&          mean,
				  Double_t&          wide,
				  Double_t&          meanErr,
				  Double_t&          wideErr);

TGraph*       spillOver          (TH2D*              matrix,
				  Double_t&          result,
				  Bool_t             verbose = false);

void          SetHistGraphColors (TH1*               hist,
				  TGraph*            g,
				  Color_t            color,
				  Style_t            mstyle);

TGraphErrors* ProduceGraph       (TH1D*              hist,
				  TH2D*              matrix,
				  const char*        name,
				  Bool_t             verbose);

Float_t       GetWeight1D        (Float_t            x,
				  TH1D*              th1);

Float_t       GetWeight2D        (Float_t            x,
				  Float_t            y,
				  TH2D*              th2);

Bool_t        InCircumference    (Float_t            z,
				  Float_t            x,
				  const Float_t      cz,
				  const Float_t      cx,
				  const Float_t      radius);

Bool_t        InEllipse          (Float_t            xx,
				  Float_t            yy,
				  Float_t            cx,
				  Float_t            cy,
				  Float_t            rx,
				  Float_t            ry,
				  Float_t            theta);

Bool_t        InParabola         (const Float_t      xx,
				  const Float_t      yy,
				  const Float_t      x0,
				  const Float_t      y0,
				  const Float_t      factor,
				  const Float_t      alphaDeg);

void          SetupBW            ();

void          SetupColor         ();

TLegend*      MyTLegend          (Double_t           x1,
				  Double_t           y1,
				  Double_t           x2,
				  Double_t           y2);

TString       EntriesTitle       (Float_t            width,
				  UInt_t             precision,
				  TString            units);


////////////////////////////////////////////////////////////////////////////////
//
// Auxiliary functions from met2012.C
//
////////////////////////////////////////////////////////////////////////////////
void          TAxisConfigure     (TAxis*             axis,
				  TString            title);

void          TCanvasSave        (TCanvas*           canvas,
				  TString            cname);

void          TH2DrawComponents  (TH2*               h2,
				  TString            cname,
				  TString            xtitle,
				  TString            ytitle,
				  TString            format);

void          TH1Cosmetics       (TH1*               h1,
				  Color_t            fcolor = 0,
				  Style_t            fstyle = 0,
				  Color_t            lcolor = kBlack,
				  Style_t            lstyle = kSolid,
				  Width_t            lwidth = 1,
				  Color_t            mcolor = kBlack,
				  Size_t             msize  = 1.0,
				  Style_t            mstyle = kFullCircle);


#endif
