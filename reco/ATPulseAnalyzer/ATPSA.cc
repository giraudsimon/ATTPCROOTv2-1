#include "ATPSA.hh"

// FairRoot classes
#include "FairRuntimeDb.h"
#include "FairRun.h"

// STL
#include <algorithm>

using std::max_element;
using std::min_element;
using std::distance;

ClassImp(ATPSA)

ATPSA::ATPSA()
{
  fLogger = FairLogger::GetLogger();

  FairRun *run = FairRun::Instance();
  if (!run)
    fLogger -> Fatal(MESSAGE_ORIGIN, "No analysis run!");

  FairRuntimeDb *db = run -> GetRuntimeDb();
  if (!db)
    fLogger -> Fatal(MESSAGE_ORIGIN, "No runtime database!");

  fPar = (ATDigiPar *) db -> getContainer("ATDigiPar");
  if (!fPar)
    fLogger -> Fatal(MESSAGE_ORIGIN, "ATDigiPar not found!!");

  fPadPlaneX = fPar -> GetPadPlaneX();
  fPadSizeX = fPar -> GetPadSizeX();
  fPadSizeZ = fPar -> GetPadSizeZ();
  fPadRows = fPar -> GetPadRows();
  fPadLayers = fPar -> GetPadLayers();

  fNumTbs = fPar -> GetNumTbs();
  fTBTime = fPar -> GetTBTime();
  fDriftVelocity = fPar -> GetDriftVelocity();
  fMaxDriftLength = fPar -> GetDriftLength();

  fBackGroundSuppression = kFALSE;
  fBackGroundInterp = kFALSE;
  fIsPeakFinder = kFALSE;
  fIsMaxFinder = kFALSE;
  fIsBaseCorr = kFALSE;
  fIsTimeCorr = kFALSE;

  fIsGainCalibrated   = kFALSE;
  fIsJitterCalibrated = kFALSE;
  fCalibration        = new ATCalibration();

  fBField = fPar->GetBField();
  fEField = fPar->GetEField();
  fTiltAng = fPar->GetTiltAngle();
  fTB0  =  fPar->GetTB0();
  fZk   = fPar->GetZPadPlane();

  fEntTB   = (Int_t) fPar->GetTBEntrance();

  fThreshold = -1;
  fThresholdlow = -1;
  //fAuxChannels.clear();
  fThetaPad =-103.0*TMath::Pi()/180.0;

  fIniTB = 0;
  fEndTB = 512;


  std::cout<<" ==== Parameters for Pulse Shape Analysis Task ==== "<<std::endl;
  std::cout<<" ==== Magnetic Field : "<<fBField<<" T "<<std::endl;
  std::cout<<" ==== Electric Field : "<<fEField<<" V/cm "<<std::endl;
  std::cout<<" ==== Sampling Rate : "<<fTBTime<<" ns "<<std::endl;
  std::cout<<" ==== Tilting Angle : "<<fTiltAng<<" deg "<<std::endl;
  std::cout<<" ==== Drift Velocity : "<<fDriftVelocity<<" cm/us "<<std::endl;
  std::cout<<" ==== TB0 : "<<fTB0<<std::endl;

}

ATPSA::~ATPSA()
{
}

void ATPSA::SetBaseCorrection(Bool_t value)                  { fIsBaseCorr=value;}
void ATPSA::SetTimeCorrection(Bool_t value)                  { fIsTimeCorr=value;}
//void ATPSA::SetAuxChannel(std::vector<Int_t> AuxCh)        { fAuxChannels = AuxCh;}

void ATPSA::SetSimulatedEvent(TClonesArray* MCSimPointArray) 
{
  fMCSimPointArray = MCSimPointArray;
}

void
ATPSA::SetThreshold(Int_t threshold)
{
  fThreshold = threshold;
}

void
ATPSA::SetThresholdLow(Int_t thresholdlow)
{
  fThresholdlow = thresholdlow;
}

Double_t
ATPSA::CalculateX(Double_t row)
{
  return (row + 0.5)*fPadSizeX - fPadPlaneX/2.;
}

Double_t
ATPSA::CalculateZ(Double_t peakIdx)
{
    // DEPRECATED
    return (fNumTbs - peakIdx)*fTBTime*fDriftVelocity/100.;
}

Double_t
ATPSA::CalculateZGeo(Double_t peakIdx)
{

    // This function must be consistent with the re-calibrations done before.
    return fZk - (fEntTB - peakIdx)*fTBTime*fDriftVelocity/100.;


}

Double_t
ATPSA::CalculateY(Double_t layer)
{
  return (layer + 0.5)*fPadSizeZ;
}

void
ATPSA::SetBackGroundSuppression(){
  fBackGroundSuppression = kTRUE;
}

void
ATPSA::SetBackGroundInterpolation(){
  fBackGroundInterp = kTRUE;
}

void
ATPSA::SetPeakFinder(){
  fIsPeakFinder= kTRUE;

}

void
ATPSA::SetMaxFinder(){
  fIsMaxFinder= kTRUE;

}

Double_t
ATPSA::CalculateXCorr(Double_t xvalue,Int_t Tbx)  //TODO: Yes, this can be done with one stupid function but try walking on my shoes...
{


  Double_t xcorr = xvalue - fLorentzVector.X()*(Tbx-fTB0)*fTBTime*1E-2; // Convert from ns to us and cm to mm
  return xcorr;


}
Double_t
ATPSA::CalculateYCorr(Double_t yvalue,Int_t Tby)
{
  Double_t ycorr = yvalue + fLorentzVector.Y()*(Tby-fTB0)*fTBTime*1E-2;
  return ycorr;
}
Double_t
ATPSA::CalculateZCorr(Double_t zvalue,Int_t Tbz)
{
  Double_t zcorr = fLorentzVector.Z()*(Tbz-fTB0)*fTBTime*1E-2;
  return zcorr;
}

void
ATPSA::CalcLorentzVector(){

       //fDriftVelocity*=-1;// TODO: Check sign of the Vd
       Double_t ot = (fBField/fEField)*fDriftVelocity*1E4;
       Double_t front = fDriftVelocity / (1 + ot*ot);
       Double_t TiltRad = fTiltAng*TMath::Pi()/180.0;

       Double_t x = front*ot*TMath::Sin(TiltRad);
       Double_t y = front*ot*ot*TMath::Cos(TiltRad)*TMath::Sin(TiltRad);
       Double_t z = front*(1+ ot*ot*TMath::Cos(TiltRad)*TMath::Cos(TiltRad));

       fLorentzVector.SetXYZ(x,y,z);

      //  std::cout<<" Vdx : "<<x<<std::endl;
      //   std::cout<<" Vdy : "<<y<<std::endl;
      //  std::cout<<" Vdz : "<<z<<std::endl;



}

TVector3 ATPSA::RotateDetector(Double_t x,Double_t y,Double_t z,Int_t tb)
{

        // DEPRECATED because of timebucket calibration (-271.0)
        TVector3 posRot;
        TVector3 posDet;

        posRot.SetX(x*TMath::Cos(fThetaPad) - y*TMath::Sin(fThetaPad));
        posRot.SetY(x*TMath::Sin(fThetaPad) + y*TMath::Cos(fThetaPad));
        posRot.SetZ(  (-271.0+tb)*fTBTime*fDriftVelocity/100. + fZk  );

        Double_t TiltAng = -fTiltAng*TMath::Pi()/180.0;

        posDet.SetX(posRot.X());
        posDet.SetY( -(fZk-posRot.Z())*TMath::Sin(TiltAng)   + posRot.Y()*TMath::Cos(TiltAng)  );
        posDet.SetZ( posRot.Z()*TMath::Cos(TiltAng) - posRot.Y()*TMath::Sin(TiltAng)  );

        return posDet;

}

void ATPSA::SetGainCalibration(TString gainFile)        { fCalibration -> SetGainFile(gainFile);}
void ATPSA::SetJitterCalibration(TString jitterFile)    { fCalibration -> SetJitterFile(jitterFile);}
void ATPSA::SetMeanK(Int_t value)                       { fMeanK = value;}
void ATPSA::SetStddevMulThresh(Double_t value)          { fStdDev = value;}


void ATPSA::SetTBLimits(std::pair<Int_t,Int_t> limits)
{
    if(limits.first>=limits.second)
    {
      std::cout<<" Warning ATPSA::SetTBLimits -  Wrong Time Bucket limits. Setting default limits (0,512) ... "<<"\n";
      fIniTB = 0;
      fEndTB = 512;

    }else{
      fIniTB = limits.first;
      fEndTB = limits.second;

    }


}

void ATPSA::TrackMCPoints(std::multimap<Int_t,std::size_t>& map, ATHit* hit)
{

      //Find every simulated point ID for each valid pad
      std::pair<MCMapIterator, MCMapIterator> result = map.equal_range(hit->GetHitPadNum());
      for (MCMapIterator it = result.first; it != result.second; it++)
      {
        //std::cout<<fMCSimPointArray->GetEntries()<<"\n";

	int count = std::distance(result.first, result.second);
    
		//if(count>1){
               // std::cout<<" Count "<<count<<"\n";
		
		        if(fMCSimPointArray!=0){
				AtTpcPoint* MCPoint = (AtTpcPoint*) fMCSimPointArray->At(it->second);
                                ATHit::MCSimPoint mcpoint(it->second,MCPoint->GetTrackID(),MCPoint->GetEIni(),MCPoint->GetEnergyLoss(),MCPoint->GetAIni(),MCPoint->GetMassNum(),MCPoint->GetAtomicNum());
				hit->SetMCSimPoint(mcpoint);
				//std::cout << " Pad Num : "<<hit->GetHitPadNum()<<" MC Point ID : "<<it->second << std::endl;
				//std::cout << " Track ID : "<<MCPoint->GetTrackID()<<" Energy (MeV) : "<<MCPoint->GetEIni()<<" Angle (deg) : "<<MCPoint->GetAIni()<<"\n";
				//std::cout << " Mass Number : "<<MCPoint->GetMassNum()<<" Atomic Number "<<MCPoint->GetAtomicNum()<<"\n";
			}
              //  }
      }

}
