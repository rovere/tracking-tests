#include "TrackingTests/TrackNtuplizer/interface/TrackNtuplizer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/eventSetupGetImplementation.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrackReco/interface/TrackResiduals.h"
#include "TMath.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <functional>
#include <sstream>
#include <iostream>
using namespace edm;
using namespace reco;

typedef math::XYZTLorentzVectorD LorentzVector;
typedef math::XYZVector Vector;
typedef math::XYZPoint Point;

using namespace std;

TrackNtuplizer::TrackNtuplizer(const edm::ParameterSet& iPara)
{
  m_outFile = 0;
  m_source = "generalTracks";
  m_associatorName = "TrackAssociatorByHits";
  beamspot_ = iPara.getParameter<edm::InputTag>( "beamspot" );
  if(iPara.exists("outfile"))m_outFile = new TFile(iPara.getParameter<string>("outfile").c_str(),"RECREATE");
  if(iPara.exists("source"))m_source = iPara.getParameter<string>("source");
  if(iPara.exists("associator"))m_associatorName = iPara.getParameter<string>("associator");
  if(iPara.exists("simSource"))m_simSource = iPara.getParameter<InputTag>("simSource");

  m_mvaSource = InputTag(m_source,"MVAVals");

  m_recoTracks = new TTree("recoTracks","",1);
  m_recoTracks->Branch("fake",&m_tvFake,"fake/F");
  m_recoTracks->Branch("duplicate",&m_tvDupl,"duplicate/F");
  m_recoTracks->Branch("iter",&m_tvIter,"iter/F");
  m_recoTracks->Branch("ndof",&m_tvNdof,"ndof/F");
  m_recoTracks->Branch("pt",&m_tvPt,"pt/F");
  m_recoTracks->Branch("phi",&m_tvPhi,"phi/F");
  m_recoTracks->Branch("nlayers",&m_tvNlayers,"nlayers/F");
  m_recoTracks->Branch("nlayers3D",&m_tvNlayers3D,"nlayers3D/F");
  m_recoTracks->Branch("nlayerslost",&m_tvNlayersLost,"nlayerslost/F");
  m_recoTracks->Branch("chi2n",&m_tvChi2n,"chi2n/F");
  m_recoTracks->Branch("chi2n_no1Dmod",&m_tvChi2n_no1Dmod,"chi2n_no1Dmod/F");
  m_recoTracks->Branch("eta",&m_tvEta,"eta/F");
  m_recoTracks->Branch("relpterr",&m_tvRelPtErr,"relpterr/F");
  m_recoTracks->Branch("nhits",&m_tvNhits,"nhits/F");
  m_recoTracks->Branch("lostin",&m_tvLostIn,"lostin/F");
  m_recoTracks->Branch("lostout",&m_tvLostOut,"lostout/F");
  m_recoTracks->Branch("minlost",&m_tvMinLost,"minlost/F");
  m_recoTracks->Branch("lostmidfrac",&m_tvLostMidFrac,"lostmidfrac/F");
  m_recoTracks->Branch("dz",&m_tvDz,"dz/F");
  m_recoTracks->Branch("d0",&m_tvD0,"d0/F");
  m_recoTracks->Branch("mvaval",&m_tvMvaVal,"mvaval/F");
  m_recoTracks->Branch("loose",&m_tvLoose,"loose/F");
  m_recoTracks->Branch("highPurity",&m_tvhighPurity,"highPurity");

  m_simTracks = new TTree("simTracks","",1);
  m_simTracks->Branch("matched",&m_tvMatched,"matched/F");
  m_simTracks->Branch("iter",&m_tvIter,"iter/F");
  m_simTracks->Branch("ndof",&m_tvNdof,"ndof/F");
  m_simTracks->Branch("pt",&m_tvPt,"pt/F");
  m_simTracks->Branch("phi",&m_tvPhi,"phi/F");
  m_simTracks->Branch("nlayers",&m_tvNlayers,"nlayers/F");
  m_simTracks->Branch("nlayers3D",&m_tvNlayers3D,"nlayers3D/F");
  m_simTracks->Branch("nlayerslost",&m_tvNlayersLost,"nlayerslost/F");
  m_simTracks->Branch("chi2n",&m_tvChi2n,"chi2n/F");
  m_simTracks->Branch("chi2n_no1Dmod",&m_tvChi2n_no1Dmod,"chi2n_no1Dmod/F");
  m_simTracks->Branch("eta",&m_tvEta,"eta/F");
  m_simTracks->Branch("relpterr",&m_tvRelPtErr,"relpterr/F");
  m_simTracks->Branch("nhits",&m_tvNhits,"nhits/F");
  m_simTracks->Branch("lostin",&m_tvLostIn,"lostin/F");
  m_simTracks->Branch("lostout",&m_tvLostOut,"lostout/F");
  m_simTracks->Branch("minlost",&m_tvMinLost,"minlost/F");
  m_simTracks->Branch("lostmidfrac",&m_tvLostMidFrac,"lostmidfrac/F");
  m_simTracks->Branch("dz",&m_tvDz,"dz/F");
  m_simTracks->Branch("d0",&m_tvD0,"d0/F");
  m_simTracks->Branch("mvaval",&m_tvMvaVal,"mvaval/F");
  m_simTracks->Branch("signal",&m_tvSignal,"signal/F");
  m_simTracks->Branch("stable",&m_tvStable,"stable/F");
  m_simTracks->Branch("loose",&m_tvLoose,"loose/F");
  m_simTracks->Branch("highPurity",&m_tvhighPurity,"highPurity");

}
//------------------------------------------------------------
//------------------------------------------------------------
TrackNtuplizer::~TrackNtuplizer()
{
  /* no op */
}
//------------------------------------------------------------
//------------------------------------------------------------
void TrackNtuplizer::beginJob()
{
}

//------------------------------------------------------------
//------------------------------------------------------------
void TrackNtuplizer::endJob()
{
  m_outFile->cd();

  m_recoTracks->Write();
  m_simTracks->Write();

  m_outFile->Close();



}
//------------------------------------------------------------
//------------------------------------------------------------
void TrackNtuplizer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  edm::ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get(m_associatorName,theAssociator);
  m_associator = theAssociator.product();
}
//------------------------------------------------------------
//------------------------------------------------------------
void TrackNtuplizer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  /* no op */
}
//------------------------------------------------------------
//------------------------------------------------------------
void TrackNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::BeamSpot> hBsp;
  iEvent.getByLabel(beamspot_, hBsp);
  reco::BeamSpot vertexBeamSpot;
  vertexBeamSpot = *hBsp;

  iSetup.get<IdealMagneticFieldRecord>().get(m_magfield);

  edm::Handle<View<Track> >handle;
  iEvent.getByLabel(m_source,handle);
  edm::Handle<View<Track> >handleHighPurity;
  iEvent.getByLabel("selectHighPurity",handleHighPurity);
  edm::Handle<TrackingParticleCollection>  simTPhandle;
  iEvent.getByLabel(m_simSource,simTPhandle);
  const TrackingParticleCollection simTracks = *(simTPhandle.product());

  edm::Handle<edm::ValueMap<float> > trackMVAStore;
  iEvent.getByLabel(m_mvaSource,trackMVAStore);

  reco::RecoToSimCollection recSimColl;
  reco::SimToRecoCollection simRecColl;

  recSimColl = m_associator->associateRecoToSim(handle,simTPhandle,&iEvent,&iSetup);
  simRecColl = m_associator->associateSimToReco(handle,simTPhandle,&iEvent,&iSetup);


  for(int i = 0; i < (int)handle->size(); i++){
    Track tk = (handle->at(i));
    m_tvFake = 1;
    m_tvDupl = 0;

    resetVars();
    setVars(tk);
    m_tvDz = tk.dz(vertexBeamSpot.position());
    m_tvD0 = tk.dxy(vertexBeamSpot.position());
    RefToBase<Track> trackRef1(handle,i);
    vector<pair<TrackingParticleRef, double> > tp1;
    if(recSimColl.find(trackRef1) != recSimColl.end())tp1 = recSimColl[trackRef1];
    if(tp1.size() > 0){
      m_tvFake = 0;
      int numAssocRecoTracks = 0;
      if(simRecColl.find(tp1[0].first) != simRecColl.end()) numAssocRecoTracks = simRecColl[tp1[0].first].size();
      if (numAssocRecoTracks >= 2) m_tvDupl = 1;
    }

    m_tvMvaVal = -99999;
    m_tvMvaVal = (*trackMVAStore)[trackRef1];

    m_recoTracks->Fill();
  }

  cout<<simTracks.size()<<" "<<simRecColl.size()<<" "<<recSimColl.size()<<" "<<trackMVAStore->size()<<endl;

  for(int i = 0; i < (int)simTracks.size(); i++){
    m_tvMatched = 0;
    TrackingParticleRef tpr(simTPhandle,i);
    TrackingParticle* tp = const_cast<TrackingParticle*>(tpr.get());
    TrackingParticle::Vector momentumTP;
    TrackingParticle::Point vertexTP;
    double dxySim(0);
    double dzSim(0);

    m_tvStable = 1;
    m_tvSignal = 0;
    if(tp->eventId().bunchCrossing() == 0 && tp->eventId().event() == 0)m_tvSignal = 1;

    momentumTP = tp->momentum();
    vertexTP = tp->vertex();
    //Calcualte the impact parameters w.r.t. PCA
    dxySim = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
    dzSim = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y())/sqrt(momentumTP.perp2()) * momentumTP.z()/sqrt(momentumTP.perp2());

    resetVars();
    std::vector<std::pair<RefToBase<Track>, double> > rt;
    if(simRecColl.numberOfAssociations(tpr) > 0) rt = (std::vector<std::pair<RefToBase<Track>, double> >) simRecColl[tpr];
    //cout<<i<<" "<<rt.size()<<endl;

    if (rt.size()!=0) {
      m_tvMatched = 1;
      reco::Track trk = *(rt.begin()->first.get());
      setVars(trk);
      m_tvMvaVal = (*trackMVAStore)[rt.begin()->first];
    }

    m_tvDz = dzSim;
    m_tvD0 = dxySim;
    m_simTracks->Fill();
  }

}
//------------------------------------------------------------
//------------------------------------------------------------
void TrackNtuplizer::setVars(reco::Track tk){
    m_tvNdof = tk.ndof();
    m_tvNlayers = tk.hitPattern().trackerLayersWithMeasurement();
    m_tvNlayers3D = tk.hitPattern().pixelLayersWithMeasurement() + tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
    m_tvNlayersLost = tk.hitPattern().trackerLayersWithoutMeasurement();

    float chi2n =  tk.normalizedChi2();
    float chi2n_no1Dmod = chi2n;

    if(tk.quality(TrackBase::TrackQuality::tight))m_tvLoose = 1;
    if(tk.quality(TrackBase::TrackQuality::highPurity))m_tvhighPurity = 1;

    int count1dhits = 0;
    for (trackingRecHit_iterator ith = tk.recHitsBegin(), edh = tk.recHitsEnd(); ith != edh; ++ith) {
      const TrackingRecHit * hit = ith->get();
      if (hit->isValid()) {
  if (typeid(*hit) == typeid(SiStripRecHit1D)) ++count1dhits;
      }
    }
    if (count1dhits > 0) {
      float chi2 = tk.chi2();
      float ndof = tk.ndof();
      chi2n = (chi2+count1dhits)/float(ndof+count1dhits);
    }
    m_tvChi2n = chi2n;
    m_tvChi2n_no1Dmod = chi2n_no1Dmod;
    m_tvEta = tk.eta();
    m_tvPt = tk.pt();
    m_tvPhi = tk.phi();
    m_tvRelPtErr = float(tk.ptError())/std::max(float(tk.pt()),0.000001f);
    m_tvNhits = tk.numberOfValidHits();
    m_tvLostIn = tk.trackerExpectedHitsInner().numberOfLostTrackerHits();
    m_tvLostOut = tk.trackerExpectedHitsOuter().numberOfLostTrackerHits();
    m_tvMinLost = std::min(m_tvLostIn,m_tvLostOut);
    m_tvLostMidFrac = float(tk.numberOfLostHits()) / float(tk.numberOfValidHits() + tk.numberOfLostHits());

    m_tvDz = 0;
    m_tvD0 = 0;

    TString algoName(tk.algoName());
    //cout<<algoName<<endl;
    algoName.ReplaceAll("iter","");
    m_tvIter = algoName.Atoi();

}
//------------------------------------------------------------
//------------------------------------------------------------
void TrackNtuplizer::resetVars(){
  m_tvNdof = -99;
  m_tvNlayers = -99;
  m_tvNlayers3D = -99;
  m_tvNlayersLost = -99;
  m_tvChi2n = -99;
  m_tvChi2n_no1Dmod = -99;
  m_tvEta = -99;
  m_tvPt = -99;
  m_tvPhi = -99;
  m_tvRelPtErr = -99;
  m_tvNhits = -99;
  m_tvLostIn = -99;
  m_tvLostOut = -99;
  m_tvMinLost = -99;
  m_tvLostMidFrac = -99;

  m_tvLoose = 0;
  m_tvhighPurity = 0;

  m_tvDz = 0;
  m_tvD0 = 0;

  m_tvIter = -99;
  m_tvMvaVal = -99;

}
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------

DEFINE_FWK_MODULE(TrackNtuplizer);
