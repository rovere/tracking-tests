// -*- C++ -*-
//
// Package:    BXFilter/BunchFilter
// Class:      BunchFilter
//
/**\class BunchFilter BunchFilter.cc BXFilter/BunchFilter/plugins/BunchFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Rovere
//         Created:  Tue, 03 Nov 2015 14:42:47 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class BunchFilter : public edm::EDFilter {
 public:
  explicit BunchFilter(const edm::ParameterSet&);
  ~BunchFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  virtual void beginJob() override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;
  // virtual void endLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;

  // ----------member data ---------------------------
  bool filter_;
  int lowest_bx_;
  int highest_bx_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BunchFilter::BunchFilter(const edm::ParameterSet& iConfig)
    : filter_(iConfig.getParameter<bool>("filter")),
      lowest_bx_(iConfig.getParameter<int>("lowest_bx")),
      highest_bx_(iConfig.getParameter<int>("highest_bx")) {
  // now do what ever initialization is needed
}

BunchFilter::~BunchFilter() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool BunchFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  if (!filter_) return false;

  int current_bx = iEvent.eventAuxiliary().bunchCrossing();
  if (current_bx >= lowest_bx_ and current_bx <= highest_bx_) return true;
  return false;
}

// ------------ method called once each job just before starting event loop
// ------------
void BunchFilter::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void BunchFilter::endJob() {}

// ------------ method called when starting to processes a run  ------------
/*
void
BunchFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
BunchFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block
// ------------
/*
void
BunchFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup
const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block
// ------------
/*
void
BunchFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup
const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void BunchFilter::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
// define this as a plug-in
DEFINE_FWK_MODULE(BunchFilter);
