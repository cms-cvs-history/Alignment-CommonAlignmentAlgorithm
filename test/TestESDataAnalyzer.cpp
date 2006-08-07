
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "FWCore/Utilities/interface/Exception.h"

//
// class decleration
//

class TestESDataAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TestESDataAnalyzer(const edm::ParameterSet&);
      ~TestESDataAnalyzer();


      virtual void analyze(const edm::Event&, const edm::EventSetup&);
   private:

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
TestESDataAnalyzer::TestESDataAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


TestESDataAnalyzer::~TestESDataAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TestESDataAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   ESHandle<TrackerGeometry> pData;

   //edm::LogInfo("TestAnalyzer") << "Retrieving tracker geometry from ES";
   iSetup.get<TrackerDigiGeometryRecord>().get(pData);
   //edm::LogInfo("TestAnalyzer") << "DONE";

   
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestESDataAnalyzer)
