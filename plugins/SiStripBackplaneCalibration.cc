/// \class SiStripBackplaneCalibration
///
/// Calibration of backplane corrections values for strip deconvolution mode. 
///
/// Note that not all algorithms support this...
///
///  \author    : Gero Flucke
///  date       : November 2012
///  $Revision: 1.6 $
///  $Date: 2012/10/25 11:07:37 $
///  (last update by $Author: flucke $)

#include "Alignment/CommonAlignmentAlgorithm/interface/IntegratedCalibrationBase.h"

//#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripLatency.h"
#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "CondFormats/SiStripObjects/interface/SiStripConfObject.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

//#include "TTree.h"
//#include "TFile.h"
//#include "TString.h"

// #include <iostream>
#include <vector>
#include <map>
//#include <sstream>
//#include <cstdio>

class SiStripBackplaneCalibration : public IntegratedCalibrationBase
{
public:
  /// Constructor
  explicit SiStripBackplaneCalibration(const edm::ParameterSet &cfg);
  
  /// Destructor
  virtual ~SiStripBackplaneCalibration();

  /// How many parameters does this calibration define?
  virtual unsigned int numParameters() const;

  // /// Return all derivatives,
  // /// default implementation uses other derivatives(..) method,
  // /// but can be overwritten in derived class for efficiency.
  // virtual std::vector<double> derivatives(const TransientTrackingRecHit &hit,
  // 					  const TrajectoryStateOnSurface &tsos,
  // 					  const edm::EventSetup &setup,
  // 					  const EventInfo &eventInfo) const;

  /// Return non-zero derivatives for x- and y-measurements with their indices by reference.
  /// Return value is their number.
  virtual unsigned int derivatives(std::vector<ValuesIndexPair> &outDerivInds,
				   const TransientTrackingRecHit &hit,
				   const TrajectoryStateOnSurface &tsos,
				   const edm::EventSetup &setup,
				   const EventInfo &eventInfo) const;

  /// Setting the determined parameter identified by index,
  /// returns false if out-of-bounds, true otherwise.
  virtual bool setParameter(unsigned int index, double value);

  /// Setting the determined parameter uncertainty identified by index,
  /// returns false if out-of-bounds, true otherwise.
  virtual bool setParameterError(unsigned int index, double error);

  /// Return current value of parameter identified by index.
  /// Returns 0. if index out-of-bounds.
  virtual double getParameter(unsigned int index) const;

  /// Return current value of parameter identified by index.
  /// Returns 0. if index out-of-bounds or if errors undetermined.
  virtual double getParameterError(unsigned int index) const;

  // /// Call at beginning of job:
  // virtual void beginOfJob(const AlignableTracker *tracker,
  // 			  const AlignableMuon *muon,
  // 			  const AlignableExtras *extras);

  /// Called at end of a the job of the AlignmentProducer.
  /// Write out determined parameters.
  virtual void endOfJob();

private:
  /// If called the first time, fill 'siStripConfObject_' //LorentzAngleInput_',
  /// later check that it has not changed.
  bool checkBackplaneInput(const edm::EventSetup &setup, const EventInfo &eventInfo);
  /// Input LorentzAngle values:
  /// - either from EventSetup of first call to derivatives(..)
  /// - NOT YET: or created from files of passed by configuration (i.e. from parallel processing)
  const SiStripConfObject* getConfObjInput();

  /// Determined parameter value for this detId (detId not treated => 0.)
  /// and the given run.
  double getParameterForDetId(unsigned int detId, edm::RunNumber_t run) const;
  /// Index of parameter for given detId (detId not treated => < 0)
  /// and the given run.
  int getParameterIndexFromDetId(unsigned int detId, edm::RunNumber_t run) const;
  /// Total number of IOVs.
  unsigned int numIovs() const;
  /// First run of iov (0 if iovNum not treated).
  edm::RunNumber_t firstRunOfIOV(unsigned int iovNum) const;

  // void writeTree(const SiStripLorentzAngle *lorentzAngle, const char *treeName) const;
  // SiStripLorentzAngle* createFromTree(const char *fileName, const char *treeName) const;

  const bool saveToDB_;
  const std::string recordNameDBwrite_;
  const std::string outFileName_;
  const std::vector<std::string> mergeFileNames_;

  edm::ESWatcher<SiStripConfObjectRcd>   watchBackPlaneRcd_;

  // const AlignableTracker *alignableTracker_;
  SiStripConfObject *siStripConfObjectInput_;
  // std::map<SiStripDetId::ModuleGeometry, double> backPlaneFractionMap_;

  std::vector<double> parameters_;
  std::vector<double> paramUncertainties_;
};

//======================================================================
//======================================================================
//======================================================================

SiStripBackplaneCalibration::SiStripBackplaneCalibration(const edm::ParameterSet &cfg)
  : IntegratedCalibrationBase(cfg),
    saveToDB_(cfg.getParameter<bool>("saveToDB")),
    recordNameDBwrite_(cfg.getParameter<std::string>("recordNameDBwrite")),
    outFileName_(cfg.getParameter<std::string>("treeFile")),
    mergeFileNames_(cfg.getParameter<std::vector<std::string> >("mergeTreeFiles")),
    //    alignableTracker_(0),
    siStripConfObjectInput_(0)
{
  // FIXME: Which granularity, leading to how many parameters?
  parameters_.resize(14, 0.); // currently one parameters per module type, start value 0.
  paramUncertainties_.resize(parameters_.size(), 0.); // dito for errors

  edm::LogInfo("Alignment") << "@SUB=SiStripBackplaneCalibration" << "Created with name "
                            << this->name() 
			    << "',\n" << this->numParameters() << " parameters to be determined."
                            << "\nsaveToDB = " << saveToDB_
                            << "\n outFileName = " << outFileName_
                            << "\n N(merge files) = " << mergeFileNames_.size();
  if (mergeFileNames_.size()) {
    edm::LogInfo("Alignment") << "@SUB=SiStripBackplaneCalibration"
                              << "First file to merge: " << mergeFileNames_[0];
  }
}
  
//======================================================================
SiStripBackplaneCalibration::~SiStripBackplaneCalibration()
{
  delete siStripConfObjectInput_;
}

//======================================================================
unsigned int SiStripBackplaneCalibration::numParameters() const
{
  return parameters_.size();
}

//======================================================================
unsigned int
SiStripBackplaneCalibration::derivatives(std::vector<ValuesIndexPair> &outDerivInds,
					 const TransientTrackingRecHit &hit,
					 const TrajectoryStateOnSurface &tsos,
					 const edm::EventSetup &setup,
					 const EventInfo &eventInfo) const
{
  // ugly const-cast:
  // But it is either only first initialisation or throwing an exception...
  const_cast<SiStripBackplaneCalibration*>(this)->checkBackplaneInput(setup, eventInfo);

  outDerivInds.clear();

  edm::ESHandle<SiStripLatency> latency;  
  setup.get<SiStripLatencyRcd>().get(latency);
  const int16_t mode = latency->singleReadOutMode();
  if (mode == 0) { // deconvolution mode: here we need the corrections

    if (hit.det()) { // otherwise 'constraint hit' or whatever
      const int index = this->getParameterIndexFromDetId(hit.det()->geographicalId(),
							 eventInfo.eventId_.run());
      if (index >= 0) { // otherwise not treated
        edm::ESHandle<MagneticField> magneticField;
        setup.get<IdealMagneticFieldRecord>().get(magneticField);
        const GlobalVector bField(magneticField->inTesla(hit.det()->surface().position()));
        const LocalVector bFieldLocal(hit.det()->surface().toLocal(bField));

        const double dZ = hit.det()->surface().bounds().thickness(); // it's a float only...
	const double tanPsi = tsos.localParameters().mixedFormatVector()[1]; //float...

	edm::ESHandle<SiStripLorentzAngle> lorentzAngleHandle;
	// Shouldn't 'deconvolution' come automatically?
	setup.get<SiStripLorentzAngleRcd>().get("deconvolution", lorentzAngleHandle);
	// Yes, mobility (= LA/By) stored in object called LA...
	const double mobility = lorentzAngleHandle->getLorentzAngle(hit.det()->geographicalId());
        // shift due to dead back plane has two parts:
	// 1) Lorentz Angle correction formula gets reduced thickness dz*(1-bp)
	// 2) 'Direct' effect is shift of effective module position in local z by bp*dz/2
	//   (see GF's presentation in alignment meeting 25.10.2012,
	//    https://indico.cern.ch/conferenceDisplay.py?confId=174266#2012-10-25)
        const double xDerivative = 0.5 * dZ * (mobility * bFieldLocal.y() - tanPsi);
	// std::cout << "derivative is " << xDerivative << " for index " << index 
	// 	  << std::endl;
	const Values derivs(xDerivative, 0.); // yDerivative = 0.
	outDerivInds.push_back(ValuesIndexPair(derivs, index));
      } 
    } else {
      edm::LogWarning("Alignment") << "@SUB=SiStripBackplaneCalibration::derivatives"
                                   << "Hit without GeomDet, skip!";
    }
  } else if (mode != 1) { // 1 == peak mode: no correction
    edm::LogWarning("Alignment") << "@SUB=SiStripBackplaneCalibration::derivatives2"
                                 << "Unknown readout mode " << mode << " !=1 && !=0, "
				 << "treat as peak, i.e. no non-zero derivatives.";
  }
  
  return outDerivInds.size();
}

//======================================================================
bool SiStripBackplaneCalibration::setParameter(unsigned int index, double value)
{
  if (index >= parameters_.size()) {
    return false;
  } else {
    parameters_[index] = value;
    return true;
  }
}

//======================================================================
bool SiStripBackplaneCalibration::setParameterError(unsigned int index, double error)
{
  if (index >= paramUncertainties_.size()) {
    return false;
  } else {
    paramUncertainties_[index] = error;
    return true;
  }
}

//======================================================================
double SiStripBackplaneCalibration::getParameter(unsigned int index) const
{
  //   if (index >= parameters_.size()) {
  //     return 0.;
  //   } else {
  //     return parameters_[index];
  //   }
  return (index >= parameters_.size() ? 0. : parameters_[index]);
}

//======================================================================
double SiStripBackplaneCalibration::getParameterError(unsigned int index) const
{
  //   if (index >= paramUncertainties_.size()) {
  //     return 0.;
  //   } else {
  //     return paramUncertainties_[index];
  //   }
  return (index >= paramUncertainties_.size() ? 0. : paramUncertainties_[index]);
}

// //======================================================================
// void SiStripBackplaneCalibration::beginOfJob(const AlignableTracker *tracker,
//                                              const AlignableMuon */*muon*/,
//                                              const AlignableExtras */*extras*/)
// {
//   alignableTracker_ = tracker;
// }


//======================================================================
void SiStripBackplaneCalibration::endOfJob()
{
  // loginfo output
  std::ostringstream out;
  out << "Parameter results\n";
  for (unsigned int iPar = 0; iPar < parameters_.size(); ++iPar) {
    out << iPar << ": " << parameters_[iPar] << " +- " << paramUncertainties_[iPar] << "\n";
  }
  edm::LogInfo("Alignment") << "@SUB=SiStripBackplaneCalibration::endOfJob" << out.str();

  // FIXME:
  // - check that input values have been identical
  // - write to TTrees
  // - write DB object


  // // now write 'input' tree
  const SiStripConfObject *input = this->getConfObjInput(); // never NULL

  if (input) {
    SiStripConfObject *output = new SiStripConfObject(*input);
//     const std::string modeS(readoutModeName_ == "peak" ? "Peak" : "Deco");
    output->update("Do_not_exist", 100.); 
    std::string shift("shift_IB1Deco");
    output->update(shift, 0.1); 
    shift = "shift_IB2Deco";
    output->update(shift, 0.15); 
    shift = "shift_OB1Deco";
    output->update(shift, 0.2); 
    shift = "shift_OB2Deco";
    output->update(shift, 0.25); 
    std::stringstream si("Input:\n");
    input->printSummary(si);

    std::stringstream so("Output:\n");
    output->printSummary(so);

    std::cout << si.str() << "\n" << so.str() << std::endl;      

  } else {
    edm::LogWarning("Alignment") << "@SUB=SiStripBackplaneCalibration::endOfJob"
                                 << "No input values!";
  }
  // const std::string treeName(this->name() + '_' + readoutModeName_ + '_');
  // this->writeTree(input, (treeName + "input").c_str());

  // if (input->getLorentzAngles().empty()) {
  //   edm::LogError("Alignment") << "@SUB=SiStripBackplaneCalibration::endOfJob"
  // 			       << "Input Lorentz angle map is empty ('"
  // 			       << readoutModeName_ << "' mode), skip writing output!";
  //   return;
  // }

  // for (unsigned int iIOV = 0; iIOV < this->numIovs(); ++iIOV) {
  //   cond::Time_t firstRunOfIOV = this->firstRunOfIOV(iIOV);
  //   SiStripLorentzAngle *output = new SiStripLorentzAngle;
  //   // Loop on map of values from input and add (possible) parameter results
  //   for (auto iterIdValue = input->getLorentzAngles().begin();
  // 	 iterIdValue != input->getLorentzAngles().end(); ++iterIdValue) {
  //     // type of (*iterIdValue) is pair<unsigned int, float>
  //     const unsigned int detId = iterIdValue->first; // key of map is DetId
  //     // Nasty: putLorentzAngle(..) takes float by reference - not even const reference!
  //     float value = iterIdValue->second + this->getParameterForDetId(detId, firstRunOfIOV);
  //     output->putLorentzAngle(detId, value); // put result in output
  //   }

  //   // Write this even for mille jobs?
  //   this->writeTree(output, (treeName + Form("result_%lld", firstRunOfIOV)).c_str());

  //   if (saveToDB_) { // If requested, write out to DB 
  //     edm::Service<cond::service::PoolDBOutputService> dbService;
  //     if (dbService.isAvailable()) {
  // 	dbService->writeOne(output, firstRunOfIOV, recordNameDBwrite_.c_str());
  // 	// no 'delete output;': writeOne(..) took over ownership
  //     } else {
  // 	delete output;
  // 	edm::LogError("BadConfig") << "@SUB=SiStripBackplaneCalibration::endOfJob"
  // 				   << "No PoolDBOutputService available, but saveToDB true!";
  //     }
  //   } else {
  //     delete output;
  //   }
  // } // end loop on IOVs
}

//======================================================================
bool SiStripBackplaneCalibration::checkBackplaneInput(const edm::EventSetup &setup,
						      const EventInfo &eventInfo)
{
  if (watchBackPlaneRcd_.check(setup)) { // new IOV of input
    edm::ESHandle<SiStripConfObject> stripConfObjHandle;
    setup.get<SiStripConfObjectRcd>().get(stripConfObjHandle);
    if (!siStripConfObjectInput_) { // nothing yet, just copy
      siStripConfObjectInput_ = new SiStripConfObject(*stripConfObjHandle);
    } else if (stripConfObjHandle->parameters != siStripConfObjectInput_->parameters) {
      // new IOV: check that BP corrections values are identical to previous
      // by comparing internal maps
      // FIXME: what if only cross talk or peak mode values differ?
      throw cms::Exception("BadInput")
	<< "SiStripBackplaneCalibration::checkBackplaneInput:\n"
	<< "Content of SiStripConfObject changed at run " << eventInfo.eventId_.run()
	<< ", but algorithm expects constant input!\n";
      return false; // not reached...
    }
  }
  
  return true;
}


// //======================================================================
// void SiStripBackplaneCalibration::checkBackPlaneFractionMap(const edm::EventSetup &setup)
// {
//   // Filled and valid? => Just return!
//   // FIXME: Why called twice?? Should we better do
//   // if (!(backPlaneFractionMap_.empty() || watchBackPlaneRcd_.check(setup))) return;
//   // or simply 
//   // if (!watchBackPlaneRcd_.check(setup))) return;
//   if (!backPlaneFractionMap_.empty() && !watchBackPlaneRcd_.check(setup)) return;

//   // All module types, see StripCPE constructor:
//   std::vector<std::pair<std::string, SiStripDetId::ModuleGeometry> > nameTypes;
//   nameTypes.push_back(std::make_pair("IB1",SiStripDetId::IB1));// 0
//   nameTypes.push_back(std::make_pair("IB2",SiStripDetId::IB2));// 1
//   nameTypes.push_back(std::make_pair("OB1",SiStripDetId::OB1));// 2
//   nameTypes.push_back(std::make_pair("OB2",SiStripDetId::OB2));// 3
//   nameTypes.push_back(std::make_pair("W1A",SiStripDetId::W1A));// 4
//   nameTypes.push_back(std::make_pair("W2A",SiStripDetId::W2A));// 5
//   nameTypes.push_back(std::make_pair("W3A",SiStripDetId::W3A));// 6
//   nameTypes.push_back(std::make_pair("W1B",SiStripDetId::W1B));// 7
//   nameTypes.push_back(std::make_pair("W2B",SiStripDetId::W2B));// 8
//   nameTypes.push_back(std::make_pair("W3B",SiStripDetId::W3B));// 9
//   nameTypes.push_back(std::make_pair("W4" ,SiStripDetId::W4 ));//10
//   nameTypes.push_back(std::make_pair("W5" ,SiStripDetId::W5 ));//11
//   nameTypes.push_back(std::make_pair("W6" ,SiStripDetId::W6 ));//12
//   nameTypes.push_back(std::make_pair("W7" ,SiStripDetId::W7 ));//13

//   edm::ESHandle<SiStripConfObject> stripConfObj;
//   setup.get<SiStripConfObjectRcd>().get(stripConfObj);

//   backPlaneFractionMap_.clear(); // Just to be sure!
//   for (auto nameTypeIt = nameTypes.begin(); nameTypeIt != nameTypes.end(); ++nameTypeIt) {
//     // See StripCPE constructor:
//     const std::string modeS(readoutModeName_ == "peak" ? "Peak" : "Deco");
//     const std::string shiftS("shift_" + nameTypeIt->first + modeS);
//     if (stripConfObj->isParameter(shiftS)) {
//       backPlaneFractionMap_[nameTypeIt->second] = stripConfObj->get<double>(shiftS);
//       std::cout << "backPlaneFraction for " << nameTypeIt->first << ": " 
//                 << backPlaneFractionMap_[nameTypeIt->second] << std::endl;
//     } else {
//       std::cout << "No " << shiftS << " in SiStripConfObject!?!";
//     }
//   }
// }

//======================================================================
const SiStripConfObject* SiStripBackplaneCalibration::getConfObjInput()
{
  // FIXME: For parallel processing in Millepede II, create SiStripConfObject
  // FIXME: from info stored in files of parallel jobs and check that they are identical.
  // FIXME: If this job has run on data, still check that LA is identical to the ones
  // FIXME: from mergeFileNames_.
  return siStripConfObjectInput_;

//   const std::string treeName(((this->name() + '_') += readoutModeName_) += "_input");
//   for (auto iFile = mergeFileNames_.begin(); iFile != mergeFileNames_.end(); ++iFile) {
//     SiStripLorentzAngle* la = this->createFromTree(iFile->c_str(), treeName.c_str());
//     // siStripLorentzAngleInput_ could be non-null from previous file of this loop
//     // or from checkLorentzAngleInput(..) when running on data in this job as well
//     if (!siStripLorentzAngleInput_ || siStripLorentzAngleInput_->getLorentzAngles().empty()) {
//       delete siStripLorentzAngleInput_; // NULL or empty
//       siStripLorentzAngleInput_ = la;
//     } else {
//       // FIXME: about comparison of maps see comments in checkLorentzAngleInput
//       if (la && !la->getLorentzAngles().empty() && // single job might not have got events
//           la->getLorentzAngles() != siStripLorentzAngleInput_->getLorentzAngles()) {
//         // Throw exception instead of error?
//         edm::LogError("NoInput") << "@SUB=SiStripBackplaneCalibration::getLorentzAnglesInput"
//                                  << "Different input values from tree " << treeName
//                                  << " in file " << *iFile << ".";
        
//       }
//       delete la;
//     }
//   }
  
//   if (!siStripLorentzAngleInput_) { // no files nor ran on events
//     siStripLorentzAngleInput_ = new SiStripLorentzAngle;
//     edm::LogError("NoInput") << "@SUB=SiStripBackplaneCalibration::getLorentzAnglesInput"
// 			     << "No input, create an empty one ('" << readoutModeName_ << "' mode)!";
//   } else if (siStripLorentzAngleInput_->getLorentzAngles().empty()) {
//     edm::LogError("NoInput") << "@SUB=SiStripBackplaneCalibration::getLorentzAnglesInput"
// 			     << "Empty result ('" << readoutModeName_ << "' mode)!";
//   }

//   return siStripLorentzAngleInput_;
}

//======================================================================
double SiStripBackplaneCalibration::getParameterForDetId(unsigned int detId,
							 edm::RunNumber_t run) const
{
  const int index = this->getParameterIndexFromDetId(detId, run);

  return (index < 0 ? 0. : parameters_[index]);
}

//======================================================================
int SiStripBackplaneCalibration::getParameterIndexFromDetId(unsigned int detId,
							    edm::RunNumber_t run) const
{
  // Return the index of the parameter that is used for this DetId.
  // If this DetId is not treated, return values < 0.
  
  // FIXME: Extend to configurable granularity? 
  //        Including treatment of run dependence?
  const SiStripDetId id(detId);
  if (id.det() == DetId::Tracker) {
    switch (id.subDetector()) {
    case SiStripDetId::TIB: // no break!
    case SiStripDetId::TID: // dito
    case SiStripDetId::TOB:
    case SiStripDetId::TEC:
      // Fine, we are in strip.
      // Note that according to DataFormats/SiStripDetId/interface/SiStripDetId.h
      // the enums of real geometries start 1.
      return id.moduleGeometry() - 1;
    default:
      ; // nothing
    }
  }

  // Either not in tracker or it should be BPIX or FPIX - not treated here!
  return -1;
}

//======================================================================
unsigned int SiStripBackplaneCalibration::numIovs() const
{
  // FIXME: Needed to include treatment of run dependence!
  return 1; 
}

//======================================================================
edm::RunNumber_t SiStripBackplaneCalibration::firstRunOfIOV(unsigned int iovNum) const
{
  // FIXME: Needed to include treatment of run dependence!
  if (iovNum < this->numIovs()) return 1;
  else return 0;
}


// //======================================================================
// void SiStripBackplaneCalibration::writeTree(const SiStripLorentzAngle *lorentzAngle,
// 					       const char *treeName) const
// {
//   if (!lorentzAngle) return;

//   TFile* file = TFile::Open(outFileName_.c_str(), "UPDATE");
//   if (!file) {
//     edm::LogError("BadConfig") << "@SUB=SiStripBackplaneCalibration::writeTree"
// 			       << "Could not open file '" << outFileName_ << "'.";
//     return;
//   }

//   TTree *tree = new TTree(treeName, treeName);
//   unsigned int id = 0;
//   float value = 0.;
//   tree->Branch("detId", &id, "detId/i");
//   tree->Branch("value", &value, "value/F");

//   for (auto iterIdValue = lorentzAngle->getLorentzAngles().begin();
//        iterIdValue != lorentzAngle->getLorentzAngles().end(); ++iterIdValue) {
//     // type of (*iterIdValue) is pair<unsigned int, float>
//     id = iterIdValue->first; // key of map is DetId
//     value = iterIdValue->second;
//     tree->Fill();
//   }
//   tree->Write();
//   delete file; // tree vanishes with the file... (?)

// }

// //======================================================================
// SiStripLorentzAngle* 
// SiStripBackplaneCalibration::createFromTree(const char *fileName, const char *treeName) const
// {
//   // Check for file existence on your own to work around
//   // https://hypernews.cern.ch/HyperNews/CMS/get/swDevelopment/2715.html:
//   TFile* file = 0;
//   FILE* testFile = fopen(fileName,"r");
//   if (testFile) {
//     fclose(testFile);
//     file = TFile::Open(fileName, "READ");
//   } // else not existing, see error below

//   TTree *tree = 0;
//   if (file) file->GetObject(treeName, tree);

//   SiStripLorentzAngle *result = 0;
//   if (tree) {
//     unsigned int id = 0;
//     float value = 0.;
//     tree->SetBranchAddress("detId", &id);
//     tree->SetBranchAddress("value", &value);

//     result = new SiStripLorentzAngle;
//     const Long64_t nEntries = tree->GetEntries();
//     for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
//       tree->GetEntry(iEntry);
//       result->putLorentzAngle(id, value);
//     }
//   } else { // Warning only since could be parallel job on no events.
//     edm::LogWarning("Alignment") << "@SUB=SiStripBackplaneCalibration::createFromTree"
//                                  << "Could not get TTree '" << treeName << "' from file '"
//                                  << fileName << (file ? "'." : "' (file does not exist).");
//   }

//   delete file; // tree will vanish with file
//   return result;
// }


//======================================================================
//======================================================================
// Plugin definition

#include "Alignment/CommonAlignmentAlgorithm/interface/IntegratedCalibrationPluginFactory.h"

DEFINE_EDM_PLUGIN(IntegratedCalibrationPluginFactory,
		   SiStripBackplaneCalibration, "SiStripBackplaneCalibration");
