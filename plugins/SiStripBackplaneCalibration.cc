/// \class SiStripBackplaneCalibration
///
/// Calibration of back plane corrections for the strip tracker, integrated in the
/// alignment algorithms. Note that not all algorithms support this...
///
/// Usally use one instance for deco mode data since peak mode should give the
//  reference - but the strip local reco foresees also calibration parameters
//  for peak data like for the Lorentz angle.
///
///  \author    : Gero Flucke
///  date       : November 2012
///  $Revision: 1.1.2.3 $
///  $Date: 2013/04/22 08:26:47 $
///  (last update by $Author: flucke $)

#include "Alignment/CommonAlignmentAlgorithm/interface/IntegratedCalibrationBase.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/TkModuleGroupSelector.h"
#include "Alignment/CommonAlignmentAlgorithm/plugins/SiStripReadoutModeEnums.h"

#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripLatency.h"
#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
// #include "CalibTracker/Records/interface/SiStripDependentRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripBackPlaneCorrection.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "TTree.h"
#include "TFile.h"
#include "TString.h"

// #include <iostream>
#include <boost/assign/list_of.hpp>
#include <vector>
#include <map>
#include <sstream>
#include <cstdio>
#include <functional>

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
  /// If called the first time, fill 'siStripBackPlaneCorrInput_',
  /// later check that Backplane has not changed.
  bool checkBackPlaneCorrectionInput(const edm::EventSetup &setup, const EventInfo &eventInfo);
  /// Input BackPlaneCorrection values:
  /// - either from EventSetup of first call to derivatives(..)
  /// - or created from files of passed by configuration (i.e. from parallel processing)
  const SiStripBackPlaneCorrection* getBackPlaneCorrectionInput();

  /// Determined parameter value for this detId (detId not treated => 0.)
  /// and the given run.
  double getParameterForDetId(unsigned int detId, edm::RunNumber_t run) const;

  void writeTree(const SiStripBackPlaneCorrection *backPlaneCorr,
		 const std::map<unsigned int,float> &errors, const char *treeName) const;
  SiStripBackPlaneCorrection* createFromTree(const char *fileName, const char *treeName) const;

  const std::string readoutModeName_;
  int16_t readoutMode_;
  const bool saveToDB_;
  const std::string recordNameDBwrite_;
  const std::string outFileName_;
  const std::vector<std::string> mergeFileNames_;

  edm::ESWatcher<SiStripBackPlaneCorrectionRcd> watchBackPlaneCorrRcd_;

  // const AlignableTracker *alignableTracker_;
  SiStripBackPlaneCorrection *siStripBackPlaneCorrInput_;

  std::vector<double> parameters_;
  std::vector<double> paramUncertainties_;

  TkModuleGroupSelector BPmoduleselector;
};

//======================================================================
//======================================================================
//======================================================================

SiStripBackplaneCalibration::SiStripBackplaneCalibration(const edm::ParameterSet &cfg)
  : IntegratedCalibrationBase(cfg),
    readoutModeName_(cfg.getParameter<std::string>("readoutMode")),
    saveToDB_(cfg.getParameter<bool>("saveToDB")),
    recordNameDBwrite_(cfg.getParameter<std::string>("recordNameDBwrite")),
    outFileName_(cfg.getParameter<std::string>("treeFile")),
    mergeFileNames_(cfg.getParameter<std::vector<std::string> >("mergeTreeFiles")),
    //    alignableTracker_(0),
    siStripBackPlaneCorrInput_(0),
    BPmoduleselector(
                     cfg.getParameter<edm::VParameterSet>("BackplaneGranularity")
                     )
{
  //specify the sub-detectors for which the LA is determined
  const std::vector<int> sdets = boost::assign::list_of(SiStripDetId::TIB)(SiStripDetId::TOB); //no TEC,TID
  BPmoduleselector.SetSubDets(sdets);
  
  // SiStripLatency::singleReadOutMode() returns
  // 1: all in peak, 0: all in deco, -1: mixed state
  // (in principle one could treat even mixed state APV by APV...)
  if (readoutModeName_ == "peak") {
    readoutMode_ = kPeakMode;
  } else if (readoutModeName_ == "deconvolution") {
    readoutMode_ = kDeconvolutionMode;
  } else {
    throw cms::Exception("BadConfig")
	  << "SiStripBackplaneCalibration:\n" << "Unknown mode '" 
	  << readoutModeName_ << "', should be 'peak' or 'deconvolution' .\n";
  }

  // FIXME: Which granularity, leading to how many parameters?
  parameters_.resize(2, 0.); // currently two parameters (TIB, TOB), start value 0.
  paramUncertainties_.resize(2, 0.); // dito for errors

  edm::LogInfo("Alignment") << "@SUB=SiStripBackplaneCalibration" << "Created with name "
                            << this->name() << " for readout mode '" << readoutModeName_
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
  //  std::cout << "Destroy SiStripBackplaneCalibration named " << this->name() << std::endl;
  delete siStripBackPlaneCorrInput_;
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
  const_cast<SiStripBackplaneCalibration*>(this)->checkBackPlaneCorrectionInput(setup, eventInfo);

  outDerivInds.clear();

  edm::ESHandle<SiStripLatency> latency;  
  setup.get<SiStripLatencyRcd>().get(latency);
  const int16_t mode = latency->singleReadOutMode();
  // std::cout << "SiStripBackplaneCalibration in mode '" << readoutModeName_ << "' finds mode "
  //           << mode << std::endl;
  //
  // std::vector<SiStripLatency::Latency> latencies = const_cast<SiStripLatency*>(latency.product())->allUniqueLatencyAndModes();
  // for (auto i = latencies.begin(); i != latencies.end(); ++i) {
  //   std::cout << static_cast<int>(i->latency) << ", mode "
  //             << static_cast<int>(i->mode) << ", id_apv " << i->detIdAndApv
  //             << std::endl;
  // }

  if(mode == readoutMode_) {
    if (hit.det()) { // otherwise 'constraint hit' or whatever
      
      const int index = BPmoduleselector.getParameterIndexFromDetId(hit.det()->geographicalId(),
                                                                    eventInfo.eventId_.run());
      if (index >= 0) { // otherwise not treated
        edm::ESHandle<MagneticField> magneticField;
        setup.get<IdealMagneticFieldRecord>().get(magneticField);
        const GlobalVector bField(magneticField->inTesla(hit.det()->surface().position()));
        const LocalVector bFieldLocal(hit.det()->surface().toLocal(bField));
        //std::cout << "SiStripBackplaneCalibration derivatives " << readoutModeName_ << std::endl;
        const double dZ = hit.det()->surface().bounds().thickness(); // it's a float only...
	const double tanPsi = tsos.localParameters().mixedFormatVector()[1]; //float...

	edm::ESHandle<SiStripLorentzAngle> lorentzAngleHandle;
	setup.get<SiStripLorentzAngleRcd>().get(readoutModeName_, lorentzAngleHandle);
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
      edm::LogWarning("Alignment") << "@SUB=SiStripBackplaneCalibration::derivatives1"
                                   << "Hit without GeomDet, skip!";
    }
  } else if (mode != kDeconvolutionMode && mode != kPeakMode) {
    // warn only if unknown/mixed mode  
    edm::LogWarning("Alignment") << "@SUB=SiStripBackplaneCalibration::derivatives2"
                                 << "Readout mode is " << mode << ", but looking for "
                                 << readoutMode_ << " (" << readoutModeName_ << ").";
  }
  
  return outDerivInds.size();
}

//======================================================================
bool SiStripBackplaneCalibration::setParameter(unsigned int index, double value)
{
  if (index >= parameters_.size()) {
    return false;
  } else {
    parameters_.at(index) = value;
    return true;
  }
}

//======================================================================
bool SiStripBackplaneCalibration::setParameterError(unsigned int index, double error)
{
  if (index >= paramUncertainties_.size()) {
    return false;
  } else {
    paramUncertainties_.at(index) = error;
    return true;
  }
}

//======================================================================
double SiStripBackplaneCalibration::getParameter(unsigned int index) const
{
  //   if (index >= parameters_.size()) {
  //     return 0.;
  //   } else {
  //     return parameters_.at(index);
  //   }
  return (index >= parameters_.size() ? 0. : parameters_.at(index));
}

//======================================================================
double SiStripBackplaneCalibration::getParameterError(unsigned int index) const
{
  //   if (index >= paramUncertainties_.size()) {
  //     return 0.;
  //   } else {
  //     return paramUncertainties_[index];
  //   }
  return (index >= paramUncertainties_.size() ? 0. : paramUncertainties_.at(index));
}

// //======================================================================
// void SiStripBackplaneCalibration::beginOfJob(const AlignableTracker *tracker,
//                                                 const AlignableMuon */*muon*/,
//                                                 const AlignableExtras */*extras*/)
// {
//   alignableTracker_ = tracker;
// }


//======================================================================
void SiStripBackplaneCalibration::endOfJob()
{
  // loginfo output
  std::ostringstream out;
  out << "Parameter results for readout mode '" << readoutModeName_ << "'\n";
  for (unsigned int iPar = 0; iPar < parameters_.size(); ++iPar) {
    out << iPar << ": " << parameters_[iPar] << " +- " << paramUncertainties_[iPar] << "\n";
  }
  edm::LogInfo("Alignment") << "@SUB=SiStripBackplaneCalibration::endOfJob" << out.str();

  std::map<unsigned int, float> errors;	// Array of errors for each detId

  // now write 'input' tree
  const SiStripBackPlaneCorrection *input = this->getBackPlaneCorrectionInput(); // never NULL
  const std::string treeName(this->name() + '_' + readoutModeName_ + '_');
  this->writeTree(input, errors, (treeName + "input").c_str()); // empty errors for input...

  if (input->getBackPlaneCorrections().empty()) {
    edm::LogError("Alignment") << "@SUB=SiStripBackplaneCalibration::endOfJob"
			       << "Input back plane correction map is empty ('"
			       << readoutModeName_ << "' mode), skip writing output!";
    return;
  }

  const unsigned int nonZeroParamsOrErrors =   // Any determined value?
    count_if (parameters_.begin(), parameters_.end(), std::bind2nd(std::not_equal_to<double>(),0.))
    + count_if(paramUncertainties_.begin(), paramUncertainties_.end(),
               std::bind2nd(std::not_equal_to<double>(), 0.));

  for (unsigned int iIOV = 0; iIOV < BPmoduleselector.numIovs(); ++iIOV) {
    cond::Time_t firstRunOfIOV = BPmoduleselector.firstRunOfIOV(iIOV);
    SiStripBackPlaneCorrection *output = new SiStripBackPlaneCorrection;
    // Loop on map of values from input and add (possible) parameter results
    for (auto iterIdValue = input->getBackPlaneCorrections().begin();
	 iterIdValue != input->getBackPlaneCorrections().end(); ++iterIdValue) {
      // type of (*iterIdValue) is pair<unsigned int, float>
      const unsigned int detId = iterIdValue->first; // key of map is DetId
      const float value = iterIdValue->second + this->getParameterForDetId(detId, firstRunOfIOV);
      output->putBackPlaneCorrection(detId, value); // put result in output
      int parameterIndex = BPmoduleselector.getParameterIndexFromDetId(detId, firstRunOfIOV);
      errors[detId] = this->getParameterError(parameterIndex);
    }

    if (saveToDB_ || nonZeroParamsOrErrors != 0) { // Skip writing mille jobs...
      this->writeTree(output, errors, (treeName + Form("result_%lld", firstRunOfIOV)).c_str());
    }

    if (saveToDB_) { // If requested, write out to DB 
      edm::Service<cond::service::PoolDBOutputService> dbService;
      if (dbService.isAvailable()) {
	dbService->writeOne(output, firstRunOfIOV, recordNameDBwrite_.c_str());
	// no 'delete output;': writeOne(..) took over ownership
      } else {
	delete output;
	edm::LogError("BadConfig") << "@SUB=SiStripBackplaneCalibration::endOfJob"
				   << "No PoolDBOutputService available, but saveToDB true!";
      }
    } else {
      delete output;
    }
  } // end loop on IOVs
}

//======================================================================
bool SiStripBackplaneCalibration::checkBackPlaneCorrectionInput(const edm::EventSetup &setup,
								const EventInfo &eventInfo)
{
  edm::ESHandle<SiStripBackPlaneCorrection> backPlaneCorrHandle;
  if (!siStripBackPlaneCorrInput_) {
    setup.get<SiStripBackPlaneCorrectionRcd>().get(readoutModeName_, backPlaneCorrHandle);
    siStripBackPlaneCorrInput_ = new SiStripBackPlaneCorrection(*backPlaneCorrHandle);
    // FIXME: Should we call 'watchBackPlaneCorrRcd_.check(setup)' as well?
    //        Otherwise could be that next check has to check via following 'else', though
    //        no new IOV has started... (to be checked)
  } else {
    if (watchBackPlaneCorrRcd_.check(setup)) { // new IOV of input - but how to check peak vs deco?
      setup.get<SiStripBackPlaneCorrectionRcd>().get(readoutModeName_, backPlaneCorrHandle);
      if (backPlaneCorrHandle->getBackPlaneCorrections() // but only bad if non-identical values
	  != siStripBackPlaneCorrInput_->getBackPlaneCorrections()) { // (comparing maps)
	// Maps are containers sorted by key, but comparison problems may arise from
	// 'floating point comparison' problems (FIXME?)
	throw cms::Exception("BadInput")
	  << "SiStripBackplaneCalibration::checkBackPlaneCorrectionInput:\n"
	  << "Content of SiStripBackPlaneCorrection changed at run " << eventInfo.eventId_.run()
	  << ", but algorithm expects constant input!\n";
	return false; // not reached...
      }
    }
  }
  
  return true;
}

//======================================================================
const SiStripBackPlaneCorrection* SiStripBackplaneCalibration::getBackPlaneCorrectionInput()
{
  // For parallel processing in Millepede II, create SiStripBackPlaneCorrection
  // from info stored in files of parallel jobs and check that they are identical.
  // If this job has run on events, still check that back plane corrections are
  // identical to the ones from mergeFileNames_.
  const std::string treeName(((this->name() + '_') += readoutModeName_) += "_input");
  for (auto iFile = mergeFileNames_.begin(); iFile != mergeFileNames_.end(); ++iFile) {
    SiStripBackPlaneCorrection* bpCorr = this->createFromTree(iFile->c_str(), treeName.c_str());
    // siStripBackPlaneCorrInput_ could be non-null from previous file of this loop
    // or from checkBackPlaneCorrectionInput(..) when running on data in this job as well
    if (!siStripBackPlaneCorrInput_ || siStripBackPlaneCorrInput_->getBackPlaneCorrections().empty()) {
      delete siStripBackPlaneCorrInput_; // NULL or empty
      siStripBackPlaneCorrInput_ = bpCorr;
    } else {
      // about comparison of maps see comments in checkBackPlaneCorrectionInput
      if (bpCorr && !bpCorr->getBackPlaneCorrections().empty() && // single job might not have got events
          bpCorr->getBackPlaneCorrections() != siStripBackPlaneCorrInput_->getBackPlaneCorrections()) {
        // Throw exception instead of error?
        edm::LogError("NoInput") << "@SUB=SiStripBackplaneCalibration::getBackPlaneCorrectionInput"
                                 << "Different input values from tree " << treeName
                                 << " in file " << *iFile << ".";
        
      }
      delete bpCorr;
    }
  }

  if (!siStripBackPlaneCorrInput_) { // no files nor ran on events
    siStripBackPlaneCorrInput_ = new SiStripBackPlaneCorrection;
    edm::LogError("NoInput") << "@SUB=SiStripBackplaneCalibration::getBackPlaneCorrectionInput"
			     << "No input, create an empty one ('" << readoutModeName_ << "' mode)!";
  } else if (siStripBackPlaneCorrInput_->getBackPlaneCorrections().empty()) {
    edm::LogError("NoInput") << "@SUB=SiStripBackplaneCalibration::getBackPlaneCorrectionInput"
			     << "Empty result ('" << readoutModeName_ << "' mode)!";
  }

  return siStripBackPlaneCorrInput_;
}

//======================================================================
double SiStripBackplaneCalibration::getParameterForDetId(unsigned int detId,
							 edm::RunNumber_t run) const
{
  const int index = BPmoduleselector.getParameterIndexFromDetId(detId, run);

  return (index < 0 ? 0. : parameters_.at(index));
}

//======================================================================
void SiStripBackplaneCalibration::writeTree(const SiStripBackPlaneCorrection *backPlaneCorrection,
					    const std::map<unsigned int, float> &errors,
					    const char *treeName) const
{
  if (!backPlaneCorrection) return;

  TFile* file = TFile::Open(outFileName_.c_str(), "UPDATE");
  if (!file) {
    edm::LogError("BadConfig") << "@SUB=SiStripBackplaneCalibration::writeTree"
			       << "Could not open file '" << outFileName_ << "'.";
    return;
  }

  TTree *tree = new TTree(treeName, treeName);
  unsigned int id = 0;
  float value = 0.;
  float error = 0.;
  tree->Branch("detId", &id, "detId/i");
  tree->Branch("value", &value, "value/F");
  tree->Branch("error", &error, "error/F");

  for (auto iterIdValue = backPlaneCorrection->getBackPlaneCorrections().begin();
       iterIdValue != backPlaneCorrection->getBackPlaneCorrections().end(); ++iterIdValue) {
    // type of (*iterIdValue) is pair<unsigned int, float>
    id = iterIdValue->first; // key of map is DetId
    value = iterIdValue->second;
    auto idErrPairIter = errors.find(id); // find error for this id - if none, fill 0. in tree
    error = (idErrPairIter != errors.end() ? idErrPairIter->second : 0.f);
    tree->Fill();
  }
  tree->Write();
  delete file; // tree vanishes with the file... (?)

}

//======================================================================
SiStripBackPlaneCorrection* 
SiStripBackplaneCalibration::createFromTree(const char *fileName, const char *treeName) const
{
  // Check for file existence on your own to work around
  // https://hypernews.cern.ch/HyperNews/CMS/get/swDevelopment/2715.html:
  TFile* file = 0;
  FILE* testFile = fopen(fileName,"r");
  if (testFile) {
    fclose(testFile);
    file = TFile::Open(fileName, "READ");
  } // else not existing, see error below

  TTree *tree = 0;
  if (file) file->GetObject(treeName, tree);

  SiStripBackPlaneCorrection *result = 0;
  if (tree) {
    unsigned int id = 0;
    float value = 0.;
    tree->SetBranchAddress("detId", &id);
    tree->SetBranchAddress("value", &value);

    result = new SiStripBackPlaneCorrection;
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
      tree->GetEntry(iEntry);
      result->putBackPlaneCorrection(id, value);
    }
  } else { // Warning only since could be parallel job on no events.
    edm::LogWarning("Alignment") << "@SUB=SiStripBackplaneCalibration::createFromTree"
                                 << "Could not get TTree '" << treeName << "' from file '"
                                 << fileName << (file ? "'." : "' (file does not exist).");
  }

  delete file; // tree will vanish with file
  return result;
}


//======================================================================
//======================================================================
// Plugin definition

#include "Alignment/CommonAlignmentAlgorithm/interface/IntegratedCalibrationPluginFactory.h"

DEFINE_EDM_PLUGIN(IntegratedCalibrationPluginFactory,
		   SiStripBackplaneCalibration, "SiStripBackplaneCalibration");
