/// \class SiPixelLorentzAngleCalibration
///
/// Calibration of Lorentz angle for the pixel,
/// integrated in the alignment algorithms.
///
/// Note that not all algorithms support this...
///
///  \author    : Gero Flucke
///  date       : September 2012
///  $Revision: 1.4.2.5 $
///  $Date: 2013/04/23 08:13:27 $
///  (last update by $Author: jbehr $)

#include "Alignment/CommonAlignmentAlgorithm/interface/IntegratedCalibrationBase.h"

#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterSelector.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include <vector>
#include <map>
#include <sstream>
#include <list>
#include <utility>
#include <set>

#include <functional>

class SiPixelLorentzAngleCalibration : public IntegratedCalibrationBase
{
public:
  /// Constructor
  explicit SiPixelLorentzAngleCalibration(const edm::ParameterSet &cfg);
  
  /// Destructor
  virtual ~SiPixelLorentzAngleCalibration();

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
  /// Return 0. if index out-of-bounds.
  virtual double getParameter(unsigned int index) const;

  /// Return current value of parameter identified by index.
  /// Returns 0. if index out-of-bounds or if errors undetermined.
  virtual double getParameterError(unsigned int index) const;

  // /// Call at beginning of job:
  virtual void beginOfJob(AlignableTracker *tracker,
  			  AlignableMuon *muon,
  			  AlignableExtras *extras);



  /// Called at end of a the job of the AlignmentProducer.
  /// Write out determined parameters.
  virtual void endOfJob();

private:
  /// If called the first time, fill 'siPixelLorentzAngleInput_',
  /// later check that LorentzAngle has not changed.
  bool checkLorentzAngleInput(const edm::EventSetup &setup, const EventInfo &eventInfo);
  /// Input LorentzAngle values:
  /// - either from EventSetup of first call to derivatives(..)
  /// - or created from files of passed by configuration (i.e. from parallel processing)
  const SiPixelLorentzAngle* getLorentzAnglesInput();

  /// Determined parameter value for this detId (detId not treated => 0.)
  /// and the given run.
  double getParameterForDetId(unsigned int detId, edm::RunNumber_t run) const;
  /// Index of parameter for given detId (detId not treated => < 0)
  /// and the given run.
  int getParameterIndexFromDetId(unsigned int detId, edm::RunNumber_t run) const;
  int getParameterIndexFromDetId_old(unsigned int detId, edm::RunNumber_t run) const;


  /// Total number of IOVs.
  unsigned int numIovs() const;
  /// First run of iov (0 if iovNum not treated).
  edm::RunNumber_t firstRunOfIOV(unsigned int iovNum) const;

  void writeTree(const SiPixelLorentzAngle *lorentzAngle,
		 const std::map<unsigned int,float>& errors, const char *treeName) const;
  SiPixelLorentzAngle* createFromTree(const char *fileName, const char *treeName) const;

  const bool saveToDB_;
  const std::string recordNameDBwrite_;
  const std::string outFileName_;
  const std::vector<std::string> mergeFileNames_;

  
  const edm::VParameterSet myLAGranularityConfig_;
  
  edm::ESWatcher<SiPixelLorentzAngleRcd> watchLorentzAngleRcd_;

  // const AlignableTracker *alignableTracker_;
  SiPixelLorentzAngle *siPixelLorentzAngleInput_;
  std::vector<double> parameters_;
  std::vector<double> paramUncertainties_;
  
  std::set<edm::RunNumber_t> globalRunRange_;
  std::vector<unsigned int> firstId_;//1:1 mapping with LAassignment
  std::vector<std::pair<std::list<Alignable*>, std::vector<edm::RunNumber_t> > > LAassignment_;

 
 
};

//======================================================================
//======================================================================
//======================================================================

SiPixelLorentzAngleCalibration::SiPixelLorentzAngleCalibration(const edm::ParameterSet &cfg)
  : IntegratedCalibrationBase(cfg),
    saveToDB_(cfg.getParameter<bool>("saveToDB")),
    recordNameDBwrite_(cfg.getParameter<std::string>("recordNameDBwrite")),
    outFileName_(cfg.getParameter<std::string>("treeFile")),
    mergeFileNames_(cfg.getParameter<std::vector<std::string> >("mergeTreeFiles")),
    myLAGranularityConfig_(cfg.getParameter<edm::VParameterSet>("LorentzAngleGranularity")),
    //    alignableTracker_(0),
    siPixelLorentzAngleInput_(0)
{
  // FIXME: Which granularity, leading to how many parameters?
//  parameters_.resize(2, 0.); // currently two parameters (ring1-4, 5-8), start value 0.
//  paramUncertainties_.resize(2, 0.); // dito for errors
//  parameters_.resize(349, 0.); // BPIX:6, 58 IOVs, start value 0.
//  paramUncertainties_.resize(349, 0.); // dito for errors
//  parameters_.resize(1072, 0.); // currently two parameters BPIX:6, FPIX:2, 134 IOVs, start value 0.
//  paramUncertainties_.resize(1072, 0.); // dito for errors
  // parameters_.resize(1274, 0.); // BPIX:24, FPIX:2, 49 IOVs, start value 0.
  // paramUncertainties_.resize(1274, 0.); // dito for errors

  edm::LogInfo("Alignment") << "@SUB=SiPixelLorentzAngleCalibration" << "Created with name "
                            << this->name() << "',\n" << this->numParameters() << " parameters to be determined,"
                            << "\nsaveToDB = " << saveToDB_
                            << "\n outFileName = " << outFileName_
                            << "\n N(merge files) = " << mergeFileNames_.size();
  if (mergeFileNames_.size()) {
    edm::LogInfo("Alignment") << "@SUB=SiPixelLorentzAngleCalibration"
                              << "First file to merge: " << mergeFileNames_[0];
  }
}
  
//======================================================================
SiPixelLorentzAngleCalibration::~SiPixelLorentzAngleCalibration()
{
  //  std::cout << "Destroy SiPixelLorentzAngleCalibration named " << this->name() << std::endl;
  delete siPixelLorentzAngleInput_;
}

//======================================================================
unsigned int SiPixelLorentzAngleCalibration::numParameters() const
{
  return parameters_.size();
}

//======================================================================
unsigned int
SiPixelLorentzAngleCalibration::derivatives(std::vector<ValuesIndexPair> &outDerivInds,
					    const TransientTrackingRecHit &hit,
					    const TrajectoryStateOnSurface &tsos,
					    const edm::EventSetup &setup,
					    const EventInfo &eventInfo) const
{
  // ugly const-cast:
  // But it is either only first initialisation or throwing an exception...
  const_cast<SiPixelLorentzAngleCalibration*>(this)->checkLorentzAngleInput(setup, eventInfo);

  outDerivInds.clear();

  if (hit.det()) { // otherwise 'constraint hit' or whatever
    
    const int index = this->getParameterIndexFromDetId(hit.det()->geographicalId(),
						       eventInfo.eventId_.run());
    if (index >= 0) { // otherwise not treated
      edm::ESHandle<MagneticField> magneticField;
      setup.get<IdealMagneticFieldRecord>().get(magneticField);
      const GlobalVector bField(magneticField->inTesla(hit.det()->surface().position()));
      const LocalVector bFieldLocal(hit.det()->surface().toLocal(bField));
      const double dZ = hit.det()->surface().bounds().thickness(); // it is a float only...
      // shift due to LA: dx = tan(LA) * dz/2 = mobility * B_y * dz/2,
      // '-' since we have derivative of the residual r = trk -hit
      const double xDerivative = bFieldLocal.y() * dZ * -0.5; // parameter is mobility!

//      const PXBDetId bdet(hit.det()->geographicalId());
//      if(bdet.subdetId()==PixelSubdetector::PixelBarrel){
//	printf("TPB: layer: %d ring: %d B: %.5f\n",bdet.layer(),bdet.module(),bFieldLocal.y());
//      }
//      const PXFDetId fdet(hit.det()->geographicalId());
//      if(fdet.subdetId()==PixelSubdetector::PixelEndcap){
//	printf("TPF: side: %d B: %.5f\n",fdet.side(),bFieldLocal.y());
//      }

      if (xDerivative) { // If field is zero, this is zero: do not return it
	const Values derivs(xDerivative, 0.); // yDerivative = 0.
	outDerivInds.push_back(ValuesIndexPair(derivs, index));
      }
    }
  } else {
    edm::LogWarning("Alignment") << "@SUB=SiPixelLorentzAngleCalibration::derivatives2"
				 << "Hit without GeomDet, skip!";
  }
  
  return outDerivInds.size();
}

//======================================================================
bool SiPixelLorentzAngleCalibration::setParameter(unsigned int index, double value)
{
  if (index >= parameters_.size()) {
    return false;
  } else {
    parameters_.at(index) = value;
    return true;
  }
}

//======================================================================
bool SiPixelLorentzAngleCalibration::setParameterError(unsigned int index, double error)
{
  if (index >= paramUncertainties_.size()) {
    return false;
  } else {
    paramUncertainties_.at(index) = error;
    return true;
  }
}

//======================================================================
double SiPixelLorentzAngleCalibration::getParameter(unsigned int index) const
{
  //   if (index >= parameters_.size()) {
  //     return 0.;
  //   } else {
  //     return parameters_[index];
  //   }
  return (index >= parameters_.size() ? 0. : parameters_.at(index));
}

//======================================================================
double SiPixelLorentzAngleCalibration::getParameterError(unsigned int index) const
{
  //   if (index >= paramUncertainties_.size()) {
  //     return 0.;
  //   } else {
  //     return paramUncertainties_[index];
  //   }
  return (index >= paramUncertainties_.size() ? 0. : paramUncertainties_.at(index));
}




//======================================================================
void SiPixelLorentzAngleCalibration::beginOfJob(AlignableTracker *aliTracker,
                                                AlignableMuon *aliMuon,
                                                AlignableExtras *aliExtras)
{
  unsigned int nparameters = 0;
  unsigned int Id = 0;
  //loop over all LA groups
  for(edm::VParameterSet::const_iterator pset = myLAGranularityConfig_.begin();
      pset != myLAGranularityConfig_.end();
      ++pset) {
    const std::vector<edm::RunNumber_t> range = pset->getParameter<std::vector<edm::RunNumber_t> >("RunRange");

    bool split = false;
    if((*pset).exists("split")) {
      split = pset->getParameter<bool>("split");
    }
    const size_t npar = (*pset).getParameterNames().size();
    
    if (npar >= 3 && !(*pset).exists("split")) {
      throw cms::Exception("BadConfig")
        << " >= 3 parameters specified in PSet BUT split parameter was not found! Maybe a typo?";
    }
    
    AlignmentParameterSelector selector(aliTracker,aliMuon, aliExtras);
    selector.clear();
    selector.addSelections((*pset).getParameter<edm::ParameterSet> ("levels"));

    const std::vector<Alignable*> &alis = selector.selectedAlignables();
    

    std::list<Alignable*> selected_alis;
    for(std::vector<Alignable*>::const_iterator it = alis.begin(); it != alis.end(); it++) {
      if((*it)->alignableObjectId() == align::AlignableDetUnit || (*it)->alignableObjectId() == align::AlignableDet) {
        if(split) {
          LAassignment_.push_back(std::make_pair(std::list<Alignable*>(1,(*it)), range));
          firstId_.push_back(Id);
          Id += range.size();
          nparameters += range.size();
        } else {
          selected_alis.push_back((*it)); //throw out HLS?
        }
      }
      const std::vector<Alignable*> &aliDaughts = (*it)->deepComponents();
      if(aliDaughts.size() > 0) {
        for (std::vector<Alignable*>::const_iterator iD = aliDaughts.begin();
             iD != aliDaughts.end(); ++iD) {
          
          if((*iD)->alignableObjectId() == align::AlignableDetUnit || (*iD)->alignableObjectId() == align::AlignableDet) {
            bool found = false;
            for(unsigned int iAlignableGroup = 0; iAlignableGroup < firstId_.size(); iAlignableGroup++) {
              const std::list<Alignable*> &a = LAassignment_.at(iAlignableGroup).first;
              
              for(std::list<Alignable*>::const_iterator iAli = a.begin();
                  iAli != a.end(); iAli++) {
                if( (*iAli)->alignableObjectId() == (*iD)->alignableObjectId()
                    && (*iAli)->id() == (*iD)->id()) {
                  found = true;
                  if(!split) {
                    throw cms::Exception("BadConfig")
                      << " Alignable already selected.";
                  }
                  break;
                }
              }
            }
            if(split) {
              if(!found) {
                LAassignment_.push_back(std::make_pair(std::list<Alignable*>(1,(*iD)), range));
                firstId_.push_back(Id);
                Id += range.size();
                nparameters += range.size();
              }
            } else {
              selected_alis.push_back((*iD));
            }
            
          



          }
        }
      }

    }
 
   
    //FIXME: add some checks whether the content of range makes sense?



  
    if(!split) {
      LAassignment_.push_back(std::make_pair(selected_alis, range));
      firstId_.push_back(Id);
      Id += range.size();
      nparameters += range.size();
    }
    
    edm::RunNumber_t firstRun = 0; 
    for(std::vector<edm::RunNumber_t>::const_iterator iRun = range.begin(); 
        iRun != range.end(); ++iRun)  {
      globalRunRange_.insert((*iRun));
      if((*iRun) > firstRun) {
        firstRun = (*iRun);
      } else {
        throw cms::Exception("BadConfig")
          << "@SUB=SiPixelLorentzAngleCalibration::beginOfJob:"
          << " Run range not sorted.";
      }
    }
  }
  

  parameters_.resize(nparameters, 0.);
  paramUncertainties_.resize(nparameters, 0.);
  
 
  //test that only pixel modules have been selected
  for(unsigned int iAlignableGroup = 0; iAlignableGroup < firstId_.size(); iAlignableGroup++) {
    const std::list<Alignable*> &alis = LAassignment_.at(iAlignableGroup).first;
    
    for(std::list<Alignable*>::const_iterator iAli = alis.begin();
        iAli != alis.end(); iAli++) {
      
      bool pixelmodule = true;
      if((*iAli)->alignableObjectId() == align::AlignableDetUnit || (*iAli)->alignableObjectId() == align::AlignableDet) {
        const PXBDetId id((*iAli)->id());
        if (id.det() != DetId::Tracker) pixelmodule = false;
        if (id.det() == DetId::Tracker && ( id.subdetId() == PixelSubdetector::PixelBarrel || id.subdetId() == PixelSubdetector::PixelEndcap)) {
          //nothing
        } else {
          pixelmodule = false;
        }
      } else {
        pixelmodule = false;
      }

      if(!pixelmodule) {
        throw cms::Exception("BadConfig") 
          << " Non-pixel modules have been selected for the LA determination in the pixel detector";
      }
    }
  }
  
 
}


//======================================================================
void SiPixelLorentzAngleCalibration::endOfJob()
{
  // loginfo output
  std::ostringstream out;
  out << "Parameter results\n";
  for (unsigned int iPar = 0; iPar < parameters_.size(); ++iPar) {
    out << iPar << ": " << parameters_[iPar] << " +- " << paramUncertainties_[iPar] << "\n";
  }
  edm::LogInfo("Alignment") << "@SUB=SiPixelLorentzAngleCalibration::endOfJob" << out.str();

  std::map<unsigned int,float> errors;	  // Array of errors for each detId

  // now write 'input' tree
  const SiPixelLorentzAngle *input = this->getLorentzAnglesInput(); // never NULL
  const std::string treeName(this->name() + '_');
  this->writeTree(input, errors, (treeName + "input").c_str()); // empty errors for input...

  if (input->getLorentzAngles().empty()) {
    edm::LogError("Alignment") << "@SUB=SiPixelLorentzAngleCalibration::endOfJob"
			       << "Input Lorentz angle map is empty, skip writing output!";
    return;
  }

  const unsigned int nonZeroParamsOrErrors =   // Any determined value?
    count_if (parameters_.begin(), parameters_.end(), std::bind2nd(std::not_equal_to<double>(),0.))
    + count_if(paramUncertainties_.begin(), paramUncertainties_.end(),
               std::bind2nd(std::not_equal_to<double>(), 0.));

  for (unsigned int iIOV = 0; iIOV < this->numIovs(); ++iIOV) {
//  for (unsigned int iIOV = 0; iIOV < 1; ++iIOV) {   // For writing out the modified values
    cond::Time_t firstRunOfIOV = this->firstRunOfIOV(iIOV);
    SiPixelLorentzAngle *output = new SiPixelLorentzAngle;
    // Loop on map of values from input and add (possible) parameter results
    for (auto iterIdValue = input->getLorentzAngles().begin();
	 iterIdValue != input->getLorentzAngles().end(); ++iterIdValue) {
      // type of (*iterIdValue) is pair<unsigned int, float>
      const unsigned int detId = iterIdValue->first; // key of map is DetId
      // Nasty: putLorentzAngle(..) takes float by reference - not even const reference!
      float value = iterIdValue->second + this->getParameterForDetId(detId, firstRunOfIOV);
//      float value = iterIdValue->second + this->getParameterForDetId(detId, firstRunOfIOV) + 0.02;  // Added 0.02 for LA in BPIX study
      output->putLorentzAngle(detId, value); // put result in output
      int parameterIndex = this->getParameterIndexFromDetId(detId, firstRunOfIOV);
      errors[detId]= this->getParameterError(parameterIndex);
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
	edm::LogError("BadConfig") << "@SUB=SiPixelLorentzAngleCalibration::endOfJob"
				   << "No PoolDBOutputService available, but saveToDB true!";
      }
    } else {
      delete output;
    }
  } // end loop on IOVs

}

//======================================================================
bool SiPixelLorentzAngleCalibration::checkLorentzAngleInput(const edm::EventSetup &setup,
							    const EventInfo &eventInfo)
{
  edm::ESHandle<SiPixelLorentzAngle> lorentzAngleHandle;
  if (!siPixelLorentzAngleInput_) {
    setup.get<SiPixelLorentzAngleRcd>().get(lorentzAngleHandle);
    siPixelLorentzAngleInput_ = new SiPixelLorentzAngle(*lorentzAngleHandle);
  } else {
    if (watchLorentzAngleRcd_.check(setup)) { // new IOV of input
      setup.get<SiPixelLorentzAngleRcd>().get(lorentzAngleHandle);
      if (lorentzAngleHandle->getLorentzAngles() // but only bad if non-identical values
	  != siPixelLorentzAngleInput_->getLorentzAngles()) { // (comparing maps)
	// Maps are containers sorted by key, but comparison problems may arise from
	// 'floating point comparison' problems (FIXME?)
	throw cms::Exception("BadInput")
	  << "SiPixelLorentzAngleCalibration::checkLorentzAngleInput:\n"
	  << "Content of SiPixelLorentzAngle changed at run " << eventInfo.eventId_.run()
	  << ", but algorithm expects constant input!\n";
	return false; // not reached...
      }
    }
  }
  
  return true;
}

//======================================================================
const SiPixelLorentzAngle* SiPixelLorentzAngleCalibration::getLorentzAnglesInput()
{
  // For parallel processing in Millepede II, create SiPixelLorentzAngle
  // from info stored in files of parallel jobs and check that they are identical.
  // If this job has run on data, still check that LA is identical to the ones
  // from mergeFileNames_.
  const std::string treeName(this->name() + "_input");
  for (auto iFile = mergeFileNames_.begin(); iFile != mergeFileNames_.end(); ++iFile) {
    SiPixelLorentzAngle* la = this->createFromTree(iFile->c_str(), treeName.c_str());
    // siPixelLorentzAngleInput_ could be non-null from previous file of this loop
    // or from checkLorentzAngleInput(..) when running on data in this job as well
    if (!siPixelLorentzAngleInput_ || siPixelLorentzAngleInput_->getLorentzAngles().empty()) {
      delete siPixelLorentzAngleInput_; // NULL or empty
      siPixelLorentzAngleInput_ = la;
    } else {
      // FIXME: about comparison of maps see comments in checkLorentzAngleInput
      if (la && !la->getLorentzAngles().empty() && // single job might not have got events
          la->getLorentzAngles() != siPixelLorentzAngleInput_->getLorentzAngles()) {
        // Throw exception instead of error?
        edm::LogError("NoInput") << "@SUB=SiPixelLorentzAngleCalibration::getLorentzAnglesInput"
                                 << "Different input values from tree " << treeName
                                 << " in file " << *iFile << ".";
        
      }
      delete la;
    }
  }

  if (!siPixelLorentzAngleInput_) { // no files nor ran on events
    siPixelLorentzAngleInput_ = new SiPixelLorentzAngle;
    edm::LogError("NoInput") << "@SUB=SiPixelLorentzAngleCalibration::getLorentzAnglesInput"
                             << "No input, create an empty one!";
  } else if (siPixelLorentzAngleInput_->getLorentzAngles().empty()) {
    edm::LogError("NoInput") << "@SUB=SiPixelLorentzAngleCalibration::getLorentzAnglesInput"
                             << "Empty result!";
  }

  return siPixelLorentzAngleInput_;
}

//======================================================================
double SiPixelLorentzAngleCalibration::getParameterForDetId(unsigned int detId,
							    edm::RunNumber_t run) const
{
  const int index = this->getParameterIndexFromDetId(detId, run);
//  const PXBDetId id(detId);
//  if (id.det() == DetId::Tracker && ( id.subdetId() == PixelSubdetector::PixelBarrel || id.subdetId() == PixelSubdetector::PixelEndcap)) {
//    edm::LogInfo("Alignment") << "@SUB=SiPixelLorentzAngleCalibration::getParameterForDetId" 
//      << "SubDetId: " << id.subdetId()
//      << " Layer: " << id.layer() 
//      << " Module: " << id.module() 
//      << " Run: " << run
//      << " Index: " << index;
//  }
  return (index < 0 ? 0. : parameters_.at(index));
}

//======================================================================
int SiPixelLorentzAngleCalibration::getParameterIndexFromDetId(unsigned int detId,
							       edm::RunNumber_t run) const
{
  // Return the index of the parameter that is used for this DetId.
  // If this DetId is not treated, return values < 0.
  
  const PXBDetId temp_id(detId);
  // const unsigned int nLayers=3;
  // const unsigned int nRings=8;
  int index = -1;
  if (temp_id.det() != DetId::Tracker || ( temp_id.subdetId() != PixelSubdetector::PixelBarrel && temp_id.subdetId() != PixelSubdetector::PixelEndcap)) return -1;
  

  int nmatched = 0;
  bool matched = false;
  for(unsigned int iAlignableGroup = 0; iAlignableGroup < firstId_.size(); iAlignableGroup++) {
    const unsigned int id0 = firstId_.at(iAlignableGroup);
    
    const std::list<Alignable*> &alis = LAassignment_.at(iAlignableGroup).first;
    const std::vector<edm::RunNumber_t> &runs = LAassignment_.at(iAlignableGroup).second;
    
  
    for(std::list<Alignable*>::const_iterator iAli = alis.begin();
        iAli != alis.end(); iAli++) {
      if( ((*iAli)->alignableObjectId() == align::AlignableDetUnit
           || (*iAli)->alignableObjectId() == align::AlignableDet)
          && (*iAli)->id() == detId) {
        nmatched++;
        break;
      }
    }
    if(!matched && nmatched == 1) {
      matched = true;
      int iovNum = -1;
      for(std::vector<edm::RunNumber_t>::const_iterator itRun = runs.begin();
          itRun != runs.end(); itRun++) {
        const edm::RunNumber_t r = (*itRun);
        if(run >= r) iovNum++;
      }

      if(iovNum == -1) {
        throw cms::Exception("BadRunRange")
	  << "@SUB=SiPixelLorentzAngleCalibration::getParameterIndexFromDetId:\n"
	  << "Bad run range for detid ('"<< detId <<"').";
      }
      index = id0 + iovNum;
      // const PXBDetId temp_id(detId);
      // const unsigned int nLayers=3;
      // const unsigned int nRings=8;
      // int index_old = -1;
      // if(temp_id.subdetId() == PixelSubdetector::PixelBarrel) {
      //   const PXBDetId id(detId);
      //   index_old = iovNum*(nLayers*nRings+2)+(id.layer()-1)*(nRings)+(id.module()-1);
      // } else if(temp_id.subdetId() == PixelSubdetector::PixelEndcap) { 
      //   const PXFDetId id(detId);
      //   index_old = iovNum*(nLayers*nRings+2)+nLayers*nRings+(id.side()-1);
      // }

      // const int index_old = getParameterIndexFromDetId_old(detId,run);
      // std::cout << "debug indices " << index << " " << index_old << std::endl;


    }
  }

  if(nmatched >= 2) {
    throw cms::Exception("BadConfig")
      << "@SUB=SiPixelLorentzAngleCalibration::getParameterIndexFromDetId:\n"
      << "Detid ('"<< detId <<"') is assigned to >= 2 groups of modules (in fact it is matched to '"<< nmatched <<"' groups).";
  }
  

  return index;
}

//======================================================================
int SiPixelLorentzAngleCalibration::getParameterIndexFromDetId_old(unsigned int detId,
                                                                   edm::RunNumber_t run) const
{
  //int index = -1;
  // Return the index of the parameter that is used for this DetId.
  // If this DetId is not treated, return values < 0.
  
  // FIXME: Extend to configurable granularity? 
  //        Including treatment of run dependence?
  const PXBDetId temp_id(detId);
  const unsigned int nLayers=3;
  const unsigned int nRings=8;
  if (temp_id.det() != DetId::Tracker || ( temp_id.subdetId() != PixelSubdetector::PixelBarrel && temp_id.subdetId() != PixelSubdetector::PixelEndcap)) return -1;


  int iovNum=-1;
  for(unsigned int iov=0; iov<this->numIovs(); iov++) {
    if(run >= this->firstRunOfIOV(iov)) iovNum=iov;
  }
  if(iovNum<0) {
    edm::LogWarning("Alignment")
      << "@SUB=SiPixelLorentzAngleCalibration::getParameterIndexFromDetId"
      << "IOV not determined for current run: " << run << " => skip!";
    return -1;
  }

  if(temp_id.subdetId() == PixelSubdetector::PixelBarrel) {
    const PXBDetId id(detId);

    if(id.layer()<1 || id.layer()>nLayers) {
      edm::LogWarning("Alignment")
	<< "@SUB=SiPixelLorentzAngleCalibration::getParameterIndexFromDetId"
	<< "Layer should be 1-3, but is " << id.layer() << "=> skip!";
      return -1;
    }

//    if (id.module() >= 1 && id.module() <= 4) return iovNum*nLayers*2+(id.layer()-1)*2+0; else    // Just BPIX 3layers * 2Zhalves
//    if (id.module() >= 5 && id.module() <= 8) return iovNum*nLayers*2+(id.layer()-1)*2+1; else

//    if (id.module() >= 1 && id.module() <= 4) return iovNum*(nLayers*2+2)+(id.layer()-1)*2+0; else  // BPIX 3layers * 2Zhalves and FPIX 2sides
//    if (id.module() >= 5 && id.module() <= 8) return iovNum*(nLayers*2+2)+(id.layer()-1)*2+1; else

    if(id.module() >= 1 && id.module() <= 8)return iovNum*(nLayers*nRings+2)+(id.layer()-1)*(nRings)+(id.module()-1); // BPIX 3layers * 8Rings and FPIX 2sides

    edm::LogWarning("Alignment")
      << "@SUB=SiPixelLorentzAngleCalibration::getParameterIndexFromDetId"
      << "Module should be 1-8, but is " << id.module() << " => skip!";
//  } else if(id.subdetId() == PixelSubdetector::PixelEndcap) return 348;   // Just BPIX 3layers * 2Zhalves

//  } else if(temp_id.subdetId() == PixelSubdetector::PixelEndcap) {  // BPIX 3layers * 2Zhalves and FPIX 2sides
//    const PXFDetId id(detId);
//    return iovNum*(nLayers*2+2)+(nLayers-1)*2+2+id.side();
//  }

  } else if(temp_id.subdetId() == PixelSubdetector::PixelEndcap) {  // BPIX 3layers * 8Rings and FPIX 2sides
    const PXFDetId id(detId);
    return iovNum*(nLayers*nRings+2)+nLayers*nRings+(id.side()-1);
  }

  return -1;

}

//======================================================================
unsigned int SiPixelLorentzAngleCalibration::numIovs() const
{
  // FIXME: Needed to include treatment of run dependence!
//  return 1; 
//  return 58; 
//  return 134;
//  return 119;	  // 100 pb-1
  return globalRunRange_.size();
}

//======================================================================
edm::RunNumber_t SiPixelLorentzAngleCalibration::firstRunOfIOV(unsigned int iovNum) const
{
  // FIXME: Needed to include treatment of run dependence!
//  if (iovNum < this->numIovs()) return 1;
//  else return 0;
//  unsigned int runNumbers[6] = {1,189147,190782,191718,193093,194896}; // 6 IOVs
//  unsigned int runNumbers[58] = {1,190707,190895,191086,191226,191271,191700,191811,191845,193336,193621,194050,194108,194120,194199,194225,194315,194428,194455,194480,194627,194691,194711,194789,194912,195016,195113,195163,195266,195304,195378,195397,195399,195540,195552,195655,195757,195774,195915,195937,195950,196096,196197,196218,196250,196353,196364,196433,196438,196452,196453,196531,198212,198230,198271,198487,198955,198969};    // 58 IOVs
//  unsigned int runNumbers[134] = {1,190707,190895,191086,191226,191271,191700,191811,191845,193336,193621,194050,194108,194120,194199,194225,194315,194428,194455,194480,194627,194691,194711,194789,194912,195016,195113,195163,195266,195304,195378,195397,195399,195540,195552,195655,195757,195774,195915,195937,195950,196096,196197,196218,196250,196353,196364,196433,196438,196452,196453,196531,198212,198230,198271,198487,198955,198969,199021,199319,199409,199435,199569,199608,199703,199754,199812,199833,199864,199876,199960,200041,200075,200091,200190,200244,200369,200473,200519,200525,200600,200786,200991,201097,201168,201191,201278,201410,201424,201602,201624,201625,201671,201707,201802,201824,202016,202045,202060,202084,202087,202178,202237,202272,202304,202328,202472,202478,202504,202973,203002,203894,203912,203987,204113,204250,204541,204553,204563,204564,204577,204601,205158,205215,205238,205311,205344,205595,205620,205667,205694,205774,205826,205921};  // 119 IOVs
//  unsigned int runNumbers[119] = {1,190738,191056,191226,191271,191721,191834,193336,193621,194050,194108,194150,194210,194314,194424,194429,194479,194533,194691,194711,194789,194912,195016,195113,195163,195266,195304,195378,195397,195399,195540,195552,195655,195758,195774,195915,195947,195950,196197,196218,196250,196362,196364,196437,196452,196453,198063,198230,198271,198487,198955,198969,199021,199319,199409,199435,199571,199608,199703,199754,199812,199833,199864,199876,199960,200041,200075,200091,200190,200244,200381,200473,200519,200532,200600,200991,201097,201168,201191,201278,201613,201657,201671,201707,201802,202013,202016,202045,202060,202084,202087,202178,202237,202272,202314,202389,202478,202504,202973,203002,203894,203912,203987,204113,204250,204554,204564,204577,204601,205158,205217,205238,205311,205344,205617,205667,205718,205826,205921}; // 119 IOVs
  // unsigned int runNumbers[49] = {1,191226,191834,193998,194115,194315,194480,194711,194912,195147,195304,195398,195645,195775,195948,196203,196349,196438,196531,198272,198955,199021,199409,199569,199699,199832,199876,200049,200190,200473,200600,201164,201196,201611,201671,201817,202054,202088,202272,202477,202972,203894,204101,204554,204577,205193,205339,205666,205781}; // 49 IOVs
  // if (iovNum < this->numIovs()) return runNumbers[iovNum];
  // else return 0;

  

  edm::RunNumber_t r = 0;

  for(std::set<edm::RunNumber_t>::const_iterator itRun = globalRunRange_.begin();
      itRun != globalRunRange_.end(); itRun++) {

    if(std::distance(globalRunRange_.begin(), itRun ) == iovNum) {
      r = (*itRun);
          
    }
  }
      
  return r;

}


//======================================================================
void SiPixelLorentzAngleCalibration::writeTree(const SiPixelLorentzAngle *lorentzAngle,
					       const std::map<unsigned int,float> &errors, 
					       const char *treeName) const
{
  if (!lorentzAngle) return;

  TFile* file = TFile::Open(outFileName_.c_str(), "UPDATE");
  if (!file) {
    edm::LogError("BadConfig") << "@SUB=SiPixelLorentzAngleCalibration::writeTree"
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

  for (auto iterIdValue = lorentzAngle->getLorentzAngles().begin();
       iterIdValue != lorentzAngle->getLorentzAngles().end(); ++iterIdValue) {
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
SiPixelLorentzAngle* 
SiPixelLorentzAngleCalibration::createFromTree(const char *fileName, const char *treeName) const
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

  SiPixelLorentzAngle *result = 0;
  if (tree) {
    unsigned int id = 0;
    float value = 0.;
    tree->SetBranchAddress("detId", &id);
    tree->SetBranchAddress("value", &value);

    result = new SiPixelLorentzAngle;
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
      tree->GetEntry(iEntry);
      result->putLorentzAngle(id, value);
    }
  } else { // Warning only since could be parallel job on no events.
    edm::LogWarning("Alignment") << "@SUB=SiPixelLorentzAngleCalibration::createFromTree"
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
		   SiPixelLorentzAngleCalibration, "SiPixelLorentzAngleCalibration");
