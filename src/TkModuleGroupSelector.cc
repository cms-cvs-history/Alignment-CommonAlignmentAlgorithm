/**
 * \file TkModuleGroupSelector.cc
 *
 *  \author Joerg Behr
 *  \date May 2013
 *  $Revision: 1.1.2.7 $
 *  $Date: 2013/05/17 15:09:21 $
 *  (last update by $Author: jbehr $)
 */

#include "Alignment/CommonAlignmentAlgorithm/interface/TkModuleGroupSelector.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterSelector.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <map>
#include <set>

//============================================================================
TkModuleGroupSelector::TkModuleGroupSelector(const edm::VParameterSet &cfg) : myGranularityConfig_(cfg),
                                                                              nparameters_(0),
                                                                              globalReferenceRun_(0)
{
  //FIXME: take list of subdetids from which modules are taken
}

//============================================================================
void TkModuleGroupSelector::fillDetIdMap(const unsigned int detid, const unsigned int groupid)
{
  //only add new entries 
  if(mapDetIdGroupId_.find(detid) == mapDetIdGroupId_.end()) {
    mapDetIdGroupId_.insert(std::pair<unsigned int, unsigned int>(detid, groupid));
  } else {
    throw cms::Exception("BadConfig")
      << "@SUB=TkModuleGroupSelector:fillDetIdMap:"
      << " Module with det ID " << detid << " already selected.";
  }
}

//============================================================================
void TkModuleGroupSelector::setSubDets(const std::vector<int> &sdets)
{
  subdetids_ = sdets;
}

//============================================================================
void TkModuleGroupSelector::setReferenceRun(const edm::ParameterSet &cfg)
{
  //extract the reference run range if defined
  if(cfg.exists("ReferenceRun")) {
    globalReferenceRun_ = cfg.getParameter<edm::RunNumber_t>("ReferenceRun");
  }
}

//============================================================================
const bool TkModuleGroupSelector::testSplitOption(const edm::ParameterSet &pset) const
{
  bool split = false;
  if(pset.exists("split")) {
    split = pset.getParameter<bool>("split");
  }
  return split;
}

//============================================================================
bool TkModuleGroupSelector::createGroup(
                                        const bool split, 
                                        unsigned int &Id, 
                                        const std::vector<edm::RunNumber_t> &range, 
                                        Alignable* iD, 
                                        const std::list<Alignable*> &selected_alis,
                                        const edm::RunNumber_t refrun
                                        )
{
  bool modules_selected = false;

  if(iD != NULL && selected_alis.size() == 0) {
    if(split) {
      referenceRun_.push_back(refrun);
      firstId_.push_back(Id);
      runRange_.push_back(range);
      this->fillDetIdMap(iD->id(), firstId_.size()-1);
      modules_selected = true;
      Id += range.size();
      nparameters_ += range.size();
    }
  } else {
    //iD == NULL
    if(!split) {
      referenceRun_.push_back(refrun);
      firstId_.push_back(Id);
      runRange_.push_back(range);
      for(std::list<Alignable*>::const_iterator it = selected_alis.begin();
          it != selected_alis.end(); it++) {
        this->fillDetIdMap((*it)->id(), firstId_.size()-1);
        modules_selected = true;
      }
      Id += range.size();
      nparameters_ += range.size();
    }
  }
  return modules_selected;
}

//============================================================================
void TkModuleGroupSelector::verifyParameterNames(const edm::ParameterSet &pset, unsigned int psetnr) const
{
  std::vector<std::string> parameterNames = pset.getParameterNames();
  for ( std::vector<std::string>::const_iterator iParam = parameterNames.begin(); 
        iParam != parameterNames.end(); ++iParam) {
    const std::string name = (*iParam);
    if(
       name != "levels" && name != "RunRange"
       && name != "split" && name != "ReferenceRun"
       ) {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector::verifyParameterNames:"
        << " Unknown parameter name '" << name << "' in PSet number " << psetnr << ". Maybe a typo?";
    }
  }
}


//============================================================================
void TkModuleGroupSelector::createModuleGroups(AlignableTracker *aliTracker,
                                               AlignableMuon *aliMuon,
                                               AlignableExtras *aliExtras)
{
  std::set<edm::RunNumber_t> localRunRange;
  nparameters_ = 0;
  unsigned int Id = 0;
  unsigned int psetnr = 0;
  //loop over all LA groups
  for(edm::VParameterSet::const_iterator pset = myGranularityConfig_.begin();
      pset != myGranularityConfig_.end();
      ++pset) {

    //test for unknown parameters
    this->verifyParameterNames((*pset),psetnr);
    psetnr++;

    bool modules_selected = false; //track whether at all a module has been selected in this group
    const std::vector<edm::RunNumber_t> range = pset->getParameter<std::vector<edm::RunNumber_t> >("RunRange");
    if(range.size() == 0) {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector::createModuleGroups:\n"
        << "Run range array empty!";
    }
    const bool split = this->testSplitOption((*pset));

    edm::RunNumber_t refrun = 0;
    if((*pset).exists("ReferenceRun")) {
      refrun = (*pset).getParameter<edm::RunNumber_t>("ReferenceRun");
    } else {
      refrun = globalReferenceRun_;
    }


    AlignmentParameterSelector selector(aliTracker,aliMuon, aliExtras);
    selector.clear();
    selector.addSelections((*pset).getParameter<edm::ParameterSet> ("levels"));

    const std::vector<Alignable*> &alis = selector.selectedAlignables();
    std::list<Alignable*> selected_alis;
    for(std::vector<Alignable*>::const_iterator it = alis.begin(); it != alis.end(); it++) {
      const std::vector<Alignable*> &aliDaughts = (*it)->deepComponents();
      for (std::vector<Alignable*>::const_iterator iD = aliDaughts.begin();
           iD != aliDaughts.end(); ++iD) {
        if((*iD)->alignableObjectId() == align::AlignableDetUnit || (*iD)->alignableObjectId() == align::AlignableDet) {
          if(split) {
            modules_selected = this->createGroup(split, Id, range, (*iD), std::list<Alignable*>(), refrun);//last parameter is a empty dummy list
          } else {
            selected_alis.push_back((*iD));
          }
        }
      }
    }
    
    if(!split) {
      modules_selected = this->createGroup(split, Id, range, NULL, selected_alis,refrun);
    }
        
    edm::RunNumber_t firstRun = 0; 
    for(std::vector<edm::RunNumber_t>::const_iterator iRun = range.begin(); 
        iRun != range.end(); ++iRun)  {
      localRunRange.insert((*iRun));
      if((*iRun) > firstRun) {
        firstRun = (*iRun);
      } else {
        throw cms::Exception("BadConfig")
          << "@SUB=TkModuleGroupSelector::createModuleGroups:"
          << " Run range not sorted.";
      }
    }

    if(!modules_selected) {
      throw cms::Exception("BadConfig") 
        << "@SUB=TkModuleGroupSelector:createModuleGroups:"
        << " No module was selected in the module group selector in group " << (firstId_.size()-1)<< ".";
    }
  }

  //copy local set into the global vector of run boundaries
  for(std::set<edm::RunNumber_t>::const_iterator itRun = localRunRange.begin();
      itRun != localRunRange.end(); itRun++) {
    globalRunRange_.push_back((*itRun));
  }
}

//============================================================================
unsigned int TkModuleGroupSelector::getNumberOfParameters() const
{
  return nparameters_;
}

//============================================================================
unsigned int TkModuleGroupSelector::numIovs() const
{
  return globalRunRange_.size();
}

//============================================================================
edm::RunNumber_t TkModuleGroupSelector::firstRunOfIOV(unsigned int iovNum) const
{
  return iovNum < this->numIovs() ? globalRunRange_.at(iovNum) : 0;
}

//======================================================================
int TkModuleGroupSelector::getParameterIndexFromDetId(unsigned int detId,
                                                      edm::RunNumber_t run) const
{
  // Return the index of the parameter that is used for this DetId.
  // If this DetId is not treated, return values < 0.

  const DetId temp_id(detId);

  int index = -1;

  bool sel = false;
  for(std::vector<int>::const_iterator itSubDets = subdetids_.begin();
      itSubDets != subdetids_.end();
      itSubDets++) {
    if (temp_id.det() == DetId::Tracker && temp_id.subdetId() == (*itSubDets)) {
      sel = true;
      break;
    }
  }

  if (temp_id.det() != DetId::Tracker || !sel) return -1;
  
  std::map<unsigned int, unsigned int>::const_iterator it = mapDetIdGroupId_.find(detId);
  if(it != mapDetIdGroupId_.end()) {
    const unsigned int iAlignableGroup = (*it).second;
    const std::vector<edm::RunNumber_t> &runs = runRange_.at(iAlignableGroup);
    const unsigned int id0 = firstId_.at(iAlignableGroup);
    const edm::RunNumber_t refrun = referenceRun_.at(iAlignableGroup);
    // assuming runs is never empty (checked in createModuleGroups(..))
    if (runs[0] > run) {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector::getParameterIndexFromDetId:\n"
        << "Run " << run << " not foreseen for detid ('"<< detId <<"').";
    }
    unsigned int iovNum = 0;
    for ( ; iovNum < runs.size(); ++iovNum) {
      if (run >= runs[iovNum]) break;
    }
    //test whether the iov contains the reference run
    if(refrun > 0) { //if > 0 a reference run number has been provided
      if(iovNum+1 == runs.size()) {
        if(refrun >= runs[iovNum])
          return -1;
      } else if( (iovNum+1) < runs.size()) {
        if(refrun >= runs[iovNum] && refrun < runs[iovNum+1]) {
          return -1;
        }
      }
    }
    index = id0 + iovNum;
  }
  return index;
}
