/**
 * \file TkModuleGroupSelector.cc
 *
 *  \author Joerg Behr
 *  \date May 2013
 *  $Revision: 1.1.2.5 $
 *  $Date: 2013/05/17 13:20:19 $
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
                                                                              nparameters_(0)
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
void TkModuleGroupSelector::setReferenceRunRange(const edm::ParameterSet &cfg)
{
  //extract the reference run range if defined
  if(cfg.exists("ReferenceRunRange")) {
    referenceRunRange_ = cfg.getParameter<std::vector<edm::RunNumber_t> >("ReferenceRunRange");
  }
  if(referenceRunRange_.size() > 0) {
    if(referenceRunRange_.size() != 2 || referenceRunRange_.at(0) >= referenceRunRange_.at(1)) {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector::setReferenceRunRange:\n"
        << " Excactly two ordered run numbers have to be provided for the reference run range.";
    }
  }
}

//============================================================================
const std::vector<edm::RunNumber_t>& TkModuleGroupSelector::getReferenceRunRange() const
{
  return referenceRunRange_;
}

//============================================================================
const bool TkModuleGroupSelector::testSplitOption(const edm::ParameterSet &pset) const
{
  bool split = false;
  if(pset.exists("split")) {
    split = pset.getParameter<bool>("split");
  }
  const size_t npar = pset.getParameterNames().size();
  
  if (npar >= 3 && !pset.exists("split")) {
    throw cms::Exception("BadConfig")
      << "@SUB=TkModuleGroupSelector:createModuleGroups:"
      << " >= 3 parameters specified in PSet BUT split parameter was not found! Maybe a typo?";
  }
  return split;
}

//============================================================================
void TkModuleGroupSelector::testGlobalRunRangeOrder() const
{
  edm::RunNumber_t firstRun = 0; 
  for(std::vector<edm::RunNumber_t>::const_iterator iRun = globalRunRange_.begin(); 
      iRun != globalRunRange_.end(); ++iRun)  {
    if((*iRun) > firstRun) {
      firstRun = (*iRun);
    } else {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector:createModuleGroups:"
        << " Global run range vector not sorted.";
    }
  }
}

//============================================================================
bool TkModuleGroupSelector::createGroup(
                                        const bool split, 
                                        unsigned int &Id, 
                                        const std::vector<edm::RunNumber_t> &range, 
                                        Alignable* iD, 
                                        const std::list<Alignable*> &selected_alis
                                        )
{
  bool modules_selected = false;

  if(iD != NULL && selected_alis.size() == 0) {
    if(split) {
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
void TkModuleGroupSelector::createModuleGroups(AlignableTracker *aliTracker,
                                               AlignableMuon *aliMuon,
                                               AlignableExtras *aliExtras)
{
  std::set<edm::RunNumber_t> localRunRange;
  nparameters_ = 0;
  unsigned int Id = 0;
  //loop over all LA groups
  for(edm::VParameterSet::const_iterator pset = myGranularityConfig_.begin();
      pset != myGranularityConfig_.end();
      ++pset) {
    bool modules_selected = false; //track whether at all a module has been selected in this group
    const std::vector<edm::RunNumber_t> range = pset->getParameter<std::vector<edm::RunNumber_t> >("RunRange");
    if(range.size() == 0) {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector::createModuleGroups:\n"
        << "Run range array empty!";
    }
    const bool split = this->testSplitOption((*pset));
    
    AlignmentParameterSelector selector(aliTracker,aliMuon, aliExtras);
    selector.clear();
    selector.addSelections((*pset).getParameter<edm::ParameterSet> ("levels"));

    const std::vector<Alignable*> &alis = selector.selectedAlignables();
    std::list<Alignable*> selected_alis;
    for(std::vector<Alignable*>::const_iterator it = alis.begin(); it != alis.end(); it++) {
      const std::vector<Alignable*> &aliDaughts = (*it)->deepComponents();
      if(aliDaughts.size() > 0) {
        for (std::vector<Alignable*>::const_iterator iD = aliDaughts.begin();
             iD != aliDaughts.end(); ++iD) {
          if((*iD)->alignableObjectId() == align::AlignableDetUnit || (*iD)->alignableObjectId() == align::AlignableDet) {
            if(split) {
              modules_selected = this->createGroup(split, Id, range, (*iD), std::list<Alignable*>());//last parameter is a empty dummy list
            } else {
              selected_alis.push_back((*iD));
            }
          }
        }
      }
    }
    
    if(!split) {
      modules_selected = this->createGroup(split, Id, range, NULL, selected_alis);
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

  // test the order of the runs
  this->testGlobalRunRangeOrder();

 
  
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
  
  // Check whether run lies inside the reference run range. Size of vector has been checked in setReferenceRunRange(...)
  if(referenceRunRange_.size() == 2 && run >= referenceRunRange_.at(0) && run <= referenceRunRange_.at(1)) {
    return -1;
  }

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
    index = id0 + iovNum;
  }
  return index;
}
