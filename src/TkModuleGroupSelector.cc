/**
 * \file TkModuleGroupSelector.cc
 *
 *  \author Joerg Behr
 *  \date May 2013
 *  $Revision: 1.1.2.1 $
 *  $Date: 2013/05/10 12:54:27 $
 *  (last update by $Author: jbehr $)
 */

#include "Alignment/CommonAlignmentAlgorithm/interface/TkModuleGroupSelector.h"

//============================================================================
TkModuleGroupSelector::TkModuleGroupSelector(const edm::VParameterSet &cfg) : myGranularityConfig_(cfg),
                                                                              nparameters_(0)
{
  //FIXME: take list of subdetids from which modules are taken
}

//============================================================================
void TkModuleGroupSelector::SetSubDets(std::vector<int> sdets)
{
  subdetids_ = sdets;
}

//============================================================================
void TkModuleGroupSelector::CreateModuleGroups(AlignableTracker *aliTracker,
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
    const std::vector<edm::RunNumber_t> range = pset->getParameter<std::vector<edm::RunNumber_t> >("RunRange");

    bool split = false;
    if((*pset).exists("split")) {
      split = pset->getParameter<bool>("split");
    }
    const size_t npar = (*pset).getParameterNames().size();
    
    if (npar >= 3 && !(*pset).exists("split")) {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector:CreateModuleGroups:"
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
          assignment_.push_back(std::make_pair(std::list<Alignable*>(1,(*it)), range));
          firstId_.push_back(Id);
          Id += range.size();
          nparameters_ += range.size();
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
              const std::list<Alignable*> &a = assignment_.at(iAlignableGroup).first;
              
              for(std::list<Alignable*>::const_iterator iAli = a.begin();
                  iAli != a.end(); iAli++) {
                if( (*iAli)->alignableObjectId() == (*iD)->alignableObjectId()
                    && (*iAli)->id() == (*iD)->id()) {
                  found = true;
                  if(!split) {
                    throw cms::Exception("BadConfig")
                      << "@SUB=TkModuleGroupSelector:CreateModuleGroups:"
                      << " Module already selected.";
                  }
                  break;
                }
              }
            }
            if(split) {
              if(!found) {
                assignment_.push_back(std::make_pair(std::list<Alignable*>(1,(*iD)), range));
                firstId_.push_back(Id);
                Id += range.size();
                nparameters_ += range.size();
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
      assignment_.push_back(std::make_pair(selected_alis, range));
      firstId_.push_back(Id);
      Id += range.size();
      nparameters_ += range.size();
    }
    
   
    edm::RunNumber_t firstRun = 0; 
    for(std::vector<edm::RunNumber_t>::const_iterator iRun = range.begin(); 
        iRun != range.end(); ++iRun)  {
      localRunRange.insert((*iRun));
      if((*iRun) > firstRun) {
        firstRun = (*iRun);
      } else {
        throw cms::Exception("BadConfig")
          << "@SUB=TkModuleGroupSelector::CreateModuleGroups:"
          << " Run range not sorted.";
      }
    }
  }
 
  for(unsigned int iAlignableGroup = 0; iAlignableGroup < firstId_.size(); iAlignableGroup++) {
    const std::list<Alignable*> &alis = assignment_.at(iAlignableGroup).first;
    int nselectedmodules = 0;
    for(std::list<Alignable*>::const_iterator iAli = alis.begin();
        iAli != alis.end(); iAli++) {
      nselectedmodules++;
      bool selectedmodule = true;
      if((*iAli)->alignableObjectId() == align::AlignableDetUnit || (*iAli)->alignableObjectId() == align::AlignableDet) {
        const DetId id((*iAli)->id());
        if (id.det() != DetId::Tracker) selectedmodule = false;
        bool sel = false;
        for(std::vector<int>::const_iterator itSubDets = subdetids_.begin();
            itSubDets != subdetids_.end();
            itSubDets++) {
          if (id.det() == DetId::Tracker && id.subdetId() == (*itSubDets)) {
            sel = true;
            break;
          }
        }
        if(!sel) selectedmodule = false;
      } else {
        selectedmodule = false;
      }

      if(!selectedmodule) {
        throw cms::Exception("BadConfig") 
          << "@SUB=TkModuleGroupSelector:CreateModuleGroups:"
          << " Badly selected module in PSet number " << iAlignableGroup;//FIXME: improve text of exception
      }
    }

    
    //test whether at all a pixel module was selected in this group
    if(nselectedmodules == 0) {
      throw cms::Exception("BadConfig") 
        << "@SUB=TkModuleGroupSelector:CreateModuleGroups:"
        << " No module was selected in PSet number " <<iAlignableGroup << " of the module group selector.";//FIXME: improve text of exception
    }
  }
  
 
  //test whether at all a module was selected
  if(nparameters_ == 0) {
    throw cms::Exception("BadConfig") 
      << "@SUB=TkModuleGroupSelector:CreateModuleGroups:"
      << " No module was selected in the module group selector.";
  }

  //copy local set into the global vector of run boundaries
  for(std::set<edm::RunNumber_t>::const_iterator itRun = localRunRange.begin();
      itRun != localRunRange.end(); itRun++) {
    globalRunRange_.push_back((*itRun));
  }
  
  edm::RunNumber_t firstRun = 0; 
  for(std::vector<edm::RunNumber_t>::const_iterator iRun = globalRunRange_.begin(); 
      iRun != globalRunRange_.end(); ++iRun)  {
    if((*iRun) > firstRun) {
      firstRun = (*iRun);
    } else {
      throw cms::Exception("BadConfig")
        << "@SUB=TkModuleGroupSelector:CreateModuleGroups:"
        << " Global run range vector not sorted.";
    }
  }
  
}

//============================================================================
unsigned int TkModuleGroupSelector::GetNumberOfParameters() const
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
  
  int nmatched = 0;
  bool matched = false;
  for(unsigned int iAlignableGroup = 0; iAlignableGroup < firstId_.size(); iAlignableGroup++) {
    const unsigned int id0 = firstId_.at(iAlignableGroup);
    
    const std::list<Alignable*> &alis = assignment_.at(iAlignableGroup).first;
    const std::vector<edm::RunNumber_t> &runs = assignment_.at(iAlignableGroup).second;
    
  
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
	  << "@SUB=TkModuleGroupSelector::getParameterIndexFromDetId:\n"
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

      // //const int index_old = getParameterIndexFromDetId_old(detId,run);
      // std::cout << "debug indices " << index << " " << index_old << std::endl;


    }
  }

  if(nmatched >= 2) {
    throw cms::Exception("BadConfig")
      << "@SUB=TkModuleGroupSelector::getParameterIndexFromDetId:\n"
      << "Detid ('"<< detId <<"') is assigned to >= 2 groups of modules (in fact it is matched to '"<< nmatched <<"' groups).";
  }
  

  return index;
}
