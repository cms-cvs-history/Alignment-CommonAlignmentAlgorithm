#ifndef Alignment_CommonAlignmentAlgorithm_TkModuleGroupSelector_h
#define Alignment_CommonAlignmentAlgorithm_TkModuleGroupSelector_h

/**
 * \file TkModuleGroupSelector.cc
 *
 * Class provides an algorithm which assigns
 * runrange-dependent (IOV-dependent)
 * indices to groups of tracker modules.
 *
 *  \author Joerg Behr
 *  \date May 2013
 *  $Revision: 1.1.2.1 $
 *  $Date: 2013/05/10 12:54:22 $
 *  (last update by $Author: jbehr $)
 *
 */


#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h"

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

#include <vector>
#include <utility>
#include <string>
#include <map>
#include <sstream>
#include <set>
#include <functional>

class AlignableTracker;
class AlignableMuon;
class AlignableExtras;

namespace edm { class EventSetup; class ParameterSet; } 

class TkModuleGroupSelector
{
public:
  /// Constructor
  //FIXME: take list of subdetids from which modules are taken
  explicit TkModuleGroupSelector(const edm::VParameterSet &cfg);
  
  /// Destructor
  ~TkModuleGroupSelector() {};

  // Reads and parses the configurations and
  // constructs the run-dependent module groups.
  void CreateModuleGroups(AlignableTracker *aliTracker,
                          AlignableMuon *aliMuon,
                          AlignableExtras *aliExtras);
  
  // Specify the sub-detectors for which modules are grouped together.
  // Modules belonging to other sub-detectors are ignored.
  void SetSubDets(std::vector<int> sdets); //FIXME: move somehow to constructor?

  // Returns the number of parameters.
  unsigned int GetNumberOfParameters() const;

  /// Total number of IOVs.
  unsigned int numIovs() const;

  /// First run of iov (0 if iovNum not treated).
  edm::RunNumber_t firstRunOfIOV(unsigned int iovNum) const;

  /// Index of parameter for given detId (detId not treated => < 0)
  /// and the given run.
  int getParameterIndexFromDetId(unsigned int detId, edm::RunNumber_t run) const;
  
 private:
  const edm::VParameterSet myGranularityConfig_;
  std::vector<edm::RunNumber_t> globalRunRange_;
  std::vector<unsigned int> firstId_;//1:1 mapping with assignment_
  std::vector<std::pair<std::list<Alignable*>, std::vector<edm::RunNumber_t> > > assignment_;
  unsigned int nparameters_;
  std::vector<int> subdetids_;
};

#endif
