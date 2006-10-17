#ifndef Alignment_CommonAlignmentAlgorithm_AlignmentParameterBuilder_h
#define Alignment_CommonAlignmentAlgorithm_AlignmentParameterBuilder_h

#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"

#include <string>
#include <vector>

/// Build Alignment Parameter Structure 

class AlignmentParameterBuilder 
{
public:

  /// Constructor
  AlignmentParameterBuilder(AlignableTracker* alignableTracker);

  /// destructor 
  ~AlignmentParameterBuilder() {};

  /// Add predefined selection of alignables defined by a string 
  void addSelection(const std::string &name, const std::vector<bool> &sel);

  /// Add arbitrary selection of Alignables 
  void add (const std::vector<Alignable*>& alignables, const std::vector<bool> &sel);

  /// Add all level 1 objects (Dets) 
  void addAllDets(const std::vector<bool> &sel);

  /// Add all level 2 (Rod or equivalent) objects 
  void addAllRods(const std::vector<bool> &sel);

  /// Add all level 3 (Layer or equivalent) objects 
  void addAllLayers(const std::vector<bool> &sel);

  /// Add all level 4 (Halfbarrel etc) objects 
  void addAllComponents(const std::vector<bool> &sel);

  /// Add all alignables 
  void addAllAlignables(const std::vector<bool> &sel);

  /// Get list of alignables for which AlignmentParameters are built 
  std::vector<Alignable*> alignables() { return theAlignables; };

  /// Remove n Alignables from list 
  void fixAlignables( int n );

  
private:

  // data members

  /// Vector of alignables 
  std::vector<Alignable*> theAlignables;

  /// Alignable tracker   
  AlignableTracker* theAlignableTracker;

  /// Alignable id converter
  TrackerAlignableId* theTrackerAlignableId;

  bool theOnlyDS;
  bool theOnlySS;
  int  theMinLayer,theMaxLayer;
  bool theSelLayers;

};

#endif
