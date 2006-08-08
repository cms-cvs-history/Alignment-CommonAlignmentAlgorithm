#ifndef Alignment_CommonAlignmentAlgorithm_TrackerAlignmentProducer_h
#define Alignment_CommonAlignmentAlgorithm_TrackerAlignmentProducer_h

//
// Package:    Alignment/CommonAlignmentAlgorithm
// Class:      AlignmentProducer
// 
//
// Description: calls alignment algorithms
//
//
// Original Author:  Frederic Ronga


#include <vector>

// Framework
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

// Alignment
#include "RecoTracker/TrackProducer/interface/TrackProducerBase.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h"
#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"

#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterBuilder.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterStore.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"

class AlignmentProducer :  public TrackProducerBase, public edm::ESProducerLooper
{

 public:

  /// Constructor
  AlignmentProducer( const edm::ParameterSet& iConfig );
  
  /// Destructor
  ~AlignmentProducer();

  // Define return type
  typedef boost::shared_ptr<TrackerGeometry> ReturnType;

  /// Produce the geometry
  virtual ReturnType produce( const TrackerDigiGeometryRecord& iRecord );

  /// Dummy implementation (job done in duringLoop)
  virtual void produce(edm::Event&, const edm::EventSetup&) {};

  /// Called at beginning of job
  virtual void beginOfJob(const edm::EventSetup&);

  /// Called at end of job
  virtual void endOfJob();

  /// Called at beginning of loop
  virtual void startingNewLoop( unsigned int iLoop );

  /// Called at end of loop
  virtual Status endOfLoop( const edm::EventSetup&, unsigned int iLoop );

  /// Called at each event 
  virtual Status duringLoop( const edm::Event&, const edm::EventSetup& );

 private:

  typedef std::vector<Alignable*> Alignables;

  AlignmentAlgorithmBase* theAlignmentAlgo;
  AlignmentParameterBuilder* theAlignmentParameterBuilder;
  AlignmentParameterStore* theAlignmentParameterStore;


  AlignableTracker* theAlignableTracker;
  TrackProducerAlgorithm theRefitterAlgo;
  ReturnType theTracker;

  unsigned int theMaxLoops;     // Number of loops to loop

  std::string theSrc; // For local use (bypass TrackerProducerBase)

  std::string stParameterSelector;
  std::string stAlignableSelector;

};

#endif
