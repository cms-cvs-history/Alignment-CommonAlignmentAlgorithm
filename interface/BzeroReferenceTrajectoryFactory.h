#ifndef Alignment_CommonAlignmentAlgorithm_BzeroReferenceTrajectoryFactory_h
#define Alignment_CommonAlignmentAlgorithm_BzeroReferenceTrajectoryFactory_h

#include "Alignment/CommonAlignmentAlgorithm/interface/TrajectoryFactoryBase.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/BzeroReferenceTrajectory.h"

/// A factory that produces instances of class BzeroReferenceTrajectory from a given TrajTrackPairCollection.


class BzeroReferenceTrajectoryFactory : public TrajectoryFactoryBase
{

public:

  BzeroReferenceTrajectoryFactory( const edm::ParameterSet & config );
  virtual ~BzeroReferenceTrajectoryFactory( void );

  /// Produce the reference trajectories.
  virtual const ReferenceTrajectoryCollection trajectories( const edm::EventSetup & setup,
							    const ConstTrajTrackPairCollection & tracks ) const;

  virtual BzeroReferenceTrajectoryFactory* clone( void ) const { return new BzeroReferenceTrajectoryFactory( *this ); }

protected:

  double theMass;
  double theMomentumEstimate;
};


#endif
