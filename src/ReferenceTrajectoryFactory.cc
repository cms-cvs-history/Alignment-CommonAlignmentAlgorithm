
#include "Alignment/CommonAlignmentAlgorithm/interface/ReferenceTrajectoryFactory.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include <iostream>


ReferenceTrajectoryFactory::ReferenceTrajectoryFactory( const edm::ParameterSet & config ) :
  TrajectoryFactoryBase( config )
{
  theMass = config.getParameter< double >( "ParticleMass" );
}

 
ReferenceTrajectoryFactory::~ReferenceTrajectoryFactory( void ) {}


const ReferenceTrajectoryFactory::ReferenceTrajectoryCollection
ReferenceTrajectoryFactory::trajectories( const edm::EventSetup & setup,
					  const TrajTrackPairCollection & tracks ) const
{
  ReferenceTrajectoryCollection trajectories;

  edm::ESHandle< MagneticField > magneticField;
  setup.get< IdealMagneticFieldRecord >().get( magneticField );

  TrajTrackPairCollection::const_iterator itTracks = tracks.begin();

  while ( itTracks != tracks.end() )
  { 
    TrajectoryInput input = innermostStateAndRecHits( *itTracks );
    // set the flag for reversing the RecHits to false, since they are already in the correct order.
    trajectories.push_back( ReferenceTrajectoryPtr( new ReferenceTrajectory( input.first, input.second, 
									     false, magneticField.product(),
									     theMaterialEffects, theMass ) ) );
    ++itTracks;
  }

  return trajectories;
}
