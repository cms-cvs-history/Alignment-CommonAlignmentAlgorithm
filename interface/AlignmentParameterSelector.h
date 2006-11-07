#ifndef ALIGNMENTPARAMETERSELECTOR_H
#define ALIGNMENTPARAMETERSELECTOR_H

/** \class AlignmentParameterSelector
 *  \author Gero Flucke (selection by strings taken from AlignableParameterBuilder)
 *
 *  Selecting Alignable's of the tracker by predefined strings,
 *  additional constraints on eta, phi, r or z are possible.
 *  Furthermore stores the 'selection' of selected AlignmentParameters.
 *
 *  $Date: 2006/11/03 16:46:57 $
 *  $Revision: 1.1 $
 *  (last update by $Author: flucke $)
 */

#include <vector>
#include <string>

class Alignable;
class AlignableTracker;
namespace edm {
  class ParameterSet;
}

class AlignmentParameterSelector {
 public:
  /// Constructor
  explicit AlignmentParameterSelector(AlignableTracker *aliTracker);

  /// Destructor
  virtual ~AlignmentParameterSelector();

  /// vector of alignables selected so far
  const std::vector<Alignable*>& selectedAlignables() const;
  /// vector of selected parameters for alignables, parallel to selectedAlignables()
  const std::vector<std::vector<bool> >& selectedParameters() const;
  /// remove all selected Alignables and geometrical restrictions
  void clear();
  /// remove all geometrical restrictions
  void clearGeometryCuts();

  /// Add several selections defined by the PSet which must contain a vstring like e.g.
  /// vstring alignableParamSelector = { "PixelHalfBarrelLadders,111000,pixelSelection",
  ///                                    "BarrelDSRods,111000",
  ///                                    "BarrelSSRods,101000"}
  /// The '111000' is decoded into vector<bool> and stored.
  /// Returns number of added selections or -1 if problems (then also an error is logged)
  /// If a string contains a third, comma separated part (e.g. ',pixelSelection'),
  /// a further PSet of that name is expected to select eta/z/phi/r-ranges
  unsigned int addSelections(const edm::ParameterSet &pSet);
  /// set geometrical restrictions to be applied on all following selections
  /// (slices defined by vdouble 'etaRanges', 'phiRanges', 'zRanges' and 'rRanges',
  /// empty array means no restriction)
  void setGeometryCuts(const edm::ParameterSet &pSet);
  /// add Alignables corresponding to predefined name, taking into account geometrical restrictions,
  /// returns number of added alignables
    unsigned int addSelection(const std::string &name, const std::vector<bool> &paramSel);
  /// as addSelection with one argument, but overwriting geometrical restrictions
      unsigned int addSelection(const std::string &name, const std::vector<bool> &paramSel, 
				const edm::ParameterSet &pSet);

  /// true if geometrical restrictions in eta, phi, r, z not satisfied
  bool outsideRanges(const Alignable *alignable) const;
  /// true if ranges.size() is even and ranges[i] <= value < ranges[i+1] for any even i
  /// ( => false if ranges.empty() == true), if(isPhi==true) takes into account phi periodicity
  bool insideRanges(double value, const std::vector<double> &ranges, bool isPhi = false) const;
  /// Decomposing input string 's' into parts separated by 'delimiter'
  std::vector<std::string> decompose(const std::string &s, std::string::value_type delimiter) const;
  /// Decoding string to select local rigid body parameters into vector<bool>,
  /// e.g. "101001" will mean x, z and gamma, but not y, alpha and beta 
  /// for (Composite)RigidBodyAlignmentParameters
  std::vector<bool> decodeParamSel(const std::string &selString) const;

 protected:
  /// adding alignables which fulfil geometrical restrictions and special switches 
   unsigned int add(const std::vector<Alignable*> &alignables, const std::vector<bool> &paramSel);
  /// some helper methods
  unsigned int addAllDets(const std::vector<bool> &paramSel);
  unsigned int addAllRods(const std::vector<bool> &paramSel);
  unsigned int addAllLayers(const std::vector<bool> &paramSel);
  unsigned int addAllAlignables(const std::vector<bool> &paramSel);

 private:
  AlignableTracker        *theTracker;
  std::vector<Alignable*>  theSelectedAlignables;
  std::vector<std::vector<bool> > theSelectedParameters;

  /// geometrical restrictions in eta, phi, r, z to be applied for next addSelection
  std::vector<double> theRangesEta;
  std::vector<double> theRangesPhi;
  std::vector<double> theRangesR;
  std::vector<double> theRangesZ;

  // further switches used in add(...)
  bool theOnlyDS;
  bool theOnlySS;
  bool theSelLayers;
  int  theMinLayer;
  int  theMaxLayer;

};
#endif