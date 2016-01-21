// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief ALEPH jet rates and event shapes at LEP 1 and 2
  class MC_EETOJETS : public Analysis {
  public:

    MC_EETOJETS()
      : Analysis("MC_EETOJETS")
    {
    }


  public:

    void init() {
      const FinalState fs;
      addProjection(fs, "FS");
      FastJets durhamjets(fs, FastJets::DURHAM, 0.7);
      durhamjets.useInvisibles(true);
      addProjection(durhamjets, "DurhamJets");

      const Thrust thrust(fs);
      addProjection(thrust, "Thrust");
      addProjection(Sphericity(fs), "Sphericity");
      addProjection(ParisiTensor(fs), "Parisi");
      addProjection(Hemispheres(thrust), "Hemispheres");

      _h_thrust = bookHisto1D("Thrust", 50, 0.0, 1.0);
      _h_heavyjetmass = bookHisto1D("HeavyJetMass", 50, 0.0, 1.0);
      _h_totaljetbroadening = bookHisto1D("TotalJetBroadening", 50, 0.0, 1.0);
      _h_widejetbroadening = bookHisto1D("WideJetBroadening", 50, 0.0, 1.0);
      _h_cparameter = bookHisto1D("CParameter", 50, 0.0, 1.0);
      _h_thrustmajor = bookHisto1D("ThrustMajor", 50, 0.0, 1.0);
      _h_thrustminor = bookHisto1D("ThrustMinor", 50, 0.0, 1.0);
      _h_jetmassdifference = bookHisto1D("JetMassDifference", 50, 0.0, 1.0);
      _h_aplanarity = bookHisto1D("Aplanarity", 50, 0.0, 1.0);
      _h_planarity  = bookHisto1D("Planarity", 50, 0.0, 1.0);
      _h_oblateness = bookHisto1D("Oblateness", 50, 0.0, 1.0);
      _h_sphericity = bookHisto1D("Sphericity", 50, 0.0, 1.0);

      _h_y_Durham[0] = bookHisto1D("y12", 100, 0.0, 3.0);
      _h_y_Durham[1] = bookHisto1D("y23", 100, 1.0, 11.0);
      _h_y_Durham[2] = bookHisto1D("y34", 100, 2.0, 15.0);
      _h_y_Durham[3] = bookHisto1D("y45", 100, 3.0, 20.0);
      _h_y_Durham[4] = bookHisto1D("y56", 100, 3.5, 15.0);

      _h_R_Durham[0] = bookHisto1D("R1", 50, -1.0, 0);
      for (size_t j = 1; j < 5; ++j) { 
        _h_R_Durham[j] = bookHisto1D("R"+to_str(j+1), 50, -9.0-double(j), 0.0-double(j));
      }
      _h_R_Durham[5] = bookHisto1D("R6", 50, -13.0, -5.0);

    }

    void analyze(const Event& e) {

      const double weight = e.weight();

      const Thrust& thrust = applyProjection<Thrust>(e, "Thrust");
      const Sphericity& sphericity = applyProjection<Sphericity>(e, "Sphericity");

      _h_thrust->fill(1.0 - thrust.thrust(),weight);
      _h_thrustmajor->fill(thrust.thrustMajor(),weight);
      _h_thrustminor->fill(thrust.thrustMinor(),weight);
      _h_oblateness->fill(thrust.oblateness(),weight);

      const Hemispheres& hemi = applyProjection<Hemispheres>(e, "Hemispheres");
      _h_heavyjetmass->fill(hemi.scaledM2high(),weight);
      _h_jetmassdifference->fill(hemi.scaledM2diff(),weight);
      _h_totaljetbroadening->fill(hemi.Bsum(),weight);
      _h_widejetbroadening->fill(hemi.Bmax(),weight);

      const ParisiTensor& parisi = applyProjection<ParisiTensor>(e, "Parisi");
      _h_cparameter->fill(parisi.C(),weight);

      _h_aplanarity->fill(sphericity.aplanarity(),weight);
      _h_planarity->fill(sphericity.planarity(),weight);
      _h_sphericity->fill(sphericity.sphericity(),weight);

      // Jet rates
      const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
      if (durjet.clusterSeq()) {
	double logynm1=0.;
	double logyn;
	for (size_t i=0; i<5; ++i) {
	  double yn = durjet.clusterSeq()->exclusive_ymerge_max(i+1);
	  if (yn<=0.0) continue;
	  logyn = -log(yn);
	  if (_h_y_Durham[i]) {
	    _h_y_Durham[i]->fill(logyn, weight);
	  }
	  for (size_t j = 0; j < _h_R_Durham[i]->numBins(); ++j) {
	    double val   = _h_R_Durham[i]->bin(j).xMin();
	    double width = _h_R_Durham[i]->bin(j).xWidth();
	    if(-val<=logynm1) break;
	    if(-val<logyn) {
	      _h_R_Durham[i]->fill(val+0.5*width, weight*width);
	    }
	  }
	  logynm1 = logyn;
	}
	for (size_t j = 0; j < _h_R_Durham[5]->numBins(); ++j) {
	  double val   = _h_R_Durham[5]->bin(j).xMin();
	  double width = _h_R_Durham[5]->bin(j).xWidth();
	  if(-val<=logynm1) break;
	  _h_R_Durham[5]->fill(val+0.5*width, weight*width);
	}
      }

    }

    void finalize() {

      normalize(_h_thrust);
      normalize(_h_heavyjetmass);
      normalize(_h_totaljetbroadening);
      normalize(_h_widejetbroadening);
      normalize(_h_cparameter);
      normalize(_h_thrustmajor);
      normalize(_h_thrustminor);
      normalize(_h_jetmassdifference);
      normalize(_h_aplanarity);
      normalize(_h_planarity);
      normalize(_h_oblateness);
      normalize(_h_sphericity);

      for (size_t n=0; n<6; ++n) {
	scale(_h_R_Durham[n], 1./sumOfWeights());
      }

      for (size_t n = 0; n < 5; ++n) {
	if (_h_y_Durham[n]) {
	  scale(_h_y_Durham[n], 1.0/sumOfWeights());
	}
      }

    }

  private:

    Histo1DPtr _h_thrust;
    Histo1DPtr _h_heavyjetmass;
    Histo1DPtr _h_totaljetbroadening;
    Histo1DPtr _h_widejetbroadening;
    Histo1DPtr _h_cparameter;
    Histo1DPtr _h_thrustmajor;
    Histo1DPtr _h_thrustminor;
    Histo1DPtr _h_jetmassdifference;
    Histo1DPtr _h_aplanarity;
    Histo1DPtr _h_planarity;
    Histo1DPtr _h_oblateness;
    Histo1DPtr _h_sphericity;

    Histo1DPtr _h_R_Durham[6];
    Histo1DPtr _h_y_Durham[5];

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_EETOJETS);

}
