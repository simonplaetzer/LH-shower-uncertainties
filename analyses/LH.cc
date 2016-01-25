// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include <cmath>

namespace Rivet {

  class LH : public Analysis {
  public:

    typedef std::vector<Histo1DPtr> HistVec;
    typedef std::map<size_t, Histo1DPtr> HistMap;
    typedef std::map<std::pair<size_t, size_t>, Histo1DPtr> Hist2Map;

    LH(const string& name, PdgId particle, std::size_t jets = 1) : Analysis(name),
      _target_particle(particle), _mass(100), _njet(jets),
      _nbins(100),
      _h_jet_pT(_njet), _h_jet_mass(_njet), _h_jet_y(_njet), _h_jet_y_plus(_njet), _h_jet_y_minus(_njet),
      _h_X_jet_dR(_njet), _h_X_jet_dphi(_njet), _h_X_jet_dy(_njet)
    {
      
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Cuts
      Cut cut = Cuts::open();


      FinalState fs;
      addProjection(fs,"FS");

      // X Definition
      IdentifiedFinalState xfinder(cut);
      xfinder.acceptId(_target_particle);
      addProjection(xfinder, "Xfinder");

      //Jet Definition
      VetoedFinalState nonXFS;
      nonXFS.addVetoOnThisFinalState(xfinder);
      FastJets jetpro(nonXFS, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "jets");

      // X inclusive observables
      _h_X_pT = bookHisto1D("X_pT", logspace(_nbins, 1.0, 2.5*_mass/GeV));
      _h_X_pT_peak = bookHisto1D("X_pT_peak", _nbins, 0.0, 50.0);
      _h_X_y = bookHisto1D("X_y", _nbins, -6, 6);

      // Jet observables
      for(size_t jetA = 0; jetA < _njet; ++jetA) {
        const std::string jetA_str = "jet"+to_str(jetA+1);
        // Jet inclusive
        _h_jet_pT[jetA] = bookHisto1D(jetA_str+"_pT", logspace(_nbins, _jetptcut, 2.5*_mass/GeV));
        _h_jet_mass[jetA] = bookHisto1D(jetA_str+"_mass", logspace(_nbins, 1.0, 100));
        
        _h_jet_y[jetA] = bookHisto1D(jetA_str+"_y", _nbins, -6.0, 6.0);
        _h_jet_y_plus[jetA].reset(new Histo1D(_nbins, 0, 5)); // Unnamed Histrogram
        _h_jet_y_minus[jetA].reset(new Histo1D(_nbins, 0, 5)); // Unnamed Histrogram
        
        _h_X_jet_dR[jetA] = bookHisto1D("X_"+jetA_str+"_dR", _nbins, 0.0, 8.0);
        _h_X_jet_dy[jetA] = bookHisto1D("X_"+jetA_str+"_dy", _nbins, -6.0, 6.0);
        _h_X_jet_dphi[jetA] = bookHisto1D("X_"+jetA_str+"_dphi", _nbins, 0.0, PI);
        
        
        // Inter-jet observables
        for (size_t jetB = jetA+1; jetB < _njet; ++jetB) {
            const std::string jetB_str = "jet"+to_str(jetB+1);
            const std::string jetA_jetB_str = jetA_str+"_"+jetB_str;
            const std::pair<size_t, size_t> ij = std::make_pair(jetA, jetB);
            
            _h_jet_jet_dR.insert(make_pair(ij, bookHisto1D(jetA_jetB_str+"_dR", _nbins, 0.0, 8.0)));
            _h_jet_jet_dy.insert(make_pair(ij, bookHisto1D(jetA_jetB_str+"_dy", _nbins, -6.0, 6.0)));
            _h_jet_jet_dphi.insert(make_pair(ij, bookHisto1D(jetA_jetB_str+"_dphi", _nbins, 0.0, PI)));
          }
        
      }


      // Exclusive observables
      for (size_t exc_jets = 0; exc_jets <= _njet; ++exc_jets) {
        const std::string exc_str = to_str(exc_jets)+"jet_exclusive";
        // X properties
        _h_exc_X_pT.insert(make_pair(exc_jets,  bookHisto1D(exc_str+"_X_pT", logspace(_nbins, 1.0, 2.5*_mass/GeV))));
        _h_exc_X_pT_peak.insert(make_pair(exc_jets,  bookHisto1D(exc_str+"_X_pT_peak", _nbins, 0.0, 50.0)));
        _h_exc_X_y.insert(make_pair(exc_jets,  bookHisto1D(exc_str+"_X_y", _nbins, -6, 6)));
        
        // Jet properties
        for (size_t jetA = 0; jetA < exc_jets; ++jetA){
          const std::string jetA_str = "jet"+to_str(jetA+1);
          const std::pair<size_t, size_t> exc_j = std::make_pair(exc_jets, jetA);

          _h_exc_jet_pT.insert(make_pair(exc_j, bookHisto1D(exc_str+"_"+jetA_str+"_pT", logspace(_nbins, _jetptcut, 2.5*_mass/GeV))));
          _h_exc_jet_y.insert(make_pair(exc_j, bookHisto1D(exc_str+"_"+jetA_str+"_y", _nbins, -6.0, 6.0)));
          _h_exc_jet_mass.insert(make_pair(exc_j, bookHisto1D(exc_str+"_"+jetA_str+"_mass", logspace(_nbins, 1.0, 100))));
        }
      }
      
      _h_jet_multi_exclusive = bookHisto1D("jet_multi_exclusive", _njet+1, -0.5, _njet+0.5);
      _h_jet_multi_inclusive = bookHisto1D("jet_multi_inclusive", _njet+1, -0.5, _njet+0.5);
      _h_jet_multi_ratio = bookScatter2D("jet_multi_ratio");
      _h_jet_HT = bookHisto1D("jet_HT", logspace(_nbins, _jetptcut, 2.5*_mass/GeV));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      Particles xMatches = applyProjection<IdentifiedFinalState>(event, "Xfinder").particles();
      if (xMatches.size() != 1) {
        vetoEvent;
      }

      // Analysis Objects
      FourMomentum xMom(xMatches[0].momentum());
      
      const Jets& jets = applyProjection<FastJets>(event, "jets").jetsByPt(_jetptcut);

      // X inclusive observables
      _h_X_pT->fill(xMom.pT()/GeV, weight);
      _h_X_pT_peak->fill(xMom.pT()/GeV, weight);
      _h_X_y->fill(xMom.rapidity(), weight);

      // Inclusive
      for (size_t jetA = 0; jetA < _njet; ++jetA) {
        if (jets.size() < jetA + 1){ break; } // skip for where there are less jets than we're looking at
        // X-jetA distances
        _h_X_jet_dR[jetA] -> fill(deltaR(xMom, jets[jetA].momentum()), weight);
        _h_X_jet_dy[jetA] -> fill((xMom.rapidity()-jets[jetA].rapidity()), weight);
        _h_X_jet_dphi[jetA] -> fill(deltaPhi(xMom, jets[jetA].momentum()), weight);
        
        
        // JetA distributions
        _h_jet_pT[jetA] -> fill(jets[jetA].pT(), weight);
        double jetA_mass2 = jets[jetA].mass2();
        _h_jet_mass[jetA] -> fill((jetA_mass2 < 0.0 ? 0.0 : sqrt(jetA_mass2)/GeV), weight);
        
        
        const double jetA_rap = jets[jetA].rapidity();
        _h_jet_y[jetA] -> fill(jetA_rap, weight);
        (jetA_rap > 0.0 ? _h_jet_y_plus : _h_jet_y_minus)[jetA]->fill(fabs(jetA_rap), weight);

        
        // jetA-jetB distances
        for (size_t jetB = jetA+1; jetB < _njet; ++jetB) {
          if (jets.size() < jetB + 1){ break; } // cannot book on jetB
          const std::pair<size_t, size_t> ij = std::make_pair(jetA, jetB);
          _h_jet_jet_dR[ij] -> fill(deltaR(jets[jetA].momentum(), jets[jetB].momentum()), weight);
          _h_jet_jet_dy[ij] -> fill((jets[jetA].rapidity() - jets[jetB].rapidity()), weight);
          _h_jet_jet_dphi[ij] -> fill(deltaPhi(jets[jetA].momentum(), jets[jetB].momentum()), weight);
        }
      }


      // Exclusive
      if (jets.size() <= _njet) {
        // X distributions
        _h_exc_X_pT[jets.size()]->fill(xMom.pT()/GeV, weight);
        _h_exc_X_pT_peak[jets.size()]->fill(xMom.pT()/GeV, weight);
        _h_exc_X_y[jets.size()]-> fill(xMom.rapidity(), weight);
        
        // Jet distributions
        for (size_t jetA = 0; jetA < jets.size(); ++jetA) {
          const std::pair<size_t, size_t> exc_j = std::make_pair(jets.size(), jetA);
          _h_exc_jet_pT[exc_j]->fill(jets[jetA].pT()/GeV, weight);
          _h_exc_jet_y[exc_j]->fill(jets[jetA].rapidity(), weight);
          double jetA_mass2 = jets[jetA].mass2();
          _h_exc_jet_mass[exc_j]->fill((jetA_mass2 < 0.0 ? 0.0 : sqrt(jetA_mass2)/GeV), weight);
        }        
      }
      
      double jet_HT = 0.0;
      foreach (const Jet& jet, jets) {
        jet_HT += jet.pT();
      }
      _h_jet_HT->fill(jet_HT, weight);

      // Jet Multiplicities
      _h_jet_multi_exclusive -> fill (jets.size(), weight);
      for(size_t i = 0; i < jets.size(); ++i) {
       _h_jet_multi_inclusive -> fill(i, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double norm = 1.0; // crossSection()/picobarn;
      
      // X inclusive observables
      normalize(_h_X_pT, norm);
      normalize(_h_X_pT_peak, norm);
      normalize(_h_X_y, norm);
      
      // Jet observables
      foreach (HistVec::value_type& h, _h_jet_pT) normalize(h, norm);
      foreach (HistVec::value_type& h, _h_jet_mass) { normalize(h, norm); }
      foreach (HistVec::value_type& h, _h_jet_y) { normalize(h, norm); }
      for (size_t jetA = 0; jetA < _njet; ++jetA) {
          divide(*_h_jet_y_plus[jetA], *_h_jet_y_minus[jetA], bookScatter2D("jet"+to_str(jetA+1)+"_y_pmratio"));
      }
      
      normalize(_h_jet_HT, norm);
      normalize(_h_jet_multi_exclusive, norm);
      normalize(_h_jet_multi_inclusive, norm);
      
      // Fill inclusive jet multi ratio
      int Nbins = _h_jet_multi_inclusive->numBins();
      for (int i = 0; i < Nbins-1; ++i) {
        _h_jet_multi_ratio->addPoint(i+1, 0, 0.5, 0);
        if (_h_jet_multi_inclusive->bin(i).sumW() > 0.0) {
          const double ratio = _h_jet_multi_inclusive->bin(i+1).sumW()/_h_jet_multi_inclusive->bin(i).sumW();
          const double relerr_i = _h_jet_multi_inclusive->bin(i).relErr();
          const double relerr_j = _h_jet_multi_inclusive->bin(i+1).relErr();
          const double err = ratio * (relerr_i + relerr_j);
          _h_jet_multi_ratio->point(i).setY(ratio, err);
        }
      }
      
      foreach (HistVec::value_type& h, _h_X_jet_dR) { normalize(h, norm); }
      foreach (HistVec::value_type& h, _h_X_jet_dphi) { normalize(h, norm); }
      foreach (HistVec::value_type& h, _h_X_jet_dy) { normalize(h, norm); }
      
      //// Inter-jet distances
      foreach (Hist2Map::value_type& h, _h_jet_jet_dR) { normalize(h.second, norm); }
      foreach (Hist2Map::value_type& h, _h_jet_jet_dphi) { normalize(h.second, norm); }
      foreach (Hist2Map::value_type& h, _h_jet_jet_dy) { normalize(h.second, norm); }

      // X + m jet exclusive observables
      //// X distributions
      foreach (HistMap::value_type& h, _h_exc_X_pT) { normalize(h.second, norm); }
      foreach (HistMap::value_type& h, _h_exc_X_pT_peak) { normalize(h.second, norm); }
      foreach (HistMap::value_type& h, _h_exc_X_y) { normalize(h.second, norm); }
      
      //// Jet distributions for 1 < n <= m
      foreach (Hist2Map::value_type& h, _h_exc_jet_pT) { normalize(h.second, norm); }
      foreach (Hist2Map::value_type& h, _h_exc_jet_mass) { normalize(h.second, norm); }
      foreach (Hist2Map::value_type& h, _h_exc_jet_y) { normalize(h.second, norm); }
    }

    //@}
  protected:
    PdgId _target_particle;
    double _mass;
    std::size_t _njet;
    int _nbins;
    double _jetptcut;

  private:

    // Data members like post-cuts event weight counters go here


    /// @name Histograms
    //@{

    // X inclusive observables
    Histo1DPtr _h_X_pT;
    Histo1DPtr _h_X_pT_peak;
    Histo1DPtr _h_X_y;

    // Jet observables
    HistVec _h_jet_pT; 
    HistVec _h_jet_mass; 
    HistVec _h_jet_y;
    HistVec _h_jet_y_plus;
    HistVec _h_jet_y_minus;

    Histo1DPtr _h_jet_HT;
    Histo1DPtr _h_jet_multi_exclusive;
    Histo1DPtr _h_jet_multi_inclusive;
    Scatter2DPtr _h_jet_multi_ratio;

    // X + >0 jet inclusive observables
    //// (X,nth jet)-distances
    HistVec _h_X_jet_dR;
    HistVec _h_X_jet_dphi;
    HistVec _h_X_jet_dy;
    
    //// Inter-jet distances
    Hist2Map _h_jet_jet_dR;
    Hist2Map _h_jet_jet_dphi;
    Hist2Map _h_jet_jet_dy;

    // X + m jet exclusive observables
    //// X distributions
    HistMap _h_exc_X_pT;
    HistMap _h_exc_X_pT_peak;
    HistMap _h_exc_X_y;

    //// Jet distributions for 1 < n <= m
    Hist2Map _h_exc_jet_pT;
    Hist2Map _h_exc_jet_mass;
    Hist2Map _h_exc_jet_y;
    
    //@}

  };

struct LH_H : public LH {
  LH_H() : LH("LH_H", PID::HIGGS, 1) {
    _mass = 125*GeV;
    _jetptcut = 40*GeV;
  }
};

struct LH_Z : public LH {
  LH_Z() : LH("LH_Z", PID::Z0, 1) {
    _mass = 91*GeV;
    _jetptcut = 40*GeV;
  }
};

struct LH_H_Heavy : public LH {
  LH_H_Heavy() : LH("LH_H_Heavy", PID::HIGGS, 1) {
    _mass = 500*GeV;
    _jetptcut = 40*GeV;
  }
};

struct LH_Z_Heavy : public LH {
  LH_Z_Heavy() : LH("LH_Z_Heavy", PID::Z0, 1) {
    _mass = 500*GeV;
    _jetptcut = 40*GeV;
  }
};


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LH_H);
  DECLARE_RIVET_PLUGIN(LH_Z);
  DECLARE_RIVET_PLUGIN(LH_H_Heavy);
  DECLARE_RIVET_PLUGIN(LH_Z_Heavy);

}
