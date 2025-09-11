// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Barak Schmookler

#pragma once

#include <Math/Vector4D.h>
#include <Math/LorentzRotation.h>
#include <Math/LorentzVector.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/Boost.h>

using ROOT::Math::PxPyPzEVector;

namespace eicrecon {

using ROOT::Math::LorentzRotation;

inline LorentzRotation determine_boost(PxPyPzEVector ei, PxPyPzEVector pi) {

  using ROOT::Math::Boost;
  using ROOT::Math::RotationX;
  using ROOT::Math::RotationY;

  // Step 1: Find the needed boosts and rotations from the incoming lepton and hadron beams
  // (note, this will give you a perfect boost, in principle you will not know the beam momenta exactly and should use an average)

  PxPyPzEVector eo = ei;
  PxPyPzEVector po = pi;

  // Define the Boost to make beams back-to-back
  const auto cmBoost = (ei + pi).BoostToCM();
  const Boost boost_to_cm(cmBoost);

  // Boost to COM frame
  pi = boost_to_cm(pi);
  ei = boost_to_cm(ei);

  // This will boost beams from a center of momentum frame back to (nearly) their original energies
  PxPyPzEVector eh(0, 0, -1*sqrt(pow(eo.E(),2)-pow(eo.M(),2)), eo.E());
  PxPyPzEVector ph(0, 0,    sqrt(pow(po.E(),2)-pow(po.M(),2)), po.E());

  const auto hoBoost = (eh + ph).BoostToCM();
  const Boost boost_to_headon(-hoBoost);

  // Boost and rotate the incoming beams to find the proper rotations TLorentzVector

  // Rotate to head-on
  RotationY rotAboutY(-1.0 * atan2(pi.Px(), pi.Pz())); // Rotate to remove x component of beams
  RotationX rotAboutX(+1.0 * atan2(pi.Py(), pi.Pz())); // Rotate to remove y component of beams

  // final matrix: P' = [BtoH][RX][RY][BtoCM]P <-- Matrix multi. goes from R to L
  LorentzRotation tf;
  tf *= boost_to_headon;
  tf *= rotAboutX;
  tf *= rotAboutY;
  tf *= boost_to_cm;

  return tf;
}

inline PxPyPzEVector apply_boost(const LorentzRotation& tf, PxPyPzEVector part) {

  // Step 2: Apply boosts and rotations to any particle 4-vector
  // (here too, choices will have to be made as to what the 4-vector is for reconstructed particles)

  // Boost and rotate particle 4-momenta into the headon frame
  tf(part);
  return part;
}

} // namespace eicrecon
