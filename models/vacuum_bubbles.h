#ifndef VACUUM_BUBBLES_H
#define VACUUM_BUBBLES_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): CosmoLattice contributors, Year: 2026

#include "CosmoInterface/cosmointerface.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace TempLat
{
  /**
   * A scalar-field benchmark for simultaneous vacuum-bubble nucleation.
   *
   * All quantities in this model are expressed in units of M, with M=1.  The
   * constructor fills the scalar field with the t=0 bubble profile.  The
   * normal CosmoLattice initializer is then asked for homogeneous initial
   * conditions, so it leaves this profile untouched and adds no Gaussian
   * fluctuations.
   */
  struct ModelPars : public TempLat::DefaultModelPars
  {
    static constexpr size_t NScalars = 1;
    static constexpr size_t NPotTerms = 1;

    using FloatType = double;
  };

#define MODELNAME vacuum_bubbles

  template <class R> using Model = MakeModel(R, ModelPars);

  class MODELNAME : public Model<MODELNAME>
  {
  private:
    FloatType delta, lambda, phiB, Rc, l0;
    int Nb;
    std::string bubbleSeed;

    template <typename S> static S getModelParameter(ParameterParser &parser, const std::string &name, S def)
    {
      // Accept the unprefixed names used by the production input and the
      // bubble_* aliases convenient for command-line smoke runs.
      const auto &params = parser.getParams();
      const auto alias = "bubble_" + name;
      if (params.count(alias) != 0) return parser.get<S>(alias);
      if (params.count(name) != 0) return parser.get<S>(name);
      return parser.get<S>(name, def);
    }

    struct BubbleCenter
    {
      FloatType x[3];
    };

    static FloatType minimumImageDistance(FloatType x, FloatType y, FloatType side)
    {
      const FloatType separation = std::abs(x - y);
      return std::min(separation, side - separation);
    }

    static FloatType centerDistance(const BubbleCenter &a, const BubbleCenter &b, FloatType side)
    {
      FloatType distanceSquared = 0;
      for (int i = 0; i < 3; ++i) {
        const FloatType separation = minimumImageDistance(a.x[i], b.x[i], side);
        distanceSquared += separation * separation;
      }
      return std::sqrt(distanceSquared);
    }

    std::vector<BubbleCenter> makeBubbleCenters(FloatType side) const
    {
      if (Nb <= 0) throw(RunParametersInconsistent("Nb must be positive for vacuum_bubbles."));
      if (side <= 0) throw(RunParametersInconsistent("The lattice side length must be positive for vacuum_bubbles."));

      // A seed sequence built from the string gives deterministic placement
      // without depending on std::hash's implementation-defined result.
      std::vector<std::uint32_t> seedWords;
      seedWords.reserve((bubbleSeed.size() + 3) / 4);
      for (std::size_t i = 0; i < bubbleSeed.size(); i += 4) {
        std::uint32_t word = 0;
        for (std::size_t j = 0; j < 4 && i + j < bubbleSeed.size(); ++j)
          word |= static_cast<std::uint32_t>(static_cast<unsigned char>(bubbleSeed[i + j])) << (8 * j);
        seedWords.push_back(word);
      }
      if (seedWords.empty()) seedWords.push_back(0);
      std::seed_seq seedSequence(seedWords.begin(), seedWords.end());
      std::mt19937_64 generator(seedSequence);
      std::uniform_real_distribution<FloatType> uniform(0, side);

      std::vector<BubbleCenter> centers;
      centers.reserve(static_cast<size_t>(Nb));
      constexpr std::size_t maxAttempts = 10000000;
      std::size_t attempts = 0;
      const FloatType minimumSeparation = 2 * Rc;
      if (Nb > 1 && minimumSeparation >= side * std::sqrt(FloatType(3)) / 2)
        throw(RunParametersInconsistent(
            "The requested center separation is too large to place two bubbles in the periodic box."));

      while (centers.size() < static_cast<size_t>(Nb)) {
        if (++attempts > maxAttempts) {
          throw(RunParametersInconsistent(
              "Could not place all vacuum bubbles after " + std::to_string(maxAttempts) +
              " attempts. Reduce Nb or Rc, increase lSide, or use a different bubble_seed."));
        }

        BubbleCenter candidate{{uniform(generator), uniform(generator), uniform(generator)}};
        bool accepted = true;
        for (const auto &center : centers) {
          if (centerDistance(candidate, center, side) <= minimumSeparation) {
            accepted = false;
            break;
          }
        }
        if (accepted) centers.push_back(candidate);
      }
      return centers;
    }

    template <typename ToolBox> void initializeBubbleProfile(ToolBox toolBox, FloatType side)
    {
      const auto centers = makeBubbleCenters(side);
      const FloatType minimumSeparation = 2 * Rc;
      FloatType measuredMinimumSeparation = side * std::sqrt(FloatType(3)) / 2;
      for (std::size_t i = 0; i < centers.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j)
          measuredMinimumSeparation =
              std::min(measuredMinimumSeparation, centerDistance(centers[i], centers[j], side));
      }

      // Outside Rc + 20 l0, 1 - tanh((r - Rc)/l0) is zero at double
      // precision.  Binning the centers and evaluating only this compact
      // numerical support avoids an O(N^3 Nb) initialization at production
      // resolution while retaining the stated tanh profile.  The bin width is
      // chosen at least as large as the support whenever the box permits it,
      // keeping the neighbor stencil small and fully periodic.
      const FloatType profileCutoff = Rc + 20 * l0;
      const std::size_t nCellsPerDimension =
          std::max<std::size_t>(1, static_cast<std::size_t>(std::floor(side / profileCutoff)));
      const FloatType cellWidth = side / nCellsPerDimension;
      const std::size_t cellRadius =
          static_cast<std::size_t>(std::ceil(profileCutoff / cellWidth));
      const std::size_t nCells = nCellsPerDimension * nCellsPerDimension * nCellsPerDimension;
      auto cellIndex = [nCellsPerDimension](std::size_t x, std::size_t y, std::size_t z) {
        return (x * nCellsPerDimension + y) * nCellsPerDimension + z;
      };
      auto periodicCell = [nCellsPerDimension](ptrdiff_t index) {
        const ptrdiff_t n = static_cast<ptrdiff_t>(nCellsPerDimension);
        index %= n;
        return static_cast<std::size_t>(index < 0 ? index + n : index);
      };

      std::vector<std::vector<std::size_t>> centerBins(nCells);
      for (std::size_t i = 0; i < centers.size(); ++i) {
        const auto &center = centers[i];
        const auto cx = std::min(nCellsPerDimension - 1, static_cast<std::size_t>(center.x[0] / cellWidth));
        const auto cy = std::min(nCellsPerDimension - 1, static_cast<std::size_t>(center.x[1] / cellWidth));
        const auto cz = std::min(nCellsPerDimension - 1, static_cast<std::size_t>(center.x[2] / cellWidth));
        centerBins[cellIndex(cx, cy, cz)].push_back(i);
      }

      // Precompute unique neighboring cells once.  The uniqueness matters for
      // small boxes where periodic wrapping maps several offsets to one cell.
      std::vector<std::vector<std::size_t>> nearbyBins(nCells);
      for (std::size_t cx = 0; cx < nCellsPerDimension; ++cx) {
        for (std::size_t cy = 0; cy < nCellsPerDimension; ++cy) {
          for (std::size_t cz = 0; cz < nCellsPerDimension; ++cz) {
            auto &bins = nearbyBins[cellIndex(cx, cy, cz)];
            for (ptrdiff_t ix = -static_cast<ptrdiff_t>(cellRadius); ix <= static_cast<ptrdiff_t>(cellRadius); ++ix) {
              for (ptrdiff_t iy = -static_cast<ptrdiff_t>(cellRadius); iy <= static_cast<ptrdiff_t>(cellRadius); ++iy) {
                for (ptrdiff_t iz = -static_cast<ptrdiff_t>(cellRadius);
                     iz <= static_cast<ptrdiff_t>(cellRadius); ++iz) {
                  const auto candidate = cellIndex(periodicCell(static_cast<ptrdiff_t>(cx) + ix),
                                                   periodicCell(static_cast<ptrdiff_t>(cy) + iy),
                                                   periodicCell(static_cast<ptrdiff_t>(cz) + iz));
                  if (std::find(bins.begin(), bins.end(), candidate) == bins.end()) bins.push_back(candidate);
                }
              }
            }
          }
        }
      }

      auto &iterator = toolBox->itX();
      iterator.begin();
      for (; iterator.end(); ++iterator) {
        const auto &latticeCoordinate = iterator.getVec();
        FloatType point[3];
        for (int i = 0; i < 3; ++i) {
          point[i] = std::fmod(latticeCoordinate[i] * dx + side / 2, side);
          if (point[i] < 0) point[i] += side;
        }

        const auto px = std::min(nCellsPerDimension - 1, static_cast<std::size_t>(point[0] / cellWidth));
        const auto py = std::min(nCellsPerDimension - 1, static_cast<std::size_t>(point[1] / cellWidth));
        const auto pz = std::min(nCellsPerDimension - 1, static_cast<std::size_t>(point[2] / cellWidth));

        FloatType fieldSquared = 0;
        for (const auto bin : nearbyBins[cellIndex(px, py, pz)]) {
          for (const auto centerIndex : centerBins[bin]) {
            const auto &center = centers[centerIndex];
            FloatType distanceSquared = 0;
            for (int i = 0; i < 3; ++i) {
              const FloatType separation = minimumImageDistance(point[i], center.x[i], side);
              distanceSquared += separation * separation;
            }
            if (distanceSquared >= profileCutoff * profileCutoff) continue;
            const FloatType radius = std::sqrt(distanceSquared);
            const FloatType profile = phiB / 2 * (1 - std::tanh((radius - Rc) / l0));
            fieldSquared += profile * profile;
          }
        }
        fldS(0_c).getSet(iterator()) = std::sqrt(fieldSquared) / fStar;
        piS(0_c).getSet(iterator()) = 0;
      }
      fldS(0_c).setGhostsAreStale();
      piS(0_c).setGhostsAreStale();

      const FloatType rStar = side / std::cbrt(static_cast<FloatType>(Nb));
      const FloatType gammaStar = rStar / (2 * Rc);
      if (toolBox->amIRoot()) {
        say << "vacuum_bubbles: L = " << side << ", R* = " << rStar << ", gamma* = " << gammaStar
            << ", minimum center separation = " << measuredMinimumSeparation << " (required > " << minimumSeparation
            << ")\n";
      }
    }

  public:
    MODELNAME(ParameterParser &parser, RunParameters<FloatType> &runPar, auto toolBox)
        : Model<MODELNAME>(parser, runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL))
    {
      if (runPar.expansion)
        throw(RunParametersInconsistent("vacuum_bubbles requires expansion=false."));

      // The model writes its own inhomogeneous profile.  The normal scalar
      // initializer must therefore take the homogeneous path, which leaves
      // this profile intact and does not add Gaussian fluctuations.
      const auto &params = parser.getParams();
      const auto icType = params.count("ICtype_S") != 0 ? parser.get<std::string>("ICtype_S") : std::string();
      if (icType != "Homogeneous" && icType != "4")
        throw(RunParametersInconsistent("vacuum_bubbles requires ICtype_S=Homogeneous."));

      delta = getModelParameter(parser, "delta", FloatType(-1.632));
      lambda = getModelParameter(parser, "lambda", FloatType(0.5));
      phiB = getModelParameter(parser, "phi_b", FloatType(2.45));
      Rc = getModelParameter(parser, "Rc", FloatType(7.15));
      l0 = getModelParameter(parser, "l0", FloatType(1.71));
      Nb = getModelParameter(parser, "Nb", 512);
      bubbleSeed = getModelParameter(parser, "bubble_seed", std::string("vacuum-bubbles-row3"));

      if (lambda <= 0) throw(RunParametersInconsistent("lambda must be positive for vacuum_bubbles."));
      if (phiB <= 0) throw(RunParametersInconsistent("phi_b must be positive for vacuum_bubbles."));
      if (Rc <= 0) throw(RunParametersInconsistent("Rc must be positive for vacuum_bubbles."));
      if (l0 <= 0) throw(RunParametersInconsistent("l0 must be positive for vacuum_bubbles."));
      if (runPar.lSide <= 0) throw(RunParametersInconsistent("lSide must be positive for vacuum_bubbles."));

      // M=1 and alpha=0.  Keeping fStar=1 leaves the scalar spectrum in M
      // units, so the output is converted to P_phi/phi_b^2 in the README.
      fldS0(0_c) = 0;
      piS0(0_c) = 0;
      alpha = 0;
      fStar = 1;
      omegaStar = 1;

      // Compute the homogeneous mass and potential before writing the
      // inhomogeneous bubble profile; the helper temporarily uses the origin.
      setInitialPotentialAndMassesFromPotential();
      initializeBubbleProfile(toolBox, runPar.lSide);
    }

    auto potentialTerms(Tag<0>) const
    {
      return FloatType(0.5) * pow<2>(fldS(0_c)) + delta / 3 * pow<3>(fldS(0_c)) +
             lambda / 4 * pow<4>(fldS(0_c));
    }

    auto potDeriv(Tag<0>) const
    {
      return fldS(0_c) + delta * pow<2>(fldS(0_c)) + lambda * pow<3>(fldS(0_c));
    }

    auto potDeriv2(Tag<0>) const
    {
      return 1 + 2 * delta * fldS(0_c) + 3 * lambda * pow<2>(fldS(0_c));
    }
  };
} // namespace TempLat

#endif // VACUUM_BUBBLES_H
