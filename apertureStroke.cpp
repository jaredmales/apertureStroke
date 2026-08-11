#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <mx/ao/analysis/aoSystem.hpp>
#include <mx/ao/sim/turbAtmosphere.hpp>
#include <mx/app/application.hpp>
#include <mx/improc/imageMasks.hpp>
#include <mx/improc/milkImage.hpp>
#include <mx/ioutils/fits/fitsFile.hpp>
#include <mx/math/constants.hpp>
#include <mx/sigproc/basisUtils2D.hpp>

#include "basisFitters.hpp"
#include "strokeUtils.hpp"

namespace pss = aperture_stroke;

template<typename _realT>
class apertureStroke : public mx::app::application
{
  public:
    using realT = _realT;
    using imageT = mx::improc::eigenImage<realT>;
    using cubeT = mx::improc::eigenCube<realT>;
    using aoSystemT = mx::AO::analysis::aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>>;
    using turbulenceT = mx::AO::sim::turbAtmosphere<aoSystemT, mx::verbose::d>;

    std::string pupilFile;
    std::string pupilName;
    int pupilArraySize {182};
    realT pupilDiameterPixels {-1};
    realT pupilDiameterMeters {25.4};
    realT pupilCentralObscuration {0};
    realT pupilThreshold {0.5};

    int nTrials {10000};
    int simulationThreads {0};
    int turbulenceOversize {3};
    realT wavelength {0.8e-6};
    realT seeingArcsec {-1};
    realT seeingWavelength {500e-9};
    realT outerScale {-1};
    int subharmonicLevel {1};
    bool outerSubHarmonics {true};
    bool psdSubtractPiston {true};
    bool psdSubtractTipTilt {false};

    std::string basisType {"zernike"};
    std::string basisFile;
    std::string basisName;
    std::string prefixFile;
    std::string prefixName {"prefix"};
    int nModes {-1};
    std::string modeCutoffs {"default"};
    std::string fitMethod {"pinv"};
    realT pinvAlpha {1e-3};
    realT pinvMaxCondition {1e6};

    int maxFourierM {-1};
    realT histogramMinimum {0};
    realT histogramMaximum {10};
    realT histogramBinWidth {0.025};

    std::string outputDirectory {"output"};
    std::string outputLabel;
    bool writePupil {true};
    bool writeBasis {false};
    bool writePhaseScreens {true};
    bool printTiming {false};
    bool configError {false};

    apertureStroke()
    {
        m_helpWidth = 150;
        m_helpLOColWidth = 31;
        m_helpCFColWidth = 31;
    }

    void setupConfig() override
    {
        config.add("pupilFile", "", "pupil.file", mx::app::argType::Required,
                   "pupil", "file", false, "string",
                   "Path to a FITS pupil mask. If empty, a circular pupil is constructed.");
        config.add("pupilName", "", "pupil.name", mx::app::argType::Required,
                   "pupil", "name", false, "string",
                   "Short pupil name used in output filenames.");
        config.add("pupilArraySize", "", "pupil.arraySize", mx::app::argType::Required,
                   "pupil", "arraySize", false, "int",
                   "Square array size for a constructed circular pupil.");
        config.add("pupilDiameterPixels", "", "pupil.diameterPixels", mx::app::argType::Required,
                   "pupil", "diameterPixels", false, "real",
                   "Outer pupil diameter in pixels; defaults to the array size.");
        config.add("pupilDiameterMeters", "", "pupil.diameterMeters", mx::app::argType::Required,
                   "pupil", "diameterMeters", false, "real",
                   "Physical outer pupil diameter in meters.");
        config.add("pupilCentralObscuration", "", "pupil.centralObscuration", mx::app::argType::Required,
                   "pupil", "centralObscuration", false, "real",
                   "Central-obscuration diameter divided by outer diameter for a circular pupil.");
        config.add("pupilThreshold", "", "pupil.threshold", mx::app::argType::Required,
                   "pupil", "threshold", false, "real",
                   "FITS pupil values at or above this value are included in the binary mask.");

        config.add("nTrials", "", "simulation.trials", mx::app::argType::Required,
                   "simulation", "trials", false, "int", "Number of independent phase screens.");
        config.add("simulationThreads", "", "simulation.threads", mx::app::argType::Required,
                   "simulation", "threads", false, "int",
                   "Independent trial workers; 0 uses the OpenMP worker limit.");
        config.add("turbulenceOversize", "", "simulation.oversize", mx::app::argType::Required,
                   "simulation", "oversize", false, "int",
                   "Integer grid oversize; 1 adds one wavefront width on each side.");
        config.add("wavelength", "", "simulation.wavelength", mx::app::argType::Required,
                   "simulation", "wavelength", false, "real",
                   "Phase-screen wavelength in meters.");
        config.add("seeingArcsec", "", "atmosphere.seeing", mx::app::argType::Required,
                   "atmosphere", "seeing", false, "real",
                   "Seeing FWHM in arcseconds; a negative value keeps the mxlib LCO default.");
        config.add("seeingWavelength", "", "atmosphere.seeingWavelength", mx::app::argType::Required,
                   "atmosphere", "seeingWavelength", false, "real",
                   "Wavelength in meters at which seeing is specified.");
        config.add("outerScale", "", "atmosphere.outerScale", mx::app::argType::Required,
                   "atmosphere", "outerScale", false, "real",
                   "Outer scale in meters; 0 is Kolmogorov/infinite and a negative value keeps the mxlib default.");
        config.add("subharmonicLevel", "", "atmosphere.subharmonicLevel", mx::app::argType::Required,
                   "atmosphere", "subharmonicLevel", false, "int",
                   "Subharmonic level passed to turbAtmosphere.");
        config.add("outerSubHarmonics", "", "atmosphere.outerSubHarmonics", mx::app::argType::Required,
                   "atmosphere", "outerSubHarmonics", false, "bool",
                   "Whether turbAtmosphere includes outer subharmonics.");
        config.add("psdSubtractPiston", "", "atmosphere.psdSubtractPiston", mx::app::argType::Required,
                   "atmosphere", "psdSubtractPiston", false, "bool",
                   "Apply mxlib's piston-removal transfer function to the generated turbulence PSD.");
        config.add("psdSubtractTipTilt", "", "atmosphere.psdSubtractTipTilt", mx::app::argType::Required,
                   "atmosphere", "psdSubtractTipTilt", false, "bool",
                   "Apply mxlib's tip/tilt-removal transfer function to the generated turbulence PSD.");

        config.add("basisType", "", "basis.type", mx::app::argType::Required,
                   "basis", "type", false, "string", "Primary basis type: zernike or file.");
        config.add("basisFile", "", "basis.file", mx::app::argType::Required,
                   "basis", "file", false, "string", "FITS cube used when basis.type=file.");
        config.add("basisName", "", "basis.name", mx::app::argType::Required,
                   "basis", "name", false, "string", "Primary basis name used in output filenames.");
        config.add("prefixFile", "", "basis.prefixFile", mx::app::argType::Required,
                   "basis", "prefixFile", false, "string",
                   "Optional FITS cube prepended to every primary-mode cutoff, such as segment PTT modes.");
        config.add("prefixName", "", "basis.prefixName", mx::app::argType::Required,
                   "basis", "prefixName", false, "string", "Prefix basis name used in output filenames.");
        config.add("nModes", "", "basis.modes", mx::app::argType::Required,
                   "basis", "modes", false, "int",
                   "Number of primary modes; negative means 1000 Zernikes or all modes from a FITS cube.");
        config.add("modeCutoffs", "", "basis.cutoffs", mx::app::argType::Required,
                   "basis", "cutoffs", false, "string",
                   "Primary-mode cutoffs: single, default, or a comma-separated list.");
        config.add("fitMethod", "", "basis.fit", mx::app::argType::Required,
                   "basis", "fit", false, "string", "Fit method: sequential, projection, lsq, or pinv.");
        config.add("pinvAlpha", "", "basis.pinvAlpha", mx::app::argType::Required,
                   "basis", "pinvAlpha", false, "real", "Pseudo-inverse Tikhonov regularization.");
        config.add("pinvMaxCondition", "", "basis.pinvMaxCondition", mx::app::argType::Required,
                   "basis", "pinvMaxCondition", false, "real", "Pseudo-inverse maximum condition number.");

        config.add("maxFourierM", "", "analysis.maxFourierM", mx::app::argType::Required,
                   "analysis", "maxFourierM", false, "int",
                   "Maximum modified-Fourier m index; -1 disables Fourier measurements and 0 uses half the pupil array size.");
        config.add("histogramMinimum", "", "analysis.histogramMinimum", mx::app::argType::Required,
                   "analysis", "histogramMinimum", false, "real", "Histogram minimum in microns surface.");
        config.add("histogramMaximum", "", "analysis.histogramMaximum", mx::app::argType::Required,
                   "analysis", "histogramMaximum", false, "real", "Histogram maximum in microns surface.");
        config.add("histogramBinWidth", "", "analysis.histogramBinWidth", mx::app::argType::Required,
                   "analysis", "histogramBinWidth", false, "real", "Histogram bin width in microns surface.");

        config.add("outputDirectory", "", "output.directory", mx::app::argType::Required,
                   "output", "directory", false, "string", "Output directory, created if necessary.");
        config.add("outputLabel", "", "output.label", mx::app::argType::Required,
                   "output", "label", false, "string", "Label used in output filenames.");
        config.add("writePupil", "", "output.writePupil", mx::app::argType::Required,
                   "output", "writePupil", false, "bool", "Write the binary pupil as FITS.");
        config.add("writeBasis", "", "output.writeBasis", mx::app::argType::Required,
                   "output", "writeBasis", false, "bool", "Write the normalized combined basis as FITS.");
        config.add("writePhaseScreens", "", "output.writePhaseScreens", mx::app::argType::Required,
                   "output", "writePhaseScreens", false, "bool", "Write the first input and final residual screens.");
        config.add("printTiming", "", "output.printTiming", mx::app::argType::Required,
                   "output", "printTiming", false, "bool", "Print aggregate wall-clock timing by simulation stage.");
    }

    void loadConfig() override
    {
        config(pupilFile, "pupilFile");
        config(pupilName, "pupilName");
        config(pupilArraySize, "pupilArraySize");
        config(pupilDiameterPixels, "pupilDiameterPixels");
        config(pupilDiameterMeters, "pupilDiameterMeters");
        config(pupilCentralObscuration, "pupilCentralObscuration");
        config(pupilThreshold, "pupilThreshold");

        config(nTrials, "nTrials");
        config(simulationThreads, "simulationThreads");
        config(turbulenceOversize, "turbulenceOversize");
        config(wavelength, "wavelength");
        config(seeingArcsec, "seeingArcsec");
        config(seeingWavelength, "seeingWavelength");
        config(outerScale, "outerScale");
        config(subharmonicLevel, "subharmonicLevel");
        config(outerSubHarmonics, "outerSubHarmonics");
        config(psdSubtractPiston, "psdSubtractPiston");
        config(psdSubtractTipTilt, "psdSubtractTipTilt");

        config(basisType, "basisType");
        config(basisFile, "basisFile");
        config(basisName, "basisName");
        config(prefixFile, "prefixFile");
        config(prefixName, "prefixName");
        config(nModes, "nModes");
        config(modeCutoffs, "modeCutoffs");
        config(fitMethod, "fitMethod");
        config(pinvAlpha, "pinvAlpha");
        config(pinvMaxCondition, "pinvMaxCondition");

        config(maxFourierM, "maxFourierM");
        config(histogramMinimum, "histogramMinimum");
        config(histogramMaximum, "histogramMaximum");
        config(histogramBinWidth, "histogramBinWidth");

        config(outputDirectory, "outputDirectory");
        config(outputLabel, "outputLabel");
        config(writePupil, "writePupil");
        config(writeBasis, "writeBasis");
        config(writePhaseScreens, "writePhaseScreens");
        config(printTiming, "printTiming");
    }

    void checkConfig() override
    {
        std::vector<std::string> unusedSections;
        config.unusedSections(unusedSections);
        if(unusedSections.empty())
        {
            return;
        }

        configError = true;
        std::cerr << "unrecognized command-line option or configuration key";
        bool haveNamedSection = false;
        for(const std::string & section : unusedSections)
        {
            if(!section.empty())
            {
                std::cerr << (haveNamedSection ? ", " : " in section(s): ") << section;
                haveNamedSection = true;
            }
        }
        std::cerr << '\n';
    }

    int execute() override
    {
        if(configError)
        {
            return -1;
        }

        if(validateScalarConfig() < 0)
        {
            return -1;
        }

        mx::fits::fitsFile<realT> fits;
        imageT pupil;
        if(loadPupil(pupil, fits) < 0)
        {
            return -1;
        }

        int wavefrontSize = pupil.rows();
        if(pupilDiameterPixels <= 0)
        {
            pupilDiameterPixels = std::min(pupil.rows(), pupil.cols());
        }
        if(pupilDiameterPixels > std::min(pupil.rows(), pupil.cols()))
        {
            std::cerr << "pupil diameter in pixels exceeds the pupil array\n";
            return -1;
        }

        if(pupilName.empty())
        {
            pupilName = pupilFile.empty() ? "circle" : std::filesystem::path(pupilFile).stem().string();
        }
        if(outputLabel.empty())
        {
            outputLabel = pupilName;
        }
        pupilName = pss::cleanLabel(pupilName);
        outputLabel = pss::cleanLabel(outputLabel);
        prefixName = pss::cleanLabel(prefixName);

        std::filesystem::path outputPath(outputDirectory);
        std::error_code fsError;
        std::filesystem::create_directories(outputPath, fsError);
        if(fsError || !std::filesystem::is_directory(outputPath))
        {
            std::cerr << "could not create output directory " << outputPath
                      << ": " << fsError.message() << '\n';
            return -1;
        }

        if(writePupil && writeFits(fits, outputPath / std::format("pupil_{}.fits", outputLabel), pupil) < 0)
        {
            return -1;
        }

        cubeT primaryModes;
        if(loadPrimaryBasis(primaryModes, pupil, fits) < 0)
        {
            return -1;
        }

        cubeT prefixModes;
        if(!prefixFile.empty() && readCube(prefixModes, prefixFile, pupil, fits) < 0)
        {
            return -1;
        }

        int nPrefixModes = prefixFile.empty() ? 0 : prefixModes.planes();
        cubeT combinedModes;
        combineModes(combinedModes, prefixModes, primaryModes, pupil.rows(), pupil.cols());
        int nFitModes = combinedModes.planes();

        std::vector<int> cutoffs;
        if(pss::parseModeCutoffs(cutoffs, modeCutoffs, nModes) < 0)
        {
            return -1;
        }

        std::vector<int> fitModeCounts(cutoffs.size());
        for(size_t n = 0; n < cutoffs.size(); ++n)
        {
            fitModeCounts[n] = nPrefixModes + cutoffs[n];
        }

        if(nFitModes > 0 && mx::sigproc::basisNormalize(combinedModes, pupil) < 0)
        {
            std::cerr << "error normalizing basis modes on the pupil\n";
            return -1;
        }

        if(writeBasis && nFitModes > 0 &&
           writeFits(fits, outputPath / std::format("basis_{}.fits", outputLabel), combinedModes) < 0)
        {
            return -1;
        }

        if(nTrials == 0)
        {
            std::cerr << "nTrials = 0; wrote setup products only to " << outputPath << '\n';
            return 0;
        }

        pss::BasisFit fitType;
        if(pss::parseBasisFit(fitType, fitMethod) < 0)
        {
            return -1;
        }

        pss::BasisLeastSquaresFitter<realT> leastSquaresFitter;
        pss::BasisPseudoInverseFitter<realT> pseudoInverseFitter;
        pss::BasisProjectionFitter<realT> projectionFitter;
        if(fitType == pss::BasisFit::projection &&
           projectionFitter.setup(combinedModes, pupil, nFitModes) < 0)
        {
            return -1;
        }
        if(fitType == pss::BasisFit::leastSquares &&
           leastSquaresFitter.setup(combinedModes, pupil, nFitModes, fitModeCounts) < 0)
        {
            return -1;
        }
        if(fitType == pss::BasisFit::pseudoInverse &&
           pseudoInverseFitter.setup(combinedModes,
                                     pupil,
                                     nFitModes,
                                     pinvAlpha,
                                     pinvMaxCondition,
                                     fitModeCounts) < 0)
        {
            return -1;
        }

        aoSystemT aoSystem;
        aoSystem.atm.loadLCO();
        aoSystem.lam_wfs(wavelength);
        aoSystem.lam_sci(wavelength);
        realT gridDiameterMeters = pupilDiameterMeters * wavefrontSize / pupilDiameterPixels;
        aoSystem.D(gridDiameterMeters);

        if(seeingArcsec > 0)
        {
            realT seeingRadians = seeingArcsec * mx::math::pi<realT>() / (180.0 * 3600.0);
            realT r0 = 0.98 * seeingWavelength / seeingRadians;
            aoSystem.atm.r_0(r0, seeingWavelength);
        }
        if(outerScale >= 0)
        {
            aoSystem.atm.L_0(std::vector<realT>(aoSystem.atm.n_layers(), outerScale));
        }

        realT singleLayerR0 = aoSystem.atm.r_0(seeingWavelength);
        realT effectiveSeeing = aoSystem.atm.fwhm0(seeingWavelength) *
                                180.0 * 3600.0 / mx::math::pi<realT>();
        realT effectiveOuterScale = aoSystem.atm.L_0(0);
        aoSystem.atm.setSingleLayer(singleLayerR0,
                                    seeingWavelength,
                                    effectiveOuterScale,
                                    0,
                                    0,
                                    0,
                                    0);
        aoSystem.psd.subPiston(psdSubtractPiston);
        aoSystem.psd.subTipTilt(psdSubtractTipTilt);

        int oversizePixels = turbulenceOversize * wavefrontSize;
        int nWorkers = simulationThreads;
#ifdef _OPENMP
        if(nWorkers == 0)
        {
            nWorkers = omp_get_max_threads();
        }
#else
        if(nWorkers == 0)
        {
            nWorkers = 1;
        }
        if(nWorkers > 1)
        {
            std::cerr << "simulation.threads requires an OpenMP build\n";
            return -1;
        }
#endif
        nWorkers = std::min(nWorkers, nTrials);

        auto setupTurbulence = [&](turbulenceT & turbulence)
        {
            turbulence.retain(true);
            turbulence.forceGen(true);
            turbulence.outerSubHarmonics(outerSubHarmonics);
            turbulence.setup(wavefrontSize,
                             oversizePixels,
                             &aoSystem,
                             static_cast<uint32_t>(subharmonicLevel));
            turbulence.setLayers(wavefrontSize + 2 * oversizePixels);
        };

        std::vector<std::unique_ptr<turbulenceT>> turbulences;
        std::vector<std::unique_ptr<mx::improc::milkImage<realT>>> wavefronts;
        turbulences.reserve(nWorkers);
        wavefronts.reserve(nWorkers);
        for(int worker = 0; worker < nWorkers; ++worker)
        {
            turbulences.emplace_back(std::make_unique<turbulenceT>());
            setupTurbulence(*turbulences.back());
            wavefronts.emplace_back(std::make_unique<mx::improc::milkImage<realT>>(
                std::format("phase_{}", worker), wavefrontSize, wavefrontSize));
        }

        std::vector<imageT> fourierModes;
        pss::FourierAmplitudeMeasurer<realT> fourierMeasurer;
        if(maxFourierM >= 0)
        {
            if(maxFourierM == 0)
            {
                maxFourierM = std::min(pupil.rows(), pupil.cols()) / 2;
            }

            if(pss::makeP2VNormalizedFourierModes(fourierModes, pupil, maxFourierM) < 0)
            {
                return -1;
            }

            if(fourierMeasurer.setup(fourierModes, pupil) < 0)
            {
                return -1;
            }
        }

        std::vector<realT> basisModeP2V = pss::basisP2Vs(combinedModes, pupil);
        std::vector<CutoffStats> statistics(cutoffs.size());
        for(size_t c = 0; c < statistics.size(); ++c)
        {
            statistics[c].nModes = cutoffs[c];
            statistics[c].nFitModes = fitModeCounts[c];
            statistics[c].p2v.resize(nTrials);
            statistics[c].p2pDifference.resize(nTrials);
            statistics[c].modeP2V.resize(statistics[c].nFitModes,
                                         std::vector<float>(nTrials));
            statistics[c].fourierP2V.resize(fourierModes.size(),
                                            std::vector<float>(nTrials));
        }

        std::string atmosphereTag = effectiveOuterScale == 0
            ? "L0inf"
            : "L0" + pss::scalarTag(effectiveOuterScale) + "m";
        atmosphereTag += std::format("_sh{}{}",
                                     subharmonicLevel,
                                     outerSubHarmonics ? "" : "_nooutersh");
        atmosphereTag += std::format("_1layer_psdP{}_psdTT{}",
                                     psdSubtractPiston ? 1 : 0,
                                     psdSubtractTipTilt ? 1 : 0);
        std::string commonTag = std::format("{}trials_{}over_seeing{:.3f}arcsec_{}",
                                            nTrials,
                                            turbulenceOversize,
                                            effectiveSeeing,
                                            atmosphereTag);

        writeRunParameters(outputPath / std::format("run_{}_{}.dat", outputLabel, commonTag),
                           pupil,
                           nPrefixModes,
                           fitType,
                           gridDiameterMeters,
                           nWorkers,
                           effectiveSeeing,
                           effectiveOuterScale,
                           cutoffs);

        std::cerr << "pupil: " << pupilName
                  << " pixels: " << pupil.sum()
                  << " aperture: " << pupilDiameterMeters << " m / " << pupilDiameterPixels << " px"
                  << " grid diameter: " << gridDiameterMeters << " m\n"
                  << "basis: " << basisName << " primary modes: " << nModes
                  << " prefix modes: " << nPrefixModes
                  << " fit: " << pss::fitName(fitType) << '\n'
                  << "seeing: " << effectiveSeeing << " arcsec at " << seeingWavelength << " m"
                  << " L0: " << effectiveOuterScale
                  << " shLevel: " << subharmonicLevel
                  << " outerSubHarmonics: " << std::boolalpha << outerSubHarmonics
                  << " psdSubtractPiston: " << psdSubtractPiston
                  << " psdSubtractTipTilt: " << psdSubtractTipTilt << std::noboolalpha << '\n'
                  << "wfSz: " << turbulences.front()->wfSz()
                  << " layers: " << turbulences.front()->nLayers()
                  << " oversize: " << turbulenceOversize
                  << " buffSz: " << turbulences.front()->buffSz()
                  << " layerSz: " << wavefrontSize + 2 * oversizePixels
                  << " workers: " << nWorkers
                  << " output: " << outputPath << '\n';

        imageT firstInputScreen;
        imageT firstResidualScreen;
        realT phaseToSurfaceMicrons = wavelength * 1e6 / mx::math::two_pi<realT>() * 0.5;
        double generationSeconds = 0;
        double projectionSeconds = 0;
        double subtractionSeconds = 0;
        double metricSeconds = 0;
        double fourierSeconds = 0;
        std::atomic<bool> trialError {false};

        // Every worker owns mutable turbulence and image state; statistics are indexed by trial.
        #pragma omp parallel for num_threads(nWorkers) schedule(static) reduction(+ : generationSeconds, projectionSeconds, subtractionSeconds, metricSeconds, fourierSeconds)
        for(int trial = 0; trial < nTrials; ++trial)
        {
            if(trialError.load(std::memory_order_relaxed))
            {
                continue;
            }

            int worker = 0;
#ifdef _OPENMP
            worker = omp_get_thread_num();
#endif
            turbulenceT & turbulence = *turbulences[worker];
            mx::improc::milkImage<realT> & wavefront = *wavefronts[worker];
            imageT inputScreen;
            imageT residual;
            auto stageStart = std::chrono::steady_clock::now();
            turbulence.genLayers();
            turbulence.shift(wavefront, 0.0);
            inputScreen = wavefront() * pupil;
            generationSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - stageStart).count();

            std::vector<realT> projectedAmplitudes;
            int projectedModesSubtracted = 0;
            stageStart = std::chrono::steady_clock::now();
            if(fitType == pss::BasisFit::projection &&
               projectionFitter.project(projectedAmplitudes, residual, inputScreen) < 0)
            {
                trialError.store(true, std::memory_order_relaxed);
                continue;
            }
            projectionSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - stageStart).count();

            for(size_t c = 0; c < statistics.size(); ++c)
            {
                CutoffStats & stats = statistics[c];
                std::vector<realT> amplitudes;

                stageStart = std::chrono::steady_clock::now();
                if(fitType == pss::BasisFit::projection)
                {
                    if(projectionFitter.subtractRange(residual,
                                                      projectedAmplitudes,
                                                      projectedModesSubtracted,
                                                      stats.nFitModes) < 0)
                    {
                        trialError.store(true, std::memory_order_relaxed);
                        break;
                    }
                    projectedModesSubtracted = stats.nFitModes;
                    amplitudes.assign(projectedAmplitudes.begin(),
                                      projectedAmplitudes.begin() + stats.nFitModes);
                }
                else if(fitType == pss::BasisFit::leastSquares)
                {
                    residual = inputScreen;
                    if(leastSquaresFitter.subtract(amplitudes, residual, c) < 0)
                    {
                        trialError.store(true, std::memory_order_relaxed);
                        break;
                    }
                }
                else if(fitType == pss::BasisFit::pseudoInverse)
                {
                    residual = inputScreen;
                    if(pseudoInverseFitter.subtract(amplitudes, residual, c) < 0)
                    {
                        trialError.store(true, std::memory_order_relaxed);
                        break;
                    }
                }
                else if(stats.nFitModes > 0)
                {
                    residual = inputScreen;
                    if(mx::sigproc::basisAmplitudes(amplitudes,
                                                    residual,
                                                    combinedModes,
                                                    pupil,
                                                    true,
                                                    1,
                                                    stats.nFitModes) < 0)
                    {
                        std::cerr << "sequential basis fit failed\n";
                        trialError.store(true, std::memory_order_relaxed);
                        break;
                    }
                }
                else
                {
                    residual = inputScreen;
                    pss::subtractPupilMean(residual, pupil);
                }
                subtractionSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - stageStart).count();

                for(size_t n = 0; n < amplitudes.size(); ++n)
                {
                    stats.modeP2V[n][trial] = std::abs(amplitudes[n]) *
                                              basisModeP2V[n] *
                                              phaseToSurfaceMicrons;
                }

                residual *= pupil;
                stageStart = std::chrono::steady_clock::now();
                stats.p2v[trial] = pss::pupilP2V(residual, pupil) * phaseToSurfaceMicrons;
                stats.p2pDifference[trial] = pss::maxAbsPixelDiff(residual, pupil) * phaseToSurfaceMicrons;
                metricSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - stageStart).count();

                if(!fourierModes.empty())
                {
                    stageStart = std::chrono::steady_clock::now();
                    std::vector<realT> fourierAmplitudes;
                    if(fourierMeasurer.measure(fourierAmplitudes, residual) < 0)
                    {
                        trialError.store(true, std::memory_order_relaxed);
                        break;
                    }
                    for(size_t n = 0; n < fourierAmplitudes.size(); ++n)
                    {
                        stats.fourierP2V[n][trial] = fourierAmplitudes[n] * phaseToSurfaceMicrons;
                    }
                    fourierSeconds +=
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - stageStart).count();
                }
            }

            if(trial == 0 && writePhaseScreens && !trialError.load(std::memory_order_relaxed))
            {
                firstInputScreen = inputScreen;
                firstResidualScreen = residual;
            }
        }

        if(trialError.load(std::memory_order_relaxed))
        {
            return -1;
        }

        if(writePhaseScreens &&
           writeFits(fits,
                     outputPath / std::format("phase_input_{}_{}.fits", outputLabel, commonTag),
                     firstInputScreen) < 0)
        {
            return -1;
        }
        if(writePhaseScreens &&
           writeFits(fits,
                     outputPath / std::format("phase_residual_{}_{}modes_{}.fits",
                                              outputLabel,
                                              statistics.back().nModes,
                                              commonTag),
                     firstResidualScreen) < 0)
        {
                return -1;
        }

        if(printTiming)
        {
            std::cerr << "timing s/trial: generation " << generationSeconds / nTrials
                      << " projection " << projectionSeconds / nTrials
                      << " subtraction " << subtractionSeconds / nTrials
                      << " metrics " << metricSeconds / nTrials
                      << " fourier " << fourierSeconds / nTrials << '\n';
        }

        for(const CutoffStats & stats : statistics)
        {
            std::string tag = std::format("{}modes_{}", stats.nModes, commonTag);
            pss::writeHistogram<realT>(outputPath / std::format("p2v_{}_{}.dat", outputLabel, tag),
                                       stats.p2v,
                                       histogramMinimum,
                                       histogramMaximum,
                                       histogramBinWidth);
            pss::writeHistogram<realT>(outputPath / std::format("p2pdiff_{}_{}.dat", outputLabel, tag),
                                       stats.p2pDifference,
                                       histogramMinimum,
                                       histogramMaximum,
                                       histogramBinWidth);
            pss::writeResidualStats<realT>(outputPath / std::format("p2v_stats_{}_{}.dat", outputLabel, tag),
                                           stats.p2v,
                                           stats.p2pDifference);
            if(!stats.fourierP2V.empty())
            {
                pss::writeFourierStats<realT>(outputPath / std::format("fourier_p2v_{}_{}.dat", outputLabel, tag),
                                              stats.fourierP2V);
            }

            if(nPrefixModes > 0)
            {
                pss::writeModeStats<realT>(outputPath / std::format("{}_p2v_{}_{}.dat",
                                                                    prefixName,
                                                                    outputLabel,
                                                                    tag),
                                           stats.modeP2V,
                                           0,
                                           nPrefixModes);
            }
            pss::writeModeStats<realT>(outputPath / std::format("{}_p2v_{}_{}.dat",
                                                                basisName,
                                                                outputLabel,
                                                                tag),
                                       stats.modeP2V,
                                       nPrefixModes,
                                       stats.nModes);
        }

        return 0;
    }

  private:
    struct CutoffStats
    {
        int nModes {0};
        int nFitModes {0};
        std::vector<float> p2v;
        std::vector<float> p2pDifference;
        std::vector<std::vector<float>> modeP2V;
        std::vector<std::vector<float>> fourierP2V;
    };

    int validateScalarConfig() const
    {
        if(nTrials < 0 || simulationThreads < 0 || turbulenceOversize < 0 || subharmonicLevel < 0)
        {
            std::cerr << "trials, threads, oversize, and subharmonic level must be non-negative\n";
            return -1;
        }
        if(pupilDiameterMeters <= 0 || wavelength <= 0 || seeingWavelength <= 0)
        {
            std::cerr << "pupil diameter and wavelengths must be positive\n";
            return -1;
        }
        if(seeingArcsec == 0)
        {
            std::cerr << "seeing must be positive or negative to use the mxlib default\n";
            return -1;
        }
        if(pupilCentralObscuration < 0 || pupilCentralObscuration >= 1)
        {
            std::cerr << "central obscuration must be in [0,1)\n";
            return -1;
        }
        if(!pupilFile.empty() && pupilCentralObscuration != 0)
        {
            std::cerr << "central obscuration applies only to a constructed circular pupil\n";
            return -1;
        }
        if(maxFourierM < -1 || histogramMaximum <= histogramMinimum || histogramBinWidth <= 0)
        {
            std::cerr << "invalid analysis limits\n";
            return -1;
        }
        return 0;
    }

    int loadPupil(imageT & pupil, mx::fits::fitsFile<realT> & fits)
    {
        if(!pupilFile.empty())
        {
            imageT input;
            if(fits.read(input, pupilFile) != mx::error_t::noerror)
            {
                std::cerr << "could not read pupil FITS file: " << pupilFile << '\n';
                return -1;
            }
            if(input.rows() != input.cols())
            {
                std::cerr << "pupil FITS image must be square\n";
                return -1;
            }

            pupil.resize(input.rows(), input.cols());
            for(int cc = 0; cc < input.cols(); ++cc)
            {
                for(int rr = 0; rr < input.rows(); ++rr)
                {
                    pupil(rr, cc) = input(rr, cc) >= pupilThreshold ? 1 : 0;
                }
            }
        }
        else
        {
            if(pupilArraySize <= 0)
            {
                std::cerr << "pupil array size must be positive\n";
                return -1;
            }
            if(pupilDiameterPixels <= 0)
            {
                pupilDiameterPixels = pupilArraySize;
            }
            if(pupilDiameterPixels > pupilArraySize)
            {
                std::cerr << "circular pupil diameter exceeds its array\n";
                return -1;
            }

            if(pss::makeCircularPupil(pupil,
                                      pupilArraySize,
                                      pupilDiameterPixels,
                                      pupilCentralObscuration) < 0)
            {
                std::cerr << "error constructing circular pupil\n";
                return -1;
            }
        }

        if(pupil.sum() == 0)
        {
            std::cerr << "pupil has no included pixels\n";
            return -1;
        }
        return 0;
    }

    int loadPrimaryBasis(cubeT & modes,
                         const imageT & pupil,
                         mx::fits::fitsFile<realT> & fits)
    {
        if(basisType == "zernike" || basisType == "zernikes")
        {
            if(nModes < 0)
            {
                nModes = 1000;
            }
            if(nModes < 0)
            {
                return -1;
            }
            basisName = basisName.empty() ? "zernike" : pss::cleanLabel(basisName);
            return pss::makeZernikeBasis(modes,
                                         nModes,
                                         pupil.rows(),
                                         pupil.cols(),
                                         static_cast<realT>(0.5 * pupilDiameterPixels));
        }

        if(basisType == "file" || basisType == "fits")
        {
            if(basisFile.empty())
            {
                std::cerr << "basis.file is required when basis.type=file\n";
                return -1;
            }
            cubeT input;
            if(readCube(input, basisFile, pupil, fits) < 0)
            {
                return -1;
            }
            if(nModes < 0)
            {
                nModes = input.planes();
            }
            if(nModes < 0 || nModes > input.planes())
            {
                std::cerr << "requested " << nModes << " primary modes but "
                          << input.planes() << " are available\n";
                return -1;
            }

            modes.resize(pupil.rows(), pupil.cols(), nModes);
            for(int n = 0; n < nModes; ++n)
            {
                modes.image(n) = input.image(n);
            }
            basisName = basisName.empty()
                ? pss::cleanLabel(std::filesystem::path(basisFile).stem().string())
                : pss::cleanLabel(basisName);
            return 0;
        }

        std::cerr << "unknown basis type: " << basisType << " (expected zernike or file)\n";
        return -1;
    }

    int readCube(cubeT & cube,
                 const std::string & path,
                 const imageT & pupil,
                 mx::fits::fitsFile<realT> & fits) const
    {
        if(fits.read(cube, path) != mx::error_t::noerror)
        {
            std::cerr << "could not read basis FITS cube: " << path << '\n';
            return -1;
        }
        if(cube.rows() != pupil.rows() || cube.cols() != pupil.cols())
        {
            std::cerr << "basis cube " << path << " is " << cube.rows() << 'x' << cube.cols()
                      << " but the pupil is " << pupil.rows() << 'x' << pupil.cols() << '\n';
            return -1;
        }
        return 0;
    }

    void combineModes(cubeT & combined,
                      const cubeT & prefix,
                      const cubeT & primary,
                      int rows,
                      int cols) const
    {
        int nPrefix = prefixFile.empty() ? 0 : prefix.planes();
        combined.resize(rows, cols, nPrefix + primary.planes());
        for(int n = 0; n < nPrefix; ++n)
        {
            combined.image(n) = prefix.image(n);
        }
        for(int n = 0; n < primary.planes(); ++n)
        {
            combined.image(nPrefix + n) = primary.image(n);
        }
    }

    template<typename fitsDataT>
    int writeFits(mx::fits::fitsFile<realT> & fits,
                  const std::filesystem::path & path,
                  const fitsDataT & data) const
    {
        if(fits.write(path.string(), data) != mx::error_t::noerror)
        {
            std::cerr << "could not write FITS file: " << path << '\n';
            return -1;
        }
        return 0;
    }

    void writeRunParameters(const std::filesystem::path & path,
                            const imageT & pupil,
                            int nPrefixModes,
                            pss::BasisFit fitType,
                            realT gridDiameterMeters,
                            int nWorkers,
                            realT effectiveSeeing,
                            realT effectiveOuterScale,
                            const std::vector<int> & cutoffs) const
    {
        std::ofstream output(path);
        output << "# apertureStroke run parameters\n"
               << "pupil.name " << pupilName << '\n'
               << "pupil.file " << (pupilFile.empty() ? "circle" : pupilFile) << '\n'
               << "pupil.rows " << pupil.rows() << '\n'
               << "pupil.pixels " << pupil.sum() << '\n'
               << "pupil.diameterPixels " << pupilDiameterPixels << '\n'
               << "pupil.diameterMeters " << pupilDiameterMeters << '\n'
               << "pupil.centralObscuration " << pupilCentralObscuration << '\n'
               << "simulation.gridDiameterMeters " << gridDiameterMeters << '\n'
               << "simulation.wavelength " << wavelength << '\n'
               << "simulation.trials " << nTrials << '\n'
               << "simulation.threads " << simulationThreads << '\n'
               << "simulation.workers " << nWorkers << '\n'
               << "simulation.oversize " << turbulenceOversize << '\n'
               << "atmosphere.seeing " << effectiveSeeing << '\n'
               << "atmosphere.seeingWavelength " << seeingWavelength << '\n'
               << "atmosphere.outerScale " << effectiveOuterScale << '\n'
               << "atmosphere.layers 1\n"
               << "atmosphere.subharmonicLevel " << subharmonicLevel << '\n'
               << "atmosphere.outerSubHarmonics " << std::boolalpha << outerSubHarmonics << '\n'
               << "atmosphere.psdSubtractPiston " << psdSubtractPiston << '\n'
               << "atmosphere.psdSubtractTipTilt " << psdSubtractTipTilt << '\n'
               << "basis.type " << basisType << '\n'
               << "basis.name " << basisName << '\n'
               << "basis.primaryModes " << nModes << '\n'
               << "basis.prefixModes " << nPrefixModes << '\n'
               << "basis.fit " << pss::fitName(fitType) << '\n'
               << "basis.pinvAlpha " << pinvAlpha << '\n'
               << "basis.pinvMaxCondition " << pinvMaxCondition << '\n'
               << "basis.cutoffs";
        for(int cutoff : cutoffs)
        {
            output << ' ' << cutoff;
        }
        output << '\n';
    }
};

int main(int argc, char ** argv)
{
    apertureStroke<float> app;
    return app.main(argc, argv);
}
