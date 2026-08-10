#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include <mx/app/application.hpp>
#include <mx/improc/eigenCube.hpp>
#include <mx/improc/eigenImage.hpp>
#include <mx/math/constants.hpp>
#include <mx/sigproc/basisUtils2D.hpp>
#include <mx/sigproc/zernike.hpp>

#include "strokeUtils.hpp"

namespace ast = aperture_stroke;

template<typename _realT>
class zernikeTheory : public mx::app::application
{
  public:
    using realT = _realT;
    using imageT = mx::improc::eigenImage<realT>;
    using cubeT = mx::improc::eigenCube<realT>;

    std::string pupilFile;
    std::string pupilName {"circle"};
    int pupilArraySize {182};
    realT pupilDiameterPixels {-1};
    realT pupilDiameterMeters {25.4};
    realT pupilCentralObscuration {0};

    realT wavelength {0.8e-6};
    realT seeingArcsec {0.79};
    realT seeingWavelength {500e-9};
    realT outerScale {0};

    std::string basisType {"zernike"};
    int nModes {1000};

    realT r0 {-1};
    realT r0Wavelength {-1};
    std::string outputFile;

    std::string outputDirectory {"output"};
    std::string outputLabel;

    zernikeTheory()
    {
        m_helpWidth = 150;
        m_helpLOColWidth = 31;
        m_helpCFColWidth = 31;
    }

    void setupConfig() override
    {
        config.add("pupilFile", "", "pupil.file", mx::app::argType::Required,
                   "pupil", "file", false, "string",
                   "Must be empty: this utility implements unobstructed circular-aperture theory only.");
        config.add("pupilName", "", "pupil.name", mx::app::argType::Required,
                   "pupil", "name", false, "string", "Pupil name used in output metadata.");
        config.add("pupilArraySize", "", "pupil.arraySize", mx::app::argType::Required,
                   "pupil", "arraySize", false, "int", "Square pupil-array size in pixels.");
        config.add("pupilDiameterPixels", "", "pupil.diameterPixels", mx::app::argType::Required,
                   "pupil", "diameterPixels", false, "real",
                   "Circular-pupil diameter in pixels; negative uses the array size.");
        config.add("pupilDiameterMeters", "", "pupil.diameterMeters", mx::app::argType::Required,
                   "pupil", "diameterMeters", false, "real", "Circular-aperture diameter in meters.");
        config.add("pupilCentralObscuration", "", "pupil.centralObscuration", mx::app::argType::Required,
                   "pupil", "centralObscuration", false, "real",
                   "Must be zero: Noll theory here is for an unobstructed circle.");

        config.add("wavelength", "", "simulation.wavelength", mx::app::argType::Required,
                   "simulation", "wavelength", false, "real",
                   "Phase wavelength in meters, used for r0 scaling and conversion to surface.");
        config.add("seeingArcsec", "", "atmosphere.seeing", mx::app::argType::Required,
                   "atmosphere", "seeing", false, "real",
                   "Seeing FWHM in arcseconds, used when theory.r0 is negative.");
        config.add("seeingWavelength", "", "atmosphere.seeingWavelength", mx::app::argType::Required,
                   "atmosphere", "seeingWavelength", false, "real",
                   "Wavelength in meters at which seeing is specified.");
        config.add("outerScale", "", "atmosphere.outerScale", mx::app::argType::Required,
                   "atmosphere", "outerScale", false, "real",
                   "Must be zero: this utility implements Kolmogorov theory.");

        config.add("basisType", "", "basis.type", mx::app::argType::Required,
                   "basis", "type", false, "string", "Must be zernike.");
        config.add("nModes", "", "basis.modes", mx::app::argType::Required,
                   "basis", "modes", false, "int", "Number of Zernike modes beginning with Noll j=2.");

        config.add("r0", "", "theory.r0", mx::app::argType::Required,
                   "theory", "r0", false, "real",
                   "Optional Fried parameter in meters; negative derives it from seeing.");
        config.add("r0Wavelength", "", "theory.r0Wavelength", mx::app::argType::Required,
                   "theory", "r0Wavelength", false, "real",
                   "Wavelength for theory.r0; negative uses atmosphere.seeingWavelength.");
        config.add("outputFile", "", "theory.outputFile", mx::app::argType::Required,
                   "theory", "outputFile", false, "string",
                   "Output table path; empty uses output.directory/zernike_theory_<label>.dat and '-' uses stdout.");

        config.add("outputDirectory", "", "output.directory", mx::app::argType::Required,
                   "output", "directory", false, "string", "Default output directory, created if necessary.");
        config.add("outputLabel", "", "output.label", mx::app::argType::Required,
                   "output", "label", false, "string", "Label used in the default output filename.");
    }

    void loadConfig() override
    {
        config(pupilFile, "pupilFile");
        config(pupilName, "pupilName");
        config(pupilArraySize, "pupilArraySize");
        config(pupilDiameterPixels, "pupilDiameterPixels");
        config(pupilDiameterMeters, "pupilDiameterMeters");
        config(pupilCentralObscuration, "pupilCentralObscuration");

        config(wavelength, "wavelength");
        config(seeingArcsec, "seeingArcsec");
        config(seeingWavelength, "seeingWavelength");
        config(outerScale, "outerScale");

        config(basisType, "basisType");
        config(nModes, "nModes");

        config(r0, "r0");
        config(r0Wavelength, "r0Wavelength");
        config(outputFile, "outputFile");

        config(outputDirectory, "outputDirectory");
        config(outputLabel, "outputLabel");
    }

    int execute() override
    {
        if(validateConfig() < 0)
        {
            return -1;
        }

        if(pupilDiameterPixels <= 0)
        {
            pupilDiameterPixels = pupilArraySize;
        }

        imageT pupil;
        if(ast::makeCircularPupil(pupil,
                                  pupilArraySize,
                                  pupilDiameterPixels,
                                  static_cast<realT>(0)) < 0)
        {
            std::cerr << "error constructing circular pupil\n";
            return -1;
        }

        cubeT modes;
        if(ast::makeZernikeBasis(modes,
                                 nModes,
                                 pupil.rows(),
                                 pupil.cols(),
                                 static_cast<realT>(0.5 * pupilDiameterPixels)) < 0)
        {
            std::cerr << "error constructing Zernike basis\n";
            return -1;
        }
        if(nModes > 0 && mx::sigproc::basisNormalize(modes, pupil) < 0)
        {
            std::cerr << "error normalizing Zernike basis\n";
            return -1;
        }

        std::vector<realT> p2vOverRms = ast::basisP2Vs(modes, pupil);

        realT referenceWavelength = r0Wavelength > 0 ? r0Wavelength : seeingWavelength;
        realT referenceR0 = r0;
        if(referenceR0 < 0)
        {
            realT seeingRadians = seeingArcsec * mx::math::pi<realT>() / (180.0 * 3600.0);
            referenceR0 = 0.98 * seeingWavelength / seeingRadians;
            referenceWavelength = seeingWavelength;
        }

        realT r0AtWavelength = referenceR0 *
                               std::pow(wavelength / referenceWavelength,
                                        static_cast<realT>(6.0 / 5.0));
        realT phaseToSurfaceMicrons = wavelength * 1e6 /
                                      mx::math::two_pi<realT>() *
                                      static_cast<realT>(0.5);

        std::unique_ptr<std::ofstream> fileStream;
        std::ostream * output = &std::cout;
        std::filesystem::path outputPath;
        if(outputFile != "-" && outputFile != "stdout")
        {
            if(outputLabel.empty())
            {
                outputLabel = pupilName.empty() ? "circle" : pupilName;
            }
            outputLabel = ast::cleanLabel(outputLabel);
            outputPath = outputFile.empty()
                ? std::filesystem::path(outputDirectory) /
                  std::format("zernike_theory_{}.dat", outputLabel)
                : std::filesystem::path(outputFile);

            if(outputPath.has_parent_path())
            {
                std::error_code fsError;
                std::filesystem::create_directories(outputPath.parent_path(), fsError);
                if(fsError)
                {
                    std::cerr << "could not create output directory: " << fsError.message() << '\n';
                    return -1;
                }
            }

            fileStream = std::make_unique<std::ofstream>(outputPath);
            if(!*fileStream)
            {
                std::cerr << "could not open theory output: " << outputPath << '\n';
                return -1;
            }
            output = fileStream.get();
        }

        if(writePrediction(*output,
                           p2vOverRms,
                           referenceR0,
                           referenceWavelength,
                           r0AtWavelength,
                           phaseToSurfaceMicrons,
                           pupil.sum()) < 0)
        {
            return -1;
        }

        if(fileStream)
        {
            std::cerr << "wrote " << outputPath << '\n';
        }
        return 0;
    }

  private:
    int validateConfig() const
    {
        if(!pupilFile.empty())
        {
            std::cerr << "zernikeTheory accepts only a constructed circular pupil\n";
            return -1;
        }
        if(pupilCentralObscuration != 0)
        {
            std::cerr << "zernikeTheory requires an unobstructed circular pupil\n";
            return -1;
        }
        if(pupilArraySize <= 0 || pupilDiameterMeters <= 0 ||
           pupilDiameterPixels > pupilArraySize || wavelength <= 0 ||
           seeingWavelength <= 0 || nModes < 0)
        {
            std::cerr << "invalid pupil, wavelength, or mode configuration\n";
            return -1;
        }
        if(outerScale != 0)
        {
            std::cerr << "zernikeTheory implements only L0=0 Kolmogorov theory\n";
            return -1;
        }
        if(basisType != "zernike" && basisType != "zernikes")
        {
            std::cerr << "zernikeTheory requires basis.type=zernike\n";
            return -1;
        }
        if(r0 == 0 || (r0 < 0 && seeingArcsec <= 0) ||
           (r0 > 0 && r0Wavelength == 0))
        {
            std::cerr << "specify positive seeing or positive theory.r0\n";
            return -1;
        }
        return 0;
    }

    int writePrediction(std::ostream & output,
                        const std::vector<realT> & p2vOverRms,
                        realT referenceR0,
                        realT referenceWavelength,
                        realT r0AtWavelength,
                        realT phaseToSurfaceMicrons,
                        realT pupilPixels) const
    {
        output << "# unobstructed circular-aperture Kolmogorov prediction\n"
               << "# diameter_m " << pupilDiameterMeters << '\n'
               << "# pupil_array_size " << pupilArraySize << '\n'
               << "# pupil_diameter_pixels " << pupilDiameterPixels << '\n'
               << "# pupil_pixels " << pupilPixels << '\n'
               << "# wavelength_m " << wavelength << '\n'
               << "# r0_reference_m " << referenceR0 << '\n'
               << "# r0_reference_wavelength_m " << referenceWavelength << '\n'
               << "# r0_at_wavelength_m " << r0AtWavelength << '\n'
               << "# phase_to_surface_um_per_rad " << phaseToSurfaceMicrons << '\n'
               << "# surface values are wavefront/2\n"
               << "# mean and rms predict the absolute modal P2V and match columns 2 and 3 of zernike_p2v tables\n"
               << "# rms is the true root-mean-square; stddev is RMS about the mean for a half-normal distribution\n"
               << "# p2v_over_rms is measured on the discrete mxlib circular pupil after basisNormalize\n"
               << "#mode mean rms stddev noll_j variance_coeff variance_rad2 coeff_rms_rad p2v_over_rms\n";

        for(size_t n = 0; n < p2vOverRms.size(); ++n)
        {
            unsigned nollIndex = static_cast<unsigned>(n + 2);
            realT varianceCoefficient =
                mx::sigproc::zernikeModeDOCDiffKolmogorov<realT>(nollIndex);
            if(varianceCoefficient <= 0)
            {
                std::cerr << "non-positive Kolmogorov variance coefficient for Noll j="
                          << nollIndex << '\n';
                return -1;
            }

            realT variance = varianceCoefficient *
                             std::pow(pupilDiameterMeters / r0AtWavelength,
                                      mx::math::five_thirds<realT>());
            realT coefficientRms = std::sqrt(variance);
            realT p2vSigma = coefficientRms * p2vOverRms[n] * phaseToSurfaceMicrons;
            realT mean = std::sqrt(static_cast<realT>(2) / mx::math::pi<realT>()) * p2vSigma;
            realT stddev = std::sqrt(static_cast<realT>(1) -
                                     static_cast<realT>(2) / mx::math::pi<realT>()) * p2vSigma;

            output << n + 1 << ' ' << mean << ' ' << p2vSigma << ' ' << stddev << ' '
                   << nollIndex << ' ' << varianceCoefficient << ' ' << variance << ' '
                   << coefficientRms << ' ' << p2vOverRms[n] << '\n';
        }

        return output ? 0 : -1;
    }
};

int main(int argc, char ** argv)
{
    zernikeTheory<float> app;
    return app.main(argc, argv);
}
