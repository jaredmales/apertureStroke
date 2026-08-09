#ifndef apertureStroke_strokeUtils_hpp
#define apertureStroke_strokeUtils_hpp

#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <mx/improc/eigenCube.hpp>
#include <mx/improc/eigenImage.hpp>
#include <mx/improc/imageMasks.hpp>
#include <mx/math/histogramUniform.hpp>
#include <mx/math/vectorUtils.hpp>
#include <mx/sigproc/fourierModes.hpp>
#include <mx/sigproc/zernike.hpp>

namespace aperture_stroke
{

template<typename realT>
int makeCircularPupil(mx::improc::eigenImage<realT> & pupil,
                      int arraySize,
                      realT diameterPixels,
                      realT centralObscuration = 0)
{
    if(arraySize <= 0 || diameterPixels <= 0 || diameterPixels > arraySize)
    {
        return -1;
    }
    if(centralObscuration < 0 || centralObscuration >= 1)
    {
        return -1;
    }

    pupil.resize(arraySize, arraySize);
    pupil.setZero();

    realT outerRadius = 0.5 * diameterPixels;
    mx::improc::maskCircle(pupil, outerRadius, static_cast<realT>(1));
    if(centralObscuration > 0)
    {
        mx::improc::maskCircle(pupil,
                               centralObscuration * outerRadius,
                               static_cast<realT>(0));
    }

    return pupil.sum() > 0 ? 0 : -1;
}

enum class BasisFit
{
    sequential,
    leastSquares,
    pseudoInverse
};

inline int parseBasisFit(BasisFit & fit, const std::string & value)
{
    if(value == "sequential" || value == "seq")
    {
        fit = BasisFit::sequential;
        return 0;
    }

    if(value == "lsq" || value == "least_squares" || value == "least-squares")
    {
        fit = BasisFit::leastSquares;
        return 0;
    }

    if(value == "pinv" || value == "pseudo_inverse" || value == "pseudo-inverse")
    {
        fit = BasisFit::pseudoInverse;
        return 0;
    }

    std::cerr << "unknown basis fit method: " << value
              << " (expected sequential, lsq, or pinv)\n";
    return -1;
}

inline std::string fitName(BasisFit fit)
{
    switch(fit)
    {
        case BasisFit::sequential:
            return "sequential";
        case BasisFit::leastSquares:
            return "lsq";
        case BasisFit::pseudoInverse:
            return "pinv";
    }

    return "unknown";
}

inline std::vector<int> defaultModeCutoffs(int maxModes)
{
    std::vector<int> cutoffs {0};

    if(maxModes >= 20)
    {
        cutoffs.push_back(20);
    }

    for(int n = 100; n <= maxModes; n += 100)
    {
        cutoffs.push_back(n);
    }

    if(cutoffs.back() != maxModes)
    {
        cutoffs.push_back(maxModes);
    }

    std::sort(cutoffs.begin(), cutoffs.end());
    cutoffs.erase(std::unique(cutoffs.begin(), cutoffs.end()), cutoffs.end());
    return cutoffs;
}

inline int parseModeCutoffs(std::vector<int> & cutoffs,
                            const std::string & value,
                            int maxModes)
{
    cutoffs.clear();

    if(value == "single")
    {
        cutoffs.push_back(maxModes);
        return 0;
    }

    if(value == "sweep" || value == "all" || value == "default")
    {
        cutoffs = defaultModeCutoffs(maxModes);
        return 0;
    }

    size_t start = 0;
    while(start <= value.size())
    {
        size_t end = value.find(',', start);
        std::string token = value.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if(token.empty())
        {
            std::cerr << "empty mode cutoff in list: " << value << '\n';
            return -1;
        }

        int cutoff = 0;
        try
        {
            cutoff = std::stoi(token);
        }
        catch(const std::exception &)
        {
            std::cerr << "invalid mode cutoff: " << token << '\n';
            return -1;
        }

        if(cutoff < 0 || cutoff > maxModes)
        {
            std::cerr << "mode cutoff " << cutoff << " is outside [0," << maxModes << "]\n";
            return -1;
        }
        cutoffs.push_back(cutoff);

        if(end == std::string::npos)
        {
            break;
        }
        start = end + 1;
    }

    std::sort(cutoffs.begin(), cutoffs.end());
    cutoffs.erase(std::unique(cutoffs.begin(), cutoffs.end()), cutoffs.end());
    return cutoffs.empty() ? -1 : 0;
}

inline std::string scalarTag(double value)
{
    std::string tag = std::format("{:.4g}", value);
    for(char & c : tag)
    {
        if(c == '.')
        {
            c = 'p';
        }
        else if(c == '-' || c == '+')
        {
            c = c == '-' ? 'm' : 'p';
        }
    }
    return tag;
}

inline std::string cleanLabel(std::string value)
{
    for(char & c : value)
    {
        if(!std::isalnum(static_cast<unsigned char>(c)) && c != '-' && c != '_')
        {
            c = '_';
        }
    }
    return value.empty() ? "pupil" : value;
}

template<typename realT>
realT maxAbsPixelDiff(const mx::improc::eigenImage<realT> & image,
                      const mx::improc::eigenImage<realT> & mask)
{
    realT maxDiff = 0;
    for(int cc = 0; cc < image.cols(); ++cc)
    {
        for(int rr = 0; rr < image.rows(); ++rr)
        {
            if(mask(rr, cc) == 0)
            {
                continue;
            }

            if(rr + 1 < image.rows() && mask(rr + 1, cc) != 0)
            {
                maxDiff = std::max(maxDiff, std::abs(image(rr + 1, cc) - image(rr, cc)));
            }
            if(cc + 1 < image.cols() && mask(rr, cc + 1) != 0)
            {
                maxDiff = std::max(maxDiff, std::abs(image(rr, cc + 1) - image(rr, cc)));
            }
        }
    }
    return maxDiff;
}

template<typename imageT, typename realT>
realT pupilP2V(const imageT & image,
               const mx::improc::eigenImage<realT> & mask)
{
    realT minValue = 0;
    realT maxValue = 0;
    bool havePixel = false;

    for(int cc = 0; cc < image.cols(); ++cc)
    {
        for(int rr = 0; rr < image.rows(); ++rr)
        {
            if(mask(rr, cc) == 0)
            {
                continue;
            }

            if(!havePixel)
            {
                minValue = image(rr, cc);
                maxValue = image(rr, cc);
                havePixel = true;
            }
            else
            {
                minValue = std::min(minValue, static_cast<realT>(image(rr, cc)));
                maxValue = std::max(maxValue, static_cast<realT>(image(rr, cc)));
            }
        }
    }

    return havePixel ? maxValue - minValue : 0;
}

template<typename realT>
void subtractPupilMean(mx::improc::eigenImage<realT> & image,
                       const mx::improc::eigenImage<realT> & mask)
{
    realT sum = 0;
    size_t nPixels = 0;
    for(int cc = 0; cc < image.cols(); ++cc)
    {
        for(int rr = 0; rr < image.rows(); ++rr)
        {
            if(mask(rr, cc) != 0)
            {
                sum += image(rr, cc);
                ++nPixels;
            }
        }
    }

    if(nPixels == 0)
    {
        return;
    }

    realT mean = sum / static_cast<realT>(nPixels);
    for(int cc = 0; cc < image.cols(); ++cc)
    {
        for(int rr = 0; rr < image.rows(); ++rr)
        {
            if(mask(rr, cc) != 0)
            {
                image(rr, cc) -= mean;
            }
        }
    }
}

template<typename realT>
std::vector<realT> basisP2Vs(const mx::improc::eigenCube<realT> & modes,
                             const mx::improc::eigenImage<realT> & mask)
{
    std::vector<realT> p2vs(modes.planes());
    for(int n = 0; n < modes.planes(); ++n)
    {
        p2vs[n] = pupilP2V(modes.image(n), mask);
    }
    return p2vs;
}

template<typename realT>
int makeZernikeBasis(mx::improc::eigenCube<realT> & modes,
                     int nModes,
                     int rows,
                     int cols,
                     realT radiusPixels)
{
    modes.resize(rows, cols, nModes);
    if(nModes == 0)
    {
        return 0;
    }
    return mx::sigproc::zernikeBasis<mx::improc::eigenCube<realT>, double>(modes, radiusPixels);
}

template<typename realT>
int makeP2VNormalizedFourierModes(std::vector<mx::improc::eigenImage<realT>> & modes,
                                  const mx::improc::eigenImage<realT> & mask,
                                  int maxM)
{
    modes.resize(maxM);
    for(int m = 1; m <= maxM; ++m)
    {
        mx::improc::eigenImage<realT> & mode = modes[m - 1];
        mode.resize(mask.rows(), mask.cols());

        mx::improc::eigenMap<realT> modeMap(mode.data(), mode.rows(), mode.cols());
        if(mx::sigproc::makeModifiedFourierMode(modeMap, static_cast<realT>(m), static_cast<realT>(0), 1) < 0)
        {
            return -1;
        }

        mode *= mask;
        realT modeP2V = pupilP2V(mode, mask);
        if(modeP2V == 0)
        {
            std::cerr << "Fourier mode m=" << m << " has zero P2V on the pupil\n";
            return -1;
        }
        mode /= modeP2V;
    }
    return 0;
}

template<typename realT>
realT fourierModeP2VAmplitude(const mx::improc::eigenImage<realT> & image,
                              const mx::improc::eigenImage<realT> & mode,
                              const mx::improc::eigenImage<realT> & mask)
{
    realT numerator = 0;
    realT denominator = 0;
    for(int cc = 0; cc < image.cols(); ++cc)
    {
        for(int rr = 0; rr < image.rows(); ++rr)
        {
            if(mask(rr, cc) != 0)
            {
                numerator += image(rr, cc) * mode(rr, cc);
                denominator += mode(rr, cc) * mode(rr, cc);
            }
        }
    }
    return denominator == 0 ? 0 : std::abs(numerator / denominator);
}

template<typename realT>
void writeModeStats(const std::filesystem::path & fileName,
                    const std::vector<std::vector<float>> & values,
                    size_t first,
                    size_t count)
{
    std::ofstream output(fileName);
    output << "#mode mean rms max\n";

    for(size_t n = 0; n < count; ++n)
    {
        const std::vector<float> & samples = values[first + n];
        realT mean = mx::math::vectorMean(samples);
        realT rms = samples.size() > 1 ? std::sqrt(mx::math::vectorVariance(samples, mean)) : 0;
        realT maxValue = *std::max_element(samples.begin(), samples.end());
        output << n + 1 << ' ' << mean << ' ' << rms << ' ' << maxValue << '\n';
    }
}

template<typename realT>
void writeFourierStats(const std::filesystem::path & fileName,
                       const std::vector<std::vector<float>> & values)
{
    std::ofstream output(fileName);
    output << "#m mean rms max\n";

    for(size_t n = 0; n < values.size(); ++n)
    {
        realT mean = mx::math::vectorMean(values[n]);
        realT rms = values[n].size() > 1 ? std::sqrt(mx::math::vectorVariance(values[n], mean)) : 0;
        realT maxValue = *std::max_element(values[n].begin(), values[n].end());
        output << n + 1 << ' ' << mean << ' ' << rms << ' ' << maxValue << '\n';
    }
}

template<typename realT>
void writeResidualStats(const std::filesystem::path & fileName,
                        const std::vector<float> & p2v,
                        const std::vector<float> & p2pDifference)
{
    std::ofstream output(fileName);
    output << "#metric mean rms max\n";

    realT mean = mx::math::vectorMean(p2v);
    realT rms = p2v.size() > 1 ? std::sqrt(mx::math::vectorVariance(p2v, mean)) : 0;
    output << "p2v " << mean << ' ' << rms << ' '
           << *std::max_element(p2v.begin(), p2v.end()) << '\n';

    mean = mx::math::vectorMean(p2pDifference);
    rms = p2pDifference.size() > 1 ? std::sqrt(mx::math::vectorVariance(p2pDifference, mean)) : 0;
    output << "p2pdiff " << mean << ' ' << rms << ' '
           << *std::max_element(p2pDifference.begin(), p2pDifference.end()) << '\n';
}

template<typename realT>
void writeHistogram(const std::filesystem::path & fileName,
                    const std::vector<float> & values,
                    realT minimum,
                    realT maximum,
                    realT binWidth)
{
    std::ofstream output(fileName);
    realT rms = std::sqrt(mx::math::vectorVariance(values, 0));
    realT maxValue = *std::max_element(values.begin(), values.end());
    mx::math::histogramUniform<realT> histogram(minimum, maximum, binWidth, values, true);

    output << "#rms: " << rms << '\n';
    output << "#max: " << maxValue << '\n';
    for(int n = 0; n < histogram.bins(); ++n)
    {
        output << histogram.binMid(n) << ' ' << histogram.freq(n) << '\n';
    }
}

} // namespace aperture_stroke

#endif
