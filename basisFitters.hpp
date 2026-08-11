#ifndef apertureStroke_basisFitters_hpp
#define apertureStroke_basisFitters_hpp

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <mx/improc/eigenCube.hpp>
#include <mx/improc/eigenImage.hpp>
#include <mx/math/eigenLapack.hpp>

namespace aperture_stroke
{

template<typename realT>
class BasisProjectionFitter
{
  public:
    int setup(const mx::improc::eigenCube<realT> & modes,
              const mx::improc::eigenImage<realT> & mask,
              int nModes)
    {
        m_modes = &modes;
        m_coords.clear();
        m_modeNorms.assign(nModes, 0);
        m_weightSum = 0;

        if(nModes < 0 || nModes > modes.planes())
        {
            std::cerr << "modal projection requested " << nModes
                      << " modes but only " << modes.planes() << " are available\n";
            return -1;
        }
        if(modes.rows() != mask.rows() || modes.cols() != mask.cols())
        {
            std::cerr << "modal projection basis and pupil dimensions do not match\n";
            return -1;
        }

        for(int cc = 0; cc < mask.cols(); ++cc)
        {
            for(int rr = 0; rr < mask.rows(); ++rr)
            {
                if(mask(rr, cc) != 0)
                {
                    m_coords.push_back({rr, cc, mask(rr, cc)});
                    m_weightSum += mask(rr, cc);
                }
            }
        }

        if(m_coords.empty() || m_weightSum <= 0)
        {
            std::cerr << "modal projection has an empty pupil\n";
            return -1;
        }

        for(int n = 0; n < nModes; ++n)
        {
            double norm = 0;
            for(const Coordinate & coord : m_coords)
            {
                double modeValue = modes.image(n)(coord.row, coord.col);
                norm += coord.weight * modeValue * modeValue;
            }
            if(norm <= 0)
            {
                std::cerr << "modal projection mode " << n << " has zero norm on the pupil\n";
                return -1;
            }
            m_modeNorms[n] = norm;
        }

        return 0;
    }

    int project(std::vector<realT> & amplitudes,
                mx::improc::eigenImage<realT> & residual,
                const mx::improc::eigenImage<realT> & image) const
    {
        if(m_modes == nullptr || image.rows() != m_modes->rows() || image.cols() != m_modes->cols())
        {
            std::cerr << "modal projection image dimensions do not match the basis\n";
            return -1;
        }

        residual = image;
        double mean = 0;
        for(const Coordinate & coord : m_coords)
        {
            mean += coord.weight * image(coord.row, coord.col);
        }
        mean /= m_weightSum;

        for(const Coordinate & coord : m_coords)
        {
            residual(coord.row, coord.col) -= static_cast<realT>(mean);
        }

        amplitudes.assign(m_modeNorms.size(), 0);
        for(size_t n = 0; n < m_modeNorms.size(); ++n)
        {
            double innerProduct = 0;
            for(const Coordinate & coord : m_coords)
            {
                innerProduct += coord.weight *
                                residual(coord.row, coord.col) *
                                m_modes->image(static_cast<int>(n))(coord.row, coord.col);
            }
            amplitudes[n] = static_cast<realT>(innerProduct / m_modeNorms[n]);
        }

        return 0;
    }

    int subtractRange(mx::improc::eigenImage<realT> & residual,
                      const std::vector<realT> & amplitudes,
                      int firstMode,
                      int endMode) const
    {
        if(firstMode < 0 || endMode < firstMode ||
           static_cast<size_t>(endMode) > amplitudes.size() ||
           static_cast<size_t>(endMode) > m_modeNorms.size())
        {
            std::cerr << "invalid modal projection subtraction range ["
                      << firstMode << ',' << endMode << ")\n";
            return -1;
        }

        for(int n = firstMode; n < endMode; ++n)
        {
            for(const Coordinate & coord : m_coords)
            {
                residual(coord.row, coord.col) -=
                    amplitudes[n] * m_modes->image(n)(coord.row, coord.col);
            }
        }

        return 0;
    }

  private:
    struct Coordinate
    {
        int row;
        int col;
        realT weight;
    };

    const mx::improc::eigenCube<realT> * m_modes {nullptr};
    std::vector<Coordinate> m_coords;
    std::vector<double> m_modeNorms;
    double m_weightSum {0};
};

template<typename realT>
class FourierAmplitudeMeasurer
{
  public:
    int setup(const std::vector<mx::improc::eigenImage<realT>> & modes,
              const mx::improc::eigenImage<realT> & mask)
    {
        m_rows = mask.rows();
        m_cols = mask.cols();
        m_coords.clear();
        if(modes.empty())
        {
            m_basis.resize(0, 0);
            m_modeNorms.resize(0);
            return 0;
        }

        for(const auto & mode : modes)
        {
            if(mode.rows() != m_rows || mode.cols() != m_cols)
            {
                std::cerr << "Fourier mode and pupil dimensions do not match\n";
                return -1;
            }
        }

        for(int cc = 0; cc < m_cols; ++cc)
        {
            for(int rr = 0; rr < m_rows; ++rr)
            {
                if(mask(rr, cc) != 0)
                {
                    m_coords.push_back({rr, cc});
                }
            }
        }

        if(m_coords.empty())
        {
            std::cerr << "Fourier amplitude measurement has an empty pupil\n";
            return -1;
        }

        m_basis.resize(static_cast<int>(m_coords.size()), static_cast<int>(modes.size()));
        m_modeNorms.resize(static_cast<int>(modes.size()));
        for(size_t n = 0; n < modes.size(); ++n)
        {
            double norm = 0;
            for(size_t p = 0; p < m_coords.size(); ++p)
            {
                double modeValue = modes[n](m_coords[p].first, m_coords[p].second);
                m_basis(static_cast<int>(p), static_cast<int>(n)) = modeValue;
                norm += modeValue * modeValue;
            }
            if(norm <= 0)
            {
                std::cerr << "Fourier mode " << n << " has zero norm on the pupil\n";
                return -1;
            }
            m_modeNorms(static_cast<int>(n)) = norm;
        }

        return 0;
    }

    int measure(std::vector<realT> & amplitudes,
                const mx::improc::eigenImage<realT> & image) const
    {
        if(image.rows() != m_rows || image.cols() != m_cols)
        {
            std::cerr << "Fourier amplitude image dimensions do not match the pupil\n";
            return -1;
        }

        Eigen::VectorXd pupilValues(static_cast<int>(m_coords.size()));
        for(size_t p = 0; p < m_coords.size(); ++p)
        {
            pupilValues(static_cast<int>(p)) = image(m_coords[p].first, m_coords[p].second);
        }

        Eigen::VectorXd coefficients = m_basis.transpose() * pupilValues;
        coefficients.array() /= m_modeNorms.array();

        amplitudes.resize(static_cast<size_t>(coefficients.size()));
        for(int n = 0; n < coefficients.size(); ++n)
        {
            amplitudes[static_cast<size_t>(n)] = static_cast<realT>(std::abs(coefficients(n)));
        }

        return 0;
    }

  private:
    int m_rows {0};
    int m_cols {0};
    std::vector<std::pair<int, int>> m_coords;
    Eigen::MatrixXd m_basis;
    Eigen::VectorXd m_modeNorms;
};

template<typename realT>
class BasisLeastSquaresFitter
{
  public:
    int setup(const mx::improc::eigenCube<realT> & modes,
              const mx::improc::eigenImage<realT> & mask,
              int nModes,
              const std::vector<int> & fitModeCounts)
    {
        m_nModes = nModes;
        m_coords.clear();
        m_modeCounts = fitModeCounts;
        m_solvers.clear();

        for(int cc = 0; cc < mask.cols(); ++cc)
        {
            for(int rr = 0; rr < mask.rows(); ++rr)
            {
                if(mask(rr, cc) != 0)
                {
                    m_coords.push_back({rr, cc});
                }
            }
        }

        if(m_coords.empty())
        {
            std::cerr << "least-squares basis fit has an empty pupil\n";
            return -1;
        }

        if(m_nModes == 0)
        {
            m_solvers.resize(m_modeCounts.size());
            return 0;
        }

        m_basis.resize(static_cast<int>(m_coords.size()), m_nModes);
        for(int n = 0; n < m_nModes; ++n)
        {
            for(size_t p = 0; p < m_coords.size(); ++p)
            {
                m_basis(static_cast<int>(p), n) = modes.image(n)(m_coords[p].first, m_coords[p].second);
            }
        }

        m_solvers.reserve(m_modeCounts.size());
        for(int nModesFit : m_modeCounts)
        {
            if(nModesFit < 0 || nModesFit > m_nModes)
            {
                std::cerr << "least-squares fit requested " << nModesFit
                          << " modes but only " << m_nModes << " are available\n";
                return -1;
            }

            m_solvers.emplace_back();
            if(nModesFit == 0)
            {
                continue;
            }

            Eigen::MatrixXd gram =
                (m_basis.leftCols(nModesFit).transpose() * m_basis.leftCols(nModesFit)) /
                static_cast<double>(m_coords.size());
            m_solvers.back().compute(gram);
            if(m_solvers.back().info() != Eigen::Success)
            {
                std::cerr << "least-squares Gram factorization failed for "
                          << nModesFit << " modes\n";
                return -1;
            }
        }

        return 0;
    }

    int subtract(std::vector<realT> & amplitudes,
                 mx::improc::eigenImage<realT> & image,
                 size_t cutoffIndex) const
    {
        int nModesFit = m_modeCounts[cutoffIndex];
        amplitudes.assign(nModesFit, 0);

        Eigen::VectorXd pupilValues(static_cast<int>(m_coords.size()));
        double mean = 0;
        for(size_t p = 0; p < m_coords.size(); ++p)
        {
            pupilValues(static_cast<int>(p)) = image(m_coords[p].first, m_coords[p].second);
            mean += pupilValues(static_cast<int>(p));
        }

        mean /= static_cast<double>(m_coords.size());
        pupilValues.array() -= mean;

        if(nModesFit > 0)
        {
            Eigen::VectorXd rhs =
                (m_basis.leftCols(nModesFit).transpose() * pupilValues) /
                static_cast<double>(m_coords.size());
            Eigen::VectorXd coefficients = m_solvers[cutoffIndex].solve(rhs);
            if(m_solvers[cutoffIndex].info() != Eigen::Success)
            {
                std::cerr << "least-squares basis solve failed\n";
                return -1;
            }

            pupilValues.noalias() -= m_basis.leftCols(nModesFit) * coefficients;
            for(int n = 0; n < nModesFit; ++n)
            {
                amplitudes[n] = static_cast<realT>(coefficients(n));
            }
        }

        for(size_t p = 0; p < m_coords.size(); ++p)
        {
            image(m_coords[p].first, m_coords[p].second) =
                static_cast<realT>(pupilValues(static_cast<int>(p)));
        }

        return 0;
    }

  private:
    int m_nModes {0};
    std::vector<int> m_modeCounts;
    std::vector<std::pair<int, int>> m_coords;
    Eigen::MatrixXd m_basis;
    std::vector<Eigen::LDLT<Eigen::MatrixXd>> m_solvers;
};

template<typename realT>
class BasisPseudoInverseFitter
{
  public:
    int setup(const mx::improc::eigenCube<realT> & modes,
              const mx::improc::eigenImage<realT> & mask,
              int nModes,
              double alpha,
              double maxCondition,
              const std::vector<int> & fitModeCounts)
    {
        m_nModes = nModes;
        m_coords.clear();
        m_modeCounts = fitModeCounts;
        m_gramPInvs.clear();

        for(int cc = 0; cc < mask.cols(); ++cc)
        {
            for(int rr = 0; rr < mask.rows(); ++rr)
            {
                if(mask(rr, cc) != 0)
                {
                    m_coords.push_back({rr, cc});
                }
            }
        }

        if(m_coords.empty())
        {
            std::cerr << "pseudo-inverse basis fit has an empty pupil\n";
            return -1;
        }

        if(m_nModes == 0)
        {
            m_gramPInvs.resize(m_modeCounts.size());
            return 0;
        }

        m_basis.resize(static_cast<int>(m_coords.size()), m_nModes);
        for(int n = 0; n < m_nModes; ++n)
        {
            for(size_t p = 0; p < m_coords.size(); ++p)
            {
                m_basis(static_cast<int>(p), n) = modes.image(n)(m_coords[p].first, m_coords[p].second);
            }
        }

        m_gramPInvs.reserve(m_modeCounts.size());
        for(int nModesFit : m_modeCounts)
        {
            if(nModesFit < 0 || nModesFit > m_nModes)
            {
                std::cerr << "pseudo-inverse fit requested " << nModesFit
                          << " modes but only " << m_nModes << " are available\n";
                return -1;
            }

            m_gramPInvs.emplace_back();
            if(nModesFit == 0)
            {
                continue;
            }

            Eigen::Array<double, -1, -1> gram =
                ((m_basis.leftCols(nModesFit).transpose() * m_basis.leftCols(nModesFit)) /
                 static_cast<double>(m_coords.size())).array();
            Eigen::Array<double, -1, -1> U;
            Eigen::Array<double, -1, -1> S;
            Eigen::Array<double, -1, -1> VT;
            double condition = 0;
            int nRejected = 0;
            double effectiveMaxCondition = maxCondition;

            int rv = mx::math::eigenPseudoInverse(m_gramPInvs.back(),
                                                  condition,
                                                  nRejected,
                                                  U,
                                                  S,
                                                  VT,
                                                  gram,
                                                  effectiveMaxCondition,
                                                  alpha);
            if(rv != 0)
            {
                std::cerr << "pseudo-inverse basis fit failed with rv = " << rv
                          << " for " << nModesFit << " modes\n";
                return -1;
            }

            std::cerr << "pinv Gram modes: " << nModesFit
                      << " condition: " << condition
                      << " rejected: " << nRejected
                      << " requestedMaxCondition: " << maxCondition
                      << " effectiveMaxCondition: " << effectiveMaxCondition
                      << " alpha: " << alpha << '\n';
        }

        return 0;
    }

    int subtract(std::vector<realT> & amplitudes,
                 mx::improc::eigenImage<realT> & image,
                 size_t cutoffIndex) const
    {
        int nModesFit = m_modeCounts[cutoffIndex];
        amplitudes.assign(nModesFit, 0);

        Eigen::VectorXd pupilValues(static_cast<int>(m_coords.size()));
        double mean = 0;
        for(size_t p = 0; p < m_coords.size(); ++p)
        {
            pupilValues(static_cast<int>(p)) = image(m_coords[p].first, m_coords[p].second);
            mean += pupilValues(static_cast<int>(p));
        }

        mean /= static_cast<double>(m_coords.size());
        pupilValues.array() -= mean;

        if(nModesFit > 0)
        {
            Eigen::VectorXd rhs =
                (m_basis.leftCols(nModesFit).transpose() * pupilValues) /
                static_cast<double>(m_coords.size());
            Eigen::VectorXd coefficients = m_gramPInvs[cutoffIndex].matrix() * rhs;

            pupilValues.noalias() -= m_basis.leftCols(nModesFit) * coefficients;
            for(int n = 0; n < nModesFit; ++n)
            {
                amplitudes[n] = static_cast<realT>(coefficients(n));
            }
        }

        for(size_t p = 0; p < m_coords.size(); ++p)
        {
            image(m_coords[p].first, m_coords[p].second) =
                static_cast<realT>(pupilValues(static_cast<int>(p)));
        }

        return 0;
    }

  private:
    int m_nModes {0};
    std::vector<int> m_modeCounts;
    std::vector<std::pair<int, int>> m_coords;
    Eigen::MatrixXd m_basis;
    std::vector<Eigen::Array<double, -1, -1>> m_gramPInvs;
};

} // namespace aperture_stroke

#endif
