#ifndef apertureStroke_basisFitters_hpp
#define apertureStroke_basisFitters_hpp

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
