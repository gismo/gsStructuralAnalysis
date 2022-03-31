 /** @file gsStructuralAnalysisUtils.h

    @brief Class providing Structural Analysis Utilities

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once


namespace gismo
{

template<class T>
class gsStructuralAnalysisOutputBase
{
public:
    gsStructuralAnalysisOutputBase(const std::string name)
    :
    m_fname(name),
    m_precision(5)
    {

    }

    gsStructuralAnalysisOutputBase(const std::string name, const gsMatrix<T> & points)
    :
    gsStructuralAnalysisOutputBase(name)
    {
        m_points = points;
    }

    void clean()
    {
        std::ofstream file;
        file.open(m_fname, std::ofstream::out | std::ofstream::trunc);
        file.close();
    }

    void _initPointHeader(const std::string name, std::vector<std::string> headers)
    {
        std::ofstream file;
        file.open(name,std::ofstream::out | std::ofstream::app);
        for (index_t p=0; p!=m_points.cols(); p++)
            for (index_t h=0; h!=headers.size(); h++)
                file<< "point"<<p<<"-"<<headers[h]<< ",";
        file.close();
    }

    void _writePointData(const std::string name, const gsMatrix<T> & pointData)
    {
        std::ofstream file;
        file.open(name,std::ofstream::out | std::ofstream::app);
        for (index_t p=0; p!=pointData.cols(); p++)
            for (index_t c=0; c!=pointData.rows(); c++)
                file  << pointData(c,p) << ",";
        file.close();
    }

    void writePointCoords(const std::string name)
    {
        std::ofstream file;
        file.open(name,std::ofstream::out);
        file  << std::setprecision(m_precision);
        _initPointHeader(file);

        for (index_t j=0; j!=m_points.cols(); j++)
            for (index_t i=0; i!=m_points.rows(); i++)
                file<<m_points(i,j)<<",";
        file.close();
    }

    void setPrecision(index_t precision) { m_precision = precision; }

protected:
    gsMatrix<T> m_points;
    std::string m_fname;
    std::ofstream m_file;
    index_t m_precision;
};

template<class T>
class gsALMOutput : public gsStructuralAnalysisOutputBase<T>
{
    typedef gsStructuralAnalysisOutputBase<T> Base;

public:
    gsALMOutput(const std::string name)
    :
    Base(name)
    {

    }

    gsALMOutput(const std::string name, const gsMatrix<T> & points)
    :
    Base(name,points)
    {

    }

    void init(std::vector<std::string> headers)
    {
        m_nheaders = headers.size();
        std::ofstream file;
        file.open(m_fname,std::ofstream::out);
        file  << std::setprecision(m_precision)
              << "Deformation norm" << ",";
        file.close();

        Base::_initPointHeader(m_fname, headers);

        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << "Lambda" << ","
              << "Indicator"
              << "\n";
        file.close();
    }

    void add(const gsMatrix<T> & pointSolutions, const gsVector<T> & solutionVector, const T lambda, const T indicator)
    {
        GISMO_ASSERT(pointSolutions.rows()==m_nheaders,"Number of solutions per point is different from the defined number of headers");
        std::ofstream file;
        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << std::setprecision(m_precision)
                << solutionVector.norm() << ",";
        file.close();

        Base::_writePointData(m_fname,pointSolutions);

        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << lambda << ","
                << indicator << ","
                << "\n";
        file.close();
    }

    void add(const gsMatrix<T> & pointSolutions, const T lambda)
    {
        GISMO_ASSERT(pointSolutions.rows()==m_nheaders,"Number of solutions per point is different from the defined number of headers");
        std::ofstream file;
        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << std::setprecision(m_precision)
                << "NA" << ",";
        file.close();

        Base::_writePointData(m_fname,pointSolutions);

        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << lambda << ","
                << "NA" << ","
                << "\n";
        file.close();
    }

protected:
    using Base::m_points;
    using Base::m_fname;
    using Base::m_precision;
    index_t m_nheaders;

};

template<class T>
class gsStaticOutput : public gsStructuralAnalysisOutputBase<T>
{
    typedef gsStructuralAnalysisOutputBase<T> Base;

public:
    gsStaticOutput(const std::string name)
    :
    Base(name)
    {

    }

    gsStaticOutput(const std::string name, const gsMatrix<T> & points)
    :
    Base(name,points)
    {

    }

    void init(std::vector<std::string> headers)
    {
        m_nheaders = headers.size();
        std::ofstream file;
        file.open(m_fname,std::ofstream::out);
        file  << std::setprecision(m_precision);
        file.close();

        Base::_initPointHeader(m_fname, headers);

        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file<< "\n";
        file.close();
    }

    void add(const gsMatrix<T> & pointSolutions)
    {
        GISMO_ASSERT(pointSolutions.rows()==m_nheaders,"Number of solutions per point is different from the defined number of headers");
        std::ofstream file;
        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << std::setprecision(m_precision);
        file.close();

        Base::_writePointData(m_fname,pointSolutions);

        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << "\n";
        file.close();
    }

protected:
    using Base::m_points;
    using Base::m_fname;
    using Base::m_precision;
    index_t m_nheaders;

};

} // namespace gismo
