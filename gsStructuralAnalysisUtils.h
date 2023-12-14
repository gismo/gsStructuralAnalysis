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
class gsStructuralAnalysisOutput
{
public:
    gsStructuralAnalysisOutput(const std::string name)
    :
    m_fname(name),
    m_precision(5),
    m_nPointHeaders(0),
    m_nOtherHeaders(0)
    {

    }

    gsStructuralAnalysisOutput(const std::string name, const gsMatrix<T> & points)
    :
    gsStructuralAnalysisOutput(name)
    {
        m_points = points;
    }

    void clean()
    {
        std::ofstream file;
        file.open(m_fname, std::ofstream::out | std::ofstream::trunc);
        file.close();
    }

    void _initPointHeader(const std::string name, std::vector<std::string> pointHeaders)
    {
        std::ofstream file;
        file.open(name,std::ofstream::out | std::ofstream::app);
        for (index_t p=0; p!=m_points.cols(); p++)
            for (size_t h=0; h!=pointHeaders.size(); h++)
                file<< "point"<<p<<"-"<<pointHeaders[h]<< ",";
        file.close();
    }

    void _initOtherHeader(const std::string name, std::vector<std::string> otherHeaders)
    {
        std::ofstream file;
        file.open(name,std::ofstream::out | std::ofstream::app);
        for (size_t h=0; h!=otherHeaders.size(); h++)
            file<<otherHeaders[h]<< ",";
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

    void _writeOtherData(const std::string name, const gsVector<T> & otherData)
    {
        std::ofstream file;
        file.open(name,std::ofstream::out | std::ofstream::app);
        for (index_t c=0; c!=otherData.size(); c++)
            file  << otherData.at(c) << ",";
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

    void init(std::vector<std::string> pointHeaders, std::vector<std::string> otherHeaders)
    {
        m_nPointHeaders = pointHeaders.size();
        m_nOtherHeaders = otherHeaders.size();
        std::ofstream file;
        file.open(m_fname,std::ofstream::out);
        file  << std::setprecision(m_precision);
        file.close();

        this->_initPointHeader(m_fname, pointHeaders);
        this->_initOtherHeader(m_fname, otherHeaders);

        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << "\n";
        file.close();
    }

    void init(std::vector<std::string> pointHeaders)
    {
        std::vector<std::string> empty;
        this->init(pointHeaders,empty);
    }

    void add(const gsMatrix<T> & pointSolutions, const gsVector<T> & otherData)
    {
        GISMO_ASSERT(pointSolutions.rows()==m_nPointHeaders,"Number of solutions per point is different from the defined number of headers. "<<pointSolutions.rows()<<" = pointSolutions.rows()==m_nPointHeaders = "<<m_nPointHeaders);
        GISMO_ASSERT(otherData.size()==m_nOtherHeaders,"Number of other data is different from the defined number of headers");
        std::ofstream file;
        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << std::setprecision(m_precision);
        file.close();

        this->_writePointData(m_fname,pointSolutions);
        this->_writeOtherData(m_fname,otherData);

        file.open(m_fname,std::ofstream::out | std::ofstream::app);
        file  << "\n";
        file.close();
    }


    void add(const gsMatrix<T> & pointSolutions)
    {
        gsVector<T> empty;
        this->add(pointSolutions,empty);
    }

    void setPrecision(index_t precision) { m_precision = precision; }

protected:
    gsMatrix<T> m_points;
    std::string m_fname;
    std::ofstream m_file;
    index_t m_precision;
    index_t m_nPointHeaders, m_nOtherHeaders;

};

} // namespace gismo
