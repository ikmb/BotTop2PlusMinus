#include <string>
#include <map>
#include <vector>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <iostream>

#include "ParsedTextfile.h"


using namespace std;

void CParsedTextfile::initVariables()
{
    m_strValidSeparator="\t";
    m_inputStream = NULL;
}

CParsedTextfile::CParsedTextfile()
{
    initVariables();
}

CParsedTextfile::CParsedTextfile(const CParsedTextfile& orig)
{
    initVariables();
    *this=orig;
}

CParsedTextfile::CParsedTextfile(const char* filename,const char* sep)
{
    initVariables();
    m_strValidSeparator=sep;
    if(!Open(filename))
        throw(string("Error opening ")+filename);
}

CParsedTextfile::~CParsedTextfile()
{
    Close();
}

CParsedTextfile& CParsedTextfile::operator=(const CParsedTextfile& orig)
{
    mapHeader=orig.mapHeader;
    m_strValidSeparator=orig.m_strValidSeparator;
    m_inputStream = orig.m_inputStream;
    return *this;
}


bool CParsedTextfile::GetHeader()
{
    try
    {
        string strLine;
        if(!getline(*m_inputStream,strLine))
            return false;

        m_actLine = strsplit(strLine,m_strValidSeparator);
        int i=0;
        mapHeader.clear();

        for(vector<string>::iterator iter=m_actLine.begin();iter != m_actLine.end();iter++, i++)
            mapHeader[*iter]=i;
        return true;
    }
    catch(...)
    {
            return false;	
    }
}

bool CParsedTextfile::Open(const char* input_file)
{
    try
    {
        if(m_inputStream == NULL)
            m_inputStream = new ifstream();
        if(m_inputStream->is_open())
        {
            cerr << "input file already open. Close before opening a new one." << endl;
            return false;
        }
        m_inputStream->open(input_file);
        if(m_inputStream->is_open())
            return true;
        return false;
    }
    catch(...)
    {
            return false;
    }
}

bool CParsedTextfile::Next()
{
    string strLine;
    if(getline(*m_inputStream,strLine))
    {
        m_actLine = strsplit(strLine,m_strValidSeparator);
        /*if(m_actLine.size() != mapHeader.size())
        {
            cerr << "warning: entry has " << m_actLine.size() << " items, header has " << mapHeader.size() << endl;
        }*/
        return true;
    }
    return false;
}

bool CParsedTextfile::ProceedTo(const char* entry)
{
    string strLine;
    while(getline(*m_inputStream,strLine))
    {
        if(strLine.compare(entry)==0)
            return true;
    }
    return false;
}

string CParsedTextfile::operator[](const char* strCol)
{
    if(mapHeader.find(strCol)==mapHeader.end())
            throw string("Invalid header column name \"")+strCol+'\"';
    return m_actLine[mapHeader[string(strCol)]];
}

string CParsedTextfile::operator[](unsigned int idx)
{
    if(idx >= m_actLine.size())
            throw "index out of bounds";
    return m_actLine[idx];
}

vector<string> CParsedTextfile::strsplit(const string& str, const string& delimiters = "\t")
{
    vector<string> tokens;
    // Skip delimiters at beginning.
    string::size_type lastPos = 0;
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        if(str.length() > pos && delimiters.find(str[pos]) != string::npos )
            lastPos=pos+1;
        else
            lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
    return tokens;
}


