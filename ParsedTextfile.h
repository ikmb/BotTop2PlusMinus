#ifndef ParsedTextfile_H_
#define ParsedTextfile_H_

using namespace std;

typedef vector<string> ParsedTextfile_OneLine;

class CParsedTextfile
{
public:
	CParsedTextfile(); 
	CParsedTextfile(const CParsedTextfile&); 
	CParsedTextfile(const char*,const char* sep = "\t");
	
	virtual ~CParsedTextfile();

        string operator[](unsigned int idx);
	string operator[](const char*);
	string operator[](string& v){return (*this)[v.c_str()];}
	CParsedTextfile& operator=(const CParsedTextfile&);
	
	bool Next();
        bool ProceedTo(const char*);
        bool isOpen(){return m_inputStream->is_open();}

        bool GetHeader();


private:
	
	map<string,int>  	mapHeader;
	string                  m_strValidSeparator;
        ifstream*               m_inputStream;
        vector<string>          m_actLine;


        vector<string> strsplit(const string& str, const string& delimiters);

	
	bool Open(const char*);
        void Close(){if(m_inputStream->is_open()){m_inputStream->close();delete m_inputStream; m_inputStream=NULL;}}
	
        void initVariables();

};

#endif /*ParsedTextfile_H_*/


