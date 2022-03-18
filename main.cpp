/* 
 * File:   main.cpp
 * Author: mwittig
 *
 * Created on December 22, 2010, 9:56 AM
 */

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <iterator>
#include <time.h>


#include "CTwoBit.h"
#include "ParsedTextfile.h"

#define FLANK_LENGTH 60
#define BASE 3

using namespace std;

vector<string> GetDirectoryListing(const char* path);
vector<string> GetDirectoryListing(const char* ,const char* , bool);
vector<string> strsplit(const string& str, const string& delimiters);
vector<string> strsplit(const string& str);
char getBotTop(char five, char three);
char GetComplNuc(const char& source);
string GetComplSequence(const string strSequence, bool ucase =true);
vector<string> GetParsedLine(string str,const string& delimiters );
void FindApproxSubs_Schwer(unsigned int k, string & P, string & T, vector<unsigned int> & Out);
void FindApproxSubs_Myers(unsigned int k, string & P, string & T, vector<unsigned int> & Out);
unsigned int AlignmentScore(string ,string,unsigned int k = 2);
char getBotTop(const char* chrom,uint32_t pos, CTwoBit& twobit);
char getPlusMinus(CParsedTextfile& manifest, CTwoBit& twobit);
void invalidSnpReport(string nuc,string snp,char strand);
pair<string,string>     getAlleles(const string&);

/*
 * 
 */

//
// /home/mwittig/coding/cpp/BotTop2PlusMinus/dist/Release/GNU-Linux-x86/bottop2plusminus
//
int main(int argc, char** argv) {

    try
    {
        if(argc==1)
        {
            cout << "run like: " << argv[0] << ' ' << "twoBitFile.2bit" << ' ' << "Manifest.csv" << endl;
            cout << "          output goes to stdout, errors to stderr" << endl;
            return EXIT_FAILURE;
        }
        if(argc != 3)
        {
            return  system(argv[0]);
        }
        CTwoBit twoBit(argv[1]);


        // read SNP annotation
        //
        //
        //const char* SNP_Annotation = "/home/mwittig/mnts/schleuse on flacker/Immuno_BeadChip_11419691_B.csv";
        //const char* SNP_Annotation = "/data/tmp/tst_Immuno_BeadChip_11419691_B.csv";
        CParsedTextfile manifest(argv[2],",");
        if(!manifest.isOpen())
            throw(string("Could not open ")+argv[2]);
        if(!manifest.ProceedTo("[Assay]"))
            throw(string("Couldn't find data start"));
        if(!manifest.GetHeader())
           throw(string("Couldn't read header")); 

        unsigned int seekg_column = 0;
        cout << "probe_set_id\tprobe_set_id\tChr\tMapInfo\tStrand\tallele_A\tallele_B\tflank\tseekg_column\tIlmnStrand\tSNP\tstrandSeq\tfwdSeq\tsql" << endl;
        while(manifest.Next())
        {
            char strand = getPlusMinus(manifest, twoBit);
            

            string chrom = string("chr")+manifest["Chr"];
            if(chrom.compare("chrMT") == 0)
                chrom="chrM";
            if(chrom.compare("chrXY") == 0)
                chrom="chrX";
            
            uint32_t pos = atol(manifest["MapInfo"].c_str());

            string fivePrime  = twoBit.getSequence(chrom.c_str(),pos-16,pos-1);
            string threePrime = twoBit.getSequence(chrom.c_str(),pos+1,pos+16);
            string origSeq    = fivePrime+'['+twoBit.getSequence(chrom.c_str(),pos,pos)+']'+threePrime;
            char   plus_strand_bot_top = getBotTop(chrom.c_str(),pos,twoBit);
            pair<string,string> alleles = getAlleles(manifest["SNP"]);
            
            
            if(strand == '-')
            {
                string strHlp = fivePrime;
                fivePrime  = GetComplSequence(threePrime,true);
                threePrime = GetComplSequence(strHlp,true);
            }
            
            //invalidSnpReport(origSeq,manifest["SNP"],strand);
            
            ostringstream sql("");
            sql << "update annotation set strand = \"" << strand << "\", flank = \"" << fivePrime+manifest["SNP"]+threePrime << "\" where probe_set_id = \"" << manifest["Name"] << "\";";
            cout << manifest["Name"] << '\t' 
                 << manifest["Name"] << '\t' 
                 << manifest["Chr"] << '\t' 
                 << manifest["MapInfo"] << '\t'
                 << strand << '\t'
                 << alleles.first << '\t'
                 << alleles.second << '\t'
                 << fivePrime+manifest["SNP"]+threePrime << '\t'
                 << seekg_column++ << '\t'
                 << manifest["IlmnStrand"] << '\t'
                 << manifest["SNP"] << '\t'
                 << origSeq << '\t'
                 << sql.str() << endl;
        }
       
        return (EXIT_SUCCESS);
    }
    catch(const string& err)
    {
        cerr << "Error:" << err << endl;
    }
    catch(...)
    {
        cerr << "unknown error" << endl;
    }
    return(EXIT_FAILURE);
}

char getPlusMinus(CParsedTextfile& manifest, CTwoBit& twobit)
{
    if(manifest["IlmnStrand"].compare("PLUS")==0)
        return '+';
    if(manifest["IlmnStrand"].compare("MINUS")==0)
        return '-';

    string act_chrom = manifest["Chr"];
    if(act_chrom.compare("XY")==0)
    {
        cerr << "unknown chromosome XY treated as chromosmome X" << endl;
        act_chrom="X";
        //cerr << "unknown chromosome XY ignored" << endl;
    }
    uint32_t position = atol(manifest["MapInfo"].c_str());

    char Allele = toupper(twobit.getSequence((string("chr")+act_chrom).c_str(),position,position)[0]);

    if(
           ( (manifest["SNP"].compare("[A/T]") == 0 || manifest["SNP"].compare("[T/A]") == 0) && (Allele == 'C' || Allele == 'G') ) ||
           ( (manifest["SNP"].compare("[C/G]") == 0 || manifest["SNP"].compare("[G/C]") == 0) && (Allele == 'T' || Allele == 'A') )
       )
    {
        cerr << manifest["Name"] << " with " << manifest["SNP"] <<  " at " << manifest["Chr"] << ':' << position << " shows inconsistency as genomic base is " << Allele << endl;
        return 'u';
    }
    if(manifest["SNP"].compare("[A/C]") == 0 || manifest["SNP"].compare("[C/A]") == 0)
    {
        switch (Allele)
        {
            case 'A':
            case 'C':
                return '+';
            case 'T':
            case 'G':
                return '-';
            default:
                return 'u';
        }
    }
    if(manifest["SNP"].compare("[A/G]") == 0 || manifest["SNP"].compare("[G/A]") == 0)
    {
        switch (Allele)
        {
            case 'A':
            case 'G':
                return '+';
            case 'T':
            case 'C':
                return '-';
            default:
                return 'u';
        }
    }
    if(manifest["SNP"].compare("[T/C]") == 0 || manifest["SNP"].compare("[C/T]") == 0)
    {
        switch (Allele)
        {
            case 'T':
            case 'C':
                return '+';
            case 'A':
            case 'G':
                return '-';
            default:
                return 'u';
        }
    }
    if(manifest["SNP"].compare("[T/G]") == 0 || manifest["SNP"].compare("[G/T]") == 0)
    {
        switch (Allele)
        {
            case 'T':
            case 'G':
                return '+';
            case 'A':
            case 'C':
                return '-';
            default:
                return 'u';
        }
    }
    char BOTTOP = getBotTop((string("chr")+manifest["Chr"]).c_str(),atol(manifest["MapInfo"].c_str()),twobit);
    if(BOTTOP == 'u')
        return 'u';
    if(BOTTOP == 'b' && manifest["IlmnStrand"].compare("BOT")==0)
        return '+';
    if(BOTTOP == 'b' && manifest["IlmnStrand"].compare("TOP")==0)
        return '-';
    if(BOTTOP == 't' && manifest["IlmnStrand"].compare("BOT")==0)
        return '-';
    if(BOTTOP == 't' && manifest["IlmnStrand"].compare("TOP")==0)
        return '+';
    return 'u';

}



char getBotTop(const char* chrom,uint32_t pos, CTwoBit& twobit)
{
    uint32_t    Pos5 = pos-1;
    uint32_t    Pos3 = pos+1;
    char cReturn = 'u';

    
    char ref_allele = twobit.getSequence(chrom,pos,pos)[0];
    switch(ref_allele)
    {
        case 'A':
            return 't';
        case 'T':
            return 'b';
    }
    
    while(Pos5!=0 && Pos3 < twobit.getChromSize(chrom) && Pos3-Pos5 < 100  && cReturn == 'u')
    {
        string seq5 = twobit.getSequence(chrom,Pos5,Pos5);
        string seq3 = twobit.getSequence(chrom,Pos3,Pos3);
        if(seq5.size() == 0 || seq3.size() == 0)
            return 'u';
        cReturn = getBotTop(seq5[0],seq3[0]);
        Pos3++;Pos5--;
    }
    return cReturn;
}

char getBotTop(char five, char three)
{
    five = toupper(five);
    three= toupper(three);

    char return_value = 'u';
    switch(five)
    {
        case 'A':
        case 'T':
            switch(three)
            {
                case 'C':
                case 'G':
                    return_value = 't';
                    break;
                default:
                    break;
            }
            break;
        case 'C':
        case 'G':
            switch(three)
            {
                case 'A':
                case 'T':
                    return_value = 'b';
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    return return_value;
}

vector<string> GetDirectoryListing(const char* path)
{
    DIR *dirp;
    struct dirent *entry;
    vector<string> strReturn;

    if(dirp = opendir(path))
    {
        while(entry = readdir(dirp))
            strReturn.push_back(entry->d_name);
        closedir(dirp);
    }
    return strReturn;
}

vector<string> GetDirectoryListing(const char* path,const char* contains, bool fullpath = true)
{
    vector<string> strReturn;
    vector<string> strhelp = GetDirectoryListing(path);
    for(vector<string>::iterator iter = strhelp.begin();iter != strhelp.end();iter++)
        if(iter->find(contains)!=string::npos)
            strReturn.push_back(string(path)+'/'+(*iter));
        
    return strReturn;
}

vector<string> strsplit(const string& val)
{
    return strsplit(val,"\t");
}

vector<string> strsplit(const string& str, const string& delimiters = "\t")
{
    vector<string> tokens;
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
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

string GetComplSequence(const string strSequence, bool ucase)
{
	try
    {
        int intPos = 0;
    	size_t intSize = strSequence.length();
    	char* buffer=(char*)malloc(intSize+1);
        if(!buffer)
            throw(string("can't allocate enough memory."));
        buffer[intSize]=0;
        while(--intSize != -1)
        {
            if(ucase)
                buffer[intPos++]=toupper(GetComplNuc(strSequence.at(intSize)));
            else
                buffer[intPos++]=GetComplNuc(strSequence.at(intSize));
        }
       	string strReturn(buffer);
        free(buffer);
        return strReturn;
    }
    catch(string strError)
    {
        throw(string("error in CMyTools::GetSeqFromFile, ")+strError);
    }
    catch(...)
    {
        throw("unknown error in CMyTools::GetSeqFromFile\n");
    }

}


char GetComplNuc(const char& source)
{
	switch(source)
	{
	case 'A':
		return 'T';
	case 'a':
		return 't';
	case 'C':
		return 'G';
	case 'c':
		return 'g';
	case 'G':
		return 'C';
	case 'g':
		return 'c';
	case 'T':
	case 'U':
		return 'A';
	case 't':
	case 'u':
		return 'a';
	case 'M':
		return 'K';
	case 'm':
		return 'k';
	case 'R':
		return 'Y';
	case 'r':
		return 'y';
	case 'Y':
		return 'R';
	case 'y':
		return 'r';
	case 'K':
		return 'M';
	case 'k':
		return 'm';
	case 'V':
		return 'B';
	case 'v':
		return 'b';
	case 'H':
		return 'D';
	case 'h':
		return 'd';
	case 'D':
		return 'H';
	case 'd':
		return 'h';
	case 'B':
		return 'V';
	case 'b':
		return 'v';
	case '(':
         return ')';
    case ')':
         return '(';
	case '[':
         return ']';
    case ']':
         return '[';
	}
	return source;
}

vector<string> GetParsedLine(string str,const string& delimiters = "\t")
{
    vector<string> tokens;
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
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


void FindApproxSubs_Schwer(unsigned int k, string & P, string & T, vector<unsigned int> & Out)
{
  // Aufgabe: Finde alle Positionen in T, an denen ein Substring von T mit P
  // annähernd übereinstimmt (mit Abstand <= k), und gebe diese in Out aus.

  unsigned int n = (unsigned int )T.length();
  unsigned int m = (unsigned int )P.length();
  unsigned int *C = new unsigned int[m];
  unsigned int last = k + 1;
  unsigned int pC,nC, i, pos;


  for(i = 0; i < m; i++) { C[i] = i + 1; }

  // Searching
  for(pos = 0; pos < n; pos++) {
    pC = 0;
    nC = 0;
    for(i = 0; i < last; i++) {
      if(P[i] == T[pos]) {
	nC = pC;
      } else {
	if(pC < nC) { nC = pC; }
	if(C[i] < nC) { nC = C[i]; }
	nC++;
      }
      pC = C[i];
      C[i] = nC;
    }
    while(C[last - 1] > k) { last--; }
    if(last == m) {
      Out.push_back(pos + 1);
    } else { last++; }
  }

}

void FindApproxSubs_Myers(unsigned int k, string & P, string & T, vector<unsigned int> & Out)
{
  // Aufgabe: Finde alle Positionen in T, an denen ein Substring von T mit P
  // annähernd übereinstimmt (mit Abstand <= k), und gebe diese in Out aus.

  unsigned int n = (unsigned int )T.length();
  unsigned int m = (unsigned int )P.length();
  unsigned int score = m;
  unsigned int i, j, pos, HN, HP, X, D0;
  unsigned int stop;

  unsigned int max_index = int((m-1) / BASE) + 1;
  unsigned int** B = new unsigned int*[256];
  for(i = 0; i < 256; i++) {
    B[i] = new unsigned int[max_index];
  }
  unsigned int *VP = new unsigned int [max_index];
  unsigned int *VN = new unsigned int [max_index];

  for(j = 0; j < m; j += BASE) {
    unsigned int index = (int) j / BASE;
    VP[index] = 0xFFFFFFFF;
    VN[index] = 0;
    for(i = 0; i < 256; i++) { B[i][index] = 0; }
    unsigned int eins = 1;
    for(i = j; i < m && i < j + BASE; i++) {
      B[P[i]][index] |= eins;
      eins<<=1;
    }
  }

  unsigned int VN_UEBERTRAG, VP_UEBERTRAG, temp;
  // scanning the text
  for(pos = 0; pos < n; pos++) {

    VN_UEBERTRAG = 0;
    VP_UEBERTRAG = 0;
    for(j = 0; j < max_index ; j++) {
      if(VN_UEBERTRAG > 0) { VN[j] |= 0x1; }
      if(VP_UEBERTRAG > 0) { VP[j] |= 0x1; }
      X = B[T[pos]][j] | VN[j];
      D0 = ((VP[j] + (X & VP[j])) ^ VP[j]) | X;
      HN = VP[j] & D0;
      HP = VN[j] |~ (VP[j] | D0);
      X = HP << 1;
      VN_UEBERTRAG = HP & 1<<BASE-1;
      VN[j] = X & D0;
      temp = HN << 1;
      VP_UEBERTRAG = temp & 1<<BASE-1 |~ VN_UEBERTRAG;
      VP[j] = temp |~ (X | D0);
    }

    stop = 1<<m-(max_index-1)*BASE-1;
    if(HP & stop) { score++; }
    else { if(HN & stop) { score--; } }
    if(score <= k) { Out.push_back(pos + 1); }
  }
}

unsigned int AlignmentScore(string P,string T, unsigned int k)
{
    vector<unsigned int> Out;

    FindApproxSubs_Schwer(k, P, T, Out);
    if(Out.size()==0)
        return 0;
    return Out.at(0);
}

void invalidSnpReport(string nuc,string snp,char strand)
{
    if(strand == '+' && snp.find(nuc) == string::npos)
        cerr << snp << '\t' << nuc << '\t' << strand << "\tdoes not fit ..." << endl;
    else if(strand == '-' && GetComplSequence(snp).find(nuc) == string::npos)
        cerr << snp << '\t' << nuc << '\t' << strand << "\tdoes not fit ..." << endl;
}

pair<string,string> getAlleles(const string& theSNP)
{
    string snp = theSNP.substr(1,theSNP.size()-2);
    vector<string>      parsed = strsplit(snp,"/");
    if(parsed.size() != 2)
        throw(string("Invalid SNP found \"")+theSNP+'\"');
    return pair<string,string>(parsed[0],parsed[1]);
}



