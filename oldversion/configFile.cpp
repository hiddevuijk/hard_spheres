

#include "configFile.h"

using namespace std;


Config::Config( const char * fname,
			const char separator_char ,
			const char comment_char )
{

	// read file
	ifstream file(fname);

	// check if file exists
	if ( !file.is_open() )
		throw "fuck this ";

	string line, lhs, rhs;
	size_t cmnt_index, sep_index;

	//for(; getline(file,line); ) {
	while( getline(file,line) ) {
		

		// check for comment, if # is inline, remove rhs
		cmnt_index = line.find(comment_char);
		if(cmnt_index != string::npos) {
			line.erase(cmnt_index,string::npos);
		}
		
		if( line.size() == 0) continue;

		// fine the separator char
		sep_index = line.find(separator_char);
	
		// if not found skip this line
		if( sep_index != string::npos ) {
			// check if valid	
			lhs = line.substr(0,sep_index);
			rhs = line.substr(sep_index+1,line.size()-1);
			//remove leading and trailing white spaces
			configFunc::trim(lhs);
			configFunc::trim(rhs);		
			// check if not empty!!!
			// check for double keys!!!
			dict[lhs] = rhs;	
		}
	}


}

void configFunc::trim(string& str)
{

	unsigned int N = str.size();
	// check if N > 0

	// remove leading white space 
	while(true) {
		if( str[0] == ' ') {
			str.erase(0,1);
		} else {
			break;
		}
	}

	// remove trailing white space
	for( unsigned int i=N-1;i>0;--i) {
		if(str[i] == '\0') continue;
		if(str[i] == ' ' ) {
			str.erase(str.end()-1);
		} else {
			break;
		}
	}
}

