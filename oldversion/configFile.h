
#ifndef GUARD_CONFIGFILE_H
#define GUARD_CONFIGFILE_H

#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

class Config {
public:
	Config( const char * fname,
			const char separator = '=',
			const char comment = '#');


	template<class T>
	T get_parameter(const std::string& parameterName) const
	{
		// check if parameter name is in dict
		T val;
		std::stringstream ss(dict.find(parameterName)->second);
		ss >> val;
		
		return val;	

	}
private:
	std::map< std::string, std::string > dict;


};

namespace configFunc {
	void trim(std::string& str);
};



#endif // GUARD_CONFIGFILE_H

