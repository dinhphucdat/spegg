#ifndef INPUT_FILE_CHECKER_H
#define INPUT_FILE_CHECKER_H

#include <algorithm>
#include <istream>
#include <fstream>

/**
 * @brief Counts the number of lines in a file
 * 
 * @param filename name of the file
 * @return int number of lines in the file
 */
int count_lines_in_file(const char *filename)
	{
	int answer;
	std::ifstream inFile(filename);
	answer = std::count(std::istreambuf_iterator<char>(inFile),
		   	    std::istreambuf_iterator<char>(),'\n');
	return answer;
	}

#endif
