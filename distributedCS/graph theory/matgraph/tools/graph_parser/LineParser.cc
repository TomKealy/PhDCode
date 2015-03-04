#include "LineParser.hxx"

/**
 * @file LineParser.cc
 * @brief Supporting code for the LineParser class
 */

void LineParser::parse() {
  theWords.clear();
  const string white(" \t\n");  // white space characters space, tab, newline
  
  long beg, end;

  beg = theLine.find_first_not_of(white); // start of first word


  while (beg != (long) string::npos) {
    end = theLine.find_first_of(white,beg); // end of current word

    if (end == (long) string::npos) {
      end = theLine.length();
    }

    // Add word to vector
    theWords.push_back(theLine.substr(beg,end-beg));

    // Find start of next word
    beg = theLine.find_first_not_of(white, end);
  }
}

string LineParser::operator[](int k) const {
  return theWords[k];
}
