#ifndef _LINEPARSER_
#define _LINEPARSER_

/**
 * @file LineParser.h
 * @brief Declaration of the LineParser class
 */

#include <vector>
#include <string>
using namespace std;

/// A device for breaking a line of text into words (white space
/// delimited only)
/**
 * A LineParser is a device for breaking a line of text (in a string)
 * into individual words. The assumption is that words are white-space
 * delimited.
 */

class LineParser {
private:
  /// This holds the line of text
  string theLine;
  /// This holds the individual words
  vector<string> theWords;
  /// Does the actual grunt work of parsing into words.
  void parse();
public:
  /// Default constructor
  /**
   * Creates an empty LineParser.
   */
  LineParser() {
    clear();
  }

  /// Primary constructor
  /**
   * Takes a line of text and breaks it into individual words.
   * @param text the line of text to parse
   */
  LineParser(const string& text) {
    clear();
    theLine = text;
    parse();
  }

  /// Copy constructor
  LineParser(const LineParser& that) {
    clear();
    theLine = that.theLine;
    parse();
  }

  /// Erase the data in this LineParser
  /**
   * There is little reason to use this method.
   */
  void clear() {
    theLine.clear();
    theWords.clear();
  }

  /// How many words are held in this LineParser?
  /**
   * @return the number of words
   */
  long size() const {
    return theWords.size();
  }

  /// Get a given word (starting with word zero)
  /**
   * @param k index of the word we want (between 0 and <code>size()-1</code>).
   * @return the k'th word
   */
  string operator[](int k) const;

  /// What is the line upon which this LineParser was built?
  /**
   * @return the line of text we parsed
   */
  string get_line() const { return theLine; }
  
};




#endif
