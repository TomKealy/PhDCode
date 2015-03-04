#include <iostream>
#include <fstream>
#include <map>
#include "LineParser.hxx"
using namespace std; 

int main() {
  const char* OUTFILE = "parsed_graph.m"; 
  ofstream fout(OUTFILE);
  if (!fout) {
    cerr << "Unable to open file " << OUTFILE << " for output" << endl;
    exit(1);
  }
  
  string theLine;           // line of text read from cin
  long nv = 0;              // number of vertices read
  map<string, long> dict;   // name to vertex number table
  map<long, string> undict; // vertex number to name table

  fout << "function parsed_graph(g)" << endl;
  fout << "edges = [" << endl;
  while(getline(cin,theLine)) {
    LineParser LP(theLine); 

    // if no words on this line, skip it
    if (LP.size() == 0) continue;

    // Check if first character of first word is #
    if ( *(LP[0].c_str()) == '#') continue;

    if (dict.count(LP[0]) == 0) {
      ++nv;
      dict[LP[0]] = nv;
      undict[nv] = LP[0];
    }
    
    if (LP.size() > 1) {
      if (dict.count(LP[1]) == 0) {
	++nv;
	dict[LP[1]] = nv;
	undict[nv] = LP[1];
      }
      long u,v;
      u = dict[LP[0]];
      v = dict[LP[1]];
      if (u != v) {
	fout << u << "\t" << v << endl;
      }
    }
  }
  fout << "];" << endl;
  fout << "resize(g,0);" << endl;
  fout << "resize(g," << nv << ");" << endl;
  fout << "add(g,edges);" << endl;

  fout << "vlabels = {" << endl;
  for(int v=1; v<=nv; ++v) {
    fout << "'" << undict[v] << "'" << endl;
  }
  fout << "};" << endl;
  fout << "label(g,vlabels)" << endl;

  return 0;

}

      
