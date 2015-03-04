#!/usr/bin/env /usr/bin/python
"""
A program for converting graffle files into SGF (simple graph format) files.
It reads and parses the XML in a .graffle file and looks for vertices and edges.
(Unfortunately, text blocks look like vertices to graffle, so they become isolated
vertices.) 

Usage: graffle2gf.py < graph.graffle > outputfile

The output is a two-column matrix (two numbers per line) as described in the
matgraph documentation "simple graph format". The output can be saved in a matrix,
call it A, and then we can build a graph from A in matgraph with the command
sgf(g,A)

The embedding of g is (essentially) the same as the embedding inside graffle. 
However, the scale is probably not right. To rescale the embedding in matgraph,
use the "scale" command like this:
scale(g,0.01)
and the embedding will be shrunk by a factor of 100 (which gives pretty good results). 
"""

# read in packages we need
from xml.dom import minidom   # for parsing XML
import sys                    # for reading from the standard input


def converter(doc):
    """
    This is a function for converting various kinds of objects we see 
    inside a graffle document.
    """
    if doc.nodeName == "#text":
	return str(doc.data)
    elif doc.nodeName == "string":
	return str(doc.firstChild.data)
    elif doc.nodeName == 'integer':
	return int(doc.firstChild.data)
    elif doc.nodeName == 'real':
	return float(doc.firstChild.data)
    elif doc.nodeName == 'dict':
	return convert_dict(doc)
    elif doc.nodeName == 'array':
	return convert_list(doc)
    elif doc.nodeName == 'plist':
	return convert_list(doc)
    else:
	return 'unknown:' + doc.nodeName

def convert_dict(doc):
    """
    This takes an XML node that represents a dictionary and converts it 
    into a Python dictionary object
    """
    ans = {}
    if doc.nodeName != 'dict':
	print 'Not a "dict" node!'
	return ans
    list = doc.childNodes
    find_key = True
    for k in range(0,len(list)):
	if find_key:
	    if list[k].nodeName != 'key':
		continue
	    else:
		key = converter(list[k].firstChild);
		find_key = False
		continue
	else:
	    if list[k].nodeName == '#text':
		continue
	    else:
		value = converter(list[k]);
		ans[key] = value
		find_key = True
    return ans


def convert_list(doc):
    """
    This takes a list-like XML object and converts it into a
    Python list
    """
    list = doc.childNodes
    ans = []
    for k in list:
	if k.nodeName != "#text":
	    ans.append(converter(k))
    return ans


def print_dict(d):
    """
    Print out a Python dictionary for viewing. This is for debugging purposes.
    """
    for k in d.keys():
	print str(k) + " --> " + str(d[k])


def print_graph(doc):
    """
    From a parsed XML file (from Graffle), print out the graph structure in 
    SGF format.
    """
    vdict, indexer = find_vertices(doc)
    edges = find_edges(doc,vdict)
    nv = len(vdict)
    ne = len(edges)
    # first row is nv,ne
    print nv,"\t", ne
    # next rows are the edges
    for e in edges:
	print e[0],"\t",e[1]
    # last set of rows are the x,y coordinates of the vertices
    for idx in range(1,nv+1):
	v = indexer[idx]
	x = vdict[v][1]
	y = -vdict[v][2]
	print x,"\t",y


def find_vertices(doc):
    """
    Find all entries that correspond to vertices in the doc. Returns a dictionary 
    that maps the Graffle ID number for each vertex into a tuple (i,x,y) where i
    is an index number (starting from 1) and x,y are the coordinates of the vertex,
    and a second dictionary that maps index numbers 1,2,...,n back to graffle ID
    numbers
    """
    d1 = {}  # v --> (i,x,y) where v is graffle's object number
    d2 = {}  # i --> v
    n = 0
    dlist = doc.getElementsByTagName('dict')
    for d in dlist:
	dd = convert_dict(d)
	if dd.has_key('Class'):
	    if dd['Class'] == 'ShapedGraphic':
		id = dd['ID']
		n = n+1
		# get the coordinates and convert them into integers
		xylist = dd['Bounds'].split(',')
		for k in range(0,4):
		    xylist[k] = float(xylist[k].strip('{} '))
		x = int(xylist[0])
		y = int(xylist[1])
		d1[id] = (n,x,y)
		d2[n] = id
    return d1,d2

def find_edges(doc,vdict):
    """ 
    Find all the edges in the graph. vdict is the first output of
    find_vertices(doc).  returns a list of tuples for the edges.
    """
    ans = []
    dlist = doc.getElementsByTagName('dict')
    for d in dlist:
	dd = convert_dict(d)
	if dd.has_key('Class'):
	    if dd['Class'] == 'LineGraphic':
		u = dd['Head']['ID']
		v = dd['Tail']['ID']
		uu = vdict[u][0]
		vv = vdict[v][0]
		if uu != vv:
		    ans.append( (uu,vv) )
    return ans

# read the standard input and find all the graph element entries
doc = minidom.parse(sys.stdin)
print_graph(doc)
