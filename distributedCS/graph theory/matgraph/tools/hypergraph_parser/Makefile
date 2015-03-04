OBJS = main.o LineParser.o
CXXFLAGS = -ansi -O

hypergraph_parser: $(OBJS)
	g++ $(OBJS) -o hypergraph_parser

.phony: clean
clean:
	rm -f $(OBJS) 

.phony: very-clean
very-clean:
	rm -f $(OBJS) hypergraph_parser