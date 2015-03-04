OBJS = main.o LineParser.o
CXXFLAGS = -ansi -O

graph_parser: $(OBJS)
	g++ $(OBJS) -o graph_parser

.phony: clean
clean:
	rm -f $(OBJS) parsed_graph.m

.phony: very-clean
very-clean:
	rm -f $(OBJS) parsed_graph.m graph_parser