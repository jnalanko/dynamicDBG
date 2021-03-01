.PHONY: fdbg-recsplit-build fdbg-recsplit-add fdbg-recsplit-del fdbg-recsplit-query
all: fdbg-recsplit-build fdbg-recsplit-add fdbg-recsplit-del fdbg-recsplit-query

test: BuildDataStructure.cpp formatutil.cpp TestUtil.cpp FDBG.cpp
	g++ -O0 -g  -std=c++17 BuildDataStructure.cpp -o bds_test -l boost_log -l boost_system -l pthread -I ../sux

fdbg-recsplit-build:
	g++ -O3 -std=c++17 fdbg-build.cpp -o fdbg-recsplit-build -l boost_log -l boost_system -l pthread -I ../sux -DNDEBUG -g -march=native

fdbg-recsplit-add:
	g++ -O3 -std=c++17 fdbg-add.cpp -o fdbg-recsplit-add -l boost_log -l boost_system -l pthread -I ../sux -DNDEBUG -g -march=native

fdbg-recsplit-del:
	g++ -O3 -std=c++17 fdbg-del.cpp -o fdbg-recsplit-del -l boost_log -l boost_system -l pthread -I ../sux -DNDEBUG -g -march=native

fdbg-recsplit-query:
	g++ -O3 -std=c++17 fdbg-query.cpp -o fdbg-recsplit-query -l boost_log -l boost_system -l pthread -I ../sux -DNDEBUG -g -march=native
