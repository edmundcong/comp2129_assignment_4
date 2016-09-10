make
time .pagerank 4 < tests/test12.in
gprof -l ./pagerank > pageranl.stats
less pagerank.stats
