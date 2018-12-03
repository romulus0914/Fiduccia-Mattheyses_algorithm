GXX = g++
GXXFLAG = -std=c++11

FLAG = -Wall -O3

SRC_BASIC = two-way_min-cut_partition.cpp
SRC_ADVANCE = k-way_min-cut_partition.cpp
EXE_BASIC = two-way_min-cut_partition
EXE_ADVANCE = k-way_min-cut_partition

all: basic advance

basic: $(SRC_BASIC)
	$(GXX) $(GXXFLAG) $(FLAG) -o $(EXE_BASIC) $(SRC_BASIC)

advance: $(SRC_ADVANCE)
	$(GXX) $(GXXFLAG) $(FLAG) -o $(EXE_ADVANCE) $(SRC_ADVANCE)

clean:
	rm -rf $(EXE_BASIC) $(EXE_ADVANCE)
