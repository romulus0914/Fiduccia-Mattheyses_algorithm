### HOW TO COMPILE
In this directory, enter following command,
```
make
```
It will generate two executable files, *two-way_min-cut_parition* and *k-way_min-cut_partition*.
If you want to remove executable file, enter following command,
```
make clean
```

### HOW TO RUN
Usage:
```
./<executable_files> <path/to/aux_file>
```
e.g.
```
./two-way_min-cut_partition publiccases/publiccases_basic/ISPD98_ibm01.aux
./k-way_min-cut_partition publiccases/publiccases_advance/ISPD98_ibm01.aux
```

### HOW TO VERIFY
In the directory of publiccases/publiccases_basic or publiccases/publiccases_advance, enter following command,
```
./Verify <aux_file>
```
e.g.
```    
./Verify ISPD98_ibm01.aux
```
