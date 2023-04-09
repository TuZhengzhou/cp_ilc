libsnark merkle circuit example

The example shows how to generate proof for one merkle path on one merkle tree with depth 3.

1/ init 
 ```
 git submodule update --init --recursive
 ```
2/ compile
 ```
 mkdir build; cd build; cmake ..; make
 ```
 You can find the "merkle" binary under the merkle folder.

3/ setup
```
./merkle setup
```

4/ prove
```
./merkle prove [data1] [data2] [data3] [data4] [data5] [data6] [data7] [data8] [index]
```
Record down the root information, which is used on verify.

for example:
run ./merkle prove 1 2 3 4 5 6 7 8 7
will output: root is f171e00bb40c83de1f09c64e2cc4558e3c327aa9e8525a467c83576071bc1045

5/ verify
```
./merkle verify [root]
```
