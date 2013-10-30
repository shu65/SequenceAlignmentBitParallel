SequenceAlignmentBitParallel
============================
SequenceAlignmentBitParallel calculates the score of global alignment by using Bit Parallel algorithm [1].

Status
----------
[![Build Status](https://travis-ci.org/shu65/SequenceAlignmentBitParallel.png?branch=master)](https://travis-ci.org/shu65/SequenceAlignmentBitParallel)
[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/shu65/sequencealignmentbitparallel/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

Build and Run Test
------------------
    cd SequenceAlignmentBitParallel/
    g++ src/PopCount.cpp src/SequenceAlignmentBitParallel.cpp test/SequenceAlignmentBitParallelTest.cpp -lgtest -lgtest_main -lpthread -o SequenceAlignmentBitParallel_test
    ./SequenceAlignmentBitParallel_test

Reference
---------
1. Benson, G., Hernandez, Y., & Loving, J. (2013). A Bit-Parallel, General Integer-Scoring Sequence Alignment Algorithm. Combinatorial Pattern Matching, 1020166, 50–61. 

