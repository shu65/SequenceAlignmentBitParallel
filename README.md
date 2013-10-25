SequenceAlignmentBitParallel
============================
SequenceAlignmentBitParallel calculates the score of global alignment by using Bit Parallel algorithm [1].

Build and Run Test
------------------
    cd SequenceAlignmentBitParallel/
    g++ src/PopCount.cpp src/SequenceAlignmentBitParallel.cpp test/SequenceAlignmentBitParallelTest.cpp -lgtest -lgtest_main -lpthread -o SequenceAlignmentBitParallel_test
    ./SequenceAlignmentBitParallel_test

Reference
---------
1. Benson, G., Hernandez, Y., & Loving, J. (2013). A Bit-Parallel, General Integer-Scoring Sequence Alignment Algorithm. Combinatorial Pattern Matching, 1020166, 50–61. 