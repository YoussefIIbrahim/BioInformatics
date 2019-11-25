import sys

import io
import unittest
import unittest.mock
from unittest import TestCase
import NeedlemanWunschAlgorithm as nw

class TestNeedlemanWunsch(TestCase):

    # def test_initScript(self):
    #     self.assertRaises(AssertionError, nw.InitScript())

    def test_readFile(self):
        nw.readFile('config.txt', 0)
        nw.readFile('sequenceA.txt', 1)
        nw.readFile('sequenceB.txt', 1)

    def test_zeros(self):
        actualList = nw.zeros(3, 3)
        result = False
        if len(actualList) > 0:
            result = all(elem == actualList[0] for elem in actualList)

        assert result

    def test_match_score(self):
        gap_penalty = -5
        match = 10
        mismatch = -10

        alpha = 'A'
        beta = 'B'
        returnValue = nw.match_score(alpha, beta, gap_penalty, match, mismatch)
        assert returnValue == mismatch

        beta = '-'
        returnValue = nw.match_score(alpha, beta, gap_penalty, match, mismatch)
        assert returnValue == gap_penalty

        alpha = '-'
        beta = 'U'
        returnValue = nw.match_score(alpha, beta, gap_penalty, match, mismatch)
        assert returnValue == gap_penalty

        alpha = 'U'
        beta = '-'
        returnValue = nw.match_score(alpha, beta, gap_penalty, match, mismatch)
        assert returnValue == gap_penalty

        beta = 'A'
        alpha = 'A'
        returnValue = nw.match_score(alpha, beta, gap_penalty, match, mismatch)
        assert returnValue == match


    def test_needleman_wunsch(self):
        SEQ1 = 'MARS'
        SEQ2 = 'SMART'
        GAP_PENALTY = -2
        MATCH = 5
        MISMATCH = -5
        MAX_PATHS = 20
        output = nw.needleman_wunsch(SEQ1, SEQ2, GAP_PENALTY, MATCH, MISMATCH, MAX_PATHS)
        outputString = 'Score: \n ' + str(output[len(output) - 1]) + '\n Alignments: \n'
        for i in range(len(output)):
            if i != len(output) - 1:
                for j in range(len(output[i]) - 1):
                    outputString += '(' + output[i][j] + ', ' + output[i][j + 1] + ') \n'
        print(outputString)
        assert outputString == "Score: \n 9\n Alignments: \n(-MAR-S, SMART-) \n(-MARS-, SMAR-T) \n"

