"""
Smith-Waterman algorithm for sequence alignment
Author: Yuhang Wang
Date: Sep-10-2014
Copyright: CC-BY-4.0 URL: http://creativecommons.org/licenses/by/4.0/
"""
#--------------------------------------------------------------------------------
# Compatibility with Python 3
#--------------------------------------------------------------------------------
from __future__ import print_function, division
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Module dependencies
#--------------------------------------------------------------------------------
import unittest
import numpy
#--------------------------------------------------------------------------------




#--------------------------------------------------------------------------------
# Smith-Waterman algorithm (wrapped in a python class)
#--------------------------------------------------------------------------------
class SmithWaterman(object):
  def __init__(self,seq1, seq2, match_score=2, mismatch_score=-1):
    self.seq1 = seq1
    self.seq2 = seq2
    self.match_score = match_score
    self.mismatch_score = mismatch_score
    self.numCols = len(seq1) # number of columns of the scoring matrix
    self.numRows = len(seq2) # number of rows of the scoring matrix
    # note: the score matrix always has one more row/column than the size of 
    #   the seq1/seq2
    self.scoreMatrix = numpy.zeros((self.numRows+1, self.numCols+1), dtype=numpy.int8)
  
  def _similarity_function(self, nucleotide1, nucleotide2):
    if nucleotide1 == nucleotide2:
      return self.match_score
    else:
      return self.mismatch_score

  def _alternative_score_1(self, nucleotide1, nucleotide2):
    """
    The 1st alternative score for scoreMatirx[i,j] (always set to "0")
    """
    return 0

  def _alternative_score_2(self, row_index, col_index, nucleotide1, nucleotide2):
    """
    The 2nd alternative score for scoreMatirx[i,j]
    """
    try:
      assert row_index > 0
      assert col_index > 0
    except:
      msg = []
      msg.append("\nINTERNAL ERROR: when computing _alternative_score_2, index out of bound")
      msg.append("Hint: must keep row_index > 0 an col_index >0,")
      msg.append("where these indices are for Smith-Waterman score matrix")
      print("".join(msg))
      raise IndexError

    score  = self.scoreMatrix[row_index-1,col_index-1] 
    score += self._similarity_function(nucleotide1, nucleotide2)
    return score

  def _alternative_score_3(self, row_index, col_index):
    """
    The 3rd alternative score for scoreMatrix[i,j]
    """
    _k = 1
    _tmp = []
    gap_penalty = -row_index
    while self.scoreMatrix[row_index-_k, col_index] - gap_penalty > 0 and row_index - _k >=0 :
      score = self.scoreMatrix[row_index-_k, col_index] + gap_penalty
      _tmp.append(score)
      _k += 1
    return max(_tmp)

  def _build_scoreMatrix(self):
    """
    Swith-Waterman alignment algorithm 
    is implemented here!
    """
    for _i in range(1, self.numRows):
      for _j in range(1, self.numCols):
        nucleotide1 = self.seq1[_i-1]
        score1 = _alternative_score_1(nucleotide1,nucleotide2)
	
    return



  def __repr__(self):
    _show_str = []
    _show_str.append("seq1: {0}\n".format(self.seq1))
    _show_str.append("seq2: {0}\n".format(self.seq2))
    _show_str.append("Score(match)={0}\n".format(self.match_score))
    _show_str.append("Score(mismatch)={0}\n".format(self.mismatch_score))
    _show_str.append("Score matrix H:\n {0}\n".format(self.scoreMatrix))
    output = "".join(_show_str)
    print(output)
    return output




#--------------------------------------------------------------------------------
# Testing (You can ignore the following part)
#--------------------------------------------------------------------------------
class TestSmithWaterman(unittest.TestCase):
  def setUp(self):
    self.seq1 = "AT"
    self.seq2 = "AC"
    self.nucleotide1 = "A"
    self.nucleotide2 = "G"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.match_score = self.SW.match_score
    self.mismatch_score = self.SW.mismatch_score
    self.scoreMatrix = self.SW.scoreMatrix
    self.expected_score1 = 0

  def test_alternative_core_1(self):
    score = self.SW._alternative_score_1(self.nucleotide1, self.nucleotide2)
    expected_score = self.expected_score1
    self.assertEqual(score, expected_score)

  def test_alternative_score_2(self):
    row_index = 1
    col_index = 1

    # case: nucleotide1 != nucleotide2
    score = self.SW._alternative_score_2(row_index, col_index, self.nucleotide1, self.nucleotide2)
    expected_score = self.mismatch_score
    self.assertEqual(score, expected_score)

    # case: nucleotide1 == nucleotide2
    score = self.SW._alternative_score_2(row_index, col_index, self.nucleotide1, self.nucleotide1)
    expected_score = self.match_score
    self.assertEqual(score, expected_score)

    # test fail case
    row_index = 0
    col_index = 1
    self.assertRaises(IndexError, 
         self.SW._alternative_score_2, row_index, col_index, self.nucleotide1, self.nucleotide2)

  def test_alternative_score_3(self):
    row_index = 1
    col_index = 1
    score = self.SW._alternative_score_3(row_index, col_index)
    expected_score = -row_index
    self.assertEqual(score, expected_score)
    
    
  def test_repr(self):
    _show_str = []
    _show_str.append("PARAMETERS:")
    _show_str.append("seq1: {0}\n".format(self.seq1))
    _show_str.append("seq2: {0}\n".format(self.seq2))
    _show_str.append("Score(match)={0}\n".format(self.match_score))
    _show_str.append("Score(mismatch)={0}\n".format(self.mismatch_score))
    _show_str.append("Score matrix H: {0}\n".format(self.scoreMatrix))
    show_str = "".join(_show_str)
    self.assertTrue(self.SW.__repr__(), show_str)


#--------------------------------------------------------------------------------
# If this file is run by itself, execute the followings
#--------------------------------------------------------------------------------
if __name__ == '__main__':
  unittest.main()

