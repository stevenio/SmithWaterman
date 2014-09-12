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
    self.numRows = len(seq1)+1 # number of rows of the scoring matrix
    self.numCols = len(seq2)+1 # number of columns of the scoring matrix
    # note: the score matrix always has one more row/column than the size of 
    #   the seq1/seq2
    self.scoreMatrix = numpy.zeros((self.numRows, self.numCols), dtype=numpy.int8)
    self.predecessorIndexMatrix = numpy.zeros((self.numRows, self.numCols, 2), dtype=numpy.int8)
    self.aligned_seq1 = []
    self.aligned_seq2 = []
    self.max_score_location = (0,0) # indices for matrix element with max score
    self.max_score = 0
    self.gap_symbol = '-'

  def _gap_penalty(self,k):
    return -k

  def _similarity_function(self, nucleotide1, nucleotide2):
    if nucleotide1 == nucleotide2:
      return self.match_score
    else:
      return self.mismatch_score

  def _score_diag(self, row_index, col_index, nucleotide1, nucleotide2):
    """
    The 2nd alternative score for scoreMatirx[i,j]
    """
    try:
      assert row_index > 0
      assert col_index > 0
    except:
      msg = []
      msg.append("\nINTERNAL ERROR: when computing _score_diag, index out of bound")
      msg.append("Hint: must keep row_index > 0 an col_index >0,")
      msg.append("where these indices are for Smith-Waterman score matrix")
      raise IndexError

    score  = self.scoreMatrix[row_index-1,col_index-1] 
    score += self._similarity_function(nucleotide1, nucleotide2)
    return score

  def _score_vertical(self, row_index, col_index):
    """
    The 3rd alternative score for scoreMatrix[i,j]
    i.e. Vertical score
    """
    _k = 1
    _tmp = [0]
    while row_index - _k >= 0 :
      gap_penalty = self._gap_penalty(_k)
      score = self.scoreMatrix[row_index-_k, col_index] + gap_penalty
      if score < 0:
        break
      else:
	_tmp.append(score)
	_k += 1
    return max(_tmp)
  
  def _score_horizontal(self, row_index, col_index):
    """
    The 4h alternative score for scoreMatrix[i,j]
    i.e. Horizontal score
    """
    _k = 1
    _tmp = [0]
    while col_index - _k >= 0:
      gap_penalty = self._gap_penalty(_k)
      score = self.scoreMatrix[row_index,col_index-_k] + gap_penalty
      if score < 0:
        break
      else:
	_tmp.append(score)
	_k += 1
    return max(_tmp)
  
  def _calc_predecessor_indices(self, row_index, col_index, scores, offsets):
    """
    Find the indices for the current cell's predessor
    Strategy: start from the diagonal neighour and compare agains its 
       vertical and horizontal neighbours.
    Inputs: row_index :: row index of current cell
            col_index :: col index of current cell
	    scores :: array of scores for the 3 neighbours
	    offsets :: index offsets of these 3 neighbours relative to the current cell
    """
    # use the first elment in scores to initialize outputs
    _i = 0
    max_score = scores[_i]
    row_offset = offsets[_i][0]
    col_offset = offsets[_i][1]
    output_row_index = row_index + row_offset
    output_col_index = col_index + col_offset
    # predecessor cell's score
    pre_cell_score = self.scoreMatrix[output_row_index, output_col_index]
    for _i in range(1,len(scores)):
      tmp_rowId = row_index + offsets[_i][0]
      tmp_colId = col_index + offsets[_i][1]
      # alternative predecessor cell's score
      tmp_cell_score = self.scoreMatrix[tmp_rowId, tmp_colId]

      criterion1 = max_score < scores[_i] 
      # if there is a tie in score contribution, the use the predecessor whose own score is higher
      criterion2 = (max_score == scores[_i] and tmp_cell_score > pre_cell_score)
      if criterion1 or criterion2 :
        max_score = scores[_i]
	row_offset = offsets[_i][0]
	col_offset = offsets[_i][1]
	output_row_index = row_index + row_offset
	output_col_index = col_index + col_offset
	# predecessor cell's score
	pre_cell_score = self.scoreMatrix[output_row_index, output_col_index]

    return (output_row_index, output_col_index, max_score)


  def _build_scoreMatrix(self):
    """
    Build Swith-Waterman alignment matrix
    """
    for _i in range(1, self.numRows):
      for _j in range(1, self.numCols):
        nucleotide1 = self.seq1[_i-1]
	assert _j-1 < self.numCols -1
	assert _j-1 >= 0
	nucleotide2 = self.seq2[_j-1]
	#print("nt1: {0}; nt2: {1}".format(nucleotide1, nucleotide2))
	score0 = self._score_diag(_i, _j, nucleotide1, nucleotide2)
	score1 = self._score_vertical(_i, _j)
	score2 = self._score_horizontal(_i, _j)

	scores = (score0, score1, score2)
	offsets = [(-1,-1), (-1,0), (0,-1)]
	pre_row_index, pre_col_index, max_score = self._calc_predecessor_indices(_i,_j,scores,offsets)
	self.scoreMatrix[_i,_j] =  max_score
	self.predecessorIndexMatrix[_i,_j,0] = pre_row_index
	self.predecessorIndexMatrix[_i,_j,1] = pre_col_index

	# update global max score and its location
	if max_score > self.max_score:
	  self.max_score = max_score
	  self.max_score_location = (_i,_j)

    return

  def align(self):
    self._build_scoreMatrix()
    row_index, col_index = self.max_score_location
    pre_row_index, pre_col_index = self.predecessorIndexMatrix[row_index,col_index,:]
    

    while self.scoreMatrix[row_index,col_index] > 0:
      delta_row_index = row_index - pre_row_index
      delta_col_index = col_index - pre_col_index
      Ttype = '' # transition type
      if delta_row_index == 1 and delta_col_index == 1:
	# note: the row_index/col_index has an offset by 1 compared with
	#   the given sequences
	self.aligned_seq1.append(self.seq1[row_index-1])
	self.aligned_seq2.append(self.seq2[col_index-1])
	Ttype = 'Diagonal'
      elif delta_row_index == 1 and delta_col_index == 0:
	# a gap in seq2
	self.aligned_seq1.append(self.seq1[row_index-1])
	self.aligned_seq2.append(self.gap_symbol)
	Ttype = "Vertical (gap in seq2)"
      elif delta_row_index == 0 and delta_col_index == 1:
	# a gap in seq1
	self.aligned_seq1.append(self.gap_symbol)
	self.aligned_seq2.append(self.seq2[col_index-1])
	Ttype = "Horizontal (gap in seq1)"

      print("({0},{1}) --> ({2},{3})".format(row_index, col_index,pre_row_index,pre_col_index), Ttype)
      # Now, update indices
      row_index, col_index = self.predecessorIndexMatrix[row_index,col_index,:]
      pre_row_index, pre_col_index = self.predecessorIndexMatrix[row_index,col_index,:]
    aligned_seq1 = "".join(reversed(self.aligned_seq1))
    aligned_seq2 = "".join(reversed(self.aligned_seq2))
    return [aligned_seq1, aligned_seq2]

    # if the location of the cell with max score is not the lower right corner
    # then just append the unused part to the end
    #if row_index < self.numRows - 1:
    #  self.aligned_seq1.append(self.seq1[row_index+1:])
    #if col_index < self.numCols - 1:
    #  self.aligned_seq2.append(self.seq2[col_index+1:])


  def __repr__(self):
    _show_str = []
    _show_str.append("seq1: {0}\n".format(self.seq1))
    _show_str.append("seq2: {0}\n".format(self.seq2))
    _show_str.append("Score(match)={0}\n".format(self.match_score))
    _show_str.append("Score(mismatch)={0}\n".format(self.mismatch_score))
    _show_str.append("Score matrix H:\n {0}\n".format(self.scoreMatrix))
    output = "".join(_show_str)
    return output

#--------------------------------------------------------------------------------
# If this file is run by itself, execute the followings
#--------------------------------------------------------------------------------
if __name__ == '__main__':
  ## test 1
  seq1 = "ACCGGCCAACTA"
  seq2 = "ACCGTGTCACTA"
  sw = SmithWaterman(seq1, seq2)
  aligned_seq1, aligned_seq2 = sw.align()
  print("\nResults:")
  print("score matrix:")
  print(sw.scoreMatrix)
  print("sequence 1: {0}".format(aligned_seq1))
  print("sequence 2: {0}".format(aligned_seq2))



