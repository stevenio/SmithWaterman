import unittest
from SmithWaterman import *
#--------------------------------------------------------------------------------
# Testing (You can ignore the following part)
#--------------------------------------------------------------------------------
# abstract test template
class ATC(unittest.TestCase):
  def setUp(self):
    self.SW = SmithWaterman('','')

  def run_diag(self, function_toBeTested,current_cell_indices,neighbour_indices, neighbour_value,expected_score):
    row_index, col_index = current_cell_indices[:]
    neighbour_row_index, neighbour_col_index = neighbour_indices[:]
    self.SW.scoreMatrix[neighbour_row_index,neighbour_col_index] = neighbour_value
    nucleotide1 = self.seq1[row_index-1]
    nucleotide2 = self.seq2[col_index-1]
    score = function_toBeTested(row_index,col_index,nucleotide1,nucleotide2)
    self.assertEqual(score,expected_score)

  def run_off_diag(self, function_toBeTested,current_cell_indices,neighbour_indices, neighbour_value,expected_score):
    row_index, col_index = current_cell_indices[:]
    neighbour_row_index, neighbour_col_index = neighbour_indices[:]
    self.SW.scoreMatrix[neighbour_row_index,neighbour_col_index] = neighbour_value
    score = function_toBeTested(row_index,col_index)
    self.assertEqual(score,expected_score)

#================================================================================
class TestScoreDiag(ATC):
  def setUp(self):
    self.seq1 = 'AG'
    self.seq2 = 'AC'
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.test_function = self.SW._score_diag

  def test_1_1(self):
    current_cell_indices = (1,1)
    neighbour_indices = (0,0)
    neighbour_value = 0
    expected_score = 2
    self.run_diag(self.test_function,current_cell_indices, neighbour_indices, neighbour_value, expected_score)
  
  def test_2_2(self):
    current_cell_indices = (2,2)
    neighbour_indices = (1,1)
    neighbour_value = 2
    expected_score = 1
    self.run_diag(self.test_function,current_cell_indices, neighbour_indices, neighbour_value, expected_score)

  def test_2_1(self):
    current_cell_indices = (2,1)
    neighbour_indices = (1,0)
    neighbour_value = 0
    expected_score = -1
    self.run_diag(self.test_function,current_cell_indices, neighbour_indices, neighbour_value, expected_score)

  def test_1_2(self):
    current_cell_indices = (1,2)
    neighbour_indices = (0,1)
    neighbour_value = 0
    expected_score = -1
    self.run_diag(self.test_function,current_cell_indices, neighbour_indices, neighbour_value, expected_score)


#================================================================================
class TestAlternativeScoreDiag2(ATC):
  def setUp(self):
    self.seq1 = "AGC"
    self.seq2 = "ACA"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.function_toBeTested = self.SW._score_diag

  def test_1_3(self):
    current_cell_indices = (1,3)
    neighbour_indices = (0,2)
    neighbour_value = 0
    expected_score = 2
    self.run_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)
    
  def test_1_2(self):
    current_cell_indices = (1,2)
    neighbour_indices = (0,1)
    neighbour_value = 0
    expected_score = -1
    self.run_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)
    
  def test_2_2(self):
    current_cell_indices = (2,2)
    neighbour_indices = (1,1)
    neighbour_value = 2 
    expected_score = 1
    self.run_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)
    
  def test_2_3(self):
    current_cell_indices = (2,3)
    neighbour_indices = (1,2)
    neighbour_value = 1 
    expected_score = 0
    self.run_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)

#================================================================================
class TestAlternativeScoreVertical(ATC):
  def setUp(self):
    self.seq1 = "AG"
    self.seq2 = "AC"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.function_toBeTested = self.SW._score_vertical

  def test_1_2(self):
    current_cell_indices = (1,2)
    neighbour_indices = (0,2)
    neighbour_value = 0
    expected_score = 0
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)

  def test_2_2(self):
    current_cell_indices = (2,2)
    neighbour_indices = (1,2)
    neighbour_value = 1
    expected_score = 0
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)

#================================================================================
class TestAlternativeScoreVertical2(ATC):
  def setUp(self):
    self.seq1 = "AGC"
    self.seq2 = "ACA"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.function_toBeTested = self.SW._score_vertical

  def test_1_3(self):
    current_cell_indices = (1,3)
    neighbour_indices = (0,3)
    neighbour_value = 0
    expected_score= 0
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)

  def test_2_3(self):
    current_cell_indices = (2,3)
    neighbour_indices = (1,3)
    neighbour_value = 2
    expected_score= 1
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)


#================================================================================
class TestAlternativeScoreHorizontal(ATC):
  def setUp(self):
    self.seq1 = "AG"
    self.seq2 = "AC"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.function_toBeTested = self.SW._score_horizontal

  def test_1_2(self):
    current_cell_indices = (1,2)
    neighbour_indices = (1,1)
    neighbour_value = 2
    expected_score= 1
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)
    
  def test_2_2(self):
    current_cell_indices = (2,2)
    neighbour_indices = (1,2)
    neighbour_value = 1
    expected_score= 0
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)
    

#================================================================================
class TestAlternativeScoreHorizontal2(ATC):
  def setUp(self):
    self.seq1 = "ACA"
    self.seq2 = "AGC"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.function_toBeTested = self.SW._score_horizontal

  def test_1_3(self):
    current_cell_indices = (1,3)
    neighbour_indices = (1,2)
    neighbour_value = 1
    expected_score= 0
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)

  def test_3_3(self):
    current_cell_indices = (3,3)
    neighbour_indices = (3,2)
    neighbour_value = 3
    expected_score= 2
    self.run_off_diag(self.function_toBeTested,current_cell_indices,neighbour_indices,neighbour_value,expected_score)


#================================================================================
class TestScoreMatrix(ATC):
  def setUp(self):
    self.seq1 = "AGC"
    self.seq2 = "ACA"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.function_toBeTested = self.SW._build_scoreMatrix
    self.SW._build_scoreMatrix()
  def test_build_scoreMatrix(self):		
    expected_scoreMatrix = numpy.array([[0,0,0,0],
                                        [0,2,1,2],
	    			        [0,1,1,1],
	    			        [0,0,3,2]
					])
    scoreMatrix = self.SW.scoreMatrix
    print("scoreMatrix",scoreMatrix)
    self.assertTrue((scoreMatrix == expected_scoreMatrix).all())

#================================================================================
class TestScoreMatrix2(ATC):
  def setUp(self):
    self.seq1 = "AGCA"
    self.seq2 = "ACAC"
    self.SW = SmithWaterman(self.seq1, self.seq2)
    self.function_toBeTested = self.SW._build_scoreMatrix
    self.SW._build_scoreMatrix()
  def test_build_scoreMatrix(self):		
    expected_scoreMatrix = numpy.array([[0,0,0,0,0],
                                        [0,2,1,2,1],
	    			        [0,1,1,1,1],
	    			        [0,0,3,2,3],
	    			        [0,2,2,5,4],
					])
    scoreMatrix = self.SW.scoreMatrix
    print("scoreMatrix",scoreMatrix)
    self.assertTrue((scoreMatrix == expected_scoreMatrix).all())

#================================================================================
class TestAlign(unittest.TestCase):
  def setUp(self):
    self.seq1 = "AGCACACA"
    self.seq2 = "ACACACTA"
    self.SW = SmithWaterman(self.seq1, self.seq2)

  def test_align(self):
    expected_seq1 = "AGCACAC-A"
    expected_seq2 = "A-CACACTA"
    output_seq1, output_seq2 = self.SW.align()
    print(self.SW.scoreMatrix)
    print(output_seq1)
    print(output_seq2)
    self.assertEqual(output_seq2, expected_seq2)
    self.assertEqual(output_seq1, expected_seq1)

#================================================================================
class TestAlign2(unittest.TestCase):
  def setUp(self):
    self.seq1 = "TTACCGGCCAACTAA"
    self.seq2 = "ACCGTGTCACTAC"
    self.SW = SmithWaterman(self.seq1, self.seq2)

  def test_align(self):
    expected_seq1 = "ACCG-GCCAACTA"
    expected_seq2 = "ACCGTGTCA-CTA"
    output_seq1, output_seq2 = self.SW.align()
    print(self.SW.scoreMatrix)
    print("output seq1:",output_seq1)
    print("output seq2:",output_seq2)
    self.assertEqual(output_seq2, expected_seq2)
    self.assertEqual(output_seq1, expected_seq1)

#--------------------------------------------------------------------------------
# If this file is run by itself, execute the followings
#--------------------------------------------------------------------------------
if __name__ == '__main__':
  unittest.main(verbosity=1)

