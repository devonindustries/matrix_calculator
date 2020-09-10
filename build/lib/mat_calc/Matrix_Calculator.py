import copy

class Matrix():
	'''
	Create an instance of the Matrix class.

	When creating an instance of the matix class, be sure to pass in the dimensions as either a list or a tuple in the form (M, N)
	where the matrix you are creating has dimensions M x N. 
	'''

	# ----------
	# DEFAULT METHODS
	# ----------

	def __init__(self, dimensions, values):
		# Dimensions are of the classic form mxn for a matrix with
		# m rows and n columns, and the values are given as a list from
		# the top left hand corner of the matrix to the bottom right.

		# This is then converted into a matrix in the form of a list, with
		# each column as a list within the matrix

		# Define each attribute before we begin assigning them so that an incomplete matrix is still valid

		self.dimensions = [] # The dimensions of the matrix stored as a list
		self.values = [] # A list containing the values of the matrix in order
		self.entries = 0 # The number of values in the matrix
		self.rows = [] # A list containing the values of the matrix row by row in independent lists
		self.det = None # The determinant of the matrix, if it exists

		# Store the dimensions of the matrix and cast them as integers
		self.dimensions = [int(x) for x in dimensions.split("x")]

		# Check that the dimensions are valid
		if len(self.dimensions) != 2:
			print("Dimensions are invalid. Please make sure you enter in the form #x#")
			return

		self.values = [float(x) for x in values.split(" ")]
		self.entries = self.dimensions[0] * self.dimensions[1]

		# Check that the numer of values is correct
		if len(self.values) != self.entries:
			print("Entered the wrong number of values. Please make sure you entered the correct number of values in the form # # ... #")
			return

		# Store each row as its own list to make matrix operations easier later on
		self.rows = [[self.values[self.dimensions[1] * y + x] for x in range(self.dimensions[1])] for y in range(self.dimensions[0])]

		# Calculate the determinant of the matrix
		self.calculateDet()

	def __str__(self):
		# Return the matrix in value form as a string
		return_string = ""
		for i in self.rows:
			return_string += str(i)
			return_string += "\n"
		return return_string

	def __add__(self, other):
		# Define matrix addition

		# First check that the addition can be carried out
		if self.dimensions == other.dimensions and self.entries == other.entries:
			sum_matrix = self # Declare the matrix to return as a copy of one of the current matrices
			for i in range(self.dimensions[0]):
				for j in range(self.dimensions[1]):
					sum_matrix.rows[i][j] += other.rows[i][j] # Loop through and take the sum of each element
			sum_matrix.refresh() # Now we update the matrix attributes
			return sum_matrix # Return the new matrix
		else:
			return None # If the operation did not succeed, return nothing

	def __sub__(self, other):
		# Define matrix subtraction
		m1 = copy.deepcopy(self)
		m2 = copy.deepcopy(other)
		return m1 + -m2 

	def __pos__(self):
		# Define matrix positive function
		for i in range(self.dimensions[0]):
			for j in range(self.dimensions[1]):
				self.rows[i][j] = abs(self.rows[i][j])
		self.refresh()
		return self

	def __neg__(self):
		# Define matrix negation
		for i in range(self.dimensions[0]):
			for j in range(self.dimensions[1]):
				self.rows[i][j] = (-1) * self.rows[i][j]
		self.refresh()
		return self

	def __mul__(self, other):
		# Define matrix multiplication
		if self.dimensions[1] == other.dimensions[0]:
			# First we calculate our new matrix dimensions, and create a new instance of the matix class
			new_dim = [self.dimensions[0], other.dimensions[1]]
			new_values = new_dim[0] * new_dim[1]
			new_matrix = Matrix('x'.join([str(x) for x in new_dim]),' '.join(["0" for x in range(new_values)]))

			# Now we begin the matrix multiplication
			for i in range(new_matrix.dimensions[0]):
				for j in range(new_matrix.dimensions[1]):
					 # Get the row and column from each matrix respectively
					mRow = self.rows[i]
					mCol = [other.rows[x][j] for x in range(other.dimensions[0])] 

					# Take the sum of the products of each term and assign it to each index in the new matrix
					new_matrix.rows[i][j] = sum([mRow[x] * mCol[x] for x in range(len(mRow))])

			# Refresh the new matrix so that it is accurate
			new_matrix.refresh()

			# Return our new matrix
			return new_matrix

		else:
			# Multiplication for this matrix is not defined
			return None

	# ----------
	# USER-DEFINED METHODS
	# ----------

	def updateValues(self):
		# A self-defined methd that updates the values of the matrix based on the rows
		self.values.clear() # Wipe the current values list

		for i in self.rows:
			self.values.extend(float(x) for x in i) # Extend the values list for each value in each sublist

	def transpose(self):
		# Find the transpose of a matrix
		
		# First we will flip the rows for the transpose matrix
		new_rows = [[] for x in range(self.dimensions[1])]

		for i in range(self.dimensions[1]): # Loop through each column in the new matrix
			for j in range(self.dimensions[0]): # Loop through each row in the new matrix
				new_rows[i].append(self.rows[j][i]) # Flip the index of each value to form the transpose
		self.rows = new_rows

		# Next we need to reverse the dimensions of the matrix
		self.dimensions = self.dimensions[::-1]

		# Finally we can update the matrix attributes
		self.refresh()

	def scale(self, scalar = 1):
		# Multiply all of the entries in the matrix by the value of scalar
		for i in range(self.dimensions[0]):
			for j in range(self.dimensions[1]):
				self.rows[i][j] = scalar * self.rows[i][j]

		# Refresh the matrix so that it is up to date
		self.refresh()

	def calculateDet(self):
		# Calculate the detemrinant of the matrix

		# ----------
		def evaluateDeterminant(sub_matrix):
			# Evaluate the determinant of a matrix

			# First we need to calculate the dimensions of the matrix that we have
			rowcols = len(sub_matrix)
			
			# Next, if the matrix is a 2x2, we can just calculte the determinant to speed up the process
			if rowcols == 2:
				return sub_matrix[0][0] * sub_matrix[1][1] - sub_matrix[0][1] * sub_matrix[1][0]
			if rowcols == 1:
				return sub_matrix[0][0]

			# Finally we can begin the recursion, and calculate the determinant of the matrix
			det = 0

			for i in range(rowcols):
				next_list = [x for x in range(rowcols)] # Find the list corresponding to the iteration of rows
				next_list.pop(i) # Remove the required index from the list
				next_matrix = [[sub_matrix[y][x] for x in range(1, rowcols)] for y in next_list] # Iterate through the new list to find our inner matrix
				det += (-1)**i*sub_matrix[i][0]*evaluateDeterminant(next_matrix) # Add the new matrix value to the determinant

			return det
		# ----------

		# First we need to check that the matrix is square
		if self.isSquare():
			self.det = evaluateDeterminant(self.rows)
		else:
			self.det = None

	def inverse(self):
		# Calculate the inverse of a matrix

		# First check that the matrix is indeed invertible
		if not(self.isInvertible()):
			return None

		# Now we go through the process of finding the inverse using determinants

		# First we need to make a new matrix with the dimensions of the current one
		new_matrix = Matrix('x'.join([str(self.dimensions[0]), str(self.dimensions[1])]), ' '.join(["0" for x in range(self.dimensions[0]*self.dimensions[1])]))

		# Loop through each term in the matrix
		for i in range(new_matrix.dimensions[0]):
			for j in range(new_matrix.dimensions[1]):
				coefficient = (-1)**(i+j) # Determines whether or not the leading coefficient is negative or positive

				# Now we need a list of the rows that we are looping through, discarding the ith row
				# and a list of the columns that we are looping through, discarding the jth column
				# in order to form our new matrix that we will calculate the determinant of

				# To do this, we will first find the rows that we need to loop through to
				# calculate the determinant at the (i, j)th position of that matrix
				rows = [x for x in range(new_matrix.dimensions[0])]
				cols = [y for y in range(new_matrix.dimensions[1])]
				rows.pop(i)
				cols.pop(j)

				# Next, define the matrix of which we will base our determinant
				det_matrix = Matrix('x'.join([str(len(rows)), str(len(cols))]), ' '.join(["0" for x in range((len(rows)) * (len(cols)))]))

				# Now we can append our values to the new determinant matrix
				for x in range(len(rows)):
					for y in range(len(cols)):
						det_matrix.rows[x][y] = self.rows[rows[x]][cols[y]]

				# Refresh our new matrix so that we can get its determinant
				det_matrix.refresh()

				# Assign the determinant value
				new_matrix.rows[i][j] = coefficient * (1 / self.det) * det_matrix.det

		# Transpose the matrix and round it
		new_matrix.transpose()

		# Return the resulting inverse
		return new_matrix

	def invert(self):
		# Invert the instance of the matrix if and only if it is invertible
		if self.inverse():
			self = self.inverse()


	def refresh(self):
		# Refresh any of the attributes of the matrix instance
		self.updateValues()
		self.calculateDet()

	# ----------
	# PROPERTY CHECKS
	# ----------
	# FIRST LIST: Properties checked by internal methods
	# SECOND LIST: Properties checked externally by the user

	def isSquare(self):
		return self.dimensions[0] == self.dimensions[1]

	def isInvertible(self):
		return self.det != 0 and self.isSquare()

	def isIdentity(self):
		# Check if a matrix is an identity matrix
		if self.dimensions[0] == self.dimensions[1]:
			for i in range(self.dimensions[0]):
				for j in range(self.dimensions[1]):
					if i == j and self.rows[i][j] != 1:
						return False
					elif i != j and self.rows[i][j] != 0:
						return False
					else:
						pass
		else:
			return False

		return True

	# ----------

	def isSingular(self):
		return self.det == 0