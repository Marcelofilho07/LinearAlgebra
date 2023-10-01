#include <array>
#include <iostream>
#include <vector>

void PrintMatrix(const std::vector<std::vector<float>>& matrix)
{
	std::cout <<" --------------PRINTING-MATRIX----------------" << std::endl;
	
	for (const std::vector<float>& r : matrix)
	{
		for(const float c : r)
		{
			std::cout << c << " ";
		}
		std::cout << std::endl;
	}
	std::cout <<" ---------------------------------------------" << std::endl;
}

float operator/ (const std::vector<std::vector<float>>& matrixA, const std::vector<std::vector<float>>& matrixB)
{
	return matrixA[0][0] / matrixB[0][0];
}


std::vector<std::vector<float>> operator- (const std::vector<std::vector<float>>& matrixA, const std::vector<std::vector<float>>& matrixB)
{
	std::vector<float> newVec;
	newVec.reserve(matrixA.size());
	for(int i = 0; i < matrixA.size(); ++i)
	{
		newVec.push_back(matrixA[0][i] - matrixB[0][i]);
	}
	
	return std::vector<std::vector<float>>{newVec};
}

std::vector<std::vector<float>> MultiMat(std::vector<std::vector<float>> A, std::vector<std::vector<float>> B)
{
	std::vector<std::vector<float>> returnMatrix;
	
	const unsigned long long rows = A.size();
	const unsigned long long innerRows = A[0].size();
	const unsigned long long cols = B[0].size();
	for(int i = 0; i < rows; ++i)
	{
		returnMatrix.emplace_back();
		for(int j = 0; j < cols; ++j)
		{
			returnMatrix[i].push_back(0);
			for(int k = 0; k < innerRows; ++k)
			{
				returnMatrix[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return returnMatrix;
}

std::vector<std::vector<float>> GenerateTranspose(std::vector<std::vector<float>>& matrix)
{
	std::vector<std::vector<float>> transpose;
	int m = matrix[0].size();
	int n = matrix.size();
	for(int i = 0; i < m; ++i)
	{
		std::vector<float> newTransposedRow;
		for(int j = 0; j < n; ++j)
		{
			newTransposedRow.push_back(matrix[j][i]);
		}

		transpose.push_back(newTransposedRow);
	}

	return transpose;
}

std::vector<std::vector<float>> GeneratePermutationMatrix(int firstRow, int secondRow, int size)
{
	std::vector<std::vector<float>> MatrixToReturn;

	for(int i = 0; i < size; ++i)
	{
		std::vector<float> newRow;
		if(i == firstRow)
		{
			for(int j = 0; j < size; ++j)
			{
				if(secondRow == j)
				{
					newRow.push_back(1.f);
				}
				else
				{
					newRow.push_back(0.f);
				}
			}
		}
		else if (i == secondRow)
		{
			for(int j = 0; j < size; ++j)
			{
				if(firstRow == j)
				{
					newRow.push_back(1.f);
				}
				else
				{
					newRow.push_back(0.f);
				}
			}
		}
		else
		{
			for(int j = 0; j < size; ++j)
			{
				if(i == j)
				{
					newRow.push_back(1.f);
				}
				else
				{
					newRow.push_back(0.f);
				}
			}
		}

		MatrixToReturn.push_back(newRow);
	}

	return MatrixToReturn;
}

std::vector<std::vector<float>> GenerateIdentityMatrix(int size)
{
	return GeneratePermutationMatrix(-1,-1, size);
}

void BackSubstitution(std::vector<std::vector<float>>& matrixA, std::vector<std::vector<float>>& matrixB, bool showSteps = false)
{
	int m = matrixA.size();
	int n = matrixA[0].size();

	for (int i = m - 1; i > 0; --i)
	{
		for(int j = (i-1); j >= 0; --j)
		{
			std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(n)};
			
			float divMatrix = matrixA[j][i] / matrixA[i][i];
			
			identityMatrix[j][i] = -divMatrix;
			
			matrixA = MultiMat(identityMatrix, matrixA);
			matrixB = MultiMat(identityMatrix, matrixB);
		}

		
		std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(n)};
			
		identityMatrix[i][i] = 1.f/matrixA[i][i];
		
		matrixA = MultiMat(identityMatrix, matrixA);
		matrixB = MultiMat(identityMatrix, matrixB);
	}

		PrintMatrix(matrixA);
		PrintMatrix(matrixB);
}

void LinearCombination(std::vector<std::vector<float>>& matrix, std::vector<std::vector<float>>& b, bool showSteps = false)
{
	/*
	 * Ax = b
	 * Linear combination of the columns
	 *
	 * Can I solve Ax = b for every b?
	 * Do the linear comb of the columns fill 3D space?
	 * Yes, if good matrix ;) (non-singular, Invertible)
	 * PA = LU
	 *
	 * U x = L B
	 * Steps:
	 * 1) 
	 * x + 2y + z = 2
	 * 3x + 8y + z = 12
	 * 4y + z = 2
	 *
	 * 1 2 1 x 2
	 * 3 8 1 y 12
	 * 0 4 1 z 2
	 *
	 * 3x3
	 * 1x0
	 * 2x0
	 * 2x1
	 * 
	 */

	constexpr int row = 3;
	constexpr int col = 3;
	//std::vector<std::vector<float>> matrix {std::vector<float>{1.f,2.f,1.f},std::vector<float>{3.f,8.f,1.f},std::vector<float>{0.f,4.f,1.f}};
	//std::vector<std::vector<float>> matrix {std::vector<float>{2.f,1.f,3.f},std::vector<float>{15.f,2.f,0.f},std::vector<float>{1.f,3.f,1.f}};
	
	//std::vector<std::vector<float>> b {std::vector<float>{2.f},std::vector<float>{12.f},std::vector<float>{2.f}};
	
	PrintMatrix(matrix);

	std::vector<std::vector<float>> LMatrix {std::vector<float>{1.f,0.f,0.f},std::vector<float>{0.f,1.f,0.f},std::vector<float>{0.f,0.f,1.f}};
	float tol = 0.00001f;
	for (int i = 0; i < col-1; ++i)
	{
		if(std::abs(matrix[i][i]) < tol)
		{
			int rowToSwap = i;
			float maxScalar = matrix[i][i];
			for(int j = (i+1); j < row; ++j)
			{
				if(std::abs(matrix[j][i]) > maxScalar)
				{
					rowToSwap = j;
				}
				maxScalar = std::max(maxScalar, std::abs(matrix[j][i]));
			}
		
			matrix = MultiMat(GeneratePermutationMatrix(i,rowToSwap,col), matrix);
		}
		
		for(int j = (i+1); j < row; ++j)
		{
			std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(row)};
			
			float divMatrix = matrix[j][i] / matrix[i][i];
			
			identityMatrix[j][i] = -divMatrix;
			LMatrix[j][i] = divMatrix;
			
			matrix = MultiMat(identityMatrix, matrix);
			b = MultiMat(identityMatrix, b);
		}
	}
	
	PrintMatrix(matrix);
	PrintMatrix(b);
}

void GaussEliminationPartialPivot(std::vector<std::vector<float>>& matrix)
{
	const unsigned int row = matrix.size();
	const unsigned int col = matrix[0].size();
	
	std::vector<std::vector<float>> b {std::vector<float>{2.f},std::vector<float>{12.f},std::vector<float>{2.f}};
	
	PrintMatrix(matrix);

	std::vector<std::vector<float>> LMatrix {std::vector<float>{1.f,0.f,0.f},std::vector<float>{0.f,1.f,0.f},std::vector<float>{0.f,0.f,1.f}};
	float tol = 0.00001f;
	for (unsigned int i = 0; i < col-1; ++i)
	{
		int rowToSwap = i;
		float maxScalar = matrix[i][i];
		for(unsigned int j = (i+1); j < row; ++j)
		{
			if(std::abs(matrix[j][i]) > maxScalar)
			{
				rowToSwap = j;
			}
			maxScalar = std::max(maxScalar, std::abs(matrix[j][i]));
		}
	
		matrix = MultiMat(GeneratePermutationMatrix(i,rowToSwap,col), matrix);
		
		for(unsigned int j = (i+1); j < row; ++j)
		{
			std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(row)};

			const float divMatrix = matrix[j][i] / matrix[i][i];
			
			identityMatrix[j][i] = -divMatrix;
			LMatrix[j][i] = divMatrix;
			
			matrix = MultiMat(identityMatrix, matrix);
			b = MultiMat(identityMatrix, b);
		}
	}
	
	PrintMatrix(matrix);
}

void GaussEliminationCompletePivot(std::vector<std::vector<float>>& matrix)
{
	const unsigned int row = matrix.size();
	const unsigned int col = matrix[0].size();
	
	std::vector<std::vector<float>> b {std::vector<float>{2.f},std::vector<float>{12.f},std::vector<float>{2.f}};
	
	PrintMatrix(matrix);

	std::vector<std::vector<float>> LMatrix {std::vector<float>{1.f,0.f,0.f},std::vector<float>{0.f,1.f,0.f},std::vector<float>{0.f,0.f,1.f}};
	float tol = 0.00001f;
	for (unsigned int i = 0; i < col-1; ++i)
	{
		int rowToSwap = i;
		int colToSwap = i;
		float maxScalar = matrix[i][i];
		for(unsigned int j = (i+1); j < row; ++j)
		{
			if(std::abs(matrix[j][i]) > maxScalar)
			{
				rowToSwap = j;
			}
			maxScalar = std::max(maxScalar, std::abs(matrix[j][i]));
		}
		for(unsigned int j = (i+1); j < col; ++j)
		{
			if(std::abs(matrix[i][j]) > maxScalar)
			{
				colToSwap = j;
			}
			maxScalar = std::max(maxScalar, std::abs(matrix[i][j]));
		}

		if (colToSwap > rowToSwap)
		{
			matrix = MultiMat(matrix,GeneratePermutationMatrix(i,colToSwap,col));
			//Save col change for back row substitution.
		}
		else
		{
			matrix = MultiMat(GeneratePermutationMatrix(i,rowToSwap,col), matrix);
		}
		
		for(unsigned int j = (i+1); j < row; ++j)
		{
			std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(row)};

			const float divMatrix = matrix[j][i] / matrix[i][i];
			
			identityMatrix[j][i] = -divMatrix;
			LMatrix[j][i] = divMatrix;
			
			matrix = MultiMat(identityMatrix, matrix);
			b = MultiMat(identityMatrix, b);
		}
	}
	
	PrintMatrix(matrix);
}

void EchelonForm(std::vector<std::vector<float>>& matrixA, std::vector<std::vector<float>>& matrixB, bool showSteps = false)
{
	const unsigned int row = matrixA.size();
	const unsigned int col = matrixA[0].size();
	
	
	PrintMatrix(matrixA);
	std::vector<std::pair<int, int>> pivots;
	std::vector<std::vector<float>> LMatrix {std::vector<float>{1.f,0.f,0.f},std::vector<float>{0.f,1.f,0.f},std::vector<float>{0.f,0.f,1.f}};

	for (unsigned int i = 0, k = 0; i < row-1  && k < col-1; ++i, ++k)
	{
		while(matrixA[i][k] == 0.f && k < col)
		{
			++k;
			if(matrixA[i][k] == 0.f)
			{
				int rowToSwap = i;
				float maxScalar = matrixA[i][i];
				for(unsigned int j = (i+1); j < row; ++j)
				{
					if(std::abs(matrixA[j][i]) > maxScalar)
					{
						rowToSwap = j;
					}
					maxScalar = std::max(maxScalar, std::abs(matrixA[j][i]));
				}
	
				matrixA = MultiMat(GeneratePermutationMatrix(i,rowToSwap,row), matrixA);
			}
		}
		if(matrixA[i][k] != 0.f)
		{
			pivots.emplace_back(i,k);
		
			for(unsigned int j = (i+1); j < row; ++j)
			{
				std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(row)};

				const float divMatrix = matrixA[j][k] / matrixA[i][k];
		
				identityMatrix[j][i] = -divMatrix;
				LMatrix[j][i] = divMatrix;
		
				matrixA = MultiMat(identityMatrix, matrixA);
				matrixB = MultiMat(identityMatrix, matrixB);
			}
		}
	}

	if(showSteps)
	{
		PrintMatrix(matrixA);
	}

	for(auto p : pivots)
	{
		if(p.first > 0)
		{
			for(int j = (p.first+1); j >= 0; --j)
			{
				std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(row)};
	
				const float divMatrix = matrixA[j][p.second] / matrixA[p.first][p.second];
		
				identityMatrix[j][p.first] = -divMatrix;
				LMatrix[j][p.first] = divMatrix;
		
				matrixA = MultiMat(identityMatrix, matrixA);
				matrixB= MultiMat(identityMatrix, matrixB);
			}
		}
		
		std::vector<std::vector<float>> identityMatrix {GenerateIdentityMatrix(row)};
	
		identityMatrix[p.first][p.first] = 1/matrixA[p.first][p.second];
		matrixA = MultiMat(identityMatrix, matrixA);
		matrixB = MultiMat(identityMatrix, matrixB);
		
	}
	
	PrintMatrix(matrixA);
	PrintMatrix(matrixB);
}


void LeastSquares(std::vector<std::vector<float>>& matrixA, std::vector<std::vector<float>>& matrixB)
{
	//AT A X^ = AT b

	std::vector<std::vector<float>> matrixATA = MultiMat(GenerateTranspose(matrixA), matrixA);
	std::vector<std::vector<float>> matrixATB = MultiMat(GenerateTranspose(matrixA), matrixB);
	PrintMatrix(matrixATA);
	PrintMatrix(matrixATB);

	
	EchelonForm(matrixATA, matrixATB, true);

	//this prints x^
	BackSubstitution(matrixATA,matrixATB);

	//this prints P
	PrintMatrix(MultiMat(matrixA, matrixATB));
}

void GramSchmidt(std::vector<std::vector<float>>& matrixA)
{
	//B = b - (AT b / AT A) A
	//C = c - (AT c / AT A) A - (BT c/ BT B) B
	// Work in progress

	std::vector<std::vector<float>> orthoVectors;
	std::vector<std::vector<float>> matrixQ;
	for(auto vec : matrixA)
	{
		std::vector<float> newVec = vec;
		for(auto ortho : orthoVectors)
		{
			std::vector<std::vector<float>> A = {ortho};
			std::vector<std::vector<float>> b = {vec};
			std::vector<std::vector<float>> ATbATA = {std::vector<float>{(MultiMat(GenerateTranspose(A), b) / MultiMat(GenerateTranspose(A), A))}};

			newVec = (std::vector<std::vector<float>>{newVec}  - MultiMat(ATbATA, A))[0];
		}
		
		orthoVectors.push_back(newVec);
	}

	PrintMatrix(orthoVectors);
}

int main()
{
	//std::vector<std::vector<float>> matrix {std::vector<float>{2.f,1.f,-1.f},std::vector<float>{-3.f,-1.f,2.f},std::vector<float>{-2.f,1.f,2.f}};
	std::vector<std::vector<float>> matrix {std::vector<float>{1.f, 2.f,2.f,2.f},std::vector<float>{2.f,4.f,6.f, 8.f},std::vector<float>{3.f,6.f,8.f, 10.f}};
	std::cout << "Hello World!" << std::endl;
	//std::vector<std::vector<float>> matrixA {std::vector<float>{1,1},std::vector<float>{1,2},std::vector<float>{1,3}};
	//std::vector<std::vector<float>> matrixB {std::vector<float>{1},std::vector<float>{2},std::vector<float>{2}};
	std::vector<std::vector<float>> matrixA {std::vector<float>{1,1},std::vector<float>{1,0},std::vector<float>{1,2}};
	//LinearCombination(matrix);
	// GaussEliminationPartialPivot(matrix);
	// matrix.clear();
	// matrix.push_back(std::vector<float>{2.f,1.f,-1.f});
	// matrix.push_back(std::vector<float>{-3.f,-1.f,2.f});
	// matrix.push_back(std::vector<float>{-2.f,1.f,2.f});
	// GaussEliminationCompletePivot(matrix);
	//EchelonForm(matrix);
	//matrix = GenerateTranspose(matrixA);

	//LeastSquares(matrixA,matrixB);

	GramSchmidt(matrixA);
	
    return 0;
    
}
