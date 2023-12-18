#include <iostream>
#include <vector>
#include <assert.h>
using namespace std;

std::vector<int> vec_matrix_mult(const vector<int>& vec, const vector<vector<int>>& matrix, const bool vec_left) {
  size_t rows = matrix.size();
  size_t cols = matrix[0].size();
  std::vector<int> result;  // 初始化结果向量为零向量
  if (vec_left) {
    // 向量在矩阵左侧

    assert(vec.size() == matrix.size() && vec.size() > (size_t)0);

    result.resize(matrix[0].size());
    for (size_t j = 0; j < rows; ++j) {
      for (size_t i = 0; i < cols; ++i) {
        result[i] = result[i] + vec[j] * matrix[j][i];
      }
    }
  } else {
    // 向量在矩阵右侧

    assert(matrix.size() > (size_t)0 && vec.size() == matrix[0].size());

    result.resize(matrix.size());
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        result[i] = result[i] + matrix[i][j] * vec[j];
      }
    }
  }

  return result;
}

int main()
{
	std::vector<int> vec1(2,1);
	std::vector<int> vec2(3,1);
	std::vector<std::vector<int>> matrix = std::vector<std::vector<int>>(2, std::vector<int>(3,1));

   	std::vector<int> mult1 = vec_matrix_mult(vec1, matrix, true);
	std::vector<int> mult2 = vec_matrix_mult(vec2, matrix, false);
	printf("%d, %d, %d\n", mult1[0], mult1[1], mult1[2]);
	printf("%d, %d\n", mult2[0], mult2[1]);
	return 0;
}