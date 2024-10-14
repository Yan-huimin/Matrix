#ifndef YHM_MATRIX_H
#define YHM_MATRIX_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

class yhm_Matrix
{
private:
    //行数与列数
    size_t yrow;
    size_t ycol;
    //存储矩阵数值
    std::vector<std::vector<double>> yMaxtrix;
public:
    //构造函数
    yhm_Matrix(size_t rows, size_t cols);

    yhm_Matrix();

    //析构函数
    ~yhm_Matrix();

    void setElement(size_t x, size_t y, double value);

    double getElement(size_t x, size_t y) const;

    size_t getRow() const;

    size_t getCol() const;

    //输出当前的矩阵
    void yhm_Print_Matrix();

    //求解该矩阵的转置矩阵
    yhm_Matrix yhm_Transpose_Matrix();

    //矩阵行列式求解
    double yhm_Cal_Det();

    // 返回以x列为中心点的代数余子式，返回类型为矩阵
    yhm_Matrix getAlgebraic_cofactor(size_t x, size_t y);
};


#endif