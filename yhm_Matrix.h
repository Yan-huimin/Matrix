#ifndef YHM_MATRIX_H
#define YHM_MATRIX_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <stdio.h>

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

    // 初始化为单位矩阵的构造函数
    yhm_Matrix(size_t rows, size_t cols, bool is_Identy_Matrix);

    // 重载构造函数
    yhm_Matrix();

    //析构函数
    ~yhm_Matrix();

    // 为矩阵某一位置的元素赋值
    void setElement(size_t x, size_t y, double value);

    // 获取矩阵某一位置的元素值
    double getElement(size_t x, size_t y) const;

    // 获取矩阵当前的行数
    size_t getRow() const;

    // 获取矩阵当前的列数
    size_t getCol() const;

    //输出当前的矩阵
    void yhm_Print_Matrix();

    //求解该矩阵的转置矩阵
    yhm_Matrix yhm_Transpose_Matrix();

    //矩阵行列式求解
    double yhm_Cal_Det();

    // 返回以x列为中心点的代数余子式，返回类型为矩阵
    yhm_Matrix getAlgebraic_cofactor(size_t x, size_t y);

    // 返回当前矩阵的逆矩阵
    yhm_Matrix getInverse();

    // 求解当前矩阵的余子式矩阵
    yhm_Matrix getConfactor_Matrix();

    // 矩阵加法
    yhm_Matrix yhm_Add(yhm_Matrix b);

    // 矩阵减法
    yhm_Matrix yhm_Sub(yhm_Matrix b);

    // 矩阵乘法
    yhm_Matrix yhm_Mul(yhm_Matrix b);

    // 矩阵的秩
    int yhm_GetRank();

    // 初等行变换
    void yhm_RowTransform(size_t x, size_t y);
};


#endif