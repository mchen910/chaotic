/**
 * @file dataframe.cpp
 * @author Matthew Chen
 * @brief A really simple dataframe with entries all of one type.
 * @version 0.1
 * @date 2024-01-02
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "dataframe.h"


namespace DES
{

template <typename T>
DataFrame<T>::DataFrame(size_t rows, size_t cols)
{
    this->rows = rows;
    this->cols = cols;

    this->entries = std::vector<std::vector<T>>(rows);
    for (int i = 0; i < rows; i++) entries[i] = std::vector<T>(cols, 0);
}


template <typename T>
std::vector<T> DataFrame<T>::getRow(size_t row)
{
    if (row >= rows)
        throw std::out_of_range("Row is out of range");
    
    return std::vector<T>(entries.at(row));
}


/**
 * @brief Get a copy of a particular row in a DataFrame.
 * 
 * @tparam T 
 * @param row 
 * @return std::vector<T> 
 */
template <typename T>
std::vector<T> DataFrame<T>::operator[] (size_t row) 
{
    return getRow(row);    // return a copy
}


template <typename T>
std::vector<T> DataFrame<T>::getCol(size_t col)
{
    if (col >= cols) 
        throw std::out_of_range("Column is out of range");

    std::vector<T> column;
    for (int i = 0; i < rows; i++)
        column.push_back(entries.at(i).at(col));
    
    return column;
}


template <typename T>
void DataFrame<T>::addRow(std::vector<T>& row)
{
    if (row.size() != cols)
        throw std::runtime_error("Row is too big");

    rows++;
    entries.push_back(row);
}



}