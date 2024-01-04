#ifndef DATAFRAME_H
#define DATAFRAME_H

#include <cstddef>
#include <iostream>
#include <tuple>
#include <vector>


namespace DES 
{

template <typename T>
class DataFrame
{

private:
    size_t _rows;
    size_t _cols;
    std::vector<std::vector<T>> entries;

public:
    DataFrame<T>                (size_t rows, size_t cols);

    std::vector<T>  getRow      (size_t row);
    std::vector<T>  getCol      (size_t col);
    std::vector<T>  operator[]  (size_t row);

    void            addRow      (std::vector<T>& row);
    size_t          getNumRows  ();
    size_t          getNumCols  ();


    friend std::ostream& operator<<(std::ostream& os, DataFrame& df)
    {
        // inline for template reasons
        for (size_t i = 0; i < df.getNumRows(); i++) {
            std::vector<T> row = df.getRow(i);
            for (size_t j = 0; j < df.getNumCols(); j++) {
                os << row.at(j) << '\t';
            }
            os << '\n';
        }
        return os;
    }
    
};



template <typename T>
DataFrame<T>::DataFrame(size_t rows, size_t cols)
{
    this->_rows = rows;
    this->_cols = cols;

    this->entries = std::vector<std::vector<T>>(rows);
    for (int i = 0; i < rows; i++) entries[i] = std::vector<T>(cols, 0);
}


template <typename T>
std::vector<T> DataFrame<T>::getRow(size_t row)
{
    if (row >= _rows)
        throw std::out_of_range("Row is out of range");
    
    return entries.at(row);
}


template <typename T>
std::vector<T> DataFrame<T>::getCol(size_t col) 
{
    if (col >= _cols) 
        throw std::out_of_range("Column is out of range");

    std::vector<T> column;
    for (int i = 0; i < _rows; i++)
        column.push_back(entries.at(i).at(col));
    
    return column;
}


template <typename T>
std::vector<T> DataFrame<T>::operator[] (size_t row) 
{
    return getRow(row);    // return a copy
}


template <typename T>
void DataFrame<T>::addRow(std::vector<T>& row)
{
    if (row.size() != _cols)
        throw std::runtime_error("Row is too big");

    _rows++;
    entries.push_back(row);
}


template <typename T>
size_t DataFrame<T>::getNumRows() 
{
    return _rows;
}


template <typename T>
size_t DataFrame<T>::getNumCols() 
{
    return _cols;
}


}

#endif