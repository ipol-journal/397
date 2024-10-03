/*
  Copyright (c) 2022 Samy Blusseau <samy.blusseau@minesparis.psl.eu>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public
  License (GNU AGPL) along with this program. If not, see
  <http://www.gnu.org/licenses/>.
*/

#ifndef MYSPARSEMATRIX_H
#define MYSPARSEMATRIX_H

#include <cassert>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>



template< typename Type >
class MySparseMatrix{
  unsigned int n_columns;
  std::vector< std::map< unsigned int , Type > > columns;
 public:
  MySparseMatrix(  ) {
    n_columns = 0;
    columns.clear();
    columns.resize(n_columns);
  }
  MySparseMatrix( unsigned int i_n_columns )  {
    n_columns = i_n_columns;
    columns.clear();
    columns.resize(n_columns);
  }
  unsigned int nColumns() const { return n_columns; }
  const std::map< unsigned int , Type > & get_column(unsigned int col) const
    {
      return columns[col];
    }
  Type & operator() (unsigned int row , unsigned int col) { return columns[col][row]; }
  Type operator() (unsigned int row , unsigned int col) const  {
    typename std::map< unsigned int , Type >::const_iterator it = columns[col].find(row);
    if( it == columns[col].end() ) return 0.0;
    return it->second;
  }
};

#endif // MYSPARSEMATRIX_H

