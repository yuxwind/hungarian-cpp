/********************************************************************
 ********************************************************************
 ** C++ class implementation of the Hungarian algorithm by David Schwarz, 2012
 **
 **
 ** C++ Hungarian Class derived from libhungarian by Cyrill Stachniss, 2004
 **
 **
 ** Solving the Minimum Assignment Problem using the 
 ** Hungarian Method.
 **
 ** ** This file may be freely copied and distributed! **
 **
 **
 ** This file is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied 
 ** warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 ** PURPOSE.  
 **
 ********************************************************************
 ********************************************************************/

#include <iostream>
#include <vector>

#ifndef HUNGARIAN_H
#define HUNGARIAN_H

typedef enum {
	HUNGARIAN_MODE_MINIMIZE_COST,
	HUNGARIAN_MODE_MAXIMIZE_UTIL,
} MODE;

typedef enum {
	HUNGARIAN_NOT_ASSIGNED,
	HUNGARIAN_ASSIGNED,
} ASSIGN;


class Hungarian
{

public:
	/** This method initialize the hungarian_problem structure and init 
	 *  the  cost matrices (missing lines or columns are filled with 0).
	 *  It returns the size of the quadratic(!) assignment matrix. **/

	Hungarian();
	Hungarian(std::vector< std::vector<int> >, int, int, MODE);

	int init( std::vector< std::vector<int> > cost_matrix, 
			   int rows, 
			   int cols, 
			   MODE mode);

	/** This method computes the optimal assignment. **/
	bool solve();

	/** Print the computed optimal assignment. **/
	void print_assignment();

	/** Print the cost matrix. **/
	void print_cost();

	/** Print cost matrix and assignment matrix. **/
	void print_status();

protected:
	bool check_solution(std::vector<int> row_dec, std::vector<int> col_inc, std::vector<int> col_mate);
	bool assign_solution(std::vector<int> row_dec, std::vector<int> col_inc, std::vector<int> col_mate);

private:

	int m_totalcost;
	int m_rows;
	int m_cols;
	std::vector< std::vector<int>> m_cost;
	std::vector< std::vector<int>> m_assignment;   

};

#endif



