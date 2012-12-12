/********************************************************************
 ********************************************************************
 ** C++ and STL class implementation of the Hungarian algorithm by David Schwarz, 2012
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


#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "hungarian.h"

#define verbose (0)

#define hungarian_test_alloc(X) do {if ((void *)(X) == NULL) fprintf(stderr, "Out of memory in %s, (%s, line %d).\n", __FUNCTION__, __FILE__, __LINE__); } while (0)

Hungarian::Hungarian()
{

}

Hungarian::Hungarian(std::vector< std::vector<int> > cost_matrix, int rows, int cols, MODE mode)
{
  int i,j, org_cols, org_rows;
  int max_cost;
  max_cost = 0;
  
  org_cols = cols;
  org_rows = rows;

  // is the number of cols  not equal to number of rows ? 
  // if yes, expand with 0-cols / 0-cols
  rows = std::max(cols, rows);
  cols = rows;
  
  m_rows = rows;
  m_cols = cols;

  m_cost.resize(rows, std::vector<int>(cols,0));
  m_assignment.resize(rows, std::vector<int>(cols,0));

  for(i=0; i<m_rows; i++) 
  {
	  for(j=0; j<m_cols; j++) 
	  {
		  m_cost[i][j] =  (i < org_rows && j < org_cols) ? cost_matrix[i][j] : 0;
		  m_assignment[i][j] = 0;

		  if (max_cost < m_cost[i][j])
		  {
			  max_cost = m_cost[i][j];
		  }
	  }
  }


  if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL) 
  {
	  for(i=0; i<m_rows; i++) 
	  {
		  for(j=0; j<m_cols; j++) 
		  {
			  m_cost[i][j] =  max_cost - m_cost[i][j];
		  }
	  }
  }
  else if (mode == HUNGARIAN_MODE_MINIMIZE_COST) {
    // nothing to do
  }
  else 
    fprintf(stderr,"%s: unknown mode. Mode was set to HUNGARIAN_MODE_MINIMIZE_COST !\n", __FUNCTION__);
}


void hungarian_print_matrix(std::vector< std::vector<int> > C, int rows, int cols) 
{
	int i,j;
	fprintf(stderr , "\n");
	for(i=0; i<rows; i++) 
	{
		fprintf(stderr, " [");
		for(j=0; j<cols; j++) 
		{
		fprintf(stderr, "%5d ",C[i][j]);
		}
		fprintf(stderr, "]\n");
	}
	fprintf(stderr, "\n");
}

void Hungarian::print_assignment() {
  hungarian_print_matrix(m_assignment, m_rows, m_cols) ;
}

void Hungarian::print_cost() {
  hungarian_print_matrix(m_cost, m_rows, m_cols) ;
}

void Hungarian::print_status() {
  
  fprintf(stderr,"cost:\n");
  print_cost();

  fprintf(stderr,"assignment:\n");
  print_assignment();
  
}

int Hungarian::init(std::vector< std::vector<int> > cost_matrix, int rows, int cols, MODE mode) 
{

  int i,j, org_cols, org_rows;
  int max_cost;
  max_cost = 0;
  
  org_cols = cols;
  org_rows = rows;

  // is the number of cols  not equal to number of rows ? 
  // if yes, expand with 0-cols / 0-cols
  rows = std::max(cols, rows);
  cols = rows;
  
  m_rows = rows;
  m_cols = cols;

  m_cost.resize(rows, std::vector<int>(cols,0));
  m_assignment.resize(rows, std::vector<int>(cols,0));

  for(i=0; i<m_rows; i++) {
    for(j=0; j<m_cols; j++) {
      m_cost[i][j] =  (i < org_rows && j < org_cols) ? cost_matrix[i][j] : 0;
      m_assignment[i][j] = 0;

      if (max_cost < m_cost[i][j])
	   max_cost = m_cost[i][j];
    }
  }


  if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL) {
    for(i=0; i<m_rows; i++) {
      for(j=0; j<m_cols; j++) {
	m_cost[i][j] =  max_cost - m_cost[i][j];
      }
    }
  }
  else if (mode == HUNGARIAN_MODE_MINIMIZE_COST) {
    // nothing to do
  }
  else 
    fprintf(stderr,"%s: unknown mode. Mode was set to HUNGARIAN_MODE_MINIMIZE_COST !\n", __FUNCTION__);
  
  return rows;
}

bool Hungarian::check_solution(std::vector<int> row_dec, std::vector<int> col_inc, std::vector<int> col_mate)
{
	int k, l, m, n;

	m = m_rows;
	n = m_cols;
  // Begin doublecheck the solution 23
	for (k=0;k<m;k++)
		for (l=0;l<n;l++)
			if (m_cost[k][l]<row_dec[k]-col_inc[l])
			return false;

	for (k=0;k<m;k++)
	{
		l=col_mate[k];
		if (l<0 || m_cost[k][l]!=row_dec[k]-col_inc[l])
			return false;
	}
	k=0;
	
	for (l=0;l<n;l++)
	{
		if (col_inc[l])
		{
			k++;
		}
	}
	if (k>m)
	{
		return false;
	}

	//everything checks out, then
	return true;
  // End doublecheck the solution 23
}
bool Hungarian::assign_solution(std::vector<int> row_dec, std::vector<int> col_inc, std::vector<int> col_mate)
{
	  // End Hungarian algorithm 18
	int i, j, k, l, m, n;

	m = m_rows;
	n = m_cols;

	for (i=0;i<m;++i)
	{
		m_assignment[i][col_mate[i]]=HUNGARIAN_ASSIGNED;
		/*TRACE("%d - %d\n", i, col_mate[i]);*/
	}
	for (k=0;k<m;++k)
	{
		for (l=0;l<n;++l)
		{
		/*TRACE("%d ",m_cost[k][l]-row_dec[k]+col_inc[l]);*/
			m_cost[k][l]=m_cost[k][l]-row_dec[k]+col_inc[l];
		}
		/*TRACE("\n");*/
	}
	for (i=0;i<m;i++)
	{
		m_totalcost+=row_dec[i];
	}
	for (i=0;i<n;i++)
	{
		m_totalcost-=col_inc[i];
	}
	if (verbose)
		fprintf(stderr, "Cost is %d\n",m_totalcost);

	return true;

}

bool Hungarian::solve()
{
	int i, j, m, n, k, l, s, t, q, unmatched, cost;

	m = m_rows;
	n = m_cols;

	int INF = std::numeric_limits<int>::max();

	std::vector<int> col_mate(m);
	std::vector<int> row_mate(n);
	std::vector<int> unchosen_row(m);
	std::vector<int> parent_row(n);
	std::vector<int> row_dec(m);
	std::vector<int> col_inc(n);
	std::vector<int> slack_row(m);
	std::vector<int> slack(n);

	cost=0;

	for (i=0;i<m_rows;i++) 
	{
		col_mate[i]=0;
		unchosen_row[i]=0;
		row_dec[i]=0;
		slack_row[i]=0;
	}

	for (j=0;j<m_cols;j++) 
	{
		row_mate[j]=0;
		parent_row[j] = 0;
		col_inc[j]=0;
		slack[j]=0;
	}

	for (i=0;i<m_rows;++i)
	{
		for (j=0;j<m_cols;++j)
		{
			m_assignment[i][j]=HUNGARIAN_NOT_ASSIGNED;
		}
	}

  // Begin subtract column minima in order to start with lots of zeroes 12
	if (verbose)
	{
		fprintf(stderr, "Using heuristic\n");
	}

	for (l=0;l<n;l++)
	{
		s=m_cost[0][l];

		for (k=1;k<m;k++) 
		{
			if (m_cost[k][l] < s)
			{
				s=m_cost[k][l];
			}
			cost += s;
		}

		if (s!=0)
		{
			for (k=0;k<m;k++)
			{
				m_cost[k][l]-=s;
			}
		}
	}
  // End subtract column minima in order to start with lots of zeroes 12

  // Begin initial state 16
	t=0;
	for (l=0;l<n;l++)
	{
		row_mate[l]= -1;
		parent_row[l]= -1;
		col_inc[l]=0;
		slack[l]=INF;
	}
	for (k=0;k<m;k++)
	{
		bool row_done = false;
		s=m_cost[k][0];

		for (l=1;l<n;l++)
		{
			if (m_cost[k][l]<s)
			{
				s=m_cost[k][l];
			}
			row_dec[k]=s;
		}
		for (l=0;l<n;l++)
		{
			if (s==m_cost[k][l] && row_mate[l]<0)
				{
					col_mate[k]=l;
					row_mate[l]=k;
					if (verbose)
					{
						fprintf(stderr, "matching col %d==row %d\n",l,k);
					}
					row_done = true;
					break;
				}
		}

		if(!row_done)
		{
			col_mate[k]= -1;

			if (verbose)
			{
				fprintf(stderr, "node %d: unmatched row %d\n",t,k);
			}
		
			unchosen_row[t++]=k;
		}

	}
  // End initial state 16

	bool checked = false;
  // Begin Hungarian algorithm 18
	if (t==0)
	{
		checked = check_solution(row_dec, col_inc, col_mate);
		if (checked)
		{
			//finish assignment, wrap up and done.
			bool assign = assign_solution(row_dec, col_inc, col_mate);
			return true;
		}
		else
		{
			if(verbose)
			{
				fprintf(stderr, "Could not solve. Error.\n");
			}
			return false;
		}
	}
	unmatched=t;

	while (1)
	{
		if (verbose)
		{
			fprintf(stderr, "Matched %d rows.\n",m-t);
		}
		q=0;
		bool try_matching;
		while (1)
		{
			while (q<t)
			{
			// Begin explore node q of the forest 19
				
				k=unchosen_row[q];
				s=row_dec[k];
				for (l=0;l<n;l++)
				{
					if (slack[l])
					{
						int del;
						del=m_cost[k][l]-s+col_inc[l];
						if (del<slack[l])
						{
							if (del==0)
							{
								if (row_mate[l]<0)
								{
									goto breakthru;
								}
								slack[l]=0;
								parent_row[l]=k;
								if (verbose)
									fprintf(stderr, "node %d: row %d==col %d--row %d\n",
									t,row_mate[l],l,k);
									unchosen_row[t++]=row_mate[l];
							}
							else
							{
								slack[l]=del;
								slack_row[l]=k;
							}
						}
					}
				}
			// End explore node q of the forest 19
				q++;	
			}
 
	  // Begin introduce a new zero into the matrix 21
		s=INF;
		for (l=0;l<n;l++)
		{
			if (slack[l] && slack[l]<s)
			{
				s=slack[l];
			}
		}
		for (q=0;q<t;q++)
		{
			row_dec[unchosen_row[q]]+=s;
		}
		for (l=0;l<n;l++)
		{
			if (slack[l])
			{
				slack[l]-=s;
				if (slack[l]==0)
				{
					// Begin look at a new zero 22
					k=slack_row[l];
					if (verbose)
					{
						fprintf(stderr, 
						"Decreasing uncovered elements by %d produces zero at [%d,%d]\n",
						s,k,l);
					}
					if (row_mate[l]<0)
					{
						for (j=l+1;j<n;j++)
							if (slack[j]==0)
							{
								col_inc[j]+=s;
							}

						goto breakthru;
						try_matching = true;
						break;
					}
					else
					{
						parent_row[l]=k;
						if (verbose)
						{
							fprintf(stderr, "node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
							unchosen_row[t++]=row_mate[l];
						}
					}
		// End look at a new zero 22
				}
			}
			else
			{
				col_inc[l]+=s;
			}
		}
	// End introduce a new zero into the matrix 21
		if(try_matching)
		{
			break;
		}


	}

    breakthru:
      // Begin update the matching 20
		if (verbose)
		{
			fprintf(stderr, "Breakthrough at node %d of %d!\n",q,t);
		}
		while (1)
		{
			j=col_mate[k];
			col_mate[k]=l;
			row_mate[l]=k;
			if (verbose)
			{
				fprintf(stderr, "rematching col %d==row %d\n",l,k);
			}
			if (j<0)
			{
				break;
			}
			k=parent_row[j];
			l=j;
		}
		// End update the matching 20
		if (--unmatched==0)
		{
			checked = check_solution(row_dec, col_inc, col_mate);
			if (checked)
			{
				//finish assignment, wrap up and done.
				bool assign = assign_solution(row_dec, col_inc, col_mate);
				return true;
			}
			else
			{
				if(verbose)
				{
					fprintf(stderr, "Could not solve. Error.\n");
				}
				return false;
			}
		}
		
		// Begin get ready for another stage 17
			t=0;
			for (l=0;l<n;l++)
			{
				parent_row[l]= -1;
				slack[l]=INF;
			}
			for (k=0;k<m;k++)
			{
				if (col_mate[k]<0)
				{
				if (verbose)
					fprintf(stderr, "node %d: unmatched row %d\n",t,k);
				unchosen_row[t++]=k;
				}
			}
		// End get ready for another stage 17
	}// back to while loop


}


