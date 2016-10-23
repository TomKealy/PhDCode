#include <mex.h>
#include <iostream>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>

//   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

struct Group
{
	std::vector<int> Data;

	bool operator < (const Group & Group_) const
	{   
		const int Size = std::min(Group_.Data.size(),Data.size());

		for(int i = 0 ; i < Size ; i++)
			if(Data[i] != Group_.Data[i])
				return Data[i] < Group_.Data[i];

		return Data.size() < Group_.Data.size(); 
	}

		
	
	Group operator + (const Group & Group_) const
	{   
		    Group Tmp;
		    
		   
		    set_intersection ( Data.begin(), Data.end(), 
			               Group_.Data.begin(), Group_.Data.end(), 
			               inserter(Tmp.Data, Tmp.Data.begin()) );
			
			
		    return Tmp;


	}
    
};




void mexFunction(int nlhs, mxArray * plhs[],int nrhs, const mxArray * prhs[]) 
{
    
	if ( nrhs != 1 )
	{
       		mexErrMsgTxt("One argument needed.");
       		return;
   	}

   	if ( nlhs > 2 )
   	{
       		mexErrMsgTxt("At most 2 outputs needed.");
       		return;
   	}   

   	const mwSize * dimsM  = mxGetDimensions(prhs[0]);
	
   	if ( mxIsLogical(prhs[0]) )
   	{

       	std::set<Group> G;
        std::set<Group> IntersectionClosure;

       	bool* M  = (bool*)mxGetPr(prhs[0]);

       	for (int c = 1; c <= dimsM[1]; c++) 
		{

			Group Tmp;
			
			for (int r = 1 ; r <= dimsM[0] ; r++)
			{
				if (not  *M ) //test if there is a one for the row r and the column c
				{
					Tmp.Data.push_back(r);
                    
				}
                                M++;
			} 

			G.insert(Tmp);

		}	
      
       	std::set<Group>::iterator it_group;
    	std::set<Group>::iterator it_intersectionclosure;

    	for (it_group = G.begin(); it_group != G.end(); it_group++) 
		{

            IntersectionClosure.insert( (*it_group) );

            std::vector<Group> Tmp;

            for (it_intersectionclosure = IntersectionClosure.begin(); it_intersectionclosure != IntersectionClosure.end(); it_intersectionclosure++) 
            {
                Tmp.push_back( (*it_group) + (*it_intersectionclosure) );
            }

            IntersectionClosure.insert(Tmp.begin(),Tmp.end());

    	}
	//===========================================================================================================
	
	std::set<Group> IntersectionClosureWithoutInclusions;
	std::set<Group>::iterator it_a;
    	std::set<Group>::iterator it_b;
	
	bool is_b_included_in_a;
	
	for (it_a = IntersectionClosure.begin(); it_a != IntersectionClosure.end(); it_a++)

	  if ( it_a->Data.size() > 0 ) // We do not consider the empty vector 
	
            {
		is_b_included_in_a = false; 
	    
		for (it_b = IntersectionClosure.begin(); it_b != IntersectionClosure.end(); it_b++)    
		    
		    if (  
			  ( it_b->Data.size() > 0  )                 & // We do not consider the empty  vector
			  ( it_a->Data.size() >  it_b->Data.size() ) & // No comparison with himself
		          std::includes( it_a->Data.begin(), it_a->Data.end(), it_b->Data.begin(), it_b->Data.end() )  
			)
			
			  {
			      is_b_included_in_a = true;
			      break;
			
			  }
		    
		if ( not is_b_included_in_a )
		  IntersectionClosureWithoutInclusions.insert( (*it_a) );

            }
	
	
	
	//===========================================================================================================
       
        const mwSize dimsOutput1[2] = { dimsM[0], IntersectionClosureWithoutInclusions.size() };
	
       	plhs[0] = mxCreateNumericArray(2,dimsOutput1,mxLOGICAL_CLASS,mxREAL);
	
        bool * Output1  = (bool*)mxGetPr(plhs[0]);
	
        //We initialize the matrix with true - we indeed want the matrix of nonzero pattern in the end
        std::fill(Output1,Output1+dimsM[0]*IntersectionClosureWithoutInclusions.size(),false);
	
        for (it_intersectionclosure = IntersectionClosureWithoutInclusions.begin(); it_intersectionclosure != IntersectionClosureWithoutInclusions.end(); it_intersectionclosure++) 
	       {
            
               for (int i=0; i < it_intersectionclosure->Data.size(); i++)
               {
                   *(Output1 + it_intersectionclosure->Data[i] - 1) = true;
               }
	       
               Output1 += dimsM[0];
	       
       	       }
     
       	return;
	
   	}
   	mexErrMsgTxt("Argument must be 'logical'.");


}
