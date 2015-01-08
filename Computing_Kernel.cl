#include "Settings_and_Parameters.h"

////////////////////////////////////////////////////////////////////////////////
// Laplacation operator definition, to calculate diffusive fluxes
////////////////////////////////////////////////////////////////////////////////

float LaplacianXY(__global float* C, int row, int column)
{
	float retval;
	float dx = dX;
	float dy = dY;
	
	int current = row * Grid_Width + column;
	int left    = row * Grid_Width + column-1;
	int right   = row * Grid_Width + column+1;
	int top     = (row-1) * Grid_Width + column;
	int bottom  = (row+1) * Grid_Width + column;
	
	retval = ( C[left] + C[right]  - 2 * C[current] )/dx/dx +
    ( C[top]  + C[bottom] - 2 * C[current] )/dy/dy;
    
	return retval;
}

////////////////////////////////////////////////////////////////////////////////
// Gradient operators definition, to calculate advective fluxes
////////////////////////////////////////////////////////////////////////////////

float GradientX(__global float* C, int row, int column)
{
	float retval;
	float dx = dX;
	
	int left=row * Grid_Width + column-1;
	int right=row * Grid_Width + column+1;
    
	retval = - ( C[left] - C[right] )/dx/2.0 ;
    
	return retval;
}

float GradientY(__global float* C, int row, int column)
{
	float retval;
	float dy =dY;
	
	int top=(row-1) * Grid_Width + column;
	int bottom=(row+1) * Grid_Width + column;
    
	retval = - ( C[top] - C[bottom] )/dy/2.0 ;
    
	return retval;
}


////////////////////////////////////////////////////////////////////////////////
// Simulation kernel
////////////////////////////////////////////////////////////////////////////////

__kernel void
Mussels_PDE_Kernel_Phase1 (__global float* C, __global float* J1X, __global float* J1Y, __global float* J2)
{
    
    float beta = Beta;
    float gm;
    
    size_t current = get_global_id(0);
	int row		= floor((float)current/(float)Grid_Width);
	int column	= current%Grid_Width; 
    
    if(row > 0 && row < Grid_Height-1 && column > 0 && column < Grid_Width-1)
	{
        //Now calculating terms for the P Matrix
        gm = ( pow(C[current],2) - beta*C[current] + 1 ) * ( 3*pow(C[current],2) - 2*beta*C[current] + 1 );
        //gm = ( C[current]*C[current] - beta*C[current] + 1 ) * ( 3*C[current]*C[current] - 2*beta*C[current] + 1 );
        J1X[current] = gm * GradientX(C, row, column);
        J1Y[current] = gm * GradientY(C, row, column);
        J2[current]  = LaplacianXY(C, row, column);
	}
    else
	{
        J1X[current] = 0;
        J1Y[current] = 0;
        J2[current]  = 0;
	}
    
    barrier(CLK_GLOBAL_MEM_FENCE);
    
}

__kernel void
Mussels_PDE_Kernel_Phase2 (__global float* C, __global float* J1X, __global float* J1Y, __global float* J2)
{
    
    float J3X;
    float J3Y;
    float J4;
    
    float d = D;
    float k = k1;
    
    size_t current = get_global_id(0);
	
	int row		= floor((float)current/(float)Grid_Width);
	int column	= current%Grid_Width;
    
    if(row > 1 && row < Grid_Height-2 && column > 1 && column < Grid_Width-2)
	{
        //Now calculating terms for the P Matrix
        J3X = d * GradientX(J1X, row, column);
        J3Y = d * GradientY(J1Y, row, column);
        J4  = d * k * LaplacianXY(J2, row, column);
        C[current]=C[current]+(J3X+J3Y-J4)*dT;
        
	}
    
    // barrier(CLK_GLOBAL_MEM_FENCE);
    
    // HANDLE Boundaries
    
    if(row<=1)
	{
        C[current]=C[(row + Grid_Height - 4) * Grid_Width + column];
	}
    else if(row>=Grid_Height-2)
	{
        C[current]=C[(row - Grid_Height + 4) * Grid_Width + column];
	}
    else if(column<=1)
	{
        C[current]=C[row * Grid_Width + column + Grid_Width - 4];
	}
    else if(column>=Grid_Width-2)
	{
        C[current]=C[row * Grid_Width + column - Grid_Width + 4];
	}
    
    // barrier(CLK_GLOBAL_MEM_FENCE);
    
} // End Mussels_PDE_Kernel
