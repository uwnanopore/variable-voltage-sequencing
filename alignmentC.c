#include "mex.h"
#include "math.h"
#include "stdlib.h"

void align( double* score_matrix, long int num_rows, long int num_cols,
        double* step_penalties, double* bad_penalties, 
        double* slip_penalties, double* modes, long int num_modes,
        long int lookback, double* slip_locations,
        int use_periodic_boundaries, double* alignment_trace,
        double* cumulate_score, double* level_is_bad, 
        double* alignment_matrix, int num_step_penalties_per_row );
void selfAlign( double* score_matrix, long int num_rows, long int num_cols,
        double* step_penalties, double SA_penalty, double* alignment_trace,
        double* cumulate_score, double* alignment_matrix, 
        int num_step_penalties_per_row );
long int mod(long int a, long int b);



void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    
    /* Get data from inputs */
    double* score_matrix = mxGetPr(prhs[0]); //the matrix of match scores
    long int num_rows = mxGetM(prhs[0]); // the number of rows i.e. measured levels
    long int num_cols = mxGetN(prhs[0]); // the number of columns i.e. ref levels
    double* step_type_penalties = mxGetPr(prhs[1]); // the array of step penalties
    double* modes = mxGetPr(prhs[2]); // the list of mode directions
    long int num_modes = mxGetN(prhs[2]); // the number of modes
    long int lookback = (long int)mxGetScalar(prhs[3]); // the distance to search for a match
    char alignment_type = (char)mxGetScalar(prhs[4]); // the type of alignment, either "R" for a reference alignment or "S" for self alignment
    double* slip_locations = mxGetPr(prhs[5]); // the places to allow for slips
    int use_periodic_boundaries = (int)mxGetScalar(prhs[6]); //whether to allow an alignment to roll over
    
    /* assign transition penalties */
    double step_penalties[(2*lookback+2)*num_modes*num_rows];
    double bad_penalties[2*num_modes];
    double slip_or_SA_penalties[num_modes];
    
    
        for (int cr = 0; cr < num_rows; cr++) {
            int row_offset_step = cr*num_modes*(2*lookback+2);
            int row_offset_type = cr*num_modes*9;
    for (int cm = 0; cm < num_modes; cm++) {
        
        bad_penalties[2*cm] = step_type_penalties[6*num_modes + cm];
        bad_penalties[2*cm+1] = step_type_penalties[7*num_modes + cm];
        slip_or_SA_penalties[cm] = step_type_penalties[8*num_modes + cm];
        
        
        if (modes[cm] == 1) {
            
            step_penalties[row_offset_step + cm*(2*lookback+2) + lookback] = step_type_penalties[row_offset_type + cm];
            step_penalties[row_offset_step + cm*(2*lookback+2) + lookback + 1] = step_type_penalties[row_offset_type + 5*num_modes + cm];
            
            for (int cs = 0; cs < lookback; cs++) {
                step_penalties[row_offset_step + cm*(2*lookback+2)+cs] = step_type_penalties[row_offset_type + num_modes + cm] + (lookback-cs-1)*step_type_penalties[row_offset_type + 2*num_modes + cm];
                step_penalties[row_offset_step + cm*(2*lookback+2)+lookback+2+cs] = step_type_penalties[row_offset_type + 3*num_modes + cm] + cs*step_type_penalties[row_offset_type + 4*num_modes + cm];
            }
            
        } else if (modes[cm] == -1) {
            
            step_penalties[row_offset_step + cm*(2*lookback+2) + lookback + 1] = step_type_penalties[row_offset_type + cm];
            step_penalties[row_offset_step + cm*(2*lookback+2) + lookback] = step_type_penalties[row_offset_type + 5*num_modes + cm];
            
            for (int cs = 0; cs < lookback; cs++) {
                step_penalties[row_offset_step + cm*(2*lookback+2)+cs] = step_type_penalties[row_offset_type + 3*num_modes + cm] + (lookback-cs-1)*step_type_penalties[row_offset_type + 4*num_modes + cm];
                step_penalties[row_offset_step + cm*(2*lookback+2)+lookback+2+cs] = step_type_penalties[row_offset_type + num_modes + cm] + cs*step_type_penalties[row_offset_type + 2*num_modes + cm];
            }
            
        }
        
    }
    }
    /* check for proper number and format of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("alignment:nrhs","Seven inputs required.");
    }
    if(nlhs!=4) {
        mexErrMsgIdAndTxt("alignment:nlhs","Three outputs required.");
    }
    if(num_modes*num_rows*9 != mxGetN(prhs[1])) {
        mexErrMsgIdAndTxt("alignment:stepCountsAndModes",
                "Number of rows in step counts must equal number of modes.");
    }
    
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)num_rows, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, (mwSize)num_rows, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, (mwSize)num_rows, mxREAL);
    plhs[3] = mxCreateDoubleMatrix( (mwSize)num_rows*num_modes, (mwSize)num_cols, mxREAL );
    
    /* get a pointer to the real data in the output matrices */
    double* alignment_trace = mxGetPr(plhs[0]);
    double* bad_levels = mxGetPr(plhs[1]);
    double* cumulate_score = mxGetPr(plhs[2]);
    double* alignment_matrix = mxGetPr(plhs[3]);
    int num_step_penalties_per_row = num_modes*(2*lookback+2);
    
    if (alignment_type == 'R') {
        //Alignment to reference
        align(
                score_matrix,
                num_rows,
                num_cols,
                step_penalties,
                bad_penalties,
                slip_or_SA_penalties,
                modes,
                num_modes,
                lookback,
                slip_locations,
                use_periodic_boundaries,
                alignment_trace,
                cumulate_score,
                bad_levels,
                alignment_matrix,
                num_step_penalties_per_row
                );
        
    } else if (alignment_type == 'S') {
        //Alignment to self
        selfAlign(
                score_matrix,
                num_rows,
                num_cols,
                step_penalties,
                slip_or_SA_penalties[0],
                alignment_trace,
                cumulate_score,
                alignment_matrix,
                num_step_penalties_per_row
                );
        
    }
    
}

/* proper modulus function; % is actually remainder */
long int mod(long int a, long int b)
{
    long int r = a % b;
    return r < 0 ? r + b : r;
}


/* alignment to reference */
void align(
        double* score_matrix,
        long int num_rows,
        long int num_cols,
        double* step_penalties,
        double* bad_penalties,
        double* slip_penalties,
        double* modes,
        long int num_modes,
        long int lookback,
        double* slip_locations,
        int use_periodic_boundaries,
        double* alignment_trace,
        double* cumulate_score,
        double* level_is_bad,
        double* alignment_matrix,
        int num_step_penalties_per_row
        ) {
    
    long int *traceback_matrix;
    int *bad_level_matrix;
    long int matrix_size = num_rows*num_cols;
    double this_score;
    double last_score;
   
    bad_level_matrix = mxCalloc(num_rows*num_cols*num_modes,sizeof(int));
    traceback_matrix = mxCalloc(num_rows*num_cols*num_modes,sizeof(long int));
    
    /* initialize the alignment matrix with the first row being the same as
     * the score matrix, plus marking bad levels */
    for (long int cc = 0; cc < num_cols; cc++) {
        for (long int cm = 0; cm < num_modes; cm++) {
            
            alignment_matrix[cm*matrix_size + cc*num_rows] = score_matrix[cc*num_rows];
            bad_level_matrix[cm*matrix_size + cc*num_rows] = 0;
            traceback_matrix[cm*matrix_size + cc*num_rows] = -1;
            
            if (score_matrix[cc*num_rows] < bad_penalties[2*cm]) {
                bad_level_matrix[cm*matrix_size+cc*num_rows] = 1;
                alignment_matrix[cm*matrix_size+cc*num_rows] = bad_penalties[2*cm];
            }
        }
    }
    
    
    
    //loop over all rows
    for (long int cr = 1; cr < num_rows; cr++) {
        
        //all columns
        for (long int cc = 0; cc < num_cols; cc++) {
            
            //and all modes
            for (long int cm = 0; cm < num_modes; cm++) {
                
                
                alignment_matrix[cm*matrix_size + cc*num_rows + cr] = -INFINITY;
                
                //check if we are at a mode transition; if so, allow for
                //the jump from the last row of the last mode.
                
                
                
                //if we are at a place where we are allowing a slip
                if (slip_locations[cr]) {
                    
                    //go through all possible transition locations
                    for (long int co = 0; co < num_cols; co++) {
                        
                        //take the slip penalty
                        
                        this_score = alignment_matrix[cm*matrix_size + co*num_rows + cr-1] + slip_penalties[cm];
                        
                        //if it's the best score so far, update the current choice
                        if (this_score > alignment_matrix[cm*matrix_size + cc*num_rows + cr]) {
                            last_score = alignment_matrix[cm*matrix_size + co*num_rows + cr-1];
                            alignment_matrix[cm*matrix_size + cc*num_rows + cr] = this_score;
                            traceback_matrix[cm*matrix_size + cc*num_rows + cr] = cm*matrix_size + co*num_rows + (cr-1);
                        }
                        
                    }
                }
                
                
                //next, we calculate the normal penalties for allowed
                //transitions. This looks different depending on whether we
                //are reading forwards or backwards, and whether we are
                //using periodic boundary conditions.
                for (long int co = cc - lookback - 1 + (modes[cm] == -1); co < cc + lookback + (modes[cm] == -1); co++ ) {
                    
                    long int ix = co;
                    
                    if (use_periodic_boundaries)
                        ix = mod(ix, num_cols);
                    else if (ix >= num_cols || ix < 0)
                        continue;
                    
                    
                    this_score = alignment_matrix[cm*matrix_size+ ix*num_rows + (cr-1)] + step_penalties[cr*num_step_penalties_per_row + (2*lookback+2)*cm + co - cc + lookback + 1 - (modes[cm] == -1)];
                    
                    //if it's the best score so far, update the current choice
                    if (this_score > alignment_matrix[cm*matrix_size + cc*num_rows + cr]) {
                        alignment_matrix[cm*matrix_size + cc*num_rows + cr] = this_score;
                        traceback_matrix[cm*matrix_size + cc*num_rows + cr] = cm*matrix_size  +ix*num_rows + cr-1;
                    }
                    
                    //if we are looking at a mode other than the first, we
                    //need to look at transitions from the previous mode
                    if (cm > 0) {
                        
                        this_score = alignment_matrix[(cm-1)*matrix_size+ ix*num_rows + (cr-1)] + step_penalties[cr*num_step_penalties_per_row + (2*lookback+2)*(cm-1) + co - cc + lookback + 1 - (modes[cm] == -1)];
                        
                        if (this_score > alignment_matrix[cm*matrix_size + cc*num_rows + cr]) {
                            traceback_matrix[cm*matrix_size + cc*num_rows + cr] = (cm-1)*matrix_size + ix*num_rows + (cr-1);
                            alignment_matrix[cm*matrix_size + cc*num_rows + cr] = alignment_matrix[(cm-1)*matrix_size + ix*num_rows + (cr-1)];
                        }
                    }
                    
                    
                    
                }
                
                alignment_matrix[cm*matrix_size + cc*num_rows + cr] += score_matrix[cc*num_rows + cr];
                bad_level_matrix[cm*matrix_size + cc*num_rows + cr] = 0;
                

                    //check the bad level penalty
                if ( alignment_matrix[cm*matrix_size + cc*num_rows + cr] - alignment_matrix[traceback_matrix[cm*matrix_size + cc*num_rows + cr]] 
                                    < bad_penalties[2*cm*cr + bad_level_matrix[traceback_matrix[cm*matrix_size + cc*num_rows + cr]]] ) {
                    bad_level_matrix[cm*matrix_size + cc*num_rows + cr] = 1;
                    traceback_matrix[cm*matrix_size + cc*num_rows + cr] = cm*matrix_size  +cc*num_rows + cr-1;
                    alignment_matrix[cm*matrix_size + cc*num_rows + cr] = alignment_matrix[traceback_matrix[cm*matrix_size + cc*num_rows + cr]] + bad_penalties[2*cm*cr+bad_level_matrix[traceback_matrix[cm*matrix_size + cc*num_rows + cr]]];
                    
                }
            }
            
        }
        
    }
    
    //find the maximum score on the bottom row
    long int current_index;
    double best_score = -INFINITY;
    
    for (long int cm = 0; cm < num_modes; cm++) {
        for (long int cc = 0; cc < num_cols; cc++) {
            
            if (alignment_matrix[cm*matrix_size + cc*num_rows + num_rows-1] > best_score) {
                best_score = alignment_matrix[cm*matrix_size + cc*num_rows + num_rows-1];
                current_index = cm*matrix_size + cc*num_rows + num_rows-1;
            }
            
        }
    }
    
    //perform the traceback
    for (long int cr = num_rows-1; cr >= 0; cr--) {
        long int which_col = (long int) ( (current_index%matrix_size) / num_rows) + 1;
        alignment_trace[cr] = which_col;
        level_is_bad[cr] = bad_level_matrix[current_index];
        cumulate_score[cr] = alignment_matrix[current_index];
        current_index = traceback_matrix[current_index];
    }

    
    mxFree(traceback_matrix);
    mxFree(bad_level_matrix);
    
}








void selfAlign(
        double* score_matrix,
        long int num_rows,
        long int num_cols,
        double* step_penalties,
        double SA_penalty,
        double* alignment_trace,
        double* cumulate_score,
        double* alignment_matrix,
        int num_step_penalties_per_row
        ) {
    
    long int *traceback_matrix;
    double this_score;
    int lookout = (num_step_penalties_per_row - 2)/2;
    
    traceback_matrix = mxCalloc(num_rows*num_cols, sizeof(long int));
    
    alignment_matrix[(num_cols-1)*num_rows] = score_matrix[(num_cols-1)*num_rows];
    traceback_matrix[(num_cols-1)*num_rows] = -1;
    
    //loop over all rows
    for (long int cr = 1; cr < num_rows; cr++) {
        int row_offset = cr*num_step_penalties_per_row;
        //all columns
        for (long int cc = fmax(num_cols-1-cr, 0); cc < num_cols; cc++) {
            
            
            alignment_matrix[cc*num_rows + cr] = -INFINITY;
           
            for (long int co = fmax(num_cols - cr, 0); co < num_cols; co++ ) {
                
                
                this_score = alignment_matrix[co*num_rows + (cr-1)] +  ( (cc == num_cols-1)?(SA_penalty + fmax(step_penalties[row_offset + lookout], step_penalties[row_offset + lookout - 1])):(step_penalties[row_offset + co - cc + num_cols - 1]) );
                
                //if it's the best score so far, update the current choice
                if (this_score > alignment_matrix[cc*num_rows + cr]) {
                    alignment_matrix[cc*num_rows + cr] = this_score;
                    traceback_matrix[cc*num_rows + cr] = co*num_rows + cr-1;
                }
                
            }
            
            alignment_matrix[cc*num_rows + cr] += score_matrix[cc*num_rows + cr];
            
            
        }
        
        
    }
   
    
    //find the maximum score on the bottom row
    long int current_index;
    double best_score = -INFINITY;
    
    for (long int cc = 0; cc < num_cols; cc++) {
        
        if (alignment_matrix[cc*num_rows + num_rows-1] > best_score) {
            best_score = alignment_matrix[cc*num_rows + num_rows-1];
            current_index = cc*num_rows + num_rows-1;
        }
        
    }
    
    
    //perform the traceback
    for (long int cr = num_rows-1; cr >= 0; cr--) {
        long int rounded = (long int) (current_index / num_rows) + 1 ;
        alignment_trace[cr] = rounded;
        cumulate_score[cr] = alignment_matrix[current_index];
        current_index = traceback_matrix[current_index];
    }
    
    
    mxFree(traceback_matrix);
    
}










/* The gateway function */

