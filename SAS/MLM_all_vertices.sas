libname home '/folders/myfolders/home';    	           

%macro mlmvertex;
	class vertex electrode test_condition trial_number;
	model ERSP_value = test_condition /solution;
	random intercept / subject = vertex;
	random intercept / subject = electrode(vertex);
	covtest 'var(vertex) = 0' 0 .;
	covtest 'var(electrode(vertex)) = 0' . 0;
	estimate 'Aud-onset vs. Prestim' test_condition 1 -1; 
%mend mlmvertex;

proc sort data=home.cong_audonset out=home.cong_audonset_sorted;	
   by vertex electrode;
run;

ods output Estimates=home.mlmOutput; 
proc glimmix data=home.cong_audonset_sorted;
   by vertex;
   %mlmvertex;
run;

proc export 
data = home.mlmOutput 		
dbms = xlsx
outfile = '/folders/myfolders/pvaluestest.xlsx' replace ;
run;
