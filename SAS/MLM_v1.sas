libname home '/folders/myfolders/home';              						           /* Edit this appropriately*/ 
%let vertexnumber = 79606;  /* This is the vertex with most electrodes (65 electrodes). Other vertices with more 
							than 60 electrodes are 79610 85531 79607 79608 86621 79605 79623 80854 
							80857 80868 80869 80872 84395 84405 85540 86629 79591 79621 */

ods output Estimates=home.mlmOutput; 
proc glimmix data=home.cong_audonset(where=(Vertex=&vertexnumber));
class vertex electrode test_condition trial_number;
model ERSP_value = test_condition /solution;
random intercept / subject = vertex;
random intercept / subject = electrode(vertex);
covtest 'var(vertex) = 0' 0 .;
covtest 'var(electrode(vertex)) = 0' . 0;
estimate 'Aud-onset vs. Prestim' test_condition 1 -1; 
run;

data _null_;
set home.mlmOutput;
if Label="Aud-onset vs. Prestim" then call symputx("pval", Probt);
run;
 
%put pval = &pval;

