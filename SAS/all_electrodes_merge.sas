libname home '/folders/myfolders/home';
%let stimdata = "/folders/myfolders/cong/audonset/audonset_aggregated.xlsx";    /*Edit this path according to the file being processed*/
%let prestimdata = "/folders/myfolders/cong/prestim/prestim_aggregated.xlsx";   /*Edit this path according to the file being processed*/
%let verticesdata = "/folders/myfolders/elecsvertices_longformatcheck.xlsx" ;  /* Do Not Edit this */
%let outfilexlsx = "/folders/myfolders/cong_audonset_one_vertex_datasubset.xlsx"; /* Edit the filename as required*/
%let outfile = home.cong_audonset;              						           /* Edit this appropriately*/ 
%let vertexnumber = 79606;  /* This is the vertex with most electrodes (65 electrodes). Other vertices with more 
							than 60 electrodes are 79610 85531 79607 79608 86621 79605 79623 80854 
							80857 80868 80869 80872 84395 84405 85540 86629 79591 79621 */
							
/*Reads in data in long format for three files
  1) stimuli condition 
  2) prestimuli condition
  3) Electrode names in each of the vertices in the MNI space
  The first two files have the following fields: 
  Subject_ID	Test_Condition	Electrode_name	MergeKey  Trial_number	ERSP_values
  -Test_condition refers to: Either prestim or one of the stimuli conditions (audonset/faceonset/facemove)
  -MergeKey is formed by combining SubjectID&Electrode_name 
  
  We need to perforn a many-to-many merge between the stimuli/prestimuli conditions and the electrodes in 
  each of the vertices such that the Subject_ID+Electrode_name is used as a merge key to find a correspondence
  between the electrodes in each of the vertices in the primary auditory cortex and their corresponding 
  ERSP values. 
  
  Change the folder & file names as follows: 
  For instance, if a visual only stimuli's audonset condition is being processed, the folder would be: 
  "........../vis/audonset/audonset_aggregated.xlsx"
  If an auditory only stimuli's faceonset condition is being processed, the folder would be: 
  "........../aud/faceonset/faceonset_aggregated.xlsx"
  If a congruent condition's facemove condition is being processed, the folder would be: 
  "........../cong/facemove/facemove_aggregated.xlsx"
  If a congruent condition's facemove condition is being processed, the folder would be: 
  "........../cong/prestim/prestim_aggregated.xlsx"

  */
proc import file = &stimdata   /*Change folder&file name*/
out=home.stimdata_aggregated  dbms=xlsx replace; 
getnames = Yes ;
run;
 

proc import file = &prestimdata     /*Change folder name*/
out=home.prestimdata_aggregated  dbms=xlsx replace; 
getnames = Yes ;
run;

proc import file = &verticesdata
out=home.all_electrodes_vertices  dbms=xlsx replace; 
getnames = Yes ;
run;

/* many to many inner join between electrodes in a stimuli condition
   and all available vertices in the primary auditory cortex*/

proc sql;
create table tempstimdata_allelecs as select
A.Vertex, A.Electrode, 
B.Subject_ID, B.Test_Condition, B.Trial_number, B.ERSP_value from
home.all_electrodes_vertices as A inner join 
home.stimdata_aggregated as B on
A.Electrode = B.MergeKey;
quit;

/* many to many inner join between the prestimuli condition
   and all available vertices in the primary auditory cortex */
proc sql;
create table tempprestim_allelecs as select
A.Vertex, A.Electrode, 
B.Subject_ID, B.Test_Condition, B.Trial_number, B.ERSP_value from
home.all_electrodes_vertices as A inner join 
home.prestimdata_aggregated as B on
A.Electrode = B.MergeKey;
quit;

/*Appending stimuli and prestimuli conditions and exporting data as a SAS file
data _null_;					   
set tempstimdata_allelecs tempprestim_allelecs;
file outfile;
put Vertex Subject Condition Electrode Trial_number ERP_value;
run;
*/

data &outfile;
set tempstimdata_allelecs tempprestim_allelecs;
run;


/*Exporting just a subset of the data as csv*/
proc export 
data = home.cong_audonset (where=(Vertex=&vertexnumber))		
dbms = xlsx
outfile = &outfilexlsx replace ;
run;


