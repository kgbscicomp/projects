libname home '/folders/myfolders/home';

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
  
  Change the folder & file names in the following line: 
  line 34
  For instance, if a visual only stimuli's audonset condition is being processed, the folder would be: 
  "........../vis/audonset_aggregated.xlsx"
  If an auditory only stimuli's faceonset condition is being processed, the folder would be: 
  "........../aud/faceonset_aggregated.xlsx"
  If a congruent condition's facemove condition is being processed, the folder would be: 
  "........../cong/facemove_aggregated.xlsx"
  
  Change folder name in the following line: line 39 such that the folder name corresponds to the condition
  being processed. For instance "............/aud/prestim_aggregated.xlsx", 
  "......./vis/prestime_aggregated.xlsx", etc. 
  
  Change file names in the following lines: Lines 72 and 78. Name files in these lines appropriately 
  as per your convenience. 
  */
proc import file = "/folders/myfolders/aud/audonset_aggregated.xlsx"    /*Change folder&file name*/
out=home.stimdata_aggregated  dbms=xlsx replace; 
getnames = Yes ;
run;
 
proc import file = "/folders/myfolders/aud/prestim_aggregated.xlsx"      /*Change folder name*/
out=home.prestimdata_aggregated  dbms=xlsx replace; 
getnames = Yes ;
run;

proc import file = "/folders/myfolders/elecsvertices_longformatcheck.xlsx" 
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

/*Appending stimuli and prestimuli conditions*/
data home.aud_audonset_data;											   /*Change file name*/
set tempstimdata_allelecs tempprestim_allelecs;
run;

/*Exporting a single vertex of the data as csv: vertex 79606 has the most electrodes. 
More options with 60+ electrodes: 79610 85531 79607 79608 86621 79605 79623 80854 80857 80868 80869 80872
84395 84405 85540 86629 79591 79621 */
proc export data = home.aud_audonset_data( where=(Vertex = 79606)	)     /*Change file name*/
dbms = xlsx
outfile = "/folders/myfolders/single_vertex_data.xlsx" replace ;
run;


