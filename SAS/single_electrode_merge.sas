libname home '/folders/myfolders/home';

proc import file = "/folders/myfolders/audonset_aggregatedd.xlsx" 
out=home.audonset_aggregated  dbms=xlsx replace; 
getnames = Yes ;
run;
 
proc import file = "/folders/myfolders/prestim_aggregated.xlsx" 
out=home.prestim_aggregated  dbms=xlsx replace; 
getnames = Yes ;
run;

proc import file = "/folders/myfolders/single_electrode" 
out=home.single_vertex  dbms=xlsx replace; 
getnames = Yes ;
run;

/* many to many inner join between audonset and a single vertex */

proc sql;
create table tempaudonset_singleelectrode as select
A.Vertex, A.Electrode, 
B.Subject, B.Condition, B.Trial_number, B.ERP_value from
home.single_vertex as A inner join 
home.audonset_aggregated as B on
A.Electrode = B.MergeKeyAudonset;
quit;

/* many to many inner join between prestim and a single vertex */
proc sql;
create table tempprestim_singleelectrode as select
A.Vertex, A.Electrode, 
B.Subject, B.Condition, B.Trial_number, B.ERP_value from
home.single_vertex as A inner join 
home.prestim_aggregated as B on
A.Electrode = B.MergeKeyPrestim;
quit;

data home.tempsingleelectrodedata;
set tempaudonset_singleelectrode tempprestim_singleelectrode;
run;

proc export data = home.tempsingleelectrodedata
outfile = "/folders/myfolders/single_electrode_merged_data.csv" replace ;
run;
