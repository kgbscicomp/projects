libname home '/folders/myfolders/home';

ods excel file="/folders/myfolders/vertex_frequency.xlsx";
proc freq data= home.all_electrodes_vertices;
   tables Vertex Electrode/ out=FreqCount outexpect sparse;
   title 'Vertex and electrode frequeny table';
run;
ods excel close;

/*
proc export data = FreqCount
dbms = xlsx
outfile = "/folders/myfolders/vertex_frequency.xlsx" replace ;
run;
*/