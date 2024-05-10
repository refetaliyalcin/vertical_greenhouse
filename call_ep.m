function [result] = call_ep(idf_name,weather_name)

code=[
"import numpy as np";
"import math";
"import eppy";
"import sys";
"# pathnameto_eppy = 'c:/eppy'";
"pathnameto_eppy = '../'";
"sys.path.append(pathnameto_eppy)";
"from eppy import modeleditor"
"from eppy.modeleditor import IDF";
"iddfile = 'Energy+.idd'";
['fname1 = ''',idf_name,'.idf'''] %put file name here
"IDF.setiddname(iddfile)";
"idf1 = IDF(fname1)";
"idf1.printidf()";
['epwfile = ''',weather_name,'.epw''']
"idf = IDF(fname1, epwfile)";
['idf.run(output_directory=''',idf_name,''')']
"from eppy.results import readhtml # the eppy module with functions to read the html";
['fname = ''',idf_name,'/eplustbl.htm'''];
"filehandle = open(fname, 'r').read()";
"htables = readhtml.titletable(filehandle) # reads the tables with their titles";
];
result=pyrun(code,["htables"])
end