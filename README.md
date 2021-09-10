# STRmix 2.7 DB Search Reader
STRmix 2.7 database search report reader to CSV

## General Description
This program parses the *results.xml* files that are part of the normal output of a successful STRmix database search run.

## Running the Program
`sm2_7_db_search_reader.py`

Run the script in the experiment folder.
The script will go through each subfolder and 
identify the results.xml files and read them.
Each read file will become a new line in the resulting CSV file for analysis.




The *results.xml* file contains most of the information found in
the report. This information is slightly different as numbers are
not rounded as they are in the report.

The program first looks inside a folder you set with the
variable 'DIRECTORY'. The program looks for run folders in the
'Directory' folder. If there are 'results.xml' files it will attempt
to parse them.

Internally, the program builds a dataframe or table with all the data
from the 'results.xml' files that it finds. After exhausting its
search of the 'DIRECTORY' folder it will stop and save a CSV file
in the same folder the program is run from.