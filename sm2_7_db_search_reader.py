# This program parses the 'results.xml' files that are part of the normal
#       output of a successful STRmix database search run.
#
#       The 'results.xml' file contains most of the information found in
#       the report. This information is slightly different as numbers are
#       not rounded as they are in the report.
#
#       The program first looks inside a folder you set with the
#       variable 'DIRECTORY'. The program looks for run folders in the
#       'Directory' folder. If there are 'results.xml' files it will attempt
#       to parse them.
#
#       This program also requires the output of the SM2.7 XML Reader program.
#       It uses the output to draw the DNA template amounts that SM
#       calculates for each contributor. From this data this program
#       pulls the lowest DNA template amount to include in the output
#       to use for the charting of the data.
#
# Internally, the program builds a dataframe or table with all the data
#       from the 'results.xml' files that it finds. After exhausting its
#       search of the 'DIRECTORY' folder it will stop and save CSV files
#       in the same folder the program is run from.
#
#       The files that are saved are as follows:
#       First, a file with only results where the LR is greater than 1
#       (supports inclusion) or a True contributor with an LR less than 1
#       (to make sure all true contributors are in the file).
#
#       Second, it will produce a series of files that contain all the
#       results of the DB search. These are split into apparent and
#       true number of contributors (ANOC and TNOC) with the number
#       appended to the end (e.g. ANOC3 or TNOC2). If there are no
#       values for a specific NOC then a file is still made but is empty.
#
#       Note that this program expects to find only True and Apparent
#       number of contributor values that range from 1 to 5 when making
#       the output files. If more than that is needed you can copy and
#       paste a new output filter to include 6 or more NOC.
#

import xml.etree.ElementTree as ET
import pandas as pd
import glob


#############################################################################
# Important variables to set for each run
#       The DIRECTORY variable sets the working directory or folder to look
#       in for STRmix DB search runs.
DIRECTORY = r'DB Results'
#
# The corresponding decon data compiled by the SM2.7 XML Reader.py program.
#       This file should be in the same directory as this program but you
#       can also indicate the full path for the file as in the DIRECTORY
#       of the file as in the DIRECTORY variable above.
DECON_FILE = "LR_data.csv"
#############################################################################


# Import the corresponding decon data to retrieve the lowest template RFU
decon_data = pd.read_csv(DECON_FILE)
decon_data = decon_data.filter(items=["Sample ID", "DNA Amount 1", "DNA Amount 2", "DNA Amount 3",
                        "DNA Amount 4", "DNA Amount 5"])
ex4_data_temp = decon_data[["DNA Amount 1", "DNA Amount 2", "DNA Amount 3",
                        "DNA Amount 4", "DNA Amount 5"]].values.tolist()
temptemp = []
for x in ex4_data_temp:
    try:
        lowest = min(filter(lambda i: i > 0, x))
    except ValueError:
        print("caught")
        lowest = 0
    temptemp.append(lowest)
decon_data["Lowest DNA Amount"] = temptemp
#print(temptemp)
#print(ex4_data)


# Constants: LRCUTOFF sets the lower limit of reported LRs
LRCUTOFF = -1

mixDictionary = {
    'M1-2': ['12M', '14F', '', '', ''],
    'M1-3': ['9F', '21F', '22M', '', ''],
    'M1-4': ['3F', '5F', '7M', '16F', ''],
    'M1-5': ['26M', '27M', '36F', '37F', '38F'],
    'M2-2': ['1M', '22F', '', '', ''],
    'M2-3': ['4M', '8F', '34F', '', ''],
    'M2-4': ['12F', '17F', '23F', '25M', ''],
    'M2-5': ['21M', '28M', '35F', '39F', '40F']
}

duplicateSamplesDictionary = {
    'Mock_17': '16F',
    'Mock_8': '21M',
    'Mock_23_EC': '12F',
    'Mock_26_B': '12F',
    'Mock_4': '8F',
    'Mock_19': '3F',
    'Mock_6': '40F'
}


# These are samples with profiles that match other samples
duplicatesSet = {'Mock_4',
                 'Mock_8',
                 'Mock_17',
                 'Mock_19',
                 'Mock_22_EC',
                 'Mock_22_SF',
                 'Mock_23_EC',
                 'Mock_23_SF',
                 'Mock_24_EC',
                 'Mock_24_SF',
                 'Mock_26_A',
                 'Mock_26_B',
                 'Mock_27_A',
                 'FD_1',
                 'FD_4',
                 'FD_6',
                 'FD_15',
                 'FD_19',
                 'FD_20',
                 'FD_21B',
                 'FD_33B',
                 }


def contributorsList(mix_num, contributors):
    mixCode = mix_num + '-' + contributors
    try:
        return mixDictionary[mixCode]
    except KeyError:
        print("caught")
        return ['', '', '', '', '']



def parseResultsXMLFile(file):
    """Takes an XML file and parses it for values needed
    to construct a table of the various LR values. It
    returns a list of the values of interest in the following
    order: Case Number, Sample ID, Case Notes, Seed, Gelman Rubin,
    AfAm Total LR, Af Am HPD LR, Af Am Unified LR, ...
    for the three remaining NIST populations. If the XML file
    is an LR from previous results file it will place a '0.0' in
    the Gelman Rubin column as that data is not carried into the XML file
    from the previous run."""

    print(f"Current file: {file}")

    tree = ET.parse(file)

    sampleData = []

    lrData = []

    if tree.find('./databaseSearchResults') is not None:
        sampleData.append('DB Search')
    else:
        raise Exception(f"The file: {file} is not a DB search result file.")

    caseNumber = tree.find('.//caseNumber').text
    sampleData.append(caseNumber)

    sampleID = tree.find('.//sampleId').text
    sampleData.append(sampleID)
    mixNumber = sampleID.split('_')[0]
    trueContributorCount = str(sampleID.split('_')[1].count('-') + 1)

    if tree.find('.//caseNotes') is not None:
        caseNotes = tree.find('.//caseNotes').text
        sampleData.append(caseNotes)
    else:
        sampleData.append('')

    seed = tree.find('.//seed').text
    sampleData.append(seed)

    contributors = tree.find('.//contributors').text
    sampleData.append(contributors)
    sampleData.append(trueContributorCount)

    # add in information for contributors but catch case for single source samples
    if 'SS' in mixNumber:
        true_contrib = sampleID.split('_')[-1].split(' ')[0]
        #print(f"True contributor for single source: {true_contrib}")
        contrib_list = [true_contrib, '', '', '', '']
    else:
        contrib_list = contributorsList(mixNumber, trueContributorCount)
    sampleData.extend(contrib_list)

    print(f"Sample ID: {sampleID}")
    lowest_DNA = decon_data.loc[decon_data['Sample ID'] == sampleID]["Lowest DNA Amount"].values[0]

    if tree.find('.//stdResult') is not None:
        resultsList = [[x.attrib["caseNumber"], x.attrib["sample"], x.find('lr').text] for x in tree.findall('.//stdResult')]
        #print(resultsList)
        alreadyAppendedSet = set()
        for result in resultsList:
            temp = []
            # print(result)
            # removing duplicate entries with different names and same names
            # check to see if any of these are needed for other experiments
            if result[1] in duplicatesSet:
                continue
            elif result[1] == '40F':
                continue

            if result[1] == 'Mock_6':
                result[1] = '40F'
            # print(result[1])
            # print(result)

            # If sample name not already in set add to set
            # if in set than already seen duplicate so discard this data
            if result[1] not in alreadyAppendedSet:
                # print("new")
                alreadyAppendedSet.add(result[1])
            else:
                # print('old')
                continue

            if result[1] in contrib_list:
                #print(result[1])
                #print(contrib_list)
                true_or_non = "True"
            else:
                true_or_non = "Non-Contributor"
            #print(template_RFU_amounts)
            #lowest_template_RFU = min(filter(lambda i: i > 0, template_RFU_amounts))
            #print(lowest_template_RFU)
            if float(result[2]) > LRCUTOFF or result[1] in contrib_list:
                temp.extend(sampleData)
                temp.append(true_or_non)
                temp.append(lowest_DNA)
                temp.extend(result)
                lrData.append(temp)


    return lrData


def makeDataFrameAndExport(array):
    """Takes an array of parsed XML data and turns it into
    a Pandas Data Frame for export as a CSV file."""

    df = pd.DataFrame(array,
                      columns=['Run Type', 'Case Number', 'Sample ID', 'Case Notes',
                               'Seed', 'Contributors', 'True Contributors',
                               'Contributor 1', 'Contributor 2', 'Contributor 3',
                               'Contributor 4', 'Contributor 5', 'True/Non', "Template RFU",
                               'Profile Type', 'Profile ID', 'LR'])

    df['Seed'] = df['Seed'].astype(int)
    df['Contributors'] = df['Contributors'].astype(int)
    df['True Contributors'] = df['True Contributors'].astype(int)
    df['LR'] = df['LR'].astype(float)


    # Write all the data to CSV files.
    df.to_csv('DB_Search_LR_data.csv', index=False)

    df_TNoc1 = df[df["True Contributors"] == 1]
    df_TNoc1.to_csv('TNOC_1_DB_Search_LR_data.csv', index=False)
    df_TNoc2 = df[df["True Contributors"] == 2]
    df_TNoc2.to_csv('TNOC_2_DB_Search_LR_data.csv', index=False)
    df_TNoc3 = df[df["True Contributors"] == 3]
    df_TNoc3.to_csv('TNOC_3_DB_Search_LR_data.csv', index=False)
    df_TNoc4 = df[df["True Contributors"] == 4]
    df_TNoc4.to_csv('TNOC_4_DB_Search_LR_data.csv', index=False)
    df_TNoc5 = df[df["True Contributors"] == 5]
    df_TNoc5.to_csv('TNOC_5_DB_Search_LR_data.csv', index=False)

    df_ANoc1 = df[df["Contributors"] == 1]
    df_ANoc1.to_csv('ANOC_1_DB_Search_LR_data.csv', index=False)
    df_ANoc2 = df[df["Contributors"] == 2]
    df_ANoc2.to_csv('ANOC_2_DB_Search_LR_data.csv', index=False)
    df_ANoc3 = df[df["Contributors"] == 3]
    df_ANoc3.to_csv('ANOC_3_DB_Search_LR_data.csv', index=False)
    df_ANoc4 = df[df["Contributors"] == 4]
    df_ANoc4.to_csv('ANOC_4_DB_Search_LR_data.csv', index=False)
    df_ANoc5 = df[df["Contributors"] == 5]
    df_ANoc5.to_csv('ANOC_5_DB_Search_LR_data.csv', index=False)


def main():
    # Using Test path with subfolder names as * wildcard
    # Could specifically refer to file or use * wildcard
    lr_array = []
    for x in glob.glob(fr'{DIRECTORY}\*\results.xml', recursive=True):
        print(f"Current file: {x}")
        temp_array = (parseResultsXMLFile(x))
        lr_array.extend(temp_array)
    #print(lr_array)
    makeDataFrameAndExport(lr_array)


if __name__ == '__main__':
    main()
